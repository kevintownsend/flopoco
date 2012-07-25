/*
  An integer multiplier mess for FloPoCo

  Authors:  Bogdan Pasca, being cleaned by F de Dinechin, Kinga Illyes and Bogdan Popa

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.
*/

/*
  Interface TODO
  support shared bitheap! In this case,
  - do not compress at the end
  - do not add the round bit
  - import the heap LSB
  - export the truncation error
  - ...
*/


/*
TODO tiling

- define intermediate data struct (list of multiplier blocks)  
multiplier block:
   - a bit saying if it should go into a DSP 
   - x and y size
   - x and y position
   - cycle ?
   - pointer to the previous tile if this tile belongs to a supertile

- a function that explores and builds this structure
   recycle Bogdan's, with optim for large mults
   at least 4 versions: (full, truncated)x(altera, xilinx), please share as much code as possible

- a function that generates VHDL (adding bits to the bit heap)
 */

/* VHDL variable names:
   X, Y: inputs
   XX,YY: after swap

   if(signedIO):
   sX, sY: signs 
   pX, pY: remaining bits, sent to the positive multiplication
*/






/* For two's complement arithmetic on n bits, the representable interval is [ -2^{n-1}, 2^{n-1}-1 ]
   so the product lives in the interval [-2^{n-1}*2^{n-1}-1,  2^n]
   The value 2^n can only be obtained as the product of the two minimal negative input values
   (the weird ones, which have no symmetric)
   Example on 3 bits: input interval [-4, 3], output interval [-12, 16] and 16 can only be obtained by -4*-4.
   So the output would be representable on 2n-1 bits in two's complement, if it werent for this weird*weird case.

   So for full signed multipliers, we just keep the 2n bits, including one bit used for only for this weird*weird case.
   Current situation is: this must be managed from outside. 
   An application that knows that it will not use the full range (e.g. signal processing, poly evaluation) can ignore the MSB bit, 
   but we produce it.



   A big TODO
  
   But for truncated signed multipliers, we should hackingly round down this output to 2^n-1 to avoid carrying around a useless bit.
   This would be a kind of saturated arithmetic I guess.

   Beware, the last bit of Baugh-Wooley tinkering does not need to be added in this case. See the TODO there.

   So the interface must provide a bit that selects this behaviour.

   A possible alternative to Baugh-Wooley that solves it (tried below but it doesn't work, zut alors)
   initially complement (xor) one negative input. This cost nothing, as this xor will be merged in the tables. Big fanout, though.
   then -x=xb+1 so -xy=xb.y+y    
   in case y pos, xy = - ((xb+1)y)  = -(xb.y +y)
   in case x and y neg, xy=(xb+1)(yb+1) = xb.yb +xb +yb +1
   It is enough not to add this lsb 1 to round down the result in this case.
   As this is relevant only to the truncated case, this lsb 1 will indeed be truncated.
   let sx and sy be the signs

   unified equation:

   px = sx xor rx  (on one bit less than x)
   py = sy xor ry

   xy = -1^(sx xor sy)( px.py + px.syb + py.sxb  )   
   (there should be a +sxsy but it is truncated. However, if we add the round bit it will do the same, so the round bit should be sx.sy)
   The final negation is done by complementing again.  
   

   Note that this only applies to truncated multipliers.
   
*/
	

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntMultiplier.hpp"
#include "IntAdder.hpp"
#include "IntMultiAdder.hpp"
#include "IntAddition/NewCompressorTree.hpp"
#include "IntAddition/PopCount.hpp"
// #include "IntMultipliers/SignedIntMultiplier.hpp"
// #include "IntMultipliers/UnsignedIntMultiplier.hpp"
// #include "IntMultipliers/LogicIntMultiplier.hpp"
// #include "IntMultipliers/IntTilingMult.hpp"
// 
using namespace std;

namespace flopoco {


	IntMultiplier::MultiplierBlock::MultiplierBlock(Operator* op, int wX, int wY, int topX, int topY, bool goToDSP,int cycle, MultiplierBlock* previous, MultiplierBlock* next)  :
	op(op), wX(wX), wY(wY), topX(topX), topY(topY)
		{
			if(cycle==-1)
			cycle = op->getCurrentCycle();
		else
			cycle=cycle;
		}



	void IntMultiplier::initialize() {
		// interface redundancy
		if(wOut<0 || wXdecl<0 || wYdecl<0) {
			THROWERROR("negative input/output size");
		}

		wFull = wXdecl+wYdecl;

		if(wOut > wFull){
			THROWERROR("wOut=" << wOut << " too large for " << (signedIO?"signed":"unsigned") << " " << wXdecl << "x" << wYdecl <<  " multiplier." );
		}

		if(wOut==0){ 
			wOut=wFull;
		}

		if(wOut<min(wXdecl, wYdecl))
			REPORT(0, "wOut<min(wX, wY): it would probably be more economical to truncate the inputs first, instead of building this multiplier.");

		wTruncated = wFull - wOut;

		if(wTruncated==0)
			g=0;
		else {
			unsigned i=0;
			mpz_class ulperror=1;
			while( wFull -wOut  > intlog2(ulperror)) {
				REPORT(DEBUG,"i = "<<i<<"  ulperror "<<ulperror<<"  log:"<< intlog2(ulperror) << "  wOut= "<<wOut<< "  wFull= "<<wFull);
				i++;
				ulperror += (i+1)*(mpz_class(1)<<i);
			}
			g=wFull-i-wOut;
			REPORT(DEBUG, "ulp truncation error=" << ulperror << "    g=" << g);
		}

		maxWeight = wOut+g;
		weightShift = wFull - maxWeight;  



		// Halve number of cases by making sure wY<=wX:
		// interchange x and y in case wY>wX

		if(wYdecl> wXdecl){
			wX=wYdecl;
			wY=wXdecl;
			vhdl << tab << declare("XX", wX, true) << " <= " << yname << ";" << endl; 
			vhdl << tab << declare("YY", wY, true) << " <= " << xname << ";" << endl; 
		}
		else{
			wX=wXdecl;
			wY=wYdecl;
			vhdl << tab << declare("XX", wX, true) << " <= " << xname << ";" << endl; 
			vhdl << tab << declare("YY", wY, true) << " <= " << yname << ";" << endl; 
		}

	}




	
	// The virtual constructor
	IntMultiplier::IntMultiplier (Operator* parentOp_, BitHeap* bitHeap_, Signal* x_, Signal* y_, int wX_, int wY_, int wOut_, int lsbWeight_, bool signedIO_, float ratio_):
		Operator ( parentOp_->getTarget()), 
		wXdecl(wX_), wYdecl(wY_), wOut(wOut_), signedIO(signedIO_), ratio(ratio_),  maxError(0.0), parentOp(parentOp_), bitHeap(bitHeap_), x(x_), y(y_) {

		srcFileName="IntMultiplier";
		useDSP = (ratio>0) & getTarget()->hasHardMultipliers();
		
		ostringstream name;
		name <<"VirtualIntMultiplier";
		if(useDSP) 
			name << "UsingDSP_";
		else
			name << "LogicOnly_";
		name << wXdecl << "_" << wYdecl <<"_" << wOut << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
		setName ( name.str() );
		
		xname = x->getName();
		yname = y->getName();
		
		initialize();

	}






	IntMultiplier::IntMultiplier (Target* target, int wX_, int wY_, int wOut_, bool signedIO_, float ratio_,  map<string, double> inputDelays_):
		Operator ( target, inputDelays_ ), wXdecl(wX_), wYdecl(wY_), wOut(wOut_), signedIO(signedIO_), ratio(ratio_), maxError(0.0) {
		srcFileName="IntMultiplier";
		setCopyrightString ( "Florent de Dinechin, Kinga Illyes, Bogdan Popa, Bogdan Pasca, 2012" );
		
		// useDSP or not? 
		useDSP = (ratio>0) & target->hasHardMultipliers();


		{
			ostringstream name;
			name <<"IntMultiplier";
			if(useDSP) 
				name << "UsingDSP_";
			else
				name << "LogicOnly_";
			name << wXdecl << "_" << wYdecl <<"_" << wOut << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
			setName ( name.str() );
		}

		parentOp=this;

		xname="X";
		yname="Y";

		initialize();

		// Set up the IO signals
		addInput ( xname  , wXdecl, true );
		addInput ( yname  , wYdecl, true );
		addOutput ( "R"  , wOut, 1 , true );

		// Set up the VHDL library style
		if(signedIO)
			useStdLogicSigned();
		else
			useStdLogicUnsigned();

		// The bit heap
		bitHeap = new BitHeap(this, wOut+g);
	

		// initialize the critical path
		setCriticalPath(getMaxInputDelays ( inputDelays_ ));

		///////////////////////////////////////
		//  architectures for corner cases   //
		///////////////////////////////////////

		// The really small ones fit in two LUTs and that's as small as it gets  
		if(wX+wY <= target->lutInputs()+2) {
			vhdl << tab << "-- Ne pouvant me fier à mon raisonnement, j'ai appris par coeur le résultat de toutes les multiplications possibles" << endl;
			SmallMultTable *t = new SmallMultTable( target, wX, wY, wOut, signedIO);
			useSoftRAM(t);
			oplist.push_back(t);
			vhdl << tab << declare("XY", wX+wY) << " <= YY & XX;"<<endl;
			inPortMap(t, "X", "XY");
			outPortMap(t, "Y", "RR");
			vhdl << instance(t, "multTable");
			vhdl << tab << "R <= RR;" << endl;
			setCriticalPath(t->getOutputDelay("Y"));

			outDelayMap["R"] = getCriticalPath();
			return;
		}

		// Multiplication by 1-bit integer is simple
		if ((wY == 1)){
			if (signedIO){
				manageCriticalPath( target->localWireDelay(wX) + target->adderDelay(wX+1) );
				vhdl << tab << "R <= (" << zg(wX+1)  << " - (XX" << of(wX-1) << " & XX)) when YY(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	
			}
			else {
				manageCriticalPath( target->localWireDelay(wX) + target->lutDelay() );
				vhdl << tab << "R <= (\"0\" &  XX) when YY(0)='1' else "<<zg(wX+1,0)<<";"<<endl;	
			}
			outDelayMap["R"] = getCriticalPath();
			return;
		}


		// Multiplication by 2-bit integer is one addition
		if ((wY == 2)) {
			// No need for the following, the mult should be absorbed in the addition
			// manageCriticalPath( target->localWireDelay() + target->lutDelay() );
			vhdl << tab << declare("R0",wX+2) << " <= (";
			if (signedIO) 
				vhdl << "XX" << of(wX-1) << " & XX" << of(wX-1);  
			else  
				vhdl << "\"00\"";
			vhdl <<  " & XX) when YY(0)='1' else "<<zg(wX+2,0)<<";"<<endl;	
			vhdl << tab << declare("R1i",wX+2) << " <= ";
			if (signedIO) 
				vhdl << "(XX" << of(wX-1) << "  & XX & \"0\")";
			else  
				vhdl << "(\"0\" & XX & \"0\")";
			vhdl << " when YY(1)='1' else " << zg(wX+2,0) << ";"<<endl;	
			vhdl << tab << declare("R1",wX+2) << " <= ";
			if (signedIO) 
				vhdl << "not R1i;" <<endl;
			else  
				vhdl <<  "R1i;" <<endl;
			
			IntAdder *resultAdder = new IntAdder( target, wX+2, inDelayMap("X", target->localWireDelay() + getCriticalPath() ) );
			oplist.push_back(resultAdder);
			
			inPortMap(resultAdder, "X", "R0");
			inPortMap(resultAdder, "Y", "R1");
			inPortMapCst(resultAdder, "Cin", (signedIO? "'1'" : "'0'"));
			outPortMap( resultAdder, "R", "RAdder");
			vhdl << tab << instance(resultAdder, "ResultAdder") << endl;
			syncCycleFromSignal("RAdder");
			setCriticalPath( resultAdder->getOutputDelay("R"));
			
			vhdl << tab << "R <= RAdder"<< range(wFull-1, wFull-wOut)<<";"<<endl;	
			outDelayMap["R"] = getCriticalPath();
			return;
		} 

#if 0 // needed to restore basic functionality while stuff is being fixed up
		
		if (signedIO)
			vhdl << tab << declare("rfull", wFull+1) << " <= XX * YY; -- that's one bit more than needed"<<endl; 
		else //sign extension is necessary for using use ieee.std_logic_signed.all; 
			// for correct inference of Xilinx DSP functions
			vhdl << tab << declare("rfull", wX + wY + 2) << " <= (\"0\" & XX) * (\"0\" & YY);"<<endl;
						
		vhdl << tab << "R <= rfull"<<range(wFull-1, wFull-wOut)<<";"<<endl;	
		return;
#endif
	
		
		// Now getting more and more generic
		if(useDSP) {
			//test if the multiplication fits into one DSP
			//int wxDSP, wyDSP;
			target->getDSPWidths(wxDSP, wyDSP, signedIO);
			bool testForward, testReverse, testFit;
			testForward     = (wX<=wxDSP)&&(wY<=wyDSP);
			testReverse = (wY<=wxDSP)&&(wX<=wyDSP);
			testFit = testForward || testReverse;
		
			
			REPORT(DEBUG,"useDSP");
			if (testFit){
				if( target->worthUsingDSP(wX, wY))
					{	REPORT(DEBUG,"worthUsingDSP");
						manageCriticalPath(target->DSPMultiplierDelay());
						if (signedIO)
							vhdl << tab << declare("rfull", wFull+1) << " <= XX * YY; -- that's one bit more than needed"<<endl; 
						else //sign extension is necessary for using use ieee.std_logic_signed.all; 
							// for correct inference of Xilinx DSP functions
							vhdl << tab << declare("rfull", wX + wY + 2) << " <= (\"0\" & XX) * (\"0\" & YY);"<<endl;
						
						vhdl << tab << "R <= rfull"<<range(wFull-1, wFull-wOut)<<";"<<endl;	
						outDelayMap["R"] = getCriticalPath();
					}
				else {
					// For this target and this size, better do a logic-only implementation
				
					buildLogicOnly(0,0,wX-1,wY-1);
				}
			}
			else {
				//wx>wxDSP && wy>wyDSP
				
				buildTiling();

			}
		} 

		else {// This target has no DSP, going for a logic-only implementation
	
			buildLogicOnly(0,0,wX-1,wY-1);
		}
		
		
	}



#define BAUGHWOOLEY 1 // the other one doesn't work any better




	/**************************************************************************/
	void IntMultiplier::buildLogicOnly(int topX, int topY, int botX, int botY) {

		manageSignBeforeMult();
		buildHeapLogicOnly(topX,topY,botX,botY);
			if(g>0) {
				int weight=wFull-wOut-1- weightShift;
				bitHeap->addBit(weight, "\'1\'", "The round bit");
			}
		manageSignAfterMult();
		bitHeap -> generateCompressorVHDL();			
#if BAUGHWOOLEY
		vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;
#else
		// negate the result if needed
		vhdl << tab << declare("RR", wOut) << " <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;
		vhdl << tab << declare("sR") << " <= sX xor sY;" << endl;

		vhdl << tab << "R <= RR when sR = '0'    else ("<< zg(wOut) << " + (not RR) + '1');" << endl;
#endif

	}


	/**************************************************************************/
	void IntMultiplier::buildTiling() {
		manageSignBeforeMult();
		buildHeapTiling();
			if(g>0) {
				int weight=wFull-wOut-1- weightShift;
				bitHeap->addBit(weight, "\'1\'", "The round bit");
			}
		manageSignAfterMult();

		bitHeap -> generateCompressorVHDL();			
		vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;
		printConfiguration();
	}



	/**************************************************************************/
	void IntMultiplier::manageSignBeforeMult() {
		int weight;
		if (signedIO){
			// split X as sx and Xr, Y as sy and Xr
			// then XY =  Xr*Yr (unsigned)      + bits and pieces
			// The purpose of the two following lines is to build the heap of bits for Xr*Yr
			// The remaining bits will be added to the heap later.
			wX--;
			wY--;
			vhdl << tab << declare("sX") << " <= XX" << of(wX) << ";" << endl;
			vhdl << tab << declare("sY") << " <= YY" << of(wY) << ";" << endl;
#if BAUGHWOOLEY
			vhdl << tab << declare("pX", wX) << " <= XX" << range(wX-1, 0) << ";" << endl;
			vhdl << tab << declare("pY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;
			// reminder: wX and wY have been decremented
			vhdl << tab << "-- Baugh-Wooley tinkering" << endl;
			setCycle(0);
			setCriticalPath(initialCP);

			// first add all the bits that need no computation
			// adding sX and sY at positions wX and wY
			weight=wX - weightShift;
			bitHeap->addBit(weight, "sX");
			weight=wY - weightShift;
			bitHeap->addBit(weight, "sY");
			// adding the one at position wX+wY+1
			// TODO: do not add this bit if we saturate a truncated mult, as it doesn't belong to the output range...
			weight=wX+wY+1 - weightShift;
			bitHeap->addBit(weight, "'1'");
			// now we have all sort of bits that need some a LUT delay to be computed
			setCriticalPath(initialCP);
			manageCriticalPath( getTarget()->localWireDelay(wY) + getTarget()->lutDelay() ) ;  
			vhdl << tab << declare("sXpYb", wY) << " <= not pY when sX='1' else " << zg(wY) << ";" << endl;
			vhdl << tab << "-- Adding these bits to the heap of bits" << endl;
			for (int k=0; k<wX; k++) {
				weight=wY+k - weightShift;
				ostringstream v;
				v << "sYpXb(" << k << ")";
				bitHeap->addBit(weight, v.str());
			}
			vhdl << endl;

			setCriticalPath(initialCP);
			manageCriticalPath( getTarget()->localWireDelay(wX) + getTarget()->lutDelay() ) ;  
			vhdl << tab << declare("sYpXb", wX) << " <= not pX when sY='1' else " << zg(wX) << ";" << endl;
			for (int k=0; k<wY; k++) {
				weight = wX+k - weightShift;
				ostringstream v;
				v << "sXpYb(" << k << ")";
				bitHeap->addBit(weight, v.str());
			}
			vhdl << endl;

			setCriticalPath(initialCP);
			manageCriticalPath( getTarget()->localWireDelay() + getTarget()->lutDelay() ) ;  
			// adding sXb and sYb at position wX+wY
			weight=wX+wY - weightShift;
			bitHeap->addBit(weight, "not sX");
			bitHeap->addBit(weight, "not sY");
			// adding sXsY at position wX+wY
			weight=wX+wY - weightShift;
			bitHeap->addBit(weight, "sX and sY");

#else // SATURATED VERSION
			// TODO manage pipeline
			// reminder: wX and wY have been decremented
			vhdl << tab << "-- Managing two's complement with saturated arithmetic" << endl;
			vhdl << tab << declare("rX", wX) << " <= XX" << range(wX-1, 0) << ";" << endl;
			vhdl << tab << declare("rY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;
			// complement X and Y and send them to the positive product
			vhdl << tab << declare("pX", wX) << " <= not rX when sX='1' else rX;" << endl;
			vhdl << tab << declare("pY", wY) << " <= not rY when sY='1' else rY;" << endl;

			// adding X and Y, possibly complemented
			vhdl << tab << declare("sYpX", wX) << " <= pX when sY='1' else " << zg(wX) << ";" << endl;
			vhdl << tab << declare("sXpY", wY) << " <= pY when sX='1' else " << zg(wY) << ";" << endl;

			// adding only the relevant bits to the bit heap
			for (weight=0; weight<wX-weightShift; weight++) {
				ostringstream v;
				v << "sYpX(" << weight+weightShift << ")";
				bitHeap->addBit(weight, v.str());
			}
			vhdl << endl;

			for (weight=0; weight<wY-weightShift; weight++) {
				ostringstream v;
				v << "sXpY(" << weight+weightShift << ")";
				bitHeap->addBit(weight, v.str());
			}
			vhdl << endl;


#endif
		}	// if(signedIO)
		else	{
			vhdl << tab << declare("pX", wX) << " <= XX;" << endl;
			vhdl << tab << declare("pY", wY) << " <= YY;" << endl;
		}
	}



	void IntMultiplier::manageSignAfterMult() {
		if (signedIO){
			// restore wX and wY
			wX++;
			wY++;
		}
	}




	void IntMultiplier::buildHeapLogicOnly(MultiplierBlock *mul,int nr) {
		Target *target=getTarget();
		vhdl << tab << "-- code generated by IntMultiplier::buildHeapLogicOnly()"<< endl;

		int dx, dy;				// Here we need to split in small sub-multiplications
		int li=target->lutInputs();
 				
		dx = li>>1;
		dy = li - dx; 
		REPORT(DEBUG, "dx="<< dx << "  dy=" << dy );
		int topX=mul->gettopX();
		int topY=mul->gettopY();
		int botX=topX+wX-1;
		int botY=topY+wY-1;
	
		int wX=botX-topX;
		int wY=botY-topY;
		int chunksX =  int(ceil( ( double(wX) / (double) dx) ));
		int chunksY =  int(ceil( ( 1+ double(wY-dy) / (double) dy) ));
		int sizeXPadded=dx*chunksX; 
		int sizeYPadded=dy*chunksY;
		int padX=sizeXPadded-wX;
		int padY=sizeYPadded-wY;
				
		REPORT(DEBUG, "X split in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
		REPORT(DEBUG, " sizeXPadded="<< sizeXPadded << "  sizeYPadded="<< sizeYPadded);
		if (chunksX + chunksY > 2) { //we do more than 1 subproduct
			vhdl << tab << "-- padding to the left" << endl;
			vhdl<<tab<<declare(join("Xp",nr),sizeXPadded)<<" <= "<< zg(padX) << " & pX" << range(botX-1,topX) << ";"<<endl;
			vhdl<<tab<<declare(join("Yp",nr),sizeYPadded)<<" <= "<< zg(padY) << " & pY" << range(botY-1, topY) << ";"<<endl;	
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++)
				vhdl<<tab<<declare(join("x",nr,"_",k),dx)<<" <= "<<join("Xp",nr)<<range((k+1)*dx-1,k*dx)<<";"<<endl;
			for (int k=0; k<chunksY ; k++)
				vhdl<<tab<<declare(join("y",nr,"_",k),dy)<<" <= "<<join("Yp",nr)<<range((k+1)*dy-1, k*dy)<<";"<<endl;
					
			REPORT(DEBUG, "maxWeight=" << maxWeight <<  "    weightShift=" << weightShift);

			SmallMultTable *t = new SmallMultTable( target, dx, dy, dx+dy, false ); // unsigned
			useSoftRAM(t);
			oplist.push_back(t);


			setCycle(0);
			setCriticalPath(initialCP);
			// SmallMultTable is built to cost this:
			manageCriticalPath( getTarget()->localWireDelay(chunksX) + getTarget()->lutDelay() ) ;  
			for (int iy=0; iy<chunksY; iy++){
				vhdl << tab << "-- Partial product row number " << iy << endl;
				for (int ix=0; ix<chunksX; ix++) { 
					vhdl << tab << declare(XY(ix,iy,nr), dx+dy) << " <= y"<< nr <<"_"<< iy << " & x" << nr <<"_"<< ix << ";"<<endl;
					inPortMap(t, "X", XY(ix,iy,nr));
					outPortMap(t, "Y", PP(ix,iy,nr));
					vhdl << instance(t, PPTbl(ix,iy,nr));
					vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;
					int maxK=dx+dy; 
					if(ix == chunksX-1)
						maxK-=padX;
					if(iy == chunksY-1)
						maxK-=padY;
					for (int k=0; k<maxK; k++) {
						ostringstream s;
						s << PP(ix,iy,nr) << of(k); // right hand side
						int weight = ix*dx+iy*dy+k - weightShift+topX+topY;
						if(weight>=0) {// otherwise these bits deserve to be truncated
							REPORT(DEBUG, "adding bit " << s.str() << " at weight " << weight); 
							bitHeap->addBit(weight, s.str());
						}
					}
					vhdl << endl;
				}
			}				

		

			// And that's it, now go compress
		
		}
	 
	}
	








	/**************************************************************************/
	void IntMultiplier::buildHeapLogicOnly(int topX, int topY, int botX, int botY,int nr) {
		Target *target=getTarget();
		vhdl << tab << "-- code generated by IntMultiplier::buildHeapLogicOnly()"<< endl;

		int dx, dy;				// Here we need to split in small sub-multiplications
		int li=target->lutInputs();
 				
		dx = li>>1;
		dy = li - dx; 
		REPORT(DEBUG, "dx="<< dx << "  dy=" << dy );


		int wX=botX-topX;
		int wY=botY-topY;
		int chunksX =  int(ceil( ( double(wX) / (double) dx) ));
		int chunksY =  int(ceil( ( 1+ double(wY-dy) / (double) dy) ));
		int sizeXPadded=dx*chunksX; 
		int sizeYPadded=dy*chunksY;
		int padX=sizeXPadded-wX;
		int padY=sizeYPadded-wY;
				
		REPORT(DEBUG, "X split in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
		REPORT(DEBUG, " sizeXPadded="<< sizeXPadded << "  sizeYPadded="<< sizeYPadded);
		if (chunksX + chunksY > 2) { //we do more than 1 subproduct
			vhdl << tab << "-- padding to the left" << endl;
			vhdl<<tab<<declare(join("Xp",nr),sizeXPadded)<<" <= "<< zg(padX) << " & pX" << range(botX-1,topX) << ";"<<endl;
			vhdl<<tab<<declare(join("Yp",nr),sizeYPadded)<<" <= "<< zg(padY) << " & pY" << range(botY-1, topY) << ";"<<endl;	
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++)
				vhdl<<tab<<declare(join("x",nr,"_",k),dx)<<" <= "<<join("Xp",nr)<<range((k+1)*dx-1,k*dx)<<";"<<endl;
			for (int k=0; k<chunksY ; k++)
				vhdl<<tab<<declare(join("y",nr,"_",k),dy)<<" <= "<<join("Yp",nr)<<range((k+1)*dy-1, k*dy)<<";"<<endl;
					
			REPORT(DEBUG, "maxWeight=" << maxWeight <<  "    weightShift=" << weightShift);

			SmallMultTable *t = new SmallMultTable( target, dx, dy, dx+dy, false ); // unsigned
			useSoftRAM(t);
			oplist.push_back(t);


			setCycle(0);
			setCriticalPath(initialCP);
			// SmallMultTable is built to cost this:
			manageCriticalPath( getTarget()->localWireDelay(chunksX) + getTarget()->lutDelay() ) ;  
			for (int iy=0; iy<chunksY; iy++){
				vhdl << tab << "-- Partial product row number " << iy << endl;
				for (int ix=0; ix<chunksX; ix++) { 
					vhdl << tab << declare(XY(ix,iy,nr), dx+dy) << " <= y"<< nr <<"_"<< iy << " & x" << nr <<"_"<< ix << ";"<<endl;
					inPortMap(t, "X", XY(ix,iy,nr));
					outPortMap(t, "Y", PP(ix,iy,nr));
					vhdl << instance(t, PPTbl(ix,iy,nr));
					vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;
					int maxK=dx+dy; 
					if(ix == chunksX-1)
						maxK-=padX;
					if(iy == chunksY-1)
						maxK-=padY;
					for (int k=0; k<maxK; k++) {
						ostringstream s;
						s << PP(ix,iy,nr) << of(k); // right hand side
						int weight = ix*dx+iy*dy+k - weightShift+topX+topY;
						if(weight>=0) {// otherwise these bits deserve to be truncated
							REPORT(DEBUG, "adding bit " << s.str() << " at weight " << weight); 
							bitHeap->addBit(weight, s.str());
						}
					}
					vhdl << endl;
				}
			}				

		

			// And that's it, now go compress
		
		}
	 
	}
	

	void IntMultiplier::splitting(int horDSP, int verDSP, int wxDSP, int wyDSP,int restX, int restY)
	{
		for(int i=0;i<horDSP;i++)
			{
				vhdl << tab << declare(join("XX",horDSP-i-1), wxDSP) << " <= XX("<<wX-(i*wxDSP)-1<<" downto "<< wX-((i+1)*wxDSP)<<");"<<endl; 
				REPORT(DEBUG,join("XX",horDSP-i-1)<< " <= XX("<<wX-(i*wxDSP)-1<<" downto "<< wX-((i+1)*wxDSP)<<");");
			}
	
		for(int i=0;i<verDSP;i++)
			{
				REPORT(DEBUG,join("YY",verDSP-i-1)<< " <= YY("<<wY-((i)*wyDSP)-1<<" downto "<< wY-(((i+1)*wyDSP)) <<");");
				vhdl << tab << declare(join("YY",verDSP-i-1), wyDSP) << " <= YY("<<(wY-(i*wyDSP))-1<<" downto "<< wY-((i+1)*wyDSP) <<");"<<endl; 
			}

				int k=0;
		for(int i=0;i<horDSP;i++)
			for(int j=0;j<verDSP;j++)
				{
					REPORT(DEBUG, "i="<<i<<" j="<<j);
					REPORT(DEBUG, "XX"<<i<<"YY"<<j<<" <= "<<join("XX",i)<<" * "<<join("YY",j)<<";");
					
					DSP* dsp = new DSP();

					dsp->setTopRightCorner(wX-((i+1)*wxDSP),wY-((j+1)*wyDSP));
					dsp->setBottomLeftCorner(wX-(i*wxDSP),wY-(j*wyDSP));
				
					dsps.push_back(dsp);

					vhdl << tab << declare (join("XX",i,"YY",j),wxDSP+wyDSP)<<" <= "<<join("XX",i)<<" * "<<join("YY",j)<<";"<<endl; 
					for(int l=wxDSP+wyDSP-1;l>=0;l--)
						{
							//	REPORT(DEBUG, wFull-wOut-g);
							if (l+i*wxDSP+j*wyDSP+restX+restY>=wFull-wOut-g)
								{
									vhdl << tab << declare (join("PP_DSP_X",i,"Y",j,"_",l))<<"<=XX"<<i<<"YY"<<j<<"("<<l<<");"<<endl;
									REPORT(DEBUG, "insert here XX"<<i<<"YY"<<j<<"("<<l<<") in column "<<l+(i*wxDSP)+(j*wyDSP)+restX+restY-(wFull-wOut-g));
									bitHeap->addBit(l+(i*wxDSP)+(j*wyDSP)+restX+restY-(wFull-wOut-g),join("PP_DSP_X",i,"Y",j,"_",l));
								}
						}
					k++;		
				}	



	}



	/**************************************************************************/
	void IntMultiplier::buildHeapTiling() {
		
		//horizontal or vertical DSP?
		int horDSP1=wX/wxDSP;
		int verDSP1=wY/wyDSP;
	    int horDSP2=wX/wyDSP;
		int verDSP2=wY/wxDSP;
		int hor=horDSP1*verDSP1;
		int ver=horDSP2*verDSP2;

		int horDSP;
		int verDSP;
        int restX;
		int restY;

		if (hor>=ver)
			{	REPORT(DEBUG, "horizontal");
				horDSP=horDSP1;
				verDSP=verDSP1;
				restX=wX-horDSP*wxDSP;
				restY=wY-verDSP*wyDSP;
				splitting(horDSP,verDSP,wxDSP,wyDSP,restX,restY);
			}
		else
			{	REPORT(DEBUG, "vertical");
				horDSP=horDSP2;
				verDSP=verDSP2;
				restX=wX-horDSP*wyDSP;
				restY=wY-verDSP*wxDSP;
				splitting(horDSP,verDSP,wyDSP,wxDSP,restX,restY);
			}

	
		REPORT(DEBUG,"restX= "<< restX);
		REPORT(DEBUG,"restY= "<< restY);
		REPORT(DEBUG,"horDSP= "<< horDSP);
		REPORT(DEBUG,"verDSP= "<< verDSP);
		REPORT(DEBUG," wX="<< wX << " wxDSP="<< wxDSP <<" wY="<< wY <<" wyDSP="<<wyDSP <<" horDSP= "<<horDSP);
		
		if((restX!=0 ) || (restY!=0))
			{

				SmallMultTable *t = new SmallMultTable( getTarget(), 3, 3, 6, false ); // unsigned
				//useSoftRAM(t);
				oplist.push_back(t);
				REPORT(DEBUG,"restX= "<<restX<<"restY= "<<restY);
				if(restY>0)
					{
						REPORT(DEBUG,"first small tiling");
						buildHeapLogicOnly(0,0,wX,restY,0);
					}
				if(restX>0)
					{
						REPORT(DEBUG, "second small tiling");
						buildHeapLogicOnly(0,restY,restX,wY,1);
					}	

		
				
		
			}
	}

	

	
	IntMultiplier::~IntMultiplier() {
	}



	// Stuff for the small multiplier table 
	IntMultiplier::SmallMultTable::SmallMultTable(Target* target, int dx, int dy, int wO, bool  signedIO) : 
		Table(target, dx+dy, wO, 0, -1, true), // logic table
		dx(dx), dy(dy), signedIO(signedIO) {
		ostringstream name; 
		srcFileName="LogicIntMultiplier::SmallMultTable";
		name <<"SmallMultTable" << dy << "x" << dx << "r" << wO << (signedIO?"signed":"unsigned");
		setName(name.str());				
	};
	
	mpz_class IntMultiplier::SmallMultTable::function(int yx){
		mpz_class r;
		int y = yx>>dx;
		int x = yx -(y<<dx);
		int wF=dx+dy;
		if(signedIO) wF--;

		if(signedIO){
			if ( x >= (1 << (dx-1)))
				x -= (1 << dx);
			if ( y >= (1 << (dy-1)))
				y -= (1 << dy);
			r = x * y;
			if ( r < 0)
				r += mpz_class(1) << wF; 
		}
		else 
			r = x*y;

		if(wOut<wF){ // wOut is that of Table
			// round to nearest, but not to nearest even
			int tr=wF-wOut; // number of truncated bits 
			// adding the round bit at half-ulp position
			r += (mpz_class(1) << (tr-1));
			r = r >> tr;
		}

		return r;
		
	};
	


	string IntMultiplier::PP(int i, int j, int nr ) {
		std::ostringstream p;		
		if(nr==-1) 
			p << "PP_X" << i << "Y" << j;
		else
			p << "PP_"<<nr<<"X" << i << "Y" << j;
		return p.str();
	};

	string IntMultiplier::PPTbl(  int i, int j,int nr) {
		std::ostringstream p;		
		if(nr==-1) 
			p << "PP_X" << i << "Y" << j << "_Tbl";
		else
			p << "PP_"<<nr<<"X" << i << "Y" << j << "_Tbl";
		return p.str();
	};

	string IntMultiplier::XY( int i, int j,int nr) {
		std::ostringstream p;		
		if(nr==-1) 
			p  << "Y" << j<< "X" << i;
		else
			p  << "Y" << j<< "X" << i<<"_"<<nr;
		return p.str();	
	};




	string IntMultiplier::heap( int i, int j) {
		std::ostringstream p;
		p  << "heap_" << i << "_" << j;
		return p.str();
	};



	void IntMultiplier::exploration(MultiplierBlock *mul)
	{
	 //splits the big block, and adds the roots of the supertilings to mulBlocks
	}


    void IntMultiplier::generateVHDL(MultiplierBlock *m, int nr,int i)
	{
	stringstream s;
	int topX=m->gettopX();
	int topY=m->gettopY();
	int botX=topX+m->getwX()-1;
	int botY=topY+m->getwY()-1;
	

	vhdl << tab << declare(join("DSP",i,"_",nr), m->getwX()+m->getwY()) << " <= XX"<<range(botX,topX)<<" * YY"<<range(botY,topY)<<";"<<endl;
	s<<join("DSP",i,"_",nr);
	
	m->setSignalName(s.str());
	m->setSignalLength(m->getwX()+m->getwY());
	}


	void IntMultiplier::chaining(MultiplierBlock *root,int i)
	{
		int nr=0;
		MultiplierBlock* mul=root;
		string prevName="";
		int prevLength=0;
		int rootLength=root->getSigLength();
		string rootName=root->getSigName();
		int newLength=rootLength;
		while(mul->getNext()!=NULL)
		{
		
		generateVHDL(mul,nr,i);
		string rootName=mul->getSigName();
		int rootLength=mul->getSigLength();
		int newLength=rootLength;
        prevName="";
		prevLength=0;

		if(mul->getPrevious()!=NULL)
		{
			 prevName=mul->getPrevious()->getSigName();
			int prevLength=mul->getPrevious()->getSigLength();

			if(rootLength<prevLength)
				newLength=prevLength;
		}
		
		
		newLength=newLength+1;
		vhdl<< tab << declare(join("DSP_ch",i,"_",nr),newLength) << prevName << " + " <<rootName<<";"<<endl;
		mul->setSignalLength(newLength);
		mul->setSignalName(join("DSP_ch_",i,"_",nr));
		nr++;
		mul=mul->getNext();
		}
		

		
		if(mul->getPrevious()!=NULL)
		{
			 prevName=mul->getPrevious()->getSigName();
			int prevLength=mul->getPrevious()->getSigLength();

			if(rootLength<prevLength)
				newLength=prevLength;
		}

		vhdl<< tab << declare(join("DSP_ch",i,"_",nr),newLength) << prevName << " + " <<rootName<<";"<<endl;
		mul->setSignalLength(newLength);
		mul->setSignalName(join("DSP_ch_",i,"_",nr));

		for(int k=newLength;k>=0;k--)
		{
			int weight=mul->gettopX()+mul->gettopY()+k-weightShift;
			if (weight>=0)
			{
				stringstream st;
				st<<join("DSP_ch_",i,"_",nr)<<"("<<k<<")";
				bitHeap->addBit(weight, st.str());
			}
		}



			
	}





	void IntMultiplier::iterate()
	{
		for(int i=0; i<mulBlocks.size();i++)
			{
				if(!mulBlocks[i]->getgoToDSP())
				{
					if((wX>wxDSP && wY>wyDSP)||(wX>wyDSP && wY>wxDSP))
						exploration(mulBlocks[i]);
					else buildHeapLogicOnly(mulBlocks[i],i);
				}
				else
				chaining(mulBlocks[i],i);
			}
		
	}







	void IntMultiplier::emulate ( TestCase* tc ) {
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svR;
		
		if (! signedIO){
			svR = svX * svY;
		}

		else{ // Manage signed digits
			mpz_class big1 = (mpz_class(1) << (wXdecl));
			mpz_class big1P = (mpz_class(1) << (wXdecl-1));
			mpz_class big2 = (mpz_class(1) << (wYdecl));
			mpz_class big2P = (mpz_class(1) << (wYdecl-1));

			if ( svX >= big1P)
				svX -= big1;

			if ( svY >= big2P)
				svY -= big2;
			
			svR = svX * svY;
			if ( svR < 0){
				svR += (mpz_class(1) << wFull); 
			}

		}
		if(wTruncated==0) 
			tc->addExpectedOutput("R", svR);
		else {
			// there is truncation, so this mult should be faithful
			svR = svR >> wTruncated;
			tc->addExpectedOutput("R", svR);
#if BAUGHWOOLEY
			svR++;
			svR &= (mpz_class(1) << (wOut)) -1;
			tc->addExpectedOutput("R", svR);
#else // saturate
			if(svR<(mpz_class(1) << wOut)-1) {
				svR++;
				tc->addExpectedOutput("R", svR);
			}
#endif			
		}
	}
	


	void IntMultiplier::buildStandardTestCases(TestCaseList* tcl)
    {
		TestCase *tc;

		mpz_class x, y;

		// 1*1
		x = mpz_class(1); 
		y = mpz_class(1); 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
		
		// -1 * -1
		x = (mpz_class(1) << wXdecl) -1; 
		y = (mpz_class(1) << wYdecl) -1; 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// The product of the two max negative values overflows the signed multiplier
		x = mpz_class(1) << (wXdecl -1); 
		y = mpz_class(1) << (wYdecl -1); 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}

    void IntMultiplier::printConfiguration()
    {
        printAreaView();
        printLozengeView();
    }

    void IntMultiplier::printAreaView()
    {
       	ostringstream figureFileName;
		figureFileName << "view_area_" << "DSP"<< ".svg";
		
		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fig.open (figureFileName.str().c_str(), ios::trunc);


        //scaling factor for the whole drawing
        int scalingFactor = 5;

        //offsets for the X and Y axes
		int offsetX = 180;
		int offsetY = 120;

		//file header
		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;

	    //draw target rectangle	
        drawTargetFigure(wX, wY, offsetX, offsetY, scalingFactor, true);

        //draw DSPs
		int xT, yT, xB, yB;

		for(unsigned i=0; i<dsps.size(); i++)
		{
			dsps[i]->getCoordinates(xT, yT, xB, yB);
			drawDSP(i, xT, yT, xB, yB, offsetX, offsetY, scalingFactor, true);
		}


        //draw truncation line
		if(wFull-wOut > 0)
		{
            drawLine(wX, wY, wOut, offsetX, offsetY, scalingFactor, true);    
		}

        //draw guard line
		if(g>0)
		{
			drawLine(wX, wY, wOut+g, offsetX, offsetY, scalingFactor, true);
		}


		fig << "</svg>" << endl;

		fig.close();
        
    }



    void IntMultiplier::drawTargetFigure(int wX, int wY, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
    {
        if(isRectangle)
       	    fig << "<rect x=\"" << offsetX << "\" y=\"" << offsetY << "\" height=\"" << wY * scalingFactor << "\" width=\"" << wX * scalingFactor 
			    <<"\" style=\"fill:rgb(255, 255, 255);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;
        else
    		fig << "<polygon points=\"" << offsetX << "," << offsetY << " " 
				<< wX*scalingFactor + offsetX << "," << offsetY << " " 
				<< wX*scalingFactor + offsetX - scalingFactor*wY << "," << wY*scalingFactor + offsetY << " "
			    << offsetX - scalingFactor*wY << "," << wY*scalingFactor + offsetY 	
				<< "\" style=\"fill:rgb(255, 255, 255);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;

        REPORT(DETAILED, "wX " << wX << "   wY " << wY);
    }

    void IntMultiplier::drawLine(int wX, int wY, int wRez, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
    {
        if(isRectangle)
	    	fig << "<line x1=\"" << offsetX + scalingFactor*(wRez - wX) << "\" y1=\"" << offsetY << "\" x2=\"" << offsetX + scalingFactor*wRez
                << "\" y2=\"" << offsetY + scalingFactor*wY <<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>" << endl ;	
        else
	    	fig << "<line x1=\"" << offsetX + scalingFactor*(wRez - wY) << "\" y1=\"" << offsetY << "\" x2=\"" << offsetX + scalingFactor*(wRez - wY)
                << "\" y2=\"" << offsetY + scalingFactor*wY <<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>" << endl ;	
    
    }

    


    void IntMultiplier::printLozengeView()
    {
       	ostringstream figureFileName;
		figureFileName << "view_lozenge_" << "DSP"<< ".svg";
		
		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fig.open (figureFileName.str().c_str(), ios::trunc);


        //scaling factor for the whole drawing
        int scalingFactor = 5;

        //offsets for the X and Y axes
		int offsetX = 180 + wY*scalingFactor;
		int offsetY = 120;

		//file header
		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;

	    //draw target lozenge	
        drawTargetFigure(wX, wY, offsetX, offsetY, scalingFactor, false);

        //draw DSPs
		int xT, yT, xB, yB;

		for(unsigned i=0; i<dsps.size(); i++)
		{
			dsps[i]->getCoordinates(xT, yT, xB, yB);
			drawDSP(i, xT, yT, xB, yB, offsetX, offsetY, scalingFactor, false);
		}


        //draw truncation line
		if(wFull-wOut > 0)
		{
            drawLine(wX, wY, wOut, offsetX, offsetY, scalingFactor, false);    
		}

        //draw guard line
		if(g>0)
		{
			drawLine(wX, wY, wOut+g, offsetX, offsetY, scalingFactor, false);
		}


		fig << "</svg>" << endl;

		fig.close();
        
    }
	

	
	void IntMultiplier::drawDSP(int i, int xT, int yT, int xB, int yB, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
	{
			
        //because the X axis is opposing, all X coordinates have to be turned around
	    int turnaroundX;

        if(isRectangle)
        {
            turnaroundX = offsetX + wX * scalingFactor;
	    	fig << "<rect x=\"" << turnaroundX - xB*scalingFactor << "\" y=\"" << yT*scalingFactor + offsetY 
                << "\" height=\"" << (yB-yT)*scalingFactor << "\" width=\"" << (xB-xT)*scalingFactor
		       	<< "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;
		    fig << "<text x=\"" << (2*turnaroundX - scalingFactor*(xT+xB))/2 -12 << "\" y=\"" << ((yT+yB)*scalingFactor)/2 + offsetY + 7
				<< "\" fill=\"blue\">D" <<  xT / wxDSP <<  yT / wyDSP  << "</text>" << endl;
        }   
        else 
        {
            turnaroundX = wX * scalingFactor;
      		fig << "<polygon points=\"" << turnaroundX - 5*xB + offsetX - 5*yT << "," << 5*yT + offsetY << " "
				<< turnaroundX - 5*xT + offsetX - 5*yT << "," << 5*yT + offsetY << " " 
				<< turnaroundX - 5*xT + offsetX - 5*yB << "," << 5*yB + offsetY << " "
				<< turnaroundX - 5*xB + offsetX - 5*yB << "," << 5*yB + offsetY
                << "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;

            fig << "<text x=\"" << (2*turnaroundX - xB*5 - xT*5 + 2*offsetX)/2 - 14 - (yT*5 + yB*5)/2 
                << "\" y=\"" << ( yT*5 + offsetY + yB*5 + offsetY )/2 + 7 
				<< "\" fill=\"blue\">D" <<  xT / wxDSP <<  yT / wyDSP  << "</text>" << endl;	
	
        }
	}

 }




