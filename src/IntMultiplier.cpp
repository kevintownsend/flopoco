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
   - pointer to the previous (and the following?) tile if this tile belongs to a supertile

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
   (the weird ones, which have no symmetry)
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

			parentOp->vhdl << tab << declare("XX", wX, true) << " <= " << yname << ";" << endl; 
			parentOp->vhdl << tab << declare("YY", wY, true) << " <= " << xname << ";" << endl; 

		}
		else{
			wX=wXdecl;
			wY=wYdecl;

			parentOp->vhdl << tab << declare("XX", wX, true) << " <= " << xname << ";" << endl; 
			parentOp->vhdl << tab << declare("YY", wY, true) << " <= " << yname << ";" << endl; 

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

			parentOp->vhdl << tab << "-- Ne pouvant me fier à mon raisonnement, j'ai appris par coeur le résultat de toutes les multiplications possibles" << endl;

			SmallMultTable *t = new SmallMultTable( target, wX, wY, wOut, signedIO);
			useSoftRAM(t);
			oplist.push_back(t);

			parentOp->vhdl << tab << declare("XY", wX+wY) << " <= YY & XX;"<<endl;

			inPortMap(t, "X", "XY");
			outPortMap(t, "Y", "RR");

	
			parentOp->vhdl << instance(t, "multTable");
			parentOp->vhdl << tab << "R <= RR;" << endl;

			setCriticalPath(t->getOutputDelay("Y"));

			outDelayMap["R"] = getCriticalPath();
			return;
		}

		// Multiplication by 1-bit integer is simple
		if ((wY == 1)){
			if (signedIO){
				manageCriticalPath( target->localWireDelay(wX) + target->adderDelay(wX+1) );

				parentOp->vhdl << tab << "R <= (" << zg(wX+1)  << " - (XX" << of(wX-1) << " & XX)) when YY(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	

			}
			else {
				manageCriticalPath( target->localWireDelay(wX) + target->lutDelay() );

				parentOp->vhdl << tab << "R <= (\"0\" &  XX) when YY(0)='1' else "<<zg(wX+1,0)<<";"<<endl;	

			}
			outDelayMap["R"] = getCriticalPath();
			return;
		}


		// Multiplication by 2-bit integer is one addition
		if ((wY == 2)) {
			// No need for the following, the mult should be absorbed in the addition
			// manageCriticalPath( target->localWireDelay() + target->lutDelay() );

			parentOp->vhdl << tab << declare("R0",wX+2) << " <= (";

			if (signedIO) 

				parentOp->vhdl << "XX" << of(wX-1) << " & XX" << of(wX-1);  

			else  

				parentOp->vhdl << "\"00\"";
			parentOp->vhdl <<  " & XX) when YY(0)='1' else "<<zg(wX+2,0)<<";"<<endl;	
			parentOp->vhdl << tab << declare("R1i",wX+2) << " <= ";

			if (signedIO) 

				parentOp->vhdl << "(XX" << of(wX-1) << "  & XX & \"0\")";

			else  

				parentOp->vhdl << "(\"0\" & XX & \"0\")";
			parentOp->vhdl << " when YY(1)='1' else " << zg(wX+2,0) << ";"<<endl;	
			parentOp->vhdl << tab << declare("R1",wX+2) << " <= ";

			if (signedIO) 

				parentOp->vhdl << "not R1i;" <<endl;

			else  

				parentOp->vhdl <<  "R1i;" <<endl;

			
			IntAdder *resultAdder = new IntAdder( target, wX+2, inDelayMap("X", target->localWireDelay() + getCriticalPath() ) );
			oplist.push_back(resultAdder);
			
			inPortMap(resultAdder, "X", "R0");
			inPortMap(resultAdder, "Y", "R1");
			inPortMapCst(resultAdder, "Cin", (signedIO? "'1'" : "'0'"));
			outPortMap( resultAdder, "R", "RAdder");

			parentOp->vhdl << tab << instance(resultAdder, "ResultAdder") << endl;

			syncCycleFromSignal("RAdder");
			setCriticalPath( resultAdder->getOutputDelay("R"));

			parentOp->vhdl << tab << "R <= RAdder"<< range(wFull-1, wFull-wOut)<<";"<<endl;	

			outDelayMap["R"] = getCriticalPath();
			return;
		} 


	
		
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

							parentOp->vhdl << tab << declare("rfull", wFull+1) << " <= XX * YY; -- that's one bit more than needed"<<endl; 

						else //sign extension is necessary for using use ieee.std_logic_signed.all; 
							// for correct inference of Xilinx DSP functions

							parentOp->vhdl << tab << declare("rfull", wX + wY + 2) << " <= (\"0\" & XX) * (\"0\" & YY);"<<endl;


						parentOp->vhdl << tab << "R <= rfull"<<range(wFull-1, wFull-wOut)<<";"<<endl;	

						outDelayMap["R"] = getCriticalPath();
					}
				else {
					// For this target and this size, better do a logic-only implementation
				
					buildLogicOnly(0,0,wX,wY);
				}
			}
			else {
			
				buildTiling();

			}
		} 

		else {// This target has no DSP, going for a logic-only implementation
	
			buildLogicOnly(0,0,wX,wY);
		}
		
		
	}



#define BAUGHWOOLEY 1 // the other one doesn't work any better




	/**************************************************************************/
	void IntMultiplier::buildLogicOnly(int topX, int topY, int botX, int botY) {

		manageSignBeforeMult();
		buildHeapLogicOnly(topX,topY,botX,botY);
		//adding the round bit
			if(g>0) {
				int weight=wFull-wOut-1- weightShift;
				bitHeap->addBit(weight, "\'1\'", "The round bit");
			}
		manageSignAfterMult();
		bitHeap -> generateCompressorVHDL();			
#if BAUGHWOOLEY

		parentOp->vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;

#else
		// negate the result if needed

		parentOp->vhdl << tab << declare("RR", wOut) << " <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;
		parentOp->vhdl << tab << declare("sR") << " <= sX xor sY;" << endl;


		parentOp->vhdl << tab << "R <= RR when sR = '0'    else ("<< zg(wOut) << " + (not RR) + '1');" << endl;

#endif

	}


	/**************************************************************************/
	void IntMultiplier::buildTiling() {
		manageSignBeforeMult();
		buildHeapTiling();
		//adding the round bit
			if(g>0) {
				int weight=wFull-wOut-1- weightShift;
				bitHeap->addBit(weight, "\'1\'", "The round bit");
			}
		manageSignAfterMult();

		bitHeap -> generateCompressorVHDL();			

		parentOp->vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;

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

			parentOp->vhdl << tab << declare("sX") << " <= XX" << of(wX) << ";" << endl;
			parentOp->vhdl << tab << declare("sY") << " <= YY" << of(wY) << ";" << endl;

#if BAUGHWOOLEY

			parentOp->vhdl << tab << declare("pX", wX) << " <= XX" << range(wX-1, 0) << ";" << endl;
			parentOp->vhdl << tab << declare("pY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;

			// reminder: wX and wY have been decremented

			parentOp->vhdl << tab << "-- Baugh-Wooley tinkering" << endl;

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

			parentOp->vhdl << tab << declare("sXpYb", wY) << " <= not pY when sX='1' else " << zg(wY) << ";" << endl;
			parentOp->vhdl << tab << "-- Adding these bits to the heap of bits" << endl;

			for (int k=0; k<wX; k++) {
				weight=wY+k - weightShift;
				ostringstream v;
				v << "sYpXb(" << k << ")";
				bitHeap->addBit(weight, v.str());
			}

			parentOp->vhdl << endl;


			setCriticalPath(initialCP);
			manageCriticalPath( getTarget()->localWireDelay(wX) + getTarget()->lutDelay() ) ;  

			parentOp->vhdl << tab << declare("sYpXb", wX) << " <= not pX when sY='1' else " << zg(wX) << ";" << endl;

			for (int k=0; k<wY; k++) {
				weight = wX+k - weightShift;
				ostringstream v;
				v << "sXpYb(" << k << ")";
				bitHeap->addBit(weight, v.str());
			}

			parentOp->vhdl << endl;


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
	  vhdl << tab << declare("rY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;

			parentOp->vhdl << tab << "-- Managing two's complement with saturated arithmetic" << endl;
			parentOp->vhdl << tab << declare("rX", wX) << " <= XX" << range(wX-1, 0) << ";" << endl;
			parentOp->vhdl << tab << declare("rY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;

			  vhdl << tab << declare("pX", wX) << " <= not rX when sX='1' else rX;" << endl;
			  vhdl << tab << declare("pY", wY) << " <= not rY when sY='1' else rY;" << endl;

			parentOp->vhdl << tab << declare("pX", wX) << " <= not rX when sX='1' else rX;" << endl;
			parentOp->vhdl << tab << declare("pY", wY) << " <= not rY when sY='1' else rY;" << endl;


			// adding X and Y, possibly complemented

			parentOp->vhdl << tab << declare("sYpX", wX) << " <= pX when sY='1' else " << zg(wX) << ";" << endl;
			parentOp->vhdl << tab << declare("sXpY", wY) << " <= pY when sX='1' else " << zg(wY) << ";" << endl;


			// adding only the relevant bits to the bit heap
			for (weight=0; weight<wX-weightShift; weight++) {
				ostringstream v;
				v << "sYpX(" << weight+weightShift << ")";
				bitHeap->addBit(weight, v.str());

			parentOp->vhdl << endl;


			for (weight=0; weight<wY-weightShift; weight++) {
				ostringstream v;
				v << "sXpY(" << weight+weightShift << ")";
				bitHeap->addBit(weight, v.str());
			}

			parentOp->vhdl << endl;



#endif
		}	// if(signedIO)
		else	{

			parentOp->vhdl << tab << declare("pX", wX) << " <= XX;" << endl;
			parentOp->vhdl << tab << declare("pY", wY) << " <= YY;" << endl;

		}
	}



	void IntMultiplier::manageSignAfterMult() {
		if (signedIO){
			// restore wX and wY
			wX++;
			wY++;
		}
	}

	/**************************************************************************/
	void IntMultiplier::buildHeapLogicOnly(int topX, int topY, int botX, int botY,int uid) {
		Target *target=getTarget();
        if(uid==-1)
			uid++;
		parentOp->vhdl << tab << "-- code generated by IntMultiplier::buildHeapLogicOnly()"<< endl;


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

			parentOp->vhdl << tab << "-- padding to the left" << endl;
			parentOp->vhdl<<tab<<declare(join("Xp",uid),sizeXPadded)<<" <= "<< zg(padX) << " & pX" << range(botX-1,topX) << ";"<<endl;
			parentOp->vhdl<<tab<<declare(join("Yp",uid),sizeYPadded)<<" <= "<< zg(padY) << " & pY" << range(botY-1, topY) << ";"<<endl;	

			//SPLITTINGS
			for (int k=0; k<chunksX ; k++)
				parentOp->vhdl<<tab<<declare(join("x",uid,"_",k),dx)<<" <= "<<join("Xp",uid)<<range((k+1)*dx-1,k*dx)<<";"<<endl;
				for (int k=0; k<chunksY ; k++)
					parentOp->vhdl<<tab<<declare(join("y",uid,"_",k),dy)<<" <= "<<join("Yp",uid)<<range((k+1)*dy-1, k*dy)<<";"<<endl;
					
			REPORT(DEBUG, "maxWeight=" << maxWeight <<  "    weightShift=" << weightShift);

			SmallMultTable *t = new SmallMultTable( target, dx, dy, dx+dy, false ); // unsigned
			useSoftRAM(t);
			oplist.push_back(t);


			setCycle(0);
			setCriticalPath(initialCP);
			// SmallMultTable is built to cost this:
			manageCriticalPath( getTarget()->localWireDelay(chunksX) + getTarget()->lutDelay() ) ;  
			for (int iy=0; iy<chunksY; iy++){

				parentOp->vhdl << tab << "-- Partial product row number " << iy << endl;

				for (int ix=0; ix<chunksX; ix++) { 

					parentOp->vhdl << tab << declare(XY(ix,iy,uid), dx+dy) << " <= y"<< uid <<"_"<< iy << " & x" << uid <<"_"<< ix << ";"<<endl;

					inPortMap(t, "X", XY(ix,iy,uid));
					outPortMap(t, "Y", PP(ix,iy,uid));

					parentOp->vhdl << instance(t, PPTbl(ix,iy,uid));
					parentOp->vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;

					int maxK=dx+dy; 
					if(ix == chunksX-1)
						maxK-=padX;
					if(iy == chunksY-1)
						maxK-=padY;
					for (int k=0; k<maxK; k++) {
						ostringstream s;
						s << PP(ix,iy,uid) << of(k); // right hand side
						int weight = ix*dx+iy*dy+k - weightShift+topX+topY;
						if(weight>=0) {// otherwise these bits deserve to be truncated
							REPORT(DEBUG, "adding bit " << s.str() << " at weight " << weight); 
							bitHeap->addBit(weight, s.str());
						}
					}

					parentOp->vhdl << endl;

				}
			}				

		

			// And that's it, now go compress
		
		}
	 
	}
	


	
	void IntMultiplier::splitting(int horDSP, int verDSP, int wxDSP, int wyDSP,int restX, int restY)
	{
		int k=0;
		for(int i=0;i<horDSP;i++)
			{
			for(int j=0;j<verDSP;j++)
				{
				
				    bitHeap->addDSP(new MultiplierBlock(parentOp,wxDSP,wyDSP,wX-(i+1)*wxDSP, wY-((j+1)*wyDSP),true,weightShift));
					
					//***** this part is needed only for the plotting ********
					DSP* dsp = new DSP();
					dsp->setTopRightCorner(wX-((i+1)*wxDSP),wY-((j+1)*wyDSP));
					dsp->setBottomLeftCorner(wX-(i*wxDSP),wY-(j*wyDSP));
					dsps.push_back(dsp);
					//********************************************************
			
				}	

			}

	}



	/**************************************************************************/
	void IntMultiplier::buildHeapTiling() {
		
		//the DSPs should be arranged horizontally or vertically?
		int horDSP1=wX/wxDSP;
		int verDSP1=wY/wyDSP;
	    int horDSP2=wX/wyDSP;
		int verDSP2=wY/wxDSP;
		int hor=horDSP1*verDSP1;
		int ver=horDSP2*verDSP2;

		int horDSP;
		int verDSP;
        int restX; //the number of lsbs of the first input which 
		int restY;

		if (hor>=ver)
			{	REPORT(DEBUG, "horizontal");
				horDSP=horDSP1;
				verDSP=verDSP1;
				restX=wX-horDSP*wxDSP;

				restY=wY-verDSP*wyDSP;
				//splitting horizontal
				splitting(horDSP,verDSP,wxDSP,wyDSP,restX,restY);
			}
		else
			{	REPORT(DEBUG, "vertical");
				horDSP=horDSP2;
				verDSP=verDSP2;
				restX=wX-horDSP*wyDSP;
				restY=wY-verDSP*wxDSP;
				//splitting vertical
				splitting(horDSP,verDSP,wyDSP,wxDSP,restX,restY);
			}

		
		if((restX!=0 ) || (restY!=0))
			{

				SmallMultTable *t = new SmallMultTable( getTarget(), 3, 3, 6, false ); // unsigned
				//useSoftRAM(t);
				oplist.push_back(t);
				REPORT(DEBUG,"restX= "<<restX<<"restY= "<<restY);
				if(restY>0)
					{
						buildHeapLogicOnly(0,0,wX,restY,0);
					}
				if(restX>0)
					{
						buildHeapLogicOnly(0,restY,restX,wY,1);
					}	

		
				
		
			}

		bitHeap->iterateDSP();

	}

	

	
	IntMultiplier::~IntMultiplier() {
	}




		string IntMultiplier::PP(int i, int j, int uid ) {
		std::ostringstream p;		
		if(uid==-1) 
			p << "PP_X" << i << "Y" << j;
		else
			p << "PP_"<<uid<<"X" << i << "Y" << j;
		return p.str();
	};

	string IntMultiplier::PPTbl(  int i, int j,int uid) {
		std::ostringstream p;		
		if(uid==-1) 
			p << "PP_X" << i << "Y" << j << "_Tbl";
		else
			p << "PP_"<<uid<<"X" << i << "Y" << j << "_Tbl";
		return p.str();
	};

	string IntMultiplier::XY( int i, int j,int uid) {
		std::ostringstream p;		
		if(uid==-1) 
			p  << "Y" << j<< "X" << i;
		else
			p  << "Y" << j<< "X" << i<<"_"<<uid;
		return p.str();	
	};




	string IntMultiplier::heap( int i, int j) {
		std::ostringstream p;
		p  << "heap_" << i << "_" << j;
		return p.str();
	};




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




