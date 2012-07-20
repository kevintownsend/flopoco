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


/* VHDL variable names:
 X, Y: inputs
 XX,YY: after swap

 if(signedIO):
   sX, sY: signs 
   pX, pY: remaining bits, sent to the positive multiplication
*/

/* A big TODO
   For two's complement arithmetic on n bits, the representable interval is [ -2^{n-1}, 2^{n-1}-1 ]
   so the product lives in the interval [-2^{n-1}*2^{n-1}-1,  2^n]
   The value 2^n can only be obtained as the product of the two minimal negative input values
   (the weird ones, which have no symmetric)
   Example on 3 bits: input interval [-4, 3], output interval [-12, 16] and 16 can only be obtained by -4*-4.
   So the output would be representable on 2n-1 bits in two's complement, if it werent for this weird*weird case.

   So for full signed multipliers, we just keep the 2n bits, including one bit used for only for this weird*weird case.
   But for truncated signed multipliers, we should hackingly round down this output to 2^n-1 to avoid carrying around a useless bit.
   This would be a kind of saturated arithmetic I guess.

   Beware, the last bit of Baugh-Wooley tinkering does not need to be added in this case. See the TODO there.

   So the interface must provide a bit that selects this behaviour.

A possible alternative to Baugh-Wooley that solves it:
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


		// interface redundancy
		if(wOut<0 || wXdecl<0 || wYdecl<0) {
			THROWERROR("negative input/output size");
		}


			wFull = wXdecl+wYdecl;

		if(wOut_ > wFull){
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

		// Set up the IO signals
		addInput ( "X"  , wXdecl, true );
		addInput ( "Y"  , wYdecl, true );
		addOutput ( "R"  , wOut, 1 , true );

		// Set up the VHDL library style
		if(signedIO)
			useStdLogicSigned();
		else
			useStdLogicUnsigned();

		// The bit heap
		REPORT(DEBUG,"before bitheap");
		bitHeap = new BitHeap(this, wOut+g);
		REPORT(DEBUG,"after bitheap");
		// Halve number of cases by making sure wY<=wX:
		// interchange x and y in case wY>wX

		if(wY_> wX_){
			wX=wY_;
			wY=wX_;
			vhdl << tab << declare("XX", wX, true) << " <= Y;" << endl; 
			vhdl << tab << declare("YY", wY, true) << " <= X;" << endl; 
		}
		else{
			wX=wX_;
			wY=wY_;
			vhdl << tab << declare("XX", wX, true) << " <= X;" << endl; 
			vhdl << tab << declare("YY", wY, true) << " <= Y;" << endl; 
		}

		// initialize the critical path
		initialCP=getMaxInputDelays ( inputDelays_ );
		setCriticalPath(initialCP);

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
					
					buildLogicOnly();
				}
			}
			else {
			//wx>wxDSP && wy>wyDSP
				
				buildTiling();

			}
		} 

		else {// This target has no DSP, going for a logic-only implementation
			buildLogicOnly();
					
		}
	}






	/**************************************************************************/
	void IntMultiplier::buildLogicOnly() {

		manageSignBeforeMult();
		buildHeapLogicOnly();
		manageSignAfterMult();
		bitHeap -> generateCompressorVHDL();			
		vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;
	}


	/**************************************************************************/
	void IntMultiplier::buildTiling() {
		manageSignBeforeMult();
		buildHeapTiling();
		manageSignAfterMult();
		bitHeap -> generateCompressorVHDL();			
		vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;
	}





#define BAUGHWOOLEY 1

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
			// reminder: wX and wY have been decremented
			vhdl << tab << "-- Managing two's complement with saturated arithmetic" << endl;
			vhdl << tab << declare("rX", wX) << " <= XX" << range(wX-1, 0) << ";" << endl;
			vhdl << tab << declare("rY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;
			// complement X and Y and send them to the positive product
			vhdl << tab << declare("pX", wX) << " <= not rX when sX='1' else " << zg(wX) << ";" << endl;
			vhdl << tab << declare("pY", wX) << " <= not rY when sY='1' else " << zg(wY) << ";" << endl;

			// adding X and Y, possibly complemented
			vhdl << tab << declare("sYpX", wX) << " <= pX when sY='1' else " << zg(wX) << ";" << endl;
			vhdl << tab << declare("sXpY", wY) << " <= pY when sX='1' else " << zg(wY) << ";" << endl;

			// adding only the relevant bits to the bit heap
			for (weight=0; weight<g; weight++) {
				ostringstream v;
				v << "sXrYb(" << wX-g+weight << ")";
				bitHeap->addBit(weight, v.str());
			}
			vhdl << endl;

			vhdl << tab << declare("sXrYb", wY) << " <= not rY when sX='1' else " << zg(wY) << ";" << endl;
			vhdl << tab << "-- Adding these bits to the heap of bits" << endl;
			for (int k=0; k<wX; k++) {
				weight=wY+k - weightShift;
				ostringstream v;
				v << "sYrXb(" << k << ")";
				bitHeap->addBit(weight, v.str());
			}
			vhdl << endl;

			weight=wX+wY - weightShift;
			bitHeap->addBit(weight, "not sX");
			bitHeap->addBit(weight, "not sY");
			// adding sXsY at position wX+wY
			weight=wX+wY - weightShift;
			bitHeap->addBit(weight, "sX and sY");


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

	/**************************************************************************/
	void IntMultiplier::buildHeapLogicOnly() {
		Target *target=getTarget();
		vhdl << tab << "-- code generated by IntMultiplier::buildLogicOnly()"<< endl;

		int dx, dy;				// Here we need to split in small sub-multiplications
		int li=target->lutInputs();
 				
		dx = li>>1;
		dy = li - dx; 
		REPORT(DEBUG, "dx="<< dx << "  dy=" << dy );

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
			vhdl<<tab<<declare("Xp",sizeXPadded)<<" <= "<< zg(padX) << " & pX" << range(wX-1, 0) << ";"<<endl;
			vhdl<<tab<<declare("Yp",sizeYPadded)<<" <= "<< zg(padY) << " & pY" << range(wY-1, 0) << ";"<<endl;	
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++)
				vhdl<<tab<<declare(join("x",k),dx)<<" <= Xp"<<range((k+1)*dx-1,k*dx)<<";"<<endl;
			for (int k=0; k<chunksY ; k++)
				vhdl<<tab<<declare(join("y",k),dy)<<" <= Yp"<<range((k+1)*dy-1, k*dy)<<";"<<endl;
					
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
					vhdl << tab << declare(XY(ix,iy), dx+dy) << " <= y" << iy << " & x" << ix << ";"<<endl;
					inPortMap(t, "X", XY(ix,iy));
					outPortMap(t, "Y", PP(ix,iy));
					vhdl << instance(t, PPTbl(ix,iy));
					vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;
					int maxK=dx+dy; 
					if(ix == chunksX-1)
		 				maxK-=padX;
					if(iy == chunksY-1)
						maxK-=padY;
					for (int k=0; k<maxK; k++) {
						ostringstream s;
						s << PP(ix,iy) << of(k); // right hand side
						int weight = ix*dx+iy*dy+k - weightShift;
						if(weight>=0) {// otherwise these bits deserve to be truncated
							REPORT(DEBUG, "adding bit " << s.str() << " at weight " << weight); 
							bitHeap->addBit(weight, s.str());
						}
					}
					vhdl << endl;
				}
			}				

			if(g>0) {
				int weight=wFull-wOut-1- weightShift;
				bitHeap->addBit(weight, "\'1\'", "The round bit");
			}

			// And that's it, now go compress
		
		}
		

	}
	
	/**************************************************************************/
	void IntMultiplier::buildHeapTiling() {
		vhdl << "-- buildTiling TODO"<< endl;
		int horDSP=wX/wxDSP;
		int verDSP=wY/wyDSP;

	
		REPORT(DEBUG," wX="<< wX << " wxDSP="<< wxDSP <<" wY="<< wY <<" wyDSP="<<wyDSP <<" horDSP= "<<horDSP);
		
		for(int i=0;i<horDSP;i++)
			{
			vhdl << tab << declare(join("XX",i), wxDSP) << " <= pX("<<wFull-(wFull-((i+1)*wxDSP))-1<<" downto "<< wFull-(wFull-(i*wxDSP))<<");"<<endl; 
			}
	
		for(int i=0;i<verDSP;i++)
			{
			vhdl << tab << declare(join("YY",i), wyDSP) << " <= pY("<<wFull-(wFull-((i+1)*wyDSP))-1<<" downto "<< wFull-(wFull-(i*wyDSP)) <<");"<<endl; 
			}
		int k=0;
		for(int i=0;i<horDSP;i++)
			for(int j=0;j<verDSP;j++)
			{
				REPORT(DEBUG, "i="<<i<<" j="<<j);
				REPORT(DEBUG,  "horDSP="<<horDSP<<" verDSP="<<verDSP);

				vhdl << tab << declare (join("XXYY",k),wxDSP+wyDSP)<<" <= "<<join("XX",i)<<" * "<<join("YY",j)<<";"<<endl; 
				for(int l=wxDSP+wyDSP-1;l>=0;l--)
				{
					REPORT(DEBUG, wFull-wOut-g);
					if (l+i*wxDSP+j*wyDSP>=wFull-wOut-g)
					{
						vhdl << tab << declare (join("PP_DSP_X",i,"Y",j,"_",l))<<"<=XXYY"<<k<<"("<<l<<");"<<endl;
						REPORT(DEBUG, "insert here "<<l+(i*wxDSP)+(j*wyDSP)-(wFull-wOut-g));
						bitHeap->addBit(l+(i*wxDSP)+(j*wyDSP)-(wFull-wOut-g),join("PP_DSP_X",i,"Y",j,"_",l));
					}
				}
				k++;		
			}	

		if(g>0) {
				int weight=g-1;
				bitHeap->addBit(weight, "\'1\'", "The round bit");
			}

		
	
	int restX=wX-horDSP*wxDSP;
	int restY=wY-verDSP*wyDSP;
	REPORT(DEBUG,"restX= "<< restX);

	if((restX!=0 ) || (restY!=0))
	
		{
			vhdl << "-- You must add logic too! TODO"<< endl;
			REPORT(DEBUG, "logic needed");
			int dx, dy;				// Here we need to split in small sub-multiplications
			int li=getTarget()->lutInputs();
 			dx = li>>1;
			dy = li - dx; 
			REPORT(DEBUG, "dx="<< dx << "  dy=" << dy );
			//dx,dy=3 for the SmallMultiTable
			/*------------------------------
			 |				|
			 |	first splitting		|
			 |		type		|
			 -----------------------********|
						|second	|
						|split.	|
						|type	|
						|	|
						|	|
						---------
			*/
			
			//FIRST TYPE OF SPLITTING 
			int chunkY1=restY/dy;
			int chunkX1=wX/dx;
			REPORT(DEBUG,"chunkX1= "<<chunkX1);
			
			for (int k=0; k<chunkX1 ; k++)
				{
					vhdl<<tab<<declare(join("x1",k),dx)<<" <= pX"<<range(wX-k*dx-1, wX-(k+1)*dx)<<";"<<endl;
					REPORT(DEBUG,"x1"<<k<<" <= pX"<<wX-k*dx-1<< " downto " <<wX-(k+1)*dx);
				}
			for (int k=0; k<chunkY1 ; k++)
				{
					vhdl<<tab<<declare(join("y1",k),dy)<<" <= pY"<<range(restY-k*dy-1, restY-(k+1)*dy)<<";"<<endl;
					REPORT(DEBUG,"y1"<<k<<" <= pY"<<restY-k*dy-1<< " downto "<< restY-(k+1)*dy);
				}

			//send to SmallMultTable


			//SECOND TYPE OF SPLITTING
			int chunkY2=wY/dy;
			int chunkX2=restX/dx;

			for (int k=0; k<chunkX2 ; k++)
				{
					vhdl<<tab<<declare(join("x2",k),dx)<<" <= pX"<<range(restX-k*dx-1, restX-(k+1)*dx)<<";"<<endl;
					REPORT(DEBUG,"x2"<<k<<" <= pX"<<restX-k*dx-1<< " downto "<< restX-(k+1)*dx);
				}

			for (int k=0; k<chunkY2 ; k++)
				{
					vhdl<<tab<<declare(join("y2",k),dy)<<" <= pY"<<range(wY-k*dy-1, wY-(k+1)*dy)<<";"<<endl;
					REPORT(DEBUG,"y2"<<k<<" <= pY"<<wY-k*dy-1<< " downto "<< wY-(k+1)*dy);
				}
			
			//send to SmallMultTable

			
			
		
		
		}

		//bitHeap->generateCompressorVHDL();
		//vhdl << tab << "R <= CompressionResult(" << wOut+g-1 << " downto "<< g << ");" << endl;



	}


	
	/******************************************************************************/



	
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
	


	string IntMultiplier::PP(int i, int j) {
		std::ostringstream p;
		p << "PPX" << i << "Y" << j;
		return p.str();
	};

	string IntMultiplier::PPTbl( int i, int j) {
		std::ostringstream p;
		p << "PPX" << i << "Y" << j << "_Tbl";
		return p.str();
	};

	string IntMultiplier::XY( int i, int j) {
		std::ostringstream p;
		p  << "Y" << j<< "X" << i;
		return p.str();
	};

	string IntMultiplier::heap( int i, int j) {
		std::ostringstream p;
		p  << "heap_" << i << "_" << j;
		return p.str();
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
			svR++;
			svR &= (mpz_class(1) << (wOut)) -1;
			tc->addExpectedOutput("R", svR);
			
		}
	}
	


	void IntMultiplier::buildStandardTestCases(TestCaseList* tcl){
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


}




