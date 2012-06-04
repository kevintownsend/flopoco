/*
An integer multiplier mess for FloPoCo

Authors:  Bogdan Pasca, being cleaned by F de Dinechin

This file is part of the FloPoCo project
developed by the Arenaire team at Ecole Normale Superieure de Lyon

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
2008-2010.
  All rights reserved.
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
#include "../IntAdder.hpp"
#include "../IntMultiAdder.hpp"
#include "../IntAddition/NewCompressorTree.hpp"
// #include "IntMultipliers/SignedIntMultiplier.hpp"
// #include "IntMultipliers/UnsignedIntMultiplier.hpp"
// #include "IntMultipliers/LogicIntMultiplier.hpp"
// #include "IntMultipliers/IntTilingMult.hpp"
// 
using namespace std;

namespace flopoco {
	extern vector<Operator*> oplist;



	
	IntMultiplier::IntMultiplier (Target* target, int wX_, int wY_, int wOut_, bool signedInputs_, float ratio_, bool useDSP_,  map<string, double> inputDelays_):
		Operator ( target, inputDelays_ ), wX(wX_), wY(wY_), wOut(wOut_), signedInputs(signedInputs_), ratio(ratio_), useDSP(useDSP_) {
		srcFileName="IntMultiplier";
		setCopyrightString ( "Florent de Dinechin 2012" );
		
		// useDSP or not? 
		useDSP = useDSP & target->hasHardMultipliers();

		// interface redundancy
		if(wOut<0 || wX<0 || wY<0) {
			THROWERROR("negative input/output size");
		}
		if(wOut>wX+wY){
			THROWERROR("wOut>wX+wY ("<< wOut << ">" << wX <<"+"<<wY<<")" );
		}

		if(wOut==0) 
			wOut=wX+wY;
		wTruncated = wX+wY - wOut;

		if(wTruncated==0)
			g=0;
		else {
			int i=0;
			mpz_class ulperror=1;
			while(wOut<wX+wY-intlog2(ulperror)) {
				i++;
				ulperror += (i+1)*(1<<i);
			}
			g=wX+wY-i-wOut;
			REPORT(DEBUG, "ulp truncation error=" << ulperror << "    g=" << g);
		}
		{
			ostringstream name;
			name <<"IntMultiplier";
			if(useDSP) 
				name << "UsingDSP_";
			else
				name << "LogicOnly_";
			name << wX << "_" << wY <<"_" << wOut << "_" << (signedInputs?"signed":"unsigned") << "_uid"<<Operator::getNewUId();;
			setName ( name.str() );
		}
		// Set up the IO signals
		addInput ( "X"  , wX, true );
		addInput ( "Y"  , wY, true );
		addOutput ( "R"  , wOut, 1 , true );

		// Set up the VHDL library style
		if(signedInputs)
			useStdLogicSigned();
		else
			useStdLogicUnsigned();

		// Halve number of cases by assuming wY<=wX:
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
		setCriticalPath(getMaxInputDelays ( inputDelays_ ));

		//architectures for corner cases 
		if ((wY == 1)){
			if (signedInputs){
				manageCriticalPath( target->localWireDelay(wX) + target->adderDelay(wX+1) );
				vhdl << tab << "R <= (" << zg(wX+1)  << " - (XX" << of(wX-1) << " & XX)) when YY(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	
			}
			else {
				manageCriticalPath( target->localWireDelay(wX) + target->lutDelay() );
				vhdl << tab << "R <= (\"0\" &  XX) when YY(0)='1' else "<<zg(wX+1,0)<<";"<<endl;	
			}
				outDelayMap["R"] = getCriticalPath();
		}
		else if ((wY == 2)) {
			if(target->lutInputs()-2>=wX){
				if(signedInputs) {
					std::ostringstream o;
					o << srcFileName << " (" << uniqueName_ << "): SORRY : ultrasmall _signed_IntMultiplier not implemented";
					throw o.str();
				} // FIXME of course
				manageCriticalPath( target->localWireDelay() + target->lutDelay() );
				vhdl << tab << "R <= XX * YY;  -- fits in one look-up table"<<endl;	
				outDelayMap["R"] = getCriticalPath();
			}
			
			else { // wY=2, one addition
				// No need for the following, the mult should be absorbed in the addition
				// manageCriticalPath( target->localWireDelay() + target->lutDelay() );
				vhdl << tab << declare("R0",wX+2) << " <= (";
				if (signedInputs) 
					vhdl << "XX" << of(wX-1) << " & XX" << of(wX-1);  
				else  
					vhdl << "\"00\"";
				vhdl <<  " & XX) when YY(0)='1' else "<<zg(wX+2,0)<<";"<<endl;	
				vhdl << tab << declare("R1i",wX+2) << " <= ";
				if (signedInputs) 
					vhdl << "(XX" << of(wX-1) << "  & XX & \"0\")";
				else  
					vhdl << "(\"0\" & XX & \"0\")";
				vhdl << " when YY(1)='1' else " << zg(wX+2,0) << ";"<<endl;	
				vhdl << tab << declare("R1",wX+2) << " <= ";
				if (signedInputs) 
					vhdl << "not R1i;" <<endl;
				else  
					vhdl <<  "R1i;" <<endl;
				
				IntAdder *resultAdder = new IntAdder( target, wX+2, inDelayMap("X", target->localWireDelay() + getCriticalPath() ) );
				oplist.push_back(resultAdder);
				
				inPortMap(resultAdder, "X", "R0");
				inPortMap(resultAdder, "Y", "R1");
				inPortMapCst(resultAdder, "Cin", (signedInputs? "'1'" : "'0'"));
				outPortMap( resultAdder, "R", "RAdder");
				vhdl << tab << instance(resultAdder, "ResultAdder") << endl;
				syncCycleFromSignal("RAdder");
				setCriticalPath( resultAdder->getOutputDelay("R"));
				
				vhdl << tab << "R <= RAdder;"<<endl;	
				outDelayMap["R"] = getCriticalPath();
			}
		} 


		
		if(useDSP) {
			
			//test if the multiplication fits into one DSP
			int wxDSP, wyDSP;
			target->getDSPWidths(wxDSP, wyDSP, signedInputs);
			bool testForward, testReverse, testFit;
			testForward     = (wX<=wxDSP)&&(wY<=wyDSP);
			testReverse = (wY<=wxDSP)&&(wX<=wyDSP);
			testFit = testForward || testReverse;
			

			if (testFit){
				if(target->worthUsingDSP(wX, wY))
					{
						manageCriticalPath(target->DSPMultiplierDelay());
						if (signedInputs)
							vhdl << tab << declare("rfull", wX + wY) << " <= XX * YY;"<<endl;
						else //sign extension is necessary for using use ieee.std_logic_signed.all; 
							// for correct inference of Xilinx DSP functions
							vhdl << tab << declare("rfull", wX + wY + 2) << " <= (\"0\" & XX) * (\"0\" & YY);"<<endl;
						
						vhdl << tab << "R <= rfull"<<range(wX + wY-1, 0)<<";"<<endl;	
						outDelayMap["R"] = getCriticalPath();
					}
				else {
					// For this target and this size, better do a logic-only implementation
					buildLogicOnly();
				}
			}
			else {
				buildTiling();
			}
		} 

		else {// This target has no DSP, going for a logic-only implementation
			buildLogicOnly();			
		}
	}
	
	/**************************************************************************/
	IntMultiplier::~IntMultiplier() {
	}









	// Stuff for the small multiplier table 
	IntMultiplier::SmallMultTable::SmallMultTable(Target* target, int dx, int dy) : 
		Table(target, dx+dy, dx+dy, 0, -1, true), // logic table
		dx(dx), dy(dy) {
		ostringstream name; 
		srcFileName="LogicIntMultiplier::SmallMultTable";
				name <<"SmallMultTable_" << dy << "x" << dx;
				setName(name.str());
				
				//				outDelayMap["Y"] = target->lutDelay();
	};
	
	mpz_class IntMultiplier::SmallMultTable::function(int yx){
		mpz_class r;
		int y = yx>>dx;
		int x = yx -(y<<dx);
		r = x*y;
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









	/**************************************************************************/
	void IntMultiplier::buildLogicOnly() {
		Target *target=getTarget();
		vhdl << tab << "-- code generated by IntMultiplier::buildLogicOnly()"<< endl;
		if(signedInputs) {
			// We will use the Baugh-Wooley tinkering:
			// split X as sx and Xr, Y as sy and Xr
			// then XY =  Xr*Yr (unsigned)      + bits and pieces
			// The purpose of the two following lines is to build the heap of bits for Xr*Yr
			// The remaining bits will be added to the heap later.
			wX--;
			wY--;
		}

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
			vhdl<<tab<<declare("Xp",sizeXPadded)<<" <= "<< zg(padX) << " & XX" << range(wX-1, 0) << ";"<<endl;
			vhdl<<tab<<declare("Yp",sizeYPadded)<<" <= "<< zg(padY) << " & YY" << range(wX-1, 0) << ";"<<endl;	
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++)
				vhdl<<tab<<declare(join("x",k),dx)<<" <= Xp"<<range((k+1)*dx-1,k*dx)<<";"<<endl;
			for (int k=0; k<chunksY ; k++)
				vhdl<<tab<<declare(join("y",k),dy)<<" <= Yp"<<range((k+1)*dy-1, k*dy)<<";"<<endl;
					

			int maxWeight = wOut+g;

			//initializing the bag of weighted bits (size is conservative, not all will be used)
			REPORT(DEBUG, "maxWeight=" << maxWeight);
			vector<unsigned int> currentHeight;
			//			for(int i=0; i<sizeXPadded+sizeYPadded; i++)
			for(int i=0; i<wX+wY; i++)
				currentHeight.push_back(0);
					
			// manageCriticalPath(target->localWireDelay(chunksX) + target->lutDelay());
					
			SmallMultTable *t = new SmallMultTable( target, dx, dy );
			oplist.push_back(t);

			for (int iy=0; iy<chunksY; iy++){
				vhdl << tab << "-- Partial product row number " << iy << endl;
				for (int ix=0; ix<chunksX; ix++) { 
					vhdl << tab << declare(XY(ix,iy), dx+dy) << " <= y" << iy << " & x" << ix << ";"<<endl;
					inPortMap(t, "X", XY(ix,iy));
					outPortMap(t, "Y", PP(ix,iy));
					vhdl << instance(t, PPTbl(ix,iy));
					vhdl << tab << "-- Adding these bits to the heap of bits" << endl << tab;
					// ignore the bits that are actually from the padding. 
					// In case of truncation, we still put the full lozange: truncation is handled when copying these bits to vectors
					// This makes everything simpler, in particular signed case handling
					// but the VHDL is longer than it could be
					int maxK=dx+dy; 
					if(ix == chunksX-1)
						maxK-=padX;
					if(iy == chunksY-1)
						maxK-=padY;
					for (int k=0; k<maxK; k++) {
						int weight=ix*dx+iy*dy+k;
						int height=currentHeight[weight];
						vhdl << declare(heap(weight,height)) << " <= " << PP(ix,iy) << of(k) << ";   ";
						currentHeight[weight] ++ ;
					}
					vhdl << endl;
				}
			}

				
			if(signedInputs) {
				int weight, height;
				// reminder: wX and wY have been decremented
				vhdl << tab << "-- Baugh-Wooley tinkering" << endl;
				vhdl << tab << declare("sX") << " <= XX" << of(wX) << ";" << endl;
				vhdl << tab << declare("sY") << " <= YY" << of(wY) << ";" << endl;
				vhdl << tab << declare("rX", wX) << " <= XX" << range(wX-1, 0) << ";" << endl;
				vhdl << tab << declare("rY", wY) << " <= YY" << range(wY-1, 0) << ";" << endl;
				vhdl << tab << declare("sXrYb", wY) << " <= not rY when sX='1' else " << zg(wY) << ";" << endl;
				vhdl << tab << declare("sYrXb", wX) << " <= not rX when sY='1' else " << zg(wX) << ";" << endl;
				vhdl << tab << "-- Adding these bits to the heap of bits" << endl << tab;
				for (int k=0; k<wX; k++) {
					weight=wY+k;
					height=currentHeight[weight];
					vhdl << declare(heap(weight,height)) << " <= sYrXb" << of(k) << ";  ";
					currentHeight[weight] ++ ;
				}
				vhdl << endl << tab;
				for (int k=0; k<wY; k++) {
					weight=wX+k;
					height=currentHeight[weight];
					vhdl << declare(heap(weight,height)) << " <= sXrYb" << of(k) << ";  ";
					currentHeight[weight] ++ ;
				}
				vhdl << endl;
				// adding sX and sY at positions wX and wY
				weight=wX;
				height=currentHeight[weight];
				vhdl << declare(heap(weight,height)) << " <= sX;  ";
				currentHeight[weight] ++ ;
				weight=wY;
				height=currentHeight[weight];
				vhdl << declare(heap(weight,height)) << " <= sY;  ";
				currentHeight[weight] ++ ;
				// adding sXb and sYb at position wX+wY
				weight=wX+wY;
				height=currentHeight[weight];
				vhdl << tab << declare(heap(weight,height)) << " <= not sX;  ";
				vhdl << declare(heap(weight,height+1)) << " <= not sY;  " << endl;
				currentHeight[weight] +=2 ;
				// adding 1 at position wX+wY+1
				weight=wX+wY+1;
				height=currentHeight[weight];
				vhdl << tab << declare(heap(weight,height)) << " <= '1';" << endl;
				currentHeight[weight] ++ ;		 
				// adding sXsY at position wX+wY
				weight=wX+wY;
				height=currentHeight[weight];
				vhdl << tab << declare(heap(weight,height)) << " <= sX and sY;" << endl;
				currentHeight[weight] ++ ;		 

				// restore wX and wY
				wX++;
				wY++;
			
			}

			if(g>0) {
				vhdl << tab << "-- The round bit " << endl;
				int weight=wX+wY-wOut-1;
				int height=currentHeight[weight];
				vhdl << declare(heap(weight,height)) << " <= '1';   ";
				currentHeight[weight] ++ ;
			}
		
			vhdl << tab << "-- Now building the heap vectors" << endl;
			
			vector<unsigned> finalHeapSizes;
			for(int weight=0;	weight<maxWeight; weight++) {
				int weightShift = wX+wY -(wOut+g);
				int bitWeight = weight + weightShift; // Here is where we manage truncation
				finalHeapSizes.push_back(currentHeight[bitWeight]);
				vhdl << tab << declare(join("HeapVector", weight), currentHeight[bitWeight]) << " <= " ;
				for (unsigned h=0; h<currentHeight[bitWeight]; h++) 
					vhdl << heap(bitWeight,h) <<  " & ";
				vhdl <<  "\"\";" << endl;
				// first copy the height vector, 
			}

			// Now instantiate a compressor tree

			NewCompressorTree *ct  = new NewCompressorTree(target, finalHeapSizes);
			oplist.push_back(ct);

			for(int weight=0;	weight<maxWeight; weight++) {
				inPortMap(ct, join("X",weight),  join("HeapVector", weight));
			}
			outPortMap(ct, "R", "RR");
			vhdl << instance(ct, "Compressor");
			//		int ctOutSize=ct->wOut;
	 
			// ctOutSize = wOut+1 but we know that the result always fits on wOut bits
			vhdl << tab << "R <= " << "RR" << range(wOut+g-1, g) << ";" << endl; 
		}
		

	}






	/**************************************************************************/
	void IntMultiplier::buildTiling() {
		vhdl << "-- buildTiling TODO"<< endl;
	}





	
	/******************************************************************************/
	void IntMultiplier::emulate ( TestCase* tc ) {
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svR;
		
		if (! signedInputs){
			svR = svX * svY;
		}

		else{ // Manage signed digits
			mpz_class big1 = (mpz_class(1) << (wX));
			mpz_class big1P = (mpz_class(1) << (wX-1));
			mpz_class big2 = (mpz_class(1) << (wY));
			mpz_class big2P = (mpz_class(1) << (wY-1));

			if ( svX >= big1P)
				svX = svX-big1;

			if ( svY >= big2P)
				svY = svY -big2;
			
			svR = svX * svY;
			if ( svR < 0){
				svR += (mpz_class(1) << (wX+wY)); 
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
		x = mpz_class(1); 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
		
		// -1 * -1
		x = (mpz_class(1) << (wX-1)) -1; 
		x = (mpz_class(1) << (wY-1)) -1; 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}


}


