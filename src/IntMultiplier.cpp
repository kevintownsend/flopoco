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
Who calls whom
the constructor calls buildLogicOnly or buildTiling
(maybe these should be unified some day)
They call buildHeapTiling or buildHeapLogicOnly

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

*/






/* For two's complement arithmetic on n bits, the representable interval is [ -2^{n-1}, 2^{n-1}-1 ]
   so the product lives in the interval [-2^{n-1}*2^{n-1}-1,  2^n]
   The value 2^n can only be obtained as the product of the two minimal negative input values
   (the weird ones, which have no symmetry)
   Example on 3 bits: input interval [-4, 3], output interval [-12, 16] and 16 can only be obtained by -4*-4.
   So the output would be representable on 2n-1 bits in two's complement, if it werent for this weird*weird case.

   So for signed multipliers, we just keep the 2n bits, including one bit used for only for this weird*weird case.
   Current situation is: this must be managed from outside:
   An application that knows that it will not use the full range (e.g. signal processing, poly evaluation) can ignore the MSB bit, 
   but we produce it.



   A big TODO ?
  
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

#define vhdl parentOp->vhdl
#define declare parentOp->declare
#define inPortMap parentOp->inPortMap
#define outPortMap parentOp->outPortMap
#define instance parentOp->instance
#define manageCriticalPath parentOp->manageCriticalPath
#define oplist parentOp->getOpListR()

	int IntMultiplier::neededGuardBits(int wX, int wY, int wOut){
		int g;
		if(wX+wY==wOut)
			g=0;
		else {
			unsigned i=0;
			mpz_class ulperror=1;
			while(wX+wY - wOut  > intlog2(ulperror)) {
				// REPORT(DEBUG,"i = "<<i<<"  ulperror "<<ulperror<<"  log:"<< intlog2(ulperror) << "  wOut= "<<wOut<< "  wFull= "<< wX+wY);
				i++;
				ulperror += (i+1)*(mpz_class(1)<<i);
			}
			g=wX+wY-i-wOut;
			// REPORT(DEBUG, "ulp truncation error=" << ulperror << "    g=" << g);
		}
		return g;
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

		g = neededGuardBits(wXdecl, wYdecl, wOut);

		maxWeight = wOut+g;
		weightShift = wFull - maxWeight;  



		// Halve number of cases by making sure wY<=wX:
		// interchange x and y in case wY>wX

		if(wYdecl> wXdecl){
			wX=wYdecl;
			wY=wXdecl;
			

			vhdl << tab << declare(join("XX",multiplierUid), wX, true) << " <= " << yname << ";" << endl; 
			vhdl << tab << declare(join("YY",multiplierUid), wY, true) << " <= " << xname << ";" << endl; 

		}
		else{
			wX=wXdecl;
			wY=wYdecl;

			vhdl << tab << declare(join("XX",multiplierUid), wX, true) << " <= " << xname << ";" << endl; 
			vhdl << tab << declare(join("YY",multiplierUid), wY, true) << " <= " << yname << ";" << endl; 

		}

			

	}




	
	// The virtual constructor
	IntMultiplier::IntMultiplier (Operator* parentOp_, BitHeap* bitHeap_, Signal* x_, Signal* y_, int wX_, int wY_, int wOut_, int lsbWeight_, bool negate_, bool signedIO_, float ratio_):
		Operator ( parentOp_->getTarget()), 
		wXdecl(wX_), wYdecl(wY_), wOut(wOut_), ratio(ratio_),  maxError(0.0), parentOp(parentOp_), bitHeap(bitHeap_), x(x_), y(y_), negate(negate_), signedIO(signedIO_) {

		isOperator=false;

		multiplierUid=parentOp->getNewUId();
		srcFileName="IntMultiplier";
		useDSP = (ratio>=0) && getTarget()->hasHardMultipliers();
		
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

		plotter = new Plotter(bitHeap);

		bitHeap->setPlotter(plotter);
		//plotter->setBitHeap(bitHeap);
		
		initialize();
		fillBitHeap();


		// leave the compression to the parent op
	}







	// The constructor for a stand-alone operator
	IntMultiplier::IntMultiplier (Target* target, int wX_, int wY_, int wOut_, bool signedIO_, float ratio_, map<string, double> inputDelays_):
		Operator ( target, inputDelays_ ), wXdecl(wX_), wYdecl(wY_), wOut(wOut_), ratio(ratio_), maxError(0.0), signedIO(signedIO_) {
		
		isOperator=true;
		srcFileName="IntMultiplier";
		setCopyrightString ( "Florent de Dinechin, Kinga Illyes, Bogdan Popa, Bogdan Pasca, 2012" );
		
		// useDSP or not? 
		//useDSP = (ratio>0) && target->hasHardMultipliers();
		useDSP = (ratio>=0)&&target->hasHardMultipliers();


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
		multiplierUid=parentOp->getNewUId();
		xname=join("X",multiplierUid);
		yname=join("Y",multiplierUid);

		initialize();

		// Set up the IO signals
		addInput ( xname  , wXdecl, true );
		addInput ( yname  , wYdecl, true );

		// TODO FIXME This 1 should be 2. It breaks TestBench but not TestBenchFile. Fix TestBench first! (check in addExpectedOutput or something)
		addOutput ( join("R",multiplierUid)  , wOut, 1 , true );

		// Set up the VHDL library style
		if(signedIO)
			useStdLogicSigned();
		else
			useStdLogicUnsigned();

		// The bit heap
		bitHeap = new BitHeap(this, wOut+g);

		plotter = new Plotter(bitHeap);

		bitHeap->setPlotter(plotter);
		//plotter->setBitHeap(bitHeap);



	

		// initialize the critical path
		setCriticalPath(getMaxInputDelays ( inputDelays_ ));

		fillBitHeap();

		bitHeap -> generateCompressorVHDL();			
		vhdl << tab << join("R",multiplierUid) << " <= " << bitHeap-> getSumName() << range(wOut+g-1, g) << ";" << endl;
	}






	/// TODO FIXME all the special case in case of a shared bit heap
	void  IntMultiplier::fillBitHeap(){

		///////////////////////////////////////
		//  architectures for corner cases   //
		// TODO manage non-standalone case here //
		///////////////////////////////////////

		// The really small ones fit in two LUTs and that's as small as it gets  
		if(wX+wY <= target()->lutInputs()+2) {
			REPORT(DETAILED,"1");
			vhdl << tab << "-- Ne pouvant me fier à mon raisonnement, j'ai appris par coeur le résultat de toutes les multiplications possibles" << endl;

			SmallMultTable *t = new SmallMultTable( target(), wX, wY, wOut, negate, signedIO, signedIO);
			//useSoftRAM(t);
			oplist.push_back(t);

			vhdl << tab << declare(join("XY",multiplierUid), wX+wY) << " <= "<<join("YY",multiplierUid)<<" & "<<join("XX",multiplierUid)<<";"<<endl;

			inPortMap(t, "X", join("XY",multiplierUid));
			outPortMap(t, "Y", join("RR",multiplierUid));
			plotter->addSmallMult(0,0,wX,wY);
			bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);
			vhdl << instance(t, "multTable");
			
			for(int k=wX+wY-1; k>=0; k--)
			{
				if(k>=wX+wY-wOut-g)
				{
					stringstream s;
					s<<join("RR",multiplierUid)<<"("<<k<<")";
					bitHeap->addBit(k,s.str());
				}	
			}		
			
			//vhdl << tab << join("R",multiplierUid)<<" <= "<<join("RR",multiplierUid)<<";" << endl;

			//setCriticalPath(t->getOutputDelay(join("Y",multiplierUid)));

			//outDelayMap[join("R",multiplierUid)] = getCriticalPath();
			return;
			
		}

		// Multiplication by 1-bit integer is simple
		if ((wY == 1)){
			REPORT(DETAILED,"2");
			if (signedIO){
				manageCriticalPath( target()->localWireDelay(wX) + target()->adderDelay(wX+1) );

				vhdl << tab << join("R",multiplierUid) <<" <= (" << zg(wX+1)  << " - ("<<join("XX",multiplierUid)<< of(wX-1) << " & "<<join("XX",multiplierUid)<<")) when "<<join("YY",multiplierUid)<<"(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	

			}
			else {
				manageCriticalPath( target()->localWireDelay(wX) + target()->lutDelay() );

				vhdl << tab << join("R",multiplierUid)<<" <= (\"0\" & "<<join("XX",multiplierUid)<<") when "<< join("YY",multiplierUid)<<"(0)='1' else "<<zg(wX+1,0)<<";"<<endl;	

			}
			outDelayMap[join("R",multiplierUid)] = getCriticalPath();
			return;
		}


		// Multiplication by 2-bit integer is one addition
		if ((wY == 2)) {
			// No need for the following, the mult should be absorbed in the addition
			// manageCriticalPath( target->localWireDelay() + target->lutDelay() );
			REPORT(DETAILED,"3");
			vhdl << tab << declare(join("R0",multiplierUid),wX+2) << " <= (";

			if (signedIO) 

				vhdl << join("XX",multiplierUid) << of(wX-1) << " & "<< join("XX",multiplierUid) << of(wX-1);  

			else  

				vhdl << "\"00\"";
			vhdl <<  " & "<<join("XX",multiplierUid)<<") when "<<join("YY",multiplierUid)<<"(0)='1' else "<<zg(wX+2,0)<<";"<<endl;	
			vhdl << tab << declare(join("R1i",multiplierUid),wX+2) << " <= ";

			if (signedIO) 

				vhdl << "("<<join("XX",multiplierUid)<< of(wX-1) << "  &  " <<join("XX",multiplierUid)<<" & \"0\")";

			else  

				vhdl << "(\"0\" & "<< join("XX",multiplierUid) <<" & \"0\")";
			vhdl << " when "<<join("YY",multiplierUid)<<"(1)='1' else " << zg(wX+2,0) << ";"<<endl;	
			vhdl << tab << declare(join("R1",multiplierUid),wX+2) << " <= ";

			if (signedIO) 

				vhdl << "not "<<join("R1i",multiplierUid)<<";" <<endl;

			else  

				vhdl << join("R1i",multiplierUid)<<";"<<endl;

			
			IntAdder *resultAdder = new IntAdder( target(), wX+2, inDelayMap("X", target()->localWireDelay() + getCriticalPath() ) );
			oplist.push_back(resultAdder);
			
			inPortMap(resultAdder, join("X",multiplierUid), join("R0",multiplierUid));
			inPortMap(resultAdder, join("Y",multiplierUid), join("R1",multiplierUid));
			inPortMapCst(resultAdder, join("Cin",multiplierUid), (signedIO? "'1'" : "'0'"));
			outPortMap( resultAdder, join("R",multiplierUid), join("RAdder",multiplierUid));

			vhdl << tab << instance(resultAdder, join("ResultAdder",multiplierUid)) << endl;

			syncCycleFromSignal(join("RAdder",multiplierUid));
			setCriticalPath( resultAdder->getOutputDelay(join("R",multiplierUid)));

			vhdl << tab << join("R",multiplierUid)<<"<= "<<  join("RAdder",multiplierUid) << range(wFull-1, wFull-wOut)<<";"<<endl;	

			outDelayMap[join("R",multiplierUid)] = getCriticalPath();
			return;
		} 


	
		
		// Now getting more and more generic
		if(useDSP) {
			//test if the multiplication fits into one DSP
			//int wxDSP, wyDSP;
			target()->getDSPWidths(wxDSP, wyDSP, signedIO);
			bool testForward, testReverse, testFit;
			testForward     = (wX<=wxDSP)&&(wY<=wyDSP);
			testReverse = (wY<=wxDSP)&&(wX<=wyDSP);
			testFit = testForward || testReverse;
		
			
			REPORT(DETAILED,"useDSP");
			if (testFit){
			REPORT(DETAILED,"testfit");
			
				if( target()->worthUsingDSP(wX, wY))
					{	REPORT(DEBUG,"worthUsingDSP");
						manageCriticalPath(target()->DSPMultiplierDelay());
						if (signedIO)
						{
							vhdl << tab << declare(join("rfull",multiplierUid), wFull+1) << " <= "<<join("XX",multiplierUid)<<"  *  "<< join("YY",multiplierUid)<<"; -- that's one bit more than needed"<<endl; 
							vhdl << tab << join("R",multiplierUid)<<" <= "<< join("rfull",multiplierUid) <<range(wFull-1, wFull-wOut)<<";"<<endl;	
							outDelayMap[join("R",multiplierUid)] = getCriticalPath();
						}
						else //sign extension is necessary for using use ieee.std_logic_signed.all; 
						{	// for correct inference of Xilinx DSP functions

							//vhdl << tab << declare(join("rfull",multiplierUid), wX + wY + 2) << " <= (\"0\" & "<<join("XX",multiplierUid)<<
							//") * (\"0\" &"<<join("YY",multiplierUid)<<");"<<endl;


							int topx=wX-wxDSP;
							int topy=wY-wyDSP;

							MultiplierBlock* m = new MultiplierBlock(wxDSP,wyDSP,topx, topy,
									join("XX",multiplierUid),join("YY",multiplierUid),weightShift);
							m->setNext(NULL);		
							m->setPrevious(NULL);			
							localSplitVector.push_back(m);
							bitHeap->addDSP(m);
							bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);
						}	
						
					}
					
					
				else {
					// For this target and this size, better do a logic-only implementation
					REPORT(DETAILED,"before buildlogic only");
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
	void IntMultiplier::buildLogicOnly() {
		buildHeapLogicOnly(0,0,wX,wY);
		//adding the round bit
		if(g>0) {
			int weight=wFull-wOut-1- weightShift;
			bitHeap->addConstantOneBit(weight);
		}
	}



	/**************************************************************************/
	void IntMultiplier::buildTiling() {
		
		
			buildHeapTiling();
			//adding the round bit
			//bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);
			if(g>0) {
				int weight=wFull-wOut-1- weightShift;
				bitHeap->addConstantOneBit(weight);
			}
			
		
	}





		/**************************************************************************/
		void IntMultiplier::buildHeapLogicOnly(int topX, int topY, int botX, int botY,int uid) {
			Target *target=getTarget();
			if(uid==-1)
				uid++;
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
			if (chunksX + chunksY > 2) { //we do more than 1 subproduct // FIXME where is the else?
				
				// Padding X to the left
				vhdl << tab << declare(join("Xp",multiplierUid,"_",uid),sizeXPadded)<<" <= ";
				if(padX>0) {
					if(signedIO)	{ // sign extension
						ostringstream signbit;
						signbit << join("XX",multiplierUid) << of(botX-1);
						vhdl  << rangeAssign(sizeXPadded-1, wX, signbit.str()) << " & ";
					}
					else {
						vhdl << zg(padX) << " & ";
					}
				}
				vhdl << join("XX",multiplierUid) << range(botX-1,topX) << ";"<<endl;

				// Padding Y to the left
				vhdl << tab << declare(join("Yp",multiplierUid,"_",uid),sizeYPadded)<<" <= ";
				if(padY>0) {
					if(signedIO)	{ // sign extension		
						ostringstream signbit;
						signbit << join("YY",multiplierUid) << of(botY-1);
						vhdl  << rangeAssign(sizeYPadded-1, wY, signbit.str()) << " & ";
					}
					else {
						vhdl << zg(padY) << " & ";
					}
				}
				vhdl << join("YY",multiplierUid) << range(botY-1,topY) << ";"<<endl;

				//SPLITTINGS
				for (int k=0; k<chunksX ; k++)
					vhdl<<tab<<declare(join("x",multiplierUid,"_",uid,"_",k),dx)<<" <= "<<join("Xp",multiplierUid,"_",uid)<<range((k+1)*dx-1,k*dx)<<";"<<endl;
					
				for (int k=0; k<chunksY ; k++)
					vhdl<<tab<<declare(join("y",multiplierUid,"_",uid,"_",k),dy)<<" <= "<<join("Yp",multiplierUid,"_",uid)<<range((k+1)*dy-1, k*dy)<<";"<<endl;
					
				REPORT(DEBUG, "maxWeight=" << maxWeight <<  "    weightShift=" << weightShift);
				SmallMultTable *tUU, *tSU, *tUS, *tSS;


				tUU = new SmallMultTable( target, dx, dy, dx+dy, negate, false, false);
				//useSoftRAM(tUU);
				oplist.push_back(tUU);

				if(signedIO) { // need for 4 different tables
					
					tSU = new SmallMultTable( target, dx, dy, dx+dy, negate, true, false );
					//useSoftRAM(tSU);
					oplist.push_back(tSU);
					
					tUS = new SmallMultTable( target, dx, dy, dx+dy, negate, false, true );
					//useSoftRAM(tUS);
					oplist.push_back(tUS);
					
					tSS = new SmallMultTable( target, dx, dy, dx+dy, negate, true, true );
					//useSoftRAM(tSS);
					oplist.push_back(tSS);
				}

				setCycle(0); // TODO FIXME for the virtual multiplier case where inputs can arrive later
				setCriticalPath(initialCP);
				// SmallMultTable is built to cost this:
				manageCriticalPath( getTarget()->localWireDelay(chunksX) + getTarget()->lutDelay() ) ;  
				for (int iy=0; iy<chunksY; iy++){

					vhdl << tab << "-- Partial product row number " << iy << endl;

					for (int ix=0; ix<chunksX; ix++) { 
					
					
						SmallMultTable *t;
						if (!signedIO) {
							t=tUU;
						}
						else { // 4  cases 
							if((ix==chunksX-1) && (iy==chunksY-1))
								t=tSS;
							else if (ix==chunksX-1)
								t=tSU;
							else if (iy==chunksY-1)
								t=tUS;
							else
								t=tUU; 
						}
						
				
						if(dx*(ix+1)+dy*(iy+1)+topX+topY>wFull-wOut-g)
						{
						
						plotter->addSmallMult(dx*(ix)+topX, dy*(iy)+topY,dx,dy);

						vhdl << tab << declare(join (XY(ix,iy,uid),multiplierUid), dx+dy) 
						               << " <= y"<< multiplierUid<<"_"<<uid <<"_"<< iy << " & x"<<multiplierUid<<"_" << uid <<"_"<< ix << ";"<<endl;

						inPortMap(t, "X", join(XY(ix,iy,uid),multiplierUid));
						outPortMap(t, "Y", join(PP(ix,iy,uid),multiplierUid));
						vhdl << instance(t, join(PPTbl(ix,iy,uid),multiplierUid));
						
						vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;
						
						// Two's complement management
						// There are really 2 cases:
						// If the result is known positive, ie if tUU and !negate, nothing to do
						// If the result is in two's complement  
						//    sign extend by adding ones on weights >= the MSB of the table, so its sign is propagated.
						//    Also need to complement the sign bit
						// Note that even when negate and tUU, the result is known negative, but it may be zero, so its actual sign is not known statically
						

						bool resultSigned;  
						if(negate || (t==tSS) || (t==tUS) || (t==tSU)) 
							resultSigned = true ;

						int maxK=dx+dy; 
 						if(ix == chunksX-1)
							maxK-=padX;
						if(iy == chunksY-1)
							maxK-=padY;
						REPORT(DEBUG,  "ix=" << ix << " iy=" << iy << "  maxK=" << maxK  << "  negate=" << negate <<  "  resultSigned="  << resultSigned );
						for (int k=0; k<maxK; k++) {
							ostringstream s;
							s << join(PP(ix,iy,uid),multiplierUid) << of(k); // right hand side
							int weight = ix*dx+iy*dy+k - weightShift+topX+topY;
							if(weight>=0) {// otherwise these bits deserve to be truncated
								if(resultSigned && (k==maxK-1)) { 
									ostringstream nots;
									nots << "not " << s.str(); 
									bitHeap->addBit(weight, nots.str());
									REPORT(DEBUG, "adding bit " << nots.str() << " at weight " << weight); 
									REPORT(DEBUG,  "  adding constant ones from weight "<< weight << " to "<< bitHeap->getMaxWeight());
									for (unsigned w=weight; w<bitHeap->getMaxWeight(); w++)
										bitHeap->addConstantOneBit(w);
								}
								else { // just add the bit
									REPORT(DEBUG, "adding bit " << s.str() << " at weight " << weight); 
									bitHeap->addBit(weight, s.str());
								}
							}
						}

						vhdl << endl;
						
						}
					}
				}				

		
		
			}
	 
		}
	


	void IntMultiplier::splitting(int horDSP, int verDSP, int wxDSP, int wyDSP,int restX, int restY)
		{
	
			int i=0;
			int j=0;
			int x=wX;
			while(i<verDSP)
			{
				j=0;
				int ok=0;
				
				while((j<horDSP)&&(ok==0))
				{	REPORT(DETAILED,"j= " << j);
					//if the top right corner is under the truncation line, then will use a DSP; otherwise, will be decided according to the ratio
					//by calling the compute function
					if((wX-(j+1)*wxDSP)+(wY-(i+1)*wyDSP)>=wFull-wOut-g)
					{
						//a DSP can be used
						MultiplierBlock* m = new MultiplierBlock(wxDSP,wyDSP,wX-(j+1)*wxDSP, wY-((i+1)*wyDSP),
									join("XX",multiplierUid),join("YY",multiplierUid),weightShift);
						m->setNext(NULL);		
						m->setPrevious(NULL);			
						REPORT(DETAILED,"new DSP sure fit");
						localSplitVector.push_back(m);
						bitHeap->addDSP(m);
						
					}
					
					else
					{
						//cannot use DSPs without checking the ratio
						ok=1;
						j--;
					}
					
					j++;
					
									
				}
				
				//here we have a horizontal block, length=(from the last used DSP's top-right coordinates to the truncation line)
				
				//determination of the x coordinate
				int y= wY-(i)*wyDSP;
				x=wX;
				while((x+y>wFull-wOut-g) && (x>0))
					x--;
				
				
				//call the function only if at least 1 bit remaining
				if((wX-j*wxDSP>0))
					compute(x,wY-(i+1)*wyDSP, wX-j*wxDSP, wY-(i)*wyDSP,wxDSP,wyDSP); 
					
				i++;		
			}
			
			
				//if there are some remaining bits on the Y 
				if(restY>0)
				{
					//determination of x coordinate (top right)
					int y=restY;
					int x=wX;
					while((x+y>wFull-wOut-g)&&(x>0))
						x--;
					compute(x,0,wX,restY,wxDSP,wyDSP);
				}
		}
		
		
		
		
		void IntMultiplier::compute(int topX, int topY, int botX, int botY,int wxDSP,int wyDSP)
		{
			int height=botY-topY;
			int width=botX-topX;
			int dspArea=wxDSP*wyDSP;
			bool was=false; // tells if the while loop was executed or not
			int botx=botX;
			int topx=topX;
			int topy=topY;
			int dsp=0;//number of used DSPs
			
			//if the width is larger then a dsp width, than we have to compute the good coordinates for the dsp
			if (width>wxDSP)
				topx=botx-wxDSP;
		
			REPORT(DETAILED,"compute called, width= "<<width);
		
			while (width>wxDSP)
			{	
				REPORT(DETAILED,"width greater than DSPwidth");
				//we need to split the block
				was=true;
				float blockArea=wxDSP*height; //the area of the block that will be analyzed
				
				//computing the area of the triangle that will be lost 
				int tx=topx;
				int ty=topy;
				while(tx+ty<wFull-wOut-g)
					tx++;
			    int triangleArea=((tx-topx)*(tx-topx))/2;
			    
			    //the area which will be used
				blockArea=blockArea-triangleArea;
				
				
				//checking the use according to the ratio
				if((blockArea>=(1.0-ratio)*dspArea))
				{  
				
					if(height<wyDSP)
						topy=topY-(wyDSP-height);
					MultiplierBlock* m = new MultiplierBlock(wxDSP,wyDSP,topx,topy,
								join("XX",multiplierUid),join("YY",multiplierUid),weightShift);
					m->setNext(NULL);		
					m->setPrevious(NULL);			
					localSplitVector.push_back(m);
					bitHeap->addDSP(m);
					
				}
				else
				{
					
					if((topx<botX-dsp*wxDSP-1))
						buildHeapLogicOnly(topx, topY,(botX-dsp*wxDSP),botY,parentOp->getNewUId());	
				}
				
				dsp++;
				botx=topx-1;
				topx=topx-wxDSP;
				width=width-wxDSP;
				if(width<=wxDSP)
					topx=topX;
			}
			
			REPORT(DETAILED,"width smaller than DSPwidth");
			
			//now the remaining block is smaller (for sure!) than the DSP width
			//we do the same area / ratio checking again	
			float blockArea=width*height;
			int triangleArea=((width)*(width))/2;
			blockArea=blockArea-triangleArea;
	
			if(blockArea>=(1.0-ratio)*dspArea)
			{
			
				if(was)
					topx=topx-(wxDSP-(botx-topx+1));
				else 
					topx=topX-(wxDSP-(botX-topX+1))-1; 
						
				if(height<wyDSP)
					topy=topY-(wyDSP-height);
						
				MultiplierBlock* m = new MultiplierBlock(wxDSP,wyDSP,topx,topy,
							join("XX",multiplierUid),join("YY",multiplierUid),weightShift);
				m->setNext(NULL);		
				m->setPrevious(NULL);			
				localSplitVector.push_back(m);
				bitHeap->addDSP(m);
			}
			else
			{
				if((topx<botX-dsp*wxDSP))
					buildHeapLogicOnly(topx,topY,botX-dsp*wxDSP,botY,parentOp->getNewUId());
			}
		
		}
	


		/** builds the tiles and the logic too*/
		/**************************************************************************/
		void IntMultiplier::buildHeapTiling() {
		
			//the DSPs should be arranged horizontally or vertically?
		
			//number of horizontal/vertical DSPs used if the tiling is horizontal
			int horDSP1=wX/wxDSP;
			int verDSP1=wY/wyDSP;

			//number of horizontal/vertical DSPs used if the tiling is vertical
			int horDSP2=wX/wyDSP;
			int verDSP2=wY/wxDSP;

			//the size of the zone filled by DSPs
			int hor=horDSP1*verDSP1;
			int ver=horDSP2*verDSP2;
 
			int horDSP;
			int verDSP;
			int restX; //the number of lsbs of the first input which remains after filling with DSPs
			int restY; //the number of lsbs of the second input which remains after filling with DSPs

			if (hor>=ver)
				{
					horDSP=horDSP1;
					verDSP=verDSP1;
					restX=wX-horDSP*wxDSP;
					restY=wY-verDSP*wyDSP;
					//splitting horizontal
					splitting(horDSP,verDSP,wxDSP,wyDSP,restX,restY);
				}
			else
				{
					horDSP=horDSP2;
					verDSP=verDSP2;
					restX=wX-horDSP*wyDSP;
					restY=wY-verDSP*wxDSP;
					//splitting vertical
					splitting(horDSP,verDSP,wyDSP,wxDSP,restX,restY);
				}

	
		bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);

		}
		
		
	

	
		IntMultiplier::~IntMultiplier() {
		}


		//signal name construction

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




		IntMultiplier::SmallMultTable::SmallMultTable(Target* target, int dx, int dy, int wO, bool negate, bool  signedX, bool  signedY ) : 
			Table(target, dx+dy, wO, 0, -1, true), // logic table
			dx(dx), dy(dy), negate(negate), signedX(signedX), signedY(signedY) {
			ostringstream name; 
			srcFileName="LogicIntMultiplier::SmallMultTable";
			name <<"SmallMultTable"<< (negate?"M":"P") << dy << "x" << dx << "r" << wO << (signedX?"Xs":"Xu") << (signedY?"Ys":"Yu")  << getuid();
			setName(name.str());				
		};
	

		mpz_class IntMultiplier::SmallMultTable::function(int yx){
			mpz_class r;
			int y = yx>>dx;
			int x = yx -(y<<dx);
			int wF=dx+dy;

			if(signedX){
				if ( x >= (1 << (dx-1)))
					x -= (1 << dx);
			}
			if(signedY){
				if ( y >= (1 << (dy-1)))
					y -= (1 << dy);
			}
			// if(negate && signedX && signedY) cerr << "  y=" << y << "  x=" << x;
			r = x * y;
			// if(negate && signedX && signedY) cerr << "  r=" << r;
			if(negate)
				r=-r;
			// if(negate && signedX && signedY) cerr << "  -r=" << r;
			if ( r < 0)
				r += mpz_class(1) << wF; 
			// if(negate && signedX && signedY) cerr << "  r2C=" << r;

			if(wOut<wF){ // wOut is that of Table
				// round to nearest, but not to nearest even
				int tr=wF-wOut; // number of truncated bits 
				// adding the round bit at half-ulp position
				r += (mpz_class(1) << (tr-1));
				r = r >> tr;
			}
			
			// if(negate && signedX && signedY) cerr << "  rfinal=" << r << endl;

			return r;
		
		};



		void IntMultiplier::emulate ( TestCase* tc ) {
			mpz_class svX = tc->getInputValue(join("X",multiplierUid));
			mpz_class svY = tc->getInputValue(join("Y",multiplierUid));
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
				}
			if(negate)
				svR = -svR;

			// manage two's complement at output
			if ( svR < 0){
				svR += (mpz_class(1) << wFull); 
			}
			if(wTruncated==0) 
				tc->addExpectedOutput(join("R",multiplierUid), svR);
			else {
				// there is truncation, so this mult should be faithful
				svR = svR >> wTruncated;
				tc->addExpectedOutput(join("R",multiplierUid), svR);
				svR++;
				svR &= (mpz_class(1) << (wOut)) -1;
				tc->addExpectedOutput(join("R",multiplierUid), svR);
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
			tc->addInput(join("X",multiplierUid), x);
			tc->addInput(join("Y",multiplierUid), y);
			emulate(tc);
			tcl->add(tc);
		
			// -1 * -1
			x = (mpz_class(1) << wXdecl) -1; 
			y = (mpz_class(1) << wYdecl) -1; 
			tc = new TestCase(this); 
			tc->addInput(join("X",multiplierUid), x);
			tc->addInput(join("Y",multiplierUid), y);
			emulate(tc);
			tcl->add(tc);

			// The product of the two max negative values overflows the signed multiplier
			x = mpz_class(1) << (wXdecl -1); 
			y = mpz_class(1) << (wYdecl -1); 
			tc = new TestCase(this); 
			tc->addInput(join("X",multiplierUid), x);
			tc->addInput(join("Y",multiplierUid), y);
			emulate(tc);
			tcl->add(tc);
		}


}




