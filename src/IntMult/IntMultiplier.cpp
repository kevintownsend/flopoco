/*
  An integer multiplier mess for FloPoCo

  TODO in the virtual multiplier case, manage properly the case when the initial instant is not  0

  Authors:  Bogdan Pasca, and then F de Dinechin, Matei Istoan, Kinga Illyes and Bogdan Popa spent 12 cumulated months getting rid of the last bits of his code

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2012.
  All rights reserved.
*/

// To enable SVG plotting, #define BITHEAP_GENERATE_SVG in BitHeap.hpp


// TODO 
// Again needs a complete overhaul I'm afraid
// - Improve the case when a multiplier fits one DSP (very common)
//     In this case even for truncated mults we should have simple truncation, g=0 (also in neededGuardBits),
//     and a special management of negate and signedIO


/* 
   Who calls whom
   the constructor calls buildLogicOnly or buildTiling
   (maybe these should be unified some day)
   They call buildXilinxTiling or buildAlteraTiling or buildHeapLogicOnly

*/



/* VHDL variable names:
   X, Y: inputs
   XX,YY: after swap

*/






/* For two's complement arithmetic on n bits, the representable interval is [ -2^{n-1}, 2^{n-1}-1 ]
   so the product lives in the interval [-2^{n-1}*2^{n-1}-1,  2^n]
   The value 2^n can only be obtained as the product of the two minimal negative input values
   (the weird ones, which have no opposite)
   Example on 3 bits: input interval [-4, 3], output interval [-12, 16] and 16 can only be obtained by -4*-4.
   So the output would be representable on 2n-1 bits in two's complement, if it werent for this weird*weird case.

   So even for signed multipliers, we just keep the 2n bits, including one bit used for only for this weird*weird case.
   Current situation is: if you don't like this you manage it from outside:
   An application that knows that it will not use the full range (e.g. signal processing, poly evaluation) can ignore the MSB bit, 
   but we produce it.

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
#include "IntAddSubCmp/IntAdder.hpp"
#include "Targets/StratixII.hpp"
#include "Targets/StratixIII.hpp"
#include "Targets/StratixIV.hpp"
#include "Plotter.hpp"

using namespace std;

namespace flopoco {

#define vhdl parentOp->vhdl
#define declare parentOp->declare
#define inPortMap parentOp->inPortMap
#define outPortMap parentOp->outPortMap
#define instance parentOp->instance
#define useSoftRAM parentOp->useSoftRAM
#define manageCriticalPath parentOp->manageCriticalPath
#define getCriticalPath parentOp->getCriticalPath
#define setCycle parentOp->setCycle
#define oplist parentOp->getOpListR()


	bool IntMultiplier::tabulatedMultiplierP(Target* target, int wX, int wY){
		return (wX+wY <=  target->lutInputs()+2); 
	}

	int IntMultiplier::neededGuardBits(Target* target, int wX, int wY, int wOut)
	{
		int g;
		if((  wOut >= wX+wY) // no rounding 
			 || tabulatedMultiplierP(target, wX, wY) ) { // multiplier will be tabulated 
			g=0;
		}
		else {
			unsigned i=0;
			mpz_class ulperror=1;
			while(wX+wY - wOut  > intlog2(ulperror)) {
				// REPORT(DEBUG,"i = "<<i<<"  ulperror "<<ulperror<<"  log:"<< intlog2(ulperror) << "  wOut= "<<wOut<< "  wFullP= "<< wX+wY);
				i++;
				ulperror += (i+1)*(mpz_class(1)<<i);
			}
			g=wX+wY-i-wOut;
			// REPORT(DEBUG, "ulp truncation error=" << ulperror << "    g=" << g);
		}
		return g;
	}


	void IntMultiplier::initialize() 
	{
		if(wXdecl<0 || wYdecl<0) {
			THROWERROR("negative input size");
		}
		wFullP = wXdecl+wYdecl;

		// Halve number of cases by making sure wY<=wX:
		// interchange x and y in case wY>wX
		// After which we negate y (the smaller) by 1/ complementing it and 2/  adding it back to the bit heap

		string newxname, newyname;
		if(wYdecl> wXdecl) {
			REPORT(DEBUG, "initialize(): Swapping X and Y")
			wX=wYdecl;	 
			wY=wXdecl;	 
			newxname=yname;
			newyname=xname;
		}
		else {
			wX=wXdecl;	 
			wY=wYdecl;
			newxname=xname;
			newyname=yname;
		}

		// The larger of the two 
		vhdl << tab << declare(addUID("XX"), wX, true) << " <= " << newxname << " ;" << endl;	 
		// possibly negate the smaller
		if(!negate)
			vhdl << tab << declare(addUID("YY"), wY, true) << " <= " << newyname << " ;" << endl;	 
		else
		{
			vhdl << tab << "-- we compute -xy as x(not(y)+1)" << endl;
			vhdl << tab << declare(addUID("YY"), wY, true) << " <= not " << newyname << " ;" << endl;	 
		}
	}





	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// The virtual constructor 
	IntMultiplier::IntMultiplier (Operator* parentOp_, BitHeap* bitHeap_, 
																Signal* x, Signal* y,
																int lsbWeightInBitHeap_, bool negate_, bool signedIO_, float DSPThreshold_):
		Operator (parentOp_->getTarget()), 
		DSPThreshold(DSPThreshold_),  
		parentOp(parentOp_), 
		bitHeap(bitHeap_), 
		lsbWeightInBitHeap(lsbWeightInBitHeap_),
		negate(negate_), 
		signedIO(signedIO_) 
	{

		multiplierUid=parentOp->getNewUId();
		srcFileName="IntMultiplier";
		useDSP = (DSPThreshold>=0) &&  parentOp->getTarget()->hasHardMultipliers();
		// TODO remove the following to switch DSP back on
		useDSP=false;

		xname = x->getName();
		yname = y->getName();
		wXdecl = x->width();
		wYdecl = y->width();

		ostringstream name;
		name <<"VirtualIntMultiplier";
		if(useDSP) 
			name << "UsingDSP_";
		else
			name << "LogicOnly_";
		name << wXdecl << "_" << wYdecl <<"_" << (lsbWeightInBitHeap<0?"m":"") << abs(lsbWeightInBitHeap) << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
		setName ( name.str() );
		
		REPORT(DEBUG, "Building " << name.str() );

		initialize();
		fillBitHeap();
		// leave the compression to the parent op
	}







	// The constructor for a stand-alone operator
	IntMultiplier::IntMultiplier (Target* target_, int wX_, int wY_, int wOut_, bool signedIO_, float DSPThreshold_, map<string, double> inputDelays_, bool enableSuperTiles_):
		Operator ( target_, inputDelays_ ), 
		wXdecl(wX_), wYdecl(wY_), wOut(wOut_),
		DSPThreshold(DSPThreshold_), negate(false), signedIO(signedIO_),enableSuperTiles(enableSuperTiles_), target(target_)
	{
		srcFileName="IntMultiplier";
		setCopyrightString ( "Florent de Dinechin, Kinga Illyes, Bogdan Popa, Bogdan Pasca, 2012" );

		// useDSP or not? 
		useDSP = (DSPThreshold>=0)&&target->hasHardMultipliers();
		{
			ostringstream name;
			name <<"IntMultiplier";
			if(useDSP) 
				name << "_UsingDSP_";
			else
				name << "_LogicOnly_";
			name << wXdecl << "_" << wYdecl <<"_" << wOut << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
			setName ( name.str() );
			REPORT(DEBUG, "Building " << name.str() );
		}

		if(wOut<0)
			THROWERROR("in standalone constructor: negative wOut");

		parentOp=this;
		multiplierUid=parentOp->getNewUId();
		xname="X";
		yname="Y";
		
		initialize();
		int g = neededGuardBits(parentOp->getTarget(), wXdecl, wYdecl, wOut);
		int possibleOutputs=1;
		if(g>=0) {
			lsbWeightInBitHeap=wOut+g-wFullP; // Should be negative; # truncated bits is the opposite of this number in this case
			possibleOutputs=2;
		}
		else
			lsbWeightInBitHeap=0;
		REPORT(DEBUG,"g=" << g << "  lsbWeightInBitHeap=" << lsbWeightInBitHeap << "   possibleOutputs=" << possibleOutputs);

		// Set up the IO signals
		addInput ( xname  , wXdecl, true );
		addInput ( yname  , wYdecl, true );
		addOutput ( "R"  , wOut, 2 , true );

		// Set up the VHDL library style // TODO
#if 0
		if(signedIO)
			useNumericStd_Signed();
		else
			useNumericStd_Unsigned();
#else
		useNumericStd();
#endif
		

		// The bit heap
		bitHeap = new BitHeap(this, wOut+g, enableSuperTiles);

		// TODO CHECK ??? A bit heap is sign-agnostic. Commented out, these methods should disapear from BitHeap
		// bitHeap->setSignedIO(signedIO);

		// initialize the critical path
		setCriticalPath(getMaxInputDelays ( inputDelays_ ));

		fillBitHeap();

		// For a standalone operator, we add the rounding-by-truncation bit, 
		// The following turns truncation into rounding, except that the overhead is large for small multipliers.
		if(lsbWeightInBitHeap<0)	{
			int weight = -lsbWeightInBitHeap-1;
			if(negate)
				bitHeap->subConstantOneBit(weight);
			else
				bitHeap->addConstantOneBit(weight);
		}
		
		bitHeap -> generateCompressorVHDL();			

		vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << range(wOut+g-1, g) << ";" << endl;
	}






	void  IntMultiplier::fillBitHeap()	{
		Plotter* plotter= bitHeap->getPlotter();
		///////////////////////////////////////
		//  architectures for corner cases   //
		///////////////////////////////////////
		
		// TODO Support negate in all the corner cases
		// To manage stand-alone cases, we just build a bit-heap of max height one, so the compression will do nothing

		// This is for legacy use of wOut -- in principle it could be removed
		int wOut=wX+wY;
		if(lsbWeightInBitHeap<0) // truncation
			wOut+=lsbWeightInBitHeap;
		
		
		// The really small ones fit in one or two LUTs and that's as small as it gets  
		if(tabulatedMultiplierP(parentOp->getTarget(), wX, wY))	 {
			vhdl << tab << "-- Ne pouvant me fier a mon raisonnement, j'ai appris par coeur le résultat de toutes les multiplications possibles" << endl;
			SmallMultTable *t = new SmallMultTable(  parentOp->getTarget(), wX, wY, wOut, negate, signedIO, signedIO);
			t->addToGlobalOpList();
			//This table is either exact, or correctly rounded if wOut<wX+wY
			// FIXME the offset is probably wrong
			vhdl << tab << declare(addUID("XY"), wX+wY) << " <= "<<addUID("YY")<<" & "<<addUID("XX")<<";"<<endl;
			inPortMap(t, "X", addUID("XY"));
			outPortMap(t, "Y", addUID("RR"));
			vhdl << instance(t, "multTable");
			useSoftRAM(t);
#if 0 // commented by Florent who didn't know what to do with the g
			 plotter->addSmallMult(0,0,wX,wY);
			 plotter->plotMultiplierConfiguration(getName(), localSplitVector, wX, wY, wOut, g);
#endif
			// Copy all the output bits of the multiplier to the bit heap if they are positive
			for (int w=0; w<wFullP; w++) 	{
				int wBH = w+lsbWeightInBitHeap;
				// REPORT(FULL, "w=" << w <<  "  wBH=" << wBH);
				if(wBH >= 0) 
					bitHeap->addBit(wBH, join(addUID("RR"), of(wBH))); 
			}		
			return;
		}
		
		// Multiplication by 1-bit integer is simple
		if ((wY == 1))		{
			vhdl << tab << "-- How to obfuscate multiplication by 1 bit: first generate a trivial bit vector" << endl;
			if (signedIO){
				manageCriticalPath(  parentOp->getTarget()->localWireDelay(wX) +  parentOp->getTarget()->adderDelay(wX+1) );				
				vhdl << tab << declare(addUID("RR"), wX+wY)  << " <= (" << zg(wX+1)  << " - ("<<addUID("XX")<< of(wX-1) 
				     << " & "<<addUID("XX")<<")) when "<<addUID("YY")<<"(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	
			}
			else {
				manageCriticalPath(  parentOp->getTarget()->localWireDelay(wX) +  parentOp->getTarget()->lutDelay() );
				vhdl << tab  << declare(addUID("RR"), wX+wY) << " <= (\"0\" & "<<addUID("XX")<<") when "<< addUID("YY")<<"(0)='1' else "<<zg(wX+1,0)<<";"<<endl;	
			}
			vhdl << tab << "-- then send its relevant bits to a useless bit heap " << endl;
			for (int w=0; w<wFullP; w++) 	{
				int wBH = w+lsbWeightInBitHeap;
				if(wBH >= 0) 
					bitHeap->addBit(wBH, join(addUID("RR"), of(wBH))); 
			}		
			vhdl << tab << "-- then compress this height-1 bit heap by doing nothing" << endl;
			outDelayMap[addUID("R")] = getCriticalPath();
			return;
		}

		// Multiplication by 2-bit integer is one addition, which is delegated to BitHeap compression anyway
		// TODO this code mostly works but it is large and unoptimal (adding 0s to bit heap)
		if ((wY == 2))		{
			string x=addUID("XX");
			string y=addUID("YY");

			vhdl << tab << declare(addUID("R0"),wX+2) << " <= (";
			if (signedIO) 
				vhdl << x << of(wX-1) << " & "<< x << of(wX-1);  
			else  
				vhdl << "\"00\"";
			vhdl <<  " & "<< x <<") when "<< y <<"(0)='1'    else "<<zg(wX+2,0)<<";"<<endl;	

			vhdl << tab << declare(addUID("R1i"),wX+2) << " <= ";
			if (signedIO) 
				vhdl << "("<< x << of(wX-1) << "  &  " << x <<" & \"0\")";
			else  
				vhdl << "(\"0\" & "<< x <<" & \"0\")";
			vhdl << " when "<< y <<"(1)='1' else " << zg(wX+2,0) << ";"<<endl;	

			vhdl << tab << declare(addUID("R1"),wX+2) << " <= ";
			if (signedIO) 
				vhdl << "not "<<addUID("R1i")<<";" <<endl;
			else  
				vhdl << addUID("R1i")<<";"<<endl;


			for (int w=0; w<wFullP; w++) 	{
				int wBH = w+lsbWeightInBitHeap;
				if(wBH >= 0) {
					bitHeap->addBit(wBH, join(addUID("R0"), of(w))); 
					bitHeap->addBit(wBH, join(addUID("R1"), of(w))); 
				}	
			}	

			// carry in bit for signed inputs
			if(signedIO)
				bitHeap->addConstantOneBit(0);
			// and that's it
			return;
		}
		
		// Multiplications that fit directly into a DSP
		int dspXSize, dspYSize;
		
		parentOp->getTarget()->getDSPWidths(dspXSize, dspYSize, signedIO);
			
		//correct the DSP sizes for Altera targets
		if(parentOp->getTarget()->getVendor() == "Altera")		{
			if(dspXSize >= 18)		{
				if(signedIO)
					dspXSize = 18;
				else
					dspXSize = 17;
			}
			if(dspYSize >= 18)	{
				if(signedIO)
					dspYSize = 18;
				else
					dspYSize = 17;
			}
		}
		
		// if we are using at least SMALL_MULT_RATIO of the DSP, then just implement 
		//	the multiplication in a DSP, without passing through a bitheap
		//	the last three conditions ensure that the multiplier can actually fit in a DSP
		//if((1.0*wX*wY >= 1.0*SMALL_MULT_RATIO*dspXSize*dspYSize) && (1.0*wX*wY < 1.0*dspXSize*dspYSize) && (wX <= dspXSize) && (wY <= dspYSize))
		if(worthUsingOneDSP(wX, wY, 0, 0, dspXSize, dspYSize) && (wX <= dspXSize) && (wY <= dspYSize))	{
			ostringstream s, zerosXString, zerosYString, zerosYNegString;
			int zerosX = dspXSize - wX + (signedIO ? 0 : 1);
			int zerosY = dspYSize - wY + (signedIO ? 0 : 1);
			int startingIndex, endingIndex;

			if(zerosX<0)
				zerosX=0;
			if(zerosY<0)
				zerosY=0;
				
			//sign extension of the inputs (or zero-extension, if working with unsigned numbers)
			if(signedIO)	{
				//sign extension
				zerosXString << "(";
				for(int i=0; i<zerosX; i++)
						zerosXString << addUID("XX") << of(wX-1) << (i!=(zerosX-1) ? " & " : ")");
				zerosYString << "(";
				zerosYNegString << "(";
				for(int i=0; i<zerosY; i++)	{
					zerosYString << addUID("YY") << of(wY-1) << (i!=(zerosY-1) ? " & " : ")");
					//zerosYNegString << addUID("YY") << "_neg" << of(wY-1) << (i!=(zerosY-1) ? " & " : ")");
					zerosYNegString << addUID("YY") << of(wY-1) << (i!=(zerosY-1) ? " & " : ")");
				}
			}
			else {// ! signedIO
				//zero extension
				zerosXString << "(" << zg(zerosX) << ")";
				zerosYString << "(" << zg(zerosY) << ")";
				zerosYNegString << "(" << zg(zerosY) << ")";
			}
				
			//if negated, the product becomes -xy = not(y)*x + x
			//	TODO: this should be more efficient than negating the product at
			//	the end, as it should be implemented in a single DSP, both the multiplication and the addition
			if(negate)
				vhdl << tab << declare(join(addUID("YY"), "_neg"), wY+zerosY) << " <= " << (zerosY>0 ? join(zerosYNegString.str(), " & ") : "") << addUID("YY") << ";" << endl;

			//manage the pipeline
			manageCriticalPath(parentOp->getTarget()->DSPMultiplierDelay());
			s << "DSP_Res_" <<  getuid();

			vhdl << tab << declare(s.str(), dspXSize+dspYSize+(signedIO ? 0 : 2))
						 << " <= (" << (zerosX>0 ? join(zerosXString.str(), " & ") : "") << addUID("XX") << ")"
						 << " *";
			if(negate)
				vhdl	 <<	" (" << addUID("YY") << "_neg);" << endl;
			else
				vhdl	 << " (" << (zerosY>0 ? join(zerosYString.str(), " & ") : "") << addUID("YY") << ");" << endl;
			
			
			//manage the pipeline: TODO
			//syncCycleFromSignal(s.str());
			
			//add the bits of x*(not y) (respectively x*y, when not negated)
			if(signedIO){
				s.str("");
				s << "not DSP_mult_" << getuid();
			}

			//FIXME: the position where the bits are added in the bitheap is no longer the same as the index of the multiplier output
			for (int w=0; w<wFullP; w++)
			{
				int wBH = w+lsbWeightInBitHeap;
				if(wBH >= 0) {
					//old version
					//bitHeap->addBit(wBH, join(s.str(), of(wBH)));
					bitHeap->addBit(wBH, join(s.str(), of(w)));
				}
			}
						
#if 0// TODO
			//keep sign-extending, if necessary
			if((bitHeap->getMaxWeight()-(endingIndex-1-startingIndex) > 1) && signedIO)
				for(int w=endingIndex-1-startingIndex; w<(int)bitHeap->getMaxWeight(); w++)
					bitHeap->addConstantOneBit(w);
			
			//add the bits of x, (not x)<<2^wY, 2^wY
			if(negate)	{
				//add x
				endingIndex	  = wX;
				startingIndex = 0+(wX+wY-wOut-g);
				for(int w=startingIndex; w<endingIndex; w++)
					if(w-startingIndex >= 0)
						bitHeap->addBit(w-startingIndex, join(addUID("XX"), of(w)));
					
				//x, when added, should be sign-extended
				if(signedIO)
					for(int w=endingIndex-startingIndex; w<(int)bitHeap->getMaxWeight(); w++)
						if(w-startingIndex >= 0)
							bitHeap->addBit(w, join(addUID("XX"), of(wX-1)));
			
				if(!signedIO)	{
					//add (not x)<<2^(wY+1)
					endingIndex	  = wX+wY-(wX+wY-wOut-g);
					startingIndex = wY-(wX+wY-wOut-g);
					for(int w=startingIndex; w<endingIndex; w++)		{
						ostringstream s;						
						s << "not(" << addUID("XX") << of(w-startingIndex) << ")";
						if((w >= 0) && (w < wOut+g))
							bitHeap->addBit(w, s.str());
					}
				 	//add 2^(wY+1)
					if(wY+1-(wX+wY-wOut-g) >= 0)
						bitHeap->addConstantOneBit(wY+1-(wX+wY-wOut-g));
				}
			}
#endif // TODO replug it 
			// this should be it
			return;
		}

		// Now getting more and more generic

		// Finish the negation of the smaller input by adding X (since -yx=not(y)x +x)

		//		setCycle(0); // TODO F2D FIXME for the virtual multiplier case where inputs can arrive later
		// setCriticalPath(initialCP);
#if 0 // TODO re-plug it 
		if(negate)		{
			vhdl << tab << "-- Finish the negation of the product (-xy as x((111..111-y)+111..111)) by adding XX + 2^wY.not XX +  2^wY" << endl;

			// Adding XX
			for(int i=0; i<wX; i++)
			{
				int w = lsbWeightInBitHeap + i-truncatedBits;
				if(w>=0)
				{
					ostringstream rhs;
					if(signedIO && i==wX-1){
						rhs << addUID("not XX") << of(i);
						// sign extension
						for(unsigned int j=i; j<bitHeap->getMaxWeight(); j++)
							bitHeap->addConstantOneBit(j); 
					}
					else
						rhs << addUID("XX") << of(i);
					bitHeap->addBit(w, rhs.str());
				}
			}

			// Adding 2^(wY+1) not XX
			for(int i=0; i<wX; i++)
			{
				int w = wY + lsbWeightInBitHeap + i-truncatedBits;
				if(w>=0 && w<(int)bitHeap->getMaxWeight())
				{
					ostringstream rhs;
					rhs << "not " << addUID("XX") << of(i);
					bitHeap->addBit(w, rhs.str());
				}
			}
			int w = wY + lsbWeightInBitHeap -truncatedBits;
			if(w>=0 && w<(int)bitHeap->getMaxWeight())
				bitHeap->addConstantOneBit(wY + lsbWeightInBitHeap -truncatedBits); 

		}
#endif

		if(useDSP) 	{
			REPORT(DETAILED,"useDSP");
			parentOp->getTarget()->getDSPWidths(wxDSP, wyDSP, signedIO);
			buildTiling();
		}
		else 	{
			// This target has no DSP, going for a logic-only implementation	
			buildLogicOnly();
		}


#if GENERATE_SVG
		plotter->plotMultiplierConfiguration(getName(), localSplitVector, wX, wY, wOut, g);
#endif
	}
	



	/**************************************************************************/
	void IntMultiplier::buildLogicOnly() 
	{
		buildHeapLogicOnly(0,0,wX,wY);
	}


	/**************************************************************************/
	void IntMultiplier::buildTiling() 
	{
#if 0
		int* multiplierWidth;
		int size;	
		
		if( parentOp->getTarget()->getVendor() == "Altera")
		{
			if ( parentOp->getTarget()->getID()=="StratixII")
			{
				StratixII* t = (StratixII*) parentOp->getTarget();
				multiplierWidth = t->getDSPMultiplierWidths();
				size = t->getNrDSPMultiplier();	

			}else if( parentOp->getTarget()->getID()=="StratixIII")
			{
				StratixIII* t = (StratixIII*) parentOp->getTarget();
				multiplierWidth = t->getDSPMultiplierWidths();
				size = t->getNrDSPMultiplier();	

			}else if( parentOp->getTarget()->getID()=="StratixIV")
			{
				StratixIV* t = (StratixIV*) parentOp->getTarget();
				multiplierWidth = t->getDSPMultiplierWidths();
				size = t->getNrDSPMultiplier();	

			}else //add Altera chips here
			{
				StratixII* t = (StratixII*) parentOp->getTarget();
				multiplierWidth = t->getDSPMultiplierWidths();
				size = t->getNrDSPMultiplier();
			}

			for(int i=0; i<size; i++)
				multWidths.push_back(multiplierWidth[i]);

			buildAlteraTiling(0, 0, wX, wY, 0, signedIO, signedIO);
		}else
		{
			// Xilinx here
			if((!signedIO) && ((wX==41) && (wY==41)))
				buildFancy41x41Tiling();
			else
				buildXilinxTiling();
		}
#endif
	}


	//the fancy tiling is used only for a hardwired case 41 41 82 
	/***********************************************************************/
	void IntMultiplier::buildFancy41x41Tiling()
	{
		//THROWERROR("fancy tiling not implemented yet");

		stringstream inx,iny;
		inx<<addUID("XX");
		iny<<addUID("YY");


		int widthX, widthY,topx,topy;

		//topright dsp;
		widthX=wxDSP;
		widthY=wyDSP;
		topx=0;
		topy=0;
		MultiplierBlock* m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeightInBitHeap);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);


		//topleft dsp
		widthX=wyDSP;
		widthY=wxDSP;
		topx=wxDSP;
		topy=0;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeightInBitHeap);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//bottomleft dsp
		widthX=wxDSP;
		widthY=wyDSP;
		topx=wyDSP;
		topy=wxDSP;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeightInBitHeap);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//bottomright dsp
		widthX=wyDSP;
		widthY=wxDSP;
		topx=0;
		topy=wyDSP;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeightInBitHeap);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//logic

		buildHeapLogicOnly(wyDSP,wyDSP,wxDSP,wxDSP,parentOp->getNewUId());

	}


	/**************************************************************************/
	void IntMultiplier::buildHeapLogicOnly(int lsbX, int lsbY, int msbX, int msbY,int blockUid)
	{
		REPORT(DETAILED,"buildheaplogiconly called for "<<lsbX<<" "<<lsbY<<" "<<msbX<<" "<<msbY);
		Target *target= parentOp->getTarget();
		
		if(blockUid == -1)
			blockUid++;    /// ???????????????

		vhdl << tab << "-- code generated by IntMultiplier::buildHeapLogicOnly()"<< endl;
		vhdl << tab << "-- buildheaplogiconly called for lsbX=" << lsbX << " lsbY=" << lsbY << " msbX="<< msbX << " msbY="<< msbY << endl;

		int dx, dy;				// Here we need to split in small sub-multiplications
		int li = target->lutInputs();

		dx = li >> 1;
		dy = li - dx; 
		REPORT(DEBUG, "dx="<< dx << "  dy=" << dy );

		int wXX=wX;
		int wYY=wY;

		int wX = msbX-lsbX;
		int wY = msbY-lsbY;
		int chunksX =  int(ceil( ( double(wX) / (double) dx) ));
		int chunksY =  int(ceil( ( 1+double(wY-dy) / (double) dy) ));
		int sizeXPadded=dx*chunksX; 
		int sizeYPadded=dy*chunksY;
		REPORT(DEBUG, "sizeXpadded"<<sizeXPadded);	

		int padX=sizeXPadded-wX;
		int padY=sizeYPadded-wY;

		REPORT(DEBUG, "X split in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
		REPORT(DEBUG, " sizeXPadded="<< sizeXPadded << "  sizeYPadded="<< sizeYPadded);
		
		//we do more than 1 subproduct 
		// FIXME where is the else?
		if (chunksX + chunksY >= 2) 
		{
			// Padding X to the right
			vhdl << tab << declare(addUID("Xp", blockUid), sizeXPadded) << " <= "
					<< addUID("XX") << range(msbX-1,lsbX) << " & " <<zg(padX) << ";" << endl;
			REPORT(DETAILED,addUID("XX") << range(msbX-1,lsbX) << " & "<<zg(padX)<<";");

			// Padding Y to the right
			vhdl << tab << declare(addUID("Yp",blockUid), sizeYPadded)<<" <= "
					<< addUID("YY") << range(msbY-1,lsbY) << " & "<< zg(padY) << ";" << endl;
			REPORT(DETAILED,addUID("YY") << range(msbY-1,lsbY) << " & "<<zg(padY)<<";");
			
			//SPLITTING
			for (int k=0; k<chunksX ; k++)
			{
				vhdl << tab << declare(join(addUID("x",blockUid),"_",k),dx) << " <= "<< addUID("Xp",blockUid) << range((k+1)*dx-1,k*dx)<<";"<<endl;
				REPORT(DETAILED,join(addUID("x",blockUid),"_",k)<<" <= "<< addUID("Xp",blockUid) << range((k+1)*dx-1,k*dx)<<";");
			}	
			for (int k=0; k<chunksY ; k++)
			{
				vhdl << tab << declare(join(addUID("y",blockUid),"_",k),dy) << " <= " << addUID("Yp",blockUid) << range((k+1)*dy-1, k*dy)<<";"<<endl;
				REPORT(DETAILED,join(addUID("y",blockUid),"_",k)<<" <= "<< addUID("Yp",blockUid) << range((k+1)*dy-1,k*dy)<<";");
			}	

			SmallMultTable *tUU, *tSU, *tUS, *tSS;

			// In the negate case we will negate the bits coming out of this table
			tUU = new SmallMultTable( target, dx, dy, dx+dy, false /*negate*/, false /*signx*/, false/*signy*/);
			tUU->addToGlobalOpList();

			if(signedIO) 
			{ // need for 4 different tables

				tSU = new SmallMultTable( target, dx, dy, dx+dy, false, true, false );
				tSU->addToGlobalOpList();

				tUS = new SmallMultTable( target, dx, dy, dx+dy, false, false, true );
				tUS->addToGlobalOpList();


				tSS = new SmallMultTable( target, dx, dy, dx+dy, false, true, true );
				tSS->addToGlobalOpList();

			}

			setCycle(0); // TODO FIXME for the virtual multiplier case where inputs can arrive later
			setCriticalPath(initialCP);
			// SmallMultTable is built to cost this:
			manageCriticalPath(  parentOp->getTarget()->localWireDelay(chunksX) + parentOp->getTarget()->lutDelay() ) ;  
			for (int iy=0; iy<chunksY; iy++)
			{
				vhdl << tab << "-- Partial product row number " << iy << endl;
				for(int ix=0; ix<chunksX; ix++)				{
					SmallMultTable *t;
					
					if(!signedIO) 
						t=tUU;
					else
					{
						// 4  cases 
						if( ((ix==chunksX-1)&&(msbX==wXX)) && ((iy==chunksY-1)&&(msbY==wYY) ))
							t=tSS;
						else if ((ix==chunksX-1)&&(msbX==wXX)) 
							t=tSU;
						else if ((iy==chunksY-1)&&(msbY==wYY))
							t=tUS;
						else
							t=tUU; 
					}

					//smallMultTable needed only if it is on the left of the truncation line
					// was if(dx*(ix+1)+dy*(iy+1)+lsbX+lsbY-padX-padY > wFullP-wOut-g)
					if(dx*(ix+1)+dy*(iy+1)+lsbX+lsbY-padX-padY + lsbWeightInBitHeap >0)	{
						bitHeap->getPlotter()->addSmallMult(dx*(ix)+lsbX-padX, dy*(iy)+lsbY-padY,dx,dy);
						REPORT(FULL,XY(ix,iy,blockUid)<<" <= " << addUID("y",blockUid) <<"_"<< iy << " & " << addUID("x",blockUid) <<"_"<< ix << ";");
						
						vhdl << tab << declare(XY(ix,iy,blockUid), dx+dy) 
						     << " <= " << addUID("y",blockUid) <<"_"<< iy << " & " << addUID("x",blockUid) <<"_"<< ix << ";"<<endl;

						inPortMap(t, "X", XY(ix,iy,blockUid));
						outPortMap(t, "Y", PP(ix,iy,blockUid));
						vhdl << instance(t, PPTbl(ix,iy,blockUid));
						useSoftRAM(t);

						vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;

						// Two's complement management
						// There are really 2 cases:
						// If the result is known positive, ie if tUU and !negate, nothing to do
						// If the result is in two's complement  
						//    sign extend by adding ones on weights >= the MSB of the table, so its sign is propagated.
						//    Also need to complement the sign bit


						// The following comments are obsolete since we negate X at the beginning of the operator:
						// Note that even when negate and tUU, the result is known negative, but it may be zero, so its actual sign is not known statically
						// Note also that in this case (negate and tUU), the result overflows the dx+dy two's complement format.
						// This is why tUU is never negated above, and we negate it below. A bit messy, but probably the most resource-efficient


						bool resultSigned = false;
						
						if((t==tSS) || (t==tUS) || (t==tSU)) 
							resultSigned = true ;

						int maxK=t->wOut; // or, dx+dy + (t==tUU && negate?1:0); 
						int minK=0;
						//if(ix == chunksX-1)
						if ((ix == 0))
							minK+=padX;
						if((iy == 0))
							//maxK-=padY;
							minK+=padY;
							
						REPORT(FULL,"The bits will be added from mink="<<minK<<"	to maxk="<<maxK);
						REPORT(FULL,  "ix=" << ix << " iy=" << iy << "  maxK=" << maxK  << "  negate=" << negate <<  "  resultSigned="  << resultSigned );
						
						for (int k=minK; k<maxK; k++) {
							ostringstream s, nots;
							s << PP(ix,iy,blockUid) << of(k); // right hand side
							nots << "not " << s.str(); 

							int weight =  ix*dx+iy*dy+k  +lsbX+lsbY-padX-padY + lsbWeightInBitHeap;

							if(weight>=0) 
							{
								// otherwise these bits are just truncated
								if(resultSigned && (k==maxK-1)) 
								{ 
									// This is a sign bit: sign-extend it by 1/ complementing and 2/ adding constant 1s all the way up to maxweight
									REPORT(FULL, "adding neg bit " << nots.str() << " at weight " << weight); 
									bitHeap->addBit(weight, nots.str());
									REPORT(FULL,  "  adding constant ones from weight "<< weight << " to "<< bitHeap->getMaxWeight());
									for (unsigned w=weight; w<bitHeap->getMaxWeight(); w++) {
										// REPORT(DEBUG, "w=" << w);
										bitHeap->addConstantOneBit(w);
									}
								}
								else 
								{ 
									// just add the bit
									bitHeap->addBit(weight, s.str());
								}
							}
						}

						vhdl << endl;
					}
				}
				REPORT(FULL,"Exiting buildHeapLogicOnly");
			}	



		}

	}

	/** checks how many DSPs will be used in case of a tiling **/
	int IntMultiplier::checkTiling(int wxDSP, int wyDSP, int& horDSP, int& verDSP)
	{
		int widthOnX=wX;
		int widthOnY=wY;
		int horDS=0;
		int verDS=0;

		//**** how many dsps will be vertical*******************************/
		int hor=0;

		//if the multiplication is signed, the first DSP will have different size, will be bigger
		if( widthOnX>=wxDSP)
			{
				hor++;
				widthOnX-=wxDSP;
			}

		if(signedIO)	
			wxDSP--;

		//how many DSPs fits on the remaining part, without the first one
		horDS=int(ceil ( (double(widthOnX) / (double) wxDSP)))+hor;
		/***********************************************************************/


		//*** how many dsps will be horizontal**********************************/
		int ver=0;

		if( widthOnY>=wyDSP)
			{
				ver++;
				widthOnY-=wyDSP;
			}


		if(signedIO)	
			wyDSP--;

		verDS=int(ceil ( (double(widthOnY) / (double) wyDSP)))+ver;
		//***********************************************************************/

		horDSP=horDS;
		verDSP=verDS;

		return verDS*horDS;
	}


	void IntMultiplier::addExtraDSPs(int lsbX, int lsbY, int botx, int boty,int wxDSP,int wyDSP)
	{
#if 0
		REPORT(DEBUG, "in addExtraDSPs at sizeX=" << wxDSP << " and sizeY=" << wyDSP << ": lsbX=" << lsbX << " lsbY=" << lsbY << " msbX=" << botx << " msbY=" << boty);
		int topx=lsbX,topy=lsbY;
		
		//if the block is on the margins of the multipliers, then the coordinates have to be reduced.
		if(lsbX<0)
			topx=0;
		else
			topx=lsbX;	

		if(lsbY<0)
			topy=0;	
		else
			topy=lsbY;

		//if the truncation line splits the block, then the used block is smaller, the coordinates need to be updated
		// was		if((botx+boty > wFullP-wOut-g) && (topx+topy < wFullP-wOut-g))
		if((botx+boty > wFullP-wOut-g) && (topx+topy < wFullP-wOut-g))
		{
			int x=topx;
			int y=boty;
			while((x+y<wFullP-wOut-g)&&(x<botx))
			{
				x++;
				topx=x;
			}

			x=botx;
			y=topy;
			while((x+y<wFullP-wOut-g)&&(y<boty))
			{
				y++;
				topy=y;	
			}	
		}

		//now check against the DSPThreshold
		if(worthUsingOneDSP(topx,topy,botx,boty,wxDSP,wyDSP))
		{
			//worth using DSP
			stringstream inx,iny;
			inx<<addUID("XX");
			iny<<addUID("YY");

			topx=botx-wxDSP;
			topy=boty-wyDSP;

			MultiplierBlock* m;
			m = new MultiplierBlock(wxDSP,wyDSP,topx,topy,inx.str(),iny.str(),truncatedBits);
			m->setNext(NULL);		
			m->setPrevious(NULL);			
			localSplitVector.push_back(m);
			bitHeap->addMultiplierBlock(m);
			REPORT(DEBUG, "in addExtraDSPs, adding a multiplier block of size wxDSP=" << wxDSP << " and wyDSP=" << wyDSP 
			       << ": lsbX=" << topx << " lsbY=" << topy << " truncatedBits=" << truncatedBits << " weight=" << m->getWeight());
			
		}
		else
		{
			//build logic	
			if((topx<botx)&&(topy<boty))
				buildHeapLogicOnly(topx,topy,botx,boty,parentOp->getNewUId());	
		}
		REPORT(FULL, "Exiting addExtraDSPs");
#endif
	}



	/** 
	 * checks the area usage of 1 dsp according to a given block and ratio(threshold)
	 * 	ratio(threshold) = says what percentage of 1 DSP area is allowed to be lost
	 **/
	//FIXME: the unused area of the DSP is not necessarly an isosceles triangle
	//			in fact, it's not necessarly a triangle
	/*
	bool IntMultiplier::worthUsingOneDSP(int lsbX, int lsbY, int msbX, int msbY,int wxDSP,int wyDSP)
	{
	
		REPORT(DEBUG, "in checktreshhold at coordinates: lsbX="<<lsbX<<" lsbY="<<lsbY<<" msbX="<<msbX<<" msbY"<<msbY);
		int widthX=(msbX-lsbX);
		int widthY=(msbY-lsbY);
		double blockArea;
		double triangle=0.0;
		double dspArea=wxDSP*wyDSP;

		if((wFullP == wOut) || (lsbX+lsbY > wFullP-wOut-g))
		{
			int x = widthX>wxDSP ? wxDSP : widthX;
			int y = widthY>wyDSP ? wyDSP : widthY;
			
			blockArea = x*y;
			//checking according to ratio/area
			if(blockArea>=(1.0-ratio)*dspArea)
					return true;
			else
					return false;
		}else
		{
			// if the truncation line splits the block, we need to subtract the area of the lost corner
			// the triangle is the area which will be lost from the block area
			
			//computing the triangle's edge (45degree => area=(edge^2)/2)		
			int x1, x2;
			int y1, y2;
			
			x1 = (lsbX < msbX-wxDSP) ? lsbX : msbX-wxDSP;
			y1 = (lsbY < msbY-wyDSP) ? lsbY : msbY-wyDSP;
			while(x1+y1<wFullP-wOut-g)
				x1++;
			
			x2 = (lsbX < msbX-wxDSP) ? lsbX : msbX-wxDSP;
			y2 = (lsbY < msbY-wyDSP) ? lsbY : msbY-wyDSP;
			while(x2+y2<wFullP-wOut-g)
				y2++;	
			
			blockArea = (widthX>wxDSP ? wxDSP : widthX) * (widthY>wyDSP ? wyDSP : widthY);
			REPORT(DEBUG, "in worthUsingOneDSP: full blocArea=" << blockArea << " x="<<x1<<" topCornerX="<<((lsbX < msbX-wxDSP) ? lsbX : msbX-wxDSP) << " y="<<y2<<" topCornerY="<<((lsbY < msbY-wyDSP) ? lsbY : msbY-wyDSP));
			
			//computing the triangle's area
			if(lsbX+lsbY<=wFullP-wOut-g)
				triangle=((x1-((lsbX < msbX-wxDSP) ? lsbX : msbX-wxDSP))*(y2-((lsbY < msbY-wyDSP) ? lsbY : msbY-wyDSP)))/2;
			else
				triangle=0.0;
			
			//the final area which is used
			blockArea = blockArea-triangle;
			REPORT(DEBUG, "in worthUsingOneDSP: usable blockArea=" << blockArea);
			REPORT(DEBUG, "in worthUsingOneDSP: dspArea=" << dspArea);
			
			//checking according to ratio/area
			if(blockArea>=(1.0-ratio)*dspArea)
					return true;
			else
					return false;
		}
	}
	*/
	
	/** 
	 * Checks the area usage of a DSP, according to a given block and DSPThreshold(threshold)
	 * 		- DSPThreshold(threshold) = says what percentage of 1 DSP area is allowed to be lost
	 * Algorithm: compute the area which is used out of a DSP, and compute 
	 * the unused ratio from this. The used area can be a triangle, a trapeze or 
	 * a pentagon, in the truncated case, or a rectangle in the case of a full 
	 * multiplier.
	 **/
	/*
	  Definition of the DSP use threshold t:
	  Consider a submultiplier block, by definition smaller than (or equal to) a DSP in both dimensions
	  let r=(submultiplier area)/(DSP area); r is between 0 and 1
	  if r >= 1-t   then use a DSP for this block 
	  So: t=0 means: any submultiplier that does not fill a DSP goes to logic
        t=1 means: any submultiplier, even very small ones, go to DSP
	*/

	bool IntMultiplier::worthUsingOneDSP(int lsbX, int lsbY, int msbX, int msbY,int wxDSP,int wyDSP)
	{
#if 0
		REPORT(DEBUG, "in checktreshhold at coordinates: lsbX="<<lsbX<<" lsbY="<<lsbY<<" msbX="<<msbX<<" msbY"<<msbY);
		
		double intersectRightX, intersectRightY, intersectTopX, intersectTopY, intersectBottomX, intersectBottomY, intersectLeftX, intersectLeftY;
		double intersectX1, intersectY1, intersectX2, intersectY2;
		int aTopEdge, bTopEdge, cTopEdge, aBottomEdge, bBottomEdge, cBottomEdge, aRightEdge, bRightEdge, cRightEdge, aLeftEdge, bLeftEdge, cLeftEdge;
		int aTruncLine, bTruncLine, cTruncLine;
		double usefulDSPArea, totalDSPArea;
		
		//if the truncation line passes above the top on the DSP
		if((wFullP == wOut) || (lsbX+lsbY > wFullP-wOut-g))
		{
			REPORT(DEBUG, "in worthUsingOneDSP: full multiplication case, truncation line does not pass through this block");
			int x = max(lsbX,msbX-wxDSP);
			int y = max(lsbY,msbY-wyDSP);
			// REPORT(DEBUG, "*********** x=" << x << "  y=" << y);
			
			usefulDSPArea = (msbX-x)*(msbY-y);
			totalDSPArea = wxDSP*wyDSP;
			
			REPORT(DEBUG, "in worthUsingOneDSP: usable blockArea=" << usefulDSPArea << "   dspArea=" << totalDSPArea);
			
			//checking according to ratio/area
			if(usefulDSPArea >= (1.0-DSPThreshold)*totalDSPArea)
					return true;
			else
					return false;
		}
		
		// equations of the lines which bound the area to be paved
		// line equation: Ax + By + C = 0
		// A=y1-y2, B=x2-x1, C=x1y2-y1x2
		aTopEdge = 0;		//top edge
		bTopEdge = 1;
		cTopEdge = -(lsbY>msbY-wyDSP ? lsbY : msbY-wyDSP);
		
		aBottomEdge = 0;		//bottom edge
		bBottomEdge = 1;
		cBottomEdge = -msbY;
		
		aLeftEdge = 1;		//left edge
		bLeftEdge = 0;
		cLeftEdge = -msbX;
		
		aRightEdge = 1;		//right edge
		bRightEdge = 0;
		cRightEdge = -(lsbX>msbX-wxDSP ? lsbX : msbX-wxDSP);
		
		//equation of truncation line - given by the 2 points (wX-g, 0) and (0, wY-g)
		aTruncLine = g-wY;
		bTruncLine = g-wX;
		cTruncLine = (wX-g)*(wY-g);
		
		//first, assume that the truncated part is a triangle
		//	then, the two intersections are with the left and bottom edge
		
		//the left edge intersected with the truncation line
		intersectLeftX = 1.0 * (bLeftEdge*cTruncLine - bTruncLine*cLeftEdge) / (aLeftEdge*bTruncLine - aTruncLine*bLeftEdge);
		intersectLeftY = 1.0 * (aTruncLine*cLeftEdge - aLeftEdge*cTruncLine) / (aLeftEdge*bTruncLine - aTruncLine*bLeftEdge);
		
		//the bottom edge intersected with the truncation line
		intersectBottomX = 1.0 * (bBottomEdge*cTruncLine - bTruncLine*cBottomEdge) / (aBottomEdge*bTruncLine - aTruncLine*bBottomEdge);
		intersectBottomY = 1.0 * (aTruncLine*cBottomEdge - aBottomEdge*cTruncLine) / (aBottomEdge*bTruncLine - aTruncLine*bBottomEdge);
		
		//check to see if the intersection points are inside the target square
		//	intersectLeft should be above the top edge
		//	intersectBottom should be to the right of right edge
		
		//check intersectLeft
		if(intersectLeftX*aTopEdge + intersectLeftY*bTopEdge + cTopEdge < 0)
		{
			//then the intersection is the one between the truncation line and 
			//	the top edge
			intersectTopX = 1.0 * (bTopEdge*cTruncLine - bTruncLine*cTopEdge) / (aTopEdge*bTruncLine - aTruncLine*bTopEdge);
			intersectTopY = 1.0 * (aTruncLine*cTopEdge - aTopEdge*cTruncLine) / (aTopEdge*bTruncLine - aTruncLine*bTopEdge);
			
			intersectX1 = intersectTopX;
			intersectY1 = intersectTopY;
		}else
		{
			intersectX1 = intersectLeftX;
			intersectY1 = intersectLeftY;
		}
		
		//check intersectBottom
		if(intersectBottomX*aRightEdge + intersectBottomY*bRightEdge + cRightEdge < 0)
		{
			//then the intersection is the one between the truncation line and 
			//	the right edge
			intersectRightX = 1.0 * (bRightEdge*cTruncLine - bTruncLine*cRightEdge) / (aRightEdge*bTruncLine - aTruncLine*bRightEdge);
			intersectRightY = 1.0 * (aTruncLine*cRightEdge - aRightEdge*cTruncLine) / (aRightEdge*bTruncLine - aTruncLine*bRightEdge);
			
			intersectX2 = intersectRightX;
			intersectY2 = intersectRightY;
		}else
		{
			intersectX2 = intersectBottomX;
			intersectY2 = intersectBottomY;
		}
		
		//renormalize the intersection points' coordinates
		intersectX1 = intersectX1<lsbX ? lsbX : intersectX1;
		intersectY2 = intersectY2<lsbY ? lsbY : intersectY2;
		
		//compute the useful DSP area
		usefulDSPArea = ((1.0*msbX-intersectX1 + 1.0*msbX-intersectX2)*(intersectY2-intersectY1)/2.0) + (1.0*msbY-intersectY2)*(msbX-intersectX2);
		totalDSPArea = wxDSP*wyDSP;
		
		REPORT(DEBUG, "in worthUsingOneDSP: truncated multiplication case");
		REPORT(DEBUG, "in worthUsingOneDSP: usable blockArea=" << usefulDSPArea);
		REPORT(DEBUG, "in worthUsingOneDSP: dspArea=" << totalDSPArea);
		REPORT(DEBUG, "in worthUsingOneDSP: intersectX1=" << intersectX1 << " intersectY1=" << intersectY1 << " intersectX2=" << intersectX2 << " intersectY2=" << intersectY2);
		
		//test if the DSP DSPThreshold is satisfied
		if(usefulDSPArea >= (1.0-DSPThreshold)*totalDSPArea)
				return true;
		else
				return false;
#endif
	}

	
	void IntMultiplier::buildAlteraTiling(int blockTopX, int blockTopY, int blockBottomX, int blockBottomY, int multIndex, bool signedX, bool signedY)
	{
#if 0
		int dspSizeX,dspSizeY;
		int widthX, widthY;
		int lsbX, lsbY, msbX, msbY;
		int multWidthsSize = multWidths.size();
		int newMultIndex = multIndex;
		bool dspSizeNotFound = true;
		bool originalSignedX = signedX;
		bool originalSignedY = signedY;
		
		//set the size of the DSP widths, preliminary
		dspSizeX = multWidths[multWidthsSize-newMultIndex-1];
		if(signedX == false)
			dspSizeX--;
		dspSizeY = multWidths[multWidthsSize-newMultIndex-1];
		if(signedY == false)
			dspSizeY--;
			
		//normalize block coordinates
		if(blockTopX<0)
			blockTopX = 0;
		if(blockTopY<0)
			blockTopY = 0;
		if(blockBottomX<0)
			blockBottomX = 0;
		if(blockBottomY<0)
			blockBottomY = 0;
		
		REPORT(DEBUG, "-----------Call to buildAlteraTiling, at dspSizeX=" << dspSizeX << " and dspSizeY=" << dspSizeY << " with parameters  - blockTopX=" << blockTopX << " blockTopY=" << blockTopY << " blockBottomX=" << blockBottomX << " blockBottomY=" << blockBottomY << (signedX ? " signed" : " unsigned") << " by " << (signedY ? "signed" : "unsigned"));
		
		//if the whole DSP is outside the range of interest, skip over it altogether
		if((blockTopX+blockTopY<wFullP-wOut-g) && (blockBottomX+blockBottomY<wFullP-wOut-g))
		{
			REPORT(DEBUG, "" << tab << tab << "adding DSP skipped (out of range of interest) at coordinates lsbX=" << blockTopX << " lsbY=" << blockTopY << " msbX=" << blockBottomX << " msbY=" << blockBottomY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
			return;
		}
		
		//if the block is of width/height zero, then end this function call, as there is nothing to do
		widthX = (blockBottomX-blockTopX+1)/dspSizeX;
		widthY = (blockBottomY-blockTopY+1)/dspSizeY;
		if(((widthX==0) && (blockTopX==blockBottomX)) || ((widthY==0) && (blockTopY==blockBottomY)))
			return;
		
		REPORT(DEBUG, "" << tab << "Testing the best DSP size");
		//search for the largest DSP size that corresponds to the ratio
		while(dspSizeNotFound)
		{
			widthX = (blockBottomX-blockTopX+1)/dspSizeX;
			widthY = (blockBottomY-blockTopY+1)/dspSizeY;
			
			REPORT(DEBUG, "" << tab << tab << "at dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY << " and widthX=" << widthX << " widthY=" << widthY);
			
			if((widthX==0) && (widthY==0))
			{
				//if both DSP dimensions are larger, the DSP block might still fit the DSPThreshold
				if(worthUsingOneDSP(blockTopX, blockTopY, blockBottomX, blockBottomY, dspSizeX, dspSizeY))
				{
					//DSPThreshold fulfilled; search is over
					dspSizeNotFound = false;
				}else
				{
					//DSPThreshold not fulfilled; pass on to the next smaller DSP size
					if(newMultIndex == multWidthsSize-1)
					{
						dspSizeNotFound = false;
					}else
					{
						newMultIndex++;
					}
					dspSizeX = multWidths[multWidthsSize-newMultIndex-1];
					if(signedX == false)
						dspSizeX--;
					dspSizeY = multWidths[multWidthsSize-newMultIndex-1];
					if(signedY == false)
						dspSizeY--;
				}
			}
			else
			{
				//there is one dimension on which the DSP can fit
				if((widthX<=1) || (widthY<=1))
				{
					int tx = ((widthX==0 || widthY==0) ? blockBottomX-dspSizeX : blockTopX), ty = ((widthX==0 || widthY==0) ? blockBottomY-dspSizeY : blockTopY);
					int bx = blockBottomX, by = blockBottomY;
					
					//if(worthUsingOneDSP(tx, ty, bx, by, dspSizeX, dspSizeY))
					if(worthUsingOneDSP((tx<blockTopX ? blockTopX : tx), (ty<blockTopY ? blockTopY : ty), bx, by, dspSizeX, dspSizeY))
					{
						//DSPThreshold fulfilled; search is over
						dspSizeNotFound = false;
					}else
					{
						//DSPThreshold not fulfilled; pass on to the next smaller DSP size
						if(newMultIndex == multWidthsSize-1)
						{
							dspSizeNotFound = false;
						}else
						{
							newMultIndex++;
						}
						dspSizeX = multWidths[multWidthsSize-newMultIndex-1];
						if(signedX == false)
							dspSizeX--;
						dspSizeY = multWidths[multWidthsSize-newMultIndex-1];
						if(signedY == false)
							dspSizeY--;
					}
				}else
				{
					dspSizeNotFound = false;
				}
			}
		}
		
		REPORT(DEBUG, "" << tab << "DSP sizes set to dspSizeX=" << dspSizeX << " and dspSizeY=" << dspSizeY);
		
		if(signedX && signedY)
		{
			REPORT(DEBUG, "" << tab << "Initial call to buildAlteraTiling, with both parameters signed");
			
			//SxS multiplication
			//	just for the top-index (for both x and y)
			//the starting point
			//	what remains to be done: one SxS multiplication in the top corner
			//							 recursive call on the top line (UxS)
			//							 recursive call on the right line (SxU)
			//							 recursive call for the rest of the multiplication (UxU)
			
			
			/*
			 * First version: sign is handled in the DSPs
			 */
			
			//top corner, SxS multiplication
			REPORT(DEBUG, "" << tab << "adding DSP (signed by signed) at coordinates lsbX=" << blockBottomX-dspSizeX << " lsbY=" << blockBottomY-dspSizeY << " msbX=" << blockBottomX << " msbY=" << blockBottomY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
			addExtraDSPs(blockBottomX-dspSizeX, blockBottomY-dspSizeY, blockBottomX, blockBottomY, dspSizeX, dspSizeY);
			
			//top line, UxS multiplications where needed, UxU for the rest
			buildAlteraTiling(blockTopX, blockTopY, blockBottomX-dspSizeX, blockBottomY, newMultIndex, false, true);
			
			//right column, SxU multiplication where needed, UxU for the rest
			//FIXME: conditions for sign are not true -> needs to be signed on the bottom part
			buildAlteraTiling(blockBottomX-dspSizeX, blockTopY, blockBottomX, blockBottomY-dspSizeY, newMultIndex, true, false);
			
			
			
			/*
			 * Second version: sign is handled separately
			 */
			/*
			//handle top line
			//		for now just call the buildHeapLogicOnly function which should handle the rest of the work
			buildHeapLogicOnly(blockTopX, blockBottomY-1, blockBottomX, blockBottomY, parentOp->getNewUId());
			 
			//handle right column
			//		for now just call the buildHeapLogicOnly function which should handle the rest of the work
			buildHeapLogicOnly(blockBottomX-1, blockTopY, blockBottomX, blockBottomY-1, parentOp->getNewUId());
			 
			//handle the rest of the multiplication
			//	the remaining part is UxU
			buildAlteraTiling(blockTopX, blockTopY, blockBottomX-1, blockBottomY-1, newMultIndex, false, false);
			*/
		}else
		{
			//SxU, UxS, UxU multiplications - they all share the same structure,
			//	the only thing that changes being the size of the operands
			
			//start tiling
			lsbX = (blockBottomX-blockTopX < dspSizeX) ? blockTopX : blockBottomX-dspSizeX;
			lsbY = (blockBottomY-blockTopY < dspSizeY) ? blockTopY : blockBottomY-dspSizeY;
			msbX = blockBottomX;
			msbY = blockBottomY;
			
			REPORT(DEBUG, "" << tab << "Original block separated in widthX=" << widthX << " by widthY=" << widthY << " blocks");
			
			//handle the bits that can be processed at the current DSP size
			for(int i=0; i<(widthY>0 ? widthY : 1); i++)
			{
				for(int j=0; j<(widthX>0 ? widthX : 1); j++)
				{
					//readjust the sign-ness and the DSP sizes
					//	only the blocks on the border need to have signed components
					if((j!=0) && (signedX))
					{
						signedX = false;
						dspSizeX--;
						
						lsbX++;
					}
					if((widthX>1) && (j==0) && (originalSignedX) && (i!=0))
					{
						signedX = originalSignedX;
						dspSizeX++;
						
						lsbX--;
					}
					
					if((i!=0) && (signedY))
					{
						signedY = false;
						dspSizeY--;
						
						lsbY++;
					}
					
					//if the whole DSP is outside the range of interest
					if((lsbX+lsbY<wFullP-wOut-g) && (msbX+msbY<wFullP-wOut-g))
					{
						if((widthX!=0) && (j!=widthX-1))
						{
							msbX = lsbX;
							lsbX = lsbX-dspSizeX;
						}
						REPORT(DEBUG, "" << tab << tab << "adding DSP skipped (out of range of interest) at coordinates lsbX=" << blockTopX << " lsbY=" << blockTopY << " msbX=" << blockBottomX << " msbY=" << blockBottomY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
						continue;						
					}
					
					if(((widthX<=1) && (widthY<=1)) && ((blockBottomX-blockTopX <= dspSizeX) && (blockBottomY-blockTopY <= dspSizeY)))
					{
						//special case; the remaining block is less than or equal to a DSP
						lsbX = blockTopX;
						lsbY = blockTopY;
						msbX = blockBottomX;
						msbY = blockBottomY;
						
						REPORT(DEBUG, "" << tab << tab << "adding DSP (to cover the whole block) at coordinates lsbX=" << blockTopX << " lsbY=" << blockTopY << " msbX=" << blockBottomX << " msbY=" << blockBottomY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
						addExtraDSPs(blockTopX, blockTopY, blockBottomX, blockBottomY, dspSizeX, dspSizeY);
					}else
					{
						//the regular case; just add a new DSP
						if(worthUsingOneDSP((lsbX<0 ? blockTopX : lsbX), (lsbY<0 ? blockTopY : lsbY), (msbX<0 ? blockTopX : msbX), (msbY<0 ? blockTopY : msbY), dspSizeX, dspSizeY))
						{
							REPORT(DEBUG, "" << tab << tab << "DSPThreshold satisfied - adding DSP at coordinates lsbX=" << lsbX << " lsbY=" << lsbY << " msbX=" << msbX << " msbY=" << msbY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
							addExtraDSPs(lsbX, lsbY, msbX, msbY, dspSizeX, dspSizeY);
						}
						else
						{
							REPORT(DEBUG, "" << tab << tab << "DSPThreshold not satisfied - recursive call at coordinates lsbX=" << lsbX << " lsbY=" << lsbY << " msbX=" << msbX << " msbY=" << msbY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY << (signedX ? " signed" : " unsigned") << " by " << (signedY ? "signed" : "unsigned"));
							if(newMultIndex == multWidthsSize-1)
							{
								REPORT(DEBUG, "" << tab << tab << tab << "size cannot be decreased anymore; will add DSP at this size");
								if((lsbX+lsbY<wFullP-wOut-g) && (msbX+msbY<wFullP-wOut-g))
								{
									if((widthX!=0) && (j!=widthX-1))
									{
										msbX = lsbX;
										lsbX = lsbX-dspSizeX;
									}
									REPORT(DEBUG, "" << tab << tab << "adding DSP skipped (out of range of interest) at coordinates lsbX=" << blockTopX << " lsbY=" << blockTopY << " msbX=" << blockBottomX << " msbY=" << blockBottomY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
									continue;						
								}else
								{
									addExtraDSPs(lsbX, lsbY, msbX, msbY, dspSizeX, dspSizeY);
								}
							}else
							{
								REPORT(DEBUG, "" << tab << tab << tab << "size can be decreased still");
								buildAlteraTiling(lsbX, lsbY, msbX, msbY, newMultIndex+1, signedX, signedY);
							}
						}
						
						if((widthX!=0) && (j!=widthX-1))
						{
							msbX = lsbX;
							lsbX = lsbX-dspSizeX;
						}
					}
				}
				
				//pass on to the next column
				if((widthY!=0) && (i!=widthY-1))
				{
					msbX = blockBottomX;
					msbY = lsbY;
					lsbX = (blockBottomX-dspSizeX < blockTopX) ? blockTopX : blockBottomX-dspSizeX;
					lsbY = (lsbY-dspSizeY < blockTopY) ? blockTopY : lsbY-dspSizeY;
				}
			}
			
			REPORT(DEBUG, "" << tab << tab << tab << "last DSP added at lsbX=" << lsbX << " lsbY=" << lsbY << " msbX=" << msbX << " msbY=" << msbY);
			
			//handle the bottom leftover bits, if necessary
			if((lsbY>0) && (lsbY != blockTopY))
			{
				REPORT(DEBUG, "" << tab << tab << "handling the bottom leftover bits at coordinates lsbX=" << lsbX << " lsbY=" << blockTopY << " msbX=" << blockBottomX << " msbY=" << lsbY << (signedX ? " signed" : " unsigned") << " by " << (signedY ? "signed" : "unsigned"));
				if((lsbX+blockTopY<wFullP-wOut-g) && (blockBottomX+lsbY<wFullP-wOut-g))
				{
					REPORT(DEBUG, "" << tab << tab << tab << "adding DSP skipped (out of range of interest) at coordinates lsbX=" << lsbX << " lsbY=" << blockTopY << " msbX=" << blockBottomX << " msbY=" << lsbY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
				}else
				{
					buildAlteraTiling(lsbX, blockTopY, blockBottomX, lsbY, newMultIndex, originalSignedX, false);
				}
			}
				
			//handle the left-side leftover bits, if necessary
			if((lsbX>0) && (lsbX != blockTopX))
			{
				REPORT(DEBUG, "" << tab << tab << "handling the left-side leftover bits at coordinates lsbX=" << blockTopX << " lsbY=" << blockTopY << " msbX=" << lsbX << " msbY=" << blockBottomY << (signedX ? " signed" : " unsigned") << " by " << (signedY ? "signed" : "unsigned"));
				if((blockTopX+blockTopY<wFullP-wOut-g) && (lsbX+blockBottomY<wFullP-wOut-g))
				{
					REPORT(DEBUG, "" << tab << tab << tab << "adding DSP skipped (out of range of interest) at coordinates lsbX=" << blockTopX << " lsbY=" << blockTopY << " msbX=" << lsbX << " msbY=" << blockBottomY << " dspSizeX=" << dspSizeX << " dspSizeY=" << dspSizeY);
				}else
				{
					buildAlteraTiling(blockTopX, blockTopY, lsbX, blockBottomY, newMultIndex, false, originalSignedY);
				}
			}
			
			REPORT(DEBUG, "-----------End of call to buildAlteraTiling, at dspSizeX=" << dspSizeX << " and dspSizeY=" << dspSizeY << " with parameters  - blockTopX=" << blockTopX << " blockTopY=" << blockTopY << " blockBottomX=" << blockBottomX << " blockBottomY=" << blockBottomY << (originalSignedX ? " signed" : " unsigned") << " by " << (originalSignedY ? "signed" : "unsigned"));
		}
#endif
	}







	void IntMultiplier::buildXilinxTiling()
	{
#if 0
		int widthXX,widthX;//local wxDSP
		int widthYY,widthY;//local wyDSP
		int hor1,hor2,ver1,ver2;	
		int horizontalDSP,verticalDSP;
		int nrDSPvertical=checkTiling(wyDSP,wxDSP,hor1,ver1); //number of DSPs used in the case of vertical tiling
		int nrDSPhorizontal=checkTiling(wxDSP,wyDSP,hor2,ver2);//number of DSPs used in case of horizontal tiling
		int botx=wX;
		int boty=wY;
		int topx,topy;

		//decides if a horizontal tiling will be used or a vertical one
		if(nrDSPvertical<nrDSPhorizontal)
		{
			widthXX=wyDSP;
			widthYY=wxDSP;
			horizontalDSP=hor1;
			verticalDSP=ver1;
		}
		else
		{
			widthXX=wxDSP;
			widthYY=wyDSP;
			horizontalDSP=hor2;
			verticalDSP=ver2;
		}


		//applying the tiles
		for(int i=0;i<verticalDSP;i++)
		{
			//managing the size of a DSP according to its position if signed
			if((signedIO)&&(i!=0))
				widthY=widthYY-1;
			else
				widthY=widthYY;	

			topy=boty-widthY;
			botx=wX;

			for(int j=0;j<horizontalDSP;j++)
			{
				//managing the size of a DSP according to its position if signed
				if((signedIO)&&(j!=0))
					widthX=widthXX-1;
				else
					widthX=widthXX;

				topx=botx-widthX;

				if(botx+boty>wFullP-wOut-g)
					addExtraDSPs(topx,topy,botx,boty,widthX,widthY);
				botx=botx-widthX;			
			}

			boty=boty-widthY;
		}
		REPORT(FULL, "Exiting buildXilinxTiling");
#endif
	}


	IntMultiplier::~IntMultiplier() {
	}


	//signal name construction

	string IntMultiplier::addUID(string name, int blockUID)
	{
		ostringstream s;
		s << name << "_m" << multiplierUid;
		if (blockUID!=-1) 
			s << "b"<< blockUID;
		return s.str() ;
	};

	string IntMultiplier::PP(int i, int j, int uid ) {
		std::ostringstream p;		
		if(uid==-1) 
			p << "PP" <<  "_X" << i << "Y" << j;
		else
			p << "PP" <<uid<<"X" << i << "Y" << j;
		return  addUID(p.str());
	};

	string IntMultiplier::PPTbl(int i, int j, int uid ) {
		std::ostringstream p;		
		if(uid==-1) 
			p << addUID("PP") <<  "_X" << i << "Y" << j << "_Tbl";
		else
			p << addUID("PP") <<"_"<<uid<<"X" << i << "Y" << j << "_Tbl";
		return p.str();
	};


	string IntMultiplier::XY( int i, int j,int uid) {
		std::ostringstream p;		
		if(uid==-1) 
			p  << "Y" << j<< "X" << i;
		else
			p  << "Y" << j<< "X" << i<<"_"<<uid;
		return  addUID(p.str());	
	};






	IntMultiplier::SmallMultTable::SmallMultTable(Target* target, int dx, int dy, int wO, bool negate, bool  signedX, bool  signedY ) : 
		Table(target, dx+dy, wO, 0, -1, true), // logic table
		dx(dx), dy(dy), wO(wO), negate(negate), signedX(signedX), signedY(signedY) {
		ostringstream name; 
		srcFileName="LogicIntMultiplier::SmallMultTable";
		// No getUid() in the name: this kind of table should be added to the globalOpList 
		name <<"SmallMultTable"<< (negate?"M":"P") << dy << "x" << dx << "r" << wO << (signedX?"Xs":"Xu") << (signedY?"Ys":"Yu");
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
		//if(!negate && signedX && signedY) cerr << "  y=" << y << "  x=" << x;
		r = x * y;
		//if(!negate && signedX && signedY) cerr << "  r=" << r;
		if(negate)
			r=-r;
		//if(negate && signedX && signedY) cerr << "  -r=" << r;
		if ( r < 0)
			r += mpz_class(1) << wF; 
		//if(!negate && signedX && signedY) cerr << "  r2C=" << r;

		if(wOut<wF){ // wOut is that of Table
			// round to nearest, but not to nearest even
			int tr=wF-wOut; // number of truncated bits 
			// adding the round bit at half-ulp position
			r += (mpz_class(1) << (tr-1));
			r = r >> tr;
		}

		//if(!negate && signedX && signedY) cerr << "  rfinal=" << r << endl;

		return r;

	};


	// Interestingly, this function is called only for a standalone constructor.
	void IntMultiplier::emulate (TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svR;

		if(!signedIO)
		{
			svR = svX * svY;
		}
		else
		{
			// Manage signed digits
			mpz_class big1 = (mpz_class(1) << (wXdecl));
			mpz_class big1P = (mpz_class(1) << (wXdecl-1));
			mpz_class big2 = (mpz_class(1) << (wYdecl));
			mpz_class big2P = (mpz_class(1) << (wYdecl-1));

			if(svX >= big1P)
				svX -= big1;

			if(svY >= big2P)
				svY -= big2;

			svR = svX * svY;
		}
		
		if(negate)
			svR = -svR;

		// manage two's complement at output
		if(svR < 0)	 
			svR += (mpz_class(1) << wFullP); 

		if(wOut>=wFullP) 
			tc->addExpectedOutput("R", svR);
		else	{
			// there is truncation, so this mult should be faithful
			svR = svR >> (wFullP-wOut);
			tc->addExpectedOutput("R", svR);
			svR++;
			svR &= (mpz_class(1) << (wOut)) -1;
			tc->addExpectedOutput("R", svR);
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




}
