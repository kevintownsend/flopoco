/*
  Polynomial Function Evaluator for FloPoCo
	This version uses piecewise polynomial approximation for a trade-off between tables and multiplier hardware.
	
  Authors: Florent de Dinechin (rewrite from scratch of original code by Mioara Joldes and Bogdan Pasca, see the attic directory)

  This file is part of the FloPoCo project
	launched by the Arénaire/AriC team of Ecole Normale Superieure de Lyon
  currently developed by the Socrate team at INSA de Lyon
 
  Initial software.
  Copyright © ENS-Lyon, INSA-Lyon, INRIA, CNRS, UCBL,  
  2008-2014.
  All rights reserved.

  */

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "../utils.hpp"


#include "FixFunctionByPiecewisePoly.hpp"
#include "FixFunctionByTable.hpp"
#include "FixHornerEvaluator.hpp"

#ifdef HAVE_SOLLYA

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0

	FixFunctionByPiecewisePoly::CoeffTable::CoeffTable(Target* target, int wIn, int wOut, PiecewisePolyApprox* polyApprox_, bool addFinalRoundBit_, int finalRoundBitPos_) :
		Table(target, wIn, wOut), polyApprox(polyApprox_), addFinalRoundBit(addFinalRoundBit_), finalRoundBitPos(finalRoundBitPos_)
	{
		srcFileName="FixFunctionByPiecewisePoly::CoeffTable";
		ostringstream name; 
		name <<"CoeffTable_" << wIn << "_" << wOut << "_" << getNewUId();
		setName(name.str());
		// outDelayMap["Y"] = target->RAMDelay(); // is this useful?
	}


	mpz_class FixFunctionByPiecewisePoly::CoeffTable::function(int x){
		mpz_class z=0;
		int currentShift=0;
		for(int i=polyApprox->degree; i>=0; i--) {
			z += (polyApprox-> getCoeff(x, i)) << currentShift; // coeff of degree i from poly number x
			// REPORT(DEBUG, "i=" << i << "   z=" << unsignedBinary(z, 64));
			if(i==0 && addFinalRoundBit){ // coeff of degree 0
				z += mpz_class(1)<<(currentShift + finalRoundBitPos - polyApprox->LSB); // add the round bit
				//REPORT(DEBUG, "i=" << i << " + z=" << unsignedBinary(z, 64));
				// This may entail an overflow of z, in the case -tiny -> +tiny, e.g. first polynomial of atan
				// This is OK modulo 2^wOut (two's complement), but will break the vhdl output of Table: fix it here.
				z = z & ((mpz_class(1)<<wOut) -1);
				// REPORT(INFO, "Adding final round bit at position " << finalRoundBitPos-polyApprox->LSB);
			}
			currentShift +=  polyApprox->MSB[i] - polyApprox->LSB +1;
		}
		return z;
	}





	FixFunctionByPiecewisePoly::FixFunctionByPiecewisePoly(Target* target, string func, int lsbIn, int msbOut, int lsbOut, int degree_, bool finalRounding_, bool plainStupidVHDL, float DSPThreshold, map<string, double> inputDelays):
		Operator(target, inputDelays), degree(degree_), finalRounding(finalRounding_){

		if(finalRounding==false){
			THROWERROR("FinalRounding=false not implemented yet" );
		}

		f=new FixFunction(func, lsbIn, msbOut, lsbOut); // this will provide emulate etc.
		
		srcFileName="FixFunctionByPiecewisePoly";
		
		ostringstream name;
		name<<"FixFunctionByPiecewisePoly_"<<getNewUId(); 
		setName(name.str()); 

		setCopyrightString("Florent de Dinechin (2014)");
		addHeaderComment("-- Evaluator for " +  f-> getDescription() + "\n"); 
		REPORT(DETAILED, "Entering: FixFunctionByPiecewisePoly \"" << func << "\" " << lsbIn << " " << msbOut << " " << lsbOut << " " << degree << " " << plainStupidVHDL);
 		int wX=-lsbIn;
		addInput("X", wX);
		int outputSize = msbOut-lsbOut+1; // TODO finalRounding would impact this line
		addOutput("Y" ,outputSize , 2);
		useNumericStd();

		if(degree==0){ // This is a simple table
			REPORT(DETAILED, "Degree 0: building a simple table");
			FixFunctionByTable* table=new FixFunctionByTable(target, func, lsbIn, msbOut, lsbOut);
			addSubComponent(table);

			inPortMap(table, "X", "X");
			outPortMap(table, "Y", "YR");
			vhdl << instance(table, "simpleTable") << endl;
			vhdl << tab << "Y <= YR;" << endl; 
		}
		else{

			// Build the polynomial approximation
			REPORT(DETAILED, "Computing polynomial approximation for target precision "<< lsbOut-2);
			double targetAcc= pow(2, lsbOut-2);
			polyApprox = new PiecewisePolyApprox(func, targetAcc, degree);
			int alpha =  polyApprox-> alpha; // coeff table input size 

			// Resize its MSB to the one input by the user. This is useful for functions that "touch" a power of two and may thus overflow by one bit.
			for (int i=0; i<(1<<alpha); i++) {
				// REPORT(DEBUG, "i=" << i << "   coeff 0 before resizing: " <<   polyApprox -> poly[i] -> coeff[0] ->getBitVector());
				polyApprox -> poly[i] -> coeff[0] -> changeMSB(msbOut);
				// REPORT(DEBUG, "i=" << i << "   coeff 0 after resizing: " <<   polyApprox -> poly[i] -> coeff[0] ->getBitVector());
			}
			polyApprox -> MSB[0] = msbOut;

			// Store it in a table
			int polyTableOutputSize=0;
			for (int i=0; i<=degree; i++) {
				polyTableOutputSize += polyApprox->MSB[i] - polyApprox->LSB +1;
			} 
			REPORT(DETAILED, "Poly table input size  = " << alpha);
			REPORT(DETAILED, "Poly table output size = " << polyTableOutputSize);

			// This is where we add the final rounding bit
			FixFunctionByPiecewisePoly::CoeffTable* coeffTable = new CoeffTable(target, alpha, polyTableOutputSize, polyApprox, 
																																					finalRounding, lsbOut-1 /*position of the round bit*/) ;
			addSubComponent(coeffTable);

			vhdl << tab << declare("A", alpha)  << " <= X" << range(wX-1, wX-alpha) << ";" << endl;
			vhdl << tab << declare("Z", wX-alpha)  << " <= X" << range(wX-alpha-1, 0) << ";" << endl;

			inPortMap(coeffTable, "X", "A");
			outPortMap(coeffTable, "Y", "Coeffs");
			vhdl << instance(coeffTable, "coeffTable") << endl;

			int currentShift=0;
			for(int i=polyApprox->degree; i>=0; i--) {
				vhdl << tab << declare(join("A",i), polyApprox->MSB[i] - polyApprox->LSB +1)  << " <= Coeffs" << range(currentShift + (polyApprox->MSB[i] - polyApprox->LSB), currentShift) << ";" << endl;
				currentShift +=  polyApprox->MSB[i] - polyApprox->LSB +1;
			}

			FixHornerEvaluator* horner = new FixHornerEvaluator(target, lsbIn+alpha, msbOut, lsbOut, degree, polyApprox->MSB, polyApprox->LSB, true, true, plainStupidVHDL, DSPThreshold);		
			addSubComponent(horner);

			inPortMap(horner, "X", "Z");
			outPortMap(horner, "R", "Ys");
			for(int i=0; i<=polyApprox->degree; i++) {
				inPortMap(horner,  join("A",i),  join("A",i));
			}
			vhdl << instance(horner, "horner") << endl;
		
			vhdl << tab << "Y <= " << "std_logic_vector(Ys);" << endl;
		}
	}



	FixFunctionByPiecewisePoly::~FixFunctionByPiecewisePoly() {
		free(f);
	}
	


	void FixFunctionByPiecewisePoly::emulate(TestCase* tc){
		f->emulate(tc);
	}

	void FixFunctionByPiecewisePoly::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		tc = new TestCase(this); 
		tc->addInput("X", 0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addInput("X", (mpz_class(1) << f->wIn) -1);
		emulate(tc);
		tcl->add(tc);

	}

}
	
#endif //HAVE_SOLLYA	
