/*
  Polynomial Function Evaluator for FloPoCo

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
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


#include "FixFunctionBySimplePoly.hpp"

#ifdef HAVE_SOLLYA

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0


	FixFunctionBySimplePoly::FixFunctionBySimplePoly(Target* target, string func, int lsbIn, int msbOut, int lsbOut, bool finalRounding_, map<string, double> inputDelays):
		Operator(target, inputDelays), finalRounding(finalRounding_){
		f=new FixFunction(func, lsbIn, msbOut, lsbOut);
		
		srcFileName="FixFunctionBySimplePoly";
		
		ostringstream name;
		name<<"FixFunctionBySimplePoly_"<<getNewUId(); 
		setName(name.str()); 

		setCopyrightString("Florent de Dinechin (2014)");
		addHeaderComment("-- Evaluator for " +  f-> getDescription() + "\n"); 
		addInput("X"  , -lsbIn);
		int outputSize = msbOut-lsbOut+1;
		addOutput("Y" ,outputSize , 2);
		useNumericStd();

		// Polynomial approximation
		double targetApproxError = pow(2,lsbOut-1); 
		poly = new BasicPolyApprox(f, targetApproxError, -1);
		double approxErrorBound = poly->approxErrorBound;

		int degree = poly->degree;
		REPORT(DEBUG, "Degree is " << degree);

		// Now building the corresponding VHDL signals
		vhdl << tab << "-- With the following polynomial, approx error bound is " << approxErrorBound << " ("<< log2(approxErrorBound) << " bits)" << endl;
		for(int i=0; i<=degree; i++) {
			FixConstant* ai = poly->coeff[i];
			coeffMSB.push_back (ai->MSB);
			coeffLSB.push_back (ai->LSB);
			coeffSize.push_back(ai->MSB - ai->LSB +1);
			//			REPORT(DEBUG, " a" << i << " = " << ai->getBitVector() << "  " << printMPFR(ai->fpValue)  );
			vhdl << tab << declareFixPoint(join("A",i), true, ai->MSB, ai->LSB)
					 << " <= " << ai->getBitVector(0 /*both quotes*/);
			if(i==0) {
				// add the final round bit that transforms the final truncation into a round
				vhdl << " + \"" << unsignedBinary(mpz_class(1)<<(lsbOut-ai->LSB-1),  ai->MSB - ai->LSB +1) << "\"; -- final round bit" <<endl;

			}
			else 
				vhdl << ";" << endl;
		}
		vhdl << tab << declareFixPoint("Xs", true, 0, lsbIn) << " <= signed('0' & X);  -- sign extension of X" << endl; 

		bool plainStupidVHDL=true;

		if(plainStupidVHDL) {
			// Here we assume all the coefficients already include the proper number of guard bits
			int sigmaMSB=coeffMSB[degree];
			int sigmaLSB=coeffLSB[degree];
			vhdl << tab << declareFixPoint(join("Sigma", 0), true, sigmaMSB, sigmaLSB) 
					 << " <= " << join("A", degree)  << ";" << endl;

			for(int i=1; i<=degree; i++) {
				// When multiplying two unsigned, the MSB is the sum of the MSBs
				// But when multiplying two signed, the MSB is the sum plus one (only used in the ultrare case -max*-max)
				vhdl << tab << declareFixPoint(join("P", i), true, sigmaMSB+0+1,  sigmaLSB  + f->lsbIn /*LSB*/) 
						 <<  " <= Xs * Sigma" << i-1 << ";" << endl;
				
				sigmaMSB = coeffMSB[degree-i]+1; // +1 to absorb addition overflow
				sigmaLSB = coeffLSB[degree-i];
				resizeFixPoint(join("Ptrunc", i), join("P", i), sigmaMSB, sigmaLSB);
				
				vhdl << tab << declareFixPoint(join("Sigma", i), true, sigmaMSB, sigmaLSB)  // sign extend the coeff
						 << " <= (" << join("A", degree-i) << of(coeffSize[degree-i]-1) << "&" << join("A", degree-i) << ") + " << join("Ptrunc", i) << ";" << endl;
			}
		}


		//Building the vector of sizes for FixHornerEvaluator
		// a0 is a bit special

		resizeFixPoint("Ys", join("Sigma", degree),  msbOut, lsbOut);

		vhdl << tab << "Y <= " << "std_logic_vector(Ys);" << endl;
	}



	FixFunctionBySimplePoly::~FixFunctionBySimplePoly() {
		free(f);
	}
	


	void FixFunctionBySimplePoly::emulate(TestCase* tc){
		f->emulate(tc);
	}

	void FixFunctionBySimplePoly::buildStandardTestCases(TestCaseList* tcl){
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
