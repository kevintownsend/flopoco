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
		f=new FixFunction(func, lsbIn, lsbOut, msbOut);
		
		srcFileName="FixFunctionBySimplePoly";
		
		ostringstream name;
		name<<"FixFunctionBySimplePoly_"<<getNewUId(); 
		setName(name.str()); 

		setCopyrightString("Florent de Dinechin (2014)");
		addHeaderComment("-- Evaluator for " +  f-> getDescription() + "\n"); 


		// Polynomial approximation
		double targetApproxError = pow(2,lsbOut-1); 
		poly = new BasicPolyApprox(f, targetApproxError, -1);
		double approxErrorBound = poly->approxErrorBound;

		int degree = poly->degree;
		REPORT(DEBUG, "Degree is " << degree);

		//Building the vector of sizes for FixHornerEvaluator
		// a0 is a bit special

		// Now building the corresponding VHDL signals
		vhdl << tab << "-- With the following polynomial, approx error bound is " << approxErrorBound << " ("<< log2(approxErrorBound) << " bits)" << endl;
		for(int i=0; i<=degree; i++) {
			FixConstant* ai = poly->coeff[i];
			//			REPORT(DEBUG, " a" << i << " = " << ai->getBitVector() << "  " << printMPFR(ai->fpValue)  );
			vhdl << tab << declare(join("A",i), ai->MSB - ai->LSB +1)
					 << " <= " << ai->getBitVector(0 /*both quotes*/) << ";" 
					 << endl;
		}


#if 1
#endif
		vhdl << tab << "-- " << endl;
	}



	FixFunctionBySimplePoly::~FixFunctionBySimplePoly() {
		free(f);
	}
	


	void FixFunctionBySimplePoly::emulate(TestCase* tc){
		f->emulate(tc, true /*correctly rounded*/);
	}
}
	
#endif //HAVE_SOLLYA	
