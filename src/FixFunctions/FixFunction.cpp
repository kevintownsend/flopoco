/*
  FixFunction object for FloPoCo

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  

  All rights reserved.

*/

#include "FixFunction.hpp"
#include <sstream>

namespace flopoco{

	FixCoefficient::FixCoefficient(sollya_obj_t coeffS) {
		mpfr_t coeffM;
		
		mpfr_init2(coeffM, 1024);
		sollya_lib_get_constant(coeffM, coeffS);
		exponent = mpfr_get_z_exp (intSignificand.get_mpz_t(), coeffM);
		cout << "Coeff: "<< intSignificand << " * 2^" << exponent << endl; 
		mpfr_clear(coeffM);
	};

	FixCoefficient::~FixCoefficient(){};



	FixFunction::FixFunction(string sollyaString_, double xmin_, double xmax_){
		mpfr_init2(xmin, 53);
		mpfr_set_d(xmin, xmin_, GMP_RNDN);
		mpfr_init2(xmax, 53);
		mpfr_set_d(xmax, xmax_, GMP_RNDN);
		finishConstruction(sollyaString_);
	};
 
	FixFunction::FixFunction(string sollyaString_, mpfr_t xmin_, mpfr_t xmax_)
	{
		mpfr_init2(xmin, mpfr_get_prec (xmin_));
		mpfr_set(xmin, xmin_, GMP_RNDN);
		mpfr_init2(xmax, mpfr_get_prec (xmax_));
		mpfr_set(xmax, xmax_, GMP_RNDN);
		finishConstruction(sollyaString_);
	}

	void FixFunction::finishConstruction(string sollyaString_)
	{
		// A bit of bookkepping
		verbose=4;
 		ostringstream completeDescription;
		completeDescription << sollyaString_  << " on [" << printMPFR(xmin) << ", " << printMPFR(xmax) << "]"; 
		description = completeDescription.str();
		srcFileName = "FixFunction";  // useful only for reporting
		uniqueName_=description;  // useful only for reporting

		// Now do the parsing in Sollya
		fS= sollya_lib_parse_string(sollyaString_.c_str());

		/* If  parse error throw an exception */
		if (sollya_lib_obj_is_error(fS))
			THROWERROR("Unable to parse input function.");

		// rangeS may as well be converted to a sollya_obj_t now
		rangeS = sollya_lib_range_from_bounds(xmin, xmax);

		REPORT(DEBUG, "Function: " << description << " successfully parsed");
		
#if 0
		/* Map the input to the [0,1[ range */
		// g(y) = scale * f(y * (xmax - xmin) + xmin)
		sollya_obj_t sXW = sollya_lib_constant_from_double(xmax - xmin);
		sollya_obj_t sXMin = sollya_lib_constant_from_double(xmin);
	
		sollya_obj_t sX = SOLLYA_ADD(SOLLYA_MUL(SOLLYA_X_, sXW), sXMin);
		sollya_obj_t sG = sollya_lib_substitute(sF, sX);
	
		node = SOLLYA_MUL(sollya_lib_constant_from_double(scale), sG);
		if (node == 0)
			throw "Sollya error when performing range mapping.";
	
		// No need to free memory for Sollya subexpressions (?)
#endif
	}

	FixFunction::~FixFunction()
	{
		mpfr_clear(xmin);
		mpfr_clear(xmax);
	  sollya_lib_clear_obj(fS);
	  sollya_lib_clear_obj(rangeS); 
		sollya_lib_clear_obj(polynomialS);
	}

	string FixFunction::getDescription() const
	{
		return description;
	}

	void FixFunction::eval(mpfr_t r, mpfr_t x) const
	{
		sollya_lib_evaluate_function_at_point(r, fS, x, NULL);
	}

	double FixFunction::eval(double x) const
	{
		mpfr_t mpX, mpR;
		double r;

		mpfr_inits(mpX, mpR, NULL);
		mpfr_set_d(mpX, x, GMP_RNDN);
		sollya_lib_evaluate_function_at_point(mpR, fS, mpX, NULL);
		r = mpfr_get_d(mpR, GMP_RNDN);
		mpfr_clears(mpX, mpR, NULL);

		return r;
	}

	sollya_obj_t FixFunction::getSollyaObj() const
	{
		return fS;
	}



	int FixFunction::buildPolyApprox(double targetAccuracy) {
		int lsb;
		// This is the LSB of the constant 
		lsb = floor(log2(targetAccuracy));
		REPORT(DEBUG, "InitialLSB=" << lsb);
		// a few constant objects
		sollya_obj_t fixedS = sollya_lib_fixed();
		sollya_obj_t absoluteS = sollya_lib_absolute();

		// Accuracy and range have to be converted to sollya objects
		sollya_obj_t targetAccuracyS = sollya_lib_constant_from_double(targetAccuracy);

		// initial evaluation of the required degree
		//guessdegree(f,I,eps,w,bound ) : (function, range, constant, function, constant) → range
		sollya_obj_t degreeIntervalS = sollya_lib_guessdegree(fS, rangeS, targetAccuracyS, NULL); 
		sollya_obj_t degreeS = sollya_lib_inf(degreeIntervalS);
		sollya_lib_clear_obj(degreeIntervalS);
		sollya_lib_get_constant_as_int(&polyApproxDegree, degreeS);
		if(DEBUG <= verbose)
			sollya_lib_printf("> FixFunction::buildPolyApprox: initial degree of poly approx is %b, I repeat, %d\n", degreeS, polyApproxDegree);


		bool success=false;
		while(not success) {
			REPORT(DEBUG, "Trying to build coefficients with LSB=" << lsb);
			// Managing memory with Sollya library sucks, it would probably be safer to we just do a big string and parse it. 
			// Build the list of coefficient LSBs for fpminimax
			ostringstream s;
			s << "[|";
			for(int i=0; i<=polyApproxDegree ; i++) {
				s << -lsb;
				// TODO ici updater lsb
				if(i<polyApproxDegree) 
					s<< ",";
			}
			s << "|]";
			sollya_obj_t coeffSizeListS = sollya_lib_parse_string(s.str().c_str());
			if(DEBUG <= verbose) {
				sollya_lib_printf("> FixFunction::buildPolyApprox: f=%b,   degreeS= %b   coeffSizeList = %b  rangeS = %b,   fixedS = %b,   absoluteS = %b\n", 
													fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS);
			}
			// Tadaaa! After all this useless noise we may launch fpminimax	
			polynomialS = sollya_lib_fpminimax(fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS, NULL);
			sollya_lib_clear_obj(coeffSizeListS);
			if(DEBUG <= verbose)
				sollya_lib_printf("> FixFunction::buildPolyApprox: obtained polynomial   %b\n", polynomialS);

			// Checking its approximation error; 
			sollya_obj_t supNormAccS = sollya_lib_parse_string("1b-10"); // This is the size of the returned interval... 10^-3 should be enough for anybody
	    sollya_obj_t supNormRangeS = sollya_lib_supnorm(polynomialS, fS, rangeS, absoluteS, supNormAccS);
			sollya_lib_clear_obj(supNormAccS);

			sollya_obj_t supNormS = sollya_lib_sup(supNormRangeS);
			sollya_lib_clear_obj(supNormRangeS);
			
			sollya_lib_get_constant_as_double(& polyApproxError, supNormS);
			sollya_lib_clear_obj(supNormS);
			
			REPORT(DEBUG, "Polynomial accuracy is " << polyApproxError);

			// Now comes the test: did we success in getting an accurate enough polynomial?
			if(polyApproxError < targetAccuracy) {
				REPORT(DEBUG, "Polynomial is accurate enough, exiting");
				success=true; 
			}
			else {
				success=false;
				lsb-=1;
				REPORT(DEBUG, "Polynomial is NOT accurate enough, pushing LSB to " << lsb << " and starting over");
				// and put this polynomial to the recycle bin
				sollya_lib_clear_obj(polynomialS);
			}
		} // exit from the while loop


		// Here we assume success...
		// Converting the polynomial to FixCoefficients
		for(int i=0; i<=polyApproxDegree; i++) {
			sollya_obj_t iS = sollya_lib_constant_from_int(i);
			sollya_lib_printf("i = %d  = %b\n", i, iS);
			sollya_obj_t coeffS = sollya_lib_coeff(polynomialS, iS);
			sollya_lib_printf("c = %b \n", coeffS);
			sollya_lib_clear_obj(iS);
			polyApprox.push_back(new FixCoefficient(coeffS));
			sollya_lib_clear_obj(coeffS);
		} 
		// Please leave the memory in the state you would like to find it when entering
		sollya_lib_clear_obj(targetAccuracyS);
 		sollya_lib_clear_obj(degreeS);
 		sollya_lib_clear_obj(fixedS);
 		sollya_lib_clear_obj(absoluteS);
		return polyApproxDegree;
	}


} //namespace
