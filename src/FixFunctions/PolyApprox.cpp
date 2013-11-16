/*
TODO cleanup, this is a copy of the former FixFunction

  PolyApprox object for FloPoCo

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  

  All rights reserved.

*/

#include "PolyApprox.hpp"
#include <sstream>

namespace flopoco{


	PolyApprox::PolyApprox(string sollyaString_, double targetAccuracy, int addGuardBitsToConstant)
	{
 		ostringstream completeDescription;
		completeDescription << sollyaString_  ; // << " on [" << printMPFR(xmin) << ", " << printMPFR(xmax) << "]"; 
		description = completeDescription.str();
		srcFileName = "PolyApprox";  // useful only for reporting
		uniqueName_=description;  // useful only for reporting

		// Now do the parsing in Sollya
		fS= sollya_lib_parse_string(sollyaString_.c_str());

		/* If  parse error throw an exception */
		if (sollya_lib_obj_is_error(fS))
			THROWERROR("Unable to parse input function.");

		REPORT(DEBUG, "Function: " << description << " successfully parsed");

		// This will be the LSB of the constant (unless extended below)
		lsb = floor(log2(targetAccuracy));
		REPORT(DEBUG, "InitialLSB=" << lsb);
		// a few constant objects
		sollya_obj_t fixedS = sollya_lib_fixed();
		sollya_obj_t absoluteS = sollya_lib_absolute();
		sollya_obj_t rangeS = sollya_lib_parse_string("[0;1]");

		// Accuracy has to be converted to sollya objects
		sollya_obj_t targetAccuracyS = sollya_lib_constant_from_double(targetAccuracy);


		// initial evaluation of the required degree
		//guessdegree(f,I,eps,w,bound ) : (function, range, constant, function, constant) → range
		sollya_obj_t degreeIntervalS = sollya_lib_guessdegree(fS, rangeS, targetAccuracyS, NULL); 
		sollya_obj_t degreeS = sollya_lib_inf(degreeIntervalS);
		sollya_obj_t degreeSupS = sollya_lib_sup(degreeIntervalS);
		sollya_lib_get_constant_as_int(&degree, degreeS);
		int degreeSup;
		sollya_lib_get_constant_as_int(&degreeSup, degreeSupS);
		if(DEBUG <= verbose)
			sollya_lib_printf("> PolyApprox::buildPolyApprox: initial degree of poly approx is %b\n", degreeIntervalS);

		// A few lines to add guard bits to the constant, it will be for free in terms of evaluation
		double constCoeffAccuracy = targetAccuracy;
		if(-1 == addGuardBitsToConstant) {
			// The following assumes faithful multiplications in a Horner scheme
			constCoeffAccuracy /= degree;
		}
		else if (addGuardBitsToConstant>0) {
			// the caller provided the number of bits to add in addGuardBitsToConstants
			constCoeffAccuracy /= addGuardBitsToConstant;
		}
		int constLSB = floor(log2(constCoeffAccuracy));

		bool success=false;
		bool tryReducingLSB=true;
		while(not success) {
			REPORT(DEBUG, "Trying to build coefficients with LSB=" << lsb << "   (LSB of the constant=" << constLSB << ")");
			// Managing memory with Sollya library sucks, it would probably be safer to we just do a big string and parse it. 
			// Build the list of coefficient LSBs for fpminimax
			ostringstream s;
			s << "[|";
			for(int i=0; i<=degree ; i++) {
					s << (i==0? -constLSB :  -lsb);
					if(i<degree) s<< ",";
			}
			s << "|]";
			sollya_obj_t coeffSizeListS = sollya_lib_parse_string(s.str().c_str());
			if(DEBUG <= verbose) {
				sollya_lib_printf("> PolyApprox::buildPolyApprox: f=%b,   degreeS= %b   coeffSizeList = %b  rangeS = %b,   fixedS = %b,   absoluteS = %b\n", 
													fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS);
				sollya_lib_printf(">   fpminimax(%b, %b, %b, %b, %b, %b);\n", 
													fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS);
				
			}
			// Tadaaa! After all this useless noise we may launch fpminimax	
			polynomialS = sollya_lib_fpminimax(fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS, NULL);
			sollya_lib_clear_obj(coeffSizeListS);
			if(DEBUG <= verbose)
				sollya_lib_printf("> PolyApprox::buildPolyApprox: obtained polynomial   %b\n", polynomialS);

			// Checking its approximation error; 
			sollya_obj_t supNormAccS = sollya_lib_parse_string("1b-10"); // This is the size of the returned interval... 10^-3 should be enough for anybody
			if(DEBUG <= verbose) {
				sollya_lib_printf(">   supnorm(%b, %b, %b, %b, %b);\n", 
													polynomialS, fS, rangeS, absoluteS, supNormAccS);
			}
	    sollya_obj_t supNormRangeS = sollya_lib_supnorm(polynomialS, fS, rangeS, absoluteS, supNormAccS);
			sollya_lib_clear_obj(supNormAccS);

			sollya_obj_t supNormS = sollya_lib_sup(supNormRangeS);
			sollya_lib_clear_obj(supNormRangeS);
			
			sollya_lib_get_constant_as_double(& approxError, supNormS);
			sollya_lib_clear_obj(supNormS);
			
			REPORT(DEBUG, "Polynomial accuracy is " << approxError);

			// Now comes the test: did we success in getting an accurate enough polynomial?
			if(approxError < targetAccuracy) {
				REPORT(DEBUG, "Polynomial is accurate enough");
				success=true; 
			}
			else {
				success=false;
				// put this polynomial to the recycle bin
				sollya_lib_clear_obj(polynomialS);
				REPORT(DEBUG, "Polynomial is NOT accurate enough");

				if(tryReducingLSB) {
					lsb-=1;
					constLSB-=1;
					tryReducingLSB=false;
					REPORT(DEBUG, "  ... pushing LSB to " << lsb << " and starting over");
				}
				else { // OK, we tried increasing LSB once and it didn't work. Maybe we should increase degree?
					if (degreeSup>degree){
						// restore LSB
						lsb+=1;
						constLSB+=1;
						// and increase degree
						degree++;
						sollya_lib_clear_obj(degreeS);
						degreeS = sollya_lib_constant_from_int(degree);
					}
					else{ // guessDegree seemed sure the degree should not be larger than that, let's keep trying with the LSB
						tryReducingLSB=true;
					}
				}
			}
		} // exit from the while loop... hopefully

		// compute the MSBs
		for (int i=0; i<=degree; i++){
			sollya_obj_t iS = sollya_lib_constant_from_int(i);
			//sollya_lib_printf("i = %d  = %b\n", i, iS);
			sollya_obj_t coeffS = sollya_lib_coeff(polynomialS, iS);
			sollya_lib_printf("c = %b \n", coeffS);
			sollya_lib_clear_obj(iS);
			double coeff;
			sollya_lib_get_constant_as_double(&coeff, coeffS);
			weightMSB.push_back(ceil(log2(fabs(coeff))));
			sollya_lib_clear_obj(coeffS);
			weightLSB.push_back((i==0? constLSB :  lsb));
		}		
		ostringstream debugstring;
		debugstring << "coeff formats (MSB, LSB): ";
		for (int i=0; i<=degree; i++){
			debugstring << "(" << weightMSB[i] << ", " << weightLSB[i] << ")";
		}
		REPORT(DEBUG, debugstring.str());

 
		// Please leave the memory in the state you would like to find it when entering
	  sollya_lib_clear_obj(fixedS);
	  sollya_lib_clear_obj(absoluteS);
	  sollya_lib_clear_obj(targetAccuracyS);
		sollya_lib_clear_obj(degreeIntervalS);
	  sollya_lib_clear_obj(degreeS);
	  sollya_lib_clear_obj(degreeSupS);
	}




		
	PolyApprox::~PolyApprox()
	{
	  sollya_lib_clear_obj(fS);
	  sollya_lib_clear_obj(polynomialS);
		//	  sollya_lib_clear_obj(S);
	}

	string PolyApprox::getDescription() const
	{
		return description;
	}

	int PolyApprox::getDegree() const
	{
		return degree;
	}

	void PolyApprox::eval(mpfr_t r, mpfr_t x) const
	{
		sollya_lib_evaluate_function_at_point(r, fS, x, NULL);
	}

	double PolyApprox::eval(double x) const
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

	sollya_obj_t PolyApprox::getSollyaPolynomial() const
	{
		return polynomialS;
	}
} //namespace
