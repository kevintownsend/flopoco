/*

  A class that manages polynomial approximation for FloPoCo (and possibly later for metalibm). 

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL

  All rights reserved.

*/

#include "PiecewisePolyApprox.hpp"
#include <sstream>

namespace flopoco{


	PiecewisePolyApprox::PiecewisePolyApprox(FixFunction *f_, double targetAccuracy, int addGuardBitsToConstant): 
		f(f_)
	{
		needToFreeF = false;
		buildApproxFromTargetAccuracy(targetAccuracy,  addGuardBitsToConstant);
	}


	PiecewisePolyApprox::PiecewisePolyApprox(string sollyaString_, double targetAccuracy, int addGuardBitsToConstant)
	{
		//  parsing delegated to FixFunction
		f = new FixFunction(sollyaString_);
		needToFreeF = true;
		buildApproxFromTargetAccuracy(targetAccuracy,  addGuardBitsToConstant);
	}



	PiecewisePolyApprox::PiecewisePolyApprox(FixFunction* f_, int degree){};

		
	PiecewisePolyApprox::~PiecewisePolyApprox()
	{
		if(needToFreeF)	free(f);
	  sollya_lib_clear_obj(polynomialS);
		//	  sollya_lib_clear_obj(S);
		if(coeff.size()!=0){
			for (int i=0; i<coeff.size(); i++)
				free(coeff[i]);
		}
	}




	void PiecewisePolyApprox::buildApproxFromTargetAccuracy(double targetAccuracy, int addGuardBitsToConstant)
	{
		sollya_obj_t fS = f->getSollyaObj(); // no need to free this one
		// This will be the LSB of the constant (unless extended below)
		LSB = floor(log2(targetAccuracy));
		REPORT(DEBUG, "InitialLSB=" << LSB);
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
			sollya_lib_printf("> PiecewisePolyApprox::buildPiecewisePolyApprox: initial degree of poly approx is %b\n", degreeIntervalS);

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
		constLSB = floor(log2(constCoeffAccuracy));

		bool success=false;
		bool tryReducingLSB=true;
		while(not success) {
			REPORT(DEBUG, "Trying to build coefficients with LSB=" << LSB << "   (LSB of the constant=" << constLSB << ")");
			// Managing memory with Sollya library sucks, it would probably be safer to we just do a big string and parse it. 
			// Build the list of coefficient LSBs for fpminimax
			ostringstream s;
			s << "[|";
			for(int i=0; i<=degree ; i++) {
					s << (i==0? -constLSB :  -LSB);
					if(i<degree) s<< ",";
			}
			s << "|]";
			sollya_obj_t coeffSizeListS = sollya_lib_parse_string(s.str().c_str());
			if(DEBUG <= verbose) {
				sollya_lib_printf("> PiecewisePolyApprox::buildPiecewisePolyApprox: f=%b,   degreeS= %b   coeffSizeList = %b  rangeS = %b,   fixedS = %b,   absoluteS = %b\n", 
													fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS);
				sollya_lib_printf(">   fpminimax(%b, %b, %b, %b, %b, %b);\n", 
													fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS);
				
			}
			// Tadaaa! After all this useless noise we may launch fpminimax	
			polynomialS = sollya_lib_fpminimax(fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS, NULL);
			sollya_lib_clear_obj(coeffSizeListS);
			if(DEBUG <= verbose)
				sollya_lib_printf("> PiecewisePolyApprox::buildPiecewisePolyApprox: obtained polynomial   %b\n", polynomialS);

			// Checking its approximation error;
			sollya_obj_t supNormS; // it will end up there 
			sollya_obj_t supNormAccS = sollya_lib_parse_string("1b-10"); // This is the size of the returned interval... 10^-3 should be enough for anybody
			if(DEBUG <= verbose) {
				sollya_lib_printf(">   supnorm(%b, %b, %b, %b, %b);\n", 
													polynomialS, fS, rangeS, absoluteS, supNormAccS);
			}
	    sollya_obj_t supNormRangeS = sollya_lib_supnorm(polynomialS, fS, rangeS, absoluteS, supNormAccS);
			if(sollya_lib_obj_is_error(supNormRangeS)) {
				cout <<  ">   Sollya infnorm failed, but do not loose all hope yet: launching dirtyinfnorm:" << endl;
				sollya_obj_t pminusfS = sollya_lib_sub(polynomialS, fS);
				if(DEBUG <= verbose) {
					sollya_lib_printf(">   dirtyinfnorm(%b, %b);\n", 
														pminusfS, rangeS);
				}
				supNormS = sollya_lib_dirtyinfnorm(pminusfS, rangeS);
				sollya_lib_clear_obj(pminusfS);
				if(sollya_lib_obj_is_error(supNormS)) {
					ostringstream o; 
					o << " ERROR in " << uniqueName_ << " (" << srcFileName << "): " << "Sollya can't seem to be able to compute the infinite norm" << endl; 
					throw o.str();
				}
			}
			else{ // supnorm succeeded, we are mostly interested in the sup of the interval 
				supNormS = sollya_lib_sup(supNormRangeS);
			}
			sollya_lib_clear_obj(supNormAccS);
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
					LSB-=1;
					constLSB-=1;
					tryReducingLSB=false;
					REPORT(DEBUG, "  ... pushing LSB to " << LSB << " and starting over");
				}
				else { // OK, we tried increasing LSB once and it didn't work. Maybe we should increase degree?
					if (degreeSup>degree){
						// restore LSB
						LSB+=1;
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
 
		// Please leave the memory in the state you would like to find it when entering
	  sollya_lib_clear_obj(fixedS);
	  sollya_lib_clear_obj(absoluteS);
	  sollya_lib_clear_obj(targetAccuracyS);
		sollya_lib_clear_obj(degreeIntervalS);
	  sollya_lib_clear_obj(degreeS);
	  sollya_lib_clear_obj(degreeSupS);
	  sollya_lib_clear_obj(rangeS);

		buildFixFormatVector();
	}




	void PiecewisePolyApprox::buildFixFormatVector()
	{
		// compute the MSBs
		for (int i=0; i<=degree; i++){
			sollya_obj_t iS = sollya_lib_constant_from_int(i);
			//sollya_lib_printf("i = %d  = %b\n", i, iS);
			sollya_obj_t coeffS = sollya_lib_coeff(polynomialS, iS);
			//sollya_lib_printf(">  c%d = %b \n", i, coeffS);
			sollya_lib_clear_obj(iS);
			
			int msb,lsb;
			mpfr_t mpcoeff, mptmp;
			// First a tentative conversion to double to get an estimate of the MSB and zeroness
			double dcoeff;
			sollya_lib_get_constant_as_double(&dcoeff, coeffS);
			if(0.0==dcoeff) {
				msb=1; // for the sign
				lsb=0;
				mpfr_init2(mpcoeff, 32);
				mpfr_set_d(mpcoeff, 0.0, GMP_RNDN);
			}
			else{
				msb = floor(log2(fabs(dcoeff)));  // ex. 2.01
				msb++; // For the sign
				lsb = ((i==0? constLSB :  LSB));
				// now we may safely allocate the proper size for the mpfr_t. Add two bits for sign + rounding.
				mpfr_init2(mpcoeff, msb-lsb+3);  
				mpfr_init2(mptmp, msb-lsb+3);  
				sollya_lib_get_constant(mpcoeff, coeffS);
				// Now recompute the MSB explicitely. 
				mpfr_abs(mptmp, mpcoeff, GMP_RNDN); // exact
				mpfr_log2(mptmp, mptmp, GMP_RNDU);
				mpfr_floor(mptmp, mptmp);
				msb = mpfr_get_si(mptmp, GMP_RNDU);
				mpfr_clear(mptmp);	
				msb++;
			}
			FixConstant* fixcoeff =	new FixConstant(msb, lsb, true/*signed*/, mpcoeff);
			coeff.push_back(fixcoeff);

			sollya_lib_clear_obj(coeffS);
			mpfr_clear(mpcoeff);
	}		
		ostringstream debugstring;
		for (int i=0; i<=degree; i++){
			debugstring << endl << "coeff " << i << ": (" << coeff[i]->getMSB() << ", " << coeff[i]->getLSB() << ")   " << coeff[i]->getBitVector();
		}
		REPORT(DEBUG, debugstring.str());
	}

	int PiecewisePolyApprox::getDegree() const	{
		return degree;
	}

	sollya_obj_t PiecewisePolyApprox::getSollyaPolynomial() const	{
		return polynomialS;
	}
} //namespace
