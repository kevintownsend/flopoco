/*

  A class that manages polynomial approximation for FloPoCo (and possibly later for metalibm).

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
	then the Socrate team at INSA de Lyon

  Initial software.
  Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL

  All rights reserved.

*/

#include "BasicPolyApprox.hpp"
#include <sstream>
#include <iomanip>

namespace flopoco{


	BasicPolyApprox::BasicPolyApprox(FixFunction *f_, double targetAccuracy, int addGuardBits):
		f(f_)
	{
		needToFreeF = false;
		initialize();
		buildApproxFromTargetAccuracy(targetAccuracy,  addGuardBits);
		buildFixFormatVector();
	}


	BasicPolyApprox::BasicPolyApprox(string sollyaString_, double targetAccuracy, int addGuardBits, bool signedIn)
	{
		//  parsing delegated to FixFunction
		f = new FixFunction(sollyaString_, signedIn);
		needToFreeF = true;
		initialize();
		buildApproxFromTargetAccuracy(targetAccuracy,  addGuardBits);
		buildFixFormatVector();
	}



	BasicPolyApprox::BasicPolyApprox(sollya_obj_t fS_, double targetAccuracy, int addGuardBits, bool signedIn)
	{
		f = new FixFunction(fS_,signedIn);
		needToFreeF = true;
		initialize();
		buildApproxFromTargetAccuracy(targetAccuracy,  addGuardBits);
		buildFixFormatVector();
	}

	BasicPolyApprox::BasicPolyApprox(sollya_obj_t fS_, int degree_, int lsb_, bool signedIn):
		degree(degree_), LSB(lsb_)
	{
		f = new FixFunction(fS_, signedIn);
		needToFreeF = true;
		initialize();
		buildApproxFromDegreeAndLSBs();
		buildFixFormatVector();
	}



	BasicPolyApprox::BasicPolyApprox(int degree_, vector<int> MSB, int LSB_, vector<mpz_class> mpzCoeff):
	  degree(degree_), LSB(LSB_)
  {
		needToFreeF = false;
		for (int i=0; i<=degree; i++){
			FixConstant* fixcoeff =	new FixConstant(MSB[i], LSB, true/*signed*/, mpzCoeff[i]);
			coeff.push_back(fixcoeff);
		}
	}


	void BasicPolyApprox::initialize() {
		srcFileName="BasicPolyApprox"; // should be somehow static but this is too much to ask me
		fixedS = sollya_lib_fixed();
		absoluteS = sollya_lib_absolute();
	}



	BasicPolyApprox::~BasicPolyApprox()
	{
		// clear the Sollya constants
	  sollya_lib_clear_obj(fixedS);
	  sollya_lib_clear_obj(absoluteS);

		// clear the FixedFunction, only if we created it here
		if(needToFreeF)	free(f);

		// clear other attributes
	  sollya_lib_clear_obj(polynomialS);
		//	  sollya_lib_clear_obj(S);
		if(coeff.size()!=0){
			for (unsigned int i=0; i<coeff.size(); i++)
				free(coeff[i]);
		}
	}





	// This is a static (class) method.
	void BasicPolyApprox::guessDegree(sollya_obj_t fS, sollya_obj_t rangeS, double targetAccuracy, int* degreeInfP, int* degreeSupP) {
		// Accuracy has to be converted to sollya objects
		// a few constant objects
		if(DEBUG <= verbose)
			sollya_lib_printf("> BasicPolyApprox::guessDegree() for function %b on range %b at target accuracy %1.5e\n", fS, rangeS, targetAccuracy);
		sollya_obj_t targetAccuracyS = sollya_lib_constant_from_double(targetAccuracy);

		// initial evaluation of the required degree
		//guessdegree(f,I,eps,w,bound ) : (function, range, constant, function, constant) → range
		sollya_obj_t degreeIntervalS = sollya_lib_guessdegree(fS, rangeS, targetAccuracyS, NULL);
		sollya_obj_t degreeInfS = sollya_lib_inf(degreeIntervalS);
		sollya_obj_t degreeSupS = sollya_lib_sup(degreeIntervalS);
		sollya_lib_get_constant_as_int(degreeInfP, degreeInfS);
		sollya_lib_get_constant_as_int(degreeSupP, degreeSupS);
		if(DEBUG <= verbose)
			sollya_lib_printf("> BasicPolyApprox::guessDegree(): degree of poly approx should be in %b\n", degreeIntervalS);
	  sollya_lib_clear_obj(targetAccuracyS);
		sollya_lib_clear_obj(degreeIntervalS);
	  sollya_lib_clear_obj(degreeInfS);
	  sollya_lib_clear_obj(degreeSupS);
	}




	void BasicPolyApprox::buildApproxFromTargetAccuracy(double targetAccuracy, int addGuardBits)
	{
		// a few constant objects
		sollya_obj_t fS = f->fS; // no need to free this one
		sollya_obj_t rangeS = f->rangeS; // no need to free this one

		// calling the class method guessDegree
		int degreeSup;
		guessDegree(fS, rangeS, targetAccuracy, &degree, &degreeSup);


		// This will be the LSB of the constant (unless extended below)
		LSB = floor(log2(targetAccuracy));
		REPORT(DEBUG, "LSB without guard bits is " << LSB);

		// A few lines to add guard bits to the constant, it will be for free in terms of evaluation
		double coeffAccuracy = targetAccuracy;

		if (-1==addGuardBits) {
			double maxEvalErrorInUlps = degree; // this assumes faithful multipliers in Horner scheme
			coeffAccuracy /= maxEvalErrorInUlps;
		}
		else if (addGuardBits>0) {
			// the caller provided the number of bits to add in addGuardBitss
			for (int i=0; i<addGuardBits; i++)
				coeffAccuracy /= 2;
		}
		LSB = floor(log2(coeffAccuracy));
		REPORT(DEBUG, "Initial LSB with guard bits is " << LSB);

		sollya_obj_t degreeS = sollya_lib_constant_from_int(degree);


		// now launch fpminimax, measure the approx error, and iterate if it fails
		bool success=false;
		bool tryReducingLSB=true;
		while(not success) {

			buildApproxFromDegreeAndLSBs();

			// did we success in getting an accurate enough polynomial?
			if(approxErrorBound < targetAccuracy) {
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
					tryReducingLSB=false;
					REPORT(DEBUG, "  ... pushing LSB to " << LSB << " and starting over");
				}
				else { // OK, we tried pushing LSB once and it didn't work. Maybe we should increase degree?
					if (degreeSup>degree){
						// restore LSB
						LSB+=1;
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
	  sollya_lib_clear_obj(degreeS);
	}




	void BasicPolyApprox::buildApproxFromDegreeAndLSBs()
	{
		sollya_obj_t fS = f->fS; // no need to free this one
		sollya_obj_t rangeS = f->rangeS; // no need to free this one
		sollya_obj_t degreeS = sollya_lib_constant_from_int(degree);

		REPORT(DEBUG, "Trying to build coefficients with LSB=" << LSB);
		// Build the list of coefficient LSBs for fpminimax
		// Sollya library is a bit painful, it is safer just build a big string and parse it.
		ostringstream s;
		s << "[|";
		for(int i=0; i<=degree ; i++) {
			s << -LSB;
			if(i<degree) s<< ",";
		}
		s << "|]";
		sollya_obj_t coeffSizeListS = sollya_lib_parse_string(s.str().c_str());
		if(DEBUG <= verbose) {
			sollya_lib_printf("> BasicPolyApprox::buildApproxFromDegreeAndLSBs:    fpminimax(%b, %b, %b, %b, %b, %b);\n",
												fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS);
		}
		// Tadaaa! After all this we may launch fpminimax
		polynomialS = sollya_lib_fpminimax(fS, degreeS, coeffSizeListS, rangeS, fixedS, absoluteS, NULL);
		sollya_lib_clear_obj(coeffSizeListS);
		if(DEBUG <= verbose)
			sollya_lib_printf("> BasicPolyApprox::buildBasicPolyApprox: obtained polynomial   %b\n", polynomialS);

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

		sollya_lib_get_constant_as_double(& approxErrorBound, supNormS);
		sollya_lib_clear_obj(supNormS);

		REPORT(DEBUG, "Polynomial accuracy is " << approxErrorBound);
		// Please leave the memory in the state you would like to find it when entering
	  sollya_lib_clear_obj(degreeS);
	}







	void BasicPolyApprox::buildFixFormatVector()
	{
		// compute the MSBs
		int msb,lsb;
		for (int i=0; i<=degree; i++){
			sollya_obj_t iS = sollya_lib_constant_from_int(i);
			//sollya_lib_printf("i = %d  = %b\n", i, iS);
			sollya_obj_t coeffS = sollya_lib_coeff(polynomialS, iS);
			//sollya_lib_printf(">  c%d = %b \n", i, coeffS);
			sollya_lib_clear_obj(iS);

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
				lsb = LSB;
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

		// printing debug string out of the final vector
		ostringstream debugstring;
		debugstring << "buildFixFormatVector:";
		int maxmsb=coeff[degree]->MSB;
		int lsb0 = coeff[0]->LSB;
		for (int i=0; i<degree; i++){
			if(coeff[i]->MSB > maxmsb)
				maxmsb = coeff[i]->MSB;
		}
		int bitwidth = maxmsb - lsb0 + 3;
		for (int i=0; i<=degree; i++){
			int lsb= coeff[i]->LSB;
			int msb= coeff[i]->MSB;
			debugstring <<  endl << "> coeff" << right << setw(2) << i << ": "
									<< " (" << setw(2)<< msb << ", " << setw(3)<< lsb << ")   "
									<< setw(bitwidth + lsb0-lsb) << coeff[i]->getBitVector()
									<< "  " << setw(10) << printMPFR(coeff[i]->fpValue) ;
		}
		REPORT(DEBUG, debugstring.str());
	}

} //namespace
