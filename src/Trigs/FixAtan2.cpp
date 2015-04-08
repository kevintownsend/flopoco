/*
   An arctangent(y/x) implementation

	Author:  Florent de Dinechin, Matei Istoan

	This file is part of the FloPoCo project

	Initial software.
	Copyright © INSA Lyon, INRIA, CNRS, UCBL,
	2014.
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

#include "../utils.hpp"
#include "Operator.hpp"
#include "FixAtan2.hpp"

using namespace std;

namespace flopoco {


	// The constructor for a stand-alone operator
	FixAtan2::FixAtan2(Target* target_, int wIn_, int wOut_,  map<string, double> inputDelays_):
		Operator(target_, inputDelays_),
		target(target_), wIn(wIn_), wOut(wOut_)
	{
	}


	FixAtan2::~FixAtan2()
	{

	}


	void FixAtan2::emulate(TestCase * tc)
	{
		mpfr_t x,y,a, constPi;
		mpfr_init2(x, 10*wIn);
		mpfr_init2(y, 10*wIn);
		mpfr_init2(a, 10*wIn);
		mpfr_init2(constPi, 10*wIn);
		mpfr_const_pi( constPi, GMP_RNDN);

		mpz_class az;

		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		// interpret as signed two'ss complement
		if (1==(svX >> (wIn-1))) // sign bit
			svX -= (1<<wIn);
		if (1==(svY >> (wIn-1))) // sign bit
			svY -= (1<<wIn);
		/* Compute correct value */

		mpfr_set_z (x, svX.get_mpz_t(), GMP_RNDN); //  exact
		mpfr_set_z (y, svY.get_mpz_t(), GMP_RNDN); //  exact

		mpfr_atan2(a, y, x, GMP_RNDN); // a between -pi and pi
		mpfr_div(a, a, constPi, GMP_RNDN); // a between -1 and 1

		// Now convert a to fix point
		// Align to fix point by adding 6 -- if we just add 4 there is a 1-bit shift in case a<0
		mpfr_add_d(a, a, 6.0, GMP_RNDN);
		mpfr_mul_2si (a, a, wOut-1, GMP_RNDN); // exact scaling

		mpz_class mask = (mpz_class(1)<<wOut) -1;

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDD); // there can be a real rounding here
		az -= mpz_class(6)<<(wOut-1);
		az &= mask;
		tc->addExpectedOutput ("R", az);

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDU); // there can be a real rounding here
		az -= mpz_class(6)<<(wOut-1);
		az &= mask;
		tc->addExpectedOutput ("R", az);

		// clean up
		mpfr_clears (x,y,a, constPi, NULL);
	}





}
