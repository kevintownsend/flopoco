/*
 * An FP Power for FloPoCo
 *
 * Author : Florent de Dinechin
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <fstream>
#include <sstream>
#include <math.h>	// for NaN
#include "FPLog.hpp"
#include "FPExp.hpp"
#include "FPPow.hpp"
#include "FPMultiplier.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"


using namespace std;


namespace flopoco{

	extern vector<Operator *> oplist;


	void FPPow::compute_error(mpfr_t & r, mpfr_t &epsE, mpfr_t& epsM, mpfr_t& epsL ) {
		mpfr_t r1, one;
		mpfr_init2(r1, 1000);
		mpfr_init2(one, 16);
		mpfr_set_d(one,  1.0, GMP_RNDN);
		mpfr_add(r1,one, epsE, GMP_RNDN);
		// r1 is (1+epsE)

		mpfr_init2(r, 1000);
		mpfr_set(r, epsL, GMP_RNDN);
		mpfr_mul(r,r, epsM, GMP_RNDN);
		mpfr_add(r,r, epsM, GMP_RNDN);
		mpfr_add(r,r, epsL, GMP_RNDN);
		mpfr_add(r,r, one, GMP_RNDN);
		// here r is 1 + epsM + epsL + epxMexpL

		// multiply it by max ( y ln(x) )

		FPNumber fpMaxFloat = FPNumber(wE, wF, FPNumber::largestPositive);
		mpfr_t maxExpIn;
		mpfr_init2(maxExpIn, 1000);
		fpMaxFloat.getMPFR(maxExpIn);
		mpfr_log(maxExpIn, maxExpIn, GMP_RNDU);
		mpfr_mul(r,r,maxExpIn, GMP_RNDN);
		if(verbose) mpfr_out_str (stderr, 10, 30, maxExpIn, GMP_RNDN); cerr << " ";

		// then take the exp
		mpfr_exp(r,r, GMP_RNDN);

		// and that's it.
		mpfr_mul(r,r,r1, GMP_RNDN);


		mpfr_clears(r1, one, NULL);

	}






	FPPow::FPPow(Target* target, int wE, int wF, int logTableSize, int expTableSize, int expDegree, int expG, int logG)
		: Operator(target), wE(wE), wF(wF)
	{

		setCopyrightString("F. de Dinechin, C. Klein  (2008)");

		ostringstream o;

		o << "FPPow_" << wE << "_" << wF << "_";
		if(target->isPipelined()) 
			o << target->frequencyMHz() ;
		else
			o << "comb";
		setName(o.str());

		addFPInput("X", wE, wF);
		addFPInput("Y", wE, wF);
		addFPOutput("R", wE, wF, 2); // 2 because faithfully rounded

		
		vhdl << tab << declare("logIn", wF+wE+logG+2) << " <= X & " << rangeAssign(wE+logG-1, 0, "'0'") << " ;" << endl; 

#if 0
		FPLog* log = new FPLog(target,  wE,  wF+wE+logG, logTableSize );
		oplist.push_back(log);
		outPortMap(log, "R", "lnX");
		inPortMap(log, "X", "logIn");
		vhdl << instance(log, "log");


		FPMultiplier* mult = FPMultiplier(target,   /*X:*/ wE, wF+wE+logG,   /*Y:*/ wE, wF,  /*R: */  wE,  wF+wE+expG,  0 /*no norm*/);
		oplist.push_back(exp);
		outPortMap(exp, "R", "powBeforeround");
		inPortMap(exp, "X", "expIn");
		vhdl << instance(exp, "exp");


		FPExp* exp = new FPLog(target,  wE,  wF+wE+logG, expTableSize, expDegree);
		oplist.push_back(exp);
		outPortMap(exp, "R", "powBeforeround");
		inPortMap(exp, "X", "expIn");
		vhdl << instance(exp, "exp");

#endif
	}	

	FPPow::~FPPow()
	{
	}







	void FPPow::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("X");
		
		/* Compute correct value */
		FPNumber fpx(wE, wF);
		FPNumber fpy(wE, wF);
		fpx = svX;
		fpy = svY;
		
		mpfr_t x,y, ru,rd;
		mpfr_init2(x,  1+wF);
		mpfr_init2(y,  1+wF);
		mpfr_init2(ru, 1+wF);
		mpfr_init2(rd, 1+wF); 
		fpx.getMPFR(x);
		fpy.getMPFR(y);
		mpfr_pow(rd, x, y, GMP_RNDD);
		mpfr_pow(ru, x, y, GMP_RNDU);
		FPNumber  fprd(wE, wF, rd);
		FPNumber  fpru(wE, wF, ru);
		mpz_class svRD = fprd.getSignalValue();
		mpz_class svRU = fpru.getSignalValue();
		tc->addExpectedOutput("R", svRD);
		tc->addExpectedOutput("R", svRU);
		mpfr_clears(x, y, ru, rd, NULL);
	}
 

	// TEST FUNCTIONS


	void FPPow::buildStandardTestCases(TestCaseList* tcl){

	}



	// One test out of 8 fully random (tests NaNs etc)
	// All the remaining ones test positive numbers.
	// with special treatment for exponents 0 and -1, 
	// and for the range reduction worst case.
 
	void FPPow::buildRandomTestCases(TestCaseList* tcl, int n){
	}

}
