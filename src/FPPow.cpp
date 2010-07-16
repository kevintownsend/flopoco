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

		addConstant("wE", "positive", wE);
		addConstant("wF", "positive", wF);
		
		

		vhdl << tab << declare("logIn", 3+wF+wE+logG+wE) << " <= X & " << rangeAssign(wE+logG-1, 0, "'0'") << " ;" << endl; 

		FPLog* log = new FPLog(target,  wE,  wF+wE+logG, logTableSize );
		oplist.push_back(log);
		inPortMap(log, "X", "logIn");
		outPortMap(log, "R", "lnX");
		vhdl << instance(log, "log");

		syncCycleFromSignal("lnX");
		nextCycle();

		FPMultiplier* mult = new FPMultiplier(target,   /*X:*/ wE, wF+wE+logG,   /*Y:*/ wE, wF,  /*R: */  wE,  wF+wE+expG,  true /* norm*/);
		oplist.push_back(mult);
		inPortMap(mult, "Y", "Y");
		inPortMap(mult, "X", "lnX");
		outPortMap(mult, "R", "P");
		vhdl << instance(mult, "mult");

		syncCycleFromSignal("P");
		nextCycle();

		FPExp* exp = new FPExp(target,  wE,  wF, 0/* means default*/, 0, expG, true);
		oplist.push_back(exp);
		inPortMap(exp, "X", "P");
		outPortMap(exp, "R", "E");
		vhdl << instance(exp, "exp");

		nextCycle();

		// That was the part that takes 99% of the size. Now the exception 

		vhdl << tab  << declare("flagsE", 2) << " <= E(wE+wF+2 downto wE+wF+1);" << endl;

		vhdl << tab  << declare("flagsX", 2) << " <= X(wE+wF+2 downto wE+wF+1);" << endl;
		vhdl << tab  << declare("signX") << " <= X(wE+wF);" << endl;
		vhdl << tab  << declare("expFieldX", wE) << " <= X(wE+wF-1 downto wF);" << endl;
		vhdl << tab  << declare("fracX", wF) << " <= X(wF-1 downto 0);" << endl;

		vhdl << tab  << declare("flagsY", 2) << " <= Y(wE+wF+2 downto wE+wF+1);" << endl;
		vhdl << tab  << declare("signY") << " <= Y(wE+wF);" << endl;
		vhdl << tab  << declare("expFieldY", wE) << " <= Y(wE+wF-1 downto wF);" << endl;
		vhdl << tab  << declare("fracY", wF) << " <= Y(wF-1 downto 0);" << endl;

		vhdl << tab  << declare("infX") << " <= '1' when flagsX=\"10\" else '0';" << endl;
		vhdl << tab  << declare("infY") << " <= '1' when flagsY=\"10\" else '0';" << endl;
		vhdl << tab  << declare("infE") << " <= '1' when flagsE=\"10\" else '0';" << endl;
		vhdl << tab  << declare("nanX") << " <= '1' when flagsX=\"11\" else '0';" << endl;
		vhdl << tab  << declare("nanY") << " <= '1' when flagsY=\"11\" else '0';" << endl;
		vhdl << tab  << declare("zeroX") << " <= '1' when flagsX=\"00\" else '0';" << endl;
		vhdl << tab  << declare("zeroY") << " <= '1' when flagsY=\"00\" else '0';" << endl;
		vhdl << tab  << declare("zeroE") << " <= '1' when flagsE=\"00\" else '0';" << endl;
		vhdl << tab  << declare("normalX") << " <= '1' when flagsX=\"01\" else '0';" << endl;
		vhdl << tab  << declare("normalY") << " <= '1' when flagsY=\"01\" else '0';" << endl;

		vhdl << tab  << declare("XisOne") << " <= '1' when X = \"0100\" & " << rangeAssign(wE-2, 0, "'1'") << " & " << rangeAssign(wF-1, 0, "'0'") << " else '0';" << endl;

		vhdl << tab  << "-- x^(+/-0)=1" << endl;
		vhdl << tab  << declare("case_x_to_0") << " <= '1' when (signX='0' and normalX='1' and zeroY='1') else '0';" << endl;

		vhdl << tab  << "-- +-0^^negative =+inf" << endl;
		vhdl << tab  << declare("case_0_to_yn") << " <= '1' when (zeroX='1' and normalY='1' and signY='1')  else '0';" << endl;

		vhdl << tab  << "-- +-0^^positive =+0" << endl;
		vhdl << tab  << declare("case_0_to_yp") << " <= '1' when (zeroX='1' and normalY='1' and signY='0')  else '0';" << endl;

		vhdl << tab  << "-- -+-0^^-inf=+inf" << endl;
		vhdl << tab  << declare("case_0_to_neginf") << " <= '1' when (zeroX='1' and infY='1' and signY='1')  else '0';" << endl;

		vhdl << tab  << "--  -0^^+inf=NaN" << endl;
		vhdl << tab  << declare("case_neg0_to_posinf") << " <= '1' when (zeroX='1' and signY='1' and infY='1' and signY='0')  else '0';" << endl;

		vhdl << tab  << "-- 0^^+inf=0" << endl;
		vhdl << tab  << declare("case_pos0_to_posinf") << " <= '1' when (zeroX='1' and signY='0' and infY='1' and signY='0')  else '0';" << endl;

		vhdl << tab  << "--  1^+-inf=NaN" << endl;
		vhdl << tab  << declare("case_pos1_to_inf") << " <= '1' when (XisOne='1'  and signX='0' and infY='1') else '0';" << endl;

		vhdl << tab  << "-- 1^+-fin=1" << endl;
		vhdl << tab  << declare("case_1_to_finite") << " <= '1' when (XisOne='1' and normalY='1') else '0';" << endl;

		vhdl << tab  << "--  +-0^+-0=NaN" << endl;
		vhdl << tab  << declare("case_0_to_0") << " <= zeroX and zeroY;" << endl;

		vhdl << tab  << "--  +inf^+-0=NaN" << endl;
		vhdl << tab  << declare("case_posinf_to_0") << " <=  '1' when (infX='1' and signX='0' and zeroY='1') else '0';" << endl;

		vhdl << tab  << "-- When do we get a NaN" << endl;
		vhdl << tab  << declare("RisComputedNaN") << " <= case_0_to_0  or (signX and normalX) or case_neg0_to_posinf;" << endl;
		vhdl << tab  << declare("RisNaN") << " <= nanX or nanY or RisComputedNaN;" << endl;

		vhdl << tab  << "-- When do we get an inf" << endl;
		vhdl << tab  << declare("RisInf") << " <= case_0_to_yn or case_0_to_neginf or infE;" << endl;

		vhdl << tab  << "-- When do we get a zero" << endl;
		vhdl << tab  << declare("RisZero") << " <= zeroE or (zeroX and normalX) or (infX and signY) or case_0_to_neginf;" << endl;

		vhdl << tab  << declare("flagR", 2) << " <= " << endl
		     << tab  << tab << "     \"11\" when RisNaN='1'" << endl
		     << tab  << tab << "else \"00\" when RisZero='1'" << endl
		     << tab  << tab << "else \"10\" when RisInf='1'" << endl
		     << tab  << tab << "else \"01\";" << endl;

		vhdl << tab << declare("signR") << " <= '0';" << endl;
		vhdl << tab << "R <= flagR & signR & E " << range(wE+wF-1, 0) <<" ;" << endl; 

	}	

	FPPow::~FPPow()
	{
	}







	void FPPow::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
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
		TestCase *tc;


		tc = new TestCase(this); 
		tc->addFPInput("X", 3.0);
		tc->addFPInput("Y", 2.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 9.0);
		tc->addFPInput("Y", 0.5);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 2.0);
		tc->addFPInput("Y", 0.5);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 3423.0);
		tc->addFPInput("Y", 0.125234);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", .23432);
		tc->addFPInput("Y", 0.5342);
		emulate(tc);
		tcl->add(tc);


	}



	// One test out of 8 fully random (tests NaNs etc)
	// All the remaining ones test positive numbers.
	// with special treatment for exponents 0 and -1, 
	// and for the range reduction worst case.
 
	void FPPow::buildRandomTestCases(TestCaseList* tcl, int n){
	}

}
