/*
  Floating Point Adder with 3 inputs for FloPoCo
 
  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Authors:   Bogdan Pasca

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

  */
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "FPAdder3Input.hpp"
#include "IntCompressorTree.hpp"

using namespace std;
namespace flopoco{

	extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPAdder3Input::FPAdder3Input(Target* target, int wE, int wF, map<string, double> inputDelays) :
		Operator(target), wE(wE), wF(wF) {
		srcFileName="FPAdder3Input";
			
		ostringstream name;
		name<<"FPAdder3Input_"<<wE<<"_"<<wF<<"_uid"<<getNewUId(); 
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca (2010)");		

		sizeRightShift = intlog2(2*wF);

		/* Set up the IO signals */
		/* Inputs: 2b(Exception) + 1b(Sign) + wE bits (Exponent) + wF bits(Fraction) */
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addFPInput ("X", wE, wF);
		addFPInput ("Y", wE, wF);
		addFPInput ("Z", wE, wF);
		addFPOutput("R", wE, wF);
		
		//split the inputs
		vhdl << tab << declare("neX",wE) << "<= X"<<range(wE+wF-1,wF)<<";"<<endl;
		vhdl << tab << declare("neY",wE) << "<= Y"<<range(wE+wF-1,wF)<<";"<<endl;
		vhdl << tab << declare("neZ",wE) << "<= Z"<<range(wE+wF-1,wF)<<";"<<endl;

		//split the inputs
		vhdl << tab << declare("mX",wF) << "<= X"<<range(wF-1,0)<<";"<<endl;
		vhdl << tab << declare("mY",wF) << "<= Y"<<range(wF-1,0)<<";"<<endl;
		vhdl << tab << declare("mZ",wF) << "<= Z"<<range(wF-1,0)<<";"<<endl;
		
		vhdl << tab << declare("sX") << "<= X"<<of(wE+wF)<<";"<<endl;
		vhdl << tab << declare("sY") << "<= Y"<<of(wE+wF)<<";"<<endl;
		vhdl << tab << declare("sZ") << "<= Z"<<of(wE+wF)<<";"<<endl;
		
		vhdl << tab << declare("excX",2) << " <= X"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("excY",2) << " <= Y"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("excZ",2) << " <= Z"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		
		//shift right exponents for numbers that are 0
		vhdl << tab << declare("eX",wE) << " <= "<<zg(wE,0)<< " when excX=\"00\" else neX;"<<endl; 
		vhdl << tab << declare("eY",wE) << " <= "<<zg(wE,0)<< " when excY=\"00\" else neY;"<<endl;
		vhdl << tab << declare("eZ",wE) << " <= "<<zg(wE,0)<< " when excZ=\"00\" else neZ;"<<endl;
		
		vhdl << tab << declare("expRes1") << " <= '1' when (excX(1)='1' or excY(1)='1' or excZ(1)='1') else '0';"<<endl;
		vhdl << tab << declare("expRes0") << " <= '1' when (((excX(1)='0' and excY(1)='0' and excZ(1)='0') and (excX(0)='1' or excY(0)='1' or excZ(0)='1') ) or " << endl;
		vhdl << tab << tab << "(excX=\"10\" and sX='1' and excY=\"10\" and sY='0') or "<<endl;
		vhdl << tab << tab << "(excX=\"10\" and sX='0' and excY=\"10\" and sY='1') or "<<endl;
		vhdl << tab << tab << "(excZ=\"10\" and sZ='1' and excY=\"10\" and sY='0') or "<<endl;
		vhdl << tab << tab << "(excZ=\"10\" and sZ='0' and excY=\"10\" and sY='1') or "<<endl;
		vhdl << tab << tab << "(excX=\"10\" and sX='1' and excZ=\"10\" and sZ='0') or "<<endl;
		vhdl << tab << tab << "(excX=\"10\" and sX='0' and excZ=\"10\" and sZ='1') or "<<endl;
		vhdl << tab << tab << "(excX=\"11\" or excY=\"11\" or excZ=\"11\")) else '0';"<<endl;
		
		vhdl << tab << declare("sgn") << "<= '1' when ((excX=\"10\" and sX='1') or (excY=\"10\" and sY='1') or (excZ=\"10\" and sZ='1') ) else '0';"<<endl; //or (sX='1' and sY='1' and sX='1')  
		
		// ---------- determine max exponent -----------------
		//det max exponent difference:
		
		vhdl << tab << declare("dexy",wE+1) << " <= (\"0\" & eX) - (\"0\" & eY);"<<endl;
		vhdl << tab << declare("deyz",wE+1) << " <= (\"0\" & eY) - (\"0\" & eZ);"<<endl;
		vhdl << tab << declare("dezx",wE+1) << " <= (\"0\" & eZ) - (\"0\" & eX);"<<endl;
		
		vhdl << tab << declare("cdexy") << " <= dexy"<<of(wE)<<";"<<endl;
		vhdl << tab << declare("cdeyz") << " <= deyz"<<of(wE)<<";"<<endl;
		vhdl << tab << declare("cdezx") << " <= dezx"<<of(wE)<<";"<<endl;
		
		vhdl << tab << declare("eMaxSel",3) << " <= cdexy & cdeyz & cdezx;"<<endl;
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("eMax",wE) << " <= eX when \"001\"|\"011\","<<endl;
		vhdl << tab << tab << " eY when \"100\"|\"101\","<<endl;
		vhdl << tab << tab << " eZ when \"010\"|\"110\","<<endl;
		vhdl << tab << tab << " eX when others;"<<endl;
		
// 		vhdl << tab << "with (eMaxSel) select " <<endl;
// 		vhdl << tab << declare("semax") << " <=  sX when \"100\"|\"110\","<<endl;
// 		vhdl << tab << tab << " sY when \"010\"|\"011\","<<endl;
// 		vhdl << tab << tab << " sZ when \"001\"|\"101\","<<endl;
// 		vhdl << tab << tab << " sX when others;"<<endl;
		
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("m1") << " <=  '1' when \"010\"|\"110\", '0' when others;"<<endl;
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("m2") << " <= '1' when \"001\"|\"011\", '0' when others;"<<endl;
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("m3") << " <= '1' when \"100\"|\"101\", '0' when others;"<<endl;
		
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("alpha",2) << " <= \"00\" when \"010\"|\"110\","<<endl;
		vhdl << tab << tab << " \"01\" when \"100\"|\"101\","<<endl;
		vhdl << tab << tab << " \"10\" when \"001\"|\"011\","<<endl;
		vhdl << tab << tab << " \"11\" when others;"<<endl;
		
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("beta",2) << " <= \"00\" when \"001\"|\"011\","<<endl;
		vhdl << tab << tab << " \"01\" when \"010\"|\"110\","<<endl;
		vhdl << tab << tab << " \"10\" when \"100\"|\"101\","<<endl;
		vhdl << tab << tab << " \"11\" when others;"<<endl;
		
		vhdl << tab << "with (eMaxSel) select " <<endl;
		vhdl << tab << declare("gama",2) << " <= \"00\" when \"100\"|\"101\","<<endl;
		vhdl << tab << tab << " \"01\" when \"001\"|\"011\","<<endl;
		vhdl << tab << tab << " \"10\" when \"010\"|\"110\","<<endl;
		vhdl << tab << tab << " \"11\" when others;"<<endl;
		
		vhdl << tab << declare("mux1out",wE+1) << " <= dexy when m1='0' else dezx;"<<endl;
		vhdl << tab << declare("mux2out",wE+1) << " <= deyz when m2='0' else dexy;"<<endl;
		vhdl << tab << declare("mux3out",wE+1) << " <= dezx when m3='0' else deyz;"<<endl;
	
		//possible negative fix
		vhdl << tab << declare("nmux1out", wE+1) << " <=  not(mux1out) + '1';" <<endl;
		vhdl << tab << declare("nmux2out", wE+1) << " <=  not(mux2out) + '1';" <<endl;
		vhdl << tab << declare("nmux3out", wE+1) << " <=  not(mux3out) + '1';" <<endl;
		
		//alignment of mantissas 
		vhdl << tab << "with (alpha) select " <<endl;
		vhdl << tab << declare("amx",wE+1) << " <= mux1out when \"00\","<<endl;
		vhdl << tab << tab << "nmux1out when \"01\", "<<endl;
		vhdl << tab << tab << zg(wE+1,0)<<" when \"10\","<<endl;
		vhdl << tab << zg(wE+1,0)<<" when others;"<<endl;
		
		//alignment of mantissas 
		vhdl << tab << "with (beta) select " <<endl;
		vhdl << tab << declare("bmx",wE+1) << " <= mux2out when \"00\","<<endl;
		vhdl << tab << tab << "nmux2out when \"01\", "<<endl;
		vhdl << tab << tab << zg(wE+1,0)<<" when \"10\","<<endl;
		vhdl << tab << zg(wE+1,0)<<" when others;"<<endl;
		
		//alignment of mantissas 
		vhdl << tab << "with (gama) select " <<endl;
		vhdl << tab << declare("gmx",wE+1) << " <= mux3out when \"00\","<<endl;
		vhdl << tab << tab << "nmux3out when \"01\", "<<endl;
		vhdl << tab << tab << zg(wE+1,0)<<" when \"10\","<<endl;
		vhdl << tab << zg(wE+1,0)<<" when others;"<<endl;
		
		
		// changed sign is semax
		// SIGN 
		//sign m2, m3
		
/*		vhdl << tab << declare("sm2")<< " semax xor sY;" << endl;
		vhdl << tab << declare("sm3")<< " semax xor sZ;" << endl;*/
		
		//create the mantissas
		
		
		//shift right exponents for numbers that are 0
		
		vhdl << tab << declare("nfX",wF+1) << " <=  \"1\" & mX;"<<endl;
		vhdl << tab << declare("nfY",wF+1) << " <=  \"1\" & mY;"<<endl;
		vhdl << tab << declare("nfZ",wF+1) << " <=  \"1\" & mZ;"<<endl;
		
		vhdl << tab << declare("fX",wF+1) << " <= "<<zg(wF+1,0)<< " when excX=\"00\" else nfX;"<<endl; 
		vhdl << tab << declare("fY",wF+1) << " <= "<<zg(wF+1,0)<< " when excY=\"00\" else nfY;"<<endl;
		vhdl << tab << declare("fZ",wF+1) << " <= "<<zg(wF+1,0)<< " when excZ=\"00\" else nfZ;"<<endl;
		
		vhdl << tab << declare("sval0", sizeRightShift) << "<= amx"<<range(sizeRightShift-1,0)<<";"<<endl;
		vhdl << tab << declare("sval1", sizeRightShift) << "<= bmx"<<range(sizeRightShift-1,0)<<";"<<endl;
		vhdl << tab << declare("sval2", sizeRightShift) << "<= gmx"<<range(sizeRightShift-1,0)<<";"<<endl;
		
		vhdl << tab << declare("sout0") << " <=  '1' when amx"<<range(wE, sizeRightShift)<<">0 else '0';"<<endl;
		vhdl << tab << declare("sout1") << " <=  '1' when bmx"<<range(wE, sizeRightShift)<<">0 else '0';"<<endl;
		vhdl << tab << declare("sout2") << " <=  '1' when gmx"<<range(wE, sizeRightShift)<<">0 else '0';"<<endl;
		
		Shifter *s0 = new Shifter(target, wF+1, 2*wF, Shifter::Right, inDelayMap("X",getCriticalPath()));
		oplist.push_back(s0);
		
		inPortMap  (s0, "X", "fX");
		inPortMap  (s0, "S", "sval0");
		outPortMap (s0, "R","sfX");
		vhdl << instance(s0, "fxShifter");
		
		Shifter *s1 = new Shifter(target, wF+1, 2*wF, Shifter::Right, inDelayMap("X",getCriticalPath()));
		oplist.push_back(s1);
		
		inPortMap  (s1, "X", "fY");
		inPortMap  (s1, "S", "sval1");
		outPortMap (s1, "R","sfY");
		vhdl << instance(s1, "fyShifter");
		
		Shifter *s2 = new Shifter(target, wF+1, 2*wF, Shifter::Right, inDelayMap("X",getCriticalPath()));
		oplist.push_back(s2);
		inPortMap  (s2, "X", "fZ");
		inPortMap  (s2, "S", "sval2");
		outPortMap (s2, "R","sfZ");
		vhdl << instance(s2, "fzShifter");
		syncCycleFromSignal("sfZ");
		
		//check if these two were not shifted out
		vhdl << tab << declare("efX", 3 + 3*wF+1) << " <= (\"000\" & sfX) when sout0='0' else "<<zg(3 + 3*wF+1,0)<<";"<<endl;
		vhdl << tab << declare("efY", 3 + 3*wF+1) << " <= (\"000\" & sfY) when sout1='0' else "<<zg(3 + 3*wF+1,0)<<";"<<endl;
		vhdl << tab << declare("efZ", 3 + 3*wF+1) << " <= (\"000\" & sfZ) when sout2='0' else "<<zg(3 + 3*wF+1,0)<<";"<<endl;
		
		//sign them properly
		
		vhdl << tab << declare("sefX", 3 + 3*wF+1) << " <= (efX xor "<<rangeAssign( 3 + 3*wF,0,"sX") << ") + sX;"<<endl;
		vhdl << tab << declare("sefY", 3 + 3*wF+1) << " <= (efY xor "<<rangeAssign( 3 + 3*wF,0,"sY") << ") + sY;"<<endl;
		vhdl << tab << declare("sefZ", 3 + 3*wF+1) << " <= (efZ xor "<<rangeAssign( 3 + 3*wF,0,"sZ") << ") + sZ;"<<endl;
		
		//add them together
		REPORT(DEBUG, "HERE");
		IntCompressorTree *ct = new IntCompressorTree(target, 3 + 3*wF+1, 3);
		oplist.push_back(ct);
		
		inPortMap( ct, "X0", "sefX");
		inPortMap( ct, "X1", "sefY");
		inPortMap( ct, "X2", "sefZ");
		outPortMap( ct, "R", "addRes");
		
		vhdl << tab << instance(ct, "Adder3Op") << endl;
		syncCycleFromSignal("addRes");
		
		vhdl << tab << declare("trSign") << " <= addRes"<<of(3 + 3*wF)<< ";"<< endl;
		// get the result back to positive and remove the signals
		
		vhdl << tab << declare("posExtF",3 + 3*wF) << " <= (addRes"<<range(3 + 3*wF-1, 0) << " xor "<<rangeAssign(3 + 3*wF-1, 0, "trSign")<<") + trSign;"<<endl; 		
		
		//sticky
		vhdl << tab << declare("sticky") << " <=  '1' when posExtF"<<range(wF-2,0)<<">"<<zg(wF-1,0)<<" else '0';"<<endl;
		
		vhdl << tab << declare("posExtFsticky", 2*wF + 4) << " <= posExtF"<<range(3 + 3*wF-1, wF) << " & sticky;" << endl;
		
		
		lzocs = new LZOCShifterSticky(target, 2*wF + 4, 2*wF + 4, intlog2(2*wF + 4), false, 0, inDelayMap("I",getCriticalPath()));
		oplist.push_back(lzocs);
		inPortMap  (lzocs, "I", "posExtFsticky");
		outPortMap (lzocs, "Count","nZerosNew");
		outPortMap (lzocs, "O","shiftedFrac");
		vhdl << instance(lzocs, "LZC_component");
		syncCycleFromSignal("shiftedFrac");

		vhdl<<tab<<declare("stk")<<"<= '1' when shiftedFrac"<<range(wF,0)<<">"<<zg(wF+1,0)<<"else '0';"<<endl;
		vhdl<<tab<<declare("rnd")<<"<= shiftedFrac"<<of(wF+1)<<";"<<endl;
		vhdl<<tab<<declare("grd")<<"<= shiftedFrac"<<of(wF+2)<<";"<<endl;
		vhdl<<tab<<declare("lsb")<<"<= shiftedFrac"<<of(wF+3)<<";"<<endl;
		
		vhdl << tab << declare("tfracR", wF+1) << "<= shiftedFrac"<<range(2*wF + 2, wF + 2)<<";"<<endl;
		vhdl<<tab<<declare("addToRoundBit")<<"<= '0' when (lsb='0' and grd='1' and rnd='0' and stk='0')  else '1';"<<endl;		
		
		vhdl << tab << declare("roundedExpFrac", 1 + wE + wF+1) << "<= (\"0\" & eMax & tfracR) + addToRoundBit;" << endl;
		vhdl << tab << declare("tnexp", wE + 1) << "<= roundedExpFrac" << range(wE+wF+1, wF+1)<<";"<<endl;
		//update exponent
		vhdl << tab << declare("upexp",wE+1) << "  <=  tnexp + CONV_STD_LOGIC_VECTOR(2,"<<wE+1<<") - ("<<zg( wE+1 - intlog2(2*wF+4), 0) << " & nZerosNew);"<<endl;
		
		//compose result
		
		vhdl << tab << declare("path1_exp", wE) << " <= upexp"<<range(wE-1,0) << ";" << endl;
		vhdl << tab << declare("path1_frac", wF) << " <= roundedExpFrac"<< range(wF,1)<<";" << endl;
		
		//exception 
		
		vhdl << tab << declare("path1_exc",2) << "<= expRes1 & expRes0;"<<endl;
		
		//sign 
		
		vhdl << tab << declare("path1_sign") << "<= sgn when path1_exc>\"01\" else trSign;" << endl;
		vhdl << tab << declare("path2_sign") << " <= '1' when ((sX='1' and sX='1' and sZ='1') or path1_sign='1') else path1_sign;"<<endl;
		
		//result path 1 for now
		vhdl << tab << "R <= path1_exc & path1_sign & path1_exp & path1_frac;" << endl;
	}

	FPAdder3Input::~FPAdder3Input() {
	}


	void FPAdder3Input::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svZ = tc->getInputValue("Z");
		
		
		/* Compute correct value */
		FPNumber fpx(wE, wF), fpy(wE, wF), fpz(wE, wF);
		fpx = svX;
		fpy = svY;
		fpz = svZ;
		mpfr_t x, y, z, r_tmp, r;
		mpfr_init2(x, 5*wF);
		mpfr_init2(y, 5*wF);
		mpfr_init2(z, 5*wF);
		mpfr_init2(r_tmp, 5*wF); 
		mpfr_init2(r, 5*wF); 
		fpx.getMPFR(x);
		fpy.getMPFR(y);
		fpz.getMPFR(z);
		mpfr_add(r_tmp, x, y, GMP_RNDN);
		mpfr_add(r, r_tmp, z, GMP_RNDN);
		// Set outputs 
		FPNumber  fpr(wE, wF, r);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		// clean up
		mpfr_clears(x, y, r, NULL);
	}





// 	void FPAdder3Input::buildStandardTestCases(TestCaseList* tcl){
// 		TestCase *tc;
// 
// 		// Regression tests 
// 		tc = new TestCase(this); 
// 		tc->addFPInput("X", 1.0);
// 		tc->addFPInput("Y", -1.0);
// 		emulate(tc);
// 		tcl->add(tc);
// 
// 		tc = new TestCase(this); 
// 		tc->addFPInput("X", 1.0);
// 		tc->addFPInput("Y", FPNumber::plusDirtyZero);
// 		emulate(tc);
// 		tcl->add(tc);
// 
// 		tc = new TestCase(this); 
// 		tc->addFPInput("X", 1.0);
// 		tc->addFPInput("Y", FPNumber::minusDirtyZero);
// 		emulate(tc);
// 		tcl->add(tc);
// 
// 		tc = new TestCase(this); 
// 		tc->addFPInput("X", FPNumber::plusInfty);
// 		tc->addFPInput("Y", FPNumber::minusInfty);
// 		emulate(tc);
// 		tcl->add(tc);
// 
// 		tc = new TestCase(this); 
// 		tc->addFPInput("X", FPNumber::plusInfty);
// 		tc->addFPInput("Y", FPNumber::plusInfty);
// 		emulate(tc);
// 		tcl->add(tc);
// 
// 		tc = new TestCase(this); 
// 		tc->addFPInput("X", FPNumber::minusInfty);
// 		tc->addFPInput("Y", FPNumber::minusInfty);
// 		emulate(tc);
// 		tcl->add(tc);
// 	
// 	}



// 	TestCase* FPAdder3Input::buildRandomTestCase(int i){
// 
// 		TestCase *tc;
// 		mpz_class x,y;
// 		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);
// 		mpz_class negative  = mpz_class(1)<<(wE+wF);
// 
// 		tc = new TestCase(this); 
// 		/* Fill inputs */
// 		if ((i & 7) == 0) {// cancellation, same exponent
// 			mpz_class e = getLargeRandom(wE);
// 			x  = getLargeRandom(wF) + (e << wF) + normalExn;
// 			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
// 		}
// 		else if ((i & 7) == 1) {// cancellation, exp diff=1
// 			mpz_class e = getLargeRandom(wE);
// 			x  = getLargeRandom(wF) + (e << wF) + normalExn;
// 			e++; // may rarely lead to an overflow, who cares
// 			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
// 		}
// 		else if ((i & 7) == 2) {// cancellation, exp diff=1
// 			mpz_class e = getLargeRandom(wE);
// 			x  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
// 			e++; // may rarely lead to an overflow, who cares
// 			y  = getLargeRandom(wF) + (e << wF) + normalExn;
// 		}
// 		else if ((i & 7) == 3) {// alignment within the mantissa sizes
// 			mpz_class e = getLargeRandom(wE);
// 			x  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
// 			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow, who cares
// 			y  = getLargeRandom(wF) + (e << wF) + normalExn;
// 		}
// 		else if ((i & 7) == 4) {// subtraction, alignment within the mantissa sizes
// 			mpz_class e = getLargeRandom(wE);
// 			x  = getLargeRandom(wF) + (e << wF) + normalExn;
// 			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow
// 			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
// 		}
// 		else if ((i & 7) == 5 || (i & 7) == 6) {// addition, alignment within the mantissa sizes
// 			mpz_class e = getLargeRandom(wE);
// 			x  = getLargeRandom(wF) + (e << wF) + normalExn;
// 			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow
// 			y  = getLargeRandom(wF) + (e << wF) + normalExn;
// 		}
// 		else{ //fully random
// 			x = getLargeRandom(wE+wF+3);
// 			y = getLargeRandom(wE+wF+3);
// 		}
// 		// Random swap
// 		mpz_class swap = getLargeRandom(1);
// 		if (swap == mpz_class(0)) {
// 			tc->addInput("X", x);
// 			tc->addInput("Y", y);
// 		}
// 		else {
// 			tc->addInput("X", y);
// 			tc->addInput("Y", x);
// 		}
// 		/* Get correct outputs */
// 		emulate(tc);
// 		return tc;
// 	}

}
