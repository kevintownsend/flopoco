/*
  Floating Point Adder for FloPoCo
 
  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Authors:   Bogdan Pasca, Florent de Dinechin

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

#include "FPAdderSinglePath.hpp"

using namespace std;

//TODO +- inf for exponent => update exception

namespace flopoco{

	extern vector<Operator*> oplist;

#define DEBUGVHDL 0


	FPAdderSinglePath::FPAdderSinglePath(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR) :
		Operator(target), wEX(wEX), wFX(wFX), wEY(wEY), wFY(wFY), wER(wER), wFR(wFR) {

		srcFileName="FPAdderSinglePath";
			
		//parameter set up. For now all wEX=wEY=wER and the same holds for fractions
		wF = wFX;
		wE = wEX;
			
		ostringstream name;
		name<<"FPAdderSinglePath_"<<wE<<"_"<<wF<<"_uid"<<getNewUId(); 
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca, Florent de Dinechin (2010)");		

		sizeRightShift = intlog2(wF+3);

		/* Set up the IO signals */
		/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addFPInput ("X", wE, wF);
		addFPInput ("Y", wE, wF);
		addFPOutput("R", wE, wF);

		//=========================================================================|
		//                          Swap/Difference                                |
		// ========================================================================|
		vhdl<<"-- Exponent difference and swap  --"<<endl;
		vhdl<< tab << declare("eXmeY",wE+1) << " <= (\"0\" & X"<<range(wE+wF-1,wF)<<") - (\"0\" & Y"<<range(wE+wF-1,wF)<<");"<<endl;
		vhdl<< tab << declare("eYmeX",wE+1) << " <= (\"0\" & Y"<<range(wE+wF-1,wF)<<") - (\"0\" & X"<<range(wE+wF-1,wF)<<");"<<endl;
				
		vhdl<< tab << declare("swap")       << " <= '1' when Y"<<range(wE+wF-1,0)<<">X"<<range(wE+wF-1,0)<<"else '0';"<<endl;
		vhdl<< tab << declare("eqdiffsign") << " <= '1' when Y"<<range(wE+wF-1,0)<<"=X"<<range(wE+wF-1,0)<<"else '0';"<<endl; 

		// depending on the value of swap, assign the corresponding values to the newX and newY signals 
		vhdl<<tab<<declare("newX",wE+wF+3) << " <= Y     when swap = '1' else X;"<<endl;
		vhdl<<tab<<declare("newY",wE+wF+3) << " <= X     when swap = '1' else Y;"<<endl;
		vhdl<<tab<<declare("expDiff",wE+1) << " <= eXmeY when swap = '0' else eYmeX;"<<endl; 

		//break down the signals
		vhdl << tab << declare("expX",wE) << "<= newX"<<range(wE+wF-1,wF)<<";"<<endl;

		vhdl << tab << declare("excX",2)  << "<= newX"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("excY",2)  << "<= newY"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("signX")   << "<= newX"<<of(wE+wF)<<";"<<endl;
		vhdl << tab << declare("signY")   << "<= newY"<<of(wE+wF)<<";"<<endl;

		vhdl<<tab<<declare("shiftedOut") << " <= '1' when (expDiff >= "<<wF+2<<") else '0';"<<endl;
		//shiftVal=the number of positions that fracY must be shifted to the right				
		vhdl<<tab<<declare("shiftVal",sizeRightShift) << " <= " ;
		if (wE>sizeRightShift) {
			vhdl << "expDiff("<< sizeRightShift-1<<" downto 0)"
				  << " when shiftedOut='0'"<<endl
				  <<tab << tab << "    else CONV_STD_LOGIC_VECTOR("<<wFX+3<<","<<sizeRightShift<<") ;" << endl; 
		}		
		else if (wE==sizeRightShift) {
			vhdl<<tab<<"expDiff;" << endl ;
		}
		else 	{ //  wE< sizeRightShift
			vhdl<<tab<<"CONV_STD_LOGIC_VECTOR(0,"<<sizeRightShift-wE <<") & exponentDifference;" <<	endl;
		}
		
		vhdl<<tab<<declare("EffSub") << " <= signX xor signY;"<<endl;
		vhdl<<tab<<declare("sdsXsYExnXY",6) << " <= signX & signY & excX & excY;"<<endl; 
		vhdl<<tab<<declare("sdExnXY",4) << " <= excX & excY;"<<endl; 


		vhdl<<tab<< declare("fracY",wF+1) << " <= '1' & newY("<<wF-1<<" downto 0);"<<endl;
	
		// shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) 
		REPORT(DETAILED, "Building far path right shifter");	
		rightShifter = new Shifter(target,wFX+1,wFX+3, Shifter::Right);
		rightShifter->changeName(getName()+"_RightShifter");
		oplist.push_back(rightShifter);
		inPortMap  (rightShifter, "X", "fracY");
		inPortMap  (rightShifter, "S", "shiftVal");
		outPortMap (rightShifter, "R","shiftedFracY");
		vhdl << instance(rightShifter, "RightShifterComponent");
		syncCycleFromSignal("shiftedFracY");

		//FIXME: compute inside shifter;
		//compute sticky bit as the or of the shifted out bits during the alignment //
		vhdl<<tab<< declare("sticky") << " <= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;
		
		//pad fraction of Y [overflow][shifted frac having inplicit 1][guard][round]
		vhdl<<tab<< declare("fracYfar", wF+4)      << " <= \"0\" & shiftedFracY("<<2*wF+3<<" downto "<<wF+1<<");"<<endl;	
		vhdl<<tab<< declare("fracYfarXorOp", wF+4) << " <= fracYfar xor ("<<wF+3<<" downto 0 => EffSub);"<<endl;
		//pad fraction of X [overflow][inplicit 1][fracX][guard bits]				
		vhdl<<tab<< declare("fracXfar", wF+4)      << " <= \"01\" & (newX("<<wF-1<<" downto 0)) & \"00\";"<<endl;
		vhdl<<tab<< declare("cInAddFar")           << " <= EffSub and not sticky;"<< endl;//TODO understand why
		
		//result is always positive.
		fracAddFar = new IntAdder(target,wF+4);
		oplist.push_back(fracAddFar);
		inPortMap  (fracAddFar, "X", "fracXfar");
		inPortMap  (fracAddFar, "Y", "fracYfarXorOp");
		inPortMap  (fracAddFar, "Cin", "cInAddFar");
		outPortMap (fracAddFar, "R","fracAddResult");
		vhdl << instance(fracAddFar, "fracAdder");
		syncCycleFromSignal("fracAddResult");

		//shift in place
		vhdl << tab << declare("fracGRS",wF+5) << "<= fracAddResult & sticky; "<<endl;
		lzocs = new LZOCShifterSticky(target, wF+5, wF+5, intlog2(wF+5), false, 0);
		oplist.push_back(lzocs);
		inPortMap  (lzocs, "I", "fracGRS");
		outPortMap (lzocs, "Count","nZerosNew");
		outPortMap (lzocs, "O","shiftedFrac");
		vhdl << instance(lzocs, "LZC_component");
		syncCycleFromSignal("shiftedFrac");

		//need to decide how much to add to the exponent
		vhdl << tab << declare("expPart",wE+2) << " <= (" << zg(wE+2-lzocs->getCountWidth(),0) <<" & nZerosNew) - 1;"<<endl;
		//update exponent
		vhdl << tab << declare("updatedExp",wE+2) << " <= (\"00\" & expX) - expPart;"<<endl;

		vhdl<<tab<<declare("stk")<<"<= shiftedFrac"<<of(1)<<" or shiftedFrac"<<of(0)<<";"<<endl;
		vhdl<<tab<<declare("rnd")<<"<= shiftedFrac"<<of(2)<<";"<<endl;
		vhdl<<tab<<declare("grd")<<"<= shiftedFrac"<<of(3)<<";"<<endl;
		vhdl<<tab<<declare("lsb")<<"<= shiftedFrac"<<of(4)<<";"<<endl;
		
		//decide what to add to the guard bit
		vhdl<<tab<<declare("addToRoundBit")<<"<= '0' when (lsb='0' and grd='1' and rnd='0' and stk='0')  else '1';"<<endl;
		//concatenate exponent with fraction to absorb the possible carry out
		vhdl<<tab<<declare("expFrac",wE+2+wF+1)<<"<= updatedExp & shiftedFrac"<<range(wF+3,3)<<";"<<endl;
		//round
		vhdl<<tab<<declare("RoundedExpFrac",wE+2+wF+1)<<"<= expFrac + addToRoundBit;"<<endl;

		//possible update to exception bits
		vhdl << tab << declare("upExc",2)<<" <= RoundedExpFrac"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		
		vhdl << tab << declare("fracR",wF)<<" <= RoundedExpFrac"<<range(wF,1)<<";"<<endl;
		vhdl << tab << declare("expR",wE) <<" <= RoundedExpFrac"<<range(wF+wE,wF+1)<<";"<<endl;
		
		//exception bits: need to be updated but for not FIXME
		vhdl <<tab<<"with sdsXsYExnXY select "<<endl;
		vhdl <<tab<<declare("excRt",2) << " <= \"00\" when \"000000\"|\"010000\"|\"100000\"|\"110000\","<<endl
		<<tab<<tab<<"\"01\" when \"000101\"|\"010101\"|\"100101\"|\"110101\","<<endl
		<<tab<<tab<<"\"10\" when \"111010\"|\"001010\"|\"001000\"|\"011000\"|\"101000\"|\"111000\"|\"000010\"|\"010010\"|\"100010\"|\"110010\"|\"001001\"|\"011001\"|\"101001\"|\"111001\"|\"000110\"|\"010110\"|\"100110\"|\"110110\", "// and signX=signY) else "<<endl
		<<tab<<tab<<"\"11\" when others;"<<endl;

		
		vhdl << tab << declare("exExpExc",4) << " <= upExc&excRt;"<<endl;
		vhdl << tab << "with (exExpExc) select "<<endl;
		vhdl << tab << declare("excRt2",2) << "<= \"00\" when \"0000\"|\"0100\"|\"1000\"|\"1100\"|\"1001\"|\"1101\","<<endl
		<<tab<<tab<<"\"01\" when \"0001\","<<endl
		<<tab<<tab<<"\"10\" when \"0010\"|\"0110\"|\"0101\","<<endl
		<<tab<<tab<<"\"11\" when others;"<<endl;
		
		vhdl<<tab<<declare("excR",2) << " <= \"00\" when (eqdiffsign='1' and EffSub='1') else excRt2;"<<endl;
		
		vhdl <<tab<<declare("signR") << "<= '0' when (sdsXsYExnXY=\"100000\" or sdsXsYExnXY=\"010000\") else signX;"<<endl;

		// assign result 
		vhdl<<tab<< declare("computedR",wE+wF+3) << " <= excR & signR & expR & fracR;"<<endl;
		vhdl<<tab<<"with sdExnXY select"<<endl;
		vhdl<<tab<<"R <= newX when \"0100\"|\"1000\"|\"1001\", newY when \"0001\"|\"0010\"|\"0110\", computedR when others;"<<endl;
	}

	FPAdderSinglePath::~FPAdderSinglePath() {
	}


	void FPAdderSinglePath::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
	
		/* Compute correct value */
		FPNumber fpx(wEX, wFX), fpy(wEY, wFY);
		fpx = svX;
		fpy = svY;
		mpfr_t x, y, r;
		mpfr_init2(x, 1+wFX);
		mpfr_init2(y, 1+wFY);
		mpfr_init2(r, 1+wFR); 
		fpx.getMPFR(x);
		fpy.getMPFR(y);
		mpfr_add(r, x, y, GMP_RNDN);

		// Set outputs 
		FPNumber  fpr(wER, wFR, r);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		// clean up
		mpfr_clears(x, y, r, NULL);
	}





	void FPAdderSinglePath::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		// Regression tests 
		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", -1.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", FPNumber::plusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", FPNumber::minusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::plusInfty);
		tc->addFPInput("Y", FPNumber::minusInfty);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::plusInfty);
		tc->addFPInput("Y", FPNumber::plusInfty);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::minusInfty);
		tc->addFPInput("Y", FPNumber::minusInfty);
		emulate(tc);
		tcl->add(tc);
	
	}



	TestCase* FPAdderSinglePath::buildRandomTestCase(int i){

		TestCase *tc;
		mpz_class x,y;
		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);
		mpz_class negative  = mpz_class(1)<<(wE+wF);

		tc = new TestCase(this); 
		/* Fill inputs */
		if ((i & 7) == 0) {// cancellation, same exponent
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
		}
		else if ((i & 7) == 1) {// cancellation, exp diff=1
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			e++; // may rarely lead to an overflow, who cares
			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
		}
		else if ((i & 7) == 2) {// cancellation, exp diff=1
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
			e++; // may rarely lead to an overflow, who cares
			y  = getLargeRandom(wF) + (e << wF) + normalExn;
		}
		else if ((i & 7) == 3) {// alignment within the mantissa sizes
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow, who cares
			y  = getLargeRandom(wF) + (e << wF) + normalExn;
		}
		else if ((i & 7) == 4) {// subtraction, alignment within the mantissa sizes
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow
			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
		}
		else if ((i & 7) == 5 || (i & 7) == 6) {// addition, alignment within the mantissa sizes
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow
			y  = getLargeRandom(wF) + (e << wF) + normalExn;
		}
		else{ //fully random
			x = getLargeRandom(wE+wF+3);
			y = getLargeRandom(wE+wF+3);
		}
		// Random swap
		mpz_class swap = getLargeRandom(1);
		if (swap == mpz_class(0)) {
			tc->addInput("X", x);
			tc->addInput("Y", y);
		}
		else {
			tc->addInput("X", y);
			tc->addInput("Y", x);
		}
		/* Get correct outputs */
		emulate(tc);
		return tc;
	}

}
