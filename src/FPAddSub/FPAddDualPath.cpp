/*
  Floating Point Adder for FloPoCo

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Authors:   Bogdan Pasca, Florent de Dinechin

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

  All rights reserved
  */

// TODO rework the pipeline properly using the newer framework
// TODO move close path prenormalization up to the Swap Difference box
//   if it becomes a part of the critical path
// TODO remove pipeline stage after finalRoundAdd if slack allows

// TODO Single path adder

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include <utils.hpp>

#include "FPAddDualPath.hpp"

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0


	FPAddDualPath::FPAddDualPath(Target* target, int wE, int wF, bool sub) :
		Operator(target), wE(wE), wF(wF),  sub(sub){

		ostringstream name, synch, synch2;

		srcFileName="FPAddDualPath";

		if(sub)
			name<<"FPSub_";
		else
			name<<"FPAdd_";
		name<<wE<<"_"<<wF;
		setNameWithFreq(name.str());

		setCopyrightString("Bogdan Pasca, Florent de Dinechin (2008)");

		sizeRightShift = intlog2(wF+3);

		/* Set up the IO signals */
		/* Inputs: 2b(Exception) + 1b(Sign) + wE bits (Exponent) + wF bits(Fraction) */
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addFPInput ("X", wE, wF);
		addFPInput ("Y", wE, wF);
		addFPOutput("R", wE, wF);

		//=========================================================================|
		//                          Swap/Difference                                |
		// ========================================================================|

		vhdl<<"-- Exponent difference and swap  --"<<endl;
		vhdl<<tab<<declare("inX",wE+wF+3) << " <= X;"<<endl;
		vhdl<<tab<<declare("inY",wE+wF+3) << " <= Y;"<<endl;

		// signal which indicates whether or not the exception bits of X are greater or equal than/to the exception bits of Y
		vhdl<<tab<<declare("exceptionXSuperiorY") << " <= '1' when inX("<<wE+wF+2<<" downto "<<wE+wF+1<<") >= inY("<<wE+wF+2<<" downto "<<wE+wF+1<<") else '0';"<<endl;

		// signal which indicates whether or not the exception bits of X are equal to the exception bits of Y
		vhdl<<tab<<declare("exceptionXEqualY") << " <= '1' when inX("<<wE+wF+2<<" downto "<<wE+wF+1<<") = inY("<<wE+wF+2<<" downto "<<wE+wF+1<<") else '0';"<<endl;

		// make the difference between the exponents of X and Y; expX - expY = expX + not(expY) + 1
		// pad exponents with sign bit
		vhdl<<tab<<declare("signedExponentX",wE+1) << " <= \"0\" & inX("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		vhdl<<tab<<declare("signedExponentY",wE+1) << " <= \"0\" & inY("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		vhdl<<tab<<declare("exponentDifferenceXY",wE+1) << " <= signedExponentX - signedExponentY ;"<<endl;
		vhdl<<tab<<declare("exponentDifferenceYX",wE) << " <= signedExponentY("<<wE-1<<" downto 0) - signedExponentX("<<wE-1<<" downto 0);"<<endl;

		// SWAP when: [excX=excY and expY>expX] or [excY>excX]
		vhdl<<tab<<declare("swap") << " <= (exceptionXEqualY and exponentDifferenceXY("<<wE<<")) or (not(exceptionXSuperiorY));"<<endl;


		string pmY="inY";
		if ( sub ) {
			vhdl << tab << declare("mY",wE+wF+3)   << " <= inY" << range(wE+wF+2,wE+wF+1) << " & not(inY"<<of(wE+wF)<<") & inY" << range(wE+wF-1,0) << ";"<<endl;
			pmY = "mY";
		}

		// depending on the value of swap, assign the corresponding values to the newX and newY signals
		vhdl<<tab<<declare("newX",wE+wF+3) << " <= " << pmY << " when swap = '1' else inX;"<<endl;
		vhdl<<tab<<declare("newY",wE+wF+3) << " <= inX when swap = '1' else " << pmY << ";"<<endl;
		vhdl<<tab<<declare("exponentDifference",wE) << " <= " << "exponentDifferenceYX"
			 << " when swap = '1' else exponentDifferenceXY("<<wE-1<<" downto 0);"<<endl;

		setCriticalPath(target->adderDelay(wE+1) +  target->lutDelay() + target->lutDelay());

		// determine if the fractional part of Y was shifted out of the operation //
		vhdl<<tab<<declare("shiftedOut") << " <= ";
		if (wE>sizeRightShift){
			manageCriticalPath(target->adderDelay(wE-sizeRightShift));
			for (int i=wE-1;i>=sizeRightShift;i--)
				if (i==sizeRightShift)
					vhdl<< "exponentDifference("<<i<<")";
				else
					vhdl<< "exponentDifference("<<i<<") or ";
			vhdl<<";"<<endl;
		}
		else
			vhdl<<tab<<"'0';"<<endl;

		//shiftVal=the number of positions that fracY must be shifted to the right
		vhdl<<tab<<declare("shiftVal",sizeRightShift) << " <= " ;
		if (wE>sizeRightShift) {
			vhdl << "exponentDifference("<< sizeRightShift-1<<" downto 0)"
				  << " when shiftedOut='0'"<<endl
				  <<tab << tab << "    else CONV_STD_LOGIC_VECTOR("<<wF+3<<","<<sizeRightShift<<") ;" << endl;
		}
		else if (wE==sizeRightShift) {
			vhdl<<tab<<"exponentDifference;" << endl ;
		}
		else 	{ //  wE< sizeRightShift
			vhdl<<tab<<"CONV_STD_LOGIC_VECTOR(0,"<<sizeRightShift-wE <<") & exponentDifference;" <<	endl;
		}

		// was nextCycle();////////////////////////////////////////////////////////////////////////////////////

		// compute EffSub as (signA xor signB) at cycle 1
		manageCriticalPath(2 * target->lutDelay() + 2*target-> localWireDelay());
		vhdl<<tab<<declare("EffSub") << " <= newX("<<wE+wF<<") xor newY("<<wE+wF<<");"<<endl;

		// compute the close/far path selection signal at cycle1
		// the close path is considered only when (signA!=signB) and |exponentDifference|<=1
		vhdl<<tab<<declare("selectClosePath") << " <= EffSub when exponentDifference("<<wE-1<<" downto "<<1<<") = ("<<wE-1<<" downto "<<1<<" => '0') else '0';"<<endl;


		// sdExnXY is a concatenation of the exception bits of X and Y, after swap, so exnX > exnY
		vhdl<<tab<<declare("sdExnXY",4) << " <= newX("<<wE+wF+2<<" downto "<<wE+wF+1<<") "
			 << "& newY("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
		vhdl<<tab<<declare("pipeSignY") << " <= newY("<<wE+wF<<");"<<endl;

		double cp_at_end_of_init = getCriticalPath();

		//=========================================================================|
		//                            close path                                   |
		//=========================================================================|


		vhdl<< endl << "-- Close Path --" << endl;

		// build the fraction signals
		// padding: [sign bit][inplicit "1"][fracX][guard bit]
		vhdl<<tab<<declare("fracXClose1",wF+3) << " <= \"01\" & newX("<<wF-1<<" downto "<<0<<") & '0';"<<endl;

		// the close path is considered when the |exponentDifference|<=1, so
		// the alignment of fracY is of at most 1 position
		vhdl<<tab<<"with exponentDifference(0) select"<<endl;
		vhdl<<tab<<declare("fracYClose1",wF+3) << " <=  \"01\" & newY("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
		vhdl<<tab<<"               \"001\" & newY("<<wF-1<<" downto "<<0<<")       when others;"<<endl;

		// substract the fraction signals for the close path;

		// instanciate the box that computes X-Y and Y-X. Note that it could take its inputs before the swap (TODO ?)
		REPORT(DETAILED, "Building close path dual mantissa subtraction box");
		dualSubClose = new 	IntDualSub(target, wF + 3, 0);
		dualSubClose->changeName(getName()+"_DualSubClose");
		oplist.push_back(dualSubClose);

		inPortMap  (dualSubClose, "X", "fracXClose1");
		inPortMap  (dualSubClose, "Y", "fracYClose1");
		outPortMap (dualSubClose, "RxMy","fracRClosexMy");
		outPortMap (dualSubClose, "RyMx","fracRCloseyMx");
		vhdl << instance(dualSubClose, "DualSubO");

		syncCycleFromSignal("fracRCloseyMx", false);/////////////////////////////////////////////////////////////////

		// register the output -- TODO merge the mux with the last stage of the adder in case of sufficient slack
		nextCycle();////////////////////////////////////////////////////////////////////////////////////

		vhdl<<tab<< declare("fracSignClose") << " <= fracRClosexMy("<<wF+2<<");"<<endl;
		vhdl<<tab<< declare("fracRClose1",wF+2) << " <= fracRClosexMy("<<wF+1<<" downto 0) when fracSignClose='0' else fracRCloseyMx("<<wF+1<<" downto 0);"<<endl;

		if (wF>40) //SP does not need this level
			nextCycle();////////////////////////////////////////////////////////////////////////////////////

		//TODO check the test if significand is all zero is useful.
		vhdl<< tab << declare("resSign") << " <= '0' when selectClosePath='1' and fracRClose1 = ("<<wF+1<<" downto 0 => '0') else"<<endl;
		// else sign(x) xor (close and sign(resclose))
		vhdl<< tab << "          newX("<<wE+wF<<") xor (selectClosePath and "
			 << "fracSignClose);"<<endl;

		// LZC + Shifting. The number of leading zeros are returned together with the shifted input
		REPORT(DETAILED, "Building close path LZC + shifter");
		lzocs = new LZOCShifterSticky(target, wF+2, wF+2, intlog2(wF+2), false, 0);

		lzocs->changeName(getName()+"_LZCShifter");
		oplist.push_back(lzocs);

		inPortMap  (lzocs, "I", "fracRClose1");
		outPortMap (lzocs, "Count","nZerosNew");
		outPortMap (lzocs, "O","shiftedFrac");
		vhdl << instance(lzocs, "LZC_component");

		syncCycleFromSignal("shiftedFrac");/////////////////////////////////////////////////////////////////
		// register the output
		nextCycle();////////////////////////////////////////////////////////////////////////////////////

		// NORMALIZATION

		// shiftedFrac(0) is the round bit, shiftedFrac(1) is the parity bit,
		// shiftedFrac(wF) is the leading one, to be discarded
		// the rounding bit is computed:
		vhdl<<tab<< declare("roundClose0") << " <= shiftedFrac(0) and shiftedFrac(1);"<<endl;
		// Is the result zero?
		vhdl<<tab<< declare("resultCloseIsZero0") << " <= '1' when nZerosNew"
				<< " = CONV_STD_LOGIC_VECTOR(" << (1<< lzocs->getCountWidth())-1 // Should be wF+2 but this is a bug of LZOCShifterSticky: for all zeroes it returns this value
				<< ", " << lzocs->getCountWidth()
			 << ") else '0';" << endl;

		// add two bits in order to absorb exceptions:
		// the second 0 will become a 1 in case of overflow,
		// the first 0 will become a 1 in case of underflow (negative biased exponent)
		vhdl<<tab<< declare("exponentResultClose",wE+2) << " <= (\"00\" & "
			 << "newX("<<wE+wF-1<<" downto "<<wF<<")) "
			 <<"- (CONV_STD_LOGIC_VECTOR(0,"<<wE-lzocs->getCountWidth()+2<<") & nZerosNew);"
			 <<endl;


		nextCycle();////////////////////////////////////////////////////////////////////////////////////
		// concatenate exponent with fractional part before rounding so the possible carry propagation automatically increments the exponent
		vhdl<<tab<<declare("resultBeforeRoundClose",wE+1 + wF+1) << " <= exponentResultClose("<<wE+1<<" downto 0) & shiftedFrac("<<wF<<" downto 1);"<<endl;
		vhdl<<tab<< declare("roundClose") << " <= roundClose0;"<<endl;
		vhdl<<tab<< declare("resultCloseIsZero") << " <= resultCloseIsZero0;"<<endl;





		//=========================================================================|
		//                              far path                                   |
		//=========================================================================|


		vhdl<< endl << "-- Far Path --" << endl;
		// get back to first cycle after exp diff/swap
		setCycleFromSignal("EffSub", true);/////////////////////////////////////////////////////////////////
		setCriticalPath(cp_at_end_of_init);

		//add implicit 1 for frac1.
		vhdl<<tab<< declare("fracNewY",wF+1) << " <= '1' & newY("<<wF-1<<" downto 0);"<<endl;

		// shift right the significand of new Y with as many positions as the exponent difference suggests (alignment)
		REPORT(DETAILED, "Building far path right shifter");
		rightShifter = new Shifter(target,wF+1,wF+3, Shifter::Right);
		rightShifter->changeName(getName()+"_RightShifter");
		oplist.push_back(rightShifter);
		inPortMap  (rightShifter, "X", "fracNewY");
		inPortMap  (rightShifter, "S", "shiftVal");
		outPortMap (rightShifter, "R","shiftedFracY");
		vhdl << instance(rightShifter, "RightShifterComponent");

		syncCycleFromSignal("shiftedFracY", false);/////////////////////////////////////////////////////////////////
		// register the output
		nextCycle();////////////////////////////////////////////////////////////////////////////////////

		// compute sticky bit as the or of the shifted out bits during the alignment //
		vhdl<<tab<< declare("sticky") << " <= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;

		// one cycle only for the sticky (the far path has more time anyway)
		nextCycle();////////////////////////////////////////////////////////////////////////////////////

		//pad fraction of Y [sign][shifted frac having inplicit 1][guard bits]
		vhdl<<tab<< declare("fracYfar", wF+4) << " <= \"0\" & shiftedFracY("<<2*wF+3<<" downto "<<wF+1<<");"<<endl;

		// depending on the signs of the operands, perform addition or substraction
		// the result will be: a + (b xor operation) + operation, where operation=0=addition and operation=1=substraction
		// the operation selector is the xor between the signs of the operands
		// perform xor
		vhdl<<tab<< declare("EffSubVector", wF+4) << " <= ("<<wF+3<<" downto 0 => EffSub);"<<endl;
		vhdl<<tab<<declare("fracYfarXorOp", wF+4) << " <= fracYfar xor EffSubVector;"<<endl;
		//pad fraction of X [sign][inplicit 1][fracX][guard bits]
		vhdl<<tab<< declare("fracXfar", wF+4) << " <= \"01\" & (newX("<<wF-1<<" downto 0)) & \"00\";"<<endl;
		vhdl<<tab<< declare("cInAddFar") << " <= EffSub and not sticky;"<< endl;

		// perform carry in addition
		REPORT(DETAILED, "Building far path adder");
		fracAddFar = new IntAdder(target,wF+4);
		fracAddFar->changeName(getName()+"_fracAddFar");
		oplist.push_back(fracAddFar);
		inPortMap  (fracAddFar, "X", "fracXfar");
		inPortMap  (fracAddFar, "Y", "fracYfarXorOp");
		inPortMap  (fracAddFar, "Cin", "cInAddFar");
		outPortMap (fracAddFar, "R","fracResultfar0");
		vhdl << instance(fracAddFar, "fracAdderFar");

		syncCycleFromSignal("fracResultfar0");/////////////////////////////////////////////////////////////////
		// register the output
		nextCycle();////////////////////////////////////////////////////////////////////////////////////

		vhdl<< tab << "-- 2-bit normalisation" <<endl;
		vhdl<< tab << declare("fracResultFarNormStage", wF+4) << " <= fracResultfar0;"<<endl;

		// NORMALIZATION
		// The leading one may be at position wF+3, wF+2 or wF+1
		//
		vhdl<<tab<< declare("fracLeadingBits", 2) << " <= fracResultFarNormStage("<<wF+3<<" downto "<<wF+2<<") ;" << endl;

		vhdl<<tab<< declare("fracResultFar1",wF) << " <=" << endl ;
		vhdl<<tab<<tab<< "     fracResultFarNormStage("<<wF+0<<" downto 1)  when fracLeadingBits = \"00\" "<<endl;
		vhdl<<tab<<tab<< "else fracResultFarNormStage("<<wF+1<<" downto 2)  when fracLeadingBits = \"01\" "<<endl;
		vhdl<<tab<<tab<< "else fracResultFarNormStage("<<wF+2<<" downto 3);"<<endl;

		vhdl<<tab<< declare("fracResultRoundBit") << " <=" << endl ;
		vhdl<<tab<<tab<< "     fracResultFarNormStage(0) 	 when fracLeadingBits = \"00\" "<<endl;
		vhdl<<tab<<tab<< "else fracResultFarNormStage(1)    when fracLeadingBits = \"01\" "<<endl;
		vhdl<<tab<<tab<< "else fracResultFarNormStage(2) ;"<<endl;

		vhdl<<tab<< declare("fracResultStickyBit") << " <=" << endl ;
		vhdl<<tab<<tab<< "     sticky 	 when fracLeadingBits = \"00\" "<<endl;
		vhdl<<tab<<tab<< "else fracResultFarNormStage(0) or  sticky   when fracLeadingBits = \"01\" "<<endl;
		vhdl<<tab<<tab<< "else fracResultFarNormStage(1) or fracResultFarNormStage(0) or sticky;"<<endl;
		// round bit
		vhdl<<tab<< declare("roundFar1") <<" <= fracResultRoundBit and (fracResultStickyBit or fracResultFar1(0));"<<endl;

		//select operation mode. This depends on wether or not the exponent must be adjusted after normalization
		vhdl<<tab<<declare("expOperationSel",2) << " <= \"11\" when fracLeadingBits = \"00\" -- add -1 to exponent"<<endl;
		vhdl<<tab<<"            else   \"00\" when fracLeadingBits = \"01\" -- add 0 "<<endl;
		vhdl<<tab<<"            else   \"01\";                              -- add 1"<<endl;

		//the second operand depends on the operation selector
		vhdl<<tab<<declare("exponentUpdate",wE+2) << " <= ("<<wE+1<<" downto 1 => expOperationSel(1)) & expOperationSel(0);"<<endl;

		// the result exponent before normalization and rounding is = to the exponent of the first operand //
		vhdl<<tab<<declare("exponentResultfar0",wE+2) << "<=\"00\" & (newX("<<wF+wE-1<<" downto "<<wF<<"));"<<endl;

		vhdl<<tab<<declare("exponentResultFar1",wE+2) << " <= exponentResultfar0 + exponentUpdate;" << endl;

		nextCycle();////////////////////////////////////////////////////////////////////////////////////

		// End of normalization stage
		vhdl<<tab<<declare("resultBeforeRoundFar",wE+1 + wF+1) << " <= "
			 << "exponentResultFar1 & fracResultFar1;" << endl;
		vhdl<<tab<< declare("roundFar") << " <= roundFar1;" << endl;





		//=========================================================================|
		//                              Synchronization                            |
		//=========================================================================|
		vhdl<<endl<<"-- Synchronization of both paths --"<<endl;

		setCycleFromSignal("resultBeforeRoundFar");/////////////////////////////////////////////////////////////////
		syncCycleFromSignal("resultBeforeRoundClose");

		//synchronize the close signal
		vhdl<<tab<<declare("syncClose") << " <= selectClosePath;"<<endl;

		// select between the results of the close or far path as the result of the operation
		vhdl<<tab<< "with syncClose select"<<endl;
		vhdl<<tab<< declare("resultBeforeRound",wE+1 + wF+1) << " <= resultBeforeRoundClose when '1',"<<endl;
		vhdl<<tab<< "                     resultBeforeRoundFar   when others;"<<endl;
		vhdl<<tab<< "with syncClose select"<<endl;
		vhdl<<tab<< declare("round") << " <= roundClose when '1',"<<endl;
		vhdl<<tab<< "         roundFar   when others;"<<endl;

		vhdl<<tab<< declare("zeroFromClose") << " <= syncClose and resultCloseIsZero;" <<endl;
		double muxDelay= target->lutDelay() + target->localWireDelay() + target->ffDelay(); // estimated delay so far (one mux)
				setCriticalPath(muxDelay);
		map<string, double> fraInputDelays;
		fraInputDelays["X"] = muxDelay;
		fraInputDelays["Y"] = 0;
		fraInputDelays["Cin"] = muxDelay;

		vhdl<< endl << "-- Rounding --" << endl;

		REPORT(DETAILED, "Building final round adder");
		// finalRoundAdd will add the mantissa concatenated with exponent, two bits reserved for possible under/overflow
		finalRoundAdd = new IntAdder(target, wE + wF + 2, fraInputDelays);
		finalRoundAdd->changeName(getName()+"_finalRoundAdd");
		oplist.push_back(finalRoundAdd);

		ostringstream zero;
		zero<<"("<<1+wE+wF<<" downto 0 => '0') ";
		inPortMap   (finalRoundAdd, "X", "resultBeforeRound");
		inPortMapCst(finalRoundAdd, "Y", zero.str() );
		inPortMap   (finalRoundAdd, "Cin", "round");
		outPortMap  (finalRoundAdd, "R","resultRounded");
		vhdl << instance(finalRoundAdd, "finalRoundAdder");

		if(finalRoundAdd->getPipelineDepth() == 0) {
			manageCriticalPath(finalRoundAdd->getOutputDelay("R")); // may insert a nextCycle
		}
		else{
			setCriticalPath(finalRoundAdd->getOutputDelay("R"));
		}
		syncCycleFromSignal("resultRounded", false);

		// We neglect the delay of the rest
		vhdl<<tab<<declare("syncEffSub") << " <= EffSub;"<<endl;

		//X
		vhdl<<tab<<declare("syncX",3+wE+wF) << " <= newX;"<<endl;

		//signY
		vhdl<<tab<<declare("syncSignY") << " <= pipeSignY;"<<endl;

		// resSign comes from closer
		vhdl<<tab<<declare("syncResSign") << " <= resSign;"<<endl;

		// compute the exception bits of the result considering the possible underflow and overflow
		vhdl<<tab<< declare("UnderflowOverflow",2) << " <= resultRounded"<<range( wE+1+wF, wE+wF)<<";"<<endl;

		vhdl<<tab<< "with UnderflowOverflow select"<<endl;
		vhdl<<tab<< declare("resultNoExn",wE+wF+3) << "("<<wE+wF+2<<" downto "<<wE+wF+1<<") <=   (not zeroFromClose) & \"0\" when \"01\", -- overflow"<<endl;
		vhdl<<tab<< "                              \"00\" when \"10\" | \"11\",  -- underflow"<<endl;
		vhdl<<tab<< "                              \"0\" &  not zeroFromClose  when others; -- normal "<<endl;

		vhdl<<tab<< "resultNoExn("<<wE+wF<<" downto 0) <= syncResSign & resultRounded("<<wE+wF-1<<" downto 0);"<<endl;

		vhdl<<tab<< declare("syncExnXY", 4) << " <= sdExnXY;"<<endl;
		vhdl<<tab<< "-- Exception bits of the result" << endl;
		vhdl<<tab<< "with syncExnXY select -- remember that ExnX > ExnY "<<endl;
		vhdl<<tab<<tab<< declare("exnR",2) <<" <= resultNoExn("<<wE+wF+2<<" downto "<<wE+wF+1<<") when \"0101\","<<endl;
		vhdl<<tab<<tab<< "        \"1\" & syncEffSub          when \"1010\","<<endl;
		vhdl<<tab<<tab<< "        \"11\"                      when \"1110\","<<endl;
		vhdl<<tab<<tab<< "        syncExnXY(3 downto 2)     when others;"<<endl;
		vhdl<<tab<< "-- Sign bit of the result" << endl;
		vhdl<<tab<< "with syncExnXY select"<<endl;
		vhdl<<tab<<tab<<declare("sgnR") << " <= resultNoExn("<<wE+wF<<")         when \"0101\","<<endl;
		vhdl<<tab<< "           syncX("<<wE+wF<<") and syncSignY when \"0000\","<<endl;
		vhdl<<tab<< "           syncX("<<wE+wF<<")               when others;"<<endl;

		vhdl<<tab<< "-- Exponent and significand of the result" << endl;
		vhdl<<tab<< "with syncExnXY select  "<<endl;
		vhdl<<tab<<tab<< declare("expsigR", wE+wF) << " <= resultNoExn("<<wE+wF-1<<" downto 0)   when \"0101\" ,"<<endl;
		vhdl<<tab<<tab<< "           syncX("<<wE+wF-1<<" downto  0)        when others; -- 0100, or at least one NaN or one infty "<<endl;

		// assign result
		vhdl<<tab<< "R <= exnR & sgnR & expsigR;"<<endl;

	}

	FPAddDualPath::~FPAddDualPath() {
	}










	void FPAddDualPath::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		/* Compute correct value */
		FPNumber fpx(wE, wF, svX);
		FPNumber fpy(wE, wF, svY);
		mpfr_t x, y, r;
		mpfr_init2(x, 1+wF);
		mpfr_init2(y, 1+wF);
		mpfr_init2(r, 1+wF);
		fpx.getMPFR(x);
		fpy.getMPFR(y);
		if(sub)
			mpfr_sub(r, x, y, GMP_RNDN);
		else
			mpfr_add(r, x, y, GMP_RNDN);

		// Set outputs
		FPNumber  fpr(wE, wF, r);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		// clean up
		mpfr_clears(x, y, r, NULL);
	}





	void FPAddDualPath::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		// Regression tests
		tc = new TestCase(this);
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", -1.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);
		tc->addFPInput("X", 2.0);
		tc->addFPInput("Y", -2.0);
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



	TestCase* FPAddDualPath::buildRandomTestCase(int i){

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

	OperatorPtr FPAddDualPath::parseArguments(Target *target, vector<string> &args) {
		int wE;
		UserInterface::parseStrictlyPositiveInt(args, "wE", &wE);
		int wF;
		UserInterface::parseStrictlyPositiveInt(args, "wF", &wF);
		return new FPAddDualPath(target, wE, wF);
	}

	void FPAddDualPath::registerFactory(){
		UserInterface::add("FPAddDualPath", // name
											 "Floating-point adder with dual-path architecture. Trades a larger circuit size for a smaller latency.",
											 "BasicFloatingPoint", // categories
											 "",
											 "wE(int): exponent size in bits; \
wF(int): mantissa size in bits;",
											 "",
											 FPAddDualPath::parseArguments
											 ) ;

	}
}
