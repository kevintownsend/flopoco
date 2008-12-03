/*
 * Floating Point Adder for FloPoCo
 *
 * Author : Bogdan Pasca, Florent de Dinechin
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


// TODO only one normalization
// TODO inverter for close path significand sub is on the critical path
// TODO ABSSUB
// TODO remove pipeline stage after finalRoundAdd if slack allows
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

#include "FPAdder.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPAdder::FPAdder(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR) :
	Operator(target), wEX(wEX), wFX(wFX), wEY(wEY), wFY(wFY), wER(wER), wFR(wFR) {

	int i, j;
	ostringstream name, synch, synch2;

	setOperatorName();
	setOperatorType();
		
	//parameter set up
	wF = wFX;
	wE = wEX;
				
	sizeRightShift = int ( ceil( log2(wF+4)));	
	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	addFPInput ("X", wEX, wFX);
	addFPInput ("Y", wEY, wFY);
	addFPOutput("R", wER, wFR);

		
	//instantiate a right shifter
	rightShifter = new Shifter(target,wFX+1,wFX+3,Right);
	oplist.push_back(rightShifter);
	
	if (isSequential())
	{
		map<string, double> inputs;
		inputs["X"]=0;
		inputs["Y"]=1.5e-10;
		inputs["Cin"]=0;
		fracSubClose = new IntAdder(target, wF + 3, inputs);
		fracSubClose->Operator::setOperatorName("fracSubClose");
		oplist.push_back(fracSubClose);

		complementAdderClose = new IntAdder(target, wF + 2);
		complementAdderClose->Operator::setOperatorName("complementAdderClose");
		oplist.push_back(complementAdderClose);

		// finalRoundAdd will add the mantissa concatenated with exponent, there is one bit reserved for possible overflow 
		finalRoundAdd = new IntAdder(target, wE + wF + 2); 
 		finalRoundAdd->Operator::setOperatorName("finalRoundAdd");
		oplist.push_back(finalRoundAdd);

		fracAddFar = new IntAdder(target,wF+4);
		fracAddFar->Operator::setOperatorName("fracAddFar");
		oplist.push_back(fracAddFar);
				
		lzocs = new LZOCShifterSticky(target, wFX+2, wFX+2,0, 0);
		if(lzocs->getCountWidth() > wE){
			cout << endl << "****WARNING****: wE < log2(wF), the generated VHDL is probably broken."<<endl;
			cout << "    Try increasing wE."<<endl;
		}
		oplist.push_back(lzocs);
	}
	
	if (isSequential()){
		
		//swap/Difference pipeline depth
		swapDifferencePipelineDepth = 1;
		
		closePathDepth = 
			fracSubClose->getPipelineDepth()  
			+ 1  // sub output is registered
			+ 1  // xor stage
			+ complementAdderClose->getPipelineDepth()
			+ 1  // output registered 
			+ lzocs->getPipelineDepth()
			+ 1  //output registered
			+ 1   // exponent update
			;

		farPathDepth = 
			rightShifter->getPipelineDepth() 
			+ 1 //  XOR stage 
			+ fracAddFar->getPipelineDepth()
			+ 1 // output of adder is registered
			+ 1 // normalization stage
			;

		if(verbose){
			cout<<"depths="<<fracAddFar->getPipelineDepth()<<" "<<" "<<" "<<rightShifter->getPipelineDepth()<<endl;
		}
		
		cout<<endl<<"Close path depth = "<< closePathDepth;
		cout<<endl<<"Far path depth   = "<< farPathDepth;
											 
		maxPathDepth = max(closePathDepth, farPathDepth);							 
		
		cout<<endl<<"Max path depth   = "<< maxPathDepth<<endl;
		
		setPipelineDepth(swapDifferencePipelineDepth + maxPathDepth + finalRoundAdd->getPipelineDepth() + 1);
	}
	
	
	
	if (isSequential()){
		
		addSignal("signedExponentX",wE+1);
		addSignal("signedExponentY",wE+1);
		//		addSignal("invSignedExponentY",wE+1);
				
		addSignal("exceptionXSuperiorY");
		addSignal("exceptionXEqualY");
		
		addSignal("exponentDifferenceXY",wE+1);
		addSignal("exponentDifferenceYX",wE);
		
		addSignal("swap",1);
		
		addDelaySignalNoReset("exponentDifference0",wE,1); 
			
		addDelaySignalNoReset("newX",wE+wF+3, 1);
		addDelaySignalNoReset("newY",wE+wF+3, 1);
							
		addSignal("sdY",wE+wF+3);
		//		addSignal("sdClose",1);
		addSignal("sdExponentDifference",wE);
		
		//close path		
		
		addSignal("fracXClose1",wF+3);
		addSignal("fracYClose1",wF+3);
		addSignal("InvFracYClose1",wFX+3);
		addDelaySignalNoReset("fracRClose0",wF+3,complementAdderClose->getPipelineDepth()+2);
		
		
		addDelaySignalNoReset("fracSignClose",1, 1 + complementAdderClose->getPipelineDepth() + 1);
		addDelaySignalNoReset("fracRClose1xor", wF+2,1);
		addDelaySignalNoReset("fracRClose1",wFX+2, lzocs->getPipelineDepth());


		addDelaySignalNoReset("resSign",1,lzocs->getPipelineDepth()+1 + finalRoundAdd->getPipelineDepth()+2);

		addDelaySignalNoReset("nZerosNew",lzocs->getCountWidth(),lzocs->getPipelineDepth());
		addDelaySignalNoReset("shiftedFrac",wFX+2,2);
				
		addDelaySignalNoReset("exponentResultClose",wEX+2, 1);
		// The following magical lines will delay the signal if needed to align it with that of the close path
		addDelaySignalNoReset("roundClose0", 1, 1);
		addDelaySignalNoReset("roundClose", 1, farPathDepth - closePathDepth);
		addDelaySignalNoReset("resultBeforeRoundClose",wE+1 + wF+1, farPathDepth - closePathDepth);


		addDelaySignalNoReset("exponentConcatFrac1",wE+wF+2,1);
		
		
		addSignal("syncResSign",1);
		
		
		//Far path
		
		addDelaySignalNoReset("shiftedOut",1, rightShifter->getPipelineDepth() + fracAddFar->getPipelineDepth() +2);
		
		addSignal("shiftVal",sizeRightShift);
		
		addSignal("fracNewY",wF+1);
		addSignal("shiftedFracY", 2*wF + 4);
		
		addDelaySignalNoReset("sticky",1,  1 + fracAddFar->getPipelineDepth() +1); 
		// + 1 for the xor stage, + 1 because output is registered
		
		addSignal("fracXfar3",wF+4);
		addSignal("fracYfar3",wF+4);
		
		addDelaySignalNoReset("cInAddFar",1,1);

		addDelaySignalNoReset("fracYfar3XorOp",wF+4,1);
		addDelaySignalNoReset("fracResultfar0",wF+4,1);
		
		
		addSignal("fracResultFarNormStage",wF+4);
		addSignal("fracLeadingBits",2);

		addSignal("exponentResultfar0",wE+2);
		addSignal("expOperationSel",2);
		addSignal("exponentUpdate",wE+2);
		addDelaySignalNoReset("exponentResultFar1",wE+2,1);	
						
		addDelaySignalNoReset("fracResultFar1",wF, 1);
		addSignal("fracResultRoundBit",1);
		addSignal("fracResultStickyBit",1);
		// The following magical line will delay the signal if needed to align it with that of the close path
		addDelaySignalNoReset("resultBeforeRoundFar",wE+1 + wF+1, closePathDepth-farPathDepth);
		addDelaySignalNoReset("roundFar",1,closePathDepth-farPathDepth);

		addDelaySignalNoReset("roundFar1",1,1);
		
		int delaySDToRound = 
			maxPathDepth 
			+ finalRoundAdd->getPipelineDepth() 
			+ 1 ; // finalRoundAdd output is registered -- TODO save it if slack allows

		addDelaySignalNoReset("pipeClose",1, delaySDToRound);
		addSignal("syncClose",1);
		addSignal("resultBeforeRound",wE+1 + wF+1);				
		addSignal("UnderflowOverflow",2);
		
		addDelaySignalNoReset("sdExnXY",4, delaySDToRound);
		addSignal("syncExnXY",4);
		
		addDelaySignalNoReset("EffSub",1, delaySDToRound);
		addSignal("syncEffSub");

		addDelaySignalNoReset("pipeX",3+wE+wF, delaySDToRound);
		addSignal("syncX",3+wE+wF);
		addSignal("round",1);

		addDelaySignalNoReset("resultRounded",wE+wF+2, 1);
				
		addDelaySignalNoReset("pipeSignY",1, delaySDToRound);
		addSignal("syncSignY");
		
		addDelaySignalNoReset("resultNoExn",wE+wF+3,1);
		addSignal("finalResult",wE+wF+3);
		addSignal("exponentResultn",wE+2);
		addSignal("fractionResultn",wF+1);
					
		//		cout<<"signals pass";	
	}
	else{
		// addSignal("exceptionXSuperiorY",1);
		// addSignal("exceptionXEqualY",1);
		// addSignal("exponentDifference0",wER+1);
		// addSignal("swap",1);
		// addSignal("exponentDifference1",wER); 
		// addSignal("exponentDifference",wER);  
		// addSignal("zeroExtendedSwap",wER);
		// addSignal("newX",wEX+wFX+3);
		// addSignal("newY",wEX+wFX+3);
		// addSignal("EffSub",1);
		// addSignal("close",1);		
		// addSignal("fracXClose1",wFX+3);
		// addSignal("fracYClose1",wFX+3);
		// addSignal("fracRClose0",wFX+3);
		// addSignal("fracRClose1",wFX+2);
		// addSignal("nZeros",wOutLZC);
		// addSignal("exponentResultClose1",wEX+2);
		// addSignal("exponentResultClose",wEX+2);
		// addSignal("shiftedFrac",2*(wFX+2));
		// addSignal("resultBeforeRoundClose",wE+wF+2);
		// addSignal("exponentConcatFrac1",wE+wF+2);
		// addSignal("shiftedFracY", 2*wF + 4);
		// addSignal("fracNewY",wF+1);
		// addSignal("fracXfar3",wF+5);
		// addSignal("fracYfar3",wF+5);
		// addSignal("fracResultfar0",wF+5);
		// addSignal("fracResultFarNormStage",wF+5);
		// addSignal("exponentResultfar0",wE+1);
		// addSignal("fracResultFar1",wF+3);
		// addSignal("exponentResultFar1",wE+1);
		// addSignal("roundFar",1);
		// addSignal("fractionResultfar0",wF+2);
		// addSignal("rs",1);
		// addSignal("exponentResultfar",wE+1);
		// addSignal("fractionResultfar",wF+1);
		// addSignal("exponentResultn",wE+2);
		// addSignal("fractionResultn",wF+1);
		// addSignal("eMax",1);
		// addSignal("eMin",1);
		// addSignal("UnderflowOverflow",2);
		// addSignal("resultNoExn",wE+wF+3);
		// addSignal("xAB",4);
		// addSignal("finalResult",wE+wF+3);
		// addSignal("fractionResultCloseC",1+wF);
		// addSignal("exponentResultCloseC",wE+2);
		// addSignal("shiftedOut",1);
		// addSignal("sticky",1);
		// addSignal("borrow",1);
		// addSignal("roundClose",1);
	}
}

FPAdder::~FPAdder() {
}

void FPAdder::setOperatorName(){
	ostringstream name;
	/* The name has the format: FPAdder_wEX_wFX_wEY_wFY_wER_wFR where: 
	   wEX = width of X exponenet and 
	   wFX = width for the fractional part of X */
	name<<"FPAdder_"<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR; 
	uniqueName_ = name.str(); 
}

void FPAdder::outputVHDL(std::ostream& o, std::string name) {
  
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;

	int bias_val=int(pow(double(2),double(wEX-1)))-1;
	int i, j; 

	licence(o,"Bogdan Pasca (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);	
	
	//output VHDL components
	if (isSequential())
		lzocs->outputVHDLComponent(o);
	else{
		leadingZeroCounter->outputVHDLComponent(o);
		leftShifter->outputVHDLComponent(o);
	}	
	rightShifter->outputVHDLComponent(o);
	
	fracSubClose->outputVHDLComponent(o);
	complementAdderClose->outputVHDLComponent(o);
	finalRoundAdd->outputVHDLComponent(o);
	fracAddFar->outputVHDLComponent(o);
	// intaddFar2->outputVHDLComponent(o);
	// intaddFar3->outputVHDLComponent(o);
		
	outputVHDLSignalDeclarations(o);	  
	beginArchitecture(o);
		
	if (isSequential()){
	
		//=========================================================================| 
		//                                                                         |
		//                             SEQUENTIAL                                  |
		//                                                                         |
		//=========================================================================|

		outputVHDLRegisters(o); o<<endl;
		//=========================================================================|
		//                          Swap/Difference                                |
		// ========================================================================|
		o<<"-- Exponent difference and swap  --"<<endl;
		// signal which indicates whether or not the exception bits of X are greater or equal than/to the exception bits of Y		  
		o<<tab<<"exceptionXSuperiorY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") >= Y("<<wEY+wFY+2<<" downto "<<wEY+wF+1<<") else '0';"<<endl;
		
		// signal which indicates whether or not the exception bits of X are equal to the exception bits of Y		  
		o<<tab<<"exceptionXEqualY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") = Y("<<wEY+wFY+2<<" downto "<<wEY+wFY+1<<") else '0';"<<endl;

		// make the difference between the exponents of X and Y; expX - expY = expX + not(expY) + 1
		// pad exponents with sign bit
		o<<tab<<"signedExponentX <= \"0\" & X("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl;
		o<<tab<<"signedExponentY <= \"0\" & Y("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl;
		o<<tab<<"exponentDifferenceXY <= signedExponentX - SignedExponentY ;"<<endl;
		o<<tab<<"ExponentDifferenceYX <= SignedExponentY("<<wE-1<<" downto 0) - SignedExponentX("<<wE-1<<" downto 0);"<<endl;

		// SWAP when: [excX=excY and expY>expX] or [excY>excX]
		o<<tab<<"swap <= (exceptionXEqualY and exponentDifferenceXY("<<wE<<")) or (not(exceptionXSuperiorY));"<<endl;

		// depending on the value of swap, assign the corresponding values to the newX and newY signals 
		o<<tab<<"newX <= Y when swap = '1' else X;"<<endl;
		o<<tab<<"newY <= X when swap = '1' else Y;"<<endl;
		o<<tab<<"exponentDifference0 <= exponentDifferenceYX when swap = '1' else exponentDifferenceXY("<<wE-1<<" downto 0);"<<endl;


		// compute EffSub as (signA xor signB) at cycle 1
		o<<tab<<"EffSub <= newX_d("<<wEX+wFX<<") xor newY_d("<<wEY+wFY<<");"<<endl;
		
		// compute the close/far path selection signal at cycle1 
		// the close path is considered only when (signA!=signB) and |exponentDifference|<=1 
		o<<tab<<"pipeClose  <= EffSub when exponentDifference0_d("<<wER-1<<" downto "<<1<<") = ("<<wER-1<<" downto "<<1<<" => '0') else '0';"<<endl;
		
		// assign interface signals
		o<<tab<<"pipeX <= newX_d;"<<endl;
		o<<tab<<"sdY <= newY_d;"<<endl;

		o<<tab<<"sdExponentDifference <= exponentDifference0_d;"<<endl;
		// sdExnXY is a concatenation of the exception bits of X and Y
		o<<tab<<"sdExnXY <= newX_d("<<wE+wF+2<<" downto "<<wE+wF+1<<") & newY_d("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
		o<<tab<<"pipeSignY <= newY_d("<<wE+wF<<");"<<endl;

		// swapDifferencePipelineDepth = 1; ==for now, this section has a constant depth of 1.         

		//=========================================================================|
		//                            close path                                   |
		//=========================================================================|
		
		o << endl << "-- Close Path --" << endl;
		
		// build the fraction signals
		// padding: [sign bit][inplicit "1"][fracX][guard bit]
		o<<tab<<"fracXClose1 <= \"01\" & pipeX("<<wFX-1<<" downto "<<0<<") & '0';"<<endl;
		// the close path is considered when the |exponentDifference|<=1, so 
		// the alignment of fracY is of at most 1 position
		o<<tab<<"with sdExponentDifference(0) select"<<endl;
		o<<tab<<"fracYClose1 <=  \"01\" & sdY("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
		o<<tab<<"               \"001\" & sdY("<<wF-1<<" downto "<<0<<")       when others;"<<endl;
		
		// substract the fraction signals for the close path; 
		// substraction fX - fY = fX + not(fY)+1
		// perform inversion
		o<<tab<<"InvFracYClose1 <= not(fracYClose1);"<<endl;
		// perform addition with carry in = 1
		o<<tab<< "FracSubInClosePath: " << fracSubClose->getOperatorName() << "  -- pipelineDepth="<< fracSubClose->getPipelineDepth() << endl;
		o<<tab<< "  port map ( X => fracXClose1 , " << endl; 
		o<<tab<< "             Y => InvFracYClose1, " << endl; 
		o<<tab<< "             Cin => '1' ," << endl;
		o<<tab<< "             R => fracRClose0, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;
		
		// if substraction result is negative - do a two's compelement of the result, otherwise leave unchanged
		// this is performed as (r xor sign(r)) + r
		// get the sign of the substraction result
		o<<tab<<"fracSignClose <= fracRClose0_d("<<wF+2<<");"<<endl;
		// perform xor
		o<<tab<<"fracRClose1xor <= fracRClose0_d("<<wF+1<<" downto 0) xor ("<<wF+1<<" downto 0 => fracSignClose);"<<endl;
		// perform carry in addition 
		o<<tab<< "ComplementAdderClosePath: " << complementAdderClose->getOperatorName() 
		 << "  -- pipelineDepth="<< complementAdderClose->getPipelineDepth()<< endl;
		o<<tab<< "  port map ( X => fracRClose1xor_d , " << endl; 
		o<<tab<< "             Y => ("<<wF+1<<" downto 0 => '0'), " << endl; 
		o<<tab<< "             Cin => fracSignClose_d ," << endl;
		o<<tab<< "             R => fracRClose1, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;
		
		
		int delayFromSD = fracSubClose->getPipelineDepth() + 1 
			+ complementAdderClose->getPipelineDepth() + 1 
			+ 1;
		string closeDelayed1 = getDelaySignalName("pipeClose", delayFromSD);
		string XDelayed1 = getDelaySignalName("pipeX", delayFromSD);
				
		//TODO check the test if significand is all zero is useful. Should not. Bug in test bench probably.
		o << tab << "resSign <= '0' when "<<closeDelayed1<<"='1' and fracRClose1_d = ("<<wF+1<<" downto 0 => '0') else"<<endl;
		// else sign(x) xor (close and sign(resclose))
		o << tab << "          " << XDelayed1 << "("<<wE+wF<<") xor (" << closeDelayed1 << " and " 
		  << getDelaySignalName("fracSignClose", 1 + complementAdderClose->getPipelineDepth()+1) << ");"<<endl;
		
			
		// LZC + Shifting. The number of leading zeros are returned together with the shifted input
		o<<tab<< "LZCOCS_component: " << lzocs->getOperatorName() 
		 << "  -- pipelineDepth="<< lzocs->getPipelineDepth()<< endl;
		o<<tab<< "      port map ( I => fracRClose1_d, " << endl; 
		o<<tab<< "                 Count => nZerosNew, "<<endl;
		o<<tab<< "                 O => shiftedFrac, " <<endl; 
		o<<tab<< "                 clk => clk, " << endl;
		o<<tab<< "                 rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;		
	
		// NORMALIZATION
		
		// shiftedFrac(0) is the round bit, shiftedFrac(1) is the parity bit, 
		// shiftedFrac(wF) is the leading one, to be discarded
		// the rounding bit is computed:
		o<<tab<< "roundClose0 <= shiftedFrac_d(0) and shiftedFrac_d(1);"<<endl;
		// add two bits in order to absorb exceptions: 
		// the second 0 will become a 1 in case of overflow, the first 0 will become a 1 in case of underflow (negative biased exponent)
		o<<tab<< "exponentResultClose <= (\"00\" & "
		 << getDelaySignalName("pipeX",closePathDepth-1) <<"("<<wE+wF-1<<" downto "<<wF<<")) "
		 <<"- (CONV_STD_LOGIC_VECTOR(0,"<<wE-lzocs->getCountWidth()+2<<") & nZerosNew_d);"
		 <<endl;
		// concatenate exponent with fractional part before rounding so the possible carry propagation automatically increments the exponent 
		o<<tab<<"resultBeforeRoundClose <= exponentResultClose_d("<<wE+1<<" downto 0) & shiftedFrac_d_d("<<wF<<" downto 1);"<<endl; 
		o<<tab<< "roundClose <= roundClose0_d;"<<endl;




		
		//=========================================================================|
		//                              far path                                   |
		//=========================================================================|
		o << endl << "-- Far Path --" << endl;
		
		//input interface signals
		
		// determine if the fractional part of Y was shifted out of the operation //
		if (wE-1>sizeRightShift){
			o<<tab<<"shiftedOut <= "; 
			for (int i=wE-1;i>=sizeRightShift;i--)
				if (((wE-1)==sizeRightShift)||(i==sizeRightShift))
					o<< "sdExponentDifference("<<i<<")";
				else
					o<< "sdExponentDifference("<<i<<") or ";
			o<<";"<<endl;
		}
		else
			o<<tab<<"shiftedOut <= '0';"<<endl; 
		
		//add inplicit 1 for frac1. 
 		o<<tab<< "fracNewY <= '1' & sdY("<<wF-1<<" downto 0);"<<endl;
		
		//selectioSignal=the number of positions that fracY must be shifted to the right				
		if (wE-1>=sizeRightShift)
			o<<tab<<"shiftVal <= sdExponentDifference("<< sizeRightShift-1<<" downto 0"<<"); " << endl; 
		else
			o<<tab<<"shiftVal <= CONV_STD_LOGIC_VECTOR(0,"<<sizeRightShift-(wE-1)<<") & sdExponentDifference("<< wE-1<<" downto 0"<<"); " <<
				endl; 			
								
		// shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) //		
		o<<tab<< "right_shifter_component: " << rightShifter->getOperatorName()
 		 << "  -- pipelineDepth="<< rightShifter->getPipelineDepth()<< endl;
		o<<tab<< "      port map ( X => fracNewY, " << endl; 
		o<<tab<< "                 S => shiftVal, " << endl; 
		o<<tab<< "                 R => shiftedFracY, " <<endl; 
		o<<tab<< "                 clk => clk, " << endl;
		o<<tab<< "                 rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;		
			
		// compute sticky bit as the or of the shifted out bits during the alignment //
		o<<tab<< "sticky<= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;
		
		//pad fraction of Y [sign][shifted frac having inplicit 1][guard bits]
		o<<tab<< "fracYfar3 <= \"0\" & shiftedFracY("<<2*wF+3<<" downto "<<wF+1<<");"<<endl;	
		
		// depending on the signs of the operands, perform addition or substraction			
		// the result will be: a + (b xor operation) + operation, where operation=0=addition and operation=1=substraction
		// the operation selector is the xor between the signs of the operands
		// perform xor 
		o<<tab<<"fracYfar3XorOp <= fracYfar3 xor ("<<wF+3<<" downto 0 => "<<getDelaySignalName("EffSub",rightShifter->getPipelineDepth())<<");"<<endl;
		//pad fraction of X [sign][inplicit 1][fracX][guard bits]				
		o<<tab<< "fracXfar3 <= \"01\" & ("<<getDelaySignalName("pipeX",rightShifter->getPipelineDepth()+1)<<"("<<wF-1<<" downto 0"<<")) & \"00\";"<<endl;
		o<<tab<< "cInAddFar <= " << getDelaySignalName("EffSub",rightShifter->getPipelineDepth())<<" and not sticky;"<< endl;
		// perform carry in addition
		o<<tab<< "fracAddFar0: " << fracAddFar->getOperatorName()
		 << "  -- pipelineDepth="<< fracAddFar->getPipelineDepth()<< endl;
		o<<tab<< "  port map ( X => fracXfar3, " << endl; 
		o<<tab<< "             Y => fracYfar3XorOp_d, " << endl; 
		// below, + 1 because XOR stage
		o<<tab<< "             Cin => cInAddFar_d," << endl;
		o<<tab<< "             R => fracResultfar0, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "            );" << endl<<endl;
	
		// result fracResultfar0 of this adder is registered
		int delayFromSD2 = rightShifter->getPipelineDepth() + 1 + fracAddFar->getPipelineDepth() + 1;
		// if the second operand was shifted out of the operation, then the result of the operation becomes = to the first operand //		
		o << tab << "-- normalisation stage" <<endl; 
		o << tab << "fracResultFarNormStage <= " << getDelaySignalName("fracResultfar0",1) << " when "<< getDelaySignalName("shiftedOut",delayFromSD2)<<"='0' else "
		  <<"(\"01\" & ("<<getDelaySignalName("pipeX",delayFromSD2)<<"("<<wF-1<<" downto 0"<<"))&\"00\");"<<endl;
 		
		// NORMALIZATION 
		// The leading one may be at position wF+3, wF+2 or wF+1
		// 
		o<<tab<< "fracLeadingBits <= fracResultFarNormStage("<<wF+3<<" downto "<<wF+2<<") ;" << endl;

		o<<tab<< "fracResultFar1 <=" << endl ;
		o<<tab<<tab<< "     fracResultFarNormStage("<<wF+0<<" downto 1)  when fracLeadingBits = \"00\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage("<<wF+1<<" downto 2)  when fracLeadingBits = \"01\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage("<<wF+2<<" downto 3);"<<endl;

		o<<tab<< "fracResultRoundBit <=" << endl ;
		o<<tab<<tab<< "     fracResultFarNormStage(0) 	 when fracLeadingBits = \"00\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage(1)    when fracLeadingBits = \"01\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage(2) ;"<<endl;
		
		o<<tab<< "fracResultStickyBit <=" << endl ;
		o<<tab<<tab<< "     "<<getDelaySignalName("sticky", 1 + fracAddFar->getPipelineDepth() +1)<<" 	 when fracLeadingBits = \"00\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage(0) or  "<<getDelaySignalName("sticky", 1 + fracAddFar->getPipelineDepth() +1)<<"   when fracLeadingBits = \"01\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage(1) or fracResultFarNormStage(0) or "<<getDelaySignalName("sticky", 1 + fracAddFar->getPipelineDepth() +1)<<";"<<endl;
		// round bit
		o<<tab<< "roundFar1 <= fracResultRoundBit and (fracResultStickyBit or fracResultFar1(0));"<<endl;

		//select operation mode. This depends on wether or not the exponent must be adjusted after normalization		
		o<<tab<<"expOperationSel <= \"11\" when fracLeadingBits = \"00\" -- add -1 to exponent"<<endl;
		o<<tab<<"            else   \"00\" when fracLeadingBits = \"01\" -- add 0 "<<endl;
		o<<tab<<"            else   \"01\";                              -- add 1"<<endl;
	
		//the second operand depends on the operation selector
		o<<tab<<"exponentUpdate <= ("<<wE+1<<" downto 1 => expOperationSel(1)) & expOperationSel(0);"<<endl;
		
		// the result exponent before normalization and rounding is = to the exponent of the first operand //
		o<<tab<<"exponentResultfar0<=\"00\" & ("<<getDelaySignalName("pipeX", delayFromSD2)<<"("<<wF+wE-1<<" downto "<<wF<<"));"<<endl;
		
		o<<tab<<"exponentResultFar1 <= exponentResultfar0 + exponentUpdate;" << endl;
		
		// End of normalization stage
		o<<tab<<"resultBeforeRoundFar <= "
		 << getDelaySignalName("exponentResultFar1",1) 
		 << " & " <<getDelaySignalName("fracResultFar1",1)<<";" << endl;		
		o<<tab<< "roundFar <= " <<getDelaySignalName("roundFar1",1)<<";" << endl;
		
				
		
		
		
		//=========================================================================|
		//                              Synchronization                            |
		//=========================================================================|				
		o<<"-- Synchronization of both paths --"<<endl;

				
		//synchronize the close signal
		o<<tab<<"syncClose <= "<<getDelaySignalName("pipeClose", maxPathDepth)<<";"<<endl; 
				
		// select between the results of the close or far path as the result of the operation 
		o<<tab<< "with syncClose select"<<endl;
		o<<tab<< "resultBeforeRound <= " << getDelaySignalName("resultBeforeRoundClose",farPathDepth - closePathDepth) <<" when '1',"<<endl;
		o<<tab<< "                     " << getDelaySignalName("resultBeforeRoundFar",  closePathDepth - farPathDepth) <<"   when others;"<<endl;
		o<<tab<< "with syncClose select"<<endl;
		o<<tab<< "round <= " << getDelaySignalName("roundClose",farPathDepth - closePathDepth) <<" when '1',"<<endl;
		o<<tab<< "         " << getDelaySignalName("roundFar", closePathDepth - farPathDepth) <<"   when others;"<<endl;
		
		o << endl << "-- Rounding --" << endl;
		// perform the actual rounding //
		o<<tab<< "FinalRoundAdder: " << finalRoundAdd->getOperatorName() 
		<< "  -- pipelineDepth="<< finalRoundAdd->getPipelineDepth()<< endl;
		o<<tab<< "  port map ( X => resultBeforeRound , " << endl; 
		o<<tab<< "             Y => ("<<1+wE+wF<<" downto 0 => '0'), " << endl; 
		o<<tab<< "             Cin => round ," << endl;
		o<<tab<< "             R => resultRounded, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;


		//sign bit
		int delaySDToRound = 
			maxPathDepth 
			+ finalRoundAdd->getPipelineDepth() 
			+ 1 ; // finalRoundAdd output is registered -- TODO save it if slack allows
		o<<tab<<"syncEffSub <= "<<getDelaySignalName("EffSub", delaySDToRound)<<";"<<endl;
		
		//X
		o<<tab<<"syncX <= "<<getDelaySignalName("pipeX", delaySDToRound)<<";"<<endl;
		
		//signY
		o<<tab<<"syncSignY <= "<<getDelaySignalName("pipeSignY", delaySDToRound)<<";"<<endl;
		
		// resSign comes from closer
		o<<tab<<"syncResSign <= "<<getDelaySignalName("resSign", lzocs->getPipelineDepth()+ 2 + finalRoundAdd->getPipelineDepth()+1)<<";"<<endl;


		o<<tab<< "exponentResultn <= resultRounded_d("<<wE+1+wF<<" downto "<<wF<<"); -- with two leading bits for underflow/underflow"<<endl;
		
		
		// compute the exception bits of the result considering the possible underflow and overflow 
		o<<tab<< "UnderflowOverflow <= exponentResultn("<<wE+1<<" downto "<<wE<<");"<<endl;
		
		o<<tab<< "with UnderflowOverflow select"<<endl;
		o<<tab<< "resultNoExn("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= \"10\" when \"01\", -- overflow"<<endl;
		o<<tab<< "                               \"00\" when \"10\" | \"11\",  -- underflow"<<endl;
		o<<tab<< "                               \"01\" when others; -- normal "<<endl;
  	
		o<<tab<< "resultNoExn("<<wE+wF<<" downto 0) <= syncResSign & resultRounded_d("<<wE+wF-1<<" downto 0);"<<endl;
	 
		o<<tab<< "syncExnXY <= "<<getDelaySignalName("sdExnXY",delaySDToRound)<<";"<<endl;

		o<<tab<< "with syncExnXY select"<<endl;
		o<<tab<< "  finalResult("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= resultNoExn("<<wE+wF+2<<" downto "<<wE+wF+1<<") when \"0101\","<<endl;
		o<<tab<< "                                 \"1\" & syncEffSub              when \"1010\","<<endl;
		o<<tab<< "                                 \"11\"            		           when \"1011\","<<endl;
		o<<tab<< "                                 syncExnXY(3 downto 2)             when others;"<<endl;
		o<<tab<< "with syncExnXY select"<<endl;
		o<<tab<< "  finalResult("<<wE+wF<<") <= resultNoExn("<<wE+wF<<")             when \"0101\","<<endl;
		o<<tab<< "                syncX("<<wE+wF<<") and syncSignY when \"0000\","<<endl;
		o<<tab<< "                syncX("<<wE+wF<<")                 when others;"<<endl;
		
		o<<tab<< "with syncExnXY select"<<endl;
		o<<tab<< "  finalResult("<<wE+wF-1<<" downto 0) <= resultNoExn("<<wE+wF-1<<" downto 0) when \"0101\","<<endl;
		o<<tab<< "                                 syncX("<<wE+wF-1<<" downto 0) when others;"<<endl;
		
		// assign results 

#if DEBUGVHDL
		// probes for debugging
		o<<tab<< "rexp <= finalResult("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		o<<tab<< "rfrac <= \"1\" & finalResult("<<wF-1<<" downto 0);"<<endl;
		o<<tab<< "rsig <= finalResult("<<wE+wF<<");"<<endl;
		o<<tab<< "rexc <= finalResult("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
#endif		
			
		o<<tab<< "R <= finalResult;"<<endl;
		
				   
	} 
	else
	{  
		
	 }
	o<< "end architecture;" << endl << endl;
}




TestIOMap FPAdder::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("X"));
	tim.add(*getSignalByName("Y"));
	tim.add(*getSignalByName("R"));
	return tim;
}

void FPAdder::fillTestCase(mpz_class a[])
{
	/* Get I/O values */
	mpz_class& svX = a[0];
	mpz_class& svY = a[1];
	mpz_class& svR = a[2];

	/* Compute correct value */
	FPNumber x(wEX, wFX), y(wEY, wFY), r(wER, wFR);
	x = svX;
	y = svY;
	r = x+y;
	/* Set outputs */
	svR = r.getSignalValue();
	
}

