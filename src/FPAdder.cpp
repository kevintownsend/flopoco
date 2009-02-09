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

// TODO move close path prenormalization up to the Swap Difference box
//   if it becomes a part of the critical path
// TODO remove pipeline stage after finalRoundAdd if slack allows

// TODO clean up X propagation to remove warnings

// TODO Rename FracXFar -> FracXFar, same for Y and other 3   
// TODO remove last pipe stage from the Shifter -- beware fix LongAcc also

// TODO Single path adder

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
	
#ifdef _WIN32
	sizeRightShift = intlog2(wF+4);	
#else
	sizeRightShift = int ( ceil( log2(wF+3)));	
#endif	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	addFPInput ("X", wEX, wFX);
	addFPInput ("Y", wEY, wFY);
	addFPOutput("R", wER, wFR);

		
	// map<string, double> inputs;
	// inputs["X"]=0;
	// inputs["Y"]=1.5e-10;
	// inputs["Cin"]=0;
	

	// -----Sub-components---------------------------------
	
	// Close path
	dualSubClose = new 	IntDualSub(target, wF + 3, 0);
	dualSubClose->Operator::setOperatorName(getOperatorName()+"_DualSubClose");
	oplist.push_back(dualSubClose);
		
	lzocs = new LZOCShifterSticky(target, wFX+2, wFX+2,0, 0);
	if(lzocs->getCountWidth() > wE){
		cout << endl << "****WARNING****: wE < log2(wF), the generated VHDL is probably broken."<<endl;
		cout << "    Try increasing wE."<<endl;
	}
	oplist.push_back(lzocs);

	// Far path
	rightShifter = new Shifter(target,wFX+1,wFX+3,Right);
	oplist.push_back(rightShifter);

 	fracAddFar = new IntAdder(target,wF+4);
	fracAddFar->Operator::setOperatorName(getOperatorName()+"_fracAddFar");
	oplist.push_back(fracAddFar);

	// finalRoundAdd will add the mantissa concatenated with exponent, two bits reserved for possible under/overflow 
	finalRoundAdd = new IntAdder(target, wE + wF + 2); 
	finalRoundAdd->Operator::setOperatorName(getOperatorName()+"_finalRoundAdd");
	oplist.push_back(finalRoundAdd);


	// Pipeline setup-------------------------------------------------

	//swap/Difference pipeline depth
	swapDifferencePipelineDepth = 1;
		
	closePathDepth = 
		+ dualSubClose->getPipelineDepth()
		+ 1  // output registered 
		+ 1  // swap 
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
	
	
	maxPathDepth = max(closePathDepth, farPathDepth);							 

	if(verbose){
		cout<<endl<<"Close path depth = "<< closePathDepth<<endl;
		cout<<"Far path depth   = "<< farPathDepth<<endl;
		//cout<<endl<<"Max path depth   = "<< maxPathDepth<<endl;
		//cout<<"depths="<<fracAddFar->getPipelineDepth()<<" "<<" "<<" "<<rightShifter->getPipelineDepth()<<endl;
	}



	if(isSequential())		
		setPipelineDepth(swapDifferencePipelineDepth + maxPathDepth + finalRoundAdd->getPipelineDepth() + 1);
	else
		setPipelineDepth(0);

	delaySDToRound = 
		maxPathDepth 
		+ finalRoundAdd->getPipelineDepth() 
		+ 1 ; // finalRoundAdd output is registered -- TODO save it if slack allows
	


	// Signals-------------------------------

	addDelaySignal("fracRClosexMy",wF+3, 1);
	addDelaySignal("fracRCloseyMx",wF+3, 1);
	addDelaySignal("fracSignClose",1, 1);


	addSignal("signedExponentX",wE+1);
	addSignal("signedExponentY",wE+1);
	//		addSignal("invSignedExponentY",wE+1);
	
	addSignal("exceptionXSuperiorY");
	addSignal("exceptionXEqualY");
	
	addSignal("exponentDifferenceXY",wE+1);
	addSignal("exponentDifferenceYX",wE);
		
	addSignal("swap",1);
	
	addDelaySignal("exponentDifference",wE,1); 
			
	addDelaySignal("newY",wE+wF+3, 1);
			
	//close path		
	
	addSignal("fracXClose1",wF+3);
	addSignal("fracYClose1",wF+3);
	
	addDelaySignal("fracRClose1",wFX+2, lzocs->getPipelineDepth());


	addDelaySignal("resSign",1,lzocs->getPipelineDepth()+1 + finalRoundAdd->getPipelineDepth()+2);

	addDelaySignal("nZerosNew",lzocs->getCountWidth(),1);
	addDelaySignal("shiftedFrac",wFX+2,2);
	
	addDelaySignal("exponentResultClose",wEX+2, 1);
	addDelaySignal("resultCloseIsZero0",1, 1);
	addDelaySignal("roundClose0", 1, 1);
	// The following magical lines will delay the signal if needed to align it with that of the close path
	addDelaySignal("resultCloseIsZero", 1, farPathDepth - closePathDepth);
	addDelaySignal("roundClose", 1, farPathDepth - closePathDepth);
	addDelaySignal("resultBeforeRoundClose",wE+1 + wF+1, farPathDepth - closePathDepth);


	addSignal("syncResSign",1);
		
	//Far path
	addDelaySignal("shiftedOut",1, rightShifter->getPipelineDepth() + fracAddFar->getPipelineDepth() +3);
		
	addDelaySignal("shiftVal",sizeRightShift, 1);
	
	addSignal("fracNewY",wF+1);
	addSignal("shiftedFracY", 2*wF + 4);
	
	addDelaySignal("sticky",1,  1 + fracAddFar->getPipelineDepth() +1); 
	// + 1 for the xor stage, + 1 because output is registered
		
	addSignal("fracXfar",wF+4);
	addSignal("fracYfar",wF+4);
		
	addDelaySignal("cInAddFar",1,1);
	
	addDelaySignal("fracYfarXorOp",wF+4,1);
	addDelaySignal("fracResultfar0",wF+4,1);
		
	addSignal("fracResultFarNormStage",wF+4);
	addSignal("fracLeadingBits",2);

	addSignal("exponentResultfar0",wE+2);
	addSignal("expOperationSel",2);
	addSignal("exponentUpdate",wE+2);
	addDelaySignal("exponentResultFar1",wE+2,1);	
		
	addDelaySignal("fracResultFar1",wF, 1);
	addSignal("fracResultRoundBit",1);
	addSignal("fracResultStickyBit",1);
	// The following magical line will delay the signal if needed to align it with that of the close path
	addDelaySignal("resultBeforeRoundFar",wE+1 + wF+1, closePathDepth-farPathDepth);
	addDelaySignal("roundFar",1,closePathDepth-farPathDepth);

	addDelaySignal("roundFar1",1,1);
	
	addDelaySignal("selectClosePath",1, delaySDToRound);
	addSignal("syncClose",1);
	addSignal("resultBeforeRound",wE+1 + wF+1);				
	addSignal("UnderflowOverflow",2);
	addDelaySignal("zeroFromClose", 1, finalRoundAdd->getPipelineDepth()+1);

	addDelaySignal("sdExnXY",4, delaySDToRound);
	addSignal("syncExnXY",4);
		
	addDelaySignal("EffSub",1, delaySDToRound);
	addSignal("syncEffSub");

	addDelaySignal("newX",3+wE+wF, delaySDToRound+1);
	addSignal("syncX",3+wE+wF);
	addSignal("round",1);
		
	addDelaySignal("resultRounded",wE+wF+2, 1);
				
	addDelaySignal("pipeSignY",1, delaySDToRound);
	addSignal("syncSignY");
		
	addDelaySignal("resultNoExn",wE+wF+3,1);
	addSignal("finalResult",wE+wF+3);
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

	licence(o,"Bogdan Pasca, Florent de Dinechin (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);	
	
	//output VHDL components
	lzocs->outputVHDLComponent(o);
	rightShifter->outputVHDLComponent(o);
	
	dualSubClose->outputVHDLComponent(o);
	finalRoundAdd->outputVHDLComponent(o);
	fracAddFar->outputVHDLComponent(o);
	// intaddFar2->outputVHDLComponent(o);
	// intaddFar3->outputVHDLComponent(o);
		
	outputVHDLSignalDeclarations(o);	  
	beginArchitecture(o);
		
	

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
		o<<tab<<"exponentDifferenceXY <= signedExponentX - signedExponentY ;"<<endl;
		o<<tab<<"exponentDifferenceYX <= signedExponentY("<<wE-1<<" downto 0) - signedExponentX("<<wE-1<<" downto 0);"<<endl;

		// SWAP when: [excX=excY and expY>expX] or [excY>excX]
		o<<tab<<"swap <= (exceptionXEqualY and exponentDifferenceXY("<<wE<<")) or (not(exceptionXSuperiorY));"<<endl;

		// depending on the value of swap, assign the corresponding values to the newX and newY signals 
		o<<tab<<"newX <= Y when swap = '1' else X;"<<endl;
		o<<tab<<"newY <= X when swap = '1' else Y;"<<endl;
		o<<tab<<"exponentDifference <= exponentDifferenceYX when swap = '1' else exponentDifferenceXY("<<wE-1<<" downto 0);"<<endl;

		// determine if the fractional part of Y was shifted out of the operation //
		if (wE>sizeRightShift){
			o<<tab<<"shiftedOut <= "; 
			for (int i=wE-1;i>=sizeRightShift;i--)
				if (i==sizeRightShift)
					o << "exponentDifference("<<i<<")";
				else
					o << "exponentDifference("<<i<<") or ";
			o<<";"<<endl;
		}
		else
			o<<tab<<"shiftedOut <= '0';"<<endl; 

		//shiftVal=the number of positions that fracY must be shifted to the right				
		if (wE>sizeRightShift) {
			o<<tab<<"shiftVal <= " 
			 << "exponentDifference("<< sizeRightShift-1<<" downto 0"<<")"
			 << " when shiftedOut='0'"<<endl
			 <<tab << tab << "    else CONV_STD_LOGIC_VECTOR("<<wFX+3<<","<<sizeRightShift<<") ;" << endl; 
		}		
		else if (wE==sizeRightShift) {
			o<<tab<<"shiftVal <= exponentDifference;" << endl ;
		}
		else 	{ //  wE< sizeRightShift
			o<<tab<<"shiftVal <= CONV_STD_LOGIC_VECTOR(0,"<<sizeRightShift-wE <<") & exponentDifference;" <<	endl;
		}

		// compute EffSub as (signA xor signB) at cycle 1
		o<<tab<<"EffSub <= " << delaySignal("newX",1) << "("<<wEX+wFX<<") xor " << delaySignal("newY",1) << "("<<wEY+wFY<<");"<<endl;
		
		// compute the close/far path selection signal at cycle1 
		// the close path is considered only when (signA!=signB) and |exponentDifference|<=1 
		o<<tab<<"selectClosePath  <= EffSub when " << delaySignal("exponentDifference",1) << "("<<wER-1<<" downto "<<1<<") = ("<<wER-1<<" downto "<<1<<" => '0') else '0';"<<endl;
		

		// sdExnXY is a concatenation of the exception bits of X and Y
		o<<tab<<"sdExnXY <= " << delaySignal("newX",1) << "("<<wE+wF+2<<" downto "<<wE+wF+1<<") "
		 << "& " << delaySignal("newY",1) << "("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
		o<<tab<<"pipeSignY <= " << delaySignal("newY",1) << "("<<wE+wF<<");"<<endl;

		// swapDifferencePipelineDepth = 1; ==for now, this section has a constant depth of 1.         

		//=========================================================================|
		//                            close path                                   |
		//=========================================================================|
		
		o << endl << "-- Close Path --" << endl;
		
		// build the fraction signals
		// padding: [sign bit][inplicit "1"][fracX][guard bit]
		o<<tab<<"fracXClose1 <= \"01\" & " << delaySignal("newX",1) << "("<<wFX-1<<" downto "<<0<<") & '0';"<<endl;

		// the close path is considered when the |exponentDifference|<=1, so 
		// the alignment of fracY is of at most 1 position
		o<<tab<<"with " << delaySignal("exponentDifference",1) << "(0) select"<<endl;
		o<<tab<<"fracYClose1 <=  \"01\" & " << delaySignal("newY",1) << "("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
		o<<tab<<"               \"001\" & " << delaySignal("newY",1) << "("<<wF-1<<" downto "<<0<<")       when others;"<<endl;
		
		// substract the fraction signals for the close path; 

		// instanciate the box that computes X-Y and Y-X. Note that it could take its inputs before the swap (TODO ?)
		o<<tab<< "DualSubO: " << dualSubClose->getOperatorName() 
		 << "  -- pipelineDepth="<< dualSubClose->getPipelineDepth()<< endl;
		o<<tab<<tab<< "port map (" << endl; 
		if(isSequential()) {
			o<<tab<<tab<<tab<< "clk => clk, " << endl;
			o<<tab<<tab<<tab<< "rst => rst, " << endl;
		}
		o<<tab<<tab<<tab<< "X => fracXClose1 , " << endl; 
		o<<tab<<tab<<tab<< "Y => FracYClose1, " << endl; 
		o<<tab<<tab<<tab<< "RxMy => fracRClosexMy, " << endl; 
		o<<tab<<tab<<tab<< "RyMx => fracRCloseyMx " << endl; 
		o<<tab<<tab<<tab<< ");" << endl<<endl;
		// register the output -- TODO merge the mux with the last stage of the adder in case of sufficient slack

		o<<tab<< "fracSignClose <= " << delaySignal("fracRClosexMy",1) << "("<<wF+2<<");"<<endl;
		o<<tab<< "fracRClose1 <= " << delaySignal("fracRClosexMy",1) << "("<<wF+1<<" downto 0) when fracSignClose='0' else " << delaySignal("fracRCloseyMx",1) << "("<<wF+1<<" downto 0);"<<endl;
		int delayFromSD = dualSubClose->getPipelineDepth() + 2; 
		string closeDelayed1 = delaySignal("selectClosePath", delayFromSD);
		string XDelayed1 = delaySignal("newX", delayFromSD+1);
				
		//TODO check the test if significand is all zero is useful. 
		o << tab << "resSign <= '0' when "<<closeDelayed1<<"='1' and" << delaySignal("fracRClose1",1) << " = ("<<wF+1<<" downto 0 => '0') else"<<endl;
		// else sign(x) xor (close and sign(resclose))
		o << tab << "          " << XDelayed1 << "("<<wE+wF<<") xor (" << closeDelayed1 << " and " 
		  << delaySignal("fracSignClose", 1) << ");"<<endl;
			
		// LZC + Shifting. The number of leading zeros are returned together with the shifted input
		o<<tab<< "LZCOCS_component: " << lzocs->getOperatorName() 
		 << "  -- pipelineDepth="<< lzocs->getPipelineDepth()<< endl;
		o<<tab<<tab<< "port map (" << endl; 
		if(isSequential()) {
			o<<tab<<tab<<tab<< "clk => clk, " << endl;
			o<<tab<<tab<<tab<< "rst => rst, " << endl;
		}
		o<<tab<<tab<<tab<< "I => " << delaySignal("fracRClose1",1) << ", " << endl; 
		o<<tab<<tab<<tab<< "Count => nZerosNew, "<<endl;
		o<<tab<<tab<<tab<< "O => shiftedFrac " <<endl; 
		o<<tab<<tab<<tab<< ");" << endl<<endl;		
	
		// NORMALIZATION
		
		// shiftedFrac(0) is the round bit, shiftedFrac(1) is the parity bit, 
		// shiftedFrac(wF) is the leading one, to be discarded
		// the rounding bit is computed:
		o<<tab<< "roundClose0 <= " << delaySignal("shiftedFrac",1) << "(0) and " << delaySignal("shiftedFrac",1) << "(1);"<<endl;
		// Is the result zero? 
		o<<tab<< "resultCloseIsZero0 <= '1' when " 
		 << delaySignal("nZerosNew",1) 
		 << " = CONV_STD_LOGIC_VECTOR(" << wF+2 << ", " << lzocs->getCountWidth() 
		 << ") else '0';" << endl;

		// add two bits in order to absorb exceptions: 
		// the second 0 will become a 1 in case of overflow, 
		// the first 0 will become a 1 in case of underflow (negative biased exponent)
		o<<tab<< "exponentResultClose <= (\"00\" & "
		 << delaySignal("newX",closePathDepth-1+1) <<"("<<wE+wF-1<<" downto "<<wF<<")) "
		 <<"- (CONV_STD_LOGIC_VECTOR(0,"<<wE-lzocs->getCountWidth()+2<<") & " << delaySignal("nZerosNew",1) << ");"
		 <<endl;
		// concatenate exponent with fractional part before rounding so the possible carry propagation automatically increments the exponent 
		o<<tab<<"resultBeforeRoundClose <= " << delaySignal("exponentResultClose",1) << "("<<wE+1<<" downto 0) & " << delaySignal("shiftedFrac",2) << "("<<wF<<" downto 1);"<<endl; 
		o<<tab<< "roundClose <= " << delaySignal("roundClose0",1) << ";"<<endl;
		o<<tab<< "resultCloseIsZero <= " << delaySignal("resultCloseIsZero0",1) << ";"<<endl;




		
		//=========================================================================|
		//                              far path                                   |
		//=========================================================================|
		o << endl << "-- Far Path --" << endl;
		
		//input interface signals
		
		
		//add implicit 1 for frac1. 
 		o<<tab<< "fracNewY <= '1' & " << delaySignal("newY",1) << "("<<wF-1<<" downto 0);"<<endl;
		
								
		// shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) //		
		o<<tab<< "right_shifter_component: " << rightShifter->getOperatorName()
 		 << "  -- pipelineDepth="<< rightShifter->getPipelineDepth()<< endl;
		o<<tab<<tab<< "port map (" << endl; 
		if(isSequential()) {
			o<<tab<<tab<<tab<< "clk => clk, " << endl;
			o<<tab<<tab<<tab<< "rst => rst, " << endl;
		}
		o<<tab<<tab<<tab<< "X => fracNewY, " << endl; 
		o<<tab<<tab<<tab<< "S => " << delaySignal("shiftVal",1) << "," << endl; 
		o<<tab<<tab<<tab<< "R => shiftedFracY " <<endl; 
		o<<tab<<tab<<tab<< ");" << endl<<endl;		
			
		// compute sticky bit as the or of the shifted out bits during the alignment //
		o<<tab<< "sticky<= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;
		
		//pad fraction of Y [sign][shifted frac having inplicit 1][guard bits]
		o<<tab<< "fracYfar <= \"0\" & shiftedFracY("<<2*wF+3<<" downto "<<wF+1<<");"<<endl;	
		
		// depending on the signs of the operands, perform addition or substraction			
		// the result will be: a + (b xor operation) + operation, where operation=0=addition and operation=1=substraction
		// the operation selector is the xor between the signs of the operands
		// perform xor 
		o<<tab<<"fracYfarXorOp <= fracYfar xor ("<<wF+3<<" downto 0 => "<<delaySignal("EffSub",rightShifter->getPipelineDepth())<<");"<<endl;
		//pad fraction of X [sign][inplicit 1][fracX][guard bits]				
		o<<tab<< "fracXfar <= \"01\" & ("<<delaySignal("newX",rightShifter->getPipelineDepth()+2)<<"("<<wF-1<<" downto 0"<<")) & \"00\";"<<endl;
		o<<tab<< "cInAddFar <= " << delaySignal("EffSub",rightShifter->getPipelineDepth())<<" and not sticky;"<< endl;
		// perform carry in addition
		o<<tab<< "fracAddFar0: " << fracAddFar->getOperatorName()
		 << "  -- pipelineDepth="<< fracAddFar->getPipelineDepth()<< endl;
		o<<tab<<tab<< "port map (" << endl; 
		if(isSequential()) {
			o<<tab<<tab<<tab<< "clk => clk, " << endl;
			o<<tab<<tab<<tab<< "rst => rst, " << endl;
		}
		o<<tab<<tab<<tab<< "X => fracXfar, " << endl; 
		o<<tab<<tab<<tab<< "Y => "<<delaySignal("fracYfarXorOp",1) << ", " << endl; 
	             	// below, + 1 because XOR stage
		o<<tab<<tab<<tab<< "Cin => "<<delaySignal("cInAddFar",1) << "," << endl;
		o<<tab<<tab<<tab<< "R => fracResultfar0 " << endl; 
		o<<tab<<tab<<tab<< ");" << endl;
	

		o << tab << "-- 2-bit normalisation" <<endl; 

		o << tab << "fracResultFarNormStage <= " << delaySignal("fracResultfar0",1) << ";"<<endl;

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
		o<<tab<<tab<< "     "<<delaySignal("sticky", 1 + fracAddFar->getPipelineDepth() +1)<<" 	 when fracLeadingBits = \"00\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage(0) or  "<<delaySignal("sticky", 1 + fracAddFar->getPipelineDepth() +1)<<"   when fracLeadingBits = \"01\" "<<endl;
		o<<tab<<tab<< "else fracResultFarNormStage(1) or fracResultFarNormStage(0) or "<<delaySignal("sticky", 1 + fracAddFar->getPipelineDepth() +1)<<";"<<endl;
		// round bit
		o<<tab<< "roundFar1 <= fracResultRoundBit and (fracResultStickyBit or fracResultFar1(0));"<<endl;

		//select operation mode. This depends on wether or not the exponent must be adjusted after normalization		
		o<<tab<<"expOperationSel <= \"11\" when fracLeadingBits = \"00\" -- add -1 to exponent"<<endl;
		o<<tab<<"            else   \"00\" when fracLeadingBits = \"01\" -- add 0 "<<endl;
		o<<tab<<"            else   \"01\";                              -- add 1"<<endl;
	
		//the second operand depends on the operation selector
		o<<tab<<"exponentUpdate <= ("<<wE+1<<" downto 1 => expOperationSel(1)) & expOperationSel(0);"<<endl;
		
		// // result fracResultfar0 of this adder is registered
	   int delayFromSD2 = rightShifter->getPipelineDepth() + 1 + fracAddFar->getPipelineDepth() + 1;

		// the result exponent before normalization and rounding is = to the exponent of the first operand //
		o<<tab<<"exponentResultfar0<=\"00\" & ("<<delaySignal("newX", delayFromSD2+1)<<"("<<wF+wE-1<<" downto "<<wF<<"));"<<endl;
		
		o<<tab<<"exponentResultFar1 <= exponentResultfar0 + exponentUpdate;" << endl;
		
		// End of normalization stage
		o<<tab<<"resultBeforeRoundFar <= "
		 << delaySignal("exponentResultFar1",1) 
		 << " & " <<delaySignal("fracResultFar1",1)<<";" << endl;		
		o<<tab<< "roundFar <= " <<delaySignal("roundFar1",1)<<";" << endl;
		
				
		
		
		
		//=========================================================================|
		//                              Synchronization                            |
		//=========================================================================|				
		o<<endl<<"-- Synchronization of both paths --"<<endl;

				
		//synchronize the close signal
		o<<tab<<"syncClose <= "<<delaySignal("selectClosePath", maxPathDepth)<<";"<<endl; 
				
		// select between the results of the close or far path as the result of the operation 
		o<<tab<< "with syncClose select"<<endl;
		o<<tab<< "resultBeforeRound <= " << delaySignal("resultBeforeRoundClose",farPathDepth - closePathDepth) <<" when '1',"<<endl;
		o<<tab<< "                     " << delaySignal("resultBeforeRoundFar",  closePathDepth - farPathDepth) <<"   when others;"<<endl;
		o<<tab<< "with syncClose select"<<endl;
		o<<tab<< "round <= " << delaySignal("roundClose",farPathDepth - closePathDepth) <<" when '1',"<<endl;
		o<<tab<< "         " << delaySignal("roundFar", closePathDepth - farPathDepth) <<"   when others;"<<endl;

		o<<tab<< "zeroFromClose <= syncClose and "  
		 << delaySignal("resultCloseIsZero", farPathDepth - closePathDepth) << ";" <<endl;
		
		o << endl << "-- Rounding --" << endl;
		// perform the actual rounding //
		o<<tab<< "FinalRoundAdder: " << finalRoundAdd->getOperatorName() 
		<< "  -- pipelineDepth="<< finalRoundAdd->getPipelineDepth()<< endl;
		o<<tab<<tab<< "port map (" << endl; 
		if(isSequential()) {
			o<<tab<<tab<<tab<< "clk => clk, " << endl;
			o<<tab<<tab<<tab<< "rst => rst, " << endl;
		}
		o<<tab<<tab<<tab<< "X => resultBeforeRound , " << endl; 
		o<<tab<<tab<<tab<< "Y => ("<<1+wE+wF<<" downto 0 => '0'), " << endl; 
		o<<tab<<tab<<tab<< "Cin => round ," << endl;
		o<<tab<<tab<<tab<< "R => resultRounded   );" << endl<<endl;


		//sign bit
		o<<tab<<"syncEffSub <= "<<delaySignal("EffSub", delaySDToRound)<<";"<<endl;
		
		//X
		o<<tab<<"syncX <= "<<delaySignal("newX", delaySDToRound+1)<<";"<<endl;
		
		//signY
		o<<tab<<"syncSignY <= "<<delaySignal("pipeSignY", delaySDToRound)<<";"<<endl;
		
		// resSign comes from closer  TODO FIXME BUG HERE if some day close path is shorter than far path
		o<<tab<<"syncResSign <= "<<delaySignal("resSign", lzocs->getPipelineDepth()+ 2 + finalRoundAdd->getPipelineDepth()+1)<<";"<<endl;


		// compute the exception bits of the result considering the possible underflow and overflow 
		o<<tab<< "UnderflowOverflow <= " << delaySignal("resultRounded",1) 
		 << "("<<wE+1+wF<<" downto "<<wE+wF<<");"<<endl;

		o<<tab<< "with UnderflowOverflow select"<<endl;
		o<<tab<< "resultNoExn("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= \"10\" when \"01\", -- overflow"<<endl;
		o<<tab<< "                              \"00\" when \"10\" | \"11\",  -- underflow"<<endl;
		o<<tab<< "                              \"0\" &  not " << delaySignal("zeroFromClose", finalRoundAdd->getPipelineDepth() + 1) << "  when others; -- normal "<<endl;
  	
		o<<tab<< "resultNoExn("<<wE+wF<<" downto 0) <= syncResSign & "<<delaySignal("resultRounded",1) << "("<<wE+wF-1<<" downto 0);"<<endl;
	 
		o<<tab<< "syncExnXY <= "<<delaySignal("sdExnXY",delaySDToRound)<<";"<<endl;
		o<<tab<< "-- Exception bits of the result" << endl;
		o<<tab<< "with syncExnXY select"<<endl;
		o<<tab<< "  finalResult("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= resultNoExn("<<wE+wF+2<<" downto "<<wE+wF+1<<") when \"0101\","<<endl;
		o<<tab<< "                                 \"1\" & syncEffSub              when \"1010\","<<endl;
		o<<tab<< "                                 \"11\"            	            when \"1011\","<<endl;
		o<<tab<< "                                 syncExnXY(3 downto 2)         when others;"<<endl;
		o<<tab<< "-- Sign bit of the result" << endl;
		o<<tab<< "with syncExnXY select"<<endl;
		o<<tab<< "  finalResult("<<wE+wF<<") <= resultNoExn("<<wE+wF<<")         when \"0101\","<<endl;
		o<<tab<< "                     syncX("<<wE+wF<<") and syncSignY when \"0000\","<<endl;
		o<<tab<< "                     syncX("<<wE+wF<<")               when others;"<<endl;
		
		o<<tab<< "-- Exponent and significand of the result" << endl;
		o<<tab<< "with syncExnXY select"<<endl;
		o<<tab<< "  finalResult("<<wE+wF-1<<" downto 0) <= resultNoExn("<<wE+wF-1<<" downto 0) when \"0101\" | \"0100\" | \"0001\","<<endl;
		o<<tab<< "                                 syncX("<<wE+wF-1<<" downto "<<wE+wF-3<<") & ("<<wE+wF-4<<" downto 0 =>'0') when others;"<<endl;
		
		// assign results 

#if DEBUGVHDL
		// probes for debugging
		o<<tab<< "rexp <= finalResult("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		o<<tab<< "rfrac <= \"1\" & finalResult("<<wF-1<<" downto 0);"<<endl;
		o<<tab<< "rsig <= finalResult("<<wE+wF<<");"<<endl;
		o<<tab<< "rexc <= finalResult("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
#endif		
			
		o<<tab<< "R <= finalResult;"<<endl;
						   
		endArchitecture(o);

}








void FPAdder::emulate(TestCase * tc)
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
	FPNumber  fpr(wER, wFR, r);

	/* Set outputs */
	mpz_class svR = fpr.getSignalValue();
	tc->addExpectedOutput("R", svR);
}





void FPAdder::buildStandardTestCases(TestCaseList* tcl){
	TestCase *tc;

	// Regression tests 
	tc = new TestCase(this); 
	tc->addInput("X", 1.0);
	tc->addInput("Y", -1.0);
	emulate(tc);
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("X", 1.0);
	tc->addInput("Y", 0.0);
	emulate(tc);
	tcl->add(tc);

	
}

