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

	ostringstream name, synch, synch2;

	name<<"FPAdder_"<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR; 
	setName(name.str()); 

	setOperatorType();

	setCopyrightString("Bogdan Pasca, Florent de Dinechin (2008)");		

	//parameter set up
	wF = wFX;
	wE = wEX;
	
	sizeRightShift = intlog2(wF+4);	// TODO Check: +3 ?

	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	addFPInput ("X", wEX, wFX);
	addFPInput ("Y", wEY, wFY);
	addFPOutput("R", wER, wFR);

	//=========================================================================|
	//                          Swap/Difference                                |
	// ========================================================================|
	vhdl<<"-- Exponent difference and swap  --"<<endl;
	vhdl<<tab<<declare("inX",wE+wF+3) << " <= X;"<<endl;
	vhdl<<tab<<declare("inY",wE+wF+3) << " <= Y;"<<endl;
	// signal which indicates whether or not the exception bits of X are greater or equal than/to the exception bits of Y		  
	vhdl<<tab<<declare("exceptionXSuperiorY") << " <= '1' when inX("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") >= inY("<<wEY+wFY+2<<" downto "<<wEY+wF+1<<") else '0';"<<endl;
		
	// signal which indicates whether or not the exception bits of X are equal to the exception bits of Y		  
	vhdl<<tab<<declare("exceptionXEqualY") << " <= '1' when inX("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") = inY("<<wEY+wFY+2<<" downto "<<wEY+wFY+1<<") else '0';"<<endl;
	
	// make the difference between the exponents of X and Y; expX - expY = expX + not(expY) + 1
	// pad exponents with sign bit
	vhdl<<tab<<declare("signedExponentX",wE+1) << " <= \"0\" & inX("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl;
	vhdl<<tab<<declare("signedExponentY",wE+1) << " <= \"0\" & inY("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl;
	vhdl<<tab<<declare("exponentDifferenceXY",wE+1) << " <= signedExponentX - signedExponentY ;"<<endl;
	vhdl<<tab<<declare("exponentDifferenceYX",wE) << " <= signedExponentY("<<wE-1<<" downto 0) - signedExponentX("<<wE-1<<" downto 0);"<<endl;
	
	// SWAP when: [excX=excY and expY>expX] or [excY>excX]
	vhdl<<tab<<declare("swap") << " <= (exceptionXEqualY and exponentDifferenceXY("<<wE<<")) or (not(exceptionXSuperiorY));"<<endl;
	
	if(0) 
		nextCycle();
	// depending on the value of swap, assign the corresponding values to the newX and newY signals 
	vhdl<<tab<<declare("newX",wE+wF+3) << " <= " << use("inY") << " when "<< use("swap") << " = '1' else " << use("inX") << ";"<<endl;
	vhdl<<tab<<declare("newY",wE+wF+3) << " <= " << use("inX") << " when "<< use("swap") << " = '1' else " << use("inY") << ";"<<endl;
	vhdl<<tab<<declare("exponentDifference",wE) << " <= " << use("exponentDifferenceYX") 
		 << " when "<< use("swap") << " = '1' else " << use("exponentDifferenceXY") << "("<<wE-1<<" downto 0);"<<endl;
	
	// determine if the fractional part of Y was shifted out of the operation //
	vhdl<<tab<<declare("shiftedOut") << " <= "; 
	if (wE>sizeRightShift){
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
		vhdl << "exponentDifference("<< sizeRightShift-1<<" downto 0"<<")"
			  << " when shiftedOut='0'"<<endl
			  <<tab << tab << "    else CONV_STD_LOGIC_VECTOR("<<wFX+3<<","<<sizeRightShift<<") ;" << endl; 
	}		
	else if (wE==sizeRightShift) {
		vhdl<<tab<<"exponentDifference;" << endl ;
	}
	else 	{ //  wE< sizeRightShift
		vhdl<<tab<<"CONV_STD_LOGIC_VECTOR(0,"<<sizeRightShift-wE <<") & exponentDifference;" <<	endl;
	}

	nextCycle();////////////////////////////////////////////////////////////////////////////////////

	// compute EffSub as (signA xor signB) at cycle 1
	vhdl<<tab<<declare("EffSub") << " <= " << use("newX") << "("<<wEX+wFX<<") xor " << use("newY") << "("<<wEY+wFY<<");"<<endl;
		
	// compute the close/far path selection signal at cycle1 
	// the close path is considered only when (signA!=signB) and |exponentDifference|<=1 
	vhdl<<tab<<declare("selectClosePath") << " <= EffSub when " << use("exponentDifference") << "("<<wER-1<<" downto "<<1<<") = ("<<wER-1<<" downto "<<1<<" => '0') else '0';"<<endl;
		

	// sdExnXY is a concatenation of the exception bits of X and Y, after swap, so exnX > exnY
	vhdl<<tab<<declare("sdExnXY",4) << " <= " << use("newX") << "("<<wE+wF+2<<" downto "<<wE+wF+1<<") "
		 << "& " << use("newY") << "("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
	vhdl<<tab<<declare("pipeSignY") << " <= " << use("newY") << "("<<wE+wF<<");"<<endl;
	
	
	//=========================================================================|
	//                            close path                                   |
	//=========================================================================|
		
	vhdl<< endl << "-- Close Path --" << endl;
	// we are already at cycle 1
	
	// build the fraction signals
	// padding: [sign bit][inplicit "1"][fracX][guard bit]
	vhdl<<tab<<declare("fracXClose1",wF+3) << " <= \"01\" & " << use("newX") << "("<<wFX-1<<" downto "<<0<<") & '0';"<<endl;
	
	// the close path is considered when the |exponentDifference|<=1, so 
	// the alignment of fracY is of at most 1 position
	vhdl<<tab<<"with " << use("exponentDifference") << "(0) select"<<endl;
	vhdl<<tab<<declare("fracYClose1",wF+3) << " <=  \"01\" & " << use("newY") << "("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
	vhdl<<tab<<"               \"001\" & " << use("newY") << "("<<wF-1<<" downto "<<0<<")       when others;"<<endl;
	
	// substract the fraction signals for the close path; 
	
	// instanciate the box that computes X-Y and Y-X. Note that it could take its inputs before the swap (TODO ?)
	dualSubClose = new 	IntDualSub(target, wF + 3, 0);
	dualSubClose->changeName(getName()+"_DualSubClose");
	oplist.push_back(dualSubClose);
	
	inPortMap  (dualSubClose, "X", "fracXClose1");
	inPortMap  (dualSubClose, "Y", "fracYClose1");
	outPortMap (dualSubClose, "RxMy","fracRClosexMy");
	outPortMap (dualSubClose, "RyMx","fracRCloseyMx");
	vhdl << instance(dualSubClose, "DualSubO");
	
	syncCycleFromSignal("fracRCloseyMx");/////////////////////////////////////////////////////////////////
	
	// register the output -- TODO merge the mux with the last stage of the adder in case of sufficient slack
	nextCycle();////////////////////////////////////////////////////////////////////////////////////

	vhdl<<tab<< declare("fracSignClose") << " <= " << use("fracRClosexMy") << "("<<wF+2<<");"<<endl;
	vhdl<<tab<< declare("fracRClose1",wFX+2) << " <= " << use("fracRClosexMy") << "("<<wF+1<<" downto 0) when fracSignClose='0' else " << use("fracRCloseyMx") << "("<<wF+1<<" downto 0);"<<endl;
				
	//TODO check the test if significand is all zero is useful. 
	vhdl<< tab << declare("resSign") << " <= '0' when " << use("selectClosePath") <<"='1' and " << use("fracRClose1") << " = ("<<wF+1<<" downto 0 => '0') else"<<endl;
	// else sign(x) xor (close and sign(resclose))
	vhdl<< tab << "          " << use("newX") << "("<<wE+wF<<") xor (" << use("selectClosePath") << " and " 
		 << use("fracSignClose") << ");"<<endl;
	
	// LZC + Shifting. The number of leading zeros are returned together with the shifted input
	lzocs = new LZOCShifterSticky(target, wFX+2, wFX+2,0, 0);
	if(lzocs->getCountWidth() > wE){
		cout << endl << "****WARNING****: wE < log2(wF), the generated VHDL is probably broken."<<endl;
		cout << "    Try increasing wE."<<endl;
	}

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
	vhdl<<tab<< declare("roundClose0") << " <= " << use("shiftedFrac") << "(0) and " << use("shiftedFrac") << "(1);"<<endl;
	// Is the result zero? 
	vhdl<<tab<< declare("resultCloseIsZero0") << " <= '1' when " 
		 << use("nZerosNew") 
		 << " = CONV_STD_LOGIC_VECTOR(" << wF+2 << ", " << lzocs->getCountWidth() 
		 << ") else '0';" << endl;

	// add two bits in order to absorb exceptions: 
	// the second 0 will become a 1 in case of overflow, 
	// the first 0 will become a 1 in case of underflow (negative biased exponent)
	vhdl<<tab<< declare("exponentResultClose",wEX+2) << " <= (\"00\" & "
		 << use("newX") <<"("<<wE+wF-1<<" downto "<<wF<<")) "
		 <<"- (CONV_STD_LOGIC_VECTOR(0,"<<wE-lzocs->getCountWidth()+2<<") & " << use("nZerosNew") << ");"
		 <<endl;


	nextCycle();////////////////////////////////////////////////////////////////////////////////////
	// concatenate exponent with fractional part before rounding so the possible carry propagation automatically increments the exponent 
	vhdl<<tab<<declare("resultBeforeRoundClose",wE+1 + wF+1) << " <= " << use("exponentResultClose") << "("<<wE+1<<" downto 0) & " << use("shiftedFrac") << "("<<wF<<" downto 1);"<<endl; 
	vhdl<<tab<< declare("roundClose") << " <= " << use("roundClose0") << ";"<<endl;
	vhdl<<tab<< declare("resultCloseIsZero") << " <= " << use("resultCloseIsZero0") << ";"<<endl;
	



	
	//=========================================================================|
	//                              far path                                   |
	//=========================================================================|


	vhdl<< endl << "-- Far Path --" << endl;
	// get back to first cycle after exp diff/swap 
	setCycleFromSignal("EffSub");/////////////////////////////////////////////////////////////////
	
	//add implicit 1 for frac1. 
	vhdl<<tab<< declare("fracNewY",wF+1) << " <= '1' & " << use("newY") << "("<<wF-1<<" downto 0);"<<endl;
	
	// shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) //		
	rightShifter = new Shifter(target,wFX+1,wFX+3,Right);
	rightShifter->changeName(getName()+"_RightShifter");
	oplist.push_back(rightShifter);
	inPortMap  (rightShifter, "X", "fracNewY");
	inPortMap  (rightShifter, "S", "shiftVal");
	outPortMap (rightShifter, "R","shiftedFracY");
	vhdl << instance(rightShifter, "RightShifterComponent");

	syncCycleFromSignal("shiftedFracY");/////////////////////////////////////////////////////////////////
	// register the output 
	nextCycle();////////////////////////////////////////////////////////////////////////////////////
			
	// compute sticky bit as the or of the shifted out bits during the alignment //
	vhdl<<tab<< declare("sticky") << " <= '0' when (" << use("shiftedFracY") << "("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;

	// one cycle only for the sticky (the far path has more time anyway)
	nextCycle();////////////////////////////////////////////////////////////////////////////////////
		
	//pad fraction of Y [sign][shifted frac having inplicit 1][guard bits]
	vhdl<<tab<< declare("fracYfar", wF+4) << " <= \"0\" & " << use("shiftedFracY") << "("<<2*wF+3<<" downto "<<wF+1<<");"<<endl;	
		
	// depending on the signs of the operands, perform addition or substraction			
	// the result will be: a + (b xor operation) + operation, where operation=0=addition and operation=1=substraction
	// the operation selector is the xor between the signs of the operands
	// perform xor 
	vhdl<<tab<<declare("fracYfarXorOp", wF+4) << " <= fracYfar xor ("<<wF+3<<" downto 0 => "<< use("EffSub")<<");"<<endl;
	//pad fraction of X [sign][inplicit 1][fracX][guard bits]				
	vhdl<<tab<< declare("fracXfar", wF+4) << " <= \"01\" & ("<< use("newX")<<"("<<wF-1<<" downto 0"<<")) & \"00\";"<<endl;
	vhdl<<tab<< declare("cInAddFar") << " <= " << use("EffSub")<<" and not " << use("sticky") << ";"<< endl;

	// perform carry in addition
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
	vhdl<< tab << declare("fracResultFarNormStage", wF+4) << " <= " << use("fracResultfar0") << ";"<<endl;
	
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
	vhdl<<tab<<tab<< "     "<< use("sticky")<<" 	 when fracLeadingBits = \"00\" "<<endl;
	vhdl<<tab<<tab<< "else fracResultFarNormStage(0) or  "<< use("sticky")<<"   when fracLeadingBits = \"01\" "<<endl;
	vhdl<<tab<<tab<< "else fracResultFarNormStage(1) or fracResultFarNormStage(0) or " << use("sticky")<<";"<<endl;
	// round bit
	vhdl<<tab<< declare("roundFar1") <<" <= fracResultRoundBit and (fracResultStickyBit or fracResultFar1(0));"<<endl;

	//select operation mode. This depends on wether or not the exponent must be adjusted after normalization		
	vhdl<<tab<<declare("expOperationSel",2) << " <= \"11\" when fracLeadingBits = \"00\" -- add -1 to exponent"<<endl;
	vhdl<<tab<<"            else   \"00\" when fracLeadingBits = \"01\" -- add 0 "<<endl;
	vhdl<<tab<<"            else   \"01\";                              -- add 1"<<endl;
	
	//the second operand depends on the operation selector
	vhdl<<tab<<declare("exponentUpdate",wE+2) << " <= ("<<wE+1<<" downto 1 => expOperationSel(1)) & expOperationSel(0);"<<endl;
		
	// the result exponent before normalization and rounding is = to the exponent of the first operand //
	vhdl<<tab<<declare("exponentResultfar0",wE+2) << "<=\"00\" & ("<< use("newX")<<"("<<wF+wE-1<<" downto "<<wF<<"));"<<endl;
		
	vhdl<<tab<<declare("exponentResultFar1",wE+2) << " <= exponentResultfar0 + exponentUpdate;" << endl;
		
	nextCycle();////////////////////////////////////////////////////////////////////////////////////

	// End of normalization stage
	vhdl<<tab<<declare("resultBeforeRoundFar",wE+1 + wF+1) << " <= "
		 << use("exponentResultFar1") << " & " << use("fracResultFar1")<<";" << endl;		
	vhdl<<tab<< declare("roundFar") << " <= " << use("roundFar1")<<";" << endl;
		
				
		
		
		
	//=========================================================================|
	//                              Synchronization                            |
	//=========================================================================|				
	vhdl<<endl<<"-- Synchronization of both paths --"<<endl;

	setCycleFromSignal("resultBeforeRoundFar");/////////////////////////////////////////////////////////////////
	syncCycleFromSignal("resultBeforeRoundClose");

	//synchronize the close signal
	vhdl<<tab<<declare("syncClose") << " <= "<< use("selectClosePath")<<";"<<endl; 
				
	// select between the results of the close or far path as the result of the operation 
	vhdl<<tab<< "with syncClose select"<<endl;
	vhdl<<tab<< declare("resultBeforeRound",wE+1 + wF+1) << " <= " << use("resultBeforeRoundClose") <<" when '1',"<<endl;
	vhdl<<tab<< "                     " << use("resultBeforeRoundFar") <<"   when others;"<<endl;
	vhdl<<tab<< "with syncClose select"<<endl;
	vhdl<<tab<< declare("round") << " <= " << use("roundClose") <<" when '1',"<<endl;
	vhdl<<tab<< "         " << use("roundFar") <<"   when others;"<<endl;
 
	vhdl<<tab<< declare("zeroFromClose") << " <= syncClose and " << use("resultCloseIsZero") << ";" <<endl;
		
	vhdl<< endl << "-- Rounding --" << endl;
	// perform the actual rounding //
	// finalRoundAdd will add the mantissa concatenated with exponent, two bits reserved for possible under/overflow 
	finalRoundAdd = new IntAdder(target, wE + wF + 2); 
	finalRoundAdd->changeName(getName()+"_finalRoundAdd");
	oplist.push_back(finalRoundAdd);
	
	ostringstream zero;
	zero<<"("<<1+wE+wF<<" downto 0 => '0') ";  
	
	inPortMap   (finalRoundAdd, "X", "resultBeforeRound");
	inPortMapCst(finalRoundAdd, "Y", zero.str() );
	inPortMap   (finalRoundAdd, "Cin", "round");
	outPortMap  (finalRoundAdd, "R","resultRounded");
	vhdl << instance(finalRoundAdd, "finalRoundAdder");

	syncCycleFromSignal("resultRounded");
	nextCycle();////////////////////////////////////////////////////////////////////////////////////
	//sign bit
	vhdl<<tab<<declare("syncEffSub") << " <= "<< use("EffSub")<<";"<<endl;
  	
	//X
	vhdl<<tab<<declare("syncX",3+wE+wF) << " <= "<< use("newX")<<";"<<endl;
		
	//signY
	vhdl<<tab<<declare("syncSignY") << " <= "<< use("pipeSignY")<<";"<<endl;
		
	// resSign comes from closer
	vhdl<<tab<<declare("syncResSign") << " <= "<< use("resSign")<<";"<<endl;

	// compute the exception bits of the result considering the possible underflow and overflow 
	vhdl<<tab<< declare("UnderflowOverflow",2) << " <= " << use("resultRounded") 
		 << "("<<wE+1+wF<<" downto "<<wE+wF<<");"<<endl;

	vhdl<<tab<< "with UnderflowOverflow select"<<endl;
	vhdl<<tab<< declare("resultNoExn",wE+wF+3) << "("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= \"10\" when \"01\", -- overflow"<<endl;
	vhdl<<tab<< "                              \"00\" when \"10\" | \"11\",  -- underflow"<<endl;
	vhdl<<tab<< "                              \"0\" &  not " << use("zeroFromClose") << "  when others; -- normal "<<endl;
  	
	vhdl<<tab<< "resultNoExn("<<wE+wF<<" downto 0) <= syncResSign & "<< use("resultRounded") << "("<<wE+wF-1<<" downto 0);"<<endl;
	 
	vhdl<<tab<< declare("syncExnXY", 4) << " <= "<< use("sdExnXY")<<";"<<endl;
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

FPAdder::~FPAdder() {
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

	// Set outputs 
	FPNumber  fpr(wER, wFR, r);
	mpz_class svR = fpr.getSignalValue();
	tc->addExpectedOutput("R", svR);

	// clean up
	mpfr_clears(x, y, r, NULL);
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
	tc->addInput("Y", FPNumber::plusDirtyZero);
	emulate(tc);
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("X", 1.0);
	tc->addInput("Y", FPNumber::minusDirtyZero);
	emulate(tc);
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("X", FPNumber::plusInfty);
	tc->addInput("Y", FPNumber::minusInfty);
	emulate(tc);
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("X", FPNumber::plusInfty);
	tc->addInput("Y", FPNumber::plusInfty);
	emulate(tc);
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("X", FPNumber::minusInfty);
	tc->addInput("Y", FPNumber::minusInfty);
	emulate(tc);
	tcl->add(tc);
	
}

