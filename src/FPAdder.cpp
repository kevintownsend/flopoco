/*
 * Floating Point Adder for FloPoCo
 *
 * Author : Bogdan Pasca
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

/**
 * The FPAdder constructor
 * @param[in]		target		the target device
 * @param[in]		wEX			the the with of the exponent for the f-p number X
 * @param[in]		wFX			the the with of the fraction for the f-p number X
 * @param[in]		wEY			the the with of the exponent for the f-p number Y
 * @param[in]		wFY			the the with of the fraction for the f-p number Y
 * @param[in]		wER			the the with of the exponent for the addition result
 * @param[in]		wFR			the the with of the fraction for the addition result
 **/
FPAdder::FPAdder(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR) :
	Operator(target), wEX(wEX), wFX(wFX), wEY(wEY), wFY(wFY), wER(wER), wFR(wFR) {

	int i, j;
	ostringstream name, synch, synch2;

	/* Set up the status of the operator. Options = sequential|combinatorial */
	if (target->is_pipelined()) 
		set_sequential();
	else
		set_combinatorial(); 
		
	//parameter set up
	if (is_sequential()){
	}
	else{
		wF = wFX;
		wE = wEX;
		
	}	
		
		
	sizeRightShift = int ( ceil( log2(wF+4)));	


	/* The name has the format: FPAdder_wEX_wFX_wEY_wFY_wER_wFR where: wEX = width of X exponenet and wFX = width for the fractional part of X */
	name.str("");
	name<<"FPAdder_"<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR; 
	unique_name = name.str(); 
	
	
	
	
	
	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	add_input ("X", wEX + wFX + 3);
	add_input ("Y", wEY + wFY + 3);
	add_output("rexp", wE);
	add_output("rfrac", wF+1);
	add_output("rsig", 1 );
	add_output("rexc",2);
	
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		   
	
	// instantiate a Leading zero counter
	wOutLZC = int(ceil(log2(wFX+2+1)));
	leadingZeroCounter = new LZOC(target, wFX+2, wOutLZC);
	oplist.push_back(leadingZeroCounter);
	
	//instantiate a left shifter
	leftShifter = new Shifter(target,wFX+2,wFX+2,Left);
	oplist.push_back(leftShifter);
	
	//instantiate a right shifter
	rightShifter = new Shifter(target,wFX+1,wFX+3,Right);
	oplist.push_back(rightShifter);
	
	
	if (is_sequential()){
	}
	else{
		add_signal("exceptionXSuperiorY",1);
		add_signal("exceptionXEqualY",1);
		add_signal("exponentDifference0",wER+1);
		add_signal("swap",1);
		add_signal("exponentDifference1",wER); 
		add_signal("exponentDifference",wER);  
		add_signal("zeroExtendedSwap",wER);
		add_signal("newX",wEX+wFX+3);
		add_signal("newY",wEX+wFX+3);
		add_signal("signAB",1);
		add_signal("close",1);		
		add_signal("fracXClose1",wFX+3);
		add_signal("fracYClose1",wFX+3);
		add_signal("fracRClose0",wFX+3);
		add_signal("fracRClose1",wFX+2);
		add_signal("nZeros",wOutLZC);
		add_signal("exponentResultClose1",wEX+2);
		add_signal("exponentResultClose",wEX+2);
		add_signal("shiftedFrac",2*(wFX+2));
		add_signal("exponentConcatFrac0",wE+wF+2);
		add_signal("exponentConcatFrac1",wE+wF+2);
		add_signal("exceptionSelector",4);
		add_signal("rexc0",2);
		add_signal("rexc1",2);
		add_signal("shiftedFracY", 2*wF + 4);
		add_signal("fracNewY",wF+1);
		add_signal("fracXfar3",wF+5);
		add_signal("fracYfar3",wF+5);
		add_signal("fracResultfar0",wF+5);
		add_signal("fracResultfar0wSh",wF+5);
		add_signal("exponentResultfar0",wE+1);
		add_signal("fracResultfar1",wF+3);
		add_signal("exponentResultfar1",wE+1);
		add_signal("round",1);
		add_signal("fractionResultfar0",wF+2);
		add_signal("rs",1);
		add_signal("exponentResultfar",wE+1);
		add_signal("fractionResultfar",wF+1);
		add_signal("exponentResultn",wE+2);
		add_signal("fractionResultn",wF+1);
		add_signal("eMax",1);
		add_signal("eMin",1);
		add_signal("eTest",2);
		add_signal("nRn",wE+wF+3);
		add_signal("xAB",4);
		add_signal("nnR",wE+wF+3);
		add_signal("fractionResultCloseC",1+wF);
		add_signal("exponentResultCloseC",wE+2);
		add_signal("shiftedOut",1);
		add_signal("stiky",1);
		add_signal("borrow",1);
		add_signal("roundClose",1);
	}
	
	
	
		
}

/**
 * FPAdder destructor
 */
FPAdder::~FPAdder() {
}



/**
 * Method belonging to the Operator class overloaded by the FPAdder class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void FPAdder::output_vhdl(std::ostream& o, std::string name) {
  
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;

	int bias_val=int(pow(double(2),double(wEX-1)))-1;
	int i, j; 

	Licence(o,"Bogdan Pasca (2008)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	new_architecture(o,name);	
	leadingZeroCounter->output_vhdl_component(o);
	leftShifter->output_vhdl_component(o);
	rightShifter->output_vhdl_component(o);
	output_vhdl_signal_declarations(o);	  
	begin_architecture(o);
	
	if (is_sequential()){
		/* the sequential version */
		   
	} 
	else
	{  
		/* the combinational version */
		
	  /* ========================= Swap/Difference =============================*/
	  /* signal which indicates whether or not the exception bits of X are greater or equal than the exception bits of Y */		  
    o<<tab<<"exceptionXSuperiorY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") >= Y("<<wEY+wFY+2<<" downto "<<wEY+wF+1<<") else"<<endl;
		o<<tab<<"                       '0';"<<endl;
			
		/* signal which indicates whether or not the exception bits of X are equal to the exception bits of Y */		  
		o<<tab<<"exceptionXEqualY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") = Y("<<wEY+wFY+2<<" downto "<<wEY+wFY+1<<") else"<<endl;
		o<<tab<<"                    '0';"<<endl;

		/* make the difference between the exponents of X and Y. Pad exponents with sign bit before pexponentResultfarorming the substraction */
		o<<tab<<"exponentDifference0 <= (\"0\" & X("<<wEX+wFX-1<<" downto "<<wFX<<")) - (\"0\" & Y("<<wEY+wFY-1<<" downto "<<wFY<<"));"<<endl;

		/* swapping is performed when:    the exception bits of X and Y are equal and the exponent of Y is greater than the exponent of X (i.e. sign bit of exponent difference=1)
		                              or  the exception bits of Y are superior to the exception bits of X (EX: exc. bits of X are 00 and of Y are 01) */
		o<<tab<<"swap <= (exponentDifference0("<<wER<<") and exceptionXEqualY) or (not exceptionXSuperiorY);"<<endl;

		/* depending on the value of swap, assign the corresponding values to this section's outputs */
		o<<tab<<"newX <= Y when swap = '1' else"<<endl;
		o<<tab<<"        X;"<<endl;   
		o<<tab<<"newY <= X when swap = '1' else"<<endl;
    o<<tab<<"        Y;"<<endl;

		/* readjust for the exponents difference after the potential swap */
		o<<tab<<"exponentDifference1 <= exponentDifference0("<<wER-1<<" downto "<<0<<") xor ("<<wER-1<<" downto "<<0<<" => swap);"<<endl;
		o<<tab<<"zeroExtendedSwap    <= "<< zero_generator(wE-1,0)<<" & swap;"<<endl;
		o<<tab<<"exponentDifference  <= exponentDifference1  + zeroExtendedSwap;"<<endl;
		
		
		//==========================================================================
		//==========================================================================
		//CLOSE PATH
				
		/* check if we are found on the CLOSE or the FAR path */
		/* compute signAB as signA xor signB. The close path is considered only when signA!=signB */
		o<<tab<<"signAB <= X("<<wEX+wFX<<") xor Y("<<wEY+wFY<<");"<<endl;
		o<<tab<<"close  <= signAB when exponentDifference("<<wER-1<<" downto "<<1<<") = ("<<wER-1<<" downto "<<1<<" => '0') else"<<endl;
		o<<tab<<"          '0';"<<endl;

		/* build the fraction signals */
		/* the close pathe is considered when the |exponentDifference|<=1 */
		o<<tab<<"fracXClose1 <= \"01\" & newX("<<wFX-1<<" downto "<<0<<") & '0';"<<endl;
		o<<tab<<"with exponentDifference(0) select"<<endl;
		o<<tab<<"fracYClose1 <=  \"01\" & newY("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
		o<<tab<<"               \"001\" & newY("<<wF-1<<" downto "<<0<<")       when others;"<<endl;

		/* substract the fraction signals for the close path (for the close path the signs of the inputs are not be equal */
		o<<tab<<"fracRClose0 <= fracXClose1 - fracYClose1;"<<endl;
		
		/* if substraction result is negative - do a two's compelement of the result */
		o<<tab<<"with fracRClose0("<<wF+2<<") select"<<endl;
		o<<tab<<"fracRClose1 <= fracRClose0("<<wF+1<<" downto "<<0<<")                 when '0',"<<endl;
		o<<tab<<"               ("<<wF+1<<" downto 0 => '0') - fracRClose0("<<wF+1<<" downto 0) when others;"<<endl;
		
		/*  LZC */
		o<<tab<< "LZC_component: " << leadingZeroCounter->unique_name << endl;
		o<<tab<< "      port map ( I => fracRClose1, " << endl; 
		o<<tab<< "                 OZB => '0', " << endl; 
		o<<tab<< "                 O => nZeros " <<endl; 
		o<<tab<< "               );" << endl<<endl;		
				
		/* shift and round */
		o<<tab<< "left_shifter_component: " << leftShifter->unique_name << endl;
		o<<tab<< "      port map ( X => fracRClose1, " << endl; 
		o<<tab<< "                 S => nZeros, " << endl; 
		o<<tab<< "                 R => shiftedFrac " <<endl; 
		o<<tab<< "               );" << endl<<endl;		
		
		/* extend expoenet result before normalization and rounding with 2 bits, one for signaling underflow and for overflow */
		o<<tab<< "exponentResultClose1 <= \"00\" & newX("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl; 
		o<<tab<< "exponentResultClose <= exponentResultClose1 - ("<<zero_generator(wE+1-wOutLZC+1, 0)<<" &  nZeros);"<<endl;
		
		/* during fraction alignment, the fraction of Y is shifted at most one position to the right, so 1 extra bit is enough to perform rounding */
		o<<tab<< "roundClose <= shiftedFrac(0) and shiftedFrac(1);"<<endl;
		
		/* concatenate exponent with fractional part before rounding so the possible carry propagation automatically increments the exponent */
		o<<tab<<" exponentConcatFrac0 <= exponentResultClose("<<wE+1<<" downto 0) & shiftedFrac("<<wF<<" downto 1);"<<endl; 
		
		/* perform the actual rounding */
		o<<tab<<" exponentConcatFrac1 <= exponentConcatFrac0 + CONV_STD_LOGIC_VECTOR(roundClose,"<<1+1+wE+wF<<");"<<endl; 
		
		o<<tab<<" fractionResultCloseC <= \"1\" & exponentConcatFrac1("<<wF-1<<" downto 0);"<<endl;
		o<<tab<<" exponentResultCloseC <= exponentConcatFrac1("<<wF+wE+1<<" downto "<<wF<<");"<<endl;
		
		//=========================================================================
		//=========================================================================
		//Far path
		
		o<<tab<< "fracNewY <= '1' & newY("<<wF-1<<" downto 0);"<<endl;
		
		/* shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) */		
		o<<tab<< "right_shifter_component: " << rightShifter->unique_name << endl;
		o<<tab<< "      port map ( X => fracNewY, " << endl; 
		o<<tab<< "                 S => exponentDifference("<< sizeRightShift-1<<" downto 0"<<"), " << endl; 
		o<<tab<< "                 R => shiftedFracY " <<endl; 
		o<<tab<< "               );" << endl<<endl;		
		
		/* determine if the fractional part of Y was shifted out of the opperation */
		//TODO
		o<<tab<<" shiftedOut <= "; 
		for (int i=wE-1;i>=sizeRightShift;i--)
			if (((wE-1)==sizeRightShift)||(i==sizeRightShift))
				o<< "exponentDifference("<<i<<")";
			else
				o<< "exponentDifference("<<i<<") or ";
		o<<";"<<endl;
		
		/* compute stiky bit as the or of the shifted out bits during the alignment */
		o<<tab<< " stiky<= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;
		o<<tab<< " borrow <= stiky and signAB;"<<endl;//in the case of substraction and when the shifted out bits contain a 1, borrow:=1
				
		o<<tab<< "fracXfar3 <= \"01\" & newX("<<wF-1<<" downto 0) & \"0\" & \"0\" & \"0\";"<<endl;
  	o<<tab<< "fracYfar3 <= \"0\" & shiftedFracY("<<2*wF+3<<" downto "<<2*wF+4- (wF+4)+1<<") & stiky;"<<endl;	
		
		/* depending on the sign, perform addition or substraction */			
		o<<tab<< "with signAB select"<<endl;
    o<<tab<< "fracResultfar0 <= fracXfar3 - fracYfar3 when '1',"<<endl;
    o<<tab<< "        					fracXfar3 + fracYfar3 when others;"<<endl;
		
		/* if the second operand was shifted out of the operation, then the result of the operation becomes = to the first operand */		
		o<<tab<< "fracResultfar0wSh <= fracResultfar0 when shiftedOut='0' else fracXfar3;"<<endl;
		/* the result exponent before normalization and rounding is = to the exponent of the first operand */
		o<<tab<< "exponentResultfar0 <= \"0\" & newX("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		
		/*perform normalization and recalculation of the stiky bit */
	  o<<tab<<"fracResultfar1<=fracResultfar0wSh("<<wF+3<<" downto 2) & "
	  											<<"(fracResultfar0wSh(1) or fracResultfar0wSh(0)) when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"01\" else"<<endl;
    o<<tab<< "    					 fracResultfar0wSh("<<wF+2<<" downto 0) 	 when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"00\" else"<<endl;
    o<<tab<< "    					 fracResultfar0wSh("<<wF+4<<" downto 3) & (fracResultfar0wSh(2) or fracResultfar0wSh(1) or fracResultfar0wSh(0));"<<endl;
		/* readjust exponent after normalization */
	  o<<tab<< "exponentResultfar1 <= exponentResultfar0                                when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"01\" else"<<endl;
	  o<<tab<< "        							exponentResultfar0 - (("<<wE<<" downto 1 => '0') & \"1\") when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"00\" else"<<endl;
	  o<<tab<< "        							exponentResultfar0 + (("<<wE<<" downto 1 => '0') & \"1\");"<<endl;

		//rounding
		o<<tab<< "round <= fracResultfar1(1) and (fracResultfar1(2) or fracResultfar1(0));"<<endl;
		o<<tab<< "fractionResultfar0 <= (\"0\" & fracResultfar1("<<wF+2<<" downto 2)) + (("<<wF<<" downto 0 => '0') & round);"<<endl;

		o<<tab<< "exponentResultfar <= exponentResultfar1 + (("<<wE-1<<" downto 0 => '0') & fractionResultfar0("<<wF+1<<"));"<<endl;
		o<<tab<< "fractionResultfar <= (fractionResultfar0("<<wF+1<<") or fractionResultfar0("<<wF<<")) & fractionResultfar0("<<wF-1<<" downto 0);"<<endl;		
		
		//handle sign ==============================================================
		o << tab << "rs <= '0' when close = '1' and fracRClose1 = ("<<wF+1<<" downto 0 => '0') else"<<endl;
    o << tab << "          newX("<<wE+wF<<") xor (close and fracRClose0("<<wF+2<<"));"<<endl;
		
		/*select between the results of the close or far path as the result of the operation*/
		o<<tab<< "with close select"<<endl;
    o<<tab<< "exponentResultn <= exponentResultCloseC when '1',"<<endl;
    o<<tab<< "                   (\"0\" & exponentResultfar) when others;"<<endl;
		o<<tab<< "with close select"<<endl;
		o<<tab<< "fractionResultn <= fractionResultCloseC when '1',"<<endl;
		o<<tab<< "                   fractionResultfar when others;"<<endl;
		
		/* compute the exception bits of the result considering the possible underflow and overflow */
		o<<tab<< "eTest <= exponentResultn("<<wE+1<<" downto "<<wE<<");"<<endl;
		o<<tab<< "with eTest select"<<endl;
    o<<tab<< "nRn("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= \"10\" when \"01\","<<endl;
    o<<tab<< "                               \"00\" when \"10\" | \"11\","<<endl;
    o<<tab<< "                               \"01\" when others;"<<endl;
  	o<<tab<< "nRn("<<wE+wF<<" downto 0) <= rs & exponentResultn("<<wE-1<<" downto 0) & fractionResultn("<<wF-1<<" downto 0);"<<endl;
		
		o<<tab<< "xAB <= newX("<<wE+wF+2<<" downto "<<wE+wF+1<<") & newY("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
  
	  o<<tab<< "with xAB select"<<endl;
	  o<<tab<< "  nnR("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= nRn("<<wE+wF+2<<" downto "<<wE+wF+1<<") when \"0101\","<<endl;
	  o<<tab<< "                                 \"1\" & signAB              when \"1010\","<<endl;
	  o<<tab<< "                                 \"11\"                      when \"1011\","<<endl;
	  o<<tab<< "                                 xAB(3 downto 2)             when others;"<<endl;
	  o<<tab<< "with xAB select"<<endl;
	  o<<tab<< "  nnR("<<wE+wF<<") <= nRn("<<wE+wF<<")      when \"0101\","<<endl;
	  o<<tab<< "                newX("<<wE+wF<<")  and newY("<<wE+wF<<") when \"0000\","<<endl;
	  o<<tab<< "                newX("<<wE+wF<<")     when others;"<<endl;

	  o<<tab<< "with xAB select"<<endl;
	  o<<tab<< "  nnR("<<wE+wF-1<<" downto 0) <= nRn("<<wE+wF-1<<" downto 0) when \"0101\","<<endl;
	  o<<tab<< "                                newX("<<wE+wF-1<<" downto 0) when others;"<<endl;
			
		/* assign results */
		o<<tab<< "rexp <= nnR("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		o<<tab<< "rfrac <= \"1\" & nnR("<<wF-1<<" downto 0);"<<endl;
		o<<tab<< "rsig <= nnR("<<wE+wF<<");"<<endl;
		o<<tab<< "rexc <= nnR("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
		
	}
	o<< "end architecture;" << endl << endl;
}

 

/**
 * A zero generator method which takes as input two arguments and returns a string of zeros with quotes as stated by the second argurment
 * @param[in] n		    integer argument representing the number of zeros on the output string
 * @param[in] margins	integer argument determining the position of the quotes in the output string. The options are: -2= no quotes; -1=left quote; 0=both quotes 1=right quote
 * @return returns a string of zeros with the corresonding quotes given by margins
 **/
string FPAdder::zero_generator(int n, int margins)
{
ostringstream left,full, right, zeros;
int i;

	for (i=1; i<=n;i++)
		zeros<<"0";

	left<<"\""<<zeros.str();
	full<<left.str()<<"\"";
	right<<zeros.str()<<"\"";

	switch(margins){
		case -2: return zeros.str(); break;
		case -1: return left.str(); break;
		case  0: return full.str(); break;
		case  1: return right.str(); break;
		default: return full.str();
	}
}


TestCaseList FPAdder::generateStandardTestCases(int n)
{
	// TODO
	return TestCaseList();
}

TestCaseList FPAdder::generateRandomTestCases(int n)
{
	Signal sx = *get_signal_by_name("X");
	Signal sy = *get_signal_by_name("Y");
	Signal srexp = *get_signal_by_name("rexp");
	Signal srfra = *get_signal_by_name("rfrac");
	Signal srexc = *get_signal_by_name("rexc");
	Signal srsgn = *get_signal_by_name("rsig"); 

	TestCaseList tcl;	/* XXX: Just like Lyon's Transportion Company. :D */
	FloFP x(wEX, wFX), y(wEY, wFY), r(wER, wFR);

	for (int i = 0; i < n; i++)
	{
		x = getLargeRandom(sx.width()-2) + (mpz_class(1) << (wEX + wFX + 1));
		y = getLargeRandom(sy.width()-2) + (mpz_class(1) << (wEY + wFY + 1));
		r = x+y;

		TestCase tc;
		tc.addInput(sx, x.getSignalValue());
		tc.addInput(sy, y.getSignalValue());
		
		tc.addExpectedOutput(srexc, r.getExceptionSignalValue());
		if (r.getExceptionSignalValue() != 0)
		{
			tc.addExpectedOutput(srsgn, r.getSignSignalValue());
		}
		// Exponent and fraction are not defined for zero, inf or NaN
		if (r.getExceptionSignalValue() == 1)
		{
			tc.addExpectedOutput(srexp, r.getExponentSignalValue());
			tc.addExpectedOutput(srfra, r.getFractionSignalValue());
		}
		tcl.add(tc);
	}

	return tcl;
}



