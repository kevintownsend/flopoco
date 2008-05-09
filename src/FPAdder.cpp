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
		nBlock = (wF+5)/4; cout<<"nBlock="<<nBlock<<endl;
		wN = log2(nBlock)+3; cout<<"wN="<<wN<<endl;
		wLastBlock = wF+2-4*(nBlock-1);
		wRanges = pow(2,(wN-2))-1;
		wTree = pow(2,(wN-2))-2*(wN-3);  
		//wAddr = wN-3-i;
	}	
		
		
		

	/* The name has the format: FPAdder_wEX_wFX_wEY_wFY_wER_wFR where: wEX = width of X exponenet and wFX = width for the fractional part of X */
	name.str("");
	name<<"FPAdder_"<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR; 
	unique_name = name.str(); 
	
	
	
	
	
	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	add_input ("X", wEX + wFX + 3);
	add_input ("Y", wEY + wFY + 3);
	add_output("fr",wFX+2);
	add_output("ed",wER);  
	add_output("n", wN -2);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		   
	
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
		
		
		for (i=0;i<wN-3;i++){
			ostringstream treeName;
			treeName<<"tree"<<i;
			wAddr = wN-3-i;
			add_signal(treeName.str(),pow(2,wAddr+1)-1);
		}
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

		/* make the difference between the exponents of X and Y. Pad exponents with sign bit before performing the substraction */
		o<<tab<<"exponentDifference0 <= (\"0\" & X("<<wEX+wFX-1<<" downto "<<wFX<<")) - (\"0\" & Y("<<wEY+wFY-1<<" downto "<<wFY<<"));"<<endl;

		/* swapping is performed when:    the exception bits of X and Y are equal and the exponrnt of Y is greater than the exponent of X
		                              or  the exception bits of Y are superior to the exception bits of X (EX: exc. bits of X are 00 and of Y are 01) */
		o<<tab<<"swap <= (exponentDifference0("<<wER<<") and exceptionXEqualY) or (not exceptionXSuperiorY);"<<endl;

		/* depending on the value of swap, assign the corresponding values to this section's outputs */
		o<<tab<<"newX <= Y when swap = '1' else"<<endl;
		o<<tab<<"        X;"<<endl;   
		o<<tab<<"newY <= X when swap = '1' else"<<endl;
    o<<tab<<"        Y;"<<endl;

		/* readjust for the exponents difference after the potential swap */
		o<<tab<<"exponentDifference1 <= exponentDifference0("<<wER-1<<" downto "<<0<<") xor ("<<wER-1<<" downto "<<0<<" => swap);"<<endl;
		o<<tab<<"zeroExtendedSwap <= "<< zero_generator(wE-1,0)<<" & swap;"<<endl;
		o<<tab<<"exponentDifference  <= exponentDifference1  + zeroExtendedSwap;"<<endl;
		//==========================================================================
		
		/* check if we are found on the CLOSE or the FAR path */
		/* compute signAB as signA xor signB. The close path is considered only when signA!=signB */
		o<<tab<<"signAB <= X("<<wEX+wFX<<") xor Y("<<wEY+wFY<<");"<<endl;
		o<<tab<<"close <= signAB when exponentDifference("<<wER-1<<" downto "<<1<<") = ("<<wER-1<<" downto "<<1<<" => '0') else"<<endl;
		o<<tab<<"         '0';"<<endl;

		/* build the fraction signals */
		o<<tab<<"fracXClose1 <= \"01\" & newX("<<wFX-1<<" downto "<<0<<") & '0';"<<endl;
		o<<tab<<"with exponentDifference(0) select"<<endl;
		o<<tab<<"fracYClose1 <=  \"01\" & newY("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
		o<<tab<<"               \"001\" & newY("<<wF-1<<" downto "<<0<<")       when others;"<<endl;

		/* substract the fraction signals for the close path (for the close path the signs of the inputs must not be equal */
		o<<tab<<"fracRClose0 <= fracXClose1 - fracYClose1;"<<endl;
		
		o<<tab<<"with fracRClose0("<<wF+2<<") select"<<endl;
		o<<tab<<"fracRClose1 <= fracRClose0("<<wF+1<<" downto "<<0<<")                 when '0',"<<endl;
		o<<tab<<"               ("<<wF+1<<" downto 0 => '0') - fracRClose0("<<wF+1<<" downto 0) when others;"<<endl;
		
		/* =============================== LZC ===================================*/
/*
		//the width of the input entering the leading zero counter
		int wInput = wFX+2;
			
		int lgWInput = int(floor(log2(wFX+2)));
		int powLgWInput = pow(2,lgWInput);
		
		add_signal("f",powLgWInput);
		
		if (powLgWInput==wInput){
			//the width of the input is a power of 2
			f<=fracRClose1;
		}
		else{
			for (i=wInput-1;i>=0;i++)
				if (i==wInput-1)
					f(i+powLgWInput-wInput)<=fracRClose1(i);
				else
					f(i+powLgWInput-wInput)<=fracRClose1(i) or f(i+powLgWInput-wInput+1);
					
			//the rest of the string is assigned 1
			f(powLgWInput-wInput-1 downto 0) <= (powLgWInput-wInput-1 downto 0 => '1');
		}
		
	
		//the number of bits for the signal representing the number of leading zeros of fracRClose1
		int wN = int(ceil(log2(powLgWInput+1)));    
		
		ostringstream orBits;
				for (j=wN-1;j<=0;j++)
					if (j==(wN-1))
						orBits<<"f("<<j<<")";
					else
						orBits<<"f("<<j<<") or ";
				//negate or string and assign to output
				
		o<<tab<<"n("<<wN-1<<") <= not("<< orBits.str() <<");"<<endl;
		
		with n(wN-1) select
		n(wN-2) <= not(or( f(wN-1 downto wN/2))) when '0' else 
							 1;
		
		
		
		
		
		
		
		
		
		
		
		
		int wRange;
		int idxLow, idxHigh;
		int availableBits;
		
		//indexes on the input signal between we are working on
		idxLow  = 0;
		idxHigh = wFX + 1;
		
		abvailableBits = wFX + 2;
		
		for (i=wN-1; i>=0;i++){
			//the selected width in bits of the input for the current bit of the output
			wRange = pow(2,i);
			if (wRange <= availableBits){
				//create the or string
				ostringstream orBits;
				for (j=idxHigh;j<=idxHigh-wRange+1;j++)
					if (j==(idxHigh-wRange+1))
						orBits<<"f("<<j<<")";
					else
						orBits<<"f("<<j<<") or ";
				//negate or string and assign to output
				o<<tab<<"n("<<i<<") <= not("<< orBits.str() <<");"<<endl;
			}
			
			
		
		}
		
		
		
		//shift and round =========================================================
		//TODO
		
		
		
		//===============OUTPUTS==================================================
		o<<"n<=n0;"<<endl;
		o<<"fr<=fracRClose1;"<<endl;
		o<<"ed<=exponentDifference;"<<endl;

		
	*/					



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

