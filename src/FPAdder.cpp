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
 */
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
		wF = wFX;
		wE = wEX;
	}else{
		wF = wFX;
		wE = wEX;
	}	
		
		
	sizeRightShift = int ( ceil( log2(wF+4)));	



	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	add_input ("X", wEX + wFX + 3);
	add_input ("Y", wEY + wFY + 3);

	add_output("rexp", wE);
	add_output("rfrac", wF+1);
	add_output("rsig", 1 );
	add_output("rexc",2);

	// instantiate a Leading zero counter
	if (! is_sequential())
	{
		wOutLZC = int(ceil(log2(wFX+2+1)));
		leadingZeroCounter = new LZOC(target, wFX+2, wOutLZC);
		oplist.push_back(leadingZeroCounter);

		//instantiate a left shifter
		leftShifter = new Shifter(target,wFX+2,wFX+2,Left);
		oplist.push_back(leftShifter);
	}
		

	//instantiate a right shifter
	rightShifter = new Shifter(target,wFX+1,wFX+3,Right);
	oplist.push_back(rightShifter);
	
	if (is_sequential())
	{

		intaddClose1 = new IntAdder(target, wF + 3);
		intaddClose1->set_operator_name("intaddClose1");
		oplist.push_back(intaddClose1);

		intaddClose2 = new IntAdder(target, wF + 2);
		intaddClose2->set_operator_name("intaddClose2");
		oplist.push_back(intaddClose2);

		intaddClose3 = new IntAdder(target, wE + wF + 2);
		intaddClose3->set_operator_name("intaddClose3");
		oplist.push_back(intaddClose3);

		intaddFar1 = new IntAdder(target,wF+5);
		intaddFar1->set_operator_name("intaddFar1");
		oplist.push_back(intaddFar1);

		intaddFar2 = new IntAdder(target,wE+1);
		intaddFar2->set_operator_name("intaddFar2");
		oplist.push_back(intaddFar2);		

		intaddFar3 = new IntAdder(target,wE+1 + wF+1);
		intaddFar3->set_operator_name("intaddFar3");
		oplist.push_back(intaddFar3);		
				
		lzocs = new LZOCShifterSticky(target, wFX+2, wFX+2,0, 0);
		oplist.push_back(lzocs);
	}
	
	if (is_sequential()){
		
		//swap/Difference pipeline depth
		swapDifferencePipelineDepth = 1;
		
		closePathDepth = intaddClose1->pipeline_depth() + 
										 intaddClose2->pipeline_depth() + 
										 lzocs->pipeline_depth() + 
										 intaddClose3->pipeline_depth() + 4;
	
		farPathDepth = intaddFar1->pipeline_depth()+
									 intaddFar2->pipeline_depth()+
									 intaddFar3->pipeline_depth()+
									 rightShifter->pipeline_depth()
									 +5;
		
		cout<<"depths="<<intaddFar1->pipeline_depth()<<" "<<intaddFar2->pipeline_depth()<<" "<<intaddFar3->pipeline_depth()<<" "<<rightShifter->pipeline_depth()<<endl;
		
		
		cout<<endl<<"Close path depth = "<< closePathDepth;
		cout<<endl<<"Far path depth   = "<< farPathDepth;
											 
		maxPathDepth = max(closePathDepth, farPathDepth);							 
		
		cout<<endl<<"Max path depth   = "<< maxPathDepth<<endl;
		
		set_pipeline_depth(swapDifferencePipelineDepth + maxPathDepth + 1);
	}
	
	
	
	if (is_sequential()){
		
		add_signal("signedExponentX",wE+1);
		add_signal("signedExponentY",wE+1);
		add_signal("invSignedExponentY",wE+1);
				
		add_registered_signal("exceptionXSuperiorY",1);
		add_registered_signal("exceptionXEqualY",1);
		
		add_delay_signal("exponentDifference0",wE+1,1);
		
		add_delay_signal("swap",1,1);
		add_delay_signal("extendSwap",wE,1);
		
		add_signal("exponentDifference1",wER); 
		add_signal("exponentDifference",wER);  
			
		add_delay_signal("newX",wE+wF+3, 1);
		add_delay_signal("newY",wE+wF+3, 1);
							
		add_signal("sdX",wE+wF+3);
		add_signal("sdY",wE+wF+3);
		add_signal("sdSignAB",1);
		add_signal("sdClose",1);
		add_signal("sdExponentDifference",wE);
		
		
		//close path
		add_delay_signal("cClose",1,intaddClose1->pipeline_depth() + intaddClose2->pipeline_depth()+2);
		add_delay_signal("cX",wE+wF+3,intaddClose1->pipeline_depth() + intaddClose2->pipeline_depth());
		add_signal("cY",wE+wF+3);
		add_signal("cExponentDifference",wE);
		
		add_delay_signal("expXClose",wE, intaddClose1->pipeline_depth() + intaddClose2->pipeline_depth() + 3 + lzocs->pipeline_depth() ); 
		
		add_signal("fracXClose1",wF+3);
		add_signal("fracYClose1",wF+3);
		add_signal("InvFracYClose1",wFX+3);
		add_delay_signal("fracRClose0",wF+3,intaddClose2->pipeline_depth()+2);
		
		
		add_delay_signal("fracSignClose",1,1);
		add_delay_signal("fracRClose1xor", wF+2,1);
		add_delay_signal("fracRClose1",wFX+2, lzocs->pipeline_depth());


		add_delay_signal("crs",1,lzocs->pipeline_depth() + intaddClose3->pipeline_depth()+2);
	//	add_delay_signal("exponentResultClose1",wEX+2,leadingZeroCounter->pipeline_depth()+leftShifter->pipeline_depth());

		add_delay_signal("nZerosNew",lzocs->getCountWidth(),lzocs->pipeline_depth());
		add_delay_signal("shiftedFrac",wFX+2,2);
				
		add_signal("exponentResultClose",wEX+2);
		add_signal("exponentResultClose1",wEX);
		add_signal("roundClose",1);
		add_signal("exponentConcatFrac0",wE+wF+2);
		add_delay_signal("exponentConcatFrac1",wE+wF+2,2);
		
		
		add_signal("fractionResultCloseC",wF + 1);
		add_signal("exponentResultCloseC",wE + 2);
		add_signal("crsOut",1);
		
		
		//Far path
		add_delay_signal("fX",wE+wF+3,rightShifter->pipeline_depth() + intaddFar1->pipeline_depth() );//CHECK XXX
		add_delay_signal("fracFX",wF,rightShifter->pipeline_depth() + intaddFar1->pipeline_depth()+2);
		add_delay_signal("expFX",wE,rightShifter->pipeline_depth() + intaddFar1->pipeline_depth()+3);
		
		
		add_signal("fY",wE+wF+3);
		add_signal("fExponentDifference",wE);
		add_delay_signal("fSignAB",1,rightShifter->pipeline_depth());
		
		add_delay_signal("shiftedOut",1, rightShifter->pipeline_depth() + intaddFar1->pipeline_depth() +2);
		
		add_signal("selectionSignal",sizeRightShift);
		
		add_signal("fracNewY",wF+1);
		add_signal("shiftedFracY", 2*wF + 4);
		
		add_signal("stiky",1);
		
		add_delay_signal("fracXfar3",wF+5,2); //XXX
		add_signal("fracYfar3",wF+5);
		
		add_delay_signal("opSelector",1,2);
		add_delay_signal("fracYfar3XorOp",wF+5,2);
		add_delay_signal("fracResultfar0",wF+5,2);
		
		
		add_signal("fracResultfar0wSh",wF+5);
		add_signal("exponentResultfar0",wE+1);
		
		add_delay_signal("exponentResultfar1",wE+1,2);	
		add_delay_signal("expOperationSel",2,2);
		
		
		add_delay_signal("expSecOperand",wE+1,2);
		
		add_delay_signal("fracResultfar1",wF+3,intaddFar2->pipeline_depth()+2);
		add_signal("expConcatFrac",wE+1 + wF+1);
		add_delay_signal("expConcatFracResult",wE+1 + wF+1,2);
		add_delay_signal("round",1,intaddFar2->pipeline_depth()+2);
		
		
		add_signal("exponentResultfar",wE+1);
		add_signal("fractionResultfar",wF+1);
		
		add_delay_signal("pipeClose",1, maxPathDepth);
		add_signal("syncClose",1);
		
		add_signal("eTest",2);
		add_signal("sdxXY",4);
		
		add_delay_signal("pipeXAB",4,maxPathDepth);
		add_delay_signal("syncXAB",4,2);
		
		add_delay_signal("pipeSignAB",1,maxPathDepth);
		add_delay_signal("syncSignAB",1,2);

		add_delay_signal("pipeX",3+wE+wF,maxPathDepth);
		add_delay_signal("syncX",3+wE+wF,2);
				
		add_signal("sdSignY",1);
		add_delay_signal("pipeSignY",1,maxPathDepth);
		add_delay_signal("syncSignY",1,2);
		
		add_delay_signal("nRn",wE+wF+3,2);
		add_signal("nnR",wE+wF+3);
		add_signal("exponentResultn",wE+2);
		add_signal("fractionResultn",wF+1);
		
		if (closePathDepth == farPathDepth){
		
			add_signal("expResClose",wE+2);
			add_signal("expResFar",wE+1);
			add_signal("fracResClose",wF+1);
			add_signal("fracResFar",wF+1);
			add_signal("syncRS",1);	
			
		}else
			if (closePathDepth > farPathDepth){
				add_delay_signal("pipeExpResFar",wE+1,closePathDepth-farPathDepth);
				add_delay_signal("pipeFracResFar",wF+1,closePathDepth-farPathDepth);
						
				add_signal("expResClose",wE+2);
				add_signal("expResFar",wE+1);
				add_signal("fracResClose",wF+1);
				add_signal("fracResFar",wF+1);
				add_signal("syncRS",1);	
			}else	{
				add_delay_signal("pipeExpResClose",wE+2,farPathDepth-closePathDepth);
				add_delay_signal("pipeFracResClose",wF+1,farPathDepth-closePathDepth);
				add_delay_signal("pipeSyncRS",1,farPathDepth-closePathDepth);		
			
				add_signal("expResClose",wE+2);
				add_signal("expResFar",wE+1);
				add_signal("fracResClose",wF+1);
				add_signal("fracResFar",wF+1);
				add_signal("syncRS",1);	
						
			}
			
		cout<<"signals pass";	
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

/** Method for setting the operator name
	* @param[in] prefix the prefix that is going to be placed to the default name
	*                   of the operator
	* @param[out] postfix the postfix that is going to be placed to the default
	*                     name of the operator
	*/
void FPAdder::set_operator_name(std::string prefix, std::string postfix){
	ostringstream name;
	/* The name has the format: FPAdder_wEX_wFX_wEY_wFY_wER_wFR where: 
	   wEX = width of X exponenet and 
	   wFX = width for the fractional part of X */
	name<<"FPAdder_"<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR; 
	unique_name = name.str(); 
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
	
	//output VHDL components
	if (is_sequential())
		lzocs->output_vhdl_component(o);
	else{
		leadingZeroCounter->output_vhdl_component(o);
		leftShifter->output_vhdl_component(o);
	}	
	rightShifter->output_vhdl_component(o);
	
	intaddClose1->output_vhdl_component(o);
	intaddClose2->output_vhdl_component(o);
	intaddClose3->output_vhdl_component(o);
	intaddFar1->output_vhdl_component(o);
	intaddFar2->output_vhdl_component(o);
	intaddFar3->output_vhdl_component(o);
		
	output_vhdl_signal_declarations(o);	  
	begin_architecture(o);
		
	if (is_sequential()){
	
		//=========================================================================| 
		//                                                                         |
		//                             SEQUENTIAL                                  |
		//                                                                         |
		//=========================================================================|
	
		output_vhdl_registers(o); o<<endl;
		//=========================================================================|	
	  //                          Swap/Difference                                |
	  // ========================================================================|
	  o<<"-- Swap Difference --"<<endl;
    // signal which indicates whether or not the exception bits of X are greater or equal than the exception bits of Y		  
    o<<tab<<"exceptionXSuperiorY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") >= Y("<<wEY+wFY+2<<" downto "<<wEY+wF+1<<") else"<<endl;
		o<<tab<<"                       '0';"<<endl;
			
		// signal which indicates whether or not the exception bits of X are equal to the exception bits of Y		  
		o<<tab<<"exceptionXEqualY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") = Y("<<wEY+wFY+2<<" downto "<<wEY+wFY+1<<") else"<<endl;
		o<<tab<<"                    '0';"<<endl;

		// make the difference between the exponents of X and Y; expX - expY = expX + not(expY) + 1
		// pad exponents with sign bit
		o<<tab<<"signedExponentX <= \"0\" & X("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl;
		o<<tab<<"signedExponentY <= \"0\" & Y("<<wEX+wFX-1<<" downto "<<wFX<<");"<<endl;
		// make not(expY)	
		o<<tab<<"invSignedExponentY <= not(signedExponentY);"<<endl;
		// perform addition with carry in
		o<<tab<<"exponentDifference0 <= signedExponentX + invSignedExponentY + '1';"<<endl;
		
  	// SWAP when: [excX=excY and expY>expX] or [excY>excX]
		o<<tab<<"swap <= (exceptionXEqualY and exponentDifference0("<<wE<<")) or (not(exceptionXSuperiorY));"<<endl;

		// depending on the value of swap, assign the corresponding values to the newX and newY signals 
		o<<tab<<"newX <= Y when swap = '1' else"<<endl;
		o<<tab<<"        X;"                    <<endl;   
		o<<tab<<"newY <= X when swap = '1' else"<<endl;
    o<<tab<<"        Y;"                    <<endl;

		// readjust for the exponents difference after the potential swap
		// sign extend swap to the size of the exponent
		o<<tab<<"extendSwap <= ("<<wE-1<<" downto "<<0<<" => swap);"<<endl;
		// compute NewExponentDifference = (OldExponentDifference xor extendedSWAP) + swap.
		o<<tab<<"exponentDifference1 <= exponentDifference0_d("<<wE-1<<" downto "<<0<<") xor extendSwap_d;"<<endl;
		o<<tab<<"exponentDifference <= exponentDifference1 + swap_d;"<<endl;
	
		// compute sdSignAB as (signA xor signB)
		o<<tab<<"sdSignAB <= newX_d("<<wEX+wFX<<") xor newY_d("<<wEY+wFY<<");"<<endl;
		
		// check if we are found on the CLOSE or the FAR path 
		// the close path is considered only when (signA!=signB) and |exponentDifference|<=1 
		o<<tab<<"sdClose  <= sdSignAB when exponentDifference("<<wER-1<<" downto "<<1<<") = ("<<wER-1<<" downto "<<1<<" => '0') else"<<endl;
		o<<tab<<"           '0';"<<endl;
		
		// assign interface signals
		o<<tab<<"sdX <= newX_d;"<<endl;
		o<<tab<<"sdY <= newY_d;"<<endl;
		o<<tab<<"sdExponentDifference <= exponentDifference;"<<endl;
		// sdxXY is a concatenation of the exception bits of X and Y		
		o<<tab<<"sdxXY <= newX_d("<<wE+wF+2<<" downto "<<wE+wF+1<<") & newY_d("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
		o<<tab<<"sdSignY <= newY_d("<<wE+wF<<");"<<endl;
	
		// swapDifferencePipelineDepth = 1; ==for now, this section has a constant depth of 1.         
	
		
		//=========================================================================|
		//                            close path                                   |
		//=========================================================================|
		
		o<<"-- Close Path --"<<endl;
		// input interface signal assignment
		o<<tab<<"cX <= sdX;"<<endl;
		o<<tab<<"cY <= sdY;"<<endl;
		o<<tab<<"expXClose <= cX("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		o<<tab<<"cExponentDifference <= sdExponentDifference;"<<endl;
		o<<tab<<"cClose <= sdClose;"<<endl;	
		
		// build the fraction signals
			// padding: [sign bit][inplicit "1"][fracX][guard bit]
			o<<tab<<"fracXClose1 <= \"01\" & cX("<<wFX-1<<" downto "<<0<<") & '0';"<<endl;
			// the close path is considered when the |exponentDifference|<=1, so 
			// the alignment of fracY is of at most 1 position
			o<<tab<<"with cExponentDifference(0) select"<<endl;
			o<<tab<<"fracYClose1 <=  \"01\" & cY("<<wF-1<<" downto "<<0<<") & '0' when '0',"<<endl;
			o<<tab<<"               \"001\" & cY("<<wF-1<<" downto "<<0<<")       when others;"<<endl;

		// substract the fraction signals for the close path; 
		// substraction fX - fY = fX + not(fY)+1
			// perform inversion
			o<<tab<<"InvFracYClose1 <= not(fracYClose1);"<<endl;
			// perform addition with carry in = 1
			o<<tab<< "int_adder_componentc1: " << intaddClose1->getOperatorName() << endl;
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
			o<<tab<< "int_adder_componentc2: " << intaddClose2->getOperatorName() << endl;
			o<<tab<< "  port map ( X => fracRClose1xor_d , " << endl; 
			o<<tab<< "             Y => ("<<wF+1<<" downto 0 => '0'), " << endl; 
			o<<tab<< "             Cin => fracSignClose_d ," << endl;
			o<<tab<< "             R => fracRClose1, " << endl; 
			o<<tab<< "             clk => clk, " << endl;
			o<<tab<< "             rst => rst " << endl;
			o<<tab<< "               );" << endl<<endl;
		
		ostringstream cCloseSNPostPipeline;
		ostringstream cXSNPostPipeline;
		
		cCloseSNPostPipeline << get_delay_signal_name("cClose",intaddClose1->pipeline_depth()+intaddClose2->pipeline_depth()+2);
		cXSNPostPipeline << get_delay_signal_name("pipeX",intaddClose1->pipeline_depth()+intaddClose2->pipeline_depth()+2);
				
		o << tab << "crs <= '0' when "<<cCloseSNPostPipeline.str()<<"='1' and fracRClose1 = ("<<wF+1<<" downto 0 => '0') else"<<endl;
    o << tab << "          "<<cXSNPostPipeline.str()<<"("<<wE+wF<<") xor ("<<cCloseSNPostPipeline.str()<<" and "
		                        <<get_delay_signal_name("fracRClose0",intaddClose2->pipeline_depth()+2)<<"("<<wF+2<<"));"<<endl;
		
			
		// LZC + Shifting. The number of leading zeros are returned together with the sifted input
		o<<tab<< "LZCOCS_component: " << lzocs->getOperatorName() << endl;
		o<<tab<< "      port map ( I => fracRClose1_d, " << endl; 
		o<<tab<< "                 Count => nZerosNew, "<<endl;
		o<<tab<< "                 O => shiftedFrac, " <<endl; 
		o<<tab<< "                 clk => clk, " << endl;
		o<<tab<< "                 rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;		
	
		// NORMALIZATION
		int delayExp;
		delayExp=lzocs->pipeline_depth() + intaddClose2->pipeline_depth() + intaddClose1->pipeline_depth()+3;
		o<<tab<< "exponentResultClose1<= "<<get_delay_signal_name("expXClose",delayExp)
		                                  <<" - (CONV_STD_LOGIC_VECTOR(0,"<<wE-lzocs->getCountWidth()<<") & nZerosNew);"<<endl;
		
		// ROUNDING
		// during fraction alignment, the fraction of Y is shifted at most one position to the right, so 1 extra bit is enough to perform rounding 
			// the rounding bit is computed:
			o<<tab<< "roundClose <= shiftedFrac(0) and shiftedFrac(1);"<<endl;
			// add two bits in order to absorb exceptions 		 		
			o<<tab<< "exponentResultClose <= \"00\" & exponentResultClose1;"<<endl;
			// concatenate exponent with fractional part before rounding so the possible carry propagation automatically increments the exponent 
			o<<tab<<" exponentConcatFrac0 <= exponentResultClose("<<wE+1<<" downto 0) & shiftedFrac("<<wF<<" downto 1);"<<endl; 
			// perform the actual rounding //
			o<<tab<< "int_adder_componentc3: " << intaddClose3->getOperatorName() << endl;
			o<<tab<< "  port map ( X => exponentConcatFrac0 , " << endl; 
			o<<tab<< "             Y => ("<<1+wE+wF<<" downto 0 => '0'), " << endl; 
			o<<tab<< "             Cin => roundClose ," << endl;
			o<<tab<< "             R => exponentConcatFrac1, " << endl; 
			o<<tab<< "             clk => clk, " << endl;
			o<<tab<< "             rst => rst " << endl;
			o<<tab<< "               );" << endl<<endl;
	
		// assign output interface signals			
		o<<tab<<" fractionResultCloseC <= \"1\" & exponentConcatFrac1_d("<<wF-1<<" downto 0);"<<endl;
		o<<tab<<" exponentResultCloseC <= exponentConcatFrac1_d("<<wF+wE+1<<" downto "<<wF<<");"<<endl;
	
		o<<tab<<" crsOut <= "<<get_delay_signal_name("crs", lzocs->pipeline_depth() + intaddClose3->pipeline_depth()+2)<<";"<<endl;

		
		//=========================================================================|
		//                              far path                                   |
		//=========================================================================|
		o<<"-- Far Path --"<<endl;
		
		//input interface signals
		o<<tab<<"fX <= sdX;"<<endl;
		o<<tab<<"fracFX<=fX("<<wF-1<<" downto 0"<<");"<<endl;
		o<<tab<<"expFX<=fX("<<wF+wE-1<<" downto "<<wF<<");"<<endl;
		o<<tab<<"fY <= sdY;"<<endl;
		o<<tab<<"fExponentDifference <= sdExponentDifference;"<<endl;
		o<<tab<<"fSignAB <= sdSignAB;"<<endl;
		
		// determine if the fractional part of Y was shifted out of the opperation //
		if (wE-1>sizeRightShift){
			o<<tab<<" shiftedOut <= "; 
			for (int i=wE-1;i>=sizeRightShift;i--)
				if (((wE-1)==sizeRightShift)||(i==sizeRightShift))
					o<< "fExponentDifference("<<i<<")";
				else
					o<< "fExponentDifference("<<i<<") or ";
			o<<";"<<endl;
		}
		else
			o<<tab<<" shiftedOut <= '0';"<<endl; 
						
		//add inplicit 1 for frac1. Fêter la fin de l'année Universitaire
 		o<<tab<< "fracNewY <= '1' & fY("<<wF-1<<" downto 0);"<<endl;
		
		//selectioSignal=the number of positions that fracY must be shifted to the right				
		if (wE-1>sizeRightShift)
			o<<tab<<"selectionSignal <= fexponentDifference("<< sizeRightShift-1<<" downto 0"<<"); " << endl; 
		else
			o<<tab<<"selectionSignal <= CONV_STD_LOGIC_VECTOR(0,"<<sizeRightShift-(wE-1)<<") & fexponentDifference("<< wE-1<<" downto 0"<<"); " << endl; 			
								
		// shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) //		
		o<<tab<< "right_shifter_component: " << rightShifter->getOperatorName() << endl;
		o<<tab<< "      port map ( X => fracNewY, " << endl; 
		o<<tab<< "                 S => selectionSignal, " << endl; 
		o<<tab<< "                 R => shiftedFracY, " <<endl; 
		o<<tab<< "                 clk => clk, " << endl;
		o<<tab<< "                 rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;		
			
		// compute stiky bit as the or of the shifted out bits during the alignment //
		o<<tab<< "stiky<= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;
						
		//pad fraction of X [sign][inplicit 1][fracX][guard bits]				
		o<<tab<< "fracXfar3 <= \"01\" & "<<get_delay_signal_name("fracFX",rightShifter->pipeline_depth())<<" & \"000\";"<<endl;
  	//pad fraction of Y [sign][shifted frac having inplicit 1][guard bits]
  	o<<tab<< "fracYfar3 <= \"0\" & shiftedFracY("<<2*wF+3<<" downto "<<2*wF+4- (wF+4)+1<<") & stiky;"<<endl;	
		
		// depending on the signs of the operands, perform addition or substraction			
		// the result will be: a + (b xor operation) + operation, where operation=0=addition and operation=1=substraction
			// the operation selector is the xor between the signs of the operands
			o<<tab<<"opSelector <= "<<get_delay_signal_name("fSignAB",rightShifter->pipeline_depth())<<";"<<endl;
			// perform xor 
			o<<tab<<"fracYfar3XorOp <= fracYfar3 xor ("<<wF+4<<" downto 0 => opSelector);"<<endl;
			// perform carry in addition
			o<<tab<< "int_adder_componentf1: " << intaddFar1->getOperatorName()<< endl;
			o<<tab<< "  port map ( X => fracXfar3_d, " << endl; 
			o<<tab<< "             Y => fracYfar3XorOp_d, " << endl; 
			o<<tab<< "             Cin => opSelector_d ," << endl;
			o<<tab<< "             R => fracResultfar0, " << endl; 
			o<<tab<< "             clk => clk, " << endl;
			o<<tab<< "             rst => rst " << endl;
			o<<tab<< "            );" << endl<<endl;
				
		// if the second operand was shifted out of the operation, then the result of the operation becomes = to the first operand //		
		o<<tab<< "fracResultfar0wSh <= fracResultfar0_d when "<< get_delay_signal_name("shiftedOut",rightShifter->pipeline_depth()+intaddFar1->pipeline_depth()+2)<<"='0' else "
		                            <<"(\"01\" & "<<get_delay_signal_name("fracFX",rightShifter->pipeline_depth()+intaddFar1->pipeline_depth()+2)<<"&\"000\");"<<endl;
		

		// the result exponent before normalization and rounding is = to the exponent of the first operand //
		o<<tab<<"exponentResultfar0<=\"0\" & "<<get_delay_signal_name("expFX",rightShifter->pipeline_depth()+intaddFar1->pipeline_depth()+3)<<";"<<endl;
		
		
		//perform NORMALIZATION and recalculation of the stiky bit 
		//delay the possible exponent adjustment until the rounding phase
	  o<<tab<<"fracResultfar1<=fracResultfar0wSh("<<wF+3<<" downto 2) & "
	  											<<"(fracResultfar0wSh(1) or fracResultfar0wSh(0)) when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"01\" else"<<endl;
    o<<tab<< "    					 fracResultfar0wSh("<<wF+2<<" downto 0) 	 when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"00\" else"<<endl;
    o<<tab<< "    					 fracResultfar0wSh("<<wF+4<<" downto 3) & (fracResultfar0wSh(2) or fracResultfar0wSh(1) or fracResultfar0wSh(0));"<<endl;
					
		
		//select operation mode. This depends on wether or not the exponent must be adjusted after normalization		
		o<<tab<<"expOperationSel <= \"00\" when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"01\" else"<<endl;
		o<<tab<<"                   \"11\" when fracResultfar0wSh("<<wF+4<<" downto "<<wF+3<<") = \"00\" else"<<endl;
		o<<tab<<"                   \"01\";"<<endl;
		
		//the second operand depends on the operation selector
		o<<tab<<"expSecOperand <= (("<<wE<<" downto 1 => '0') & expOperationSel(0)) xor ("<<wE<<" downto 0 => expOperationSel(1));"<<endl;
		
		// readjust exponent after normalization //
		o<<tab<< "int_adder_componentf2: " << intaddFar2->getOperatorName() << endl;
		o<<tab<< "  port map ( X => exponentResultfar0, " << endl; 
		o<<tab<< "             Y => expSecOperand_d, " << endl; 
		o<<tab<< "             Cin => expOperationSel_d(1) ," << endl;
		o<<tab<< "             R => exponentResultfar1, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;
		
		
		//ROUNDING
		o<<tab<< "round <= fracResultfar1(1) and (fracResultfar1(2) or fracResultfar1(0));"<<endl;

		//concatenation of the exponent and the fraction. 	
		o<<tab<<"expConcatFrac <= exponentResultfar1_d & "<<get_delay_signal_name("fracResultfar1",intaddFar2->pipeline_depth()+2)<<"("<<wF+2<<" downto 2);"<<endl;
		
		// round //
		o<<tab<< "int_adder_componentf3: " << intaddFar3->getOperatorName() << endl;
		o<<tab<< "  port map ( X => expConcatFrac, " << endl; 
		o<<tab<< "             Y => (others => '0'), " << endl; 
		o<<tab<< "             Cin => "<<get_delay_signal_name("round",intaddFar2->pipeline_depth()+2) <<" ," << endl;
		o<<tab<< "             R => expConcatFracResult, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;
		
		o<<tab<< "exponentResultfar <= expConcatFracResult_d("<<wE+1+wF<<" downto "<<wF+1<<");"<<endl;
		o<<tab<< "fractionResultfar <= expConcatFracResult_d("<<wF<<" downto "<<0<<");"<<endl;
				
  	//==========================================================================
		//==========================================================================
		
		
		
		//=========================================================================|
		//                              Synchronization                            |
		//=========================================================================|				
		o<<"-- SYNCH --"<<endl;

		// Synchronization Interface Output
				
		//synchronize the close signal
		o<<tab<<"pipeClose <= sdClose;"<<endl;
		o<<tab<<"syncClose <= "<<get_delay_signal_name("pipeClose", maxPathDepth)<<";"<<endl; 
				
		//exception bits
		o<<tab<<"pipeXAB <= sdxXY;"<<endl;
		o<<tab<<"syncXAB <= "<<get_delay_signal_name("pipeXAB",maxPathDepth)<<";"<<endl;
				
		//sign bit
		o<<tab<<"pipeSignAB <= sdSignAB;"<<endl;
		o<<tab<<"syncSignAB <= "<<get_delay_signal_name("pipeSignAB",maxPathDepth)<<";"<<endl;
				
		//X
		o<<tab<<"pipeX <= sdX ;"<<endl;
		o<<tab<<"syncX <= "<<get_delay_signal_name("pipeX",maxPathDepth)<<";"<<endl;
		
		//signY
		o<<tab<<"pipeSignY <= sdSignY ;"<<endl;
		o<<tab<<"syncSignY <= "<<get_delay_signal_name("pipeSignY",maxPathDepth)<<";"<<endl;
		
		if (closePathDepth == farPathDepth){
		
			o<<tab<<"expResClose <= exponentResultCloseC;"<<endl;
			o<<tab<<"expResFar <= exponentResultfar;"<<endl;
			
			o<<tab<<"fracResClose <= fractionResultCloseC;"<<endl;
			o<<tab<<"fracResFar <= fractionResultfar;"<<endl;
			
			o<<tab<<"syncRS <= crsOut;"<<endl;
		}
		else
			if (closePathDepth > farPathDepth){	
				o<<tab<<"expResClose <= exponentResultCloseC;"<<endl;
				
				o<<tab<<"pipeExpResFar <= exponentResultfar;"<<endl;
				o<<tab<<"expResFar <= "<<get_delay_signal_name("pipeExpResFar",closePathDepth-farPathDepth) <<";"<<endl;
				
				o<<tab<<"fracResClose <= fractionResultCloseC;"<<endl;
				
				o<<tab<<"pipeFracResFar <= fractionResultfar;"<<endl;
				o<<tab<<"fracResFar <= "<<get_delay_signal_name("pipeFracResFar",closePathDepth-farPathDepth) <<";"<<endl;
				
				o<<tab<<"syncRS <= crsOut;"<<endl;
			}else{
				o<<tab<<"pipeExpResClose <= exponentResultCloseC;"<<endl;
				o<<tab<<"expResClose <= "<<get_delay_signal_name("pipeExpResClose",farPathDepth-closePathDepth) <<";"<<endl;
								
				o<<tab<<"expResFar <= exponentResultfar;"<<endl;
				
				o<<tab<<"pipeFracResClose <= fractionResultCloseC;"<<endl;
				o<<tab<<"fracResClose <= "<<get_delay_signal_name("pipeFracResClose",farPathDepth-closePathDepth) <<";"<<endl;
				
				o<<tab<<"fracResFar <= fractionResultfar;"<<endl;
				
				o<<tab<<"pipeSyncRS <= crsOut;"<<endl;
				o<<tab<<"syncRS <= "<<get_delay_signal_name("pipeSyncRS",farPathDepth-closePathDepth) <<";"<<endl;
			
			
			}
			
		
		
		
		//=========================================================================|
		//                     final result composition                            |
		//=========================================================================|				
		o<<"-- Roundings --"<<endl;
				
		// select between the results of the close or far path as the result of the operation 
		o<<tab<< "with syncClose select"<<endl;
    o<<tab<< "exponentResultn <= expResClose when '1',"<<endl;
    o<<tab<< "                   (\"0\" & expResFar) when others;"<<endl;
		o<<tab<< "with syncClose select"<<endl;
		o<<tab<< "fractionResultn <= fracResClose when '1',"<<endl;
		o<<tab<< "                   fracResFar when others;"<<endl;
				
		// compute the exception bits of the result considering the possible underflow and overflow 
		o<<tab<< "eTest <= exponentResultn("<<wE+1<<" downto "<<wE<<");"<<endl;
		
		o<<tab<< "with eTest select"<<endl;
    o<<tab<< "nRn("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= \"10\" when \"01\","<<endl;
    o<<tab<< "                               \"00\" when \"10\" | \"11\","<<endl;
    o<<tab<< "                               \"01\" when others;"<<endl;
  	
  	o<<tab<< "nRn("<<wE+wF<<" downto 0) <= syncRS & exponentResultn("<<wE-1<<" downto 0) & fractionResultn("<<wF-1<<" downto 0);"<<endl;
	 
	  o<<tab<< "with syncXAB_d select"<<endl;
	  o<<tab<< "  nnR("<<wE+wF+2<<" downto "<<wE+wF+1<<") <= nRn_d("<<wE+wF+2<<" downto "<<wE+wF+1<<") when \"0101\","<<endl;
	  o<<tab<< "                                 \"1\" & syncSignAB_d              when \"1010\","<<endl;
	  o<<tab<< "                                 \"11\"            		           when \"1011\","<<endl;
	  o<<tab<< "                                 syncXAB_d(3 downto 2)             when others;"<<endl;
	  o<<tab<< "with syncXAB_d select"<<endl;
	  o<<tab<< "  nnR("<<wE+wF<<") <= nRn_d("<<wE+wF<<")             when \"0101\","<<endl;
	  o<<tab<< "                syncX_d("<<wE+wF<<") and syncSignY_d when \"0000\","<<endl;
	  o<<tab<< "                syncX_d("<<wE+wF<<")                 when others;"<<endl;

	  o<<tab<< "with syncXAB_d select"<<endl;
	  o<<tab<< "  nnR("<<wE+wF-1<<" downto 0) <= nRn_d("<<wE+wF-1<<" downto 0) when \"0101\","<<endl;
	  o<<tab<< "                                 syncX_d("<<wE+wF-1<<" downto 0) when others;"<<endl;
			
		// assign results 
		o<<tab<< "rexp <= nnR("<<wE+wF-1<<" downto "<<wF<<");"<<endl;
		o<<tab<< "rfrac <= \"1\" & nnR("<<wF-1<<" downto 0);"<<endl;
		o<<tab<< "rsig <= nnR("<<wE+wF<<");"<<endl;
		o<<tab<< "rexc <= nnR("<<wE+wF+2<<" downto "<<wE+wF+1<<");"<<endl;
					
		











				   
	} 
	else
	{  
		//=========================================================================| 
		//                                                                         |
		//                             COMBINATORIAL                               |
		//                                                                         |
		//=========================================================================|
		
	  /* == Swap/Difference ==*/
	  /* signal which indicates whether or not the exception bits of X are greater or equal than the exception bits of Y */		  
    o<<tab<<"exceptionXSuperiorY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") >= Y("<<wEY+wFY+2<<" downto "<<wEY+wF+1<<") else"<<endl;
		o<<tab<<"                       '0';"<<endl;
			
		/* signal which indicates whether or not the exception bits of X are equal to the exception bits of Y */		  
		o<<tab<<"exceptionXEqualY <= '1' when X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<") = Y("<<wEY+wFY+2<<" downto "<<wEY+wFY+1<<") else"<<endl;
		o<<tab<<"                    '0';"<<endl;

		/* make the difference between the exponents of X and Y. Pad exponents with sign bit before prforming the substraction */
		o<<tab<<"exponentDifference0 <= (\"0\" & X("<<wEX+wFX-1<<" downto "<<wFX<<")) - (\"0\" & Y("<<wEY+wFY-1<<" downto "<<wFY<<"));"<<endl;

		/* swapping is performed when:    the exception bits of X and Y are equal and the exponent of Y is greater than the exponent of X (i.e. sign bit of exponent difference=1)
		                              or  the exception bits of Y are superior to the exception bits of X (EX: exc. bits of X are 00 and of Y are 01) */
		o<<tab<<"swap <= (exponentDifference0("<<wER<<") and exceptionXEqualY) or (not exceptionXSuperiorY);"<<endl;

		/* depending on the value of swap, assign the corresponding values to the newX and newY signals */
		o<<tab<<"newX <= Y when swap = '1' else"<<endl;
		o<<tab<<"        X;"<<endl;   
		o<<tab<<"newY <= X when swap = '1' else"<<endl;
    o<<tab<<"        Y;"<<endl;

		/* readjust for the exponents difference after the potential swap */
		o<<tab<<"exponentDifference1 <= exponentDifference0("<<wER-1<<" downto "<<0<<") xor ("<<wER-1<<" downto "<<0<<" => swap);"<<endl;
		o<<tab<<"zeroExtendedSwap    <= "<< zero_generator(wE-1,0)<<" & swap;"<<endl;
		o<<tab<<"exponentDifference  <= exponentDifference1  + zeroExtendedSwap;"<<endl;
		
			
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
		o<<tab<< "LZC_component: " << leadingZeroCounter->getOperatorName() << endl;
		o<<tab<< "      port map ( I => fracRClose1, " << endl; 
		o<<tab<< "                 OZB => '0', " << endl; 
		o<<tab<< "                 O => nZeros " <<endl; 
		o<<tab<< "               );" << endl<<endl;		
				
		/* shift and round */
		o<<tab<< "left_shifter_component: " << leftShifter->getOperatorName() << endl;
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
		o<<tab<< "right_shifter_component: " << rightShifter->getOperatorName() << endl;
		o<<tab<< "      port map ( X => fracNewY, " << endl; 
		o<<tab<< "                 S => exponentDifference("<< sizeRightShift-1<<" downto 0"<<"), " << endl; 
		o<<tab<< "                 R => shiftedFracY " <<endl; 
		o<<tab<< "               );" << endl<<endl;		
		
		/* determine if the fractional part of Y was shifted out of the opperation */
		//TODO
		if (wE-1>=sizeRightShift)
		{
			o<<tab<<" shiftedOut <= "; 
			for (int i=wE-1;i>=sizeRightShift;i--)
				if (((wE-1)==sizeRightShift)||(i==sizeRightShift))
					o<< "exponentDifference("<<i<<")";
				else
					o<< "exponentDifference("<<i<<") or ";
			o<<";"<<endl;
		}
		else
		o<<tab<<" shiftedOut <= '0'"; 
		
		
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

TestIOMap FPAdder::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("X"));
	tim.add(*get_signal_by_name("Y"));
	tim.add(*get_signal_by_name("rexc"));
	tim.add(*get_signal_by_name("rsig"));
	tim.add(*get_signal_by_name("rexp"));
	tim.add(*get_signal_by_name("rfrac"));
	return tim;
}

void FPAdder::fillTestCase(mpz_class a[])
{
	/* Get I/O values */
	mpz_class& svX = a[0];
	mpz_class& svY = a[1];
	mpz_class& svExc = a[2];
	mpz_class& svSig = a[3];
	mpz_class& svExp = a[4];
	mpz_class& svFra = a[5];

	/* Compute correct value */
	FPNumber x(wEX, wFX), y(wEY, wFY), r(wER, wFR);
	x = svX;
	y = svY;
	r = x+y;
	
	/* Set outputs */
	/* XXX: This is old-school; TestBenches are able
	 * to compare whole FP signals */
	svExc = r.getExceptionSignalValue();
	if (svExc < 3)
	{
		svSig = r.getSignSignalValue();
		if (svExc == 1)
		{
			svExp = r.getExponentSignalValue();
			svFra = r.getFractionSignalValue();
		}
	}
}

