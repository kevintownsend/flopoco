/*
 * Floating Point Multiplier for FloPoCo
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

#include "FPMultiplier.hpp"

using namespace std;
extern vector<Operator*> oplist;


/**
 * The FPMutliplier constructor
 * @param[in]		target		the target device
 * @param[in]		wEX			the the with of the exponent for the number f-p number X
 * @param[in]		wFX			the the with of the fraction for the number f-p number X
 * @param[in]		wEY			the the with of the exponent for the number f-p number Y
 * @param[in]		wFY			the the with of the fraction for the number f-p number Y
 * @param[in]		wER			the the with of the exponent for the number multiplication result
 * @param[in]		wFR			the the with of the fraction for the number multiplication result
 **/
FPMultiplier::FPMultiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm) :
	Operator(target), wEX(wEX), wFX(wFX), wEY(wEY), wFY(wFY), wER(wER), wFR(wFR) {

	int i, j;
	ostringstream name, synch, synch2;

	//set the normalized attribute  
	if (norm == 0) 
		normalized = false;
	else
		normalized = true;


	//parameters setup
	addition_chunk_width = 17;
	
	



	/** Name Setup procedure
	 * The name has the format: FPMultiplier_wEX_wFX_wEY_wFY_wER_wFR
	 * wEX = width of X exponenet
	 * wFX = width for the fractional part of X
	 **/
	name.str("");
	name <<"FPMultiplier_"; 
	name<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR;
	unique_name = name.str(); //set the attribute
	

	// Set up the IO signals
	add_input ("X", wEX + wFX + 3);// format: 2b(Exception) + 1b(Sign)
	add_input ("Y", wEY + wFY + 3);//+ wEX bits (Exponent) + wFX bits(Fraction)

	add_output ("ResultExponent"   , wER    );  
	add_output ("ResultSignificand", wFR + 1);
	add_output ("ResultException"  , 2      );
	add_output ("ResultSign"       , 1      ); 

	//Set up if we are in the sequential or combinational case
	if (target->is_pipelined()) 
		set_sequential();
	else
		set_combinatorial(); 
	   
	//Instantiate an IntMultiplier -> multiply the significands
	intmult = new IntMultiplier(target, wFX+1, wFY+1);
	oplist.push_back(intmult);

	//Setup the attibute cocerning the IntMultiplier pipeline 
	IntMultPipelineDepth = intmult->pipeline_depth();


	// Unregistered signals
	if ((normalized)||!(is_sequential()) ) {
		if (!(1 + wFR >= wFX + wFY + 2)) {
			add_signal("check_string",wFX+wFY+2 - (1+wFR) );
			add_signal("reunion_signal",2+wEX+wFR);
			add_signal("reunion_signal_out",2+wEX+wFR);
			add_signal("significand_synch_tb2",1); 
			add_signal("middle_signal",2+wEX+wFR);
			add_signal("reunion_signal_p_upd",2+wEX+wFR); 
		}
		add_signal("mult_selector",1);
		add_signal("exponent_p",1+wEX);
		add_signal("exception_upd",2);
	}	
	add_signal("significandX", wFX+1);   //The Significands for the input signals X and Y
	add_signal("significandY", wFY+1);
	add_signal("significand_product", wFX +wFY +2);   

	add_signal("exponentX", wEX);	//The signals for the exponents of X and Y
	add_signal("exponentY", wEY);
	add_signal("bias",wEX+1); //bias to be substracted from the sum of exponents  
	add_signal("exponents_sum_minus_bias_ext",wEX+1); //extended with 1 bit

	add_signal("exception_selector",4);//signal that selects the case for the exception

	if (is_sequential())  {
		// Unregistered signals
		add_signal("exponent_synch", wEX);
		add_signal("significand_synch", wFX + wFY + 2); 
		add_signal("sign_synch", 1); 
		add_signal("exception_synch",2);
		add_signal("significand_synch_tb2_out",1);	
		add_signal("middle_output_out",1);
		add_signal("sign_synch2_out",1);
		add_signal("exception_synch2_out",2);
		
		if (!normalized)
			add_signal("temp_exp_concat_fract", 1 + wEX + wFX + wFY + 2);
		
		// Registered Signals
		if (normalized){
			if (!(1+wFR >= wFX+wFY+2)) {
				add_signal("reunion_signal_p",2+wEX+wFR);   
				add_registered_signal_with_sync_reset("middle_output",1);
				add_registered_signal_with_sync_reset("exponent_synch2",wEX);
				add_registered_signal_with_sync_reset("exception_synch2",2);
				add_registered_signal_with_sync_reset("significand_synch2",wFX + wFY + 2);
				add_registered_signal_with_sync_reset("sign_synch2",1);
				add_registered_signal_with_sync_reset("mult_selector2",1);
			
				reunion_signal_width = 2+ wEX + wFR;
				reunion_signal_parts = int ( ceil( double(reunion_signal_width)/double(addition_chunk_width)));
				addition_last_chunk_width = reunion_signal_width - (reunion_signal_parts - 1)*addition_chunk_width;
							
				for (j=1; j<=reunion_signal_parts;j++)	
					for (i=1;i<=reunion_signal_parts;i++){	
						name.str("");
						name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
						if (i!=reunion_signal_parts)
							add_registered_signal_with_sync_reset(name.str(), addition_chunk_width + 1 );
						else
							add_registered_signal_with_sync_reset(name.str(), addition_last_chunk_width );	
	         }
				
				if (reunion_signal_parts>1){
					add_delay_signal("reunion_signal_level",2+wEX+wFR,reunion_signal_parts);
					add_delay_signal("significand_synch_tb2_level",1,reunion_signal_parts);		
					add_delay_signal("middle_output_level",1,reunion_signal_parts);			
					add_delay_signal("sign_synch2_level",1,reunion_signal_parts);
					add_delay_signal("exception_synch2_level",2,reunion_signal_parts);
				}
			}  
		}

		add_registered_signal_with_sync_reset("Exponents_Sum", wEX+1 );
		add_registered_signal_with_sync_reset("Exponents_Sum_Minus_Bias", wEX );

		//registers for the case when the integer multiplication requires 
		//more than two pipeline levels
		for (i=0;i<=IntMultPipelineDepth-3;i++)	{
			synch.str(""); synch2.str("");
			synch<<"Exponent_Synchronization_"<<i;
			add_registered_signal_with_sync_reset(synch.str(), wEX );
			synch2<<"Exception_Synchronization_"<<i;
			add_registered_signal_with_sync_reset(synch2.str(), 2 );	
		}

		//registers for the case when the int mult has less than two levels
		for (i=0;i<=1-IntMultPipelineDepth;i++)	{
			synch.str("");
			synch<<"Int_Synchronization_"<<i;
			add_registered_signal_with_sync_reset(synch.str(), wFX+wFY+2);		
		}	

		add_delay_signal("Sign_Synchronization", 1, max(2,IntMultPipelineDepth));

		add_registered_signal_with_sync_reset("Result_Exception",2); 
		add_registered_signal_with_sync_reset("Result_Exception_With_EA",2);
	}//end if sequential
	else
	{//if combinational
		add_signal("middle_output",1);
		add_signal("exponent_p2",wEX);
		add_signal("Exponents_Sum", wEX+1 );
		add_signal("Exponents_Sum_Minus_Bias", wEX );   
		add_signal("Result_Exception",2); 
		add_signal("Result_Exception_With_EA",2);
	}

	//set pipeline depth 
	if (!is_sequential())
		set_pipeline_depth(0);
	else
		if (normalized)
			if (1+wFR >= wFX+wFY+2)   
				set_pipeline_depth(max(2,IntMultPipelineDepth));
			else
				if (reunion_signal_parts>1)
					set_pipeline_depth(1 + max(2,IntMultPipelineDepth) + reunion_signal_parts );
				else
					set_pipeline_depth(1 + max(2,IntMultPipelineDepth));
		else
			set_pipeline_depth(max(2,IntMultPipelineDepth));


}


FPMultiplier::~FPMultiplier() {
}


void FPMultiplier::output_vhdl(std::ostream& o, std::string name) {
  
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;

	int bias_val=int(pow(double(2),double(wEX-1)))-1;
	int i, j; 

	Licence(o,"Bogdan Pasca (2008)");
	Operator::StdLibs(o);
	o<<endl<<endl;	
	output_vhdl_entity(o);

	new_architecture(o,name);
	
	intmult->output_vhdl_component(o);
	output_vhdl_signal_declarations(o);	  
	o<<endl;

	begin_architecture(o);
	
	if (is_sequential()){
		//common code for both normalized and non-normalized version
		output_vhdl_registers(o); o<<endl;
				
		o<<tab<<"exponentX <= X("<<wEX + wFX -1<<" downto "<<wFX<<");"<<endl; 
		o<<tab<<"exponentY <= Y("<<wEY + wFY -1<<" downto "<<wFY<<");"<<endl<<endl;

		//Add exponents -> put exponents sum into register
		o<<tab<<"Exponents_Sum <= (\"0\" & exponentX) + (\"0\" & exponentY);"<<endl; //wEX+1 bits
		o<<tab<<"bias <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wEX+1<<");"<<endl; 

		//substract the bias value from the exponent sum
		o<<tab<<"exponents_sum_minus_bias_ext <= Exponents_Sum_d - bias;"<<endl;//wEX + 1   
		o<<tab<<"Exponents_Sum_Minus_Bias <= exponents_sum_minus_bias_ext("<<wEX-1<<" downto 0);"<<endl; //wEX

		//build significands 	 	
		o<<tab<< "significandX <= \"1\" & X("<<(wFX-1)<<" downto 0);"<<endl;
		o<<tab<< "significandY <= \"1\" & Y("<<(wFY-1)<<" downto 0);"<<endl<<endl;

		//multiply significands 
		o<<tab<< "int_multiplier_component: " << intmult->unique_name << endl;
		o<<tab<< "      port map ( X => significandX, " << endl; //wFX+1 bits (1 bit is the hidden 1)
		o<<tab<< "                 Y => significandY, " << endl; //wFY+1 bits
		o<<tab<< "                 R => significand_product, " << endl; //wFX+wFY+2 bits
		o<<tab<< "                 clk => clk, " << endl;
		o<<tab<< "                 rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;

		//Concatenate the exception bits of the two input numbers
		o<<tab<<"exception_selector <= X("<<wEX + wFX +2<<" downto "<<wEX + wFX + 1<<") & Y("<<wEY + wFY +2<<" downto "<<wEY + wFY +1<<");"<<endl<<endl;

		//select proper exception with selector			
		o<<tab<<"process (exception_selector)"<<endl;
		o<<tab<<"begin"<<endl;
		o<<tab<<"case (exception_selector) is"<< endl;
		o<<tab<<tab<<"when \"0000\" => Result_Exception <= \"00\";   "<<endl;
		o<<tab<<tab<<"when \"0001\" => Result_Exception <= \"00\";   "<<endl;
		o<<tab<<tab<<"when \"0010\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"0011\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"0100\" => Result_Exception <= \"00\";   "<<endl;
		o<<tab<<tab<<"when \"0101\" => Result_Exception <= \"01\";   "<<endl;//to see
		o<<tab<<tab<<"when \"0110\" => Result_Exception <= \"10\";   "<<endl;
		o<<tab<<tab<<"when \"0111\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1000\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1001\" => Result_Exception <= \"10\";   "<<endl;
		o<<tab<<tab<<"when \"1010\" => Result_Exception <= \"10\";   "<<endl;
		o<<tab<<tab<<"when \"1011\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1100\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1101\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1110\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1111\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when others =>   Result_Exception <= \"11\";   "<<endl;    
		o<<tab<<"end case;"<<endl;
		o<<tab<<"end process;"<<endl<<endl;

		//use MSB of exponents_sum_minus_bias_ext to siagnal overflow		
		o<<tab<<"Result_Exception_With_EA <= Result_Exception_d when exponents_sum_minus_bias_ext("<< wEX<<")='0' else"<<endl;
		o<<tab<<"                            \"11\";"<<endl<<endl;
									
		// synchronization barrier (IntMultiplication | Exponent Addition | Exception Computation | Sign Computation
		// 0. when IntMultiplication | Exponent Addition are allready synchronized	
		// 1. when registers are added to exponent addition
		// 2. when registers are added to integer multiplication
	  
		if (IntMultPipelineDepth==2){ 
			o<<tab<<"exponent_synch    <= Exponents_Sum_Minus_Bias_d;"<<endl;
			o<<tab<<"significand_synch <= significand_product;"<<endl;
			o<<tab<<"exception_synch   <= Result_Exception_With_EA_d;"<<endl;
		}
	   else{	
			if (IntMultPipelineDepth > 2){
				//connect the first synchronization reg to the output of the last reg from exponent addition 
				o<<tab<<"Exponent_Synchronization_0  <= Exponents_Sum_Minus_Bias_d;"<<endl;
				o<<tab<<"Exception_Synchronization_0 <= Result_Exception_With_EA_d;"<<endl; 									

				for (i=1;i<=IntMultPipelineDepth-3;i++){
					o<<tab<<"Exponent_Synchronization_"<<i<<"   <= "<<"Exponent_Synchronization_"<<i-1<<"_d;"<<endl;		
					o<<tab<<"Exception_Synchronization_"<<i<<" <= "<<"Exception_Synchronization_"<<i-1<<"_d;"<<endl;
				}
					o<<tab<<"exponent_synch    <= Exponent_Synchronization_"<<	IntMultPipelineDepth-3 <<"_d;"<<endl;
					o<<tab<<"exception_synch   <= Exception_Synchronization_"<<IntMultPipelineDepth-3 <<"_d;"<<endl;
					o<<tab<<"significand_synch <= significand_product;"<<endl;			
			}
			else{
				//add registers to the output of IntMultiplier so that exponent addition and 
				//significand multiplication become synchronized
				o<<tab<<"Int_Synchronization_0 <= significand_product;"<<endl;
				for (i=1;i<=1-IntMultPipelineDepth;i++){
					synch1.str(""); synch2.str("");
					synch2<<"Int_Synchronization_"<<i;
					synch1<<"Int_Synchronization_"<<i-1<<"_d";
					o<<tab<<synch2.str()<<"<="<<synch1.str()<<";"<<endl;		
				}
				o<<tab<<"exponent_synch    <= Exponents_Sum_Minus_Bias_d;"<<endl;
				o<<tab<<"exception_synch   <= Result_Exception_With_EA_d;"<<endl;
				o<<tab<<"significand_synch <= Int_Synchronization_"<< 1-IntMultPipelineDepth<<"_d;"<<endl;
			}
		}
		
		//sign synchronization
		o<<tab<<"Sign_Synchronization <= X("<<wEX + wFX<<") xor Y("<<wEY + wFY<<");"<<endl;
		o<<tab<<"sign_synch <= "<<get_delay_signal_name("Sign_Synchronization",max(2,IntMultPipelineDepth))<<";"<<endl<<endl;
	

		if (normalized==true){	
			//normalize result
			//assign selector signal to MSB of signif. multiplic 
			o<<tab<<"mult_selector <= significand_synch ("<<wFX+wFY+1<<");"<<endl;

			o<<tab<<"exponent_p <= (\"0\" & exponent_synch ) + (CONV_STD_LOGIC_VECTOR(0 ,"<<wEX<<") & mult_selector);"<<endl;

			//signal the exception
			o<<tab<<"exception_upd <= exception_synch when exponent_p("<<wEX<<")='0' else"<<endl;
			o<<tab<<"                \"11\";"<<endl;

			//round result
			//check is rounding is needed
			if (1+wFR >= wFX+wFY+2) {
				//=>no rounding needed - possible padding
				//generate zeros of length 1+wFR - wFX+wFY+2 
				zeros1.str("");
				zeros1 << zero_generator(1+wFR - (wFX+wFY+2) , 0);

				//make the outputs
				o<<tab<<"ResultExponent <= exponent_p("<<wEX - 1<<" downto 0);"<<endl;  
				o<<tab<<"ResultSignificand <= significand_synch & "<<zeros1.str()<<" when mult_selector='1' else"<<endl;
				o<<tab<<"                     significand_synch("<<wFX+wFY<<" downto 0) & "<<zeros1.str()<<" & \"0\";"<<endl;
				o<<tab<<"ResultException <= exception_upd;"<<endl;
				o<<tab<<"ResultSign <= sign_synch;"<<endl;

			}
			else{
				//check if in the middle of two FP numbers
				//generate two xor strings
				str1.str(""); str2.str(""); //init
				str1<<"\"1"<< zero_generator( wFX+wFY+2 - (1+wFR) -1, 1);
				str2<<"\"01"<< zero_generator( wFX+wFY+2 - (1+wFR) -1-1, 1);

				zeros1.str("");
				zeros1<< zero_generator( wFX+wFY+2 - (1+wFR) , 0);

				o<<tab<<"check_string <= (significand_synch ("<< wFX+wFY+2 - (1+wFR) -1<<" downto 0) xor "<<str1.str()<<") when mult_selector='1' else"<<endl;
				o<<tab<<"               ( (\"0\" & significand_synch ("<< wFX+wFY+2 - (1+wFR) -1-1<<" downto 0)) xor "<<str2.str()<<");"<<endl;

				o<<tab<<"process(clk)"<<endl;
				o<<tab<<"begin"<<endl;
				o<<tab<<tab<<"   if (clk'event and clk='1') then "<<endl;  
				o<<tab<<tab<<"      if ( "<<zeros1.str()<<" = check_string ) then"<<endl; 
				o<<tab<<tab<<"         middle_output <= '1';"<<endl;
				o<<tab<<tab<<"      else "<<endl;
				o<<tab<<tab<<"         middle_output <= '0';"<<endl;
				o<<tab<<tab<<"      end if;"<<endl;
				o<<tab<<tab<<"   end if; "<<endl;
				o<<tab<<tab<<"end process; "<<endl;

				//propagate rest of signals 1 level
				o<<tab<<"exponent_synch2 <= exponent_p("<<wEX-1<<" downto 0);"<<endl;
				o<<tab<<"exception_synch2 <= exception_upd;"<<endl;
				o<<tab<<"significand_synch2 <= significand_synch;"<<endl;
				o<<tab<<"sign_synch2 <= sign_synch;"<<endl;
				o<<tab<<"mult_selector2 <= mult_selector;"<<endl;

				o<<tab<<"reunion_signal <= (\"0\" & exponent_synch2_d & significand_synch2_d("<< wFX+wFY<<" downto "<<wFX+wFY - (wFR) <<" )) when mult_selector2_d='1' else"<<endl;
				o<<tab<<"                 (\"0\" & exponent_synch2_d & significand_synch2_d("<< wFX+wFY-1<<" downto "<<wFX+wFY - (wFR) -1<<" ));"<<endl; 

				
				               
				o<<tab<<"significand_synch_tb2 <= significand_synch2_d("<< wFX+wFY+1 - wFR<<" ) when mult_selector2_d='1' else"<<endl;
				o<<tab<<"                         significand_synch2_d("<< wFX+wFY+1 - wFR -1<<" );"<<endl;              
				
				if (reunion_signal_parts != 1){
					for (j = 1; j <= reunion_signal_parts; j++){
						if (j == 1)
							o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <= (\"0\" & reunion_signal("<<j*addition_chunk_width-1<<" downto "<<(j-1)*addition_chunk_width<<")) + '1';"<<endl;
						else
							if (j!=reunion_signal_parts)
								o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <= (\"0\" & reunion_signal("<<j*addition_chunk_width-1<<" downto "<<(j-1)*addition_chunk_width<<"));"<<endl;	
							else
								o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <=  reunion_signal("<<(j-1)*addition_chunk_width + addition_last_chunk_width -1<<" downto "<<(j-1)*addition_chunk_width<<");"<<endl;	
					}
					
					//put tohether the pipeline
				for (j=1; j<=reunion_signal_parts-1;j++)	
					for (i=1;i<=reunion_signal_parts;i++)
						if ((i<=j)||(i>j+1))
							o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
						else
							o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<addition_chunk_width<<") ;"<<endl; 

					//put the result back together
					zeros.str("");
					for (i = reunion_signal_parts; i > 0; i--) 
						if ( (i!=1) && (i!=reunion_signal_parts) )
							zeros<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0) &" ;
						else
							if (i==reunion_signal_parts)
								zeros<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_last_chunk_width -1<<" downto 0) & " ;
							else
								zeros<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0)" ;
					
				
					o<<tab<<"reunion_signal_level <=reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_out <= "<<get_delay_signal_name("reunion_signal_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"significand_synch_tb2_level <= significand_synch_tb2;"<<endl;
					o<<tab<<"significand_synch_tb2_out <= "<<get_delay_signal_name("significand_synch_tb2_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"middle_output_level <= middle_output_d;"<<endl;
					o<<tab<<"middle_output_out <= "<<get_delay_signal_name("middle_output_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"sign_synch2_level <= sign_synch2_d;"<<endl;
					o<<tab<<"sign_synch2_out <= "<<get_delay_signal_name("sign_synch2_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"exception_synch2_level <= exception_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out <= "<<get_delay_signal_name("exception_synch2_level", reunion_signal_parts)<<";" <<endl;		
														
					o<<tab<<"reunion_signal_p <= "<< zeros.str()<<";"<< endl;
					
				
				}   
				else{
					//add 1 to the reunited signal for rounding & normalize purposes
					o<<tab<<"reunion_signal_out <= reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_p <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2+ wEX + wFR<<");"<<endl;
					o<<tab<<"middle_output_out <= middle_output_d;"<<endl;
					o<<tab<<"significand_synch_tb2_out <= significand_synch_tb2;"<<endl;
					o<<tab<<"sign_synch2_out <= sign_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out <= exception_synch2_d"<<endl;
					
				}				
				   		     
				   		     
				o<<tab<<"middle_signal <= reunion_signal_p when significand_synch_tb2_out = '1' else"<<endl;
				o<<tab<<"                reunion_signal_out;"<<endl;    		     
				   		     
				o<<tab<<"reunion_signal_p_upd <= reunion_signal_p when middle_output_out='0' else"<<endl;
				o<<tab<<"                       middle_signal;"<<endl;
				                        
				//update exception                       
				o<<tab<<"ResultException <= exception_synch2_out when reunion_signal_p_upd("<<wEX + wFR+1<<")='0' else"<<endl;
				o<<tab<<"                   \"11\";"<<endl;                                                   

				o<<tab<<"ResultExponent <= reunion_signal_p_upd("<<wEX + wFR<<" downto "<<wFR+1<<");"<<endl;  
				o<<tab<<"ResultSignificand <= \"1\" & reunion_signal_p_upd("<<wFR<<" downto 1);"<<endl;
				o<<tab<<"ResultSign <= sign_synch2_out;"<<endl;
			}
		}
		else { 
			if (1+wFR < wFX+wFY+2){
				//do rounding
				zeros1.str("");	zeros2.str("");		
				zeros1 << zero_generator(wEX + wFR + 2, -1)<<"1"<<zero_generator(wFX+wFY+2-(wFR+1)-1, 1);
				zeros2 << zero_generator(wEX + wFR + 3, -1)<<"1"<<zero_generator(wFX+wFY+2-(wFR+1)-2, 1);  

				o<<tab<<"temp_exp_concat_fract <= (\"0\" & exponent_synch & significand_synch) + "<< zeros1.str()<<" when significand_synch("<<wEX+wEY+1<<")='1' else	"<<endl;
				o<<tab<<"                           (\"0\" & exponent_synch & significand_synch) + "<<zeros2.str()<<";"<<endl;	

				o<<tab<<"ResultSignificand <= temp_exp_concat_fract ("<<wFX+wFY+1<<" downto "<< wFX+wFY+1 - (wFR+1)+1<<");"<<endl;
				o<<tab<<"ResultExponent    <= temp_exp_concat_fract ("<< wEX+wFX+wFY + 1 <<" downto  "<< wFX+wFY+2 <<");"<<endl; 
				o<<tab<<"ResultException   <= exception_synch when temp_exp_concat_fract("<<wEX+wFX+wFY + 2 <<")='0' else"<<endl;
				o<<tab<<"		                 \"11\";"<<endl;			
			}
			else{
				zeros.str("");//initialization
				zeros<< zero_generator((wFR+1) - (wFX+wFY+2), 0 );

			if (wFX+wFY+2 == wFR+1)
				o<<tab<<"ResultSignificand <= significand_synch;"<<endl;
			else 	
				o<<tab<<"ResultSignificand <= significand_synch & "<<zeros.str()<<";"<<endl;
				
			o<<tab<<"ResultExponent <= exponent_synch;"<<endl;
			o<<tab<<"ResultException <= exception_synch;"<<endl;
			}				
    
			o<<tab<< "ResultSign <= sign_synch;"<<endl;
		}//end else not normalized
       
	}//end if pipelined
	else
	{  
		//the combinational version
		//========================================================================

		o<<tab<<"exponentX <= X("<<wEX + wFX -1<<" downto "<<wFX<<");"<<endl; 
		o<<tab<<"exponentY <= Y("<<wEY + wFY -1<<" downto "<<wFY<<");"<<endl<<endl;

		//Add exponents -> put exponents sum into register
		o<<tab<<"Exponents_Sum <= (\"0\" & exponentX) + (\"0\" & exponentY);"<<endl; //wEX+1 bits
		o<<tab<<"bias <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wEX+1<<");"<<endl; 

		//substract the bias value from the exponent sum
		o<<tab<<"exponents_sum_minus_bias_ext <= Exponents_Sum - bias;"<<endl;//wEX + 1   
		o<<tab<<"Exponents_Sum_Minus_Bias <= exponents_sum_minus_bias_ext("<<wEX-1<<" downto 0);"<<endl; //wEX

		//build significands 	 	
		o<<tab<< "significandX <= \"1\" & X("<<(wFX-1)<<" downto 0);"<<endl;
		   o<<tab<< "significandY <= \"1\" & Y("<<(wFY-1)<<" downto 0);"<<endl<<endl;
			
		//multiply significands 
		o<<tab<< "int_multiplier_component: " << intmult->unique_name << endl;
		o<<tab<< "      port map ( X => significandX, " << endl; //wFX+1 bits (1 bit is the hidden 1)
		o<<tab<< "                 Y => significandY, " << endl; //wFY+1 bits
		o<<tab<< "                 R => significand_product " << endl; //wFX+wFY+2 bits
		o<<tab<< "               );" << endl<<endl;

		//Concatenate the exception bits of the two input numbers
		o<<tab<<"exception_selector <= X("<<wEX + wFX +2<<" downto "<<wEX + wFX + 1<<") & Y("<<wEY + wFY +2<<" downto "<<wEY + wFY +1<<");"<<endl<<endl;

		//select proper exception with selector			
		o<<tab<<"process (exception_selector)"<<endl;
		o<<tab<<"begin"<<endl;
		o<<tab<<"case (exception_selector) is"<< endl;
		o<<tab<<tab<<"when \"0000\" => Result_Exception <= \"00\";   "<<endl;
		o<<tab<<tab<<"when \"0001\" => Result_Exception <= \"00\";   "<<endl;
		o<<tab<<tab<<"when \"0010\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"0011\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"0100\" => Result_Exception <= \"00\";   "<<endl;
		o<<tab<<tab<<"when \"0101\" => Result_Exception <= \"01\";   "<<endl;//to see
		o<<tab<<tab<<"when \"0110\" => Result_Exception <= \"10\";   "<<endl;
		o<<tab<<tab<<"when \"0111\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1000\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1001\" => Result_Exception <= \"10\";   "<<endl;
		o<<tab<<tab<<"when \"1010\" => Result_Exception <= \"10\";   "<<endl;
		o<<tab<<tab<<"when \"1011\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1100\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1101\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1110\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when \"1111\" => Result_Exception <= \"11\";   "<<endl;
		o<<tab<<tab<<"when others =>   Result_Exception <= \"11\";   "<<endl;    
		o<<tab<<"end case;"<<endl;
		o<<tab<<"end process;"<<endl<<endl;

		//use MSB of exponents_sum_minus_bias_ext to siagnal overflow		
		o<<tab<<"Result_Exception_With_EA <= Result_Exception when exponents_sum_minus_bias_ext("<< wEX<<")='0' else"<<endl;
		o<<tab<<"				                    \"11\";"<<endl<<endl;
										
		//normalize result
		//assign selector signal to MSB of signif. multiplic 
		o<<tab<<"mult_selector <= significand_product ("<<wFX+wFY+1<<");"<<endl;

		o<<tab<<"exponent_p <= (\"0\" & Exponents_Sum_Minus_Bias ) + (CONV_STD_LOGIC_VECTOR(0 ,"<<wEX<<") & mult_selector);"<<endl;

		//signal the exception
		o<<tab<<"exception_upd <= Result_Exception_With_EA when exponent_p("<<wEX<<")='0' else"<<endl;
		o<<tab<<"                \"11\";"<<endl;

		//round result
		//check is rounding is needed
		if (1+wFR >= wFX+wFY+2) {
			//=>no rounding needed - possible padding
			//generate zeros of length 1+wFR - wFX+wFY+2 
			zeros1.str("");
			zeros1 << zero_generator(1+wFR - wFX+wFY+2 , 0);

			//make the outputs
			o<<tab<<"ResultExponent <= exponent_p("<<wEX - 1<<" downto 0);"<<endl;  
			o<<tab<<"ResultSignificand <= significand_product & "<<zeros1.str()<<" when mult_selector='1' else"<<endl;
			o<<tab<<"                     significand_product("<<wFX+wFY<<" downto 0) & "<<zeros1.str()<<" & \"0\";"<<endl;
			o<<tab<<"ResultException <= exception_upd;"<<endl;
			o<<tab<<"ResultSign <= X("<<wEX + wFX<<") xor Y("<<wEY + wFY<<");"<<endl;
		}
		else{
			//check if in the middle of two FP numbers

			//generate two xor strings
			str1.str(""); str2.str(""); //init
			str1<<"\"1"<< zero_generator( wFX+wFY+2 - (1+wFR) -1, 1);
			str2<<"\"01"<< zero_generator( wFX+wFY+2 - (1+wFR) -1-1, 1);

			zeros1.str("");
			zeros1<< zero_generator( wFX+wFY+2 - (1+wFR) , 0);

			o<<tab<<"check_string <= (significand_product ("<< wFX+wFY+2 - (1+wFR) -1<<" downto 0) xor "<<str1.str()<<") when mult_selector='1' else"<<endl;
			o<<tab<<"                ( (\"0\" & significand_product ("<< wFX+wFY+2 - (1+wFR) -1-1<<" downto 0)) xor "<<str2.str()<<");"<<endl;

			o<<tab<<"process(check_string)"<<endl;
			o<<tab<<"begin"<<endl;
			o<<tab<<tab<<"      if ( "<<zeros1.str()<<" = check_string ) then"<<endl; 
			o<<tab<<tab<<"         middle_output <= '1';"<<endl;
			o<<tab<<tab<<"      else "<<endl;
			o<<tab<<tab<<"         middle_output <= '0';"<<endl;
			o<<tab<<tab<<"      end if;"<<endl;
			o<<tab<<tab<<"end process; "<<endl;

			o<<tab<<"exponent_p2 <= exponent_p("<<wEX-1<<" downto 0);"<<endl;
				      
			o<<tab<<"reunion_signal <= (\"0\" & exponent_p2 & significand_product("<< wFX+wFY<<" downto "<<wFX+wFY - (wFR) <<" )) when mult_selector='1' else"<<endl;
			o<<tab<<"                  (\"0\" & exponent_p2 & significand_product("<< wFX+wFY-1<<" downto "<<wFX+wFY - (wFR) -1<<" ));"<<endl; 
			               
			//significand_synch_tb2 = LSB of the wFR part of the significand_product 
			o<<tab<<"significand_synch_tb2 <= significand_product("<< wFX+wFY+1 - wFR<<" ) when  mult_selector='1' else"<<endl;
			o<<tab<<"                         significand_product("<< wFX+wFY+1 - wFR -1<<" );"<<endl;              
			   		     
		
			//add 1 to the reunited signal for rounding & normalize purposes
			o<<tab<<"reunion_signal_p <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2+ wEX + wFR<<");"<<endl;
			   		     
			o<<tab<<"middle_signal <= reunion_signal_p when significand_synch_tb2 = '1' else"<<endl;
			o<<tab<<"                 reunion_signal;"<<endl;    		     
			   		     
			o<<tab<<"reunion_signal_p_upd <= reunion_signal_p when middle_output='0' else"<<endl;
			o<<tab<<"                        middle_signal;"<<endl;
			                        
			//update exception                       
			o<<tab<<"ResultException <= exception_upd when reunion_signal_p_upd("<<wEX + wFR+1<<")='0' else"<<endl;
			o<<tab<<"                   \"11\";"<<endl;                                                   

			o<<tab<<"ResultExponent <= reunion_signal_p_upd("<<wEX + wFR<<" downto "<<wFR+1<<");"<<endl;  
			o<<tab<<"ResultSignificand <= \"1\" & reunion_signal_p_upd("<<wFR<<" downto 1);"<<endl;
			o<<tab<<"ResultSign <= X("<<wEX + wFX<<") xor Y("<<wEY + wFY<<");"<<endl;
		}

	}//end the combinational part   
  
	o<< "end architecture;" << endl << endl;
}

 
/**
 * A zero generator method which takes as input two arguments and returns a string of zeros with quotes as stated by the second argurment
 * @param[in] n		    integer argument representing the number of zeros on the output string
 * @param[in] margins	integer argument determining the position of the quotes in the output string. The options are: -2= no quotes; -1=left quote; 0=both quotes 1=right quote
 * @return returns a string of zeros with the corresonding quotes given by margins
 **/
string FPMultiplier::zero_generator(int n, int margins)
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


