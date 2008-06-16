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
#include "FloFP.hpp"

using namespace std;
extern vector<Operator*> oplist;


/**
 * The FPMutliplier constructor
 * @param[in]		target		the target device
 * @param[in]		wEX			the the with of the exponent for the f-p number X
 * @param[in]		wFX			the the with of the fraction for the f-p number X
 * @param[in]		wEY			the the with of the exponent for the f-p number Y
 * @param[in]		wFY			the the with of the fraction for the f-p number Y
 * @param[in]		wER			the the with of the exponent for the multiplication result
 * @param[in]		wFR			the the with of the fraction for the multiplication result
 **/
FPMultiplier::FPMultiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm) :
	Operator(target), wEX(wEX), wFX(wFX), wEY(wEY), wFY(wFY), wER(wER), wFR(wFR) {

	int i, j;
	ostringstream name, synch, synch2;

	/* Set up the status of the operator. Options = sequential|combinatorial */
	if (target->is_pipelined()) 
		set_sequential();
	else
		set_combinatorial(); 

	/* set if operator outputs a normalized result */
	if (norm == 0) 
		normalized = false;
	else
		normalized = true;

	/* The name has the format: FPMultiplier_wEX_wFX_wEY_wFY_wER_wFR where: wEX = width of X exponenet and wFX = width for the fractional part of X */
	name.str("");
	name<<"FPMultiplier_"<<wEX<<"_"<<wFX<<"_"<<wEY<<"_"<<wFY<<"_"<<wER<<"_"<<wFR; 
	unique_name = name.str(); 
	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	add_FP_input ("X", wEX, wFX);
	add_FP_input ("Y", wEY, wFY);
	add_output ("ResultExponent"   , wER    );  
	add_output ("ResultSignificand", wFR + 1);
	add_output ("ResultException"  , 2      );
	add_output ("ResultSign"       , 1      ); 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		   
	/* Instantiate an IntMultiplier -> multiply the significands */
	intmult = new IntMultiplier(target, wFX+1, wFY+1);
	oplist.push_back(intmult);

	//Setup the attibute cocerning the IntMultiplier pipeline depth 
	IntMultPipelineDepth = intmult->pipeline_depth();
	
	/* Common signals for all versions */

	add_signal("significandX", wFX+1);   //The Significands for the input signals X and Y
	add_signal("significandY", wFY+1);
	add_signal("significand_product", wFX +wFY +2);   

	add_signal("exponentX", wEX);	//The signals for the exponents of X and Y
	add_signal("exponentY", wEY);
	add_signal("bias",wEX+2); //bias to be substracted from the sum of exponents  

	add_signal("exception_selector",4);//signal that selects the case for the exception

	
	if (is_sequential())
	{
	/* Common signals for the sequential version */
		
		/* Unregistered signals */
		add_signal("significand_synch", wFX + wFY + 2); 
		add_signal("exponent_synch", wEX+2);
		add_signal("exception_synch",2);
		add_signal("sign_synch", 1); 
		
		/* Registered signals */
		//add_registered_signal_with_sync_reset("Exponents_Sum_Pre_Bias_Substraction",  wEX+2 );
		//add_registered_signal_with_sync_reset("Exponents_Sum_Post_Bias_Substraction", wEX+2 );
		//add_registered_signal_with_sync_reset("Exception_Reg_0",2); 
		//add_registered_signal_with_sync_reset("Exception_Reg_1",2);
		
		add_registered_signal("Exponents_Sum_Pre_Bias_Substraction",  wEX+2 );
		add_registered_signal("Exponents_Sum_Post_Bias_Substraction", wEX+2 );
		add_registered_signal("Exception_Reg_0",2); 
		add_registered_signal("Exception_Reg_1",2);
		
		
		
			/* Synchronization registers */
			
			//when integer multiplication has more than 2 pipeline levels
			for (i=0;i<=IntMultPipelineDepth-3;i++)	{
				synch.str(""); synch2.str("");
				synch<<"Exponent_Synchronization_"<<i;
				
				//add_registered_signal_with_sync_reset(synch.str(), wEX+2 );
				add_registered_signal(synch.str(), wEX+2 );
				
				synch2<<"Exception_Synchronization_"<<i;
				//add_registered_signal_with_sync_reset(synch2.str(), 2 );
				add_registered_signal(synch2.str(), 2 );
			}
			//when the integer multiplication has less than 2 levels
			for (i=0;i<=1-IntMultPipelineDepth;i++)	{
				synch.str("");
				synch<<"Int_Synchronization_"<<i;
				//add_registered_signal_with_sync_reset(synch.str(), wFX+wFY+2);
				add_registered_signal(synch.str(), wFX+wFY+2);				
			}	
			add_delay_signal("Sign_Synchronization", 1, max(2,IntMultPipelineDepth));

		if (normalized){
			/* The sequential normalized version */
		
			add_signal("normalization_selector",1);
			add_signal("exponent_post_normalization",2+wEX);
			add_signal("exception_post_normalization",2);
			
			if (1+wFR < wFX+wFY+2) {
				
				//add_registered_signal_with_sync_reset("exponent_synch2",wEX+2);
				add_registered_signal("exponent_synch2",wEX+2);
				
				add_registered_signal_with_sync_reset("exception_synch2",2);
				
				//add_registered_signal_with_sync_reset("significand_synch2",wFX + wFY + 2);
				add_registered_signal("significand_synch2",wFX + wFY + 2);
				
				
				add_registered_signal_with_sync_reset("sign_synch2",1);
				add_registered_signal_with_sync_reset("normalization_selector2",1);
		
				add_signal("between_fp_numbers",1);	
				add_signal("LSB_of_result_significand_out",1);	
				add_signal("between_fp_numbers_out",1);
				add_signal("sign_synch2_out",1);
				add_signal("exception_synch2_out",2);
				add_signal("exponent_synch2_out",wER+2);
				add_signal("exponent_synch2_out_post_rounding",wER+2);
				
						
				add_signal("check_string",wFX+wFY+2 - (1+wFR) );
				add_signal("reunion_signal",2+wFR);
				add_signal("reunion_signal_out",2+wFR);
				add_signal("LSB_of_result_significand",1); 
				add_signal("between_fp_numbers_result_significand",2+wFR);
				add_signal("reunion_signal_post_addition",2+wFR); 
				add_signal("reunion_signal_post_rounding",2+wFR); 
				
				
				//parameters setup
				
				bool status = target->suggest_subadd_size(addition_chunk_width, 2+ wFR);
				if (!status)
					cout<<"Frequency report:"<<endl;
					cout<<tab<<"WARNING: Desired frequency cannot be reached !";
			
			
				reunion_signal_width = 2+ wFR;
				reunion_signal_parts = int ( ceil( double(reunion_signal_width)/double(addition_chunk_width)));
				addition_last_chunk_width = reunion_signal_width - (reunion_signal_parts - 1)*addition_chunk_width;
							
				
				
				if (reunion_signal_parts>1){
					for (j=1; j<=reunion_signal_parts;j++)	
						for (i=1;i<=reunion_signal_parts;i++){	
							name.str("");
							name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
							if (i!=reunion_signal_parts){
								//add_delay_signal_bus_(name.str(), addition_chunk_width + 1 ,1);
								add_delay_signal_bus_no_reset(name.str(), addition_chunk_width + 1 ,1);
							}
							else
							{
								//add_delay_signal_bus(name.str(), addition_last_chunk_width,1 );	
								add_delay_signal_bus_no_reset(name.str(), addition_last_chunk_width,1 );	
							
							}
		        }
				
					/*
					add_delay_signal("reunion_signal_level",2+wFR,reunion_signal_parts);
					add_delay_signal("LSB_of_result_significand_level",1,reunion_signal_parts);		
					add_delay_signal("between_fp_numbers_level",1,reunion_signal_parts);			
					add_delay_signal("sign_synch2_level",1,reunion_signal_parts);
					add_delay_signal("exception_synch2_level",2,reunion_signal_parts);
					add_delay_signal("exponent_synch2_level",2+wER,reunion_signal_parts);
					*/
					
					add_delay_signal_no_reset("reunion_signal_level",2+wFR,reunion_signal_parts);
					add_delay_signal_no_reset("LSB_of_result_significand_level",1,reunion_signal_parts);		
					add_delay_signal_no_reset("between_fp_numbers_level",1,reunion_signal_parts);			
					add_delay_signal_no_reset("sign_synch2_level",1,reunion_signal_parts);
					add_delay_signal_no_reset("exception_synch2_level",2,reunion_signal_parts);
					add_delay_signal_no_reset("exponent_synch2_level",2+wER,reunion_signal_parts);
					
					
				}
			}
		
		}else{
		/* The sequential non-normalized version */
			
			
			if ((wFX+wFY+2 > wFR+1))
			{
				/*
				add_registered_signal_with_sync_reset("exponent_synch2",wEX+2);
				add_registered_signal_with_sync_reset("exception_synch2",2);
				add_registered_signal_with_sync_reset("significand_synch2",wFX + wFY + 2);
				add_registered_signal_with_sync_reset("sign_synch2",1);
				*/
				
				add_registered_signal("exponent_synch2",wEX+2);
				add_registered_signal("exception_synch2",2);
				add_registered_signal("significand_synch2",wFX + wFY + 2);
				add_registered_signal("sign_synch2",1);
						
				add_signal("between_fp_numbers",1);	
				add_signal("LSB_of_result_significand_out",1);	
				add_signal("between_fp_numbers_out",1);
				add_signal("sign_synch2_out",1);
				add_signal("exception_synch2_out",2);
				add_signal("exponent_synch2_out",wER+2);
				add_signal("exponent_synch2_out_post_rounding",wER+2);
										
				add_signal("check_string",wFX+wFY+2 - (1+wFR) );
				
			
			
				bool status = target->suggest_subadd_size(addition_chunk_width, 3+ wFR);
				if (!status)
					cout<<"Frequency report:"<<endl;
					cout<<tab<<"WARNING: Desired frequency cannot be reached !";
						
				reunion_signal_width = 3 + wFR;
				reunion_signal_parts = int ( ceil( double(reunion_signal_width)/double(addition_chunk_width)));
				addition_last_chunk_width = reunion_signal_width - (reunion_signal_parts - 1)*addition_chunk_width;
				
				if (reunion_signal_parts>1){
					for (j=1; j<=reunion_signal_parts;j++)	
					for (i=1;i<=reunion_signal_parts;i++){	
						name.str("");
						name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
						if (i!=reunion_signal_parts)
						{
							//add_registered_signal_with_sync_reset(name.str(), addition_chunk_width + 1 );
							add_registered_signal_with_sync_reset(name.str(), addition_chunk_width + 1 );
						}
						else{
							//add_registered_signal_with_sync_reset(name.str(), addition_last_chunk_width );
							add_registered_signal(name.str(), addition_last_chunk_width );		
						}
	        }
				
					/*
					add_delay_signal("reunion_signal_level",3+wFR,reunion_signal_parts);
					add_delay_signal("LSB_of_result_significand_level",1,reunion_signal_parts);		
					add_delay_signal("between_fp_numbers_level",1,reunion_signal_parts);			
					add_delay_signal("sign_synch2_level",1,reunion_signal_parts);
					add_delay_signal("exception_synch2_level",2,reunion_signal_parts);
					add_delay_signal("exponent_synch2_level",2+wER,reunion_signal_parts);
					*/
					
					add_delay_signal_no_reset("reunion_signal_level",3+wFR,reunion_signal_parts);
					add_delay_signal_no_reset("LSB_of_result_significand_level",1,reunion_signal_parts);		
					add_delay_signal_no_reset("between_fp_numbers_level",1,reunion_signal_parts);			
					add_delay_signal_no_reset("sign_synch2_level",1,reunion_signal_parts);
					add_delay_signal_no_reset("exception_synch2_level",2,reunion_signal_parts);
					add_delay_signal_no_reset("exponent_synch2_level",2+wER,reunion_signal_parts);					
					
					
				}

				add_signal("reunion_signal"                       ,1+(1+wFR)+1);
				add_signal("reunion_signal_out"                   ,3+wFR);
				add_signal("LSB_of_result_significand"            ,1    ); 
				add_signal("between_fp_numbers_result_significand",3+wFR);
				add_signal("reunion_signal_post_addition"         ,3+wFR); 
				add_signal("reunion_signal_post_rounding"         ,3+wFR);
						
				//add_signal("temp_exp_concat_fract", 1 + wEX + wFX + wFY + 2);
			}
		}
	}
	else{ 
		/* Signals for the combinational version */
		add_signal("normalization_selector",1);
		add_signal("exponent_post_normalization",2+wER);
		add_signal("exception_post_normalization",2);
				
		add_signal("Exponents_Sum_Pre_Bias_Substraction", wER+2 );
		add_signal("Exponents_Sum_Post_Bias_Substraction", wER+2 );   
		add_signal("Exception_Reg_0",2); 
		add_signal("Exception_Reg_1",2);
		if ((1 + wFR < wFX + wFY + 2)) {
			add_signal("check_string",wFX+wFY+2 - (1+wFR) );
			add_signal("reunion_signal",2+wFR);
			add_signal("LSB_of_result_significand",1); 
			add_signal("between_fp_numbers_result_significand",2+wFR);
			add_signal("reunion_signal_post_addition",2+wFR); 
			add_signal("reunion_signal_post_rounding",2+wFR); 
			add_signal("exponent_post_rounding",wER+2);
			add_signal("between_fp_numbers",1);
		}
		
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
			if (1+wFR>=2+wFX+wFY)
				set_pipeline_depth(max(2,IntMultPipelineDepth));
			else
				if (reunion_signal_parts==1)
					set_pipeline_depth(max(2,IntMultPipelineDepth) + 1);
				else
					set_pipeline_depth(max(2,IntMultPipelineDepth) + 1 + reunion_signal_parts);
		
}

/**
 * FPMultiplier destructor
 */
FPMultiplier::~FPMultiplier() {
}



/**
 * Method belonging to the Operator class overloaded by the FPMultiplier class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
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
		/* common code for both normalized and non-normalized version */
		output_vhdl_registers(o); o<<endl;
		
		/* Exponent Handling */
			o<<tab<<"exponentX <= X("<<wEX + wFX -1<<" downto "<<wFX<<");"<<endl; 
			o<<tab<<"exponentY <= Y("<<wEY + wFY -1<<" downto "<<wFY<<");"<<endl<<endl;
			//Add exponents -> put sum into register
			o<<tab<<"Exponents_Sum_Pre_Bias_Substraction <= (\"00\" & exponentX) + (\"00\" & exponentY);"<<endl; //wEX+2 bits	[*Optimization point - pipelining*]
			o<<tab<<"bias <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wEX+2<<");"<<endl; 
			//substract the bias value from the exponents' sum
			o<<tab<<"Exponents_Sum_Post_Bias_Substraction <= Exponents_Sum_Pre_Bias_Substraction_d - bias;"<<endl;//wEX + 2   
		
		/* Significand Handling */
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

		/* Exception Handling */
			//Concatenate the exception bits of the two input numbers
			o<<tab<<"exception_selector <= X("<<wEX + wFX +2<<" downto "<<wEX + wFX + 1<<") & Y("<<wEY + wFY +2<<" downto "<<wEY + wFY +1<<");"<<endl<<endl;
			//select proper exception with selector			
			o<<tab<<"process (exception_selector)"<<endl;
			o<<tab<<"begin"<<endl;
			o<<tab<<"case (exception_selector) is"<< endl;
			o<<tab<<tab<<"when \"0000\" => Exception_Reg_0 <= \"00\";   "<<endl;
			o<<tab<<tab<<"when \"0001\" => Exception_Reg_0 <= \"00\";   "<<endl;
			o<<tab<<tab<<"when \"0010\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"0011\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"0100\" => Exception_Reg_0 <= \"00\";   "<<endl;
			o<<tab<<tab<<"when \"0101\" => Exception_Reg_0 <= \"01\";   "<<endl;
			o<<tab<<tab<<"when \"0110\" => Exception_Reg_0 <= \"10\";   "<<endl;
			o<<tab<<tab<<"when \"0111\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1000\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1001\" => Exception_Reg_0 <= \"10\";   "<<endl;
			o<<tab<<tab<<"when \"1010\" => Exception_Reg_0 <= \"10\";   "<<endl;
			o<<tab<<tab<<"when \"1011\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1100\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1101\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1110\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1111\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when others =>   Exception_Reg_0 <= \"11\";   "<<endl;    
			o<<tab<<"end case;"<<endl;
			o<<tab<<"end process;"<<endl<<endl;
			
			//propagate 1 level the computed exception in order to synchronize with the exponents sum						
			o<<tab<<"Exception_Reg_1 <= Exception_Reg_0_d;"<<endl; 						
									
		/* synchronization barrier (IntMultiplication | Exponent Addition | Exception Computation | Sign Computation */
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
		/* Cases: 		
			0. when IntMultiplication | Exponent Addition are allready synchronized	
			1. when registers are added to exponent addition
			2. when registers are added to integer multiplication
		*/	
	  
		if (IntMultPipelineDepth==2){
			// when IntMultiplication | Exponent Addition are allready synchronized	 
			o<<tab<<"exponent_synch    <= Exponents_Sum_Post_Bias_Substraction_d;"<<endl;
			o<<tab<<"significand_synch <= significand_product;"<<endl;
			o<<tab<<"exception_synch   <= Exception_Reg_1_d;"<<endl;
		}
	   else{	
			if (IntMultPipelineDepth > 2){
			/* when registers are added to exponent addition */
				//connect the first synchronization reg to the output of the last reg from exponent addition 
				o<<tab<<"Exponent_Synchronization_0  <= Exponents_Sum_Post_Bias_Substraction_d;"<<endl;
				o<<tab<<"Exception_Synchronization_0 <= Exception_Reg_1_d;"<<endl; 									

				for (i=1;i<=IntMultPipelineDepth-3;i++){
					o<<tab<<"Exponent_Synchronization_"<<i<<"   <= "<<"Exponent_Synchronization_"<<i-1<<"_d;"<<endl;		
					o<<tab<<"Exception_Synchronization_"<<i<<" <= "<<"Exception_Synchronization_"<<i-1<<"_d;"<<endl;
				}
					o<<tab<<"exponent_synch    <= Exponent_Synchronization_"<<	IntMultPipelineDepth-3 <<"_d;"<<endl;
					o<<tab<<"exception_synch   <= Exception_Synchronization_"<<IntMultPipelineDepth-3 <<"_d;"<<endl;
					o<<tab<<"significand_synch <= significand_product;"<<endl;			
			}
			else{
				//add registers to the output of IntMultiplier so that exponent addition and significand multiplication become synchronized
				o<<tab<<"Int_Synchronization_0 <= significand_product;"<<endl;
				for (i=1;i<=1-IntMultPipelineDepth;i++){
					synch1.str(""); synch2.str("");
					synch2<<"Int_Synchronization_"<<i;
					synch1<<"Int_Synchronization_"<<i-1<<"_d";
					o<<tab<<synch2.str()<<"<="<<synch1.str()<<";"<<endl;		
				}
				o<<tab<<"exponent_synch    <= Exponents_Sum_Post_Bias_Substraction_d;"<<endl;
				o<<tab<<"exception_synch   <= Exception_Reg_1_d;"<<endl;
				o<<tab<<"significand_synch <= Int_Synchronization_"<< 1-IntMultPipelineDepth<<"_d;"<<endl;
			}
		}
		//sign synchronization
		o<<tab<<"Sign_Synchronization <= X("<<wEX + wFX<<") xor Y("<<wEY + wFY<<");"<<endl;
		o<<tab<<"sign_synch <= "<<get_delay_signal_name("Sign_Synchronization",max(2,IntMultPipelineDepth))<<";"<<endl<<endl;
	
	
		if (normalized==true){	
			/* normalize result */
			
			//assign selector signal to MSB of signif. multiplic 
			o<<tab<<"normalization_selector <= significand_synch ("<<wFX+wFY+1<<");"<<endl;

			o<<tab<<"exponent_post_normalization <= (exponent_synch ) + (CONV_STD_LOGIC_VECTOR(0 ,"<<wEX+1<<") & normalization_selector);"<<endl;

			/* result rounding */
			//check is rounding is needed
			if (1+wFR >= wFX+wFY+2) {
				/* =>no rounding needed - possible padding */
				o<<tab<<"with exponent_post_normalization("<< wER+1 <<" downto "<< wER <<") select"<<endl;		
				o<<tab<<"exception_post_normalization <= exception_synch when \"00\","<<endl;
				o<<tab<<"                            \"10\"             when \"01\", "<<endl;
				o<<tab<<"                            \"00\"             when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"             when others;"<<endl;						
				
				//make the outputs
				o<<tab<<"ResultExponent     <= exponent_post_normalization("<<wEX - 1<<" downto 0);"<<endl;  
				o<<tab<<"ResultSignificand  <= significand_synch & "<<zero_generator(1+wFR - (wFX+wFY+2) , 0)<<" when normalization_selector='1' else"<<endl;
				o<<tab<<"                      significand_synch("<<wFX+wFY<<" downto 0) & "<<zero_generator(1+wFR - (wFX+wFY+2) , 0)<<" & \"0\";"<<endl;
				o<<tab<<"ResultException    <= exception_post_normalization;"<<endl;
				o<<tab<<"ResultSign         <= sign_synch;"<<endl;

			}
			else{ 
				o<<tab<<"exception_post_normalization <= exception_synch;"<<endl;
				//check if in the middle of two FP numbers
				//generate two xor strings
				str1.str(""); str2.str(""); //init
				str1<<"\"1"<< zero_generator( wFX+wFY+2 - (1+wFR) -1, 1);
				str2<<"\"01"<< zero_generator( wFX+wFY+2 - (1+wFR) -1-1, 1);

				zeros1.str("");
				zeros1<< zero_generator( wFX+wFY+2 - (1+wFR) , 0);

				o<<tab<<"check_string <= (significand_synch ("<< wFX+wFY+2 - (1+wFR) -1<<" downto 0) xor "<<str1.str()<<") when normalization_selector='1' else"<<endl;
				o<<tab<<"               ( (\"0\" & significand_synch ("<< wFX+wFY+2 - (1+wFR) -1-1<<" downto 0)) xor "<<str2.str()<<");"<<endl;

				o<<tab<<"process(clk)"<<endl;
				o<<tab<<"begin"<<endl;
				o<<tab<<tab<<"   if (clk'event and clk='1') then "<<endl;  
				o<<tab<<tab<<"      if ( "<<zeros1.str()<<" = check_string ) then"<<endl; 
				o<<tab<<tab<<"         between_fp_numbers <= '1';"<<endl;
				o<<tab<<tab<<"      else "<<endl;
				o<<tab<<tab<<"         between_fp_numbers <= '0';"<<endl;
				o<<tab<<tab<<"      end if;"<<endl;
				o<<tab<<tab<<"   end if; "<<endl;
				o<<tab<<tab<<"end process; "<<endl;

				//propagate rest of signals 1 level in order to compensate with the 1 bit delay caused by above comparisson
				o<<tab<<"exponent_synch2         <= exponent_post_normalization;"<<endl;
				o<<tab<<"exception_synch2        <= exception_post_normalization;"<<endl;
				o<<tab<<"significand_synch2      <= significand_synch;"<<endl;
				o<<tab<<"sign_synch2             <= sign_synch;"<<endl;
				o<<tab<<"normalization_selector2 <= normalization_selector;"<<endl;
				//~~~~~~~~~~~~~~~~~~~~~


				o<<tab<<"reunion_signal <= (\"0\" & significand_synch2_d("<< wFX+wFY<<" downto "<<wFX+wFY - (wFR) <<" )) when normalization_selector2_d='1' else"<<endl;
				o<<tab<<"                  (\"0\" & significand_synch2_d("<< wFX+wFY-1<<" downto "<<wFX+wFY - (wFR) -1<<" ));"<<endl; 
              
				o<<tab<<"LSB_of_result_significand <= significand_synch2_d("<< wFX+wFY+1 - wFR<<" ) when normalization_selector2_d='1' else"<<endl;
				o<<tab<<"                             significand_synch2_d("<< wFX+wFY+1 - wFR -1<<" );"<<endl;              
				
				/* pipeline the possible carry_in propagation caused by rounding */
				
				//connect the first part of the pipeline structure to the signals
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
					
					//put together the pipeline
					for (j=1; j<=reunion_signal_parts-1;j++)	
						for (i=1;i<=reunion_signal_parts;i++)
							if ((i<=j)||(i>j+1))
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
							else
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<addition_chunk_width<<") ;"<<endl; 

					//put the result back together
					ostringstream reunion_signal_post_1_addition;
					for (i = reunion_signal_parts; i > 0; i--) 
						if ( (i!=1) && (i!=reunion_signal_parts) )
							reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0) &" ;
						else
							if (i==reunion_signal_parts)
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_last_chunk_width -1<<" downto 0) & ";
							else
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0)" ;
					
				
					//propagate the rest of the signals reunion_signal_parts steps
					
					o<<tab<<"reunion_signal_level <=reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_out <= "<<get_delay_signal_name("reunion_signal_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"LSB_of_result_significand_level <= LSB_of_result_significand;"<<endl;
					o<<tab<<"LSB_of_result_significand_out <= "<<get_delay_signal_name("LSB_of_result_significand_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"between_fp_numbers_level <= between_fp_numbers;"<<endl;
					o<<tab<<"between_fp_numbers_out <= "<<get_delay_signal_name("between_fp_numbers_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"sign_synch2_level <= sign_synch2_d;"<<endl;
					o<<tab<<"sign_synch2_out <= "<<get_delay_signal_name("sign_synch2_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"exception_synch2_level <= exception_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out <= "<<get_delay_signal_name("exception_synch2_level", reunion_signal_parts)<<";" <<endl;		
														
					o<<tab<<"reunion_signal_post_addition <= "<< reunion_signal_post_1_addition.str()<<";"<< endl;
					
					o<<tab<<"exponent_synch2_level <= exponent_synch2_d;"<<endl;
					o<<tab<<"exponent_synch2_out <= "<<get_delay_signal_name("exponent_synch2_level", reunion_signal_parts)<<";" <<endl;						
				
				}   
				else{
					/* when the carry propagation on the reunion signal is not pipelined */
					//add 1 to the reunited signal for rounding & normalize purposes
					o<<tab<<"reunion_signal_out            <= reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_post_addition  <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2+ wFR<<");"<<endl;
					o<<tab<<"between_fp_numbers_out             <= between_fp_numbers;"<<endl;
					o<<tab<<"LSB_of_result_significand_out <= LSB_of_result_significand;"<<endl;
					o<<tab<<"sign_synch2_out               <= sign_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out          <= exception_synch2_d;"<<endl;
					o<<tab<<"exponent_synch2_out <= exponent_synch2_d;"<<endl;
					
				}				
				   		     
				   		     
				o<<tab<<"between_fp_numbers_result_significand <= reunion_signal_post_addition when LSB_of_result_significand_out = '1' else"<<endl;
				o<<tab<<"                                         reunion_signal_out;"<<endl;    		     
				   		     
				o<<tab<<"reunion_signal_post_rounding <= reunion_signal_post_addition when between_fp_numbers_out='0' else"<<endl;
				o<<tab<<"                       	     between_fp_numbers_result_significand;"<<endl;
				                        
				
				
				o<<tab<<"exponent_synch2_out_post_rounding <= exponent_synch2_out + (CONV_STD_LOGIC_VECTOR(0,"<<wEX+1<<") & reunion_signal_post_rounding("<<1+wFR<<"));"<<endl;    
				
				                        
				//update exception 
				o<<tab<<"with exponent_synch2_out_post_rounding("<<wER+1<<" downto "<< wER <<") select"<<endl;                       
				o<<tab<<"ResultException   <= exception_synch2_out      when \"00\","<<endl;
				o<<tab<<"                            \"10\"             when \"01\", "<<endl;
				o<<tab<<"                            \"00\"             when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"             when others;"<<endl;						
			                
				                        
				                        
				//update exception                       
				//o<<tab<<"ResultException   <= exception_synch2_out when reunion_signal_post_rounding("<<wEX + wFR+1<<")='0' else"<<endl;
				//o<<tab<<"                     \"10\";"<<endl;                                                   

				o<<tab<<"ResultExponent    <= exponent_synch2_out_post_rounding("<<wER-1<<" downto "<<0<<");"<<endl;  
				o<<tab<<"ResultSignificand <= \"1\" & reunion_signal_post_rounding("<<wFR<<" downto 1);"<<endl;
				o<<tab<<"ResultSign        <= sign_synch2_out;"<<endl;
			}
		}
		else { //TODO
		/* Non-Normalized case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			if (1+wFR < wFX+wFY+2){
				//check if in the middle of two FP numbers
				//generate two xor strings
				str1.str(""); 
				str1<<"\"1"<< zero_generator( wFX+wFY+2 - (1+wFR) -1, 1);
				
				zeros1.str("");
				zeros1<< zero_generator( wFX+wFY+2 - (1+wFR) , 0);

				o<<tab<<"check_string <= (significand_synch ("<< wFX+wFY+2 - (1+wFR) -1<<" downto 0) xor "<<str1.str()<<"); "<<endl;
				
				o<<tab<<"process(clk)"<<endl;
				o<<tab<<"begin"<<endl;
				o<<tab<<tab<<"   if (clk'event and clk='1') then "<<endl;  
				o<<tab<<tab<<"      if ( "<<zeros1.str()<<" = check_string ) then"<<endl; 
				o<<tab<<tab<<"         between_fp_numbers <= '1';"<<endl;
				o<<tab<<tab<<"      else "<<endl;
				o<<tab<<tab<<"         between_fp_numbers <= '0';"<<endl;
				o<<tab<<tab<<"      end if;"<<endl;
				o<<tab<<tab<<"   end if; "<<endl;
				o<<tab<<tab<<"end process; "<<endl;

				//propagate rest of signals 1 level in order to compensate with the 1 bit delay caused by above comparisson
				o<<tab<<"exponent_synch2         <= exponent_synch;"<<endl;
				o<<tab<<"exception_synch2        <= exception_synch;"<<endl;
				o<<tab<<"significand_synch2      <= significand_synch;"<<endl;
				o<<tab<<"sign_synch2             <= sign_synch;"<<endl;
				//~~~~~~~~~~~~~~~~~~~~~
				
				o<<tab<<"reunion_signal <= (\"0\"  & significand_synch2_d("<< wFX+wFY+1<<" downto "<<wFX+wFY+1 - (wFR+1) <<" ));"<<endl;; 
              
				o<<tab<<"LSB_of_result_significand <= significand_synch2_d("<< wFX+wFY+1 - wFR<<" );"<<endl;
				
				/* pipeline the possible carry_in propagation caused by rounding */
				
				//connect the first part of the pipeline structure to the signals
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
					
					//put together the pipeline
					for (j=1; j<=reunion_signal_parts-1;j++)	
						for (i=1;i<=reunion_signal_parts;i++)
							if ((i<=j)||(i>j+1))
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
							else
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<addition_chunk_width<<") ;"<<endl; 

					//put the result back together
					ostringstream reunion_signal_post_1_addition;
					for (i = reunion_signal_parts; i > 0; i--) 
						if ( (i!=1) && (i!=reunion_signal_parts) )
							reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0) &" ;
						else
							if (i==reunion_signal_parts)
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_last_chunk_width -1<<" downto 0) & ";
							else
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunion_signal_parts<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0)" ;
					
				
					//propagate the rest of the signals reunion_signal_parts steps
					
					o<<tab<<"reunion_signal_level <=reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_out <= "<<get_delay_signal_name("reunion_signal_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"LSB_of_result_significand_level <= LSB_of_result_significand;"<<endl;
					o<<tab<<"LSB_of_result_significand_out <= "<<get_delay_signal_name("LSB_of_result_significand_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"between_fp_numbers_level <= between_fp_numbers;"<<endl;
					o<<tab<<"between_fp_numbers_out <= "<<get_delay_signal_name("between_fp_numbers_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"sign_synch2_level <= sign_synch2_d;"<<endl;
					o<<tab<<"sign_synch2_out <= "<<get_delay_signal_name("sign_synch2_level", reunion_signal_parts)<<";" <<endl;		
					
					o<<tab<<"exception_synch2_level <= exception_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out <= "<<get_delay_signal_name("exception_synch2_level", reunion_signal_parts)<<";" <<endl;		
														
					o<<tab<<"reunion_signal_post_addition <= "<< reunion_signal_post_1_addition.str()<<";"<< endl;
					
					o<<tab<<"exponent_synch2_level <= exponent_synch2_d;"<<endl;
					o<<tab<<"exponent_synch2_out <= "<<get_delay_signal_name("exponent_synch2_level", reunion_signal_parts)<<";" <<endl;						
				}   
				else{
					/* when the carry propagation on the reunion signal is not pipelined */
					//add 1 to the reunited signal for rounding & normalize purposes
					o<<tab<<"reunion_signal_out            <= reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_post_addition  <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2+wFR<<");"<<endl;
					o<<tab<<"between_fp_numbers_out             <= between_fp_numbers;"<<endl;
					o<<tab<<"LSB_of_result_significand_out <= LSB_of_result_significand;"<<endl;
					o<<tab<<"sign_synch2_out               <= sign_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out          <= exception_synch2_d;"<<endl;
					o<<tab<<"exponent_synch2_out           <= exponent_synch2_d;"<<endl;
				}				
				
		     
				   		     
				o<<tab<<"between_fp_numbers_result_significand <= reunion_signal_post_addition when LSB_of_result_significand_out = '1' else"<<endl;
				o<<tab<<"                                         reunion_signal_out;"<<endl;    		     
				   		     
				o<<tab<<"reunion_signal_post_rounding <= reunion_signal_post_addition when between_fp_numbers_out='0' else"<<endl;
				o<<tab<<"                       	     between_fp_numbers_result_significand;"<<endl;
				
				o<<tab<<"exponent_synch2_out_post_rounding <= exponent_synch2_out + (CONV_STD_LOGIC_VECTOR(0,"<<wEX+1<<") & reunion_signal_post_rounding("<<2+wFR<<"));"<<endl;    
				
				                        
				//update exception 
				o<<tab<<"with exponent_synch2_out_post_rounding("<<wER+1<<" downto "<< wER <<") select"<<endl;                       
				o<<tab<<"ResultException   <= exception_synch2_out      when \"00\","<<endl;
				o<<tab<<"                            \"10\"             when \"01\", "<<endl;
				o<<tab<<"                            \"00\"             when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"             when others;"<<endl;						
			

				o<<tab<<"ResultExponent    <= exponent_synch2_out_post_rounding("<<wEX -1<<" downto "<<0<<");"<<endl;  
					
				o<<tab<<"ResultSignificand <= (\"11\" & reunion_signal_post_rounding("<<wFR-1<<" downto 1)) when reunion_signal_out("<<wFR+1<<" downto "<<wFR<<")=\"11\" else"<<endl;
				o<<tab<<"					  reunion_signal_post_rounding("<<wFR+1<<" downto 1);"<<endl;	
				o<<tab<<"ResultSign        <= sign_synch2_out;"<<endl;

			}
			else{
			/* No rounding is needed. If wFR+1>wFX+wFY+2 then the ResultSignificand must be padded with 0s to the right */
				if (wFX+wFY+2 == wFR+1)
					o<<tab<<"ResultSignificand <= significand_synch;"<<endl;
				else 	
					o<<tab<<"ResultSignificand <= significand_synch & "<<zero_generator((wFR+1) - (wFX+wFY+2), 0 )<<";"<<endl;
					
				o<<tab<<"with exponent_synch("<<wER+1<<" downto "<< wER <<") select"<<endl;                       
				o<<tab<<"ResultException   <= exception_synch  when \"00\","<<endl;
				o<<tab<<"                            \"10\"    when \"01\", "<<endl;
				o<<tab<<"                            \"00\"    when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"    when others;"<<endl;							
					
				o<<tab<<"ResultExponent  <= exponent_synch("<<wER-1<<" downto 0 "<< ");"<<endl;
				
				o<<tab<<"ResultSign <= sign_synch;"<<endl;
			}				
    
			
		}//end else not normalized
       
	}//end if pipelined
	else
	{  
		/* the combinational version */
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
		/* Exponent Handling */
			// Fetch exponents from inputs 
			o<<tab<<"exponentX <= X("<<wEX + wFX -1<<" downto "<<wFX<<");"<<endl; 
			o<<tab<<"exponentY <= Y("<<wEY + wFY -1<<" downto "<<wFY<<");"<<endl<<endl;
			//Add exponents and put sum into register
			o<<tab<<"Exponents_Sum_Pre_Bias_Substraction <= (\"00\" & exponentX) + (\"00\" & exponentY);"<<endl; //wEX+2 bits
			o<<tab<<"bias <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wER+2<<");"<<endl; 
			//substract the bias value from the exponent sum
			o<<tab<<"Exponents_Sum_Post_Bias_Substraction <= Exponents_Sum_Pre_Bias_Substraction - bias;"<<endl;//wEX + 2   

		/* Significand Handling */
			//fetch and build significands by adding a "1" in MSB position 
			o<<tab<< "significandX <= \"1\" & X("<<(wFX-1)<<" downto 0);"<<endl;
			o<<tab<< "significandY <= \"1\" & Y("<<(wFY-1)<<" downto 0);"<<endl<<endl;
			//multiply significands 
			o<<tab<< "int_multiplier_component: " << intmult->unique_name << endl;
			o<<tab<< "      port map ( X => significandX, " << endl; //wFX+1 bits (1 bit is the hidden 1)
			o<<tab<< "                 Y => significandY, " << endl; //wFY+1 bits
			o<<tab<< "                 R => significand_product " << endl; //wFX+wFY+2 bits
			o<<tab<< "               );" << endl<<endl;

		/* Exception bits Handling */
			//Concatenate the exception bits of the two input numbers
			o<<tab<<"exception_selector <= X("<<wEX + wFX +2<<" downto "<<wEX + wFX + 1<<") & Y("<<wEY + wFY +2<<" downto "<<wEY + wFY +1<<");"<<endl<<endl;
			//select proper exception with selector			
			o<<tab<<"process (exception_selector)"<<endl;
			o<<tab<<"begin"<<endl;
			o<<tab<<"case (exception_selector) is"<< endl;
			o<<tab<<tab<<"when \"0000\" => Exception_Reg_0 <= \"00\";   "<<endl;
			o<<tab<<tab<<"when \"0001\" => Exception_Reg_0 <= \"00\";   "<<endl;
			o<<tab<<tab<<"when \"0010\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"0011\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"0100\" => Exception_Reg_0 <= \"00\";   "<<endl;
			o<<tab<<tab<<"when \"0101\" => Exception_Reg_0 <= \"01\";   "<<endl;//to see
			o<<tab<<tab<<"when \"0110\" => Exception_Reg_0 <= \"10\";   "<<endl;
			o<<tab<<tab<<"when \"0111\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1000\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1001\" => Exception_Reg_0 <= \"10\";   "<<endl;
			o<<tab<<tab<<"when \"1010\" => Exception_Reg_0 <= \"10\";   "<<endl;
			o<<tab<<tab<<"when \"1011\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1100\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1101\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1110\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when \"1111\" => Exception_Reg_0 <= \"11\";   "<<endl;
			o<<tab<<tab<<"when others =>   Exception_Reg_0 <= \"11\";   "<<endl;    
			o<<tab<<"end case;"<<endl;
			o<<tab<<"end process;"<<endl<<endl;
			
			
		/* Normalization */
			// assign selector signal to MSB of signif. multiplication 
			o<<tab<<"normalization_selector <= significand_product ("<<wFX+wFY+1<<");"<<endl;
			// add the selector bit to the exponent in order to normalize the result 
			o<<tab<<"exponent_post_normalization <= Exponents_Sum_Post_Bias_Substraction + (CONV_STD_LOGIC_VECTOR(0 ,"<<wER+1<<") & normalization_selector);"<<endl;
			o<<tab<<"exception_post_normalization <= Exception_Reg_0;"<<endl;			

		/* Rounding */
			if (1+wFR >= wFX+wFY+2) {
				/* No rounding needed */
				/* Possible 0 padding to the length of the fractional part of the result */
				o<<tab<<"with exponent_post_normalization("<< wER+1 <<" downto "<< wER <<") select"<<endl;		
				o<<tab<<" ResultException <= exception_post_normalization when \"00\","<<endl;
				o<<tab<<"                            \"10\"                           when \"01\", "<<endl;
				o<<tab<<"                            \"00\"                           when \"11\"|\"10\","<<endl;
				o<<tab<<"                             \"11\"                          when others;"<<endl;
											
				// Assign the architecture outputs
				o<<tab<<"ResultExponent    <= exponent_post_normalization("<<wEX - 1<<" downto 0);"<<endl;  
				o<<tab<<"ResultSignificand <= significand_product & "<< zero_generator(1+wFR - (wFX+wFY+2) , 0) <<" when normalization_selector='1' else"<<endl;
				o<<tab<<"                     significand_product("<<wFX+wFY<<" downto 0) & "<<zero_generator(1+wFR - (wFX+wFY+2) , 0)<<" & \"0\";"<<endl;
				//o<<tab<<"ResultException   <= exception_post_normalization;"<<endl;
				o<<tab<<"ResultSign        <= X("<<wEX + wFX<<") xor Y("<<wEY + wFY<<");"<<endl;
			}
			else{
				/* Rounding needed */
				/* Check case when resut is in the middle of two FP numbers */

				// generate two 1000...0000 and 0100....000 strings 
				str1.str(""); str2.str(""); //init
				str1<<"\"1" << zero_generator( wFX+wFY+2 - (1+wFR) -1, 1);
				str2<<"\"01"<< zero_generator( wFX+wFY+2 - (1+wFR) -1 -1, 1);
				
				// the check string represent the LSBits of the significand product to the right of the rounding bit
				o<<tab<<"check_string <= (significand_product ("<< wFX+wFY+2 - (1+wFR) -1<<" downto 0) xor "<<str1.str()<<") when normalization_selector='1' else"<<endl;
				o<<tab<<"                ( (\"0\" & significand_product ("<< wFX+wFY+2 - (1+wFR) -1-1<<" downto 0)) xor "<<str2.str()<<");"<<endl;

				// In the middle of 2 FP numbers <=> LSBits to the right of the rounding bit are 100...00
				o<<tab<<"process(check_string)"<<endl;
				o<<tab<<"begin"<<endl;
				o<<tab<<tab<<"      if ( "<< zero_generator( wFX+wFY+2 - (1+wFR) , 0) <<" = check_string ) then"<<endl; 
				o<<tab<<tab<<"         between_fp_numbers <= '1';"<<endl;
				o<<tab<<tab<<"      else "<<endl;
				o<<tab<<tab<<"         between_fp_numbers <= '0';"<<endl;
				o<<tab<<tab<<"      end if;"<<endl;
				o<<tab<<tab<<"end process; "<<endl;

				// the signal exponent_post_normalization2 has wEX + 1 bits
				
				// the reunion signal is a concatenation of 0 & Significand      
				o<<tab<<"reunion_signal <= (\"0\" & significand_product("<< wFX+wFY<<" downto "<<wFX+wFY - (wFR) <<" )) when normalization_selector='1' else"<<endl;
				o<<tab<<"                  (\"0\" & significand_product("<< wFX+wFY-1<<" downto "<<wFX+wFY - (wFR) -1<<" ));"<<endl; 
				               
				// LSB_of_result_significand = LSB of the wFR part of the significand_product 
				o<<tab<<"LSB_of_result_significand <= significand_product("<< wFX+wFY+1 - wFR<<" ) when  normalization_selector='1' else"<<endl;
				o<<tab<<"                             significand_product("<< wFX+wFY+1 - wFR -1<<" );"<<endl;              
				   		     			
				// add 1 to the reunited signal for rounding & normalize purposes
				o<<tab<<"reunion_signal_post_addition <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2 + wFR<<");"<<endl;
				
				o<<tab<<"between_fp_numbers_result_significand <= reunion_signal_post_addition when LSB_of_result_significand = '1' else"<<endl;
				o<<tab<<"                                         reunion_signal;"<<endl;    		     
					
				o<<tab<<"reunion_signal_post_rounding <= reunion_signal_post_addition when between_fp_numbers='0' else"<<endl;
				o<<tab<<"                                between_fp_numbers_result_significand;"<<endl;
				                   
				
				o<<tab<<"exponent_post_rounding <= exponent_post_normalization + (CONV_STD_LOGIC_VECTOR(0,"<<wER+1<<") & reunion_signal_post_rounding("<<wFR+1<<"));"<<endl;
				
				                        
				//update exception        
				 
				o<<tab<<"with exponent_post_rounding("<< wER+1 <<" downto "<< wER <<") select"<<endl;		
				o<<tab<<"ResultException <= exception_post_normalization when \"00\","<<endl;
				o<<tab<<"                            \"10\"              when \"01\", "<<endl;
				o<<tab<<"                            \"00\"              when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"              when others;"<<endl;						               
                                         
				o<<tab<<"ResultExponent    <= exponent_post_rounding("<<wER -1<<" downto "<<0<<");"<<endl;  
				o<<tab<<"ResultSignificand <= \"1\" & reunion_signal_post_rounding("<<wFR<<" downto 1);"<<endl;
				o<<tab<<"ResultSign        <= X("<<wEX + wFX<<") xor Y("<<wEY + wFY<<");"<<endl;
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

TestIOMap FPMultiplier::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("X"));
	tim.add(*get_signal_by_name("Y"));
	tim.add(*get_signal_by_name("ResultException"));
	tim.add(*get_signal_by_name("ResultSign"));
	tim.add(*get_signal_by_name("ResultExponent"));
	tim.add(*get_signal_by_name("ResultSignificand"));
	return tim;
}

void FPMultiplier::fillTestCase(mpz_class a[])
{
	/* Get I/Os */
	mpz_class &svx = a[0];
	mpz_class &svy = a[1];
	mpz_class &svexc = a[2];
	mpz_class &svsgn = a[3];
	mpz_class &svexp = a[4];
	mpz_class &svfra = a[5];

	/* Compute result */
	FloFP x(wEX, wFX), y(wEY, wFY), r(wER, wFR, normalized);
	x = svx; y = svy;
	r = x * y;

	svexc = r.getExceptionSignalValue();
	svsgn = r.getSignSignalValue();
	// Exponent and fraction are not defined for zero, inf or NaN
	if (r.getExceptionSignalValue() == 1)
	{
		svexp = r.getExponentSignalValue();
		svfra = r.getFractionSignalValue();
	}
}

