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
#include "FPNumber.hpp"

using namespace std;
extern vector<Operator*> oplist;


FPMultiplier::FPMultiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm) :
	Operator(target), wEX_(wEX), wFX_(wFX), wEY_(wEY), wFY_(wFY), wER_(wER), wFR_(wFR) {

	int i, j;
	ostringstream name, synch, synch2;

	setOperatorName();
	setOperatorType();

	/* set if operator outputs a normalized_ result */
	if (norm == 0) 
		normalized_ = false;
	else
		normalized_ = true;

	
	/* Set up the IO signals */
	/* Inputs: 2b(Exception) + 1b(Sign) + wEX_ bits (Exponent) + wFX_ bits(Fraction) */
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	addFPInput ("X", wEX_, wFX_);
	addFPInput ("Y", wEY_, wFY_);

	// FIXME remove all these output ports in the normalized case. Or, add them to all the operators
	addSignal ("ResultExponent"   , wER_    );  
	addSignal ("ResultSignificand", wFR_ + 1);
	addSignal ("ResultException"  , 2      );
	addSignal ("ResultSign"       , 1      ); 
	if(normalized_) {
		addFPOutput ("R"   , wER_, wFR    );  
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		   
	/* Instantiate an intmult_iplier -> multiply the significands */
	intmult_ = new IntMultiplier(target, wFX_+1, wFY_+1);
	oplist.push_back(intmult_);

	//Setup the attibute cocerning the intmult_iplier pipeline depth 
	intMultPipelineDepth_ = intmult_->getPipelineDepth();
	
	/* Common signals for all versions */

	addSignal("significandX", wFX_+1);   //The Significands for the input signals X and Y
	addSignal("significandY", wFY_+1);
	addSignal("significand_product", wFX_ +wFY_ +2);   

	addSignal("exponentX", wEX_);	//The signals for the exponents of X and Y
	addSignal("exponentY", wEY_);
	addSignal("bias",wEX_+2); //bias to be substracted from the sum of exponents  

	addSignal("exception_selector",4);//signal that selects the case for the exception

	
	if (isSequential())
	{
	/* Common signals for the sequential version */
		
		/* Unregistered signals */
		addSignal("significand_synch", wFX_ + wFY_ + 2); 
		addSignal("exponent_synch", wEX_+2);
		addSignal("exception_synch",2);
		addSignal("sign_synch", 1); 
		
		/* Registered signals */
		//addDelaySignal("Exponents_Sum_Pre_Bias_Substraction",  wEX_+2 );
		//addDelaySignal("Exponents_Sum_Post_Bias_Substraction", wEX_+2 );
		//addDelaySignal("Exception_Reg_0",2); 
		//addDelaySignal("Exception_Reg_1",2);
		
		addDelaySignal("Exponents_Sum_Pre_Bias_Substraction",  wEX_+2 );
		addDelaySignal("Exponents_Sum_Post_Bias_Substraction", wEX_+2 );
		addDelaySignal("Exception_Reg_0",2); 
		addDelaySignal("Exception_Reg_1",2);
		
		
		
			/* Synchronization registers */
			
			//when integer multiplication has more than 2 pipeline levels
			for (i=0;i<=intMultPipelineDepth_-3;i++)	{
				synch.str(""); synch2.str("");
				synch<<"Exponent_Synchronization_"<<i;
				
				//addDelaySignal(synch.str(), wEX_+2 );
				addDelaySignal(synch.str(), wEX_+2 );
				
				synch2<<"Exception_Synchronization_"<<i;
				//addDelaySignal(synch2.str(), 2 );
				addDelaySignal(synch2.str(), 2 );
			}
			//when the integer multiplication has less than 2 levels
			for (i=0;i<=1-intMultPipelineDepth_;i++)	{
				synch.str("");
				synch<<"Int_Synchronization_"<<i;
				//addDelaySignal(synch.str(), wFX_+wFY_+2);
				addDelaySignal(synch.str(), wFX_+wFY_+2);				
			}	
			addDelaySignal("Sign_Synchronization", 1, max(2,intMultPipelineDepth_));

		if (normalized_){
			/* The sequential normalized_ version */
		
			addSignal("normalization_selector",1);
			addSignal("exponent_post_normalization",2+wEX_);
			addSignal("exception_post_normalization",2);
			
			if (1+wFR_ < wFX_+wFY_+2) {
				
				//addDelaySignal("exponent_synch2",wEX_+2);
				addDelaySignal("exponent_synch2",wEX_+2);
				
				addDelaySignal("exception_synch2",2);
				
				//addDelaySignal("significand_synch2",wFX_ + wFY_ + 2);
				addDelaySignal("significand_synch2",wFX_ + wFY_ + 2);
				
				
				addDelaySignal("sign_synch2",1);
				addDelaySignal("normalization_selector2",1);
		
				addSignal("between_fp_numbers",1);	
				addSignal("LSB_of_result_significand_out",1);	
				addSignal("between_fp_numbers_out",1);
				addSignal("sign_synch2_out",1);
				addSignal("exception_synch2_out",2);
				addSignal("exponent_synch2_out",wER_+2);
				addSignal("exponent_synch2_out_post_rounding",wER_+2);
				
						
				addSignal("check_string",wFX_+wFY_+2 - (1+wFR_) );
				addSignal("reunion_signal",2+wFR_);
				addSignal("reunion_signal_out",2+wFR_);
				addSignal("LSB_of_result_significand",1); 
				addSignal("between_fp_numbers_result_significand",2+wFR_);
				addSignal("reunion_signal_post_addition",2+wFR_); 
				addSignal("reunion_signal_post_rounding",2+wFR_); 
				
				
				//parameters setup
				
				bool status = target->suggestSubaddSize(additionChunkWidth_, 2+ wFR_);
				if (!status)
					cout<<"Frequency report:"<<endl;
					cout<<tab<<"WARNING: Desired frequency cannot be reached !";
			
			
				reunionSignalWidth_ = 2+ wFR_;
				reunionSignalParts_ = int ( ceil( double(reunionSignalWidth_)/double(additionChunkWidth_)));
				additionLastChunkWidth_ = reunionSignalWidth_ - (reunionSignalParts_ - 1)*additionChunkWidth_;
							
				
				
				if (reunionSignalParts_>1){
					for (j=1; j<=reunionSignalParts_;j++)	
						for (i=1;i<=reunionSignalParts_;i++){	
							name.str("");
							name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
							if (i!=reunionSignalParts_){
								//addDelaySignalBus(name.str(), additionChunkWidth_ + 1 ,1);
								addDelaySignalBus(name.str(), additionChunkWidth_ + 1 ,1);
							}
							else
							{
								//addDelaySignalBus(name.str(), additionLastChunkWidth_,1 );	
								addDelaySignalBus(name.str(), additionLastChunkWidth_,1 );	
							
							}
		        }
				
					/*
					addDelaySignal("reunion_signal_level",2+wFR_,reunionSignalParts_);
					addDelaySignal("LSB_of_result_significand_level",1,reunionSignalParts_);		
					addDelaySignal("between_fp_numbers_level",1,reunionSignalParts_);			
					addDelaySignal("sign_synch2_level",1,reunionSignalParts_);
					addDelaySignal("exception_synch2_level",2,reunionSignalParts_);
					addDelaySignal("exponent_synch2_level",2+wER_,reunionSignalParts_);
					*/
					
					addDelaySignal("reunion_signal_level",2+wFR_,reunionSignalParts_);
					addDelaySignal("LSB_of_result_significand_level",1,reunionSignalParts_);		
					addDelaySignal("between_fp_numbers_level",1,reunionSignalParts_);			
					addDelaySignal("sign_synch2_level",1,reunionSignalParts_);
					addDelaySignal("exception_synch2_level",2,reunionSignalParts_);
					addDelaySignal("exponent_synch2_level",2+wER_,reunionSignalParts_);
					
					
				}
			}
		
		}else{
		/* The sequential non-normalized_ version */
			
			
			if ((wFX_+wFY_+2 > wFR_+1))
			{
				/*
				addDelaySignal("exponent_synch2",wEX_+2);
				addDelaySignal("exception_synch2",2);
				addDelaySignal("significand_synch2",wFX_ + wFY_ + 2);
				addDelaySignal("sign_synch2",1);
				*/
				
				addDelaySignal("exponent_synch2",wEX_+2);
				addDelaySignal("exception_synch2",2);
				addDelaySignal("significand_synch2",wFX_ + wFY_ + 2);
				addDelaySignal("sign_synch2",1);
						
				addSignal("between_fp_numbers",1);	
				addSignal("LSB_of_result_significand_out",1);	
				addSignal("between_fp_numbers_out",1);
				addSignal("sign_synch2_out",1);
				addSignal("exception_synch2_out",2);
				addSignal("exponent_synch2_out",wER_+2);
				addSignal("exponent_synch2_out_post_rounding",wER_+2);
										
				addSignal("check_string",wFX_+wFY_+2 - (1+wFR_) );
				
			
			
				bool status = target->suggestSubaddSize(additionChunkWidth_, 3+ wFR_);
				if (!status)
					cout<<"Frequency report:"<<endl;
					cout<<tab<<"WARNING: Desired frequency cannot be reached !";
						
				reunionSignalWidth_ = 3 + wFR_;
				reunionSignalParts_ = int ( ceil( double(reunionSignalWidth_)/double(additionChunkWidth_)));
				additionLastChunkWidth_ = reunionSignalWidth_ - (reunionSignalParts_ - 1)*additionChunkWidth_;
				
				if (reunionSignalParts_>1){
					for (j=1; j<=reunionSignalParts_;j++)	
					for (i=1;i<=reunionSignalParts_;i++){	
						name.str("");
						name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
						if (i!=reunionSignalParts_)
						{
							//addDelaySignal(name.str(), additionChunkWidth_ + 1 );
							addDelaySignal(name.str(), additionChunkWidth_ + 1 );
						}
						else{
							//addDelaySignal(name.str(), additionLastChunkWidth_ );
							addDelaySignal(name.str(), additionLastChunkWidth_ );		
						}
	        }
				
					/*
					addDelaySignal("reunion_signal_level",3+wFR_,reunionSignalParts_);
					addDelaySignal("LSB_of_result_significand_level",1,reunionSignalParts_);		
					addDelaySignal("between_fp_numbers_level",1,reunionSignalParts_);			
					addDelaySignal("sign_synch2_level",1,reunionSignalParts_);
					addDelaySignal("exception_synch2_level",2,reunionSignalParts_);
					addDelaySignal("exponent_synch2_level",2+wER_,reunionSignalParts_);
					*/
					
					addDelaySignal("reunion_signal_level",3+wFR_,reunionSignalParts_);
					addDelaySignal("LSB_of_result_significand_level",1,reunionSignalParts_);		
					addDelaySignal("between_fp_numbers_level",1,reunionSignalParts_);			
					addDelaySignal("sign_synch2_level",1,reunionSignalParts_);
					addDelaySignal("exception_synch2_level",2,reunionSignalParts_);
					addDelaySignal("exponent_synch2_level",2+wER_,reunionSignalParts_);					
					
					
				}

				addSignal("reunion_signal"                       ,1+(1+wFR_)+1);
				addSignal("reunion_signal_out"                   ,3+wFR_);
				addSignal("LSB_of_result_significand"            ,1    ); 
				addSignal("between_fp_numbers_result_significand",3+wFR_);
				addSignal("reunion_signal_post_addition"         ,3+wFR_); 
				addSignal("reunion_signal_post_rounding"         ,3+wFR_);
						
				//addSignal("temp_exp_concat_fract", 1 + wEX_ + wFX_ + wFY_ + 2);
			}
		}
	}
	else{
		 
		/* Signals for the combinational version */
		addSignal("normalization_selector",1);
		addSignal("exponent_post_normalization",2+wER_);
		addSignal("exception_post_normalization",2);
				
		addSignal("Exponents_Sum_Pre_Bias_Substraction", wER_+2 );
		addSignal("Exponents_Sum_Post_Bias_Substraction", wER_+2 );   
		addSignal("Exception_Reg_0",2); 
		addSignal("Exception_Reg_1",2);
		if ((1 + wFR_ < wFX_ + wFY_ + 2)) {
			addSignal("check_string",wFX_+wFY_+2 - (1+wFR_) );
			addSignal("reunion_signal",2+wFR_);
			addSignal("LSB_of_result_significand",1); 
			addSignal("between_fp_numbers_result_significand",2+wFR_);
			addSignal("reunion_signal_post_addition",2+wFR_); 
			addSignal("reunion_signal_post_rounding",2+wFR_); 
			addSignal("exponent_post_rounding",wER_+2);
			addSignal("between_fp_numbers",1);
		}
		
	}	

	
	//set pipeline depth 
	if (!isSequential())
		setPipelineDepth(0);
	else
		if (normalized_)
			if (1+wFR_ >= wFX_+wFY_+2)   
				setPipelineDepth(max(2,intMultPipelineDepth_));
			else
				if (reunionSignalParts_>1)
					setPipelineDepth(1 + max(2,intMultPipelineDepth_) + reunionSignalParts_ );
				else
					setPipelineDepth(1 + max(2,intMultPipelineDepth_));
		else
			if (1+wFR_>=2+wFX_+wFY_)
				setPipelineDepth(max(2,intMultPipelineDepth_));
			else
				if (reunionSignalParts_==1)
					setPipelineDepth(max(2,intMultPipelineDepth_) + 1);
				else
					setPipelineDepth(max(2,intMultPipelineDepth_) + 1 + reunionSignalParts_);

}

FPMultiplier::~FPMultiplier() {
}

void FPMultiplier::setOperatorName(){
	/* The name has the format: FPMultiplier_wEX__wFX__wEY__wFY__wER__wFR_ where: wEX_ = width of X exponenet and wFX_ = width for the fractional part of X */
	ostringstream name;
	name.str("");
	name<<"FPMultiplier_"<<wEX_<<"_"<<wFX_<<"_"<<wEY_<<"_"<<wFY_<<"_"<<wER_<<"_"<<wFR_; 
	uniqueName_ = name.str(); 
}


void FPMultiplier::outputVHDL(std::ostream& o, std::string name) {
  	
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;
	
	int bias_val=int(pow(double(2),double(wEX_-1)))-1;
	int i, j; 

	licence(o,"Bogdan Pasca (2008)");
	Operator::stdLibs(o);
		
	o<<endl<<endl;	
	outputVHDLEntity(o);

	newArchitecture(o,name);
	
	intmult_->outputVHDLComponent(o);
	outputVHDLSignalDeclarations(o);	  
	o<<endl;

	beginArchitecture(o);
	
	
	if (isSequential()){
		/* common code for both normalized_ and non-normalized_ version */
		outputVHDLRegisters(o); o<<endl;
		
		/* Exponent Handling */
			o<<tab<<"exponentX <= X("<<wEX_ + wFX_ -1<<" downto "<<wFX_<<");"<<endl; 
			o<<tab<<"exponentY <= Y("<<wEY_ + wFY_ -1<<" downto "<<wFY_<<");"<<endl<<endl;
			//Add exponents -> put sum into register
			o<<tab<<"Exponents_Sum_Pre_Bias_Substraction <= (\"00\" & exponentX) + (\"00\" & exponentY);"<<endl; //wEX_+2 bits	[*Optimization point - pipelining*]
			o<<tab<<"bias <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wEX_+2<<");"<<endl; 
			//substract the bias value from the exponents' sum
			o<<tab<<"Exponents_Sum_Post_Bias_Substraction <= Exponents_Sum_Pre_Bias_Substraction_d - bias;"<<endl;//wEX_ + 2   
		
		/* Significand Handling */
			//build significands 	 	
			o<<tab<< "significandX <= \"1\" & X("<<(wFX_-1)<<" downto 0);"<<endl;
			o<<tab<< "significandY <= \"1\" & Y("<<(wFY_-1)<<" downto 0);"<<endl<<endl;
			//multiply significands 
			o<<tab<< "int_multiplier_component: " << intmult_->getOperatorName() << endl;
			o<<tab<< "      port map ( X => significandX, " << endl; //wFX_+1 bits (1 bit is the hidden 1)
			o<<tab<< "                 Y => significandY, " << endl; //wFY_+1 bits
			o<<tab<< "                 R => significand_product, " << endl; //wFX_+wFY_+2 bits
			o<<tab<< "                 clk => clk, " << endl;
			o<<tab<< "                 rst => rst " << endl;
			o<<tab<< "               );" << endl<<endl;

		/* Exception Handling */
			//Concatenate the exception bits of the two input numbers
			o<<tab<<"exception_selector <= X("<<wEX_ + wFX_ +2<<" downto "<<wEX_ + wFX_ + 1<<") & Y("<<wEY_ + wFY_ +2<<" downto "<<wEY_ + wFY_ +1<<");"<<endl<<endl;
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
									
		/* synchronization barrier (intmult_iplication | Exponent Addition | Exception Computation | Sign Computation */
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
		/* Cases: 		
			0. when intmult_iplication | Exponent Addition are allready synchronized	
			1. when registers are added to exponent addition
			2. when registers are added to integer multiplication
		*/	
	  
		if (intMultPipelineDepth_==2){
			// when intmult_iplication | Exponent Addition are allready synchronized	 
			o<<tab<<"exponent_synch    <= Exponents_Sum_Post_Bias_Substraction_d;"<<endl;
			o<<tab<<"significand_synch <= significand_product;"<<endl;
			o<<tab<<"exception_synch   <= Exception_Reg_1_d;"<<endl;
		}
	   else{	
			if (intMultPipelineDepth_ > 2){
			/* when registers are added to exponent addition */
				//connect the first synchronization reg to the output of the last reg from exponent addition 
				o<<tab<<"Exponent_Synchronization_0  <= Exponents_Sum_Post_Bias_Substraction_d;"<<endl;
				o<<tab<<"Exception_Synchronization_0 <= Exception_Reg_1_d;"<<endl; 									

				for (i=1;i<=intMultPipelineDepth_-3;i++){
					o<<tab<<"Exponent_Synchronization_"<<i<<"   <= "<<"Exponent_Synchronization_"<<i-1<<"_d;"<<endl;		
					o<<tab<<"Exception_Synchronization_"<<i<<" <= "<<"Exception_Synchronization_"<<i-1<<"_d;"<<endl;
				}
					o<<tab<<"exponent_synch    <= Exponent_Synchronization_"<<	intMultPipelineDepth_-3 <<"_d;"<<endl;
					o<<tab<<"exception_synch   <= Exception_Synchronization_"<<intMultPipelineDepth_-3 <<"_d;"<<endl;
					o<<tab<<"significand_synch <= significand_product;"<<endl;			
			}
			else{
				//add registers to the output of intmult_iplier so that exponent addition and significand multiplication become synchronized
				o<<tab<<"Int_Synchronization_0 <= significand_product;"<<endl;
				for (i=1;i<=1-intMultPipelineDepth_;i++){
					synch1.str(""); synch2.str("");
					synch2<<"Int_Synchronization_"<<i;
					synch1<<"Int_Synchronization_"<<i-1<<"_d";
					o<<tab<<synch2.str()<<"<="<<synch1.str()<<";"<<endl;		
				}
				o<<tab<<"exponent_synch    <= Exponents_Sum_Post_Bias_Substraction_d;"<<endl;
				o<<tab<<"exception_synch   <= Exception_Reg_1_d;"<<endl;
				o<<tab<<"significand_synch <= Int_Synchronization_"<< 1-intMultPipelineDepth_<<"_d;"<<endl;
			}
		}
		//sign synchronization
		o<<tab<<"Sign_Synchronization <= X("<<wEX_ + wFX_<<") xor Y("<<wEY_ + wFY_<<");"<<endl;
		o<<tab<<"sign_synch <= "<<delaySignal("Sign_Synchronization",max(2,intMultPipelineDepth_))<<";"<<endl<<endl;
	
	
		if (normalized_==true){	
			/* normalize result */
			
			//assign selector signal to MSB of signif. multiplic 
			o<<tab<<"normalization_selector <= significand_synch ("<<wFX_+wFY_+1<<");"<<endl;

			o<<tab<<"exponent_post_normalization <= (exponent_synch ) + (CONV_STD_LOGIC_VECTOR(0 ,"<<wEX_+1<<") & normalization_selector);"<<endl;

			/* result rounding */
			//check is rounding is needed
			if (1+wFR_ >= wFX_+wFY_+2) {
				/* =>no rounding needed - possible padding */
				o<<tab<<"with exponent_post_normalization("<< wER_+1 <<" downto "<< wER_ <<") select"<<endl;		
				o<<tab<<"exception_post_normalization <= exception_synch when \"00\","<<endl;
				o<<tab<<"                            \"10\"             when \"01\", "<<endl;
				o<<tab<<"                            \"00\"             when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"             when others;"<<endl;						
				
				//make the outputs
				o<<tab<<"ResultExponent     <= exponent_post_normalization("<<wEX_ - 1<<" downto 0);"<<endl;  
				o<<tab<<"ResultSignificand  <= significand_synch & "<<zeroGenerator(1+wFR_ - (wFX_+wFY_+2) , 0)<<" when normalization_selector='1' else"<<endl;
				o<<tab<<"                      significand_synch("<<wFX_+wFY_<<" downto 0) & "<<zeroGenerator(1+wFR_ - (wFX_+wFY_+2) , 0)<<" & \"0\";"<<endl;
				o<<tab<<"ResultException    <= exception_post_normalization;"<<endl;
				o<<tab<<"ResultSign         <= sign_synch;"<<endl;

			}
			else{ 
				o<<tab<<"exception_post_normalization <= exception_synch;"<<endl;
				//check if in the middle of two FP numbers
				//generate two xor strings
				str1.str(""); str2.str(""); //init
				str1<<"\"1"<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) -1, 1);
				str2<<"\"01"<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) -1-1, 1);

				zeros1.str("");
				zeros1<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) , 0);

				o<<tab<<"check_string <= (significand_synch ("<< wFX_+wFY_+2 - (1+wFR_) -1<<" downto 0) xor "<<str1.str()<<") when normalization_selector='1' else"<<endl;
				o<<tab<<"               ( (\"0\" & significand_synch ("<< wFX_+wFY_+2 - (1+wFR_) -1-1<<" downto 0)) xor "<<str2.str()<<");"<<endl;

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


				o<<tab<<"reunion_signal <= (\"0\" & significand_synch2_d("<< wFX_+wFY_<<" downto "<<wFX_+wFY_ - (wFR_) <<" )) when normalization_selector2_d='1' else"<<endl;
				o<<tab<<"                  (\"0\" & significand_synch2_d("<< wFX_+wFY_-1<<" downto "<<wFX_+wFY_ - (wFR_) -1<<" ));"<<endl; 
              
				o<<tab<<"LSB_of_result_significand <= significand_synch2_d("<< wFX_+wFY_+1 - wFR_<<" ) when normalization_selector2_d='1' else"<<endl;
				o<<tab<<"                             significand_synch2_d("<< wFX_+wFY_+1 - wFR_ -1<<" );"<<endl;              
				
				/* pipeline the possible carry_in propagation caused by rounding */
				
				//connect the first part of the pipeline structure to the signals
				if (reunionSignalParts_ != 1){
					for (j = 1; j <= reunionSignalParts_; j++){
						if (j == 1)
							o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <= (\"0\" & reunion_signal("<<j*additionChunkWidth_-1<<" downto "<<(j-1)*additionChunkWidth_<<")) + '1';"<<endl;
						else
							if (j!=reunionSignalParts_)
								o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <= (\"0\" & reunion_signal("<<j*additionChunkWidth_-1<<" downto "<<(j-1)*additionChunkWidth_<<"));"<<endl;	
							else
								o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <=  reunion_signal("<<(j-1)*additionChunkWidth_ + additionLastChunkWidth_ -1<<" downto "<<(j-1)*additionChunkWidth_<<");"<<endl;	
					}
					
					//put together the pipeline
					for (j=1; j<=reunionSignalParts_-1;j++)	
						for (i=1;i<=reunionSignalParts_;i++)
							if ((i<=j)||(i>j+1))
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
							else
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<additionChunkWidth_<<") ;"<<endl; 

					//put the result back together
					ostringstream reunion_signal_post_1_addition;
					for (i = reunionSignalParts_; i > 0; i--) 
						if ( (i!=1) && (i!=reunionSignalParts_) )
							reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunionSignalParts_<<"_Reg_"<<i<<"_d("<<additionChunkWidth_ -1<<" downto 0) &" ;
						else
							if (i==reunionSignalParts_)
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunionSignalParts_<<"_Reg_"<<i<<"_d("<<additionLastChunkWidth_ -1<<" downto 0) & ";
							else
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunionSignalParts_<<"_Reg_"<<i<<"_d("<<additionChunkWidth_ -1<<" downto 0)" ;
					
				
					//propagate the rest of the signals reunionSignalParts_ steps
					
					o<<tab<<"reunion_signal_level <=reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_out <= "<<delaySignal("reunion_signal_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"LSB_of_result_significand_level <= LSB_of_result_significand;"<<endl;
					o<<tab<<"LSB_of_result_significand_out <= "<<delaySignal("LSB_of_result_significand_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"between_fp_numbers_level <= between_fp_numbers;"<<endl;
					o<<tab<<"between_fp_numbers_out <= "<<delaySignal("between_fp_numbers_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"sign_synch2_level <= sign_synch2_d;"<<endl;
					o<<tab<<"sign_synch2_out <= "<<delaySignal("sign_synch2_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"exception_synch2_level <= exception_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out <= "<<delaySignal("exception_synch2_level", reunionSignalParts_)<<";" <<endl;		
														
					o<<tab<<"reunion_signal_post_addition <= "<< reunion_signal_post_1_addition.str()<<";"<< endl;
					
					o<<tab<<"exponent_synch2_level <= exponent_synch2_d;"<<endl;
					o<<tab<<"exponent_synch2_out <= "<<delaySignal("exponent_synch2_level", reunionSignalParts_)<<";" <<endl;						
				
				}   
				else{
					/* when the carry propagation on the reunion signal is not pipelined */
					//add 1 to the reunited signal for rounding & normalize purposes
					o<<tab<<"reunion_signal_out            <= reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_post_addition  <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2+ wFR_<<");"<<endl;
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
				                        
				
				
				o<<tab<<"exponent_synch2_out_post_rounding <= exponent_synch2_out + (CONV_STD_LOGIC_VECTOR(0,"<<wEX_+1<<") & reunion_signal_post_rounding("<<1+wFR_<<"));"<<endl;    
				
				                        
				//update exception 
				o<<tab<<"with exponent_synch2_out_post_rounding("<<wER_+1<<" downto "<< wER_ <<") select"<<endl;                       
				o<<tab<<"ResultException   <= exception_synch2_out      when \"00\","<<endl;
				o<<tab<<"                            \"10\"             when \"01\", "<<endl;
				o<<tab<<"                            \"00\"             when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"             when others;"<<endl;						
			                
				                        
				                        
				//update exception                       
				//o<<tab<<"ResultException   <= exception_synch2_out when reunion_signal_post_rounding("<<wEX_ + wFR_+1<<")='0' else"<<endl;
				//o<<tab<<"                     \"10\";"<<endl;                                                   

				o<<tab<<"ResultExponent    <= exponent_synch2_out_post_rounding("<<wER_-1<<" downto "<<0<<");"<<endl;  
				o<<tab<<"ResultSignificand <= \"1\" & reunion_signal_post_rounding("<<wFR_<<" downto 1);"<<endl;
				o<<tab<<"ResultSign        <= sign_synch2_out;"<<endl;
				o<<tab<<"R <= ResultException & ResultSign & ResultExponent & ResultSignificand("<<wFR_-1<<" downto 0);" << endl;
			}
		}
		else { //TODO
		/* Non-normalized_ case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			if (1+wFR_ < wFX_+wFY_+2){
				//check if in the middle of two FP numbers
				//generate two xor strings
				str1.str(""); 
				str1<<"\"1"<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) -1, 1);
				
				zeros1.str("");
				zeros1<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) , 0);

				o<<tab<<"check_string <= (significand_synch ("<< wFX_+wFY_+2 - (1+wFR_) -1<<" downto 0) xor "<<str1.str()<<"); "<<endl;
				
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
				
				o<<tab<<"reunion_signal <= (\"0\"  & significand_synch2_d("<< wFX_+wFY_+1<<" downto "<<wFX_+wFY_+1 - (wFR_+1) <<" ));"<<endl;; 
              
				o<<tab<<"LSB_of_result_significand <= significand_synch2_d("<< wFX_+wFY_+1 - wFR_<<" );"<<endl;
				
				/* pipeline the possible carry_in propagation caused by rounding */
				
				//connect the first part of the pipeline structure to the signals
				if (reunionSignalParts_ != 1){
					for (j = 1; j <= reunionSignalParts_; j++){
						if (j == 1)
							o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <= (\"0\" & reunion_signal("<<j*additionChunkWidth_-1<<" downto "<<(j-1)*additionChunkWidth_<<")) + '1';"<<endl;
						else
							if (j!=reunionSignalParts_)
								o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <= (\"0\" & reunion_signal("<<j*additionChunkWidth_-1<<" downto "<<(j-1)*additionChunkWidth_<<"));"<<endl;	
							else
								o<<"Last_Addition_Level_"<<1<<"_Reg_"<<j<<" <=  reunion_signal("<<(j-1)*additionChunkWidth_ + additionLastChunkWidth_ -1<<" downto "<<(j-1)*additionChunkWidth_<<");"<<endl;	
					}
					
					//put together the pipeline
					for (j=1; j<=reunionSignalParts_-1;j++)	
						for (i=1;i<=reunionSignalParts_;i++)
							if ((i<=j)||(i>j+1))
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
							else
								o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<additionChunkWidth_<<") ;"<<endl; 

					//put the result back together
					ostringstream reunion_signal_post_1_addition;
					for (i = reunionSignalParts_; i > 0; i--) 
						if ( (i!=1) && (i!=reunionSignalParts_) )
							reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunionSignalParts_<<"_Reg_"<<i<<"_d("<<additionChunkWidth_ -1<<" downto 0) &" ;
						else
							if (i==reunionSignalParts_)
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunionSignalParts_<<"_Reg_"<<i<<"_d("<<additionLastChunkWidth_ -1<<" downto 0) & ";
							else
								reunion_signal_post_1_addition<<"Last_Addition_Level_"<<reunionSignalParts_<<"_Reg_"<<i<<"_d("<<additionChunkWidth_ -1<<" downto 0)" ;
					
				
					//propagate the rest of the signals reunionSignalParts_ steps
					
					o<<tab<<"reunion_signal_level <=reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_out <= "<<delaySignal("reunion_signal_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"LSB_of_result_significand_level <= LSB_of_result_significand;"<<endl;
					o<<tab<<"LSB_of_result_significand_out <= "<<delaySignal("LSB_of_result_significand_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"between_fp_numbers_level <= between_fp_numbers;"<<endl;
					o<<tab<<"between_fp_numbers_out <= "<<delaySignal("between_fp_numbers_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"sign_synch2_level <= sign_synch2_d;"<<endl;
					o<<tab<<"sign_synch2_out <= "<<delaySignal("sign_synch2_level", reunionSignalParts_)<<";" <<endl;		
					
					o<<tab<<"exception_synch2_level <= exception_synch2_d;"<<endl;
					o<<tab<<"exception_synch2_out <= "<<delaySignal("exception_synch2_level", reunionSignalParts_)<<";" <<endl;		
														
					o<<tab<<"reunion_signal_post_addition <= "<< reunion_signal_post_1_addition.str()<<";"<< endl;
					
					o<<tab<<"exponent_synch2_level <= exponent_synch2_d;"<<endl;
					o<<tab<<"exponent_synch2_out <= "<<delaySignal("exponent_synch2_level", reunionSignalParts_)<<";" <<endl;						
				}   
				else{
					/* when the carry propagation on the reunion signal is not pipelined */
					//add 1 to the reunited signal for rounding & normalize purposes
					o<<tab<<"reunion_signal_out            <= reunion_signal;"<<endl;
					o<<tab<<"reunion_signal_post_addition  <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2+wFR_<<");"<<endl;
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
				
				o<<tab<<"exponent_synch2_out_post_rounding <= exponent_synch2_out + (CONV_STD_LOGIC_VECTOR(0,"<<wEX_+1<<") & reunion_signal_post_rounding("<<2+wFR_<<"));"<<endl;    
				
				                        
				//update exception 
				o<<tab<<"with exponent_synch2_out_post_rounding("<<wER_+1<<" downto "<< wER_ <<") select"<<endl;                       
				o<<tab<<"ResultException   <= exception_synch2_out      when \"00\","<<endl;
				o<<tab<<"                            \"10\"             when \"01\", "<<endl;
				o<<tab<<"                            \"00\"             when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"             when others;"<<endl;						
			

				o<<tab<<"ResultExponent    <= exponent_synch2_out_post_rounding("<<wEX_ -1<<" downto "<<0<<");"<<endl;  
					
				o<<tab<<"ResultSignificand <= (\"11\" & reunion_signal_post_rounding("<<wFR_-1<<" downto 1)) when reunion_signal_out("<<wFR_+1<<" downto "<<wFR_<<")=\"11\" else"<<endl;
				o<<tab<<"					  reunion_signal_post_rounding("<<wFR_+1<<" downto 1);"<<endl;	
				o<<tab<<"ResultSign        <= sign_synch2_out;"<<endl;

			}
			else{
			/* No rounding is needed. If wFR_+1>wFX_+wFY_+2 then the ResultSignificand must be padded with 0s to the right */
				if (wFX_+wFY_+2 == wFR_+1)
					o<<tab<<"ResultSignificand <= significand_synch;"<<endl;
				else 	
					o<<tab<<"ResultSignificand <= significand_synch & "<<zeroGenerator((wFR_+1) - (wFX_+wFY_+2), 0 )<<";"<<endl;
					
				o<<tab<<"with exponent_synch("<<wER_+1<<" downto "<< wER_ <<") select"<<endl;                       
				o<<tab<<"ResultException   <= exception_synch  when \"00\","<<endl;
				o<<tab<<"                            \"10\"    when \"01\", "<<endl;
				o<<tab<<"                            \"00\"    when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"    when others;"<<endl;							
					
				o<<tab<<"ResultExponent  <= exponent_synch("<<wER_-1<<" downto 0 "<< ");"<<endl;
				
				o<<tab<<"ResultSign <= sign_synch;"<<endl;
			}				
    
			
		}//end else not normalized_
       
	}//end if pipelined
	else
	{  
		/* the combinational version */
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
		/* Exponent Handling */
			// Fetch exponents from inputs 
			o<<tab<<"exponentX <= X("<<wEX_ + wFX_ -1<<" downto "<<wFX_<<");"<<endl; 
			o<<tab<<"exponentY <= Y("<<wEY_ + wFY_ -1<<" downto "<<wFY_<<");"<<endl<<endl;
			//Add exponents and put sum into register
			o<<tab<<"Exponents_Sum_Pre_Bias_Substraction <= (\"00\" & exponentX) + (\"00\" & exponentY);"<<endl; //wEX_+2 bits
			o<<tab<<"bias <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wER_+2<<");"<<endl; 
			//substract the bias value from the exponent sum
			o<<tab<<"Exponents_Sum_Post_Bias_Substraction <= Exponents_Sum_Pre_Bias_Substraction - bias;"<<endl;//wEX_ + 2   

		/* Significand Handling */
			//fetch and build significands by adding a "1" in MSB position 
			o<<tab<< "significandX <= \"1\" & X("<<(wFX_-1)<<" downto 0);"<<endl;
			o<<tab<< "significandY <= \"1\" & Y("<<(wFY_-1)<<" downto 0);"<<endl<<endl;
			//multiply significands 
			o<<tab<< "int_multiplier_component: " << intmult_->getOperatorName() << endl;
			o<<tab<< "      port map ( X => significandX, " << endl; //wFX_+1 bits (1 bit is the hidden 1)
			o<<tab<< "                 Y => significandY, " << endl; //wFY_+1 bits
			o<<tab<< "                 R => significand_product " << endl; //wFX_+wFY_+2 bits
			o<<tab<< "               );" << endl<<endl;

		/* Exception bits Handling */
			//Concatenate the exception bits of the two input numbers
			o<<tab<<"exception_selector <= X("<<wEX_ + wFX_ +2<<" downto "<<wEX_ + wFX_ + 1<<") & Y("<<wEY_ + wFY_ +2<<" downto "<<wEY_ + wFY_ +1<<");"<<endl<<endl;
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
			o<<tab<<"normalization_selector <= significand_product ("<<wFX_+wFY_+1<<");"<<endl;
			// add the selector bit to the exponent in order to normalize the result 
			o<<tab<<"exponent_post_normalization <= Exponents_Sum_Post_Bias_Substraction + (CONV_STD_LOGIC_VECTOR(0 ,"<<wER_+1<<") & normalization_selector);"<<endl;
			o<<tab<<"exception_post_normalization <= Exception_Reg_0;"<<endl;			

		/* Rounding */
			if (1+wFR_ >= wFX_+wFY_+2) {
				/* No rounding needed */
				/* Possible 0 padding to the length of the fractional part of the result */
				o<<tab<<"with exponent_post_normalization("<< wER_+1 <<" downto "<< wER_ <<") select"<<endl;		
				o<<tab<<" ResultException <= exception_post_normalization when \"00\","<<endl;
				o<<tab<<"                            \"10\"                           when \"01\", "<<endl;
				o<<tab<<"                            \"00\"                           when \"11\"|\"10\","<<endl;
				o<<tab<<"                             \"11\"                          when others;"<<endl;
											
				// Assign the architecture outputs
				o<<tab<<"ResultExponent    <= exponent_post_normalization("<<wEX_ - 1<<" downto 0);"<<endl;  
				o<<tab<<"ResultSignificand <= significand_product & "<< zeroGenerator(1+wFR_ - (wFX_+wFY_+2) , 0) <<" when normalization_selector='1' else"<<endl;
				o<<tab<<"                     significand_product("<<wFX_+wFY_<<" downto 0) & "<<zeroGenerator(1+wFR_ - (wFX_+wFY_+2) , 0)<<" & \"0\";"<<endl;
				//o<<tab<<"ResultException   <= exception_post_normalization;"<<endl;
				o<<tab<<"ResultSign        <= X("<<wEX_ + wFX_<<") xor Y("<<wEY_ + wFY_<<");"<<endl;
			}
			else{
				/* Rounding needed */
				/* Check case when resut is in the middle of two FP numbers */

				// generate two 1000...0000 and 0100....000 strings 
				str1.str(""); str2.str(""); //init
				str1<<"\"1" << zeroGenerator( wFX_+wFY_+2 - (1+wFR_) -1, 1);
				str2<<"\"01"<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) -1 -1, 1);
				
				// the check string represent the LSBits of the significand product to the right of the rounding bit
				o<<tab<<"check_string <= (significand_product ("<< wFX_+wFY_+2 - (1+wFR_) -1<<" downto 0) xor "<<str1.str()<<") when normalization_selector='1' else"<<endl;
				o<<tab<<"                ( (\"0\" & significand_product ("<< wFX_+wFY_+2 - (1+wFR_) -1-1<<" downto 0)) xor "<<str2.str()<<");"<<endl;

				// In the middle of 2 FP numbers <=> LSBits to the right of the rounding bit are 100...00
				o<<tab<<"process(check_string)"<<endl;
				o<<tab<<"begin"<<endl;
				o<<tab<<tab<<"      if ( "<< zeroGenerator( wFX_+wFY_+2 - (1+wFR_) , 0) <<" = check_string ) then"<<endl; 
				o<<tab<<tab<<"         between_fp_numbers <= '1';"<<endl;
				o<<tab<<tab<<"      else "<<endl;
				o<<tab<<tab<<"         between_fp_numbers <= '0';"<<endl;
				o<<tab<<tab<<"      end if;"<<endl;
				o<<tab<<tab<<"end process; "<<endl;

				// the signal exponent_post_normalization2 has wEX_ + 1 bits
				
				// the reunion signal is a concatenation of 0 & Significand      
				o<<tab<<"reunion_signal <= (\"0\" & significand_product("<< wFX_+wFY_<<" downto "<<wFX_+wFY_ - (wFR_) <<" )) when normalization_selector='1' else"<<endl;
				o<<tab<<"                  (\"0\" & significand_product("<< wFX_+wFY_-1<<" downto "<<wFX_+wFY_ - (wFR_) -1<<" ));"<<endl; 
				               
				// LSB_of_result_significand = LSB of the wFR_ part of the significand_product 
				o<<tab<<"LSB_of_result_significand <= significand_product("<< wFX_+wFY_+1 - wFR_<<" ) when  normalization_selector='1' else"<<endl;
				o<<tab<<"                             significand_product("<< wFX_+wFY_+1 - wFR_ -1<<" );"<<endl;              
				   		     			
				// add 1 to the reunited signal for rounding & normalize purposes
				o<<tab<<"reunion_signal_post_addition <= reunion_signal + CONV_STD_LOGIC_VECTOR(1, "<< 2 + wFR_<<");"<<endl;
				
				o<<tab<<"between_fp_numbers_result_significand <= reunion_signal_post_addition when LSB_of_result_significand = '1' else"<<endl;
				o<<tab<<"                                         reunion_signal;"<<endl;    		     
					
				o<<tab<<"reunion_signal_post_rounding <= reunion_signal_post_addition when between_fp_numbers='0' else"<<endl;
				o<<tab<<"                                between_fp_numbers_result_significand;"<<endl;
				                   
				
				o<<tab<<"exponent_post_rounding <= exponent_post_normalization + (CONV_STD_LOGIC_VECTOR(0,"<<wER_+1<<") & reunion_signal_post_rounding("<<wFR_+1<<"));"<<endl;
				
				                        
				//update exception        
				 
				o<<tab<<"with exponent_post_rounding("<< wER_+1 <<" downto "<< wER_ <<") select"<<endl;		
				o<<tab<<"ResultException <= exception_post_normalization when \"00\","<<endl;
				o<<tab<<"                            \"10\"              when \"01\", "<<endl;
				o<<tab<<"                            \"00\"              when \"11\"|\"10\","<<endl;
				o<<tab<<"                            \"11\"              when others;"<<endl;						               
                                         
				o<<tab<<"ResultExponent    <= exponent_post_rounding("<<wER_ -1<<" downto "<<0<<");"<<endl;  
				o<<tab<<"ResultSignificand <= \"1\" & reunion_signal_post_rounding("<<wFR_<<" downto 1);"<<endl;
				o<<tab<<"ResultSign        <= X("<<wEX_ + wFX_<<") xor Y("<<wEY_ + wFY_<<");"<<endl;
				o<<tab<<"R <= ResultException & ResultSign & ResultExponent & ResultSignificand("<<wFR_-1<<" downto 0);" << endl;
			}

	}//end the combinational part   
  
	o<< "end architecture;" << endl << endl;
}




// FIXME: the following is only functional for a correctly rounded multiplier.
// The old code for the non-normalized case is commented out below, just in case.
void FPMultiplier::emulate(TestCase * tc)
{
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");

	/* Compute correct value */
	FPNumber fpx(wEX_, wFX_), fpy(wEY_, wFY_);
	fpx = svX;
	fpy = svY;
	mpfr_t x, y, r;
	mpfr_init2(x, 1+wFX_);
	mpfr_init2(y, 1+wFY_);
	mpfr_init2(r, 1+wFR_); 
	fpx.getMPFR(x);
	fpy.getMPFR(y);
	mpfr_mul(r, x, y, GMP_RNDN);

	// Set outputs 
	FPNumber  fpr(wER_, wFR_, r);
	mpz_class svR = fpr.getSignalValue();
	tc->addExpectedOutput("R", svR);

	// clean up
	mpfr_clears(x, y, r, NULL);
}


#if 0
void FPMultiplier::fillTestCase(mpz_class a[])
{
	/* Get I/Os */
	mpz_class &svx = a[0];
	mpz_class &svy = a[1];
	mpz_class &svr = a[2];
	/* Compute result */
	FPNumber fpx(wEX_, wFX_), fpy(wEY_, wFY_);
	fpx = svx; fpy = svy;

	mpfr_t x, y, r;
	mpfr_init2(x, wFX_+1);
	mpfr_init2(y, wFY_+1);
	mpfr_init2(r,   wFX_ + wFY_ + 2); // r will hold an exact product
	fpx.getMPFR(x);
	fpy.getMPFR(y);
	if(normalized_) {
		mpfr_mul(r, x, y, GMP_RNDN);
		FPNumber fpr(wER_, wFR_, r);
		svr = fpr.getSignalValue();
	}
	else {
		throw std::string("fillTestCase not yet implememented for non-normalized results");
#if 0
		//This is the relevant SWW code purged from FPNumber& FPNumber::operator=(mpfr_t mp_)
		// It should be implemented here instead of polluting all of FPNumber.

/* SwW: Slipper when Wet
 * When the FPMultiplier is set not to normalise results,
 * it outputs significants and exponents which have little 
 * logic when viewed from a software side. Unfortunately,
 * in order to have a good test bench, we have to emulate
 * those behaviours in software. Whenever you see this tag
 * expect hard-to-understand or ilogic code.
 *
 * What happens is that the binary fraction comma is placed
 * after the second bit. Exponents are
 * simply added together. Sometimes the first bit
 * is 1 (the result is „overnormalised”), other times
 * it is 0 (the result is normalized). We will detect this
 * „condition” by seeing if the exponent of the normalized result
 * is higher than the sum of the exponents of the operands.
 *
 * We will do the following. We will always store numbers as
 * normalized numbers (1.mantissa), exponents like FPMultiplier
 * and when necessary (i.e. mustunnormalize==true)), we will do
 * rounding with one bit „faster” and shift the significant in
 * getFractionSignalValue().
 */

	if (!normalise && mustAddLeadingZero)
	{
		/* SwW: we need to round with one bit earlier */ 
		mpfr_mul_2si(mp, mp, wF-1, GMP_RNDN);
		mpfr_get_z(mantissa.get_mpz_t(), mp, GMP_RNDN);
		mantissa = mantissa << 1;
	}
	/* SwW: exponent is smaller in the „normalised” case */
	if (!normalise && !mustAddLeadingZero)
		exp--;

#endif
	}

	mpfr_clears(x, y, r, 0, NULL);

}
#endif // end commented out code
