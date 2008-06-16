/*
 * An integer multiplier for FloPoCo
 *
 * Authors : Bogdan Pasca and Florent de Dinechin
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
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"
#include "IntMultiplier.hpp"

using namespace std;
extern vector<Operator*> oplist;

/** 
 * The constructor of the IntMultiplier class
 * @param target argument of type Target containing the data for which this operator will be optimized
 * @param wInX integer argument representing the width in bits of the input X 
 * @param wInY integer argument representing the width in bits of the input Y
 **/
IntMultiplier:: IntMultiplier(Target* target, int wInX, int wInY) :
	Operator(target), wInX(wInX), wInY(wInY), wOut(wInX + wInY){
 
	int depth, i, j;
	ostringstream name, nameH, nameL;

	//gets information if the multiplier is or not to be pipelined
	if(target->is_pipelined())
		set_sequential();
	else
		set_combinatorial();  

	/* Name Setup procedure
	 *  The name has the format: IntMultiplier_wInX_wInY
	 *  wInX = width of the X input
	 *  wInY = width of the Y input
	 */  
	name.str("");
	name <<"IntMultiplier_"<<wInX<<"_"<<wInY;
	unique_name = name.str(); 
	
	/* Set up the IO signals
	 * X and Y have wInX and wInY bits respectively 
	 * R has wOut bits where wOut = (wInX + WInY) bits
	 */
	  add_input ("X", wInX);
	  add_input ("Y", wInY);
	  add_output("R", wOut);
	 
  
	
	if(is_sequential()) 
	{
	// set up the variables and parameters for the sequential version	
		
		
		//given the input widths and the specific target this method suggests the chunk widths for X and Y
		bool test = target->suggest_submult_size(multiplier_width_X, multiplier_width_Y, wInX, wInY);
		cout<<"Frequency report (multiplications):"<<endl;
		if (test)
			cout<<tab<<"Frequency can be reached = "<<"YES"<<" suggested sizes: X="<<multiplier_width_X<<" Y="<<multiplier_width_Y<<endl;
		else
			cout<<tab<<"WARNING: Frequency can be reached = "<<"NO"<<" suggested sizes: X="<<multiplier_width_X<<" Y="<<multiplier_width_Y<<endl;
		
		// decide in how many parts should the numbers be split depending on the recommended sizes
		partsX = ((wInX % multiplier_width_X)==0) ? (wInX / multiplier_width_X) : (wInX / multiplier_width_X+1);
		partsY = ((wInY % multiplier_width_Y)==0) ? (wInY / multiplier_width_Y) : (wInY / multiplier_width_Y+1);
		
		//reverse
		if (multiplier_width_X < multiplier_width_Y){
			i = wInY;
			wInY = wInX;
			wInX = i;
			
			i = multiplier_width_X;
			multiplier_width_X = multiplier_width_Y;
			multiplier_width_Y = i;
			
			i = partsX;
			partsX = partsY;
			partsY = i;
			
			reverse = true;
		}else 
			reverse = false;
		
		number_of_zerosX = partsX * multiplier_width_X - wInX;
		number_of_zerosY = partsY * multiplier_width_Y - wInY;
		
		// signals that will contain the inputs X and Y but padded with 0 so that 
		// their width becomes a multiple of the target multiplier width
		add_signal("i_X", partsX * multiplier_width_X);
		add_signal("i_Y", partsY * multiplier_width_Y);

		//declare the signals needed to split X 
		for (i=1;i<=partsX;i++){
			name.str("");
			name <<"X_"<<i;
			//add_registered_signal_with_sync_reset(name.str(), multiplier_width_X);
			add_registered_signal(name.str(), multiplier_width_X);  
		}

		//declare the signals needed to split Y 
		for (i=1;i<=partsY;i++){
			name.str("");
			name <<"Y_"<<i;
			//add_registered_signal_with_sync_reset(name.str(), multiplier_width_Y); 
			add_registered_signal(name.str(), multiplier_width_Y); 
		}

		//declare the registers needed to store the partial products
		for (i=1;i<=partsY;i++)
			for (j=1;j<=partsX;j++){
				name.str("");;
				name <<"Y_"<<i<<"_X_"<<j;
				//add_registered_signal_with_sync_reset(name.str(), multiplier_width_X + multiplier_width_Y);
				add_registered_signal(name.str(), multiplier_width_X + multiplier_width_Y);
			} 

		if (!((partsX==partsY) && (partsX==1)))
		{
		//this signals are not needed when simple 1blockx1block multiplications are performed
		
			//declare the high/low part decomposition signals
			for (i=1;i<=partsY;i++)
				for (j=1;j<=partsX;j++){
					nameH.str("");; nameL.str("");;
					nameH <<"Y_"<<i<<"_X_"<<j<<"_H";
					nameL <<"Y_"<<i<<"_X_"<<j<<"_L";
					add_signal(nameH.str(), multiplier_width_X); //the largest of the two
					add_signal(nameL.str(), multiplier_width_X);
				} 

			//signals to regroup all the high parts and all the low parts
			for (i=1;i<=partsY;i++)  {
				nameH.str("");; nameL.str("");;
				nameH <<"High_"<<i;
				nameL <<"Low_"<<i;
				add_signal(nameH.str(), partsX * multiplier_width_X);
				add_signal(nameL.str(), partsX * multiplier_width_X);
			}
		}  
		
		
		if (partsY>1) {
		//declare the signals which form the two trees (high and low)
			
			/* Use a state of the art adder to take care of pipelining the long additions */
			intadd2 = new IntAdder(target, 1 + partsX * multiplier_width_X);
			oplist.push_back(intadd2);
			
			
			for (i=1;i<=partsY;i++){
				 name.str("");
				 name<<"Buffer_H_"<<i;
				 //add_delay_signal(name.str(), partsX * multiplier_width_X,intadd2->pipeline_depth());
				 add_delay_signal_bus_no_reset(name.str(), partsX * multiplier_width_X,intadd2->pipeline_depth());
			}
			
			
			
			//for the low and high tree
			for (j=0; j<partsY;j++)
				for (i=1;i<=partsY-j;i++){
					nameH.str("");; nameL.str("");;
					nameL <<"L_Level_"<<j<<"_Reg_"<<i;
					nameH <<"H_Level_"<<j<<"_Reg_"<<i;
					if ((i==1)||(i==2)){
						
						add_signal(nameL.str(), 1 + partsX * multiplier_width_X);
						add_signal(nameH.str(), 1 + partsX * multiplier_width_X);
						ostringstream temp;
						if (i==1)
						{
							temp.str("");
							temp<<"L_Level_"<<j<<"_summand";
							add_signal(temp.str(),1 + partsX * multiplier_width_X);
							temp.str("");
							temp<<"H_Level_"<<j<<"_summand";
							add_signal(temp.str(),1 + partsX * multiplier_width_X);
						}
						
					}
					else{
						//add_delay_signal(nameL.str(), 1 + partsX * multiplier_width_X,intadd2->pipeline_depth());
						add_delay_signal_bus_no_reset(nameL.str(), 1 + partsX * multiplier_width_X,intadd2->pipeline_depth());
						
						//add_delay_signal(nameH.str(), 1 + partsX * multiplier_width_X,intadd2->pipeline_depth());
						add_delay_signal_bus_no_reset(nameH.str(), 1 + partsX * multiplier_width_X,intadd2->pipeline_depth());
					}
					
					
				}

			
			//one more registers for the Low part for synchronization 
			//add_delay_signal("Low1", 1 + partsX * multiplier_width_X,intadd2->pipeline_depth());
			add_delay_signal_bus_no_reset("Low1", 1 + partsX * multiplier_width_X,intadd2->pipeline_depth());
			
			//one more signal for the high part
			add_signal("High1", 1 + partsX * multiplier_width_X);

			
			//partsY! registers to hold and compute the low part of the result in the 
			//pipeline levels
			for (j=1; j<=partsY;j++)
				for (i=1;i<=j;i++){
					name.str("");;
					name<<"PartialBits_Level_"<<j<<"_Reg_"<<i;
					//add_delay_signal(name.str(), multiplier_width_Y + 1,intadd2->pipeline_depth() );
					add_delay_signal_bus_no_reset(name.str(), multiplier_width_Y + 1,intadd2->pipeline_depth() );
			 }

			/* Use a state of the art adder to take care of pipelining the last addition */
			intadd1 = new IntAdder(target, partsX * multiplier_width_X);
			oplist.push_back(intadd1);


			
			//add_delay_signal("PartialBits_Reg",partsY * multiplier_width_Y,intadd1->pipeline_depth());
			add_delay_signal_bus_no_reset("PartialBits_Reg",partsY * multiplier_width_Y,intadd1->pipeline_depth());  
			
			add_signal("temp_result",  partsX * multiplier_width_X );
			add_signal("partial_bits", partsY * multiplier_width_Y );
			add_signal("full_result",  partsX * multiplier_width_X + partsY * multiplier_width_Y);

		}//end if partsY>1
		else
		{
			
			if (partsX>1)
			{
				intadd = new IntAdder(target, partsX * multiplier_width_X);
				oplist.push_back(intadd);

				//Setup the attibute cocerning the IntMultiplier pipeline 
				IntAddPipelineDepth = intadd->pipeline_depth();

				cout<<endl<<" IntAddPipelineDepth = "<<IntAddPipelineDepth<<endl;

				add_signal("addition_term", partsX  * multiplier_width_X);

				add_signal("addition_result", partsX  * multiplier_width_X);
				
				//add_delay_signal("delayed_bits",multiplier_width_Y, IntAddPipelineDepth);
				add_delay_signal_bus_no_reset("delayed_bits",multiplier_width_Y, IntAddPipelineDepth);
			}
	
			add_signal("temp_result", partsX  * multiplier_width_X + partsY * multiplier_width_Y);
 		}

	
		// Set up the pipeline
		if (partsY==1)
			if (partsX==1)
				depth=2;
			else
				depth=2 + IntAddPipelineDepth;
		else
			depth= 2 + partsY*intadd2->pipeline_depth() + intadd1->pipeline_depth();
 
		set_pipeline_depth(depth);
  
	}//end if sequential
}



/**
 * Int Multiplier destructor
 */
IntMultiplier::~IntMultiplier() {
}





/**
 * Method belonging to the Operator class overloaded by the IntMultiplier class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void IntMultiplier::output_vhdl(std::ostream& o, std::string name) {
	
	//variables
	ostringstream zeros, names, nameH, nameL, concatH, concatL,
					  first_summand, second_summand, the_bits;

	int i, j, k;
		
	Licence(o,"Bogdan Pasca and Florent de Dinechin (2008)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	
	new_architecture(o,name);
	

	if ((partsY==1)&&(partsX>1))
		intadd->output_vhdl_component(o);//output the IntAdder component
	else
		if (!((partsX==partsY) && (partsX==1))){
			intadd1->output_vhdl_component(o);
			intadd2->output_vhdl_component(o);
		}
	
	
	output_vhdl_signal_declarations(o);
	begin_architecture(o);
 
	if (is_sequential()) { 
	
		output_vhdl_registers(o); 
		
		//pad inputs with zeros
		pad_inputs(o);
		
		//split X into chunks and assign the chunks to X_1 ...X_partsX 
		split(o,"X",partsX, multiplier_width_X);			  
		split(o,"Y",partsY, multiplier_width_Y);
		
		//make all multiplications concurrently
		do_multiplications(o);
		
		if (!((partsX==partsY) && (partsX==1)))
		{	
			//decompose the results into high and low part
			decompose_low_high(o);
			//Regroup all the high parts together and all the low parts together
			regroup_low_high(o);
		}
		
		new_line(o);
		//pass to the adder tree
		if (partsY==1) {
			if (partsX>1){
				
				o<<"addition_term <= ("<< zero_generator(multiplier_width_Y,0)<<" & Low_1("<< partsX*multiplier_width_X-1<<" downto "<< multiplier_width_Y<<"));"<<endl;
				map_adder(o, "addition_term", "High_1", "addition_result");
				
				//some of the bits of Low_1 do no participate in the addtion. These bits together with the result of the addition need to be merged to form the final result.
				//=> they need to be delayed with a number of clock cycles equal to the pipeline levels of the addition
				o<<tab<<"delayed_bits <= Low_1("<<multiplier_width_Y-1<<" downto 0 )"<<";"<<endl;
				
				//the concatenation of the addtion result and the delayed bits form the result
				o<<tab<<"temp_result <= addition_result & "<<get_delay_signal_name("delayed_bits", IntAddPipelineDepth) <<";"<<endl;
			}else{
				//when both X and Y have 1 part the result is simply the product of these two parts
				o<<tab<<"temp_result <= Y_1_X_1_d;"<<endl;
			}
			
			//output the result
			o<<tab<<"R <= temp_result("<<partsX*multiplier_width_X + partsY*multiplier_width_Y-1 -number_of_zerosX - number_of_zerosY<<" downto 0);"<<endl;
		}
		else //(the number of parts of Y is >1)
		{ 

			new_line(o);
			//connect the low part of the multiplication results to the low tree
			connect_low(o);
						
			new_line(o);
			//the high part of the multiplication needs to be delayed one clock cycle before it will be inserted into the adder tree
			high_to_buffer(o);	
			
			new_line(o);
			//the buffers need to be connected to the high tree
			connect_buffers(o);
				
			//Connect the low part tohether
			link_low_adder_structure(o);
			 
			//Connect the HIGH part tohether 
			link_high_adder_structure(o);
			
			//connect the partsY! registers which compute the low part of the result
			connect_partial_bits(o);
			
			//pipeline the addition which gives the gigh part of the multiplication result;
	  	pipeline_addition(o); //the result of the addition is called temp_result
	  
			delay_partial_bits(o);
						
			o<<tab<<"full_result <= temp_result & partial_bits;"<<endl; 
			o<<tab<<"R <= full_result("<<  partsX*multiplier_width_X + partsY*multiplier_width_Y - 1 -number_of_zerosX - number_of_zerosY<<" downto 0);"<<endl<<endl;
		  
		}
	}else{ 
		//the combinational version
		o << tab<<"R <= X * Y ;" << endl;
	}
	
	end_architecture(o);
}



/**
 * A zero generator method which takes as input two arguments and returns a string of zeros with quotes as stated by the second argurment
 * @param[in] n		    integer argument representing the number of zeros on the output string
 * @param[in] margins	integer argument determining the position of the quotes in the output string. The options are: -2= no quotes; -1=left quote; 0=both quotes 1=right quote
 * @return returns a string of zeros with the corresonding quotes given by margins
 **/
string IntMultiplier::zero_generator(int n, int margins)
{
	ostringstream left,full, right, zeros;
	int i;

	//generate the zeros without quotes
	for (i=1; i<=n;i++)
		zeros<<"0";

	//generate the 3 strings with left quote, both and right quote
	left<<"\""<<zeros.str();
	full<<left.str()<<"\"";
	right<<zeros.str()<<"\"";

	//margins determines which of the strings are outputed
	switch(margins){
		case -2: return zeros.str(); break;
		case -1: return left.str();  break;
		case  0: return full.str();  break;
		case  1: return right.str(); break;
		default: return full.str();
	}
}


/**
 * A method which pads the inputs with 0. 
 * @param[in,out] o - the stream to which the new line will be added
 **/
void IntMultiplier::pad_inputs(std::ostream& o)
{
	if (reverse == true){
		o<<tab<< "i_Y <= "<< zero_generator(number_of_zerosY,0)<<" & X;"<<endl;
		o<<tab<< "i_X <= "<< zero_generator(number_of_zerosX,0)<<" & Y;"<<endl;
		new_line(o);
	}else{
		o<<tab<< "i_X <= "<< zero_generator(number_of_zerosX,0)<<" & X;"<<endl;
		o<<tab<< "i_Y <= "<< zero_generator(number_of_zerosY,0)<<" & Y;"<<endl;
		new_line(o);
	}		
}


/**
 * A method which splits the input given by name into [parts] parts, each part having the width part_width
 * @param[in,out]	o 			- the stream to which the code will be inputed
 * @param[in]		name		- the name of the input string
 * @param[in]		parts		- the number of parts the input will be split in
 * @param[in]		part_width	- the width of the part	
 **/
void IntMultiplier::split(std::ostream& o, std::string name, int parts, int part_width)
{
int i;
ostringstream temp;
	
	for (i=1;i<=parts;i++){
		temp.str("");
		temp << name <<"_"<<i;
		o<<tab<< temp.str() << " <= i_"<<name<<"("<< i * part_width -1 <<" downto "<< (i-1) * part_width<<");"<<endl;
	}
}


/**
 * A does the multiplications of the parts X and Y were splitted in
 * @param[in,out] o - the stream to which the multiplications are added
 **/
void IntMultiplier::do_multiplications(std::ostream& o)
{
int i,j;
ostringstream names;	
	for (i=1;i<=partsY;i++)
		for (j=1;j<=partsX;j++){  
			names.str("");;
			names <<"Y_"<<i<<"_X_"<<j;
			o<<tab<< names.str() << " <= " << "Y_"<< i <<"_d * " << "X_"<<j << "_d;"<<endl;
		} 
	new_line(o);
}


/**
 * Decompose the multiplication of two chunks into high and low parts
 * @param[in,out] o - the stream to which the decomposition is written
 **/
void IntMultiplier::decompose_low_high(std::ostream& o)
{
int i,j;
ostringstream names, nameH, nameL;
	for (i=1;i<=partsY;i++)
		for (j=1;j<=partsX;j++){
			names.str("");; nameH.str("");; nameL.str("");;
			names <<"Y_"<<i<<"_X_"<<j<<"_d";
			nameH <<"Y_"<<i<<"_X_"<<j<<"_H";
			nameL <<"Y_"<<i<<"_X_"<<j<<"_L";
			o<<tab<< nameH.str() << " <= " << names.str() <<"("<< multiplier_width_X + multiplier_width_Y - 1 <<" downto " << multiplier_width_Y <<");"<<endl;
			o<<tab<< nameL.str() << " <= " <<zero_generator(multiplier_width_X - multiplier_width_Y ,0)<<" & "<< names.str() <<"("<< multiplier_width_Y - 1 <<" downto " << 0 <<");"<<endl;
			} 
}

/**
 * Regroup the decomposition result into vectors low and high vectors
 * @param[in,out] o - the stream to which the regrouping is written
 **/
void IntMultiplier::regroup_low_high(std::ostream& o)
{
int i,j;
ostringstream nameL,nameH, concatH, concatL;
	for (i=1;i<=partsY;i++){
		nameH.str(""); nameL.str(""); concatH.str(""); concatL.str("");
		nameH <<"High_"<<i;
		nameL <<"Low_"<<i;
		//create the concatenation strings for low and high parts
		for (j=partsX;j>=1;j--){
			if (j>1){
				concatL<<"Y_"<<i<<"_X_"<<j<<"_L & ";
				concatH<<"Y_"<<i<<"_X_"<<j<<"_H & ";
			}else{
			  concatL<<"Y_"<<i<<"_X_"<<j<<"_L";
			  concatH<<"Y_"<<i<<"_X_"<<j<<"_H";
			}
		}
		o<<tab<< nameL.str() << "  <= " << concatL.str() << ";"<<endl;
		o<<tab<< nameH.str() << " <= " << concatH.str() << ";"<<endl;
	}
	new_line(o);
}


/**
 * Adder mapping
 * @param[in,out] 	o 			- the stream to which the mapping is written
 * @param[in]		left_term	- the name of the left term of the addition as a string
 * @param[in]		right_term	- the name of the right term of the addition as a string
 * @param[in]		result		- the name of the result of the addition
 **/
void IntMultiplier::map_adder(std::ostream& o,std::string left_term, std::string right_term, std::string result )
{
	o<<tab<< "int_adder_component: " << intadd->unique_name << endl;
	o<<tab<< "  port map ( X => "<< left_term << ", " << endl; 
	o<<tab<< "             Y => "<< right_term<< ", " << endl; 
	o<<tab<< "             Cin => '0' ," << endl;
	o<<tab<< "             R => "<< result<< ", " << endl; 
	o<<tab<< "             clk => clk, " << endl;
	o<<tab<< "             rst => rst " << endl;
	o<<tab<< "               );" << endl<<endl;
}

/**
 * Connect the low part of the multiplication to the low part of the adder structure
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::connect_low(std::ostream& o)
{
int i;
	for(i=1;i<=partsY;i++)
		o<<tab<<"L_Level_0_Reg_"<<i<<" <= \"0\" & Low_"<<i<<";"<<endl;
}


/**
 * Connect the high part of the multiplication to a buffer
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::high_to_buffer(std::ostream& o)
{
int i;
	for(i=1;i<=partsY;i++)
		o<<tab<<"Buffer_H_"<< i << " <= High_" <<i<<";"<<endl;
}

/**
 * Connect the buffer to the high part of the adder structure
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::connect_buffers(std::ostream& o)
{
int i;
	for(i=1;i<=partsY;i++)
	{
	ostringstream temp;
	temp<<"Buffer_H_"<<i;
		o<<tab<<"H_Level_0_Reg_"<<i<<" <= \"0\" & "<<get_delay_signal_name(temp.str(),intadd2->pipeline_depth())<<";"<<endl;
	}
}

/**
 * Link together the components of the low part of the adder structure
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::link_low_adder_structure(std::ostream& o)
{
int i,j;
ostringstream first_summand, second_summand;
	for (i=0;i<partsY-1;i++) {
		for (j=1;j<=partsY-i;j++) {
			if (j==1) {
				//concatenate zeros in front of the first register, add it with the second register and place the result in the next level's first register
				first_summand.str(""); second_summand.str("");
				first_summand<<"("<< zero_generator(multiplier_width_Y ,0)<<" & L_Level_"<<i<<"_Reg_"<<j<<"("<<partsX*multiplier_width_X<<" downto "<<multiplier_width_Y<<"))"; 
				second_summand<<"L_Level_"<<i<<"_Reg_"<<j+1;
				
				o<<tab<<"L_Level_"<<i+1<<"_summand <= "<<first_summand.str()<<";"<<endl;								
				o<<tab<< "int_adder_component_low_"<<i<<": " << intadd2->unique_name << endl;
				o<<tab<< "  port map ( X => L_Level_"<<i+1<<"_summand , " << endl; 
				o<<tab<< "             Y => "<<second_summand.str()<< ", " << endl; 
				o<<tab<< "             Cin => '0'," << endl;
				o<<tab<< "             R => L_Level_"<<i+1<<"_Reg_"<<j<<", " << endl; 
				o<<tab<< "             clk => clk, " << endl;
				o<<tab<< "             rst => rst " << endl;
				o<<tab<< "               );" << endl<<endl;

				//o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j<<"<="<< first_summand.str() << " + " << second_summand.str()<<";"<< endl;
			}
			else
				if (j==2) {} //do nothing (the second register on the level is always processed with the first
					else //for the rest of the registers just propagate their contents to the next level
					{
						ostringstream tmpSig;
						tmpSig<<"L_Level_"<<i<<"_Reg_"<<j;
						o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<< get_delay_signal_name(tmpSig.str(),intadd2->pipeline_depth())<<";"<<endl;
					}	
		}
	}

	//the zeros + [the high bits of the last register of the low part] 
	o<<tab<<"Low1 <= "<<zero_generator(multiplier_width_Y,0)<<" & L_Level_"<<partsY-1<<"_Reg_1("<<partsX * multiplier_width_X<<" downto "<<multiplier_width_Y<<")"<< ";"<<endl;
}


/**
 * Link together the components of the high part of the adder structure
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::link_high_adder_structure(std::ostream& o)
{
int i,j;
ostringstream first_summand, second_summand;
	for (i=0;i<partsY-1;i++){
		for (j=1;j<=partsY-i;j++){
			if (j==1){
				//concatenate zeros in front of the first register, add it with the second register and place the result in the next level's first register
				first_summand.str("");; second_summand.str("");;
				first_summand<<"("<< zero_generator(multiplier_width_Y,0)<<" & H_Level_"<<i<<"_Reg_"<<j<<"("<<partsX*multiplier_width_X<<" downto "<<multiplier_width_Y<<"))"; 
				second_summand<<"H_Level_"<<i<<"_Reg_"<<j+1<<"";
				
				o<<tab<<"H_Level_"<<i+1<<"_summand <= "<<first_summand.str()<<";"<<endl;								
				o<<tab<< "int_adder_component_high_"<<i<<": " << intadd2->unique_name << endl;
				o<<tab<< "  port map ( X => H_Level_"<<i+1<<"_summand , " << endl; 
				o<<tab<< "             Y => "<<second_summand.str()<< ", " << endl; 
				o<<tab<< "             Cin => '0'," << endl;
				o<<tab<< "             R => H_Level_"<<i+1<<"_Reg_"<<j<<", " << endl; 
				o<<tab<< "             clk => clk, " << endl;
				o<<tab<< "             rst => rst " << endl;
				o<<tab<< "               );" << endl<<endl;
					
				//o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j<<"<="<< first_summand.str() << " + " << second_summand.str()<<";"<< endl;
			}
			else
				if (j==2) {} //do nothing (the second register on the level is always processed with the first
					else //for the rest of the registers just propagate their contents to the next level
					{
					ostringstream tmpSig;
						tmpSig<<"H_Level_"<<i<<"_Reg_"<<j;
						o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<< get_delay_signal_name(tmpSig.str(),intadd2->pipeline_depth())<<";"<<endl;
								
					//	o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<<"H_Level_"<<i<<"_Reg_"<<j<<"_d;"<<endl;
					}
		}
	}

	o<<tab<<"High1 <= H_Level_"<<partsY-1<<"_Reg_1"<<";"<<endl;
}

/**
 * The partial bits structure outputs the low part of the result by using short additions
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::connect_partial_bits(std::ostream& o)
{
int i,j;
	for (j=1; j<=partsY;j++)
		for (i=1;i<=j;i++){
			if ((j==1) and (i==1) )
				o<<tab<<"PartialBits_Level_1_Reg_1 <= \"0\" & L_Level_0_Reg_1("<< multiplier_width_Y -1 <<" downto 0);"<< endl;
			else{
				//just propagate the registers from 1 to j-1
				if (i<j)
				{
					ostringstream temp;
					temp<<"PartialBits_Level_"<<j-1<<"_Reg_"<<i;
					o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i <<"<= "<<get_delay_signal_name(temp.str(),intadd2->pipeline_depth())<<";"<<endl;
				}
				else //if i=j compute the jth register
				{
					ostringstream temp;
					temp<<"PartialBits_Level_"<<j-1<<"_Reg_"<<i-1;
					o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i<<"<= "<<"(\"0\" & L_Level_"<<j-1<<"_Reg_1("<<multiplier_width_Y-1<<" downto 0)) + "
																														<<"(\"0\" & H_Level_"<<j-2<<"_Reg_1("<<multiplier_width_Y-1<<" downto 0)) + "
																														<< get_delay_signal_name(temp.str(),intadd2->pipeline_depth())<<"("<<multiplier_width_Y<<");"<<endl; 
				}
			}	
		}
}


/**
 * The addition which gives the high part of the product between X and Y is pipelined to obtain better performance
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::pipeline_addition(std::ostream& o)
{
ostringstream temp;
temp<<"PartialBits_Level_"<< partsY <<"_Reg_"<< partsY;
		o<<tab<< "int_adder_component: " << intadd1->unique_name << endl;
		o<<tab<< "  port map ( X => "<< get_delay_signal_name("Low1",intadd2->pipeline_depth())<< "("<<partsX * multiplier_width_X -1<<" downto 0"<<") , " << endl; 
		o<<tab<< "             Y => High1("<<partsX * multiplier_width_X -1<<" downto 0"<<"), " << endl; 
		o<<tab<< "             Cin => "<<get_delay_signal_name(temp.str(),intadd2->pipeline_depth()) <<"("<< multiplier_width_Y <<")," << endl;
		o<<tab<< "             R => temp_result, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;
}

/**
 * Delay the low bits of the result with as many clock counts as the high bits addition pipeline take
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::delay_partial_bits(std::ostream& o)
{
	int i;
	ostringstream the_bits;
	
	//make the string that will gather all partial bits from the last level
	the_bits.str("");;
	for (i=partsY;i>0;i--)
	{
	ostringstream temp;
	temp<<"PartialBits_Level_"<< partsY <<"_Reg_"<<i;
		if (i!=1) 
			the_bits << get_delay_signal_name(temp.str(),intadd2->pipeline_depth())<<"("<< multiplier_width_Y-1 << " downto 0" <<")"<<" & ";
		else
			the_bits << get_delay_signal_name(temp.str(),intadd2->pipeline_depth())<<"("<< multiplier_width_Y-1 << " downto 0" <<")";
	
	}
	//setup the pipeline part that will carry the low part of the final result
	o<<tab<<"PartialBits_Reg <= "<<the_bits.str()<<";"<<endl;
	o<<tab<<"partial_bits <= "<<get_delay_signal_name("PartialBits_Reg",intadd1->pipeline_depth())<<";"<<endl;
}



TestIOMap IntMultiplier::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("X"));
	tim.add(*get_signal_by_name("Y"));
	tim.add(*get_signal_by_name("R"));
	return tim;
}

void IntMultiplier::fillTestCase(mpz_class a[])
{
	mpz_class &x = a[0];
	mpz_class &y = a[1];
	mpz_class &r = a[2];

	r = x * y;
}

