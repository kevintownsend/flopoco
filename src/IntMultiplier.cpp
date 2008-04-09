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

	//============================ Name Setup procedure=========================
	//  The name has the format: IntMultiplier_wInX_wInY
	//  wInX = width of the X input
	//  wInY = width of the Y input
	name.str("");;
	name <<"IntMultiplier_"<<wInX<<"_"<<wInY;
	unique_name = name.str(); //attribute assignment
	
	//============================ Set up the IO signals========================
	//  X and Y have wInX and wInY bits respectively 
	//  R has wOut bits where wOut = (wInX + WInY) bits
	add_input ("X", wInX);
	add_input ("Y", wInY);
	add_output("R", wOut);
	//==========================================================================
  
	
	
	if(is_sequential()) 
	{
	// set up the variables and parameters for the sequential version	
		
		
		//given the input widths and the specific target this method suggests the chunk widths for X and Y
		bool test = target->suggest_submult_size(multiplier_width_X, multiplier_width_Y, wInX, wInY);
		cout<<endl<<tab<<" -- Frequency can be reached = "<<test<<" suggested sizes: X="<<multiplier_width_X<<" Y="<<multiplier_width_Y<<" --"<<endl;
		
		
		

		// decide in how many parts should the numbers be split depending on the recommended sizes
		partsX = int(ceil( double(wInX) / double(multiplier_width_X)));
		partsY = int(ceil( double(wInY) / double(multiplier_width_Y)));

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
			add_registered_signal_with_sync_reset(name.str(), multiplier_width_X); 
		}

		//declare the signals needed to split Y 
		for (i=1;i<=partsY;i++){
			name.str("");
			name <<"Y_"<<i;
			add_registered_signal_with_sync_reset(name.str(), multiplier_width_Y); 
		}

		//declare the registers needed to store the partial products
		for (i=1;i<=partsY;i++)
			for (j=1;j<=partsX;j++){
				name.str("");;
				name <<"Y_"<<i<<"_X_"<<j;
				add_registered_signal_with_sync_reset(name.str(), multiplier_width_X + multiplier_width_Y);
			} 

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
		
		
		if (partsY>1) {
		//declare the signals which form the two trees (high and low)
			for (i=1;i<=partsY;i++){
				 name.str("");
				 name<<"Buffer_H_"<<i;
				 add_registered_signal_with_sync_reset(name.str(), partsX * multiplier_width_X);
			}
			
			//for the low and high tree
			for (j=0; j<partsY;j++)
				for (i=1;i<=partsY-j;i++){
					nameH.str("");; nameL.str("");;
					nameL <<"L_Level_"<<j<<"_Reg_"<<i;
					nameH <<"H_Level_"<<j<<"_Reg_"<<i;
					add_registered_signal_with_sync_reset(nameL.str(), partsX * multiplier_width_X);
					add_registered_signal_with_sync_reset(nameH.str(), partsX * multiplier_width_X);
				}

			
			//one more registers for the Low part for synchronization 
			add_registered_signal_with_sync_reset("Low1", partsX * multiplier_width_X);
			
			//one more signal for the high part
			add_signal("High1", partsX * multiplier_width_X);

			
			//partsY! registers to hold and compute the low part of the result in the 
			//pipeline levels
			for (j=1; j<=partsY;j++)
				for (i=1;i<=j;i++){
					name.str("");;
					name<<"PartialBits_Level_"<<j<<"_Reg_"<<i;
					add_registered_signal_with_sync_reset(name.str(), multiplier_width_Y + 1 );
			 }

			bool add_test = target->suggest_subadd_size(addition_chunk_width,partsX * multiplier_width_X);
			if (add_test==true)
				cout<<endl<<tab<<"addition chunk size = "<<addition_chunk_width<<endl;
			else
				cerr<<endl<<"Cannot reach the desired frequency for the addition";
			
			pipe_levels = int(ceil(double(partsX * multiplier_width_X)/double(addition_chunk_width))); 
			cout<<endl<<"added pipeline levels = "<<pipe_levels<<endl;
			
			for (j=1; j<=pipe_levels;j++)
				for (i=1;i<=pipe_levels;i++){
					name.str("");;
					name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
					if (i==pipe_levels)
						add_registered_signal_with_sync_reset(name.str(), partsX * multiplier_width_X - (pipe_levels-1)*(addition_chunk_width));
					else
						add_registered_signal_with_sync_reset(name.str(), addition_chunk_width + 1 );
			}

			//TODO
			for (j=1; j<=partsX;j++){ 
				name.str("");;
				name<<"PartialBits_Reg_"<<j;
				add_registered_signal_with_sync_reset(name.str(), partsY * multiplier_width_Y);
			} 
			
			add_signal("temp_result",  partsX * multiplier_width_X );
			add_signal("partial_bits", partsY * multiplier_width_Y );
			add_signal("full_result",  partsX * multiplier_width_X + partsY * multiplier_width_Y);

		}//end if partsY>1
		else
		{
			
			if (partsX>1)
			{
			//Instantiate an IntMultiplier -> multiply the significands

				intadd = new IntAdder(target, partsX * multiplier_width_X);
				oplist.push_back(intadd);

				//Setup the attibute cocerning the IntMultiplier pipeline 
				IntAddPipelineDepth = intadd->pipeline_depth();

				cout<<endl<<" IntAddPipelineDepth = "<<IntAddPipelineDepth<<endl;

				add_signal("addition_term", partsX  * multiplier_width_X);

				add_signal("addition_result", partsX  * multiplier_width_X);
				add_delay_signal("delayed_bits",multiplier_width_Y, IntAddPipelineDepth);
			}
			
						
				add_signal("temp_result", partsX  * multiplier_width_X + partsY * multiplier_width_Y );

 		}

	
		// Set up the pipeline
		//======================================================================
		if (partsY==1)
			if (partsX==1)
				depth=2;
			else
				depth=2 + IntAddPipelineDepth;
		else
			depth= 1 + 2 + partsY + pipe_levels;
 
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
	//--------------------------------------------
	
	Licence(o,"Bogdan Pasca and Florent de Dinechin (2008)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	new_architecture(o,name);
	

	if ((partsY==1)&&(partsX>1))
		intadd->output_vhdl_component(o);//output the IntAdder component
		
	output_vhdl_signal_declarations(o);
	begin_architecture(o);
 
	if (is_sequential()) { 
	
		output_vhdl_registers(o); 
		new_line(o);
		
		//pad inputs with zeros
		pad_inputs(o);
		
		//split X into chunks and assign the chunks to X_1 ...X_partsX 
		split(o,"X",partsX, multiplier_width_X);			  
		split(o,"Y",partsY, multiplier_width_Y);
		
		new_line(o);
		//make all multiplications concurrently
		do_multiplications(o);
			
		new_line(o);
		//decompose the results into high and low part
		decompose_low_high(o);

		new_line(o);		
		//Regroup all the high parts together and all the low parts together
		regroup_low_high(o);
		
		
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
			
			new_line(o);
			//connect the partsY! registers which compute the low part of the result
			connect_partial_bits(o);
			
			new_line(o);
			//pipeline the addition which gives the gigh part of the multiplication result;
	  		pipeline_addition(o); //the result of the addition is called temp_result
	  
			delay_partial_bits(o);
			
			o<<tab<<"full_result <= temp_result & partial_bits;"<<endl; 
			o<<tab<<"R <= full_result("<<  partsX*multiplier_width_X + partsY*multiplier_width_Y - 1 -number_of_zerosX - number_of_zerosY<<" downto 0);"<<endl;
		  
			o<<endl<<endl;
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
	if (reverse == false){
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
		o<<tab<<"L_Level_0_Reg_"<<i<<" <= Low_"<<i<<";"<<endl;
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
		o<<tab<<"H_Level_0_Reg_"<<i<<" <= Buffer_H_"<<i<<"_d;"<<endl;
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
				first_summand<<"("<< zero_generator(multiplier_width_Y,0)<<" & L_Level_"<<i<<"_Reg_"<<j<<"_d("<<partsX*multiplier_width_X-1<<" downto "<<multiplier_width_Y<<"))"; 
				second_summand<<"L_Level_"<<i<<"_Reg_"<<j+1<<"_d";
				o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j<<"<="<< first_summand.str() << " + " << second_summand.str()<<";"<< endl;
			}
			else
				if (j==2) {} //do nothing (the second register on the level is always processed with the first
					else //for the rest of the registers just propagate their contents to the next level
						o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<<"L_Level_"<<i<<"_Reg_"<<j<<"_d"<<";"<<endl;
		}
	}

	//the zeros + [the high bits of the last register of the low part] 
	o<<tab<<"Low1 <= "<<zero_generator(multiplier_width_Y,0)<<" & L_Level_"<<partsY-1<<"_Reg_1_d("<<partsX * multiplier_width_X-1<<" downto "<<multiplier_width_Y<<")"<< ";"<<endl;
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
				//concatenate zeros in front of the first register, add it with thesecond register and place the result in the next level's first register
				first_summand.str("");; second_summand.str("");;
				first_summand<<"("<< zero_generator(multiplier_width_Y,0)<<" & H_Level_"<<i<<"_Reg_"<<j<<"_d("<<partsX*multiplier_width_X-1<<" downto "<<multiplier_width_Y<<"))"; 
				second_summand<<"H_Level_"<<i<<"_Reg_"<<j+1<<"_d";
				o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j<<"<="<< first_summand.str() << " + " << second_summand.str()<<";"<< endl;
			}
			else
				if (j==2) {} //do nothing (the second register on the level is always processed with the first
					else //for the rest of the registers just propagate their contents to the next level
						o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<<"H_Level_"<<i<<"_Reg_"<<j<<"_d;"<<endl;
		}
	}

	o<<tab<<"High1 <= H_Level_"<<partsY-1<<"_Reg_1_d"<<";"<<endl;
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
				o<<tab<<"PartialBits_Level_1_Reg_1 <= \"0\" & L_Level_0_Reg_1_d("<< multiplier_width_Y -1 <<" downto 0);"<< endl;
			else{
				//just propagate the registers from 1 to j-1
				if (i<j)
					o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i <<"<= "<<"PartialBits_Level_"<<j-1<<"_Reg_"<< i <<"_d;"<<endl;
				else //if i=j compute the jth register
					o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i <<"<= "<<"(\"0\" & L_Level_"<<j-1<<"_Reg_1_d("<<multiplier_width_Y-1<<" downto 0)) + "<<"(\"0\" & H_Level_"<<j-2<<"_Reg_1_d("<<multiplier_width_Y-1<<" downto 0)) + "<< "PartialBits_Level_"<<j-1<<"_Reg_"<<i-1<<"_d("<< multiplier_width_Y  <<") ;"<<endl; 
			}	
		}
}


/**
 * The addition which gives the high part of the product between X and Y is pipelined to obtain better performance
 * @param[in,out] 	o 			- the stream
 **/
void IntMultiplier::pipeline_addition(std::ostream& o)
{
int i,j;
	//connect first parts of addition registers
	for (j=1; j<=pipe_levels;j++){
		if ((j==1)&&(pipe_levels!=1))
			o<<tab<<"Last_Addition_Level_1_Reg_"<<j<<" <= (\"0\" & Low1_d("<<j*addition_chunk_width -1<<" downto "<<(j-1)*addition_chunk_width<<")) + ( \"0\" & High1("<<j*addition_chunk_width -1<<" downto "<<(j-1)*addition_chunk_width<<")) + PartialBits_Level_"<< partsY <<"_Reg_"<< partsY<<"_d("<< multiplier_width_Y <<");"<<endl;  
		else
			if ((j==1)&&(pipe_levels==1))
				o<<tab<<"Last_Addition_Level_1_Reg_"<<j<<" <= Low1_d("<<partsX*multiplier_width_X-1<<" downto "<<"0) + High1("<<partsX*multiplier_width_X-1<<" downto 0) + PartialBits_Level_"<< partsY <<"_Reg_"<< partsY<<"_d("<< multiplier_width_Y <<");"<<endl;  
			else
				if (j==pipe_levels)
					o<<tab<<"Last_Addition_Level_1_Reg_"<<j<<" <= Low1_d("<<partsX*multiplier_width_X-1<<" downto "<<(j-1)*addition_chunk_width<<") + High1("<<partsX*multiplier_width_X-1<<" downto "<<(j-1)*addition_chunk_width<<");"<<endl;  
				else
					o<<tab<<"Last_Addition_Level_1_Reg_"<<j<<" <= (\"0\" & Low1_d("<<j*addition_chunk_width -1<<" downto "<<(j-1)*addition_chunk_width<<")) + (\"0\" & High1("<<j*addition_chunk_width -1<<" downto "<<(j-1)*addition_chunk_width<<"));"<<endl;  
	}

	//put tohether the pipeline
	for (j=1; j<=partsX-1;j++)
		for (i=1;i<=partsX;i++)
			if ((i<=j)||(i>j+1))
				o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
			else
				o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<addition_chunk_width<<") ;"<<endl; 
				 
	//put the result back together
	ostringstream addition_string;
	for (i=pipe_levels;i>0;i--) 
		if (i!=1)
			if (i!=pipe_levels)
				addition_string<<"Last_Addition_Level_"<<pipe_levels<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0) &" ;
			else
				addition_string<<"Last_Addition_Level_"<<pipe_levels<<"_Reg_"<<i<<"_d &";
		else
			if (pipe_levels!=1)
				addition_string<<"Last_Addition_Level_"<<pipe_levels<<"_Reg_"<<i<<"_d("<<addition_chunk_width -1<<" downto 0);" ;
			else
				addition_string<<"Last_Addition_Level_"<<pipe_levels<<"_Reg_"<<i<<"_d("<<partsX*multiplier_width_X -1<<" downto 0);" ;

	o<<tab<<"temp_result <= "<<addition_string.str()<<endl;
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
		if (i!=1) 
			the_bits << "PartialBits_Level_"<< partsY <<"_Reg_"<<i<<"_d("<< multiplier_width_Y-1 << " downto 0" <<")"<<" & ";
		else
			the_bits << "PartialBits_Level_"<< partsY <<"_Reg_"<<i<<"_d("<< multiplier_width_Y-1 << " downto 0" <<")";

	//setup the pipeline part that will carry the low part of the final result
	o<<tab<<"PartialBits_Reg_1 <= "<<the_bits.str()<<";"<<endl;
	for (i=2;i<=pipe_levels;i++)
		o<<tab<<"PartialBits_Reg_"<<i<<" <= PartialBits_Reg_"<<i-1<<"_d;"<<endl;

	o<<tab<<"partial_bits <= PartialBits_Reg_"<<pipe_levels<<"_d;"<<endl;
}









//////////////////////////TEST FUNCTIONS//////////////////////////


void IntMultiplier::add_standard_test_cases(vector<TestCase> &list){
  map<string, mpz_class> in;
  multimap<string, mpz_class> out;

  // 0*0 = 0
  in.clear();
  out.clear();
  in["X"]=0;    in["Y"]=0;    out.insert(pair<string,mpz_class>("R",in["X"] * in["Y"]));  
  add_test_case(list, in, out);

  // 42*17 = ??
  out.clear();
  in["X"]=42;    in["Y"]=17;         out.insert(pair<string,mpz_class>("R",in["X"] * in["Y"]));
  add_test_case(list, in, out);
}

void IntMultiplier::add_random_test_cases(vector<TestCase> &list, int n){
  //TODO
}
