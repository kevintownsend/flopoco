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

IntMultiplier:: IntMultiplier(Target* target, int wInX, int wInY) :
	Operator(target), wInX_(wInX), wInY_(wInY), wOut_(wInX + wInY){
 
	int depth, i, j;
	ostringstream name, nameH, nameL;

	setOperatorType();
	setOperatorName();
	
	addInput ("X", wInX_);
	addInput ("Y", wInY_);
	addOutput("R", wOut_); /* wOut_ = wInX_ + wInY_ */
  
	if(isSequential()) 
	{
	// set up the variables and parameters for the sequential version	
		//given the input widths and the specific target this method suggests the chunk widths for X and Y
		bool test = target->suggestSubmultSize(multiplierWidthX_, multiplierWidthY_, wInX_, wInY_);
		cout<<"Multiplication frequency report:"<<endl;
		if (test)
			cout<<tab<<"Frequency can be reached. Suggested multiplier sizes are: X="<<multiplierWidthX_<<" Y="<<multiplierWidthY_<<endl;
		else
			cout<<tab<<"WARNING: It might not be possible to reach frequency "<<target->frequencyMHz() << "MHz"
			         <<" with suggested sizes: X="<<multiplierWidthX_<<" Y="<<multiplierWidthY_<<endl;		

		// decide in how many parts should the numbers be split depending on the recommended sizes
		partsX_ = ((wInX_ % multiplierWidthX_)==0) ? (wInX_ / multiplierWidthX_) : (wInX_ / multiplierWidthX_+1);
		partsY_ = ((wInY_ % multiplierWidthY_)==0) ? (wInY_ / multiplierWidthY_) : (wInY_ / multiplierWidthY_+1);
		
		//reverse_
		if (multiplierWidthX_ < multiplierWidthY_){
			i = wInY_;
			wInY_ = wInX_;
			wInX_ = i;
			
			i = multiplierWidthX_;
			multiplierWidthX_ = multiplierWidthY_;
			multiplierWidthY_ = i;
			
			i = partsX_;
			partsX_ = partsY_;
			partsY_ = i;
			
			reverse_ = true;
		}else 
			reverse_ = false;
		
		numberOfZerosX_ = partsX_ * multiplierWidthX_ - wInX_;
		numberOfZerosY_ = partsY_ * multiplierWidthY_ - wInY_;
		
		// signals that will contain the inputs X and Y but padded with 0 so that 
		// their width becomes a multiple of the target multiplier width
		addSignal("i_X", partsX_ * multiplierWidthX_);
		addSignal("i_Y", partsY_ * multiplierWidthY_);

		//declare the signals needed to split X 
		for (i=1;i<=partsX_;i++){
			name.str("");
			name <<"X_"<<i;
			//addDelaySignal(name.str(), multiplierWidthX_);
			addDelaySignal(name.str(), multiplierWidthX_);  
		}

		//declare the signals needed to split Y 
		for (i=1;i<=partsY_;i++){
			name.str("");
			name <<"Y_"<<i;
			//addDelaySignal(name.str(), multiplierWidthY_); 
			addDelaySignal(name.str(), multiplierWidthY_); 
		}

		//declare the registers needed to store the partial products
		for (i=1;i<=partsY_;i++)
			for (j=1;j<=partsX_;j++){
				name.str("");;
				name <<"Y_"<<i<<"_X_"<<j;
				//addDelaySignal(name.str(), multiplierWidthX_ + multiplierWidthY_);
				addDelaySignal(name.str(), multiplierWidthX_ + multiplierWidthY_);
			} 

		if (!((partsX_==partsY_) && (partsX_==1)))
		{
		//this signals are not needed when simple 1blockx1block multiplications are performed
		
			//declare the high/low part decomposition signals
			for (i=1;i<=partsY_;i++)
				for (j=1;j<=partsX_;j++){
					nameH.str("");; nameL.str("");;
					nameH <<"Y_"<<i<<"_X_"<<j<<"_H";
					nameL <<"Y_"<<i<<"_X_"<<j<<"_L";
					addSignal(nameH.str(), multiplierWidthX_); //the largest of the two
					addSignal(nameL.str(), multiplierWidthX_);
				} 

			//signals to regroup all the high parts and all the low parts
			for (i=1;i<=partsY_;i++)  {
				nameH.str("");; nameL.str("");;
				nameH <<"High_"<<i;
				nameL <<"Low_"<<i;
				addSignal(nameH.str(), partsX_ * multiplierWidthX_);
				addSignal(nameL.str(), partsX_ * multiplierWidthX_);
			}
		}  
		
		if (partsY_>1) {
		//declare the signals which form the two trees (high and low)
			
			/* Use a state of the art adder to take care of pipelining the long additions */
			intAdd2_ = new IntAdder(target, 1 + partsX_ * multiplierWidthX_);
			oplist.push_back(intAdd2_);
			
			for (i=1;i<=partsY_;i++){
				 name.str("");
				 name<<"Buffer_H_"<<i;
				 //add_delay_signal(name.str(), partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
				 addDelaySignalBus(name.str(), partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
			}
			
			//for the low and high tree
			for (j=0; j<partsY_;j++)
				for (i=1;i<=partsY_-j;i++){
					nameH.str("");; nameL.str("");;
					nameL <<"L_Level_"<<j<<"_Reg_"<<i;
					nameH <<"H_Level_"<<j<<"_Reg_"<<i;
					if ((i==1)||(i==2)){
						
						addSignal(nameL.str(), 1 + partsX_ * multiplierWidthX_);
						addSignal(nameH.str(), 1 + partsX_ * multiplierWidthX_);
						ostringstream temp;
						if (i==1)
						{
							temp.str("");
							temp<<"L_Level_"<<j<<"_summand";
							addSignal(temp.str(),1 + partsX_ * multiplierWidthX_);
							temp.str("");
							temp<<"H_Level_"<<j<<"_summand";
							addSignal(temp.str(),1 + partsX_ * multiplierWidthX_);
						}
						
					}
					else{
						//add_delay_signal(nameL.str(), 1 + partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
						addDelaySignalBus(nameL.str(), 1 + partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
						
						//add_delay_signal(nameH.str(), 1 + partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
						addDelaySignalBus(nameH.str(), 1 + partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
					}
				}
		
			//one more registers for the Low part for synchronization 
			//add_delay_signal("Low1", 1 + partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
			addDelaySignalBus("Low1", 1 + partsX_ * multiplierWidthX_,intAdd2_->getPipelineDepth());
			
			//one more signal for the high part
			addSignal("High1", 1 + partsX_ * multiplierWidthX_);

			
			//partsY_! registers to hold and compute the low part of the result in the 
			//pipeline levels
			for (j=1; j<=partsY_;j++)
				for (i=1;i<=j;i++){
					name.str("");;
					name<<"PartialBits_Level_"<<j<<"_Reg_"<<i;
					//add_delay_signal(name.str(), multiplierWidthY_ + 1,intAdd2_->getPipelineDepth() );
					addDelaySignalBus(name.str(), multiplierWidthY_ + 1,intAdd2_->getPipelineDepth() );
			 }

			/* Use a state of the art adder to take care of pipelining the last addition */
			intAdd1_ = new IntAdder(target, partsX_ * multiplierWidthX_);
			oplist.push_back(intAdd1_);
			
			//add_delay_signal("PartialBits_Reg",partsY_ * multiplierWidthY_,intAdd1_->getPipelineDepth());
			addDelaySignalBus("PartialBits_Reg",partsY_ * multiplierWidthY_,intAdd1_->getPipelineDepth());  
			
			addSignal("temp_result",  partsX_ * multiplierWidthX_ );
			addSignal("partial_bits", partsY_ * multiplierWidthY_ );
			addSignal("full_result",  partsX_ * multiplierWidthX_ + partsY_ * multiplierWidthY_);

		}//end if partsY_>1
		else
		{
			
			if (partsX_>1)
			{
				intAdd_ = new IntAdder(target, partsX_ * multiplierWidthX_);
				oplist.push_back(intAdd_);

				//Setup the attibute cocerning the IntMultiplier pipeline 
				IntAddPipelineDepth = intAdd_->getPipelineDepth();

				cout<<endl<<" IntAddPipelineDepth = "<<IntAddPipelineDepth<<endl;

				addSignal("addition_term", partsX_  * multiplierWidthX_);

				addSignal("addition_result", partsX_  * multiplierWidthX_);
				
				//add_delay_signal("delayed_bits",multiplierWidthY_, IntAddPipelineDepth);
				addDelaySignalBus("delayed_bits",multiplierWidthY_, IntAddPipelineDepth);
			}
	
			addSignal("temp_result", partsX_  * multiplierWidthX_ + partsY_ * multiplierWidthY_);
 		}
	
		// Set up the pipeline
		if (partsY_==1)
			if (partsX_==1)
				depth=2;
			else
				depth=2 + IntAddPipelineDepth;
		else
			depth= 2 + partsY_*intAdd2_->getPipelineDepth() + intAdd1_->getPipelineDepth();
 
		setPipelineDepth(depth);
  
	}//end if sequential
}

IntMultiplier::~IntMultiplier() {
}

void IntMultiplier::setOperatorName(){	
	/* Name Setup procedure
	 *  The name has the format: IntMultiplier_wInX__wInY_
	 *  wInX_ = width of the X input
	 *  wInY_ = width of the Y input
	 */  
	ostringstream name;
	name <<"IntMultiplier_"<<wInX_<<"_"<<wInY_;
	uniqueName_ = name.str(); 
}

void IntMultiplier::outputVHDL(std::ostream& o, std::string name) {
	
	//variables
	ostringstream zeros, names, nameH, nameL, concatH, concatL,
					  first_summand, second_summand, the_bits;

	int i, j, k;
		
	licence(o,"Bogdan Pasca and Florent de Dinechin (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);

	
	if (isSequential()) { 
		if ((partsY_==1)&&(partsX_>1))
			intAdd_->outputVHDLComponent(o);//output the IntAdder component
		else
			if (!((partsX_==partsY_) && (partsX_==1))){
				intAdd1_->outputVHDLComponent(o);
				intAdd2_->outputVHDLComponent(o);
			}
	}
		
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
 
	if (isSequential()) { 
	
		outputVHDLRegisters(o); 
		
		//pad inputs with zeros
		padInputs(o);
		
		//split X into chunks and assign the chunks to X_1 ...X_partsX_ 
		split(o,"X",partsX_, multiplierWidthX_);			  
		split(o,"Y",partsY_, multiplierWidthY_);
		
		//make all multiplications concurrently
		doMultiplications(o);
		
		if (!((partsX_==partsY_) && (partsX_==1)))
		{	
			//decompose the results into high and low part
			decomposeLowHigh(o);
			//Regroup all the high parts together and all the low parts together
			regroupLowHigh(o);
		}
		
		newLine(o);
		//pass to the adder tree
		if (partsY_==1) {
			if (partsX_>1){
				
				o<<"addition_term <= ("<< zeroGenerator(multiplierWidthY_,0)<<" & Low_1("<< partsX_*multiplierWidthX_-1<<" downto "<< multiplierWidthY_<<"));"<<endl;
				mapAdder(o, "addition_term", "High_1", "addition_result");
				
				//some of the bits of Low_1 do no participate in the addtion. These bits together with the result of the addition need to be merged to form the final result.
				//=> they need to be delayed with a number of clock cycles equal to the pipeline levels of the addition
				o<<tab<<"delayed_bits <= Low_1("<<multiplierWidthY_-1<<" downto 0 )"<<";"<<endl;
				
				//the concatenation of the addtion result and the delayed bits form the result
				o<<tab<<"temp_result <= addition_result & "<<getDelaySignalName("delayed_bits", IntAddPipelineDepth) <<";"<<endl;
			}else{
				//when both X and Y have 1 part the result is simply the product of these two parts
				o<<tab<<"temp_result <= Y_1_X_1_d;"<<endl;
			}
			
			//output the result
			o<<tab<<"R <= temp_result("<<partsX_*multiplierWidthX_ + partsY_*multiplierWidthY_-1 -numberOfZerosX_ - numberOfZerosY_<<" downto 0);"<<endl;
		}
		else //(the number of parts of Y is >1)
		{ 

			newLine(o);
			//connect the low part of the multiplication results to the low tree
			connectLow(o);
						
			newLine(o);
			//the high part of the multiplication needs to be delayed one clock cycle before it will be inserted into the adder tree
			highToBuffer(o);	
			
			newLine(o);
			//the buffers need to be connected to the high tree
			connectBuffers(o);
				
			//Connect the low part tohether
			linkLowAdderStructure(o);
			 
			//Connect the HIGH part tohether 
			linkHighAdderStructure(o);
			
			//connect the partsY_! registers which compute the low part of the result
			connectPartialBits(o);
			
			//pipeline the addition which gives the gigh part of the multiplication result;
	  	pipelineAddition(o); //the result of the addition is called temp_result
	  
			delayPartialBits(o);
						
			o<<tab<<"full_result <= temp_result & partial_bits;"<<endl; 
			o<<tab<<"R <= full_result("<<  partsX_*multiplierWidthX_ + partsY_*multiplierWidthY_ - 1 -numberOfZerosX_ - numberOfZerosY_<<" downto 0);"<<endl<<endl;
		  
		}
	}else{ 
		//the combinational version
		o << tab<<"R <= X * Y ;" << endl;
	}
	
	endArchitecture(o);
}

void IntMultiplier::padInputs(std::ostream& o){
ostringstream padString;

	if (reverse_ == true){
		if (numberOfZerosY_==0)
			padString<<"";
		else 
			padString<<zeroGenerator(numberOfZerosY_,0)<<" &";
		o<<tab<< "i_Y <= "<< padString.str() <<"X;"<<endl;

		padString.str(""); //steam initialization		
		if (numberOfZerosX_==0)
			padString<<"";
		else 
			padString<<zeroGenerator(numberOfZerosX_,0)<<" &";

		o<<tab<< "i_X <= "<< padString.str()<<" Y;"<<endl;
		newLine(o);
	}else{
		if (numberOfZerosX_==0)
			padString<<"";
		else 
			padString<<zeroGenerator(numberOfZerosX_,0)<<" &";
		o<<tab<< "i_X <= "<< padString.str()<<" X;"<<endl;

		padString.str(""); //steam initialization		

		if (numberOfZerosY_==0)
			padString<<"";
		else 
			padString<<zeroGenerator(numberOfZerosY_,0)<<" &";
		o<<tab<< "i_Y <= "<< padString.str() <<" Y;"<<endl;

		newLine(o);
	}		
}

void IntMultiplier::split(std::ostream& o, std::string name, int parts, int part_width){
int i;
ostringstream temp;
	
	for (i=1;i<=parts;i++){
		temp.str("");
		temp << name <<"_"<<i;
		o<<tab<< temp.str() << " <= i_"<<name<<"("<< i * part_width -1 <<" downto "<< (i-1) * part_width<<");"<<endl;
	}
}

void IntMultiplier::doMultiplications(std::ostream& o)
{
int i,j;
ostringstream names;	
	for (i=1;i<=partsY_;i++)
		for (j=1;j<=partsX_;j++){  
			names.str("");;
			names <<"Y_"<<i<<"_X_"<<j;
			o<<tab<< names.str() << " <= " << "Y_"<< i <<"_d * " << "X_"<<j << "_d;"<<endl;
		} 
	newLine(o);
}

void IntMultiplier::decomposeLowHigh(std::ostream& o)
{
int i,j;
ostringstream names, nameH, nameL;
	for (i=1;i<=partsY_;i++)
		for (j=1;j<=partsX_;j++){
			names.str("");; nameH.str("");; nameL.str("");;
			names <<"Y_"<<i<<"_X_"<<j<<"_d";
			nameH <<"Y_"<<i<<"_X_"<<j<<"_H";
			nameL <<"Y_"<<i<<"_X_"<<j<<"_L";
			o<<tab<< nameH.str() << " <= " << names.str() <<"("<< multiplierWidthX_ + multiplierWidthY_ - 1 <<" downto " << multiplierWidthY_ <<");"<<endl;
			o<<tab<< nameL.str() << " <= " <<zeroGenerator(multiplierWidthX_ - multiplierWidthY_ ,0)<<" & "<< names.str() <<"("<< multiplierWidthY_ - 1 <<" downto " << 0 <<");"<<endl;
			} 
}

void IntMultiplier::regroupLowHigh(std::ostream& o)
{
int i,j;
ostringstream nameL,nameH, concatH, concatL;
	for (i=1;i<=partsY_;i++){
		nameH.str(""); nameL.str(""); concatH.str(""); concatL.str("");
		nameH <<"High_"<<i;
		nameL <<"Low_"<<i;
		//create the concatenation strings for low and high parts
		for (j=partsX_;j>=1;j--){
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
	newLine(o);
}

void IntMultiplier::mapAdder(std::ostream& o,std::string left_term, std::string right_term, std::string result )
{
	o<<tab<< "int_adder_component: " << intAdd_->getOperatorName() << endl;
	o<<tab<< "  port map ( X => "<< left_term << ", " << endl; 
	o<<tab<< "             Y => "<< right_term<< ", " << endl; 
	o<<tab<< "             Cin => '0' ," << endl;
	o<<tab<< "             R => "<< result<< ", " << endl; 
	o<<tab<< "             clk => clk, " << endl;
	o<<tab<< "             rst => rst " << endl;
	o<<tab<< "               );" << endl<<endl;
}

void IntMultiplier::connectLow(std::ostream& o)
{
int i;
	for(i=1;i<=partsY_;i++)
		o<<tab<<"L_Level_0_Reg_"<<i<<" <= \"0\" & Low_"<<i<<";"<<endl;
}

void IntMultiplier::highToBuffer(std::ostream& o)
{
int i;
	for(i=1;i<=partsY_;i++)
		o<<tab<<"Buffer_H_"<< i << " <= High_" <<i<<";"<<endl;
}

void IntMultiplier::connectBuffers(std::ostream& o)
{
int i;
	for(i=1;i<=partsY_;i++)
	{
	ostringstream temp;
	temp<<"Buffer_H_"<<i;
		o<<tab<<"H_Level_0_Reg_"<<i<<" <= \"0\" & "<<getDelaySignalName(temp.str(),intAdd2_->getPipelineDepth())<<";"<<endl;
	}
}

void IntMultiplier::linkLowAdderStructure(std::ostream& o)
{
int i,j;
ostringstream first_summand, second_summand;
	for (i=0;i<partsY_-1;i++) {
		for (j=1;j<=partsY_-i;j++) {
			if (j==1) {
				//concatenate zeros in front of the first register, add it with the second register and place the result in the next level's first register
				first_summand.str(""); second_summand.str("");
				first_summand<<"("<< zeroGenerator(multiplierWidthY_ ,0)<<" & L_Level_"<<i<<"_Reg_"<<j<<"("<<partsX_*multiplierWidthX_<<" downto "<<multiplierWidthY_<<"))"; 
				second_summand<<"L_Level_"<<i<<"_Reg_"<<j+1;
				
				o<<tab<<"L_Level_"<<i+1<<"_summand <= "<<first_summand.str()<<";"<<endl;								
				o<<tab<< "int_adder_component_low_"<<i<<": " << intAdd2_->getOperatorName() << endl;
				o<<tab<< "  port map ( X => L_Level_"<<i+1<<"_summand , " << endl; 
				o<<tab<< "             Y => "<<second_summand.str()<< ", " << endl; 
				o<<tab<< "             Cin => '0'," << endl;
				o<<tab<< "             R => L_Level_"<<i+1<<"_Reg_"<<j<<", " << endl; 
				o<<tab<< "             clk => clk, " << endl;
				o<<tab<< "             rst => rst " << endl;
				o<<tab<< "               );" << endl<<endl;
			}
			else
				if (j==2) {} //do nothing (the second register on the level is always processed with the first
					else //for the rest of the registers just propagate their contents to the next level
					{
						ostringstream tmpSig;
						tmpSig<<"L_Level_"<<i<<"_Reg_"<<j;
						o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<< getDelaySignalName(tmpSig.str(),intAdd2_->getPipelineDepth())<<";"<<endl;
					}	
		}
	}

	//the zeros + [the high bits of the last register of the low part] 
	o<<tab<<"Low1 <= "<<zeroGenerator(multiplierWidthY_,0)<<" & L_Level_"<<partsY_-1<<"_Reg_1("<<partsX_ * multiplierWidthX_<<" downto "<<multiplierWidthY_<<")"<< ";"<<endl;
}

void IntMultiplier::linkHighAdderStructure(std::ostream& o)
{
int i,j;
ostringstream first_summand, second_summand;
	for (i=0;i<partsY_-1;i++){
		for (j=1;j<=partsY_-i;j++){
			if (j==1){
				//concatenate zeros in front of the first register, add it with the second register and place the result in the next level's first register
				first_summand.str("");; second_summand.str("");;
				first_summand<<"("<< zeroGenerator(multiplierWidthY_,0)<<" & H_Level_"<<i<<"_Reg_"<<j<<"("<<partsX_*multiplierWidthX_<<" downto "<<multiplierWidthY_<<"))"; 
				second_summand<<"H_Level_"<<i<<"_Reg_"<<j+1<<"";
				
				o<<tab<<"H_Level_"<<i+1<<"_summand <= "<<first_summand.str()<<";"<<endl;								
				o<<tab<< "int_adder_component_high_"<<i<<": " << intAdd2_->getOperatorName() << endl;
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
						o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<< getDelaySignalName(tmpSig.str(),intAdd2_->getPipelineDepth())<<";"<<endl;
								
					//	o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<<"H_Level_"<<i<<"_Reg_"<<j<<"_d;"<<endl;
					}
		}
	}

	o<<tab<<"High1 <= H_Level_"<<partsY_-1<<"_Reg_1"<<";"<<endl;
}

void IntMultiplier::connectPartialBits(std::ostream& o)
{
int i,j;
	for (j=1; j<=partsY_;j++)
		for (i=1;i<=j;i++){
			if ((j==1) and (i==1) )
				o<<tab<<"PartialBits_Level_1_Reg_1 <= \"0\" & L_Level_0_Reg_1("<< multiplierWidthY_ -1 <<" downto 0);"<< endl;
			else{
				//just propagate the registers from 1 to j-1
				if (i<j)
				{
					ostringstream temp;
					temp<<"PartialBits_Level_"<<j-1<<"_Reg_"<<i;
					o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i <<"<= "<<getDelaySignalName(temp.str(),intAdd2_->getPipelineDepth())<<";"<<endl;
				}
				else //if i=j compute the jth register
				{
					ostringstream temp;
					temp<<"PartialBits_Level_"<<j-1<<"_Reg_"<<i-1;
					o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i<<"<= "<<"(\"0\" & L_Level_"<<j-1<<"_Reg_1("<<multiplierWidthY_-1<<" downto 0)) + "
					      <<"(\"0\" & H_Level_"<<j-2<<"_Reg_1("<<multiplierWidthY_-1<<" downto 0)) + "
					      << getDelaySignalName(temp.str(),intAdd2_->getPipelineDepth())<<"("<<multiplierWidthY_<<");"<<endl; 
				}
			}	
		}
}

void IntMultiplier::pipelineAddition(std::ostream& o)
{
ostringstream temp;
temp<<"PartialBits_Level_"<< partsY_ <<"_Reg_"<< partsY_;
		o<<tab<< "int_adder_component: " << intAdd1_->getOperatorName() << endl;
		o<<tab<< "  port map ( X => "<< getDelaySignalName("Low1",intAdd2_->getPipelineDepth())<< "("<<partsX_ * multiplierWidthX_ -1<<" downto 0"<<") , " << endl; 
		o<<tab<< "             Y => High1("<<partsX_ * multiplierWidthX_ -1<<" downto 0"<<"), " << endl; 
		o<<tab<< "             Cin => "<<getDelaySignalName(temp.str(),intAdd2_->getPipelineDepth()) <<"("<< multiplierWidthY_ <<")," << endl;
		o<<tab<< "             R => temp_result, " << endl; 
		o<<tab<< "             clk => clk, " << endl;
		o<<tab<< "             rst => rst " << endl;
		o<<tab<< "               );" << endl<<endl;
}

void IntMultiplier::delayPartialBits(std::ostream& o)
{
	int i;
	ostringstream the_bits;
	
	//make the string that will gather all partial bits from the last level
	the_bits.str("");;
	for (i=partsY_;i>0;i--)
	{
	ostringstream temp;
	temp<<"PartialBits_Level_"<< partsY_ <<"_Reg_"<<i;
		if (i!=1) 
			the_bits << getDelaySignalName(temp.str(),intAdd2_->getPipelineDepth())<<"("<< multiplierWidthY_-1 << " downto 0" <<")"<<" & ";
		else
			the_bits << getDelaySignalName(temp.str(),intAdd2_->getPipelineDepth())<<"("<< multiplierWidthY_-1 << " downto 0" <<")";
	
	}
	//setup the pipeline part that will carry the low part of the final result
	o<<tab<<"PartialBits_Reg <= "<<the_bits.str()<<";"<<endl;
	o<<tab<<"partial_bits <= "<<getDelaySignalName("PartialBits_Reg",intAdd1_->getPipelineDepth())<<";"<<endl;
}

TestIOMap IntMultiplier::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("X"));
	tim.add(*getSignalByName("Y"));
	tim.add(*getSignalByName("R"));
	return tim;
}

void IntMultiplier::fillTestCase(mpz_class a[])
{
	mpz_class &x = a[0];
	mpz_class &y = a[1];
	mpz_class &r = a[2];

	r = x * y;
}

