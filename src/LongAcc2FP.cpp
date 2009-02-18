/*
 * A post-normalization unit for the FloPoCo Long Accumulator 
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
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "utils.hpp"
#include "Operator.hpp"
#include "LongAcc2FP.hpp"
#include "IntAdder.hpp"

using namespace std;

extern vector<Operator*> oplist;

LongAcc2FP::LongAcc2FP(Target* target, int LSBA, int MSBA, int wEOut, int wFOut): 
	Operator(target), 
	LSBA_(LSBA), MSBA_(MSBA), wEOut_(wEOut), wFOut_(wFOut)
{
	ostringstream name;
	name <<"LongAcc2FP_"
		  <<(LSBA_>=0?"":"M")<<abs(LSBA_)<<"_"
		  <<(MSBA_>=0?"":"M")<<abs(MSBA_)<<"_"
		  <<wEOut_<<"_"<<wFOut_;
	setName(name.str());
	if (target->isPipelined()) 
		setSequential();
	else
		setCombinatorial();
	
	if (isSequential())
		extraPipeLevel=1;
	else
		extraPipeLevel=0;

	// Set up various architectural parameters
	sizeAcc_ = MSBA_-LSBA_+1;	
	
	//instantiate a leading zero/one counter
#ifdef _WIN32
	countWidth_ = intlog2(sizeAcc_)+1;
#else
	countWidth_ = log2(sizeAcc_)+1;
#endif

	lzocShifterSticky_ = new LZOCShifterSticky(target, sizeAcc_, wFOut_ + 1,0,-1); // target, inputNbOfBits,outputNbOfBits, computeSticky, countType 
	oplist.push_back(lzocShifterSticky_);
	
	adder_ = new IntAdder(target, wFOut_ + 1);
	oplist.push_back(adder_);

	// This operator is a sequential one
	setPipelineDepth(lzocShifterSticky_->getPipelineDepth() + adder_->getPipelineDepth()+extraPipeLevel);

	//compute the bias value
	expBias_ = (1<<(wEOut_-1)) -1; 

	countWidth_ = lzocShifterSticky_->getCountWidth();

	//inputs and outputs
	addInput  ("A", sizeAcc_);
	addInput  ("AccOverflow",1);
	addOutput ("R", 3 + wEOut_ + wFOut_);

	addDelaySignal("resultSign0",1,lzocShifterSticky_->getPipelineDepth()+adder_->getPipelineDepth()+extraPipeLevel);	
	addDelaySignal("AccOverflowFlag",1,lzocShifterSticky_->getPipelineDepth());	
	
	addSignal("nZO"    ,countWidth_);
	addDelaySignal("resFrac",wFOut_ + 1,adder_->getPipelineDepth()+extraPipeLevel);
	addSignal("notResFrac", wFOut_ + 1);
	addDelaySignal("postResFrac", wFOut_ + 1,extraPipeLevel);
	addSignal("resultFraction",wFOut_ + 1);
	
	addSignal("expBias"     ,wEOut_);
	addSignal("c2MnZO"      ,countWidth_+1);
	addSignal("expAdj"      ,countWidth_+1);
	addDelaySignal("expRAdjusted",wEOut_,adder_->getPipelineDepth()+extraPipeLevel);
	addSignal("excBits"     ,2);

	//rare case
	if (countWidth_+1 > wEOut_)	{
		addSignal("expAdjustedExt",countWidth_+1);
		addSignal("signExpAdjustedExt",1);	
		addSignal("modulusExpAdjusted",countWidth_);
		addSignal("maxExponent",countWidth_);
		addSignal("expOverflow",1);
		addSignal("expUnderflow",1);
	}
	
	addDelaySignal("excRes",2,adder_->getPipelineDepth()+extraPipeLevel);

}

LongAcc2FP::~LongAcc2FP() {
}


void LongAcc2FP::outputVHDL(ostream& o, string name) {

	licence(o,"Bogdan Pasca (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	lzocShifterSticky_->outputVHDLComponent(o);
	adder_->outputVHDLComponent(o);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	if (isSequential())	
	outputVHDLRegisters(o); 

	o<<endl;

	//the sign of the result 
	o<<tab<<"resultSign0 <= A("<<sizeAcc_-1<<");"<<endl;
	o<<tab<<"AccOverflowFlag <= AccOverflow;"<<endl;	
	
	//count the number of zeros/ones in order to determine the value of the exponent
	//Leading Zero/One counter 
	o<<tab<< "LZOCShifterSticky: " << lzocShifterSticky_->getName() << endl;
	o<<tab<< "      port map ( I => A, "                      << endl; 
	o<<tab<< "                 Count => nZO, "                    << endl; 
	o<<tab<< "                 O => resFrac, "                    << endl; 
	o<<tab<< "                 OZB => resultSign0 ";
	if (isSequential()){ 
	o<<","<<endl;
	o<<tab<< "                 clk => clk, "                  << endl;		
	o<<tab<< "                 rst => rst "                   << endl;
	}
	o<<tab<< "               );"                       << endl<< endl;		

	//the exponent bias	
	o<<tab<<"expBias <= CONV_STD_LOGIC_VECTOR("<<expBias_<<","<<wEOut_<<");"<<endl; //fixed value

	// c2MnZO is  -nZO in 2'c complement with a twist. Usually for 2's complement all bits must be inverted
	// and then a 1 must be added to the LSB. We postpone the extra 1 addition until later.
	o<<tab<<"c2MnZO <= \"1\" & not(nZO);"<<endl; //the bit inversion is done in O(1)
	
	// the 1 is added here, as a carry in bit for the substraction MSBA - nZO,
	// which is an addition in 2's complement
	o<<tab<<"expAdj <= CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<countWidth_+1<<") + c2MnZO + 1;"<<endl;
	
	if (countWidth_+1 < wEOut_){
		// this case is encountered most often.
		// by adding expBias to expAdj, we always get a positive quantity, so no need to converte back to SM
		o<<tab<<"expRAdjusted <= expBias + "<< // + sign extended expAdj
		                         "(("<<wEOut_-1<<" downto "<<countWidth_+1<<" => expAdj("<<countWidth_<<")) & expAdj);"<<endl;
		o<<tab<<"excBits <=\"01\";"<<endl;
	}
	else
		if (countWidth_+1 == wEOut_){
			o<<tab<<"expRAdjusted <= expBias + expAdj;"<<endl;
			o<<tab<<"excBits <=\"01\";"<<endl;
		}	
		else
		{
		// case encountered when we wish to put the contents of the accumulator in a floating point number 
		// with a very small exponent
			o<<tab<<"expAdjustedExt <= (("<<countWidth_<<" downto "<<wEOut_<<" =>'0') & expR) + expAdj;"<<endl;
			
			o<<tab<<"signExpAdjustedExt <= expAdjustedExt("<<countWidth_<<");"<<endl;
			o<<tab<<"modulusExpAdjusted <= (("<<countWidth_-1<<" downto 0 => signExpAdjustedExt) xor expAdjustedExt("<<countWidth_-1<<" downto 0)) + signExpAdjustedExt;"<<endl;
				
			o<<tab<<"maxExponent <= CONV_STD_LOGIC_VECTOR("<<(1<<wEOut_)-1<<" , "<<countWidth_<<");"<<endl;
			
			o<<tab<<"expOverflow <= '1' when (modulusExpAdjusted>maxExponent) and (signExpAdjustedExt='0') else '0';"<<endl;
			o<<tab<<"expUnderflow <= '1' when (modulusExpAdjusted>maxExponent) and (signExpAdjustedExt='1') else '0';"<<endl;
			
			o<<tab<<"excBits(1) <= expOverflow and  not(expUnderflow);"<<endl;
			o<<tab<<"excBits(0) <= not(expOverflow) and not(expUnderflow);"<<endl;
					
			o<<tab<<"expRAdjusted<=expAdjustedExt("<<wEOut_-1<<" downto 0);"<<endl;
		}			
	
	o<<tab<<"excRes <=  excBits when ("<<delaySignal("AccOverflowFlag",lzocShifterSticky_->getPipelineDepth())<<"='0') else \"10\";"<<endl;

	//c2 of the fraction
	o<<tab<< " notResFrac <= not("<<delaySignal("resFrac",extraPipeLevel)<<");"<<endl;

	//count the number of zeros/ones in order to determine the value of the exponent
	//Leading Zero/One counter 
	o<<tab<< "Adder: " << adder_->getName() << endl;
	o<<tab<< "      port map ( X => notResFrac, "  << endl; 
	o<<tab<< "                 Y => "<< zeroGenerator(wFOut_ + 1, 0)<<", " << endl; 
	o<<tab<< "                 Cin =>"<<delaySignal("resultSign0",lzocShifterSticky_->getPipelineDepth()+extraPipeLevel)<< ", "<< endl; 
	o<<tab<< "                 R => postResFrac "           <<endl;
	if (isSequential()){ 
	o<<","<<endl;
	o<<tab<< "                 clk => clk, "                  << endl;		
	o<<tab<< "                 rst => rst "                   << endl;
	}
	o<<tab<< "               );"                       << endl<< endl;		
		
	o<<tab<<" resultFraction <= "<<delaySignal("resFrac",adder_->getPipelineDepth()+extraPipeLevel)<<
			                   " when ("<<delaySignal("resultSign0",lzocShifterSticky_->getPipelineDepth()+adder_->getPipelineDepth()+extraPipeLevel)<<"='0')"<<
							   " else postResFrac;"<<endl;
	
	o<<tab<< "R <= "<<delaySignal("excRes",adder_->getPipelineDepth()+extraPipeLevel) <<" & "
			 <<delaySignal("resultSign0",lzocShifterSticky_->getPipelineDepth()+ adder_->getPipelineDepth()+extraPipeLevel)<<" & "
			 << delaySignal("expRAdjusted",adder_->getPipelineDepth()+extraPipeLevel)<<"("<<wEOut_-1<<" downto 0) & "
			 << "resultFraction("<<wFOut_-1<<" downto 0);"<<endl;
	
	endArchitecture(o);
}



	//determine if the value to be added to the exponent bias is negative or positive.
	// expTest _ 1 if true exponent value (unbiased) is positive
	//         | 0 if exp value negatives
//	o<<tab<<"expTest <= '1' when (CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<countWidth_<<")>nZO) else '0';"<<endl;
	
	// shift 
	/*
	o<<tab<< "left_shifter_component: " << leftShifter_->getName() << endl;
	o<<tab<< "      port map ( X => "<<delaySignal("pipeA",leadZOCounter_->getPipelineDepth())<<", "<< endl; 
	o<<tab<< "                 S => nZO, " << endl; 
	o<<tab<< "                 R => resFrac, " <<endl; 
	o<<tab<< "                 clk => clk, " << endl;
	o<<tab<< "                 rst => rst " << endl;
	o<<tab<< "               );" << endl<<endl;		
	*/
