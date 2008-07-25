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
#include "utils.hpp"
#include "Operator.hpp"
#include "LongAcc2FP.hpp"

using namespace std;

extern vector<Operator*> oplist;

LongAcc2FP::LongAcc2FP(Target* target, int MaxMSBX, int LSBA, int MSBA, int wEOut, int wFOut): 
	Operator(target), 
	MaxMSBX_(MaxMSBX), LSBA_(LSBA), MSBA_(MSBA), wEOut_(wEOut), wFOut_(wFOut)
{
	int i;
	setOperatorName();
	setSequential();

	// Set up various architectural parameters
	sizeAcc_ = MSBA_-LSBA_+1;	
	
	//instantiate a leading zero/one counter
	wOutLZOC_ = log2(sizeAcc_)+1;
	leadZOCounter_ = new LZOC(target, sizeAcc_, wOutLZOC_);
	oplist.push_back(leadZOCounter_);

	//instantiate a left shifter
	leftShifter_ = new Shifter(target,sizeAcc_,sizeAcc_,Left);
	oplist.push_back(leftShifter_);

	// This operator is a sequential one
	setPipelineDepth(leadZOCounter_->getPipelineDepth() + leftShifter_->getPipelineDepth());

	//compute the bias value
	expBias_ = (1<<(wEOut_-1)) -1; 

	//inputs and outputs
	addInput  ("A", sizeAcc_);
	addOutput ("R", 3 + wEOut_ + wFOut_);

	
	addDelaySignalNoReset("resultSign0",1,leadZOCounter_->getPipelineDepth() + leftShifter_->getPipelineDepth());	
	addDelaySignalBusNoReset("pipeA",sizeAcc_,leadZOCounter_->getPipelineDepth());
	addDelaySignalNoReset("expRAdjusted",wEOut_,leftShifter_->getPipelineDepth());
	addDelaySignalNoReset("excBits",2, leftShifter_->getPipelineDepth());
	
	addSignal("nZO",wOutLZOC_);
	addSignal("c2nZO",wOutLZOC_+1);
	addSignal("expTest",1);
	addSignal("expAdjValue",wEOut_);
	addSignal("excRes",2);
	addSignal("expR",wEOut_);
	
	addSignal("resFrac",sizeAcc_ + sizeAcc_);
	
	addSignal("posExpAdj",wOutLZOC_+1);
	addSignal("negExpAdj",wOutLZOC_+1);
	addSignal("expAdj",wOutLZOC_+1);
		
	addSignal("expAdjustedExt",wOutLZOC_+1);
	addSignal("signExpAdjustedExt",1);
	addSignal("modulusExpAdjusted",wOutLZOC_);
	addSignal("maxExponent",wOutLZOC_);
	addSignal("expOverflow",1);
	addSignal("expUnderflow",1);
}

LongAcc2FP::~LongAcc2FP() {
}

void LongAcc2FP::setOperatorName(){
	ostringstream name;
	name <<"LongAcc2FP_"
			 <<(MaxMSBX_>=0?"":"M")<<abs(MaxMSBX_)<<"_"
			 <<(LSBA_>=0?"":"M")<<abs(LSBA_)<<"_"
			 <<(MSBA_>=0?"":"M")<<abs(MSBA_)<<"_"
			 <<wEOut_<<"_"<<wFOut_;
	uniqueName_=name.str();

}

void LongAcc2FP::outputVHDL(ostream& o, string name) {
int i;
	licence(o,"Bogdan Pasca (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	leadZOCounter_->outputVHDLComponent(o);
	leftShifter_->outputVHDLComponent(o);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	outputVHDLRegisters(o); o<<endl;

	//the sign of the result 
	o<<tab<<"resultSign0 <= A("<<sizeAcc_-1<<");"<<endl;
	
	o<<tab<<"pipeA <= A;"<<endl;
	
	//count the number of zeros/ones in order to determine the size of the exponent
	//Leading Zero/One counter 
	o<<tab<< "LZOC_component: " << leadZOCounter_->getOperatorName() << endl;
	o<<tab<< "      port map ( I => A, "                      << endl; 
	o<<tab<< "                 OZB => resultSign0, "          << endl; 
	o<<tab<< "                 O => nZO, "                    << endl; 
	o<<tab<< "                 clk => clk, "                  << endl;
	o<<tab<< "                 rst => rst "                   << endl;
	o<<tab<< "               );"                       << endl<< endl;		

	//readjust exponent
	
	//compute the 2's complement value of the number of leading Z/O
	o<<tab<<"c2nZO <= ("<<wOutLZOC_<<" downto 0 => '1') - (\"0\" & nZO) + '1';"<<endl;
	
	//determine if the value to be added to the exponent bias is negative or positive.
	o<<tab<<"expTest <= '1' when (CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<wOutLZOC_<<")>nZO) else '0';"<<endl;
	
	//compute the exponent operations in both cases => may need pipelining but are generally short additions	
	//o<<tab<<"posExpAdj <= CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<wOutLZOC_+1<<")- (\"0\" & nZO); "<<endl;
	o<<tab<<"negExpAdj <= CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<wOutLZOC_+1<<") + c2nZO;"<<endl;
	
	//select the actual value which will be added to the expoent bias
	//o<<tab<<"expAdj <= posExpAdj when expTest='1' else"<<endl;
	//o<<tab<<"          negExpAdj;"<<endl;					
	
	o<<tab<<"expAdj <= negExpAdj;"<<endl;
		
	//the exponent bias	
	o<<tab<<"expR <= CONV_STD_LOGIC_VECTOR("<<expBias_<<","<<wEOut_<<");"<<endl;
	
	
	if (wOutLZOC_+1 < wEOut_){
		o<<tab<<"expRAdjusted <= expR + (("<<wEOut_-1<<" downto "<<wOutLZOC_+1<<" => expAdj("<<wOutLZOC_<<")) & expAdj);"<<endl;
		o<<tab<<"excBits <=\"01\";"<<endl;
	}
	else
		if (wOutLZOC_+1 == wEOut_){
			o<<tab<<"expRAdjusted <= expR + expAdj;"<<endl;
			o<<tab<<"excBits <=\"01\";"<<endl;
		}	
		else
		{
			o<<tab<<"expAdjustedExt <= (("<<wOutLZOC_<<" downto "<<wEOut_<<" =>'0') & expR) + expAdj;"<<endl;
			
			o<<tab<<"signExpAdjustedExt <= expAdjustedExt("<<wOutLZOC_<<");"<<endl;
			o<<tab<<"modulusExpAdjusted <= (("<<wOutLZOC_-1<<" downto 0 => signExpAdjustedExt) xor expAdjustedExt("<<wOutLZOC_-1<<" downto 0)) + signExpAdjustedExt;"<<endl;
				
			o<<tab<<"maxExponent <= CONV_STD_LOGIC_VECTOR("<<(1<<wEOut_)-1<<" , "<<wOutLZOC_<<");"<<endl;
			
			o<<tab<<"expOverflow <= '1' when (modulusExpAdjusted>maxExponent) and (signExpAdjustedExt='0') else '0';"<<endl;
			o<<tab<<"expUnderflow <= '1' when (modulusExpAdjusted>maxExponent) and (signExpAdjustedExt='1') else '0';"<<endl;
			
			o<<tab<<"excBits(1) <= expOverflow and  not(expUnderflow);"<<endl;
			o<<tab<<"excBits(0) <= not(expOverflow) and not(expUnderflow);"<<endl;
					
			o<<tab<<"expRAdjusted<=expAdjustedExt("<<wEOut_-1<<" downto 0);"<<endl;
	}			
		
	// shift 
	o<<tab<< "left_shifter_component: " << leftShifter_->getOperatorName() << endl;
	o<<tab<< "      port map ( X => "<<getDelaySignalName("pipeA",leadZOCounter_->getPipelineDepth())<<", "<< endl; 
	o<<tab<< "                 S => nZO, " << endl; 
	o<<tab<< "                 R => resFrac, " <<endl; 
	o<<tab<< "                 clk => clk, " << endl;
	o<<tab<< "                 rst => rst " << endl;
	o<<tab<< "               );" << endl<<endl;		
		
	o<<tab<<"excRes <= "<<getDelaySignalName("excBits", leftShifter_->getPipelineDepth() ) <<";"<<endl;	
	o<<tab<< "R <= excRes & "
			 <<getDelaySignalName("resultSign0",leadZOCounter_->getPipelineDepth())<<" & "
			 <<getDelaySignalName("expRAdjusted",leftShifter_->getPipelineDepth())<<"("<<wEOut_-1<<" downto 0) & ";

	if (sizeAcc_-1-wFOut_>=0)		 
			o <<"resFrac("<<sizeAcc_-2<<" downto "<<sizeAcc_-1-wFOut_<<");"<<endl;
	else
			o <<"resFrac("<<sizeAcc_-2<<" downto "<<0<<")& CONV_STD_LOGIC_VECTOR(0,"<<-(sizeAcc_-1-wFOut_)<<");"<<endl;
	
	endArchitecture(o);
		
}



