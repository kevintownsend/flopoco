/*
 * An integer adder for FloPoCo
 *
 * It may be pipelined to arbitrary frequency.
 * Also useful to derive the carry-propagate delays for the subclasses of Target
 *
 * Authors : Florent de Dinechin, Bogdan Pasca
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
#include "IntAdder.hpp"

using namespace std;

IntAdder::IntAdder(Target* target, int wIn, map<string, double> inputDelays):
Operator(target), wIn_(wIn), inputDelays_(inputDelays)
{
	setOperatorName();
	setOperatorType();

	// Set up the IO signals
	addInput ("X"  , wIn_);
	addInput ("Y"  , wIn_);
	addInput ("Cin", 1  );
	addOutput("R"  , wIn_);

	if (verbose){
		cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
		cout <<"delay for Y is   "<< inputDelays["Y"]<<endl;
		cout <<"delay for Cin is "<< inputDelays["Cin"]<<endl;
	}

	if (isSequential()){
		//compute the maximum input delay
		maxInputDelay = 0;
		map<string, double>::iterator iter;
		for (iter = inputDelays.begin(); iter!=inputDelays.end();++iter)
			if (iter->second > maxInputDelay)
				maxInputDelay = iter->second;
	
		if (verbose)
			cout << "The maximum input delay is "<<	maxInputDelay<<endl;
	
		double	objectivePeriod;
		objectivePeriod = 1/ target->frequency();
		
		if (verbose)
			cout << "Objective period is "<< objectivePeriod<<" at an objective frequency of "<<target->frequency() << endl;

		if (objectivePeriod<maxInputDelay){
			//It is the responsability of the previous components to not have a delay larger than the period
			cout << "Warning, the combinatorial delay at the input of "<<this->getOperatorName()<<"is above limit"<<endl;
			maxInputDelay = objectivePeriod;
		}
		
		if (((objectivePeriod - maxInputDelay) - target->lutDelay())<0)	{
			bufferedInputs = 1;
			maxInputDelay=0;
			//bool status = target->suggestSubaddSize(chunkSize_ ,wIn_);
			nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
			cSize = new int[nbOfChunks+1];
			cSize[nbOfChunks-1]=( ((wIn_%chunkSize_)==0)?chunkSize_:wIn_-(nbOfChunks-1)*chunkSize_);
			for (int i=0;i<=nbOfChunks-2;i++)
				cSize[i]=chunkSize_;				
		}
		else{
			int cS0; 
			bufferedInputs=0;
			int maxInAdd = ceil(((objectivePeriod - maxInputDelay) - target->lutDelay())/target->carryPropagateDelay()); 			
			cS0 = (maxInAdd<=wIn_?maxInAdd:wIn_);
			if ((wIn_-cS0)>0)
			{
				int newWIn = wIn_-cS0;
				//bool status = target->suggestSubaddSize(chunkSize_,newWIn);
				nbOfChunks = ceil( double(newWIn)/double(chunkSize_));
				cSize = new int[nbOfChunks+1];
				cSize[0] = cS0;
				cSize[nbOfChunks]=( (( (wIn_-cSize[0])%chunkSize_)==0)?chunkSize_:(wIn_-cSize[0])-(nbOfChunks-1)*chunkSize_);
				for (int i=1;i<=nbOfChunks-1;i++)
					cSize[i]=chunkSize_;				
				nbOfChunks++;			
			}
			else{
				nbOfChunks=1;
				cSize = new int[1];
				cSize[0] = cS0;
			}
		}
		
		if (verbose){
			cout<<endl<<endl<<"Buffered Inputs "<<(bufferedInputs?"yes":"no")<<endl;
			cout<<endl;
			for (int i=nbOfChunks-1;i>=0;i--)
				cout<<cSize[i]<<" ";
			cout<<endl; 
		}	
		
		for (int i=0;i<nbOfChunks;i++){
			ostringstream t; t<<"X"<<i;	
			addDelaySignalBus(t.str(),cSize[i],bufferedInputs); t.str(""); t<<"Y"<<i;
			addDelaySignalBus(t.str(),cSize[i],bufferedInputs); t.str(""); 
			if (i==0)
				addDelaySignal("Carry",1,bufferedInputs); 
		}
		
		for (int i=0; i<nbOfChunks-1;i++){
			ostringstream t;
			t<<"cin"<<i+1<<"R"<<i;
			addDelaySignal(t.str(),cSize[i]+1,1);
		}
		
		for (int i=0; i<nbOfChunks;i++){
			ostringstream t;
			t<<"R"<<i;
			addDelaySignalBus(t.str(),cSize[i],nbOfChunks-2-i);
		}	
		
		for (int i=0; i<nbOfChunks;i++){
			ostringstream t; t<<"sX"<<i;
			addDelaySignalBus(t.str(),cSize[i],i); t.str(""); t<<"sY"<<i;
			addDelaySignalBus(t.str(),cSize[i],i); t.str("");
			if (i==0)
			addSignal("cin0",1);
		}	

		setPipelineDepth(nbOfChunks-1+bufferedInputs);

		outDelayMap["R"] = target->adderDelay(cSize[nbOfChunks-1]); 
		if (verbose)
			cout<< "Last addition size is "<<cSize[nbOfChunks-1]<< " having a delay of "<<target->adderDelay(cSize[nbOfChunks-1])<<endl;

	}
	
}

IntAdder::~IntAdder() {
}

void IntAdder::setOperatorName(){
	ostringstream name;
	name << "IntAdder_" << wIn_;
	uniqueName_ = name.str();
}

void IntAdder::outputVHDL(std::ostream& o, std::string name) {
	ostringstream signame;
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007, 2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	if(isSequential()){
		//first assignments. This level might be transformed from signals to registers
		// if the (T-inputDelay)<lutDelay
		for (int i=0;i<nbOfChunks;i++){
			if (i==0)
				o << tab << "Carry <= Cin; "<<endl;
			int sum=0;
			for (int j=0;j<=i;j++) sum+=cSize[j];
			o << tab << "X"<<i<<" <= X("<<sum-1<<" downto "<<sum-cSize[i]<<");"<<endl;
			o << tab << "Y"<<i<<" <= Y("<<sum-1<<" downto "<<sum-cSize[i]<<");"<<endl;
		}
		//connect first assignments to second level of signals
		for (int i=0;i<nbOfChunks;i++){
			ostringstream xi,yi;
			xi << "X"<<i;
			yi << "Y"<<i;
			o << tab << "sX"<<i<<" <= " << delaySignal(xi.str(), bufferedInputs)<<";"<<endl;
			o << tab << "sY"<<i<<" <= " << delaySignal(yi.str(), bufferedInputs)<<";"<<endl;
			if (i==0)
			o << tab << "cin0 <= "<<delaySignal("Carry",bufferedInputs)<<";"<<endl;
		}

		//additions		
		for (int i=0;i<nbOfChunks;i++){
			ostringstream sxi,syi;
			sxi << "sX"<<i;
			syi << "sY"<<i;
			if (i==0 && nbOfChunks>1)
				o << tab << "cin"<<i+1<<"R"<<i<<" <= (\"0\" & sX"<<i<<") + (\"0\" & sY"<<i<<") + cin0;"<<endl;
			else 
				if (i<nbOfChunks-1)
					o << tab << "cin"<<i+1<<"R"<<i<<" <= ( \"0\" & " << delaySignal(sxi.str(), i) << ")"
					  << " + ( \"0\" & " << delaySignal(syi.str(),i)<< ")"
					  << " + cin"<<i<<"R"<<i-1<<"_d("<<cSize[i-1]<<");"<<endl;
		}

		//assign the partial additions which will propagate to the result
		for (int i=0;i<nbOfChunks;i++){
			ostringstream sxi,syi;
			sxi << "sX"<<i;
			syi << "sY"<<i;
			if (i<nbOfChunks-1)
				o << tab << "R"<<i<<" <= cin"<<i+1<<"R"<<i<<"_d("<<cSize[i]-1<<" downto 0);"<<endl;
			else{
				o << tab << "R"<<i<<" <= "<< delaySignal(sxi.str(), i)
				  << " + " << delaySignal(syi.str(), i);
				if (nbOfChunks>1)				
					o << " + cin"<<i<<"R"<<i-1<<"_d("<<cSize[i-1]<<");"<<endl;
				else
					o << " + cin0;"<<endl;
			}
		}
		
		//assign output by composing the result
		o << tab << "R <= ";
		for (int i=nbOfChunks-1;i>=0;i--){
			ostringstream ri;
			ri << "R"<<i;
			if (i==0)
				o << delaySignal(ri.str(), nbOfChunks-2-i)<<";"<<endl;
			else
				o << delaySignal(ri.str(),nbOfChunks-2-i)<<" & ";			
		} 
		o<<endl;

		outputVHDLRegisters(o);
	}
	else{
		o << tab << "R <= X + Y + Cin;" <<endl;
	}
	o << "end architecture;" << endl << endl;
}


void IntAdder::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");
	mpz_class svC = tc->getInputValue("Cin");

	mpz_class svR = svX + svY + svC;
	// Don't allow overflow
	mpz_clrbit(svR.get_mpz_t(),wIn_); 

	tc->addExpectedOutput("R", svR);
}


