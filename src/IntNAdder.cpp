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
#include "IntNAdder.hpp"

using namespace std;

IntNAdder::IntNAdder(Target* target, int wIn, int N, map<string, double> inputDelays):
Operator(target), wIn_(wIn), N_(N), inputDelays_(inputDelays) 
{
	ostringstream name;
	name << "IntNAdder_" << wIn_;
	setName(name.str());

	// Set up the IO signals
	for (int i=0; i<N; i++){
		name.str(""); //init a ostringstream variable
		name << "X"<<i; 
		addInput (name.str() , wIn_);
	}

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
			cout << "Warning, the combinatorial delay at the input of "<<this->getName()<<"is above limit"<<endl;
			maxInputDelay = objectivePeriod;
		}

		if (((objectivePeriod - maxInputDelay) - target->lutDelay())<0)	{
			bufferedInputs = 1;
			maxInputDelay=0;
			target->suggestSubaddSize(chunkSize_ ,wIn_);
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
				target->suggestSubaddSize(chunkSize_,newWIn);
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
		
				//=================================================
		if (N>=2){
		//split the inputs ( this should be reusable )
			for (int i=0;i<N;i++)
				for (int j=0; j<nbOfChunks; j++){
					ostringstream name;
					//the naming standard: sX j _ i _ l
					//j=the chunk index i is the input index and l is the current level
					name << "sX"<<j<<"_"<<i<<"_l"<<0;
					int low=0, high=0;
					for (int k=0;k<=j;k++)
						high+=cSize[k];
					for (int k=0;k<=j-1;k++)
						low+=cSize[k];
					vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << " \"0\" & X"<<i<<range(high-1,low)<<";"<<endl;
				}
			
			////////////////////////////////////////////////
			int currentLevel = 1;
			for (int l=1; l<N; l++){
				//perform addition round; there is only one additon per round 
				for (int j=0; j<nbOfChunks; j++){
					ostringstream dname, uname1, uname2, uname3;
					dname << "sX"<<j<<"_0_l"<<l;
					uname1 << "sX"<<j<<"_0_l"<<l-1;
					uname2 << "sX"<<j<<"_1_l"<<l-1;
					uname3 << "sX"<<j-1<<"_0_l"<<l-1;
					vhdl << tab << declare(dname.str(),cSize[j]+1) << " <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<")";
					if ((j==0)&&(l==1)) vhdl << " + Cin";
					if (j>0) vhdl << " + " << use(uname3.str())<<"("<<cSize[j-1]<<") ";
					vhdl << ";" << endl;				}
				
				//from this point we just add two numbers with internal propagations, so perform addition and then take care of the propagations in a loop-like manner
				for (int propL=2; propL<=N-l; propL++)
					for (int j=0; j<nbOfChunks; j++){
						ostringstream dname, uname;
						dname << "sX"<<j<<"_"<<propL-1<<"_l"<<l;
						uname << "sX"<<j<<"_"<<propL<<"_l"<<l-1;
						vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= " << use(uname.str()) << ";" <<endl;
					}
					nextCycle();
			currentLevel++;
			}
			////////////////////////////////////////////////
			
			for (int i=2; i<nbOfChunks+1; i++){
				for (int j=i-1; j< nbOfChunks ; j++){
					ostringstream dname, uname1, uname2;
					dname << "sX"<<j<<"_0_l"<<currentLevel;
					uname1 << "sX"<<j<<"_0_l"<<currentLevel-1;
					uname2 << "sX"<<j-1<<"_0_l"<<currentLevel-1;
					vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= ( \"0\" & " << use(uname1.str())<<range(cSize[j]-1,0) << ") + " << use(uname2.str()) << "(" << cSize[j-1] << ")" << ";" <<endl;
				}
			currentLevel++;
			if (i<nbOfChunks) nextCycle();
			}
			currentLevel--;
			vhdl << tab << "R <= ";
			int k=0;
			for (int i=nbOfChunks-1; i>=0; i--){
				ostringstream uname;
				uname << "sX"<<i<<"_0_l"<<currentLevel-k;
				vhdl << use(uname.str()) << range(cSize[i]-1,0);
				if (i > 0) vhdl << " & ";
				k++;
			}
			vhdl << ";" <<endl;
		}
	}
}

IntNAdder::~IntNAdder() {
}


void IntNAdder::emulate(TestCase* tc)
{
	mpz_class svX[N_];
	for (int i=0; i<N_; i++){
	ostringstream iName;
		iName << "X"<<i;
		svX[i] = tc->getInputValue(iName.str());
	}
	mpz_class svC =  tc->getInputValue("Cin");

	mpz_class svR = svX[0] + svC;
	mpz_clrbit(svR.get_mpz_t(),wIn_); 
	for (int i=1; i<N_; i++){
		svR = svR + svX[i];
		mpz_clrbit(svR.get_mpz_t(),wIn_); 
	}

	tc->addExpectedOutput("R", svR);

}


