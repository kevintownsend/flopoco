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


namespace flopoco{

	IntNAdder::IntNAdder(Target* target, int wIn, int N, map<string, double> inputDelays):
		Operator(target), wIn_(wIn), N_(N), inputDelays_(inputDelays) 
	{
		ostringstream name;
		name << "IntNAdder_" << wIn_<<"_"<<N;
		setName(name.str());

		setCopyrightString("Bogdan Pasca (2009)");

		// Set up the IO signals
		for (int i=0; i<N; i++)
			addInput ( join("X",i) , wIn_, true);

		addInput ("Cin", 1  );
		addOutput("R"  , wIn_);

		if (isSequential()){
			double objectivePeriod = 1 / target->frequency();	
			int lastChunkSize;
			int *cIndex;                       /**< array containing the indexes for all Chunks*/
			maxInputDelay = getMaxInputDelays (inputDelays);
			if (verbose)
				cout << "The maximum input delay is: " << maxInputDelay << endl;
	
			if ( maxInputDelay > objectivePeriod ){
				//the maximum combinatorial delay of the input is larger than the objective period, so the requested frequency might not be reached.
				cout << "WARNING: the combinatorial delay at the input of " << this->getName() << " is above objective period "<<endl;
				maxInputDelay = objectivePeriod;
			}

			if ( ((objectivePeriod - maxInputDelay) - target->adderDelay(1) ) < 0 )	{
				//if not enough slack is available for any combinatorial circuit, we register the inputs
				nextCycle();
				target->suggestSubaddSize(chunkSize_ ,wIn_);
				nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
				lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
			}
			else{
				int maxInAdd;
				//explore 2 designs and chose the best				
				target->suggestSlackSubaddSize(maxInAdd, wIn_, maxInputDelay); 
				int nbOfChunksFirstDesign = ceil(double(wIn_)/double(maxInAdd));
				int scoreFirstDesign = nbOfChunksFirstDesign - 1;
				if (verbose) cout << "Exploring first design ... score is:"<< scoreFirstDesign << endl;
			
				target->suggestSubaddSize(maxInAdd, wIn_); 
				int nbOfChunksSecondDesign = ceil(double(wIn_)/double(maxInAdd));
				int scoreSecondDesign = nbOfChunksSecondDesign;
				if (verbose) cout << "Exploring second design ... score is:"<< scoreSecondDesign << endl;
			
				if ((scoreFirstDesign > scoreSecondDesign) &&
					 (maxInputDelay <= 0)) // this expresion was added to ensure that the implemented design will have the necessary input delay
					{
						if (verbose) cout << "Implementation of the second design" << endl;
						nbOfChunks = nbOfChunksSecondDesign;
						target->suggestSubaddSize(chunkSize_, wIn_); 
						lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
					}else{
					if (verbose) cout << "Implementation of the first design" << endl;
					nbOfChunks = nbOfChunksFirstDesign;
					target->suggestSlackSubaddSize(chunkSize_, wIn_, maxInputDelay); 
					lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
				}
			}
			//the sizes of the chunks
			cSize = new int[nbOfChunks+1];
			if ( nbOfChunks > 1 ){
				for (int i=0; i<nbOfChunks-1; i++)
					cSize[i] = chunkSize_;
				cSize[nbOfChunks-1] = lastChunkSize;
			}
			else{
				nbOfChunks = 1;
				cSize = new int[1];
				cSize[0] = wIn_;
			}
			//the indexes in the inputs of the chunks
			cIndex = new int[nbOfChunks];
			cIndex[0]= cSize[0];
			for (int i=1; i < nbOfChunks; i++)
				cIndex[i] = cIndex[i-1] + cSize[i];
	
			if (verbose){
				cout << "The chunk sizes[MSB-->LSB]: "<<endl;
				for (int i=nbOfChunks-1;i>=0;i--)
					cout<<cSize[i]<<" ";
				cout<<endl; 
				cout << "The index sizes[MSB-->LSB]: "<<endl;
				for (int i=nbOfChunks-1;i>=0;i--)
					cout<<cIndex[i]<<" ";
				cout<<endl; 
			}	
		
			//=================================================
			if (N>=2){
				double delay = 0.0;
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
						if(high-1<=wIn)
							vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << " \"0\" & X"<<i<<range(high-1,low)<<";"<<endl;
						else
							{
								if(low<wIn)
									vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << zg(cSize[j]-(wIn -low)+1,0) <<" & X"<<i<<range(wIn-1,low)<<";"<<endl;
								else
									if(low==wIn)
										vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << zg(cSize[j],0) <<" & X"<<i<<of(wIn-1)<<";"<<endl;
									else
										vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << zg(cSize[j]+1,0) <<";"<<endl;
							}
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
						vhdl << ";" << endl;				
					}
				
					//from this point we just add two numbers with internal propagations, so perform addition and then take care of the propagations in a loop-like manner
					for (int propL=2; propL<=N-l; propL++)
						for (int j=0; j<nbOfChunks; j++){
							ostringstream dname, uname;
							dname << "sX"<<j<<"_"<<propL-1<<"_l"<<l;
							uname << "sX"<<j<<"_"<<propL<<"_l"<<l-1;
							vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= " << use(uname.str()) << ";" <<endl;
						}
				
					if (nbOfChunks>1) { 
						currentLevel++;
						nextCycle();
					}else{
						currentLevel++;
						delay += target->adderDelay(wIn_) + target->localWireDelay();
						if (delay > objectivePeriod){
							nextCycle();
							delay = target->adderDelay(wIn_) + target->localWireDelay();
						}
					}
				}
				////////////////////////////////////////////////
				vhdl << tab << "--final propagations " << endl; 
			
				if (verbose)
					cout << "The number of chunks is: " << nbOfChunks << endl;
			
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
				//ostringstream uname;
				//uname << "sX"<<nbOfChunks-1<<"_0_l"<<currentLevel-k;
				//vhdl << use(uname.str()) << range(cSize[nbOfChunks-1]-1,0);
				//vhdl << ";" <<endl;
			
				for (int i=nbOfChunks-1; i>=0; i--){
					ostringstream uname;
					uname << "sX"<<i<<"_0_l"<<currentLevel-k;
				
					vhdl << use(uname.str()) << range(cSize[i]-1,0);
				
				
					if (i > 0) vhdl << " & ";
					k++;
				}
				vhdl << ";" <<endl;
			}else
				if (N==1){
					//split the inputs ( this should be reusable )
					for (int i=0;i<N;i++)
						for (int j=0; j<nbOfChunks; j++){
							ostringstream name;
							//the naming standard: sX j _ i _ l
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
					for (int i=2; i<nbOfChunks+2; i++){
						for (int j=i-2; j< nbOfChunks ; j++){
							ostringstream dname, uname1, uname2;
							dname << "sX"<<j<<"_0_l"<<currentLevel;
							uname1 << "sX"<<j<<"_0_l"<<currentLevel-1;
							uname2 << "sX"<<j-1<<"_0_l"<<currentLevel-1;
							if (j>0)
								vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= ( \"0\" & " << use(uname1.str())<<range(cSize[j]-1,0) << ") + " << use(uname2.str()) << "(" << cSize[j-1] << ")" << ";" <<endl;
							else
								vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= ( \"0\" & " << use(uname1.str())<<range(cSize[j]-1,0) << ") + Cin ;" <<endl;
						}
						currentLevel++;
						if (i<nbOfChunks+1) nextCycle();
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
		}else{
			vhdl << tab << "R <= ";
			for (int i=N-1; i>=0; i--){
				vhdl << "X"<<i<< " + ";
			} 
			vhdl << " Cin ;"<<endl;
		}
	}

	IntNAdder::~IntNAdder() {
	}


	void IntNAdder::emulate(TestCase* tc)
	{
		mpz_class svX;
		mpz_class svC =  tc->getInputValue("Cin");
		mpz_class svR = svC;

		for (int i=0; i<N_; i++){
			ostringstream iName;
			iName << "X"<<i;
			svX = tc->getInputValue(iName.str());
			svR = svR + svX;
			mpz_clrbit(svR.get_mpz_t(),wIn_); 
		}
	
		tc->addExpectedOutput("R", svR);
	}


}
