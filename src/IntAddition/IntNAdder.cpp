/*
  An integer multi-operand adder for FloPoCo
 
  Authors : Bogdan Pasca
 
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "IntNAdder.hpp"

using namespace std;

namespace flopoco{

	IntNAdder::IntNAdder(Target* target, int wIn, int N, map<string, double> inputDelays):
		Operator(target, inputDelays), wIn_(wIn), N_(N), inputDelays_(inputDelays) 
	{
	
		ostringstream name;
		name << "IntNAdder_" << wIn_<<"_"<<N<<"_uid"<<Operator::getNewUId();
		srcFileName = "IntNAdder";
		setName(name.str());

		setCopyrightString("Bogdan Pasca (2009, 2010)");

		// Set up the IO signals
		for (int i=0; i<N; i++)
			addInput ( join("X",i) , wIn_, true);

//		addInput ("Cin", 1  );
		addOutput("R"  , wIn_, 1, true);

		if (isSequential()){
	
			double objectivePeriod = 1 / target->frequency();	
			int lastChunkSize;
			int *cIndex;                       /**< array containing the indexes for all Chunks*/
			maxInputDelay = getMaxInputDelays (inputDelays);
			REPORT(DEBUG, "The maximum input delay is: " << maxInputDelay);
	
			if ( maxInputDelay > objectivePeriod ){
				//the maximum combinatorial delay of the input is larger than the objective period, so the requested frequency might not be reached.
				REPORT(INFO, "WARNING: the combinatorial delay at the input of " << this->getName() << " is above objective period ");
				maxInputDelay = objectivePeriod;
				nextCycle();
				setCriticalPath(0.0);
				target->suggestSubaddSize(chunkSize_ ,wIn_);
				nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
				lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
			} else if ( ((objectivePeriod - maxInputDelay) - target->adderDelay(1) ) < 0 )	{
				//if not enough slack is available for any combinatorial circuit, we register the inputs
				nextCycle();
				setCriticalPath(0.0);
				target->suggestSubaddSize(chunkSize_ ,wIn_);
				nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
				lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
			} else{
				int maxInAdd;
				//explore 2 designs and chose the best				
				target->suggestSlackSubaddSize(maxInAdd, wIn_, maxInputDelay); 
				int nbOfChunksFirstDesign = ceil(double(wIn_)/double(maxInAdd));
				int scoreFirstDesign = nbOfChunksFirstDesign - 1;
				REPORT(DEBUG, "Exploring first design ... score is:"<< scoreFirstDesign);
				target->suggestSubaddSize(maxInAdd, wIn_); 
				int nbOfChunksSecondDesign = ceil(double(wIn_)/double(maxInAdd));
				int scoreSecondDesign = nbOfChunksSecondDesign;
				REPORT(DEBUG, "Exploring second design ... score is:"<< scoreSecondDesign);
				if (scoreFirstDesign > scoreSecondDesign){
					REPORT(DEBUG, "Implementation of the second design");
					nextCycle();//
					setCriticalPath(0.0);
					nbOfChunks = nbOfChunksSecondDesign;
					target->suggestSubaddSize(chunkSize_, wIn_); 
					lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
				}else{
					REPORT(DEBUG, "Implementation of the first design");
					setCriticalPath(maxInputDelay);
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
			
			if (wIn>1){ 
				if (N>=2){
					//split the inputs (this should be reusable)
					for (int i=0;i<N;i++)
						for (int j=0; j<nbOfChunks; j++){
							ostringstream name;
							//the naming standard: sX j _ i _ l
							//j=the chunk index i is the input index and l is the current level
							name << "sX"<<j<<"_"<<i<<"_l"<<0;
							int low=0, high=0;
							for (int k=0;k<=j  ;k++)  	high+=cSize[k];
							for (int k=0;k<=j-1;k++)	low +=cSize[k];
							if(high-1<=wIn)
								vhdl << tab << declare (name.str(),cSize[j]+1) << " <=  \"0\" & X"<<i<<range(high-1,low)<<";"<<endl;
							else{
								if(low<wIn)
									vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << zg(cSize[j]-(wIn -low)+1,0) <<" & X"<<i<<range(wIn-1,low)<<";"<<endl;
								else
									if(low==wIn)
										vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << zg(cSize[j],0) <<" & X"<<i<<of(wIn-1)<<";"<<endl;
									else
										vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << zg(cSize[j]+1,0) <<";"<<endl;
							}
						}
			
					int currentLevel = 1;
					for (int l=1; l<N; l++){
						//perform addition round; there is only one additon per round
						manageCriticalPath(target->adderDelay(cSize[0]));
						 
						for (int j=0; j<nbOfChunks; j++){
							ostringstream dname, uname1, uname2, uname3;
							dname << "sX"<<j<<"_0_l"<<l;
							uname1 << "sX"<<j<<"_0_l"<<l-1;
							uname2 << "sX"<<j<<"_1_l"<<l-1;
							uname3 << "sX"<<j-1<<"_0_l"<<l-1;
							vhdl << tab << declare(dname.str(),cSize[j]+1) << " <= (\"0\" & "<< uname1.str()<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< uname2.str()<<range(cSize[j]-1,0)<<")";
//							if ((j==0)&&(l==1)) vhdl << " + Cin";
							if (j>0) vhdl << " + " << uname3.str()<<"("<<cSize[j-1]<<") ";
							vhdl << ";" << endl;				
						}
				
						//from this point we just add two numbers with internal propagations, so perform addition and then take care of the propagations in a loop-like manner
						for (int propL=2; propL<=N-l; propL++)
							for (int j=0; j<nbOfChunks; j++){
								ostringstream dname, uname;
								dname << "sX"<<j<<"_"<<propL-1<<"_l"<<l;
								uname << "sX"<<j<<"_"<<propL<<"_l"<<l-1;
								vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= " << uname.str() << ";" <<endl;
							}
				
						currentLevel++;
					}

					REPORT(DEBUG, "The number of chunks is: " << nbOfChunks);
					
					if (nbOfChunks>1){
						vhdl << tab << "--final propagations " << endl; 
			
						for (int i=2; i<nbOfChunks+1; i++){
							manageCriticalPath(target->adderDelay(cSize[i]));
							for (int j=i-1; j< nbOfChunks ; j++){
								ostringstream dname, uname1, uname2;
								dname <<  "sX"<<j<<"_0_l"<<currentLevel;
								uname1 << "sX"<<j<<"_0_l"<<currentLevel-1;
								uname2 << "sX"<<j-1<<"_0_l"<<currentLevel-1;
								vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= ( \"0\" & " << uname1.str()<<range(cSize[j]-1,0) << ") + " 
										                                                            << uname2.str()<<of(cSize[j-1])<<";" <<endl;
							}
							currentLevel++;
						}
					}
					currentLevel--;
				
					outDelayMap["R"] = getCriticalPath();
					vhdl << tab << "R <= ";
					int k=0;

					for (int i=nbOfChunks-1; i>=0; i--){
						ostringstream uname;
						uname << "sX"<<i<<"_0_l"<<currentLevel-k;
						vhdl << uname.str() << range(cSize[i]-1,0);
						if (i > 0) vhdl << " & ";
						k++;
					}
					vhdl << ";" <<endl;
				}else{
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
								vhdl << tab << declare (name.str(),cSize[j]+1) << " <=  \"0\" & X"<<i<<range(high-1,low)<<";"<<endl;
							}
	
						int currentLevel = 1;
						for (int i=2; i<nbOfChunks+2; i++){
							manageCriticalPath(target->adderDelay(cSize[i]));
							for (int j=i-2; j< nbOfChunks ; j++){
								ostringstream dname, uname1, uname2;
								dname << "sX"<<j<<"_0_l"<<currentLevel;
								uname1 << "sX"<<j<<"_0_l"<<currentLevel-1;
								uname2 << "sX"<<j-1<<"_0_l"<<currentLevel-1;
								if (j>0)
									vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= ( \"0\" & " << uname1.str()<<range(cSize[j]-1,0) << ") + " << uname2.str() << of(cSize[j-1])<<";" <<endl;
								else
									vhdl << tab << declare(dname.str(), cSize[j]+1) << " <= ( \"0\" & " << uname1.str()<<range(cSize[j]-1,0) << ");"<<endl; // + Cin ;" <<endl;
							}
							currentLevel++;
						}
						currentLevel--;
						vhdl << tab << "R <= ";
						int k=0;
						for (int i=nbOfChunks-1; i>=0; i--){
							ostringstream uname;
							uname << "sX"<<i<<"_0_l"<<currentLevel-k;
							vhdl << uname.str() << range(cSize[i]-1,0);
							if (i > 0) vhdl << " & ";
							k++;
						}
						vhdl << ";" <<endl;
						outDelayMap["R"] = getCriticalPath();
					}
				}
				
			}else{
				vhdl << tab << "R <= ";
					for (int i=N-1; i>0; i--){
						vhdl << "X"<<i<< " + ";
					} 
					vhdl << " X0 ;"<<endl;
			}
		}else{ 
			/* if combinatorial */
			vhdl << tab << "R <= ";
			for (int i=N-1; i>0; i--){
				vhdl << "X"<<i<< " + ";
			} 
			vhdl << " X0 ;"<<endl;
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
