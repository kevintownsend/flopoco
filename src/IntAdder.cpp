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


namespace flopoco{

	IntAdder::IntAdder(Target* target, int wIn, map<string, double> inputDelays, int aType):
		Operator(target), wIn_(wIn), inputDelays_(inputDelays)
	{
		ostringstream name;
		setCopyrightString("Bogdan Pasca, Florent de Dinechin (2008-2009)");		
		name << "IntAdder_" << wIn_<<"_f"<<target->frequencyMHz();
		setName(name.str());

		// Set up the IO signals
		addInput ("X"  , wIn_, true);
		addInput ("Y"  , wIn_, true);
		addInput ("Cin", 1 );
		addOutput("R"  , wIn_, 1 , true);

		if (verbose) 
			cout << printInputDelays( inputDelays );

		if (isSequential()){
			//general variables for all versions
			double objectivePeriod = 1 / target->frequency();	
		
		
			//**********************************************************************
			//**********************************************************************
			//FASTEST AND SMALLEST PIPELINED VERSION OF CPA1
			if ( aType == 0 ){

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
				
					if (scoreFirstDesign > scoreSecondDesign){
						if (verbose) cout << "Implementation of the second design" << endl;
						nbOfChunks = nbOfChunksSecondDesign;
						target->suggestSubaddSize(chunkSize_, wIn_); 
						lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
					}else{
						if (verbose) cout << "Implementation of the first design" << endl;
						nbOfChunks = nbOfChunksFirstDesign;
						target->suggestSubaddSize(chunkSize_, wIn_); 
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

				////////////////////////////////////////////////////////////////////////
		
				for (int i=0; i < nbOfChunks; i++){
					vhdl << tab << declare (join("sum_l",0,"_idx",i), cSize[i]+1, true) << " <= " << "( \"0\" & " << use("X") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ") + "
						  << "( \"0\" & " << use("Y") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ")" ;
					if (i==0) vhdl << " + " << use("Cin");
					vhdl << ";" << endl;
				}
		
				for (int i=1; i <= nbOfChunks-1 ; i++){
					nextCycle(); ///////////////////////////////////////////////////////
					for (int j=i; j <= nbOfChunks-1; j++){
						vhdl << tab << declare(join("sum_l",i,"_idx",j), cSize[j]+1, true) << " <= " << "( \"0\" & " << use(join("sum_l",i-1,"_idx",j))<<range(cSize[j]-1,0) << ") + "
							  << use(join("sum_l",i-1,"_idx",j-1))<<of(cSize[j-1])<<";"<<endl;
					}
				}
		
				vhdl << tab << "R <= ";
				for (int i=nbOfChunks-1; i >= 1; i--){
					vhdl << use(join("sum_l",i,"_idx",i))<<range(cSize[i]-1,0)<< " & ";
				}
				vhdl << use("sum_l0_idx0")<<range(cSize[0]-1,0)<<";"<<endl;

				outDelayMap["R"] = target->adderDelay(cSize[nbOfChunks-1]); 
				if (verbose)
					cout<< "Last addition size is "<<cSize[nbOfChunks-1]<< " having a delay of "<<target->adderDelay(cSize[nbOfChunks-1])<<endl;
			}
			//**********************************************************************
			//**********************************************************************
			//SECOND BEST PIPELINED VERSION OF CPA1
			else if ( aType == 1 ){
				int selectedDesign;
				int firstChunkSize, middleChunkSize;
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
					lastChunkSize = wIn_ % chunkSize_;
				}
				else{
					//explore 2 designs and chose the best	
					int maxInAdd;
								
					int typeOfChunks = 1;
								
					target->suggestSlackSubaddSize(firstChunkSize, wIn_, maxInputDelay); 
					if (wIn_-firstChunkSize > 0){
						target->suggestSubaddSize(middleChunkSize, wIn_); 
						if ( wIn_-firstChunkSize < middleChunkSize ) 
							typeOfChunks++;
						else
							typeOfChunks+=2;
					
						lastChunkSize = ( (wIn_-firstChunkSize) % middleChunkSize == 0 ? middleChunkSize : (wIn_-firstChunkSize) % middleChunkSize );
					}
					int nbOfChunksFirstDesign = 1;
				
					if (typeOfChunks > 1){
						if (typeOfChunks > 2)
							nbOfChunksFirstDesign = 2 + (wIn_-firstChunkSize-lastChunkSize)/middleChunkSize;
						else
							nbOfChunksFirstDesign = typeOfChunks;
					}
					int scoreFirstDesign = nbOfChunksFirstDesign - 1;
					if (verbose) cout << "Exploring first design ... score is:"<< scoreFirstDesign << endl;
				
					target->suggestSubaddSize(maxInAdd, wIn_); 
					int nbOfChunksSecondDesign = ceil(double(wIn_)/double(maxInAdd));
					int scoreSecondDesign = nbOfChunksSecondDesign-1;
					if (verbose) cout << "Exploring second design ... score is:"<< scoreSecondDesign << endl;
				
					if (scoreFirstDesign > scoreSecondDesign){
						selectedDesign = 2;
						if (verbose) cout << "Implementation of the second design" << endl;
						nbOfChunks = nbOfChunksSecondDesign;
						target->suggestSubaddSize(chunkSize_, wIn_); 
						lastChunkSize = ( wIn_ % chunkSize_ ==0 ? chunkSize_ : wIn_ % chunkSize_) ;
					}else{
						selectedDesign = 1;
						if (verbose) cout << "Implementation of the first design" << endl;
						nbOfChunks = nbOfChunksFirstDesign;
						target->suggestSubaddSize(chunkSize_, wIn_); 
					}
				}
				//the sizes of the chunks
				cSize = new int[nbOfChunks+1];
				if ( nbOfChunks > 1 ){
					if (selectedDesign ==1)
						cSize[0] = firstChunkSize;
					else
						cSize[0] = chunkSize_;
					
					for (int i=1; i<nbOfChunks-1; i++)
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

				//TODO - fill-in the implementation
				for (int i=0; i < nbOfChunks; i++){
					vhdl << tab << declare( join("x",i), cSize[i],true) << " <= " << use("X") << range(cIndex[i]-1,(i>0?cIndex[i-1]:0)) << ";" << endl;
					vhdl << tab << declare( join("y",i), cSize[i],true) << " <= " << use("Y") << range(cIndex[i]-1,(i>0?cIndex[i-1]:0)) << ";" << endl;
				}
			
				for (int i=0; i < nbOfChunks; i++){
					vhdl << tab << declare(join("sum",i),cSize[i]+1,true) << " <= " << "( \"0\" & "<< use(join("x",i)) << ") + " << "( \"0\" & "<< use(join("y",i)) << ")  + ";
					if (i==0)
						vhdl << use("Cin");
					else
						vhdl << use(join("sum",i-1))<<of(cSize[i-1]);
					vhdl << ";"	<< endl;
				
					if (i < nbOfChunks-1)
						nextCycle();
				}
	
				vhdl << tab << "R <= ";
				for (int i=nbOfChunks-1; i >= 1; i--){
					vhdl << use(join("sum",i))<<range(cSize[i]-1,0)<< " & ";
				}
				vhdl << use("sum0")<<range(cSize[0]-1,0)<<";"<<endl;
				outDelayMap["R"] = target->adderDelay(cSize[nbOfChunks-1]); 
			}
			else{}
		}else{
			vhdl << tab << " R <= X + Y + Cin;" << endl;
			outDelayMap["R"] = target->adderDelay(wIn_); 
		}
	
	
	}

	IntAdder::~IntAdder() {
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


}
