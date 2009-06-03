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
	ostringstream name;
	name << "IntAdder_" << wIn_<<"_f"<<target->frequencyMHz();
	setName(name.str());

	// Set up the IO signals
	addInput ("X"  , wIn_);
	addInput ("Y"  , wIn_);
	addInput ("Cin", 1  );
	addOutput("R"  , wIn_);

	if (verbose)
		cout << printInputDelays( inputDelays );

	if (isSequential()){
		maxInputDelay = getMaxInputDelays (inputDelays);
		if (verbose)
			cout << "The maximum input delay is: " << maxInputDelay << endl;
	
		double objectivePeriod = 1 / target->frequency();
		
		if ( maxInputDelay > objectivePeriod ){
			//the maximum combinatorial delay of the input is larger than the objective period, so the requested frequency might not be reached.
			cout << "WARNING: the combinatorial delay at the input of " << this->getName() << " is above objective period "<<endl;
			maxInputDelay = objectivePeriod;
		}

		if ( ((objectivePeriod - maxInputDelay) - target->adderDelay(1) ) < 0 )	{
			nextCycle();////////////////////////////////////////////////////////
			target->suggestSubaddSize(chunkSize_ ,wIn_);
			if (verbose) cout << "The suggested subaddition chunk size is: " << chunkSize_ << endl;
			nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
			cSize = new int[nbOfChunks+1];
			for (int i=0; i <= nbOfChunks-2; i++)
				cSize[i]=chunkSize_;				
			cSize[nbOfChunks-1]=( ((wIn_%chunkSize_)==0)?chunkSize_:wIn_-(nbOfChunks-1)*chunkSize_);
		}
		else{

			int maxInAdd;
			//the size of the first addition is dependent on the input slack
			target->suggestSlackSubaddSize(maxInAdd, wIn_, maxInputDelay); 
						 			
			int firstChunkSize = ( maxInAdd <= wIn_? maxInAdd : wIn_ );
			if ( wIn_ - firstChunkSize > 0){
				int remainingBits = wIn_ - firstChunkSize;
				target->suggestSubaddSize( chunkSize_, remainingBits);
				nbOfChunks = ceil( double(remainingBits)/double(chunkSize_));
				nbOfChunks++;
				cSize = new int[nbOfChunks];
				cSize[0] = firstChunkSize;
				cSize[nbOfChunks-1]=( (( remainingBits % chunkSize_)==0) ? chunkSize_ : remainingBits - (nbOfChunks-2)*chunkSize_ );
				for (int i=1; i <= nbOfChunks-2; i++)
					cSize[i] = chunkSize_;				
							
			}
			else{
				nbOfChunks = 1;
				cSize = new int[1];
				cSize[0] = firstChunkSize;
			}
		}

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
			vhdl << tab << declare (join("sum_l",0,"_idx",i), cSize[i]+1) << " <= " << "( \"0\" & " << use("X") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ") + "
			                                                             << "( \"0\" & " << use("Y") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ")" ;
			if (i==0) vhdl << " + " << use("Cin");
			vhdl << ";" << endl;
		}
		
		for (int i=1; i <= nbOfChunks-1 ; i++){
			nextCycle(); ///////////////////////////////////////////////////////
			for (int j=i; j <= nbOfChunks-1; j++){
				vhdl << tab << declare(join("sum_l",i,"_idx",j), cSize[j]+1) << " <= " << "( \"0\" & " << use(join("sum_l",i-1,"_idx",j))<<range(cSize[j]-1,0) << ") + "
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


