/*
 * An integer adder for FloPoCo
 *
 * It may be pipelined to arbitrary frequency.
 * Also useful to derive the carry-propagate delays for the subclasses of Target
 *
 * Authors : Bogdan Pasca
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

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "LongIntAdderCmpCmpAdd.hpp"
#include "IntAdder.hpp"
// #define XILINX_OPTIMIZATION 0


using namespace std;
namespace flopoco{
extern vector<Operator*> oplist;

	LongIntAdderCmpCmpAdd::LongIntAdderCmpCmpAdd(Target* target, int wIn, map<string, double> inputDelays):
		Operator(target), wIn_(wIn), inputDelays_(inputDelays) 
	{
		srcFileName="LongIntAdderCmpCmpAdd";
		setName(join("LongIntAdderCmpCmpAdd_", wIn_));
//		int version = 1; /* this will go into the parameters */
		
		// Set up the IO signals
		for (int i=0; i<2; i++)
			addInput ( join("X",i) , wIn_);
		addInput ("Cin", 1  );
		addOutput("R"  , wIn_);

		if (verbose){
			cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
			cout <<"delay for Y is   "<< inputDelays["Y"]<<endl;
			cout <<"delay for Cin is "<< inputDelays["Cin"]<<endl;
		}

		if (isSequential()){
			
			//compute the maximum input delay
			maxInputDelay = getMaxInputDelays(inputDelays);
			if (verbose)
				cout << "The maximum input delay is "<<	maxInputDelay<<endl;
			
			cSize = new int[2000];
			REPORT(3, "-- The new version: direct mapping without 0/1 padding, IntAdders instantiated");
			double	objectivePeriod = double(1) / target->frequency();
			REPORT(2, "Objective period is "<< objectivePeriod <<" at an objective frequency of "<<target->frequency());
			target->suggestSubaddSize(chunkSize_ ,wIn_);
			REPORT(2, "The chunkSize for first two chunks is: " << chunkSize_ );
			
			if (2*chunkSize_ >= wIn_){
				cerr << "ERROR FOR NOW -- instantiate int adder, dimmension too small for LongIntAdderCmpCmpAdd" << endl;
				exit(0);
			}
			
			cSize[0] = chunkSize_;
			cSize[1] = chunkSize_;
			
			bool finished = false; /* detect when finished the first the first
			phase of the chunk selection algo */
			int width = wIn_ - 2*chunkSize_; /* remaining size to split into chunks */
			int propagationSize = 0; /* carry addition size */
			int chunkIndex = 2; /* the index of the chunk for which the size is
			to be determined at the current step */
			bool invalid = false; /* the result of the first phase of the algo */
			
			/* FIRST PHASE */
			REPORT(3, "FIRST PHASE chunk splitting");
			while (not (finished))	 {
				REPORT(2, "The width is " << width);
				propagationSize+=2;
				double delay = objectivePeriod - target->adderDelay(width)- target->adderDelay(propagationSize); //2*target->localWireDelay()  -
				REPORT(2, "The value of the delay at step " << chunkIndex << " is " << delay);
				if ((delay > 0) || (width < 4)) {
					REPORT(2, "finished -> found last chunk of size: " << width);
					cSize[chunkIndex] = width;
					finished = true;
				}else{
					REPORT(2, "Found regular chunk ");
					int cs; 
					double slack =  target->adderDelay(propagationSize) ; //+ 2*target->localWireDelay()
					REPORT(2, "slack is: " << slack);
					REPORT(2, "adderDelay of " << propagationSize << " is " << target->adderDelay(propagationSize) );
					target->suggestSlackSubaddSize( cs, width, slack);
					REPORT(2, "size of the regular chunk is : " << cs);
					width = width - cs;
					cSize[chunkIndex] = cs;
					
					if ( (cSize[chunkIndex-1]<=2) && (cSize[chunkIndex-1]<=2) && ( invalid == false) ){
						REPORT(1, "[WARNING] Register level inserted after carry-propagation chain");
						invalid = true; /* invalidate the current splitting */
					}
					chunkIndex++; /* as this is not the last pair of chunks,
					pass to the next pair */
				}
			}
			REPORT(2, "First phase return valid result: " << invalid);
			
			/* SECOND PHASE: 
			only if first phase is cannot return a valid chunk size
			decomposition */
			if (invalid){
				REPORT(2,"SECOND PHASE chunk splitting ...");
				target->suggestSubaddSize(chunkSize_ ,wIn_);
				lastChunkSize_ = (wIn_% chunkSize_ ==0 ? chunkSize_ :wIn_% chunkSize_);
				
				/* the index of the last chunk pair */
				chunkIndex = (wIn_% chunkSize_ ==0 ? ( wIn_ / chunkSize_) - 1 :  (wIn_-lastChunkSize_) / chunkSize_ ); 								
				for (int i=0; i < chunkIndex; i++)
					cSize[i] = chunkSize_;
				/* last chunk is handled separately  */
				cSize[chunkIndex] = lastChunkSize_;
			}
			
			/* VERIFICATION PHASE: check if decomposition is correct */		
			REPORT(2, "found " << chunkIndex + 1  << " chunks ");
			nbOfChunks = chunkIndex + 1; 
			int sum = 0;
			ostringstream chunks;
			for (int i=chunkIndex; i>=0; i--){
				chunks << cSize[i] << " ";
				sum+=cSize[i];
			}
			chunks << endl;
			REPORT(2, "Chunks are: " << chunks.str());
			REPORT(2, "The chunk size sum is " << sum << " and initial width was " << wIn_);
			if (sum != wIn_){
				cerr << "ERROR: check the algo" << endl; /*should never get here ... */
				exit(0);
			}
			
			//=================================================
			//split the inputs ( this should be reusable )
			vhdl << tab << "--split the inputs into chunks of bits depending on the frequency" << endl;
			for (int i=0;i<2;i++){
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
					if (j==0)
						vhdl << tab << declare (name.str(),cSize[j]+1) << " <=  \"0\" & X"<<i<<range(high-1,low)<<";"<<endl;
					else
						vhdl << tab << declare (name.str(),cSize[j]) << " <= X"<<i<<range(high-1,low)<<";"<<endl;
				}
			}	
			vhdl << tab << declare("scIn",1) << " <= Cin;"<<endl;
			
			int l=1;
			for (int j=0; j<nbOfChunks; j++){
				if (j>0){ //for all chunks greater than zero we perform this comparissons
					vhdl<<tab<<declare(join("sX",j,"_0_l",l,"_Zero"))<< " <= '1' when "<< join("sX",j,"_0_l",l-1)<< " > not("<<join("sX",j,"_1_l",l-1)<<") else '0';"<<endl;
					vhdl<<tab<<declare(join("sX",j,"_0_l",l,"_One"))<< "  <= '1' when "<< join("sX",j,"_0_l",l-1)<< " >= not("<<join("sX",j,"_1_l",l-1)<<") else '0';"<<endl;
				}else{
					//for the zero chunk we directly perform the addition
					vhdl<<tab<< "-- the carry resulting from the addition of the chunk + Cin is obtained directly" << endl;
					vhdl<<tab<<declare(join("sX",j,"_0_l",l,"_Cin"),cSize[j]+1)<< "  <= " << join("sX",j,"_0_l",l-1)<<" + "<<join("sX",j,"_1_l",l-1)<<" + scIn;"<<endl;
				}
			}
			
			vhdl << tab <<"--form the two carry string"<<endl;
			vhdl << tab << declare("carryStringZero",nbOfChunks-2) << " <= "; 
			for (int i=nbOfChunks-3; i>=0; i--) {
				vhdl << "sX"<<i+1<<"_0_l"<<l<<"_Zero" << (i>0?" & ":";") ;
			} vhdl << endl;
			
			vhdl << tab << declare("carryStringOne",  nbOfChunks-2) << "  <= "; 
			for (int i=nbOfChunks-3; i>=0; i--) {
				vhdl << "sX"<<i+1<<"_0_l"<<l<<"_One" << (i>0?" & ":";");
			} vhdl << endl;
			
			//				nextCycle();/////////////////////
			
			vhdl << tab << "--perform the short carry additions" << endl;
			ostringstream unameCin;
			unameCin  << "sX"<<0<<"_0_l"<<l<<"_Cin";
			vhdl << tab << declare("rawCarrySum",nbOfChunks-2) << " <= carryStringOne + carryStringZero + " << unameCin.str() <<of(cSize[0])<< ";" << endl;
			
			
			
			//				vhdl << tab << declare("manipulatedSum",nbOfChunks-2) << "<= carryStringOne AND ( (not(rawCarrySum) AND not(carryStringZero)) OR carryStringZero);" << endl; //strike of genious
			vhdl << tab << declare("manipulatedSum",nbOfChunks-2) << "<= (not(rawCarrySum) AND carryStringOne) OR carryStringZero;" << endl; //strike of genious
			
			//				if (invalid)
			//					nextCycle();/////////////////////
			
			vhdl << tab <<"--get the final pipe results"<<endl;
			for ( int i=0; i<nbOfChunks; i++){
				ostringstream unameZero, unameOne, unameCin;
				unameZero << "sX"<<i<<"_0_l"<<l<<"_Zero";
				unameOne  << "sX"<<i<<"_0_l"<<l<<"_One";
				unameCin  << "sX"<<0<<"_0_l"<<l<<"_Cin";
				if (i==0) 
					vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << unameCin.str()<< range(cSize[0]-1,0) <<  ";" << endl;
				else {
					if (i==1){ 
							vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << join("sX",i,"_0_l0")<<range(cSize[i]-1,0) 
						<< " + " << join("sX",i,"_1_l0")<<range(cSize[i]-1,0) 
						<< " + sX0_0_l1_Cin"<<of(cSize[i-1]) <<";"<<endl;
					}else{      
						vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << join("sX",i,"_0_l0")<<range(cSize[i]-1,0) 
						<< " + " << join("sX",i,"_1_l0")<<range(cSize[i]-1,0) 
						<< " + manipulatedSum"<<of(i-2) <<";"<<endl;
					}
				}
			}
			
			vhdl << tab << "R <= ";
			int k=0;
			for (int i=nbOfChunks-1; i>=0; i--){
				vhdl << use(join("res",i));
				if (i > 0) vhdl << " & ";
				k++;
			}
			vhdl << ";" <<endl;
			
		}else
			vhdl << tab << " R <= X0 + X1 + Cin;"<<endl;
	}

	LongIntAdderCmpCmpAdd::~LongIntAdderCmpCmpAdd() {
	}


	void LongIntAdderCmpCmpAdd::emulate(TestCase* tc)
	{
		mpz_class svX[2];
		for (int i=0; i<2; i++){
			ostringstream iName;
			iName << "X"<<i;
			svX[i] = tc->getInputValue(iName.str());
		}
		mpz_class svC =  tc->getInputValue("Cin");

		mpz_class svR = svX[0] + svC;
		mpz_clrbit(svR.get_mpz_t(),wIn_); 
		for (int i=1; i<2; i++){
			svR = svR + svX[i];
			mpz_clrbit(svR.get_mpz_t(),wIn_); 
		}

		tc->addExpectedOutput("R", svR);

	}


}
