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
#include "LongIntAdder.hpp"

using namespace std;


namespace flopoco{

	LongIntAdder::LongIntAdder(Target* target, int wIn, map<string, double> inputDelays):
		Operator(target), wIn_(wIn), inputDelays_(inputDelays) 
	{
		ostringstream name;
		name << "LongIntAdder_" << wIn_;
		setName(name.str());

		// Set up the IO signals
		for (int i=0; i<2; i++){
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
			objectivePeriod = 1/ (1.45*target->frequency());
		
			if (verbose)
				cout << "Objective period is "<< objectivePeriod<<" at an objective frequency of "<<target->frequency() << endl;

			if (objectivePeriod<maxInputDelay){
				//It is the responsability of the previous components to not have a delay larger than the period
				cout << "Warning, the combinatorial delay at the input of "<<this->getName()<<" is above limit"<<endl;
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
						target->setFrequency(target->frequency() * 1.45);
						target->suggestSubaddSize(chunkSize_,newWIn);
						target->setFrequency(target->frequency() / 1.45);
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
			//split the inputs ( this should be reusable )
			vhdl << tab << "--split the inputs into chunks of bits depending on the frequency" << endl;
			for (int i=0;i<2;i++)
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
			vhdl << tab << declare("scIn",1) << " <= " << "Cin;"<<endl;
			if (true){
				int l=1;
				for (int j=0; j<nbOfChunks; j++){
					ostringstream dnameZero, dnameOne, uname1, uname2, dnameCin;
					dnameZero << "sX"<<j<<"_0_l"<<l<<"_Zero";
					dnameOne  << "sX"<<j<<"_0_l"<<l<<"_One";
					dnameCin  << "sX"<<j<<"_0_l"<<l<<"_Cin";
				
					uname1 << "sX"<<j<<"_0_l"<<l-1;
					uname2 << "sX"<<j<<"_1_l"<<l-1;
					if (j>0){ //for all chunks greater than zero we perform this additions
						vhdl << tab << declare(dnameZero.str(),cSize[j]+1) << " <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<");"<<endl;
						if (j<nbOfChunks-1)
						vhdl << tab << declare(dnameOne.str(),cSize[j]+1) << "  <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + '1';"<<endl;
					}else{
						vhdl << tab << "-- the carry resulting from the addition of the chunk + Cin is obtained directly" << endl;
						vhdl << tab << declare(dnameCin.str(),cSize[j]+1) << "  <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + "<<use("scIn")<<";"<<endl;
					}
				}
			
			
				vhdl << tab <<"--form the two carry string"<<endl;
				vhdl << tab << declare("carryStringZero",2*(nbOfChunks-2)) << " <= "; 
				for (int i=2*(nbOfChunks-2)-1; i>=0; i--) {
					ostringstream dnameZero;
					if (i%2==0){
						dnameZero << "sX"<<(i/2)+1<<"_0_l"<<l<<"_Zero";
						vhdl << " " << use(dnameZero.str()) <<"(" << cSize[(i/2)+1] << ")";
					}else
						vhdl << " \"0\" ";
					if (i>0) vhdl << " & ";
					else     vhdl << " ; ";
				} vhdl << endl;
				vhdl << tab << declare("carryStringOne",2*(nbOfChunks-2)) << "  <= "; 
				for (int i=2*(nbOfChunks-2)-1; i>=0; i--) {
					ostringstream dnameOne;
					if (i%2==0){
						dnameOne << "sX"<<(i/2)+1<<"_0_l"<<l<<"_One";
						vhdl << " " << use(dnameOne.str()) <<"(" << cSize[(i/2)+1] << ") ";
					}else
						vhdl << " \"1\" ";
					if (i>0) vhdl << " & ";
					else     vhdl << " ; ";
				} vhdl << endl;

				nextCycle();/////////////////////

				vhdl << tab << "--perform the short carry additions" << endl;
				ostringstream unameCin;
				unameCin  << "sX"<<0<<"_0_l"<<l<<"_Cin";
				vhdl << tab << declare("rawCarrySum",2*(nbOfChunks-2)) << " <= " << use("carryStringOne") << " + " << use("carryStringZero") << " + " << use(unameCin.str()) << "(" << cSize[0] << ")" << " ;" << endl;

				//			nextCycle();/////////////////////
				vhdl << tab <<"--get the final pipe results"<<endl;
				for ( int i=0; i<nbOfChunks; i++){
					ostringstream unameZero, unameOne, unameCin;
					unameZero << "sX"<<i<<"_0_l"<<l<<"_Zero";
					unameOne  << "sX"<<i<<"_0_l"<<l<<"_One";
					unameCin  << "sX"<<0<<"_0_l"<<l<<"_Cin";
					if (i==0) vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameCin.str())<< range(cSize[0]-1,0) <<  ";" << endl;
					else {
						if (i==1) vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " + " << use(unameCin.str()) << "(" << cSize[0] << ")" << ";"<<endl;
						else      vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " + not(" << use("rawCarrySum")<<"("<<2*(i-2)+1<<"));"<<endl;
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




			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			if (false){
				//the first version
				////////////////////////////////////////////////
				//perform addition round; there is only one additon per round 
				int l=1;
				for (int j=0; j<nbOfChunks; j++){
					ostringstream dnameZero, dnameOne, uname1, uname2;
					dnameZero << "sX"<<j<<"_0_l"<<l<<"_Zero";
					dnameOne  << "sX"<<j<<"_0_l"<<l<<"_One";
					uname1 << "sX"<<j<<"_0_l"<<l-1;
					uname2 << "sX"<<j<<"_1_l"<<l-1;
					vhdl << tab << declare(dnameZero.str(),cSize[j]+1) << " <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<");"<<endl;
					vhdl << tab << declare(dnameOne.str(),cSize[j]+1) << " <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + '1';"<<endl;
				}

				nextCycle(); //////////////////////////////////////////////////////

				vhdl << tab << declare("carryStringZero",nbOfChunks) << " <= "; 
				for (int i=nbOfChunks-1; i>=0; i--) {
					ostringstream dnameZero;
					dnameZero << "sX"<<i<<"_0_l"<<l<<"_Zero";
					vhdl << " " << use(dnameZero.str()) <<"(" << cSize[i] << ")";
					if (i>0) vhdl << " & ";
					else     vhdl << " ; ";
				} vhdl << endl;
				vhdl << tab << declare("carryStringOne",nbOfChunks) << " <= "; 
				for (int i=nbOfChunks-1; i>=0; i--) {
					ostringstream dnameOne;
					dnameOne << "sX"<<i<<"_0_l"<<l<<"_One";
					vhdl << " " << use(dnameOne.str()) <<"(" << cSize[i] << ")";
					if (i>0) vhdl << " & ";
					else     vhdl << " ; ";
				} vhdl << endl;
				for (int i=0; i<nbOfChunks; i++){
					ostringstream unameZero, unameOne;;
					unameZero << "sX"<<i<<"_0_l"<<l<<"_Zero";
					unameOne << "sX"<<i<<"_0_l"<<l<<"_One";
					vhdl << tab << declare(join("s",i),2) << " <= " << " (\"0\" & "<<use("carryStringZero") << range(i,i) << ") + " 
						  << " (\"0\" & "<<use("carryStringOne")  << range(i,i) << ") + ";
					if (i==0)
						vhdl << use("scIn")<<";"<<endl;
					else
						vhdl << use(join("s",i-1))<<"("<<1<<");"<<endl;
				}
				vhdl << tab << declare("carrySum",nbOfChunks) << " <= ";
				for (int i=nbOfChunks-1;i>=0;i--){
					vhdl << use(join("s",i))<< "("<< (i==nbOfChunks-1?0:1)<<")" << (i==0?";":" & ");
				} vhdl << endl;

				nextCycle(); //////////////////////////////////////////////////////
		
				//get the final pipe results;
				for ( int i=0; i< nbOfChunks; i++){
					ostringstream unameZero, unameOne;
					unameZero << "sX"<<i<<"_0_l"<<l<<"_Zero";
					unameOne  << "sX"<<i<<"_0_l"<<l<<"_One";
					//			if (i==0)vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " when " << use("scIn") <<"='0' else " << use(unameOne.str())<<range(cSize[i]-1,0) << ";" << endl;
					//			else vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " when " << use("carrySum")<<"("<<i-1<<")"
					//			                                                    <<"='0' else " << use(unameOne.str())<<range(cSize[i]-1,0) << ";" << endl;
					if (i==0) vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " + " << use("scIn") << ";" << endl;
					else vhdl << tab << declare(join("res",i),cSize[i]) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " + " << use("carrySum")<<"("<<i-1<<");"<<endl;
				}
		
		
				vhdl << tab << "R <= ";
				int k=0;
				for (int i=nbOfChunks-1; i>=0; i--){
					vhdl << use(join("res",i));
					if (i > 0) vhdl << " & ";
					k++;
				}
				vhdl << ";" <<endl;
			}
		}else
			vhdl << tab << " R <= X0 + X1 + Cin;"<<endl;
	}

	LongIntAdder::~LongIntAdder() {
	}


	void LongIntAdder::emulate(TestCase* tc)
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
