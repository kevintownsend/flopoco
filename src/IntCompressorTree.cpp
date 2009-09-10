/*
 * An integer compressor tree for FloPoCo
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
#include "IntCompressorTree.hpp"

using namespace std;
extern vector<Operator*> oplist;

IntCompressorTree::IntCompressorTree(Target* target, int wIn, int N, map<string, double> inputDelays):
Operator(target), wIn_(wIn), N_(N), inputDelays_(inputDelays) 
{
	ostringstream name;
	name << "IntCompressorTree_" << wIn_<<"_"<<N_;
	setName(name.str());
	setCopyrightString("Bogdan Pasca (2009)");

	// Set up the IO signals
	for (int i=0; i<N; i++){
		name.str(""); //init a ostringstream variable
		name << "X"<<i; 
		addInput (name.str() , wIn_);
	}

	addOutput("R"  , wIn_);

	if (verbose){
		cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
		cout <<"delay for Y is   "<< inputDelays["Y"]<<endl;
	}

	int lutSize = target->lutInputs();
	bool processing = true;
	int nbOfInputs = N_;
	int treeLevel = 1;

	IntAdder *finalAdder = new IntAdder(target, wIn_);
	oplist.push_back(finalAdder);	

	//for homogeneous signal names
	for (int j=0; j<nbOfInputs;j++){
		name.str("");
		name << "level_" << treeLevel-1 << "_sum_"<<j;
		vhdl << tab << declare(name.str(),wIn_, true) << " <= " << use(join("X",j)) << ";" << endl;
	}
	
	while (processing){
		if (nbOfInputs == 2){
			name.str("");
			name << "level_" << treeLevel-1 << "_sum_";			
			vhdl << endl;
			inPortMap(finalAdder,"X",use(join(name.str(),0)));
			inPortMap(finalAdder,"Y",use(join(name.str(),1)));
			inPortMapCst(finalAdder,"Cin","'0'");
			outPortMap(finalAdder,"R","myR");
			vhdl << instance(finalAdder,"FinalAdder_CompressorTree") << endl;
			cout << tab << "Final Adder Instantiation " << endl;
			
			syncCycleFromSignal("myR");
			
			vhdl << tab << "R <= " << use("myR") << ";" << endl;
			processing = false;
		}else{
			int a[16];    //possibly in the future LUTS may have up to 16 inputs
			int sol[16];
			int bestSol[16];
	
			for (int i=0;i <= lutSize; i++){
				sol[i]=0;
				bestSol[i]=0;
			}
			bestSol[0]=999;
	
			for (int i=lutSize; i>=1 ; i--){
				a[i]= intlog2(double(i));
			} 
			a[0]=0;
	
			if (verbose){
				for (int i=1; i<= lutSize; i++)
					cerr << a[i] << ", "; 
				cerr << endl;
			}

			bt(1, lutSize, nbOfInputs, sol, a, nbOfInputs, bestSol);
	
			clog << "BACKTRACKING FINISHED" << endl;
			clog << endl;
			clog << " Solution = ";
			for (int i=1; i <=lutSize; i++)
				clog << bestSol[i] << ", ";
			clog << " having score " << bestSol[0] << endl;
			
			int currentlyMapped = 0;
			int currentOutput = 0;
			int currentCompressor = 0;
			for (int i=lutSize; i>=1; i--){
				cout << "mapping compressors " << i << " to " << intlog2(i) << endl;
				int sumSize = intlog2(i);			
				for (int h=1; h<=bestSol[i]; h++){
					cout << tab << "number " << h << endl;
					if (i>2){
						for (int j=0; j < wIn_; j++){ //do the compressor computation
							name.str("");
							name << "level_" << treeLevel << "_compressor_"<<currentCompressor<< "_column_" << j ;
							vhdl << tab << declare(name.str(), sumSize, true) << " <= ";
							for (int k=currentlyMapped; k<currentlyMapped+i; k++) {
								name.str("");
								name << "level_" << treeLevel-1 << "_sum_"<<k;
								vhdl << "("	<< zg(sumSize-1,0) << " & " << use(name.str())<<of(j)<<")";
								if (k < currentlyMapped+i-1)
									vhdl << " + ";
								}
								vhdl << ";" << endl;
							}
						currentCompressor++;
					}
					if (i == 2){
						name.str("");
						name << "level_" << treeLevel << "_compressor_"<<currentCompressor<< "_column_" << 0 ;
						vhdl << tab << declare(name.str(), wIn_, true) << " <= ";
						name.str("");
						name << "level_" << treeLevel-1 << "_sum_"<<currentlyMapped;
						vhdl << use(name.str()) << ";" << endl;

						name.str("");
						name << "level_" << treeLevel << "_compressor_"<<currentCompressor<< "_column_" << 1 ;
						vhdl << tab << declare(name.str(), wIn_, true) << " <= ";
						name.str("");
						name << "level_" << treeLevel-1 << "_sum_"<<currentlyMapped+1;
						vhdl << use(name.str()) << ";" << endl;
						currentCompressor++;
					}
					if (i == 1){
						name.str("");
						name << "level_" << treeLevel << "_compressor_"<<currentCompressor<< "_column_" << 0 ;
						vhdl << tab << declare(name.str(), wIn, true) << " <= ";
						name.str("");
						name << "level_" << treeLevel-1 << "_sum_"<<currentlyMapped;
						vhdl << use(name.str()) << ";" << endl;
						currentCompressor++;
					}
				
				
					//form the summs for the current compressor
					if (i>2){
						for (int m=0; m < sumSize; m++){
							name.str("");
							name << "level_" << treeLevel << "_sum_"<< currentOutput + m;
							vhdl << tab << declare (name.str(), wIn_, true) << " <= ";
							for (int k=wIn_-1-m; k >= 0; k--)	{							
								name.str("");//emptyName
								name << "level_" << treeLevel << "_compressor_"<<currentCompressor-1<< "_column_" << k ;
								vhdl << use(name.str())<<of(m);
								if (k!=0)
									vhdl << " & ";
							}
							if (m>0)
								vhdl << " & " <<  zg(m,0);
							vhdl << ";" << endl;
						}
					}else{
						if (i==2){
							name.str("");
							name << "level_" << treeLevel << "_sum_"<< currentOutput;
							vhdl << tab << declare (name.str(), wIn_, true) << " <= ";
							name.str("");//emptyName
							name << "level_" << treeLevel << "_compressor_"<<currentCompressor-1<< "_column_" << 0 ;
							vhdl << use(name.str())<<";"<<endl;
							
							name.str("");
							name << "level_" << treeLevel << "_sum_"<< currentOutput+1;
							vhdl << tab << declare (name.str(), wIn_, true) << " <= ";
							name.str("");//emptyName
							name << "level_" << treeLevel << "_compressor_"<<currentCompressor-1<< "_column_" << 1 ;
							vhdl << use(name.str())<<";"<<endl;
						}else
							if (i==1){
								name.str("");
								name << "level_" << treeLevel << "_sum_"<< currentOutput;
								vhdl << tab << declare (name.str(), wIn_, true) << " <= ";
								name.str("");//emptyName
								name << "level_" << treeLevel << "_compressor_"<<currentCompressor-1<< "_column_" << 0 ;
								vhdl << use(name.str())<<";"<<endl;
							}
					
					}
					
					currentlyMapped += i;
					currentOutput += intlog2(i);
				}
			}
			//form the sums
			
			nbOfInputs = bestSol[0];
			treeLevel++;
		}
	}		
}		


bool IntCompressorTree::solution(int k, int n, int targetSum, int * sol, int * coef){
	int val=0;
	for (int i=1; i<=n; i++)
		val+=sol[i]*i;
	if ((k==n+1) && ( val==targetSum))
		return true;
	else
		return false;
}

void IntCompressorTree::printSolution(int n, int * sol, int * coef, int *bestSol){
	cout << endl << "solution !!! :"; 
	for (int i=1; i<=n; i++)
		cout << sol[i] << ", ";

	int val=0;
	for (int i=1; i<=n; i++)
		val+=sol[i]*coef[i];
		
	cout << " score = " << val << " coeficinets=";
	
	if (val < bestSol[0]){
		for (int i=1; i<=n; i++)
			bestSol[i] = sol[i];
		bestSol[0]=val;
	}
	
	for (int i=1; i<=n ; i ++)
		cout << coef[i] << "::";	

	
}

bool IntCompressorTree::successor( int k, int sum, int * sol){
	if (sol[k] < sum) {
		sol[k]++;
		return true;
	}else
		return false;
}

bool IntCompressorTree::valid(int k, int sum, int n, int * sol, int * coef){
//	cout << "check valid " << endl;
	int val = 0;
	for (int i=1; i<= (k>n?n:k); i++)  
		val+=sol[i]*coef[i];
	
//	cout << "end check" << endl;	
	if ((val <= sum ) && (k<=n+1))
		return true;
	else
		return false;
}

void IntCompressorTree::bt(int k, int n, int sum, int* sol, int* coef, int targetSum, int * bestSol){
//	cout<<endl << "k=" << k << " === ";
//	int val=0;
//	for (int i=1; i<= (k>n?n:k); i++){
//		val+= sol[i]* coef[i];
//		cout << sol [ i ] << " : ";
//	}
//	cout << "current_sum=" << val << ", target_sum="<<targetSum;
	
	if (solution(k, n, targetSum, sol, coef))
		printSolution(n, sol, coef, bestSol);
	else{
		sol[k]=-1;
			while (successor(k, sum, sol)) { //to make more efficinet
//				cout << "doing while with k = "<<k << " and sol[k]"<<sol[k] << endl;
				if (valid(k,sum, n, sol, coef)){ 
//					cout << tab << " go to k+1" << endl;
					bt(k+1, n, sum, sol, coef, targetSum, bestSol);
				}
				
			}
		
	}
//	cout << "return from recursion: k="<<k<<endl;
}



IntCompressorTree::~IntCompressorTree() {
}


void IntCompressorTree::emulate(TestCase* tc)
{
	mpz_class svX;
	mpz_class svR = 0;

	for (int i=0; i<N_; i++){
	ostringstream iName;
		iName << "X"<<i;
		svX = tc->getInputValue(iName.str());
		svR = svR + svX;
		mpz_clrbit(svR.get_mpz_t(),wIn_); 
	}
	tc->addExpectedOutput("R", svR);

}


