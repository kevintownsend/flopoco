/*
  An integer compressor tree for FloPoCo
  
  Authors : Bogdan Pasca

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

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


namespace flopoco{


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
			addInput (name.str() , wIn_, true);
		}

		addOutput("R"  , wIn_, 1, true);

		if (verbose){
			cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
////			cout <<"delay for Y is   "<< inputDelays["Y"]<<endl;
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
			if ((nbOfInputs == 1) || (wIn == 1)){
				vhdl << tab << "R <= ";
				for (int i=N-1; i>=0; i--){
					if (i>0)
						vhdl << "X"<<i<< " + ";
					else
						vhdl << "X"<<i<< ";"<<endl;
				} 
				processing = false;
			}else if (nbOfInputs == 2){
				name.str("");
				name << "level_" << treeLevel-1 << "_sum_";			
				vhdl << endl;
				inPortMap(finalAdder,"X",use(join(name.str(),0)));
				inPortMap(finalAdder,"Y",use(join(name.str(),1)));
				inPortMapCst(finalAdder,"Cin","'0'");
				outPortMap(finalAdder,"R","myR");
				vhdl << instance(finalAdder,"FinalAdder_CompressorTree") << endl;
				if (verbose) cerr << tab << "Final Adder Instantiation " << endl;
			
				syncCycleFromSignal("myR");
			
				vhdl << tab << "R <= myR;" << endl;
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
	
				if (verbose){
					cerr << "BACKTRACKING FINISHED" << endl;
					cerr << endl;
					cerr << " Solution = ";
				
				for (int i=1; i <=lutSize; i++)
					clog << bestSol[i] << ", ";
				clog << " having score " << bestSol[0] << endl;
				}
			
				int currentlyMapped = 0;
				int currentOutput = 0;
				int currentCompressor = 0;
				for (int i=lutSize; i>=1; i--){
					if (verbose) cerr << "mapping compressors " << i << " to " << intlog2(i) << endl;
					int sumSize = intlog2(i);			
					for (int h=1; h<=bestSol[i]; h++){
						if (verbose) cerr << tab << "number " << h << endl;
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
		if (verbose){
			cerr << endl << "solution !!! :"; 
			for (int i=1; i<=n; i++)
				cout << sol[i] << ", ";
		}

		int val=0;
		for (int i=1; i<=n; i++)
			val+=sol[i]*coef[i];
		
		if (verbose) cerr << " score = " << val << " coeficinets=";
	
		if (val < bestSol[0]){
			for (int i=1; i<=n; i++)
				bestSol[i] = sol[i];
			bestSol[0]=val;
		}
	
		if (verbose)
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
		int val = 0;
		for (int i=1; i<= (k>n?n:k); i++)  
			val+=sol[i]*coef[i];
	
		if ((val <= sum ) && (k<=n+1))
			return true;
		else
			return false;
	}

	void IntCompressorTree::bt(int k, int n, int sum, int* sol, int* coef, int targetSum, int * bestSol){
	
		if (solution(k, n, targetSum, sol, coef))
			printSolution(n, sol, coef, bestSol);
		else{
			sol[k]=-1;
			while (successor(k, sum, sol)) { //TODO optimize 
				if (valid(k,sum, n, sol, coef)){ 
					bt(k+1, n, sum, sol, coef, targetSum, bestSol);
				}
			}
		}
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

	void IntCompressorTree::outputVHDL(std::ostream& o, std::string name) {
		licence(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.std_logic_unsigned.all;" << endl;
		o << "library work;" << endl;
		outputVHDLEntity(o);
		newArchitecture(o,name);
		o << buildVHDLComponentDeclarations();	
		o << buildVHDLSignalDeclarations();
		beginArchitecture(o);		
		o<<buildVHDLRegisters();
		o << vhdl.str();
		endArchitecture(o);
	}


}
