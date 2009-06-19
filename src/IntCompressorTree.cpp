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

	// Set up the IO signals
	for (int i=0; i<N; i++){
		name.str(""); //init a ostringstream variable
		name << "X"<<i; 
		addInput (name.str() , wIn_);
	}

	//addInput ("Cin", 1  );
	addOutput("R"  , wIn_);

	if (verbose){
		cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
		cout <<"delay for Y is   "<< inputDelays["Y"]<<endl;
	}

	int lutSize = target->lutInputs();
	int nbOfCompressors;
	bool processing = true;
	int nbOfInputs = N_;
	int treeLevel = 1;
	int inputsLastCompressor;
	
	IntAdder *finalAdder = new IntAdder(target, wIn_);
	oplist.push_back(finalAdder);	

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
			nbOfCompressors = ceil( double(nbOfInputs) / double(lutSize) );
			if ((nbOfCompressors * lutSize) == 	nbOfInputs)
				inputsLastCompressor = lutSize;
			else
				inputsLastCompressor = nbOfInputs - (nbOfCompressors-1) * lutSize;

				for (int i=0; i < nbOfCompressors ; i++){ //for each compressor stages in the level 
					//form the summs 
					ostringstream name;
					int sumSize = intlog2( (i==nbOfCompressors-1? inputsLastCompressor : lutSize));					

					for (int j=0; j < wIn_; j++){ //do the compressor computation
						name.str("");
						name << "level_" << treeLevel << "_compressor_"<<i<< "_column_" << j ;
						vhdl << tab << declare(name.str(), sumSize, true) << " <= ";
						if ( i < nbOfCompressors -1){
							for (int k=lutSize*i; k< lutSize*(i+1); k++) {
								name.str("");
								name << "level_" << treeLevel-1 << "_sum_"<<k;
								vhdl << "("	<< zeroGenerator(sumSize-1,0) << " & " << use(name.str())<<of(j)<<")";
								if (k < lutSize*(i+1)-1)
									vhdl << " + ";
							}
							vhdl << ";" << endl;
						}else{
							for (int k=lutSize*i; k< lutSize*i + inputsLastCompressor; k++) {
								name.str("");
								name << "level_" << treeLevel-1 << "_sum_"<<k;
								vhdl << "("	<< zeroGenerator(sumSize-1,0) << " & " << use(name.str())<<of(j)<<")";
								if (k < lutSize*i + inputsLastCompressor-1)
									vhdl << " + ";
							}
							vhdl << ";" << endl;
						}
					}
					//form the summs for the current compressor
						for (int m=0; m < sumSize; m++){
							name.str("");
							name << "level_" << treeLevel << "_sum_"<< i*intlog2(lutSize) + m;
							vhdl << tab << declare (name.str(), wIn_, true) << " <= ";
							for (int k=wIn_-1-m; k >= 0; k--)	{							
								name.str("");//emptyName
								name << "level_" << treeLevel << "_compressor_"<<i<< "_column_" << k ;
								vhdl << use(name.str())<<of(m);
								if (k!=0)
									vhdl << " & ";
							}
							if (m>0)
								vhdl << " & " <<  zeroGenerator(m,0);
							vhdl << ";" << endl;
						}
					
				}
				
			
			cout << tab << "Level " << treeLevel << endl;
			cout << tab << tab << "Inputs/level: " << nbOfInputs << endl;
			cout << tab << tab << "Compressors: " << nbOfCompressors << endl;
			cout << tab << tab << "InputsLastCompressor: " << inputsLastCompressor << endl;
			
			nbOfInputs = (nbOfCompressors-1)*intlog2(lutSize) + intlog2(inputsLastCompressor);
			treeLevel++;
			
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


