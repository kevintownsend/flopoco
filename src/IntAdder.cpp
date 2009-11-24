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
#include "IntAdder.hpp"

using namespace std;

namespace flopoco{
extern vector<Operator*> oplist;

	IntAdder::IntAdder(Target* target, int wIn, map<string, double> inputDelays, int optimizeType):
		Operator(target), wIn_(wIn), inputDelays_(inputDelays)
	{
		srcFileName="IntAdder";
		ostringstream name;
		setCopyrightString("Bogdan Pasca, Florent de Dinechin (2008-2009)");		
		name << "IntAdder_" << wIn_<<"_f"<<target->frequencyMHz();
		setName(name.str());

		// Set up the IO signals
		addInput ("X"  , wIn_, true);
		addInput ("Y"  , wIn_, true);
		addInput ("Cin", 1 );
		addOutput("R"  , wIn_, 1 , true);

		REPORT( DETAILED, "Printing input delays ... " << printInputDelays(inputDelays));

		if (isSequential()){
		
			int selectedImplementation = implementationSelector( target, wIn, inputDelays, optimizeType);
			REPORT(INFO, "The implementation is: " << (selectedImplementation==0? "Classical": (selectedImplementation==1?"Alternative":"Short-Lantency")));
		
			/* shared variables */
			objectivePeriod = 1 / target->frequency();	

			//**********************************************************************
			//**********************************************************************
			// Classical implementation of pipelined addition
			if ( selectedImplementation == CLASSICAL ){

				maxInputDelay = getMaxInputDelays (inputDelays);
				REPORT( DETAILED, "The maximum input delay is: " << maxInputDelay );
		
				if ( maxInputDelay > objectivePeriod ){
					REPORT( INFO, "WARNING: combinatorial delay at the input of " << this->getName() << " is > objective period ");
					maxInputDelay = objectivePeriod;
				}

				if ( maxInputDelay == 0 ){
					/* the non-slack version */			
					REPORT(DETAILED, "Building architecture for classical version: no slack");	
					updateParameters ( target, alpha, beta, k);
					REPORT(DETAILED, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
				}else{
					if (classicalSlackVersion == 0){
						REPORT(DETAILED, "Building architecture for classical version with slack version 0");
						updateParameters ( target, inputDelays, alpha, beta, gamma, k);
						REPORT(DETAILED, "alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" k="<<k);
					}else{
						nextCycle(); ///////////////////////////////////////////////
						REPORT(DETAILED, "Building architecture for classical version with slack version 1");
						updateParameters (target, alpha, beta, k);
						REPORT(DETAILED, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
					}
				}
				
				//the sizes of the chunks
				cSize = new int[k+1];
				if ( k > 1 ){
					if (selectedDesign == 1 )
						cSize[0] = gamma;
					else
						cSize[0] = alpha;
					
					for (int i=1; i<k-1; i++)
						cSize[i] = alpha;
					cSize[k-1] = beta;
				}
				else{
					k = 1;
					cSize = new int[1];
					cSize[0] = wIn_;
				}

				//the indexes in the inputs of the chunks
				cIndex = new int[k];
				cIndex[0]= cSize[0];
				for (int i=1; i < k; i++)
					cIndex[i] = cIndex[i-1] + cSize[i];
		
				ostringstream verb;
				verb << "The chunk sizes[MSB-->LSB]: ";
				for (int i=k-1;i>=0;i--)
					verb << cSize[i]<<" ";
				REPORT( DETAILED, verb.str());
				verb.str("");
				verb << "The index sizes[MSB-->LSB]: ";
				for (int i=k-1;i>=0;i--)
					verb<<cIndex[i]<<" ";
				verb<<endl; 
				REPORT( DEBUG, verb.str());
				
				/* The implementation */
				for (int i=0; i < k; i++){
					vhdl << tab << declare( join("x",i), cSize[i],true) << " <= " << use("X") << range(cIndex[i]-1,(i>0?cIndex[i-1]:0)) << ";" << endl;
					vhdl << tab << declare( join("y",i), cSize[i],true) << " <= " << use("Y") << range(cIndex[i]-1,(i>0?cIndex[i-1]:0)) << ";" << endl;
				}
			
				for (int i=0; i < k; i++){
					vhdl << tab << declare(join("sum",i),cSize[i]+1,true) << " <= " << "( \"0\" & "<< use(join("x",i)) << ") + " << "( \"0\" & "<< use(join("y",i)) << ")  + ";
					if (i==0) vhdl << use("Cin");
					else      vhdl << use(join("sum",i-1))<<of(cSize[i-1]);
					vhdl << ";"	<< endl;
					if (i < k-1)
						nextCycle(); //////////////////////////////////////////
				}
	
				vhdl << tab << "R <= ";
				for (int i=k-1; i >= 1; i--){
					vhdl << use(join("sum",i))<<range(cSize[i]-1,0)<< " & ";
				}
				vhdl << use("sum0")<<range(cSize[0]-1,0)<<";"<<endl;
				/* the output is asociated with the combinatorial delay caused 
				by the most-significant bits addition */
				outDelayMap["R"] = target->adderDelay(cSize[k-1]); 
			}
			
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			else if (selectedImplementation == ALTERNATIVE) {

				maxInputDelay = getMaxInputDelays (inputDelays);
				REPORT( DETAILED, "The maximum input delay is: " << maxInputDelay );
		
				if ( maxInputDelay > objectivePeriod ){
					REPORT( INFO, "WARNING: combinatorial delay at the input of " << this->getName() << " is > objective period ");
					maxInputDelay = objectivePeriod;
				}

				if ( maxInputDelay == 0 ){
					/* the non-slack version */			
					REPORT(DETAILED, "Building architecture for alternative version: no slack");	
					updateParameters ( target, alpha, beta, k);
					REPORT(DETAILED, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
				}else{
					if (alternativeSlackVersion == 0){
						REPORT(DETAILED, "Building architecture for alternative version with slack version 0");
						updateParameters ( target, inputDelays, alpha, beta, k);
						REPORT(DETAILED, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
					}else{
						nextCycle(); ///////////////////////////////////////////////
						REPORT(DETAILED, "Building architecture for alternative version with slack version 1");
						updateParameters (target, alpha, beta, k);
						REPORT(DETAILED, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
					}
				}

				//the sizes of the chunks
				cSize = new int[k+1];
				if ( k > 1 ){
					for (int i=0; i<k-1; i++)
						cSize[i] = alpha;
					cSize[k-1] = beta;
				}
				else{
					k = 1;
					cSize = new int[1];
					cSize[0] = wIn_;
				}

				//the indexes in the inputs of the chunks
				cIndex = new int[k];
				cIndex[0]= cSize[0];
				for (int i=1; i < k; i++)
					cIndex[i] = cIndex[i-1] + cSize[i];
		
				ostringstream verb;
				verb << "The chunk sizes[MSB-->LSB]: ";
				for (int i=k-1;i>=0;i--)
					verb << cSize[i]<<" ";
				REPORT( DETAILED, verb.str());
				verb.str("");
				verb << "The index sizes[MSB-->LSB]: ";
				for (int i=k-1;i>=0;i--)
					verb<<cIndex[i]<<" ";
				verb<<endl; 
				REPORT( DEBUG, verb.str());


				////////////////////////////////////////////////////////////////////////
		
				for (int i=0; i < k; i++){
					vhdl << tab << declare (join("s_sum_l",0,"_idx",i), cSize[i]+1, true) << " <= " << "( \"0\" & " << use("X") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ") + "
						  << "( \"0\" & " << use("Y") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ")" ;
					if (i==0) vhdl << " + " << use("Cin");
					vhdl << ";" << endl;
				}
				for (int i=0; i < k; i++){
					vhdl << tab << declare (join("sum_l",0,"_idx",i), cSize[i], true) << " <= " << use(join("s_sum_l",0,"_idx",i))<<range(cSize[i]-1,0) << ";" << endl;
					vhdl << tab << declare (join("c_l",0,"_idx",i), 1, true) << " <= " << use(join("s_sum_l",0,"_idx",i))<<range(cSize[i],cSize[i]) << ";" << endl;
				}			
		
				for (int i=1; i <= k-1 ; i++){
					nextCycle(); ///////////////////////////////////////////////////////
					vhdl << tab << declare(join("sum_l",i,"_idx",i), cSize[i]+1, true) << " <= " << "( \"0\" & " << use(join("sum_l",0,"_idx",i))<< ") + " 
					                                                                             << use(join("c_l",0,"_idx",i-1))<<range(0,0) ;
					if (i>1) 
						vhdl << " + " << use(join("sum_l",i-1,"_idx",i-1))<<of(cSize[i-1]);
					vhdl<<";"<<endl;
				}
		
				vhdl << tab << "R <= ";
				for (int i=k-1; i >= 1; i--){
					vhdl << use(join("sum_l",i,"_idx",i))<<range(cSize[i]-1,0)<< " & ";
				}
				vhdl << use("sum_l0_idx0")<<range(cSize[0]-1,0)<<";"<<endl;

				outDelayMap["R"] = target->adderDelay(cSize[k-1]); 
			
				
			
			}
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			
			else if (selectedImplementation == SHORTLATENCY) {
				if (shortLatencyVersion == 0){
					tryOptimizedChunkSplittingShortLatency ( target, wIn , k );
				}else{
					updateParameters (target, alpha, beta, k);
					//the sizes of the chunks
					cSize = new int[k+1];
					if ( k > 1 ){
						for (int i=0; i<k-1; i++)
							cSize[i] = alpha;
						cSize[k-1] = beta;
					}
					else{
						k = 1;
						cSize = new int[1];
						cSize[0] = wIn_;
					}
				}

				REPORT( DEBUG, " version " << shortLatencyVersion);
				
				ostringstream verb;
				verb << "Implementing ShortLatency with chunks: ";
				for (int i=k-1; i>=0; i--)
					verb << " " << cSize[i];
			
				REPORT(DETAILED, verb.str());
			
			
				//split the inputs 
				for (int i=0;i<2;i++)
					for (int j=0; j<k; j++){
						ostringstream name;
						//the naming standard: sX j _ i _ l
						//j=the chunk index i is the input index and l is the current level
						name << "sX"<<j<<"_"<<i<<"_l"<<0;
						int low=0, high=0;
						for (int k=0;k<=j;k++)
							high+=cSize[k];
						for (int k=0;k<=j-1;k++)
							low+=cSize[k];
						vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << " \"0\" & "<<(i==0?"X":"Y")<<range(high-1,low)<<";"<<endl;
					}
				vhdl << tab << declare("scIn",1) << " <= " << "Cin;"<<endl;
			
				int l=1;
				for (int j=0; j<k; j++){
					ostringstream dnameZero, dnameOne, uname1, uname2, dnameCin;
					dnameZero << "sX"<<j<<"_0_l"<<l<<"_Zero";
					dnameOne  << "sX"<<j<<"_0_l"<<l<<"_One";
					dnameCin  << "sX"<<j<<"_0_l"<<l<<"_Cin";
				
					uname1 << "sX"<<j<<"_0_l"<<l-1;
					uname2 << "sX"<<j<<"_1_l"<<l-1;
					
#ifdef XILINX_OPTIMIZATION
	// the xst synthetsies x+y and x+y+1 slower if this optimization is not used				
					bool pipe = target->isPipelined();
					target->setNotPipelined();

					IntAdder *adder = new IntAdder(target, cSize[j]+1);
					ostringstream a;
					a.str("");

					vhdl << tab << declare( uname1.str()+"ext", cSize[j]+1 ) << " <= " <<  "\"0\" & "  <<  use(uname1.str()) << range(cSize[j]-1,0) << ";" << endl;
					vhdl << tab << declare( uname2.str()+"ext", cSize[j]+1 ) << " <= " <<  "\"0\" & "  <<  use(uname2.str()) << range(cSize[j]-1,0) << ";" << endl;
			
					if (j>0){ //for all chunks greater than zero we perform this additions
							bool found = false;
							for(unsigned k=0; k<oplist.size(); k++) {
								if  ( (oplist[k]->getName()).compare(adder->getName()) ==0 ){
									REPORT(3, "found in opList ... not adding");
									found = true;
								}
							}
						if (found == false)						{
							REPORT(3, "Not found in list, adding " << adder->getName());
							oplist.push_back(adder);
						}
				
						inPortMapCst(adder, "X", uname1.str()+"ext" );
						inPortMapCst(adder, "Y", uname2.str()+"ext" );
						inPortMapCst(adder, "Cin", "'0'");
						outPortMap(adder, "R", dnameZero.str());
						a<< "Adder" << cSize[j]+1 << "Zero" << j;
						vhdl << instance( adder, a.str() );
				
						if (j<k-1){
					
							inPortMapCst(adder, "X", uname1.str()+"ext" );
							inPortMapCst(adder, "Y", uname2.str()+"ext" );
							inPortMapCst(adder, "Cin", "'1'");
							outPortMap(adder, "R", dnameOne.str());
							a.str("");
							a<< "Adder" << cSize[j]+1 << "One" << j;
							vhdl << instance( adder, a.str() );
						}
					if (pipe)
						target->setPipelined();
					else
						target->setNotPipelined();
					}else{
						vhdl << tab << "-- the carry resulting from the addition of the chunk + Cin is obtained directly" << endl;
						vhdl << tab << declare(dnameCin.str(),cSize[j]+1) << "  <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + "<<use("scIn")<<";"<<endl;
					}
#else					
					if (j>0){ //for all chunks greater than zero we perform this additions
						vhdl << tab << declare(dnameZero.str(),cSize[j]+1) << " <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<");"<<endl;
						if (j<k-1)
						vhdl << tab << declare(dnameOne.str(),cSize[j]+1) << "  <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + '1';"<<endl;
					}else{
						vhdl << tab << "-- the carry resulting from the addition of the chunk + Cin is obtained directly" << endl;
						vhdl << tab << declare(dnameCin.str(),cSize[j]+1) << "  <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + "<<use("scIn")<<";"<<endl;
					}
#endif
				}
			
				vhdl << tab <<"--form the two carry string"<<endl;
				vhdl << tab << declare("carryStringZero",2*(k-2)) << " <= "; 
				for (int i=2*(k-2)-1; i>=0; i--) {
					ostringstream dnameZero;
					if (i%2==0){
						dnameZero << "sX"<<(i/2)+1<<"_0_l"<<l<<"_Zero";
						vhdl << " " << use(dnameZero.str()) <<"(" << cSize[(i/2)+1] << ")";
					}else
						vhdl << " \"0\" ";
					if (i>0) vhdl << " & ";
					else     vhdl << " ; ";
				} vhdl << endl;
				vhdl << tab << declare("carryStringOne",2*(k-2)) << "  <= "; 
				for (int i=2*(k-2)-1; i>=0; i--) {
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
				vhdl << tab << declare("rawCarrySum",2*(k-2)) << " <= " << use("carryStringOne") << " + " << use("carryStringZero") << " + " << use(unameCin.str()) << "(" << cSize[0] << ")" << " ;" << endl;

				if (shortLatencyVersion==1)
					nextCycle();/////////////////////
				vhdl << tab <<"--get the final pipe results"<<endl;
				for ( int i=0; i<k; i++){
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
				for (int i=k-1; i>=0; i--){
					vhdl << use(join("res",i));
					if (i > 0) vhdl << " & ";
				}
				vhdl << ";" <<endl;			
				outDelayMap["R"] = target->adderDelay(cSize[k-1]); 			
			
			}
			else {
				vhdl << tab << " R <= X + Y + Cin;" << endl;
				outDelayMap["R"] = target->adderDelay(wIn_); 
			}
		// COMBINATORIAL VERSION	
		}else{
			vhdl << tab << " R <= X + Y + Cin;" << endl;
			outDelayMap["R"] = target->adderDelay(wIn_); 
		}
	
	
	}

	IntAdder::~IntAdder() {
	}


	int IntAdder::implementationSelector(Target* target, int wIn, map<string, double> inputDelays, int optimizeType){
	
		if ((target->getID()=="Virtex4")||(target->getID()=="Virtex5")||(target->getID()=="Spartan3")){
			/* XILINX TARGETS */
			if (optimizeType == LOGIC){
				int lutCostClassical     = getLutCostClassical(target, wIn, inputDelays);
				int lutCostAlternative   = getLutCostAlternative(target, wIn, inputDelays);
				int lutCostShortLantency = getLutCostShortLatency(target, wIn, inputDelays);
				
				REPORT(DETAILED, "LUT cost for the classical     mathod " << lutCostClassical );
				REPORT(DETAILED, "LUT cost for the alternative   mathod " << lutCostAlternative );
				REPORT(DETAILED, "LUT cost for the short-latency mathod " << lutCostShortLantency	 );
				
				int bestScore = min( min( lutCostClassical, lutCostAlternative), lutCostShortLantency);
				
				if (bestScore == lutCostClassical)
					return 0;
				else if (bestScore == lutCostAlternative)
					return 1;
				else if (bestScore == lutCostShortLantency)
					return 2;
				else{
					cerr << " Error in cost evaluation " << endl;
					exit(-1);
				}
			} else if (optimizeType == REGISTER){
				int regCostClassical     = getRegCostClassical(target, wIn, inputDelays);
				int regCostAlternative   = getRegCostAlternative(target, wIn, inputDelays);
				int regCostShortLantency = getRegCostShortLatency(target, wIn, inputDelays);
				
				REPORT(DETAILED, "REG cost for the classical     mathod " << regCostClassical );
				REPORT(DETAILED, "REG cost for the alternative   mathod " << regCostAlternative );
				REPORT(DETAILED, "REG cost for the short-latency mathod " << regCostShortLantency	 );
				
				int bestScore = min( min( regCostClassical, regCostAlternative), regCostShortLantency);
				
				if (bestScore == regCostClassical)
					return 0;
				else if (bestScore == regCostAlternative)
					return 1;
				else if (bestScore == regCostShortLantency)
					return 2;
				else{
					cerr << " Error in cost evaluation " << endl;
					exit(-1);
				}
			} else if (optimizeType == SLICE){
				int sliceCostClassical     = getSliceCostClassical(target, wIn, inputDelays);
				int sliceCostAlternative   = getSliceCostAlternative(target, wIn, inputDelays);
				int sliceCostShortLantency = getSliceCostShortLatency(target, wIn, inputDelays);
				
				REPORT(DETAILED, "SLICE cost for the classical     mathod " << sliceCostClassical );
				REPORT(DETAILED, "SLICE cost for the alternative   mathod " << sliceCostAlternative );
				REPORT(DETAILED, "SLICE cost for the short-latency mathod " << sliceCostShortLantency	 );
				
				int bestScore = min( min( sliceCostClassical, sliceCostAlternative), sliceCostShortLantency);
				
				if (bestScore == sliceCostClassical)
					return CLASSICAL;
				else if (bestScore == sliceCostAlternative)
					return ALTERNATIVE;
				else if (bestScore == sliceCostShortLantency)
					return SHORTLATENCY;
				else{
					cerr << " Error in cost evaluation " << endl;
					exit(-1);
				}
			}if (optimizeType == LATENCY){
				//TODO
				return SHORTLATENCY;
			}else {
				cerr << "Error in " <<  __FILE__ << "@" << __LINE__ << ": optimization parameter must be 0, 1 or 2" <<  endl;
				exit(-1);
			}
						
		}else{
			/* ALTERA TARGETS */
		}
		
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}
	
	
	int IntAdder::getLutCostClassical(Target* target, int wIn, map<string, double> inputDelays){
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k;
			updateParameters( target, alpha, beta, k);
			
			if (k == 1) {
				return wIn;
			}else if (k == 2){
				return (alpha + beta);
			}else{
				/* more than two chunk splitting */
				return ((4*k-9)*alpha + 3*beta);
			}
		}else{
			int version0, version1;
			int alpha, beta, gamma, k;			
			updateParameters( target, inputDelays, alpha, beta, gamma, k);
			REPORT( DEBUG, "LUT, Classical, Version 0: alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" k="<<k);
			
			if (k>0)
				if (k==1){			
					version0=wIn;
				}else if (k == 2){
					version0 = alpha + gamma;
				}else{
					version0 = 4*wIn - 3*alpha - 2*gamma  - beta;
				}
			else
				version0 = 16000; //infinity
			
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "LUT, Classical, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (k==1){			
				version1=wIn;
			}else if (k == 2){
				version1 = alpha + 3*beta;
			}else{
				version1 = 4*wIn - 3*alpha - beta;
			}
			
			REPORT(DETAILED, "Version 0: " << version0);
			REPORT(DETAILED, "Version 1: " << version1);
			
			if (version0 <= version1){
				classicalSlackVersion = 0;
				REPORT( DETAILED, "Classical slack version is 0 with cost " << version0 );
				return version0;
			}else{
				classicalSlackVersion = 1;
				REPORT( DETAILED, "Classical slack version is 1 with cost " << version1 );
				return version1;
			}
		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getLutCostAlternative(Target* target, int wIn, map<string, double> inputDelays){
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k;
			updateParameters( target, alpha, beta, k);
			
			if (k == 1) {
				return wIn;
			}else if (k == 2){
				return (alpha + 2*beta);
			}else{
				/* more than two chunk splitting */
				return ((4*k-8)*alpha + 3*beta + k-3);
			}
		}else{
			int version0, version1;
			int alpha, beta, k;			
			updateParameters( target, inputDelays, alpha, beta, k);
			REPORT( DEBUG, "LUT, Alternative, Version 0: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (k>0)
				if (k==1){			
					version0=wIn;
				}else if (k == 2){
					version0 = alpha + 2*beta;
				}else{
					version0 = (4*k-8)*alpha + 3*beta + k-3;
				}
			else
				version0 = 16000; //infinity
			
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "LUT, Alternative, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (k==1){			
				version1=wIn;
			}else if (k == 2){
				version1 = alpha + 2*beta;
			}else{
				version1 = (4*k-8)*alpha + 3*beta + k-3;
			}
			
			REPORT(DETAILED, "Version 0: " << version0);
			REPORT(DETAILED, "Version 1: " << version1);
			
			if (version0 <= version1){
				alternativeSlackVersion = 0;
				REPORT( DETAILED, "Alternative slack version is 0 with cost " << version0 );
				return version0;
			}else{
				alternativeSlackVersion = 1;
				REPORT( DETAILED, "Alternative slack version is 1 with cost " << version1 );
				return version1;
			}
		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getLutCostShortLatency(Target* target, int wIn, map<string, double> inputDelays){
		bool status = tryOptimizedChunkSplittingShortLatency ( target, wIn , k );

		if ( getMaxInputDelays( inputDelays ) == 0 ){
		/* no input slack problem */
			if (status == false){
				shortLatencyVersion = 1;
				int alpha, beta, k;
				updateParameters( target, alpha, beta, k);
				if (k>=3){
					return ((4*k - 6)*alpha + 3*beta + 2*k-4);
				}else
					return 16000; //+inf			
			} else{
				shortLatencyVersion = 0;
				/* the algorithm found a good way to split inputs and save 1 pipeline level */
				if (k>=3){
					return (3*wIn - 2*cSize[0] - cSize[k-1]  + 2*(k-2));
				}else
					return 16000; //+inf
			}
		}else{
			//TODO No implementation for now, so cost is +inf
			return 16000;
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getRegCostClassical(Target* target, int wIn, map<string, double> inputDelays){
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k;
			updateParameters( target, alpha, beta, k);
			
			if (k == 1) {
				return 0;
			}else{
				return ((3*k-5)*alpha + 2*beta + k - 1);
			}
		}else{
			int version0, version1;
			int alpha, beta, gamma, k;			
			updateParameters( target, inputDelays, alpha, beta, gamma, k);
			
			if (k>0)
				if (k==1){			
					version0=0;
				}else if (k == 2){
					version0 = 2*alpha + gamma + 1;
				}else{
					version0 = 3*wIn - 2*gamma  - beta + k - 1;
				}
			else
				version0 = 16000; //infinity
			
			updateParameters( target, alpha, beta, k);

			if (k==1){			
				version1=wIn;
			}else if (k == 2){
				version1 = 3*alpha + 2*beta + 2;
			}else{
				version1 = 3*wIn - beta + k;
			}
			
			REPORT(DETAILED, "Version 0: " << version0);
			REPORT(DETAILED, "Version 1: " << version1);
			
			if (version0 <= version1){
				classicalSlackVersion = 0;
				REPORT( DETAILED, "Classical slack version is 0 with cost " << version0 );
				return version0;
			}else{
				classicalSlackVersion = 1;
				REPORT( DETAILED, "Classical slack version is 1 with cost " << version1 );
				return version1;
			}

		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getRegCostAlternative(Target* target, int wIn, map<string, double> inputDelays){
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k;
			updateParameters( target, alpha, beta, k);
			
			if (k == 1) {
				return wIn;
			}else{
				return ((2*k-3)*alpha + beta + 2*k-3);
			}
		}else{
			int version0, version1;
			int alpha, beta, k;			
			updateParameters( target, inputDelays, alpha, beta, k);
			REPORT( DEBUG, "REG, Alternative, Version 0: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (k>0)
				if (k==1){			
					version0=0;
				}else{
					version0 = (2*k-3)*alpha + beta + 2*k-3;
				}
			else
				version0 = 16000; //infinity
			
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "REG, Alternative, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (k==1){			
				version1 = 0;
			}else{
				version1 = (4*k-8)*alpha + 3*beta + 2*k-3 + 2*wIn+1;
			}
			
			REPORT(DETAILED, "Version 0: " << version0);
			REPORT(DETAILED, "Version 1: " << version1);
			
			if (version0 <= version1){
				alternativeSlackVersion = 0;
				REPORT( DETAILED, "Alternative slack version is 0 with cost " << version0 );
				return version0;
			}else{
				alternativeSlackVersion = 1;
				REPORT( DETAILED, "Alternative slack version is 1 with cost " << version1 );
				return version1;
			}
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getRegCostShortLatency(Target* target, int wIn, map<string, double> inputDelays){
		bool status = tryOptimizedChunkSplittingShortLatency ( target, wIn , k );
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			if ( status == false ){
				shortLatencyVersion = 1;
				int alpha, beta, k;
				updateParameters( target, alpha, beta, k);
				if (k>=3){
					return ((k - 1)*alpha + beta + 4*k-7);
				}else
					return 16000; 		
			}else{
				shortLatencyVersion = 0;	
				if (k>=3){
					return ((k - 1)*alpha + beta + 2*k-7);;
				}else
					return 16000; 		
			}
		}else{
			/* TODO for slack */
			return 16000;
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getSliceCostClassical(Target* target, int wIn, map<string, double> inputDelays){
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k;
			updateParameters( target, alpha, beta, k);
			
			if (k == 1) {
				return int(ceil( double(wIn) / double(2) ) );
			}else{
				return int(ceil(double((4*k-7)*alpha + 3*beta + k - 1) / double(2)));
			}
		}else{
			int version0, version1;
			int alpha, beta, gamma, k;			
			updateParameters( target, inputDelays, alpha, beta, gamma, k);
			
			if (k>0)
				if (k==1){			
					version0= int(ceil(double(wIn)/double(2)));
				}else if (k == 2){
					version0 = int(ceil(double(3*alpha + gamma + 1)/double(2)));
				}else{
					version0 = int(ceil(double(4*wIn - alpha - 2*gamma -2*beta + k - 1)/double(2)));
				}
			else
				version0 = 16000; //infinity
			
			updateParameters( target, alpha, beta, k);

			if (k==1){			
				version1= int(ceil(double(wIn)/double(2)));
			}else if (k == 2){
				version1 = int(ceil(double(3*alpha + 2*beta + 2)/double(2)));
			}else{
				version1= int(ceil(double(4*wIn - alpha - beta)/double(2)));
			}
			
			REPORT(DETAILED, "Version 0: " << version0);
			REPORT(DETAILED, "Version 1: " << version1);
			
			if (version0 <= version1){
				classicalSlackVersion = 0;
				REPORT( DETAILED, "Classical slack version is 0 with cost " << version0 );
				return version0;
			}else{
				classicalSlackVersion = 1;
				REPORT( DETAILED, "Classical slack version is 1 with cost " << version1 );
				return version1;
			}

		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getSliceCostAlternative(Target* target, int wIn, map<string, double> inputDelays){
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k;
			updateParameters( target, alpha, beta, k);
			
			if (k == 1) {
				return int(ceil( double(wIn) / double(2) ) );
			} else if (k == 2 ){
				return int(ceil(double(alpha + 2*beta + 1) / double(2)));
			}else{
				return int(ceil(double((4*k-8)*alpha + 3*beta + 2*k - 5) / double(2)));
			}
		}else{
			int version0, version1;
			int alpha, beta, k;			
			updateParameters( target, inputDelays, alpha, beta, k);
			REPORT( DEBUG, "SLICE, Alternative, Version 0: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (k>0){
				if (k==1){			
					version0=0;
				}else if (k==2){
					version0 = int(ceil(double(alpha + 2*beta + 1 )/2.));
				}else{
					version0 = int(ceil(double((4*k-8)*alpha + 3*beta + 2*k-5 )/2.));
				}
			}else
				version0 = 16000; //infinity
			
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "SLICE, Alternative, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (k==1){			
				version1 = 0;
			}else if (k==2){
				version1 = int(ceil(double(alpha + 2*beta + 1 + 2*wIn + 1)/2.));
			}else{
				version1 = int(ceil(double((4*k-8)*alpha + 3*beta + 2*k-5 + 2*wIn + 1)/2.));
			}
		
			REPORT(DETAILED, "Version 0: " << version0);
			REPORT(DETAILED, "Version 1: " << version1);
		
			if (version0 <= version1){
				alternativeSlackVersion = 0;
				REPORT( DETAILED, "Alternative slack version is 0 with cost " << version0 );
				return version0;
			}else{
				alternativeSlackVersion = 1;
				REPORT( DETAILED, "Alternative slack version is 1 with cost " << version1 );
				return version1;
			}		
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getSliceCostShortLatency(Target* target, int wIn, map<string, double> inputDelays){
		bool status = tryOptimizedChunkSplittingShortLatency ( target, wIn , k );
		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			if ( status == false ){
				shortLatencyVersion = 1;
				int alpha, beta, k;
				updateParameters( target, alpha, beta, k);
				if (k>=3){
					return int(ceil(double((4*k-6)*alpha + 3*beta + 2*k - 4) / double(2)));
				}else
					return 16000; //TODO			
			}else{
				shortLatencyVersion = 0;
				if (k>=3){
					return int(ceil(double( 3*wIn - 2*cSize[0] - cSize[k-1] + 2*(k-1) + 1 ) / double(2)));
				}else
					return 16000; //TODO			
			}
		}else{
			/* TODO for slack */
			return 16000;
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	void IntAdder::updateParameters( Target* target, int &alpha, int &beta, int &k){
		
		target->suggestSubaddSize(alpha , wIn_); /* chunk size */
		if (wIn_ == alpha ){ 
			/* addition requires one chunk */
			beta = 0;
			k    = 1;	
		}else{
			beta = ( wIn_ % alpha == 0 ? alpha : wIn_ % alpha );
			k    = ( wIn_ % alpha == 0 ? wIn_ / alpha : ceil(double(wIn_) / double(alpha)));
		}
	}
	
	void IntAdder::updateParameters( Target* target, map<string, double> inputDelays, int &alpha, int &beta, int &gamma, int &k){

		int typeOfChunks = 1;
		bool status = target->suggestSlackSubaddSize(gamma , wIn_, getMaxInputDelays(inputDelays) ); // the first chunk size
		REPORT(DEBUG, "suggestSlackSubaddSize returns gamma="<<gamma<<" with status:"<< (status?"true":"false"));
		
		if ( ! status ){
			k=-1;
			alpha=0;
			beta=0;
			gamma=0;
		}else				
			if (wIn_ - gamma > 0){ //more than 1 chunk
				target->suggestSubaddSize(alpha, wIn_-gamma); 
				if ( wIn_-gamma == alpha ) 
					typeOfChunks++; //only two types of chunks
				else
					typeOfChunks+=2; //three types of chunks
		
				if (typeOfChunks==3)
					beta = ( (wIn_-gamma) % alpha == 0 ? alpha : (wIn_-gamma) % alpha );
				else
					beta = 0;
				
				if (typeOfChunks==2)
					k = 2;
				else
					k = 2 + ( wIn_ - beta - gamma ) / alpha;
			}else{
				alpha = 0;
				beta = 0;
				k=1;
			}
	}

	void IntAdder::updateParameters( Target* target, map<string, double> inputDelays, int &alpha, int &beta, int &k){
		bool status = target->suggestSlackSubaddSize(alpha , wIn_,  getMaxInputDelays(inputDelays)); /* chunk size */
		if (!status){
			k=-1;
			alpha=0;
			beta=0;
		}else
			if (wIn_ == alpha ){ 
				/* addition requires one chunk */
				beta = 0;
				k    = 1;	
			}else{
				beta = ( wIn_ % alpha == 0 ? alpha : wIn_ % alpha );
				k    = ( wIn_ % alpha == 0 ? wIn_ / alpha : ceil(double(wIn_) / double(alpha)));
			}

	}

	bool IntAdder::tryOptimizedChunkSplittingShortLatency(Target* target, int wIn, int k){
		cSize = new int[2000];

		REPORT(DETAILED, "Trying to optimize chunk splitting ...");

		target->suggestSubaddSize(chunkSize_ ,wIn_);
		REPORT(2, "The chunkSize for first two chunks would be: " << chunkSize_ );
	
		if (2*chunkSize_ >= wIn_){
			REPORT(DETAILED, "Failed ... the input size is too small to use the short latency algorithm");
			return false;
		}
			
		cSize[0] = chunkSize_;
		cSize[1] = chunkSize_;
	
		bool finished = false; /* detect when finished the first the first
			                   phase of the chunk selection algo */
		int width = wIn - 2*chunkSize_; /* remaining size to split into chunks */
		int propagationSize = 0; /* carry addition size */
		int chunkIndex = 2; /* the index of the chunk for which the size is
			                to be determined at the current step */
		bool invalid = false; /* the result of the first phase of the algo */
	
		/* FIRST PHASE */
		while (not (finished))	 {
			propagationSize+=2;
			double delay = objectivePeriod - target->adderDelay(width)- target->adderDelay(propagationSize); //2*target->localWireDelay()  -
			if ((delay > 0) || (width < 4)) {
				cSize[chunkIndex] = width;
				finished = true;
			}else{
				int cs; 
				double slack =  target->adderDelay(propagationSize) ; //+ 2*target->localWireDelay()
				target->suggestSlackSubaddSize( cs, width, slack);
				width = width - cs;
				cSize[chunkIndex] = cs;
				if ( (cSize[chunkIndex-1]==cSize[chunkIndex]) && (cSize[chunkIndex-1]==2) && ( invalid == false) ){
					invalid = true; /* invalidate the current splitting */
				}
				chunkIndex++; /* as this is not the last pair of chunks,
					          pass to the next pair */
		}
			}
			
		REPORT(DETAILED, "Optimized splitting algorithm result: " << (invalid?"failure":"succeeded"));
		
		if (invalid){
			k = -1;
			return false;
		}else{
			k=chunkIndex+1;
			return true;
		}
	}

//	int IntAdder::optimizedLUTCountShortLatency(Target *target, int wIn, int k){
//		int lcount = wIn;
//		lcount += wIn - ( cSize[0] + cSize[k-1] );
//		lcount += 2(k-2);
//		lcount += wIn  - cSize[0];
//		lcount += cSize[0] + cSize[1] + cSize[k-1];
//		 
//		lcount =  - cSize[0]   + 2(k-2)   + cSize[1] ;
//		lcount = 3*wIn - cSize[0] + cSize[1]  + 2*(k-2);   
//	
//	
//	}


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





//			//**********************************************************************
//			//**********************************************************************
//			//SECOND VERSION: CPA2
//			else if ( aType == 2 ){

//				maxInputDelay = getMaxInputDelays (inputDelays);
//				if (verbose)
//					cout << "The maximum input delay is: " << maxInputDelay <<
//					endl;
//		
//				if ( maxInputDelay > objectivePeriod ){
//					//the maximum combinatorial delay of the input is larger than the objective period, so the requested frequency might not be reached.
//					cout << "WARNING: the combinatorial delay at the input of " << this->getName() << " is above objective period "<<endl;
//					maxInputDelay = objectivePeriod;
//				}

//				if ( ((objectivePeriod - maxInputDelay) - target->adderDelay(1) ) < 0 )	{
//					//if not enough slack is available for any combinatorial circuit, we register the inputs
//					nextCycle();
//					target->suggestSubaddSize(chunkSize_ ,wIn_);
//					nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
//					lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
//				}
//				else{
//					int maxInAdd;
//					//explore 2 designs and chose the best				
//					target->suggestSlackSubaddSize(maxInAdd, wIn_, maxInputDelay); 
//					int nbOfChunksFirstDesign = ceil(double(wIn_)/double(maxInAdd));
//					int scoreFirstDesign = nbOfChunksFirstDesign - 1;
//					if (verbose) cout << "Exploring first design ... score is:"<< scoreFirstDesign << endl;
//				
//					target->suggestSubaddSize(maxInAdd, wIn_); 
//					int nbOfChunksSecondDesign = ceil(double(wIn_)/double(maxInAdd));
//					int scoreSecondDesign = nbOfChunksSecondDesign;
//					if (verbose) cout << "Exploring second design ... score is:"<< scoreSecondDesign << endl;
//				
//					if (scoreFirstDesign > scoreSecondDesign){
//						if (verbose) cout << "Implementation of the second design" << endl;
//						nbOfChunks = nbOfChunksSecondDesign;
//						target->suggestSubaddSize(chunkSize_, wIn_); 
//						lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
//					}else{
//						if (verbose) cout << "Implementation of the first design" << endl;
//						nbOfChunks = nbOfChunksFirstDesign;
//						target->suggestSubaddSize(chunkSize_, wIn_); 
//						lastChunkSize = ( wIn_ % chunkSize_ == 0 ? chunkSize_ : wIn_ % chunkSize_);
//					}
//				}
//				//the sizes of the chunks
//				cSize = new int[nbOfChunks+1];
//				if ( nbOfChunks > 1 ){
//					for (int i=0; i<nbOfChunks-1; i++)
//						cSize[i] = chunkSize_;
//					cSize[nbOfChunks-1] = lastChunkSize;
//				}
//				else{
//					nbOfChunks = 1;
//					cSize = new int[1];
//					cSize[0] = wIn_;
//				}
//				//the indexes in the inputs of the chunks
//				cIndex = new int[nbOfChunks];
//				cIndex[0]= cSize[0];
//				for (int i=1; i < nbOfChunks; i++)
//					cIndex[i] = cIndex[i-1] + cSize[i];
//		
//				if (verbose){
//					cout << "The chunk sizes[MSB-->LSB]: "<<endl;
//					for (int i=nbOfChunks-1;i>=0;i--)
//						cout<<cSize[i]<<" ";
//					cout<<endl; 
//					cout << "The index sizes[MSB-->LSB]: "<<endl;
//					for (int i=nbOfChunks-1;i>=0;i--)
//						cout<<cIndex[i]<<" ";
//					cout<<endl; 
//				}	

//				////////////////////////////////////////////////////////////////////////
//		
//				for (int i=0; i < nbOfChunks; i++){
//					vhdl << tab << declare (join("sum_l",0,"_idx",i), cSize[i]+1, true) << " <= " << "( \"0\" & " << use("X") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ") + "
//						  << "( \"0\" & " << use("Y") << range(cIndex[i]-1, (i>0?cIndex[i-1]:0)) << ")" ;
//					if (i==0) vhdl << " + " << use("Cin");
//					vhdl << ";" << endl;
//				}
//		
//				for (int i=1; i <= nbOfChunks-1 ; i++){
//					nextCycle(); ///////////////////////////////////////////////////////
//					for (int j=i; j <= nbOfChunks-1; j++){
//						vhdl << tab << declare(join("sum_l",i,"_idx",j), cSize[j]+1, true) << " <= " << "( \"0\" & " << use(join("sum_l",i-1,"_idx",j))<<range(cSize[j]-1,0) << ") + "
//							  << use(join("sum_l",i-1,"_idx",j-1))<<of(cSize[j-1])<<";"<<endl;
//					}
//				}
//		
//				vhdl << tab << "R <= ";
//				for (int i=nbOfChunks-1; i >= 1; i--){
//					vhdl << use(join("sum_l",i,"_idx",i))<<range(cSize[i]-1,0)<< " & ";
//				}
//				vhdl << use("sum_l0_idx0")<<range(cSize[0]-1,0)<<";"<<endl;

//				outDelayMap["R"] = target->adderDelay(cSize[nbOfChunks-1]); 
//				if (verbose)
//					cout<< "Last addition size is "<<cSize[nbOfChunks-1]<< " having a delay of "<<target->adderDelay(cSize[nbOfChunks-1])<<endl;
//			}

