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

	IntAdder::IntAdder(Target* target, int wIn, map<string, double> inputDelays, int optimizeType, bool srl, int implementation):
		Operator(target), wIn_(wIn), shortLatencyInputRegister(0)
	{
		ostringstream name;
		srcFileName="IntAdder";
		setCopyrightString("Bogdan Pasca, Florent de Dinechin (2008-2010)");		
		name << "IntAdder_" << wIn_<<"_f"<<target->frequencyMHz()
		     <<"_"<<(optimizeType==0? "logic": (optimizeType==1?"reg":"slice")) /*what the operator will be optimized for */
		     <<(srl? "_SRL": "_noSRL"); /* if architecture will be optimized for Shirt Registers */
		setName(name.str());

		// Set up the IO signals
		addInput ("X"  , wIn_, true);
		addInput ("Y"  , wIn_, true);
		addInput ("Cin", 1 );
		addOutput("R"  , wIn_, 1 , true);

		REPORT( INFO, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ");
		REPORT( INFO, "IntAdder wIn="<< wIn << " optimizeType="<<optimizeType << " srl="<< srl << " implementation="<<implementation);
		if (implementation == -1){
			REPORT( INFO, "Optimization focuses on: " << (optimizeType==0? "logic (LUT/ALUT)": (optimizeType==1?" register":"slice")) << " count");
		}else{
			REPORT( INFO, "Forcing implementation: " << (implementation==0? "Classical": (implementation==1?"Alternative":"Short-Latency")) << " Architecture");
		}
		
		REPORT( DEBUG, "Printing input delays:" << printInputDelays(inputDelays));
		REPORT( DEBUG, "Frequeny set to:"<< target->frequency());

		if (isSequential()){	
		
			int selectedImplementation = implementationSelector( target, wIn, inputDelays, optimizeType, srl, implementation);
			REPORT( DEBUG, "................................................................................");
			REPORT( INFO, "The implementation is: " << (selectedImplementation==0? "Classical": (selectedImplementation==1?"Alternative":"Short-Lantency")));
		
			/* shared variables */
			objectivePeriod	 = 1 / target->frequency();	

			//**********************************************************************
			//**********************************************************************
			// Classical implementation of pipelined addition
			if ( selectedImplementation == CLASSICAL ){
				setName( getName() + "_Classical");
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
					REPORT(DEBUG, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
				}else{
					if (classicalSlackVersion == 0){
						/* the slack version that does not buffer the inputs*/			
						REPORT(DETAILED, "Building architecture for classical version with slack: non-buffering");
						updateParameters ( target, inputDelays, alpha, beta, gamma, k);
						REPORT(DEBUG, "alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" k="<<k);
					}else{
						nextCycle(); /* bufferning the inputs */
						REPORT(DETAILED, "Building architecture for classical version with slack: buffering");
						updateParameters (target, alpha, beta, k);
						REPORT(DEBUG, "alpha="<<alpha<<" beta="<<beta<<" k="<<k);
					}
				}
				
				/* init the array with chunk sizes */
				cSize = new int[k+1];
				if ( k > 1 ){
					if (maxInputDelay != 0 )
						cSize[0] = (classicalSlackVersion == 0 ? gamma: alpha);
					else
						cSize[0] = alpha;
						
					for (int i=1; i<k-1; i++)
						cSize[i] = alpha;
					cSize[k-1] = beta;
				}else{
					k = 1;
					cSize = new int[1];
					cSize[0] = wIn_;
				}

				/* the indexes of the chunks */
				cIndex = new int[k];
				cIndex[0] = cSize[0];
				for (int i=1; i < k; i++)
					cIndex[i] = cIndex[i-1] + cSize[i];
		
				ostringstream verb;
				verb << "The chunk sizes[MSB-->LSB]: ";
				for (int i=k-1;i>=0;i--)
					verb << cSize[i]<<" ";
				REPORT( DEBUG, verb.str());
				verb.str("");
				verb << "The index sizes[MSB-->LSB]: ";
				for (int i=k-1;i>=0;i--)
					verb<<cIndex[i]<<" ";
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
				outDelayMap["R"] = target->adderDelay(cSize[k-1]) + (getCurrentCycle()>0?0:getMaxInputDelays(inputDelays)); 
			}
			
			//**********************************************************************
			//**********************************************************************
			// Alternative implementation of pipelined addition
			else if (selectedImplementation == ALTERNATIVE) {
				setName( getName() + "_Alternative");

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

			//**********************************************************************
			//**********************************************************************
			// Alternative implementation of pipelined addition
			else if (selectedImplementation == SHORTLATENCY) {
				setName( getName() + "_ShortLatency");
				REPORT( DEBUG, "SHORT-LATENCY version " << shortLatencyVersion);
				k = shortLatencyKValue;
				
				ostringstream verb;
				verb << "Implementing ShortLatency with chunks: ";
				for (int i=k-1; i>=0; i--)
					verb << " " << cSize[i];
			
				REPORT(DETAILED, verb.str());
			
				//split	 the inputs 
				for (int i=0;i<2;i++)
					for (int j=0; j<k; j++){	
						ostringstream name;
						//the naming standard: sX j _ i _ l
						//j = the chunk index i is the input index and l is the current level
						name << "sX"<<j<<"_"<<i<<"_l"<<0;
						int low=0, high=0;
						for (int k=0;k<=j;k++)
							high+=cSize[k];
						for (int k=0;k<=j-1;k++)
							low+=cSize[k];
						vhdl << tab << declare (name.str(),cSize[j]+1) << " <= " << " \"0\" & "<<(i==0?"X":"Y")<<range(high-1,low)<<";"<<endl;
					}
				vhdl << tab << declare("scIn",1) << " <= " << "Cin;"<<endl;
			
//				if (shortLatencyInputRegister ==1)
//					nextCycle();///////////////////
			
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

					IntAdder *adder = new IntAdder(target, cSize[j]+1, emptyDelayMap, optimizeType, srl);
					ostringstream a;
					a.str("");

					vhdl << tab << declare( uname1.str()+"ext", cSize[j]+1 ) << " <= " <<  "\"0\" & "  <<  use(uname1.str()) << range(cSize[j]-1,0) << ";" << endl;
					vhdl << tab << declare( uname2.str()+"ext", cSize[j]+1 ) << " <= " <<  "\"0\" & "  <<  use(uname2.str()) << range(cSize[j]-1,0) << ";" << endl;
			
					if (j>0){ //for all chunks greater than zero we perform this additions
							bool found = false;
							for(unsigned kk=0; kk<oplist.size(); kk++) {
								if  ( (oplist[kk]->getName()).compare(adder->getName()) ==0 ){
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
				
//						if (j<k-1){
					
							inPortMapCst(adder, "X", uname1.str()+"ext" );
							inPortMapCst(adder, "Y", uname2.str()+"ext" );
							inPortMapCst(adder, "Cin", "'1'");
							outPortMap(adder, "R", dnameOne.str());
							a.str("");
							a<< "Adder" << cSize[j]+1 << "One" << j;
							vhdl << instance( adder, a.str() );
//						}
					}else{
						vhdl << tab << "-- the carry resulting from the addition of the chunk + Cin is obtained directly" << endl;
						vhdl << tab << declare(dnameCin.str(),cSize[j]+1) << "  <= (\"0\" & "<< use(uname1.str())<<range(cSize[j]-1,0) <<") +  (\"0\" & "<< use(uname2.str())<<range(cSize[j]-1,0)<<") + "<<use("scIn")<<";"<<endl;
					}
					
					if (pipe)
						target->setPipelined();
					else
						target->setNotPipelined();

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

				if (k >2){			
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

					if (shortLatencyVersion > 1 )
						nextCycle();/////////////////////

					vhdl << tab << "--perform the short carry additions" << endl; //TODO: PIPELINE ADDITION
					ostringstream unameCin;
					unameCin  << "sX"<<0<<"_0_l"<<l<<"_Cin";
					vhdl << tab << declare("rawCarrySum",2*(k-2)) << " <= " << use("carryStringOne") << " + " << use("carryStringZero") << " + " << use(unameCin.str()) << "(" << cSize[0] << ")" << " ;" << endl;

					if (shortLatencyVersion > 2)
						nextCycle();/////////////////////
				}

				vhdl << tab <<"--get the final pipe results"<<endl;
				for ( int i=0; i<k; i++){
					ostringstream unameZero, unameOne, unameCin;
					unameZero << "sX"<<i<<"_0_l"<<l<<"_Zero";
					unameOne  << "sX"<<i<<"_0_l"<<l<<"_One";
					unameCin  << "sX"<<0<<"_0_l"<<l<<"_Cin";
					if (i==0) vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << use(unameCin.str())<< range(cSize[0]-1,0) <<  ";" << endl;
					else {
//						if (i==1) vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " + " << use(unameCin.str()) << "(" << cSize[0] << ")" << ";"<<endl;
//						else      vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " + not(" << use("rawCarrySum")<<"("<<2*(i-2)+1<<"));"<<endl;
						if (i==1) vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " when " << use(unameCin.str()) << "(" << cSize[0] << ")='0'" 
						                                                              << " else " << use(unameOne.str()) << range(cSize[i]-1,0) << ";"<<endl;
						else      vhdl << tab << declare(join("res",i),cSize[i],true) << " <= " << use(unameZero.str()) << range(cSize[i]-1,0) << " when " << use("rawCarrySum")<<"("<<2*(i-2)+1<<")='1'"
						                                                              << " else " << use(unameOne.str()) << range(cSize[i]-1,0) << ";"<<endl;
					}
				}
			
				vhdl << tab << "R <= ";
				for (int i=k-1; i>=0; i--){
					vhdl << use(join("res",i));
					if (i > 0) vhdl << " & ";
				}
				vhdl << ";" <<endl;			
				outDelayMap["R"] = target->adderDelay(cSize[k-1]); 			
			
				REPORT(DEBUG, " FINISHED SHORT-LATENCY IMPLEMENTATION" );
			}
			else {
				vhdl << tab << " R <= X + Y + Cin;" << endl;
				outDelayMap["R"] = target->adderDelay(wIn_); 
			}
		// COMBINATORIAL VERSION	
		}else{
			vhdl << tab << " R <= X + Y + Cin;" << endl;
			outDelayMap["R"] = target->adderDelay(wIn_) + getMaxInputDelays(inputDelays); 
		}
	
		REPORT(DEBUG, "OutDelay for R is " << outDelayMap["R"]);	
	
	}

	IntAdder::~IntAdder() {
	}


	int IntAdder::implementationSelector(Target* target, int wIn, map<string, double> inputDelays, int optimizeType, bool srl, int implementation){
	
		if ((target->getID()=="Virtex4")||(target->getID()=="Virtex5")||(target->getID()=="Spartan3")){
			/* XILINX TARGETS */
			if (implementation == -1){
				if (optimizeType == LOGIC){
					int lutCostClassical     = getLutCostClassical(target, wIn, inputDelays, srl);
					int lutCostAlternative   = getLutCostAlternative(target, wIn, inputDelays, srl);
					int lutCostShortLantency = getLutCostShortLatency(target, wIn, inputDelays, srl);
				
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "LUT cost for the classical     method " << lutCostClassical );
					REPORT(DETAILED, "LUT cost for the alternative   method " << lutCostAlternative );
					REPORT(DETAILED, "LUT cost for the short-latency method " << lutCostShortLantency	 );
				
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
					int regCostClassical     = getRegCostClassical(target, wIn, inputDelays, srl);
					int regCostAlternative   = getRegCostAlternative(target, wIn, inputDelays, srl);
					int regCostShortLantency = getRegCostShortLatency(target, wIn, inputDelays, srl);
				
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "REG cost for the classical     method " << regCostClassical );
					REPORT(DETAILED, "REG cost for the alternative   method " << regCostAlternative );
					REPORT(DETAILED, "REG cost for the short-latency method " << regCostShortLantency	 );
				
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
					int sliceCostClassical     = getSliceCostClassical(target, wIn, inputDelays, srl);
					int sliceCostAlternative   = getSliceCostAlternative(target, wIn, inputDelays, srl);
					int sliceCostShortLantency = getSliceCostShortLatency(target, wIn, inputDelays, srl);
				
					REPORT(DEBUG   , "................................................................................");
					REPORT(DETAILED, "SLICE cost for the classical     method " << sliceCostClassical );
					REPORT(DETAILED, "SLICE cost for the alternative   method " << sliceCostAlternative );
					REPORT(DETAILED, "SLICE cost for the short-latency method " << sliceCostShortLantency	 );
				
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
			}else if (implementation == CLASSICAL){
				if (optimizeType == LOGIC){
					int lutCostClassical = getLutCostClassical(target, wIn, inputDelays, srl);
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "LUT cost for the classical     method " << lutCostClassical );
					return CLASSICAL;				
				} else if (optimizeType == REGISTER){
					int regCostClassical = getRegCostClassical(target, wIn, inputDelays, srl);
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "REG cost for the classical     method " << regCostClassical );
					return CLASSICAL;
				} else if (optimizeType == SLICE){
					int sliceCostClassical = getSliceCostClassical(target, wIn, inputDelays, srl);
					REPORT(DEBUG   , "................................................................................");
					REPORT(DETAILED, "SLICE cost for the classical     method " << sliceCostClassical );
					return CLASSICAL;
				}
			}else if (implementation == ALTERNATIVE){
				if (optimizeType == LOGIC){
					int lutCostAlternative   = getLutCostAlternative(target, wIn, inputDelays, srl);
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "LUT cost for the alternative   method " << lutCostAlternative );
					return ALTERNATIVE;
				} else if (optimizeType == REGISTER){
					int regCostAlternative   = getRegCostAlternative(target, wIn, inputDelays, srl);
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "REG cost for the alternative   method " << regCostAlternative );
					return ALTERNATIVE;
				} else if (optimizeType == SLICE){
					int sliceCostAlternative   = getSliceCostAlternative(target, wIn, inputDelays, srl);
					REPORT(DEBUG   , "................................................................................");
					REPORT(DETAILED, "SLICE cost for the alternative   method " << sliceCostAlternative );
					return ALTERNATIVE;
				}
			}else if (implementation == SHORTLATENCY){
				if (optimizeType == LOGIC){
					int lutCostShortLantency = getLutCostShortLatency(target, wIn, inputDelays, srl);
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "LUT cost for the short-latency method " << lutCostShortLantency );
//					if (lutCostShortLantency != PINF)
						return SHORTLATENCY;
//					else
//						return CLASSICAL;
				} else if (optimizeType == REGISTER){
					int regCostShortLantency = getRegCostShortLatency(target, wIn, inputDelays, srl);
					REPORT( DEBUG, "................................................................................");
					REPORT(DETAILED, "REG cost for the short-latency method " << regCostShortLantency	 );
//					if (regCostShortLantency != PINF)
						return SHORTLATENCY;
//					else
//						return ALTERNATIVE;
				} else if (optimizeType == SLICE){
					int sliceCostShortLantency = getSliceCostShortLatency(target, wIn, inputDelays, srl);
					REPORT(DEBUG   , "................................................................................");
					REPORT(DETAILED, "SLICE cost for the short-latency method " << sliceCostShortLantency	 );
//					if (sliceCostShortLantency != PINF)
						return SHORTLATENCY;
//					else
//						return ALTERNATIVE;
				}
			}						
		}else{
			/* ALTERA TARGETS */
		}
		
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}
	
/******************************************************************************/	
	int IntAdder::getLutCostClassical(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters(target, alpha, beta, k);
			REPORT( DEBUG, "LUT, Classical, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (srl){
				if ((k == 1) || (k == 2)) {
					cost = wIn;
				}else{
					/* more than two chunk splitting */
					cost = wIn + (k-2)*alpha;
				}
			}else{
				cost = wIn;
			}
			REPORT( DETAILED, "Selected: Classical, NO-SLACK with LUT cost " << cost );
			return cost;
		}else{
			int version0, version1;
			int alpha, beta, gamma, k;			
			/* Version 0 */
			updateParameters(target, inputDelays, alpha, beta, gamma, k);
			REPORT( DEBUG, "LUT, Classical, SLACK, Version 0: alpha="<<alpha<<" beta`="<<beta<<" gamma="<<gamma<<" k="<<k);
			
			if (k>0)
				if (srl){
					if ((k==1) || (k==2)){			
						version0 = wIn;
					}else{
						version0 = 4*wIn - 3*alpha - 2*gamma  - beta;
					}
				}else{
					/* NO SRLs */
					version0 = wIn;
				}
			else
				version0 = PINF; //infinity
			
			/* Version 1 */
			updateParameters(target, alpha, beta, k);
			REPORT( DEBUG, "LUT, Classical, SLACK, Version 1 (buffered inputs): alpha="<<alpha<<" beta="<<beta<<" k="<<k << " p="<<k+1);

			if (srl){
				if (k==1){			
					version1 = wIn;
				}else if (k == 2){
					version1 = w + 2*beta;
				}else{
					version1 = 4*wIn - 3*alpha - beta;
				}
			}else{
				version1 = wIn;
			}		
			REPORT(DETAILED, "LUT, Classical, SLACK, Version 0: " << version0);
			REPORT(DETAILED, "LUT, Classical, SLACK, Version 1 (buffered inputs) " << version1);
			
			if (version0 <= version1){
				/* for equality version 0 has less pipeline levels */
				classicalSlackVersion = 0;
				REPORT( DETAILED, "Selected: Classical SLACK version is 0 with LUT cost " << version0 );
				return version0;
			}else{
				classicalSlackVersion = 1;
				REPORT( DETAILED, "Selected: Classical SLACK version is 1 (buffered inputs) with LUT cost " << version1 );
				return version1;
			}
		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

/******************************************************************************/
	int IntAdder::getLutCostAlternative(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "LUT, Alternative, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (srl){
				if (k == 1) {
					cost = wIn;
				}else if (k == 2){
					cost = (alpha + 2*beta);
				}else{
					/* more than two chunk splitting */
					cost = ((4*k-8)*alpha + 3*beta + k-3);
				}
			}else{
				cost = (2*wIn - alpha);
			}
			REPORT( DETAILED, "Selected: Alternative, NO-SLACK with LUT cost " << cost );
			return cost;
		}else{
			int version0, version1;
			int alpha, beta, k, kVersion0, kVersion1;			
			
			/* Version 0 */
			updateParameters( target, inputDelays, alpha, beta, k);
			kVersion0 = k;
			REPORT( DEBUG, "LUT, Alternative, SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (k>0)
				if (srl){
					if (k==1){			
						version0=wIn;
					}else if (k == 2){
						version0 = alpha + 2*beta;
					}else{
						version0 = (4*k-8)*alpha + 3*beta + k-3;
					}
				}else{
					version0 = 2*wIn - alpha;
				}
			else
				version0 = PINF; //infinity

			/* Version 1 */			
			updateParameters( target, alpha, beta, k);
			kVersion1 = k;
			REPORT( DEBUG, "LUT, Alternative, SLACK, Version 1 (buffered inputs): alpha="<<alpha<<" beta="<<beta<<" k="<<k << " p="<<k+1);

			if (srl){
				if (k==1){			
					version1=wIn;
				}else if (k == 2){
					version1 = alpha + 2*beta;
				}else{
					version1 = (4*k-8)*alpha + 3*beta + k-3;
				}
			}else{
				version1 = 2*wIn - alpha;
			}			
			
			REPORT(DETAILED, "LUT, Alternative, SLACK, Version 0: " << version0);
			REPORT(DETAILED, "LUT, Alternative, SLACK, Version 1 (buffered inputs): " << version1);
			
			if (version0 == version1){
				if (kVersion0 <= kVersion1){
					alternativeSlackVersion = 0;
					REPORT( DETAILED, "Alternative SLACK, version is 0 with LUT cost " << version0 << " and k="<<kVersion0 );
					return version0;
				}else{
					alternativeSlackVersion = 1;
					REPORT( DETAILED, "Alternative SLACK, version is 1 with LUT cost " << version1 << " and k="<<kVersion1 );
					return version1;
				}
			}else if (version0 < version1){
				alternativeSlackVersion = 0;
				REPORT( DETAILED, "Alternative SLACK version is 0 with LUT cost " << version0 );
				return version0;
			}else{
				alternativeSlackVersion = 1;
				REPORT( DETAILED, "Alternative SLACK version is 1 with LUT cost " << version1 );
				return version1;
			}
		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

/******************************************************************************/
	int IntAdder::getLutCostShortLatency(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		tryOptimizedChunkSplittingShortLatency ( target, wIn , k );
		REPORT( DEBUG, "LUT, Short-Latency: k="<<k);
		
		int cost;
		if ( getMaxInputDelays( inputDelays ) == 0 ){
		/* no input slack problem */
			if ( shortLatencyVersion == 0){
				cost = cSize[0] + 3*cSize[1];
			} else if ( shortLatencyVersion == 1){
				cost = 3*wIn - 2*cSize[0] + 2*(k-2);
			} else if ( shortLatencyVersion == 2){
				/* the algorithm found a good way to split inputs and save 1 pipeline level */
				cost = (3*wIn - 2*cSize[0] - cSize[k-1]  + 2*(k-2));
			}else if ( shortLatencyVersion == 3){
				if (k>=3){
					if (srl)
						cost = ((4*k - 6)*cSize[0] + 3*cSize[k-1] + 2*k-4);
					else
						cost = (3*wIn-2*cSize[0]+2*(k-2));
				}else
					cost = PINF; //+inf			
			}else
				cost = PINF; //+inf
			
			REPORT( DETAILED, "Selected: Short-Latency, "<< (shortLatencyVersion==0?"optimized splitting":"default splitting") <<" with LUT cost " << cost );
			return cost;
			
		}else{
			//TODO No implementation for now, so cost is +inf
			return PINF;
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

/******************************************************************************/
	int IntAdder::getRegCostClassical(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "REG, CLASSICAL, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			if (k == 1) {
				cost = 0;
			}else{
				if (srl)
					cost = w - beta;
				else
					cost = ((3*k*k -7*k+4)*alpha/2 + 2*(k-1)*beta + k-1);
			}
			
			REPORT( DETAILED, "Selected: Classical, NO-SLACK, with REG cost " << cost );
			return cost;

		}else{
			int version0, version1;
			int alpha, beta, gamma, k;			
			
			/* Version 0 */
			updateParameters( target, inputDelays, alpha, beta, gamma, k);
			REPORT( DEBUG, "REG, Classical, SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" k="<<k);
			
			if (k>0){
				if (k==1)			
					version0 = 0;

				if (srl){
					if (k == 2)
						version0 = wIn + alpha + 1;
					else
						version0 = 3*wIn - 2*gamma  - beta + k - 1;
				}else{
					if (k == 2){
						version0 = wIn + alpha + 1;
					}else{
						version0 = (k-1)*gamma + 2*(k-1)*beta + 3*(k*k-3*k+2)/2;
					}
				}
				
			}else
				version0 = PINF; //infinity
			
			/* Version 0 */
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "REG, Classical, SLACK, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (srl){
				if (k==1){			
					version1 = wIn;
				}else if (k == 2){
					version1 = 3*alpha + 2*beta + 2;
				}else{
					version1 = 3*wIn - beta + k - 1;
				}
			}else{
				if (k==1)
					version1 = wIn;
				else
					version1 = wIn + ( (3*k*k -7*k+4)*alpha/2 + 2*(k-1)*beta + k-1);
			}
			
			REPORT(DETAILED, "REG, Classical, SLACK, Version 0: " << version0);
			REPORT(DETAILED, "REG, Classical, SALCK, Version 1 (buffered inputs): " << version1);
		
			if (version0 <= version1){
				classicalSlackVersion = 0;
				REPORT( DETAILED, "Selected: Classical SLACK version is 0 with REG cost " << version0 );
				return version0;
			}else{
				classicalSlackVersion = 1;
				REPORT( DETAILED, "Selected: Classical SLACK version is 1 with REG cost " << version1 );
				return version1;
			}
		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

/******************************************************************************/
	int IntAdder::getRegCostAlternative(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "REG, ALTERNATIVE, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			
			if (k == 1) {
				cost = 0;
			}else{
				if (srl)
					cost = ((2*k-3)*alpha + beta + 2*k-3);
				else
					cost = ((k-1)*wIn + k*k-2*k+1);
			}

			REPORT( DETAILED, "Selected: Alternative, NO-SLACK, with REG cost " << cost );
			return cost;
		}else{
			int version0, version1;
			int alpha, beta, k;	
			/* Version 0 */		
			updateParameters( target, inputDelays, alpha, beta, k);
			REPORT( DEBUG, "REG, Alternative, SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<" k="<<k << " p="<<k+1);
			
			if (k>0)
				if (k==1){			
					version0 = 0;
				}else{
					if (srl)
						version0 = (2*k-3)*alpha + beta + 2*k-3;
					else
						version0 = ((k-1)*w + k*k-3*k+3);
				}
			else
				version0 = PINF; //infinity
			
			/* Version 1 */
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "REG, Alternative, SLACK, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (k==1){			
				version1 = 0;
			}else{
				if (srl)
					version1 = (4*k-8)*alpha + 3*beta + 2*k-3 + 2*wIn+1;
				else
					version1 = wIn + ((k-1)*w + k*k-3*k+3);
			}
			
			REPORT(DETAILED, "REG, Alternative, Version 0: " << version0);
			REPORT(DETAILED, "REG, Alternative, Version 1 (buffered inputs): " << version1);
			
			if (version0 <= version1){
				alternativeSlackVersion = 0;
				REPORT( DETAILED, "Selected: Alternative slack version is 0 with cost " << version0 );
				return version0;
			}else{
				alternativeSlackVersion = 1;
				REPORT( DETAILED, "Selected: Alternative slack version is 1 with cost " << version1 );
				return version1;
			}
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

/******************************************************************************/
	int IntAdder::getRegCostShortLatency(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		tryOptimizedChunkSplittingShortLatency ( target, wIn , k );
		REPORT( DEBUG, "REG, Short-Latency:  k="<<k);
		int cost;
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			if (( shortLatencyVersion == 0 ) || ( shortLatencyVersion == 1 )){
				cost = 0;
			}else if ( shortLatencyVersion == 2 ){
				if (k>=3){
					if (srl)
						cost = ((k - 1)*cSize[0] + cSize[k-1] + 2*k-7);
					else
						cost = (2*wIn + 3*k - 5);
				}else
					cost = PINF; 		
			}else if ( shortLatencyVersion == 3 ){
				if (k>=3){
					if (srl)
						cost = ((k - 1)*cSize[0] + cSize[k-1] + 4*k-7);
					else
						cost = (2*wIn + 3*k - 5);
				}else
					cost = PINF; 		
			} else 
				cost = PINF;
			
			REPORT( DETAILED, "Selected: Short-Latency with REG cost " << cost );
			return cost;

		}else{
			/* TODO for slack */
			return PINF;
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getSliceCostClassical(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");		
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "SLICE, CLASSICAL, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			
			if (k == 1) {
				cost = int(ceil( double(wIn) / double(2) ) );
			}else{
				if (srl)
					cost = int(ceil(double(wIn + (k-2)*alpha  + k - 1) / double(2)));
				else
					cost = int(ceil(double(wIn + ((3*(k*k-3*k+2))/2)*alpha + 2*(k-1)*beta) / double(2)));
			}
			REPORT( DETAILED, "Selected: Classical NO-SLACK, with SLICE cost " << cost );
			return cost;

		}else{
			int version0, version1;
			int alpha, beta, gamma, k;			
			
			/* Version 0 */
			updateParameters( target, inputDelays, alpha, beta, gamma, k);
			REPORT( DEBUG, "SLICE, Classical, SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<"gamma="<<gamma<< " k="<<k);
			
			if (k>0){
				if (k==1)			
					version0= int(ceil(double(wIn)/double(2)));
					
				if (srl){	
					if (k == 2){
						version0 = int(ceil(double(3*alpha + gamma + 1)/double(2)));
					}else{
						version0 = int(ceil(double(4*wIn - alpha - 2*gamma -2*beta + k - 1)/double(2)));
					}
				}else{
					if (k == 2){
						version0 = int(ceil(double(gamma + 3*alpha + 1) / double(2)));
					}else{
						version0 = int(ceil(double(wIn +(k-2)*gamma + 2*(k-1)*beta + alpha*(2*k*k-11*k+10)/2)/double(2)));
					}
				}
			}else
				version0 = PINF; //infinity

			/* Version 1 */			
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "SLICE, Classical, SLACK Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k<<" p="<<k+1);

			if (k==1)			
				version1= int(ceil(double(wIn)/double(2)));
			
			if (srl){
				if (k == 2){
					version1 = int(ceil(double(3*alpha + 2*beta + 2)/double(2)));
				}else{
					version1= int(ceil(double(4*wIn - alpha - beta)/double(2)));
				}
			}else{
				if (k == 2){
					version1 = int(ceil(double(2*wIn + ((3*(k*k-3*k+2))/2)*alpha + 2*(k-1)*beta) / double(2)));
				}else{
					version1 = int(ceil(double(3*wIn + ((3*(k*k-3*k+2))/2)*alpha + 2*(k-1)*beta) / double(2)));
				}
			}
			
			REPORT(DETAILED, "SLICE, Classical, SLACK, Version 0: " << version0);
			REPORT(DETAILED, "SLICE, Classical, SLACK, Version 1: " << version1);
			
			if (version0 <= version1){
				classicalSlackVersion = 0;
				REPORT( DETAILED, "Selected: Classical SLACK version is 0 with SLICE cost " << version0 );
				return version0;
			}else{
				classicalSlackVersion = 1;
				REPORT( DETAILED, "Selected: Classical SLACK version is 1 with SLICE cost " << version1 );
				return version1;
			}
		}
	
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getSliceCostAlternative(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		REPORT( DEBUG, "................................................................................");
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "SLICE, ALTERNATIVE, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k );			

			if (k == 1) {
				cost = int(ceil( double(wIn) / double(2) ) );
			} else if (srl){
				if (k == 2 ){
					cost = int(ceil(double(alpha + 2*beta + 1) / double(2)));
				}else{
					cost = int(ceil(double((4*k-8)*alpha + 3*beta + 2*k - 5) / double(2)));
				}
			}else{
				cost = int(ceil(double((k-1)*wIn + beta + k*k-2*k+1) / double(2)));
			}
			
			REPORT( DETAILED, "Selected: Alternative, NO-SLACK with SLICE cost " << cost );
			return cost;
		}else{
			int version0, version1;
			int alpha, beta, k;			
			updateParameters( target, inputDelays, alpha, beta, k);
			REPORT( DEBUG, "SLICE, Alternative, NO-SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<" k="<<k);
			
			if (k>0){
				if (k==1){			
					version0=0;
				}else if (srl){
					if (k==2){
						version0 = int(ceil(double(alpha + 2*beta + 1 )/2.));
					}else{
						version0 = int(ceil(double((4*k-8)*alpha + 3*beta + 2*k-5 )/2.));
					}
				}else{
					version0 = int(ceil(double(2*wIn + (k*k-4*k+3)*alpha + (k-2)*beta) / double(2)));
				}
			}else
				version0 = PINF; //infinity
			
			updateParameters( target, alpha, beta, k);
			REPORT( DEBUG, "SLICE, Alternative, SLACK, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k);

			if (k==1){			
				version1 = 0;
			}else if (srl){ 
				if (k==2){
					version1 = int(ceil(double(alpha + 2*beta + 1 + 2*wIn + 1)/2.));
				}else{
					version1 = int(ceil(double((4*k-8)*alpha + 3*beta + 2*k-5 + 2*wIn + 1)/2.));
				}
			}else{
				version1 = int(ceil(double(4*wIn + (k*k-4*k+3)*alpha + (k-2)*beta) / double(2)));
			}	
			REPORT(DETAILED, "SLICE, Alternative, SLACK, Version 0: " << version0);
			REPORT(DETAILED, "SLICE, Alternative, SLACK, Version 1: " << version1);
		
			if (version0 <= version1){
				alternativeSlackVersion = 0;
				REPORT( DETAILED, "Alternative SLACK version is 0 with SLICE cost " << version0 );
				return version0;
			}else{
				alternativeSlackVersion = 1;
				REPORT( DETAILED, "Alternative SLACK version is 1 with SLICE cost " << version1 );
				return version1;
			}		
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit(-1);
		return -1; 
	}

	int IntAdder::getSliceCostShortLatency(Target* target, int wIn, map<string, double> inputDelays, bool srl){
		shortLatencyInputRegister = -1;
		REPORT( DEBUG, "................................................................................");
		tryOptimizedChunkSplittingShortLatency ( target, wIn , k );
		REPORT( DEBUG, "SLICE, Short-Latency:  k="<<k);
		int cost;
		if ( getMaxInputDelays( inputDelays ) == 0 ){
			/* no input slack problem */
			if ( shortLatencyVersion == 0){
				cost = int(ceil(double(cSize[0] + 3*cSize[1])/double(2)));
			} else if ( shortLatencyVersion == 1){
				cost = int(ceil(double(3*wIn - 2*cSize[0] + 2*(k-2))/double(2)));
			}else if ( shortLatencyVersion == 2){
				ostringstream tmp;
				for (int y=k-1;y<=0;y--)
					tmp << cSize[y] << ":";
				
				REPORT(DEBUG, "SLICE, Short-Latency cSize=> " << tmp.str() );
			
//				shortLatencyVersion = 0;
				if (k>=3){
					if (srl)
						cost = int(ceil(double( 3*wIn - 2*cSize[0] - cSize[k-1] + 2*(k-1) + 1 ) / double(2)));
					else
						cost = int(ceil(double( 4*wIn - 2*cSize[0] - cSize[k-1] + 5*k -9 ) / double(2)));
				}else
					cost = PINF; //TODO			
			}else if ( shortLatencyVersion == 3){
				if (k>=3){
					if (srl)
						cost = int(ceil(double((4*k-6)*cSize[0] + 3*cSize[k-1] + 2*k - 4) / double(2)));
					else
						cost = int(ceil(double(4*wIn - 2*cSize[0] - cSize[k-1] +5*k-9) / double(2)));
				}else
					cost = PINF; //TODO			
			}
			REPORT( DETAILED, "Selected: Short-Latency with SLICE cost " << cost );
			return cost;
		}else{
			shortLatencyInputRegister = 1;
//			shortLatencyVersion = 1;
			/* TODO for slack */
			
			
			return PINF;
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

	void IntAdder::tryOptimizedChunkSplittingShortLatency(Target* target, int wIn, int &k){
		cSize = new int[2000];
		for (int u=0; u<2000; u++)
			cSize[u] = -1;

		int alpha0;
		double tSelect = target->lutDelay() + target->localWireDelay();

		double k1,k2; 
		target->getAdderParameters(k1,k2);
//		alpha0 = int(floor(((1.0 / target->frequency()) - tSelect - k1) / k2));
		target->suggestSlackSubaddSize(alpha0, wIn, tSelect);
		int alpha;
		target->suggestSubaddSize(alpha,wIn);
		
		double C =  ((1.0 / target->frequency()) - tSelect - 2*k1 + 2*k2)/k2;
		int U = int (floor (C));
				
		REPORT( DEBUG, "U="<<U<<" C="<<C);
		int maxW;
		if (U < 0)
			U = 0;
		
		if (U >= 0)
			if (U % 2 ==0)
				maxW = 2*alpha0 + U*(U+2)/4;
			else
				maxW = 2*alpha0 + (U+1)*(U+1)/4;
		else
			maxW = -1;
	
		REPORT( DEBUG, "alpha is " << alpha );
		REPORT( DEBUG, "Max addition size for latency 0, two chunk architecture:" << 2*alpha0);	
		REPORT( DEBUG, "Max addition size for latency 0 is:" << maxW);
		
		double C2 = ((1.0 / target->frequency()) - tSelect - k1 + k2)/(2*k2);
		int U2 = int(floor(C2));
		
		double C3 = ((1.0 / target->frequency()) - k1 + k2)/(2*k2);
		int U3 = int(floor(C3));
		
		REPORT(DEBUG, "Max addition size for latency 1 is:" << (U2+2)*alpha0 );
		REPORT(DEBUG, "Max addition size for latency 2 is:" << (U3+2)*alpha );
		
		if (wIn <= 2*alpha0){
			// TWO CHUNK ARCHITECTURE
			cSize[0]= int(floor(wIn/2));
			cSize[1]= wIn - cSize[0];
			REPORT ( DEBUG, " Chunk sizes are: " << cSize[0] << " " << cSize[1]);
			shortLatencyVersion = 0;
			k = 2; 
		}else if (wIn <= maxW){
			// LATENCY 0 Architecture
			cSize[0]=alpha0;
			cSize[1]=alpha0;			
			int tmpWIn = wIn - 2*alpha0;
			int i=2;
			while (tmpWIn > 0){
				if (tmpWIn - ( U-2*(i-2)) >= 0){
					cSize[i]=U-2*(i-2);
					tmpWIn -= cSize[i];
				}else{
					cSize[i] = tmpWIn;
					tmpWIn -= cSize[i];
				}
				i++;	
			}
			k=i;
			ostringstream tmp;
			for (int kk=i-1; kk>=0;kk--){
				tmp <<  cSize[kk] << " ";
			}	
			REPORT( DEBUG, " Chunks " << k <<"  Sizes: " << tmp.str());
			shortLatencyVersion = 1;
		} else if (( wIn > maxW ) && ( wIn <= (U2+2)*alpha0 )){
			//LATENCY 1 architecture
			if ( wIn % alpha0 == 0 )
				k = wIn / alpha0;
			else
				k = int(ceil (double(wIn) / double(alpha0)) );
			
			for (int p=0; p<k-1; p++)
				cSize[p] = alpha0;
			
			cSize[k-1] = ( wIn % alpha0 == 0 ? alpha0 : wIn % alpha0);

			ostringstream tmp;
			for (int kk=k-1; kk>=0;kk--){
				tmp <<  cSize[kk] << " ";
			}	
			REPORT( DEBUG, " Chunk " << k<< "  Sizes: " << tmp.str());
			shortLatencyVersion = 2;
		}  else if ( ( wIn > (U2+2)*alpha0 ) && (wIn <= (U3+2)*alpha)){
			if ( wIn % alpha == 0 )
				k = wIn / alpha;
			else
				k = int(ceil (double(wIn) / double(alpha)) );
			
			for (int p=0; p<k-1; p++)
				cSize[p] = alpha;
			
			cSize[k-1] = ( wIn % alpha == 0 ? alpha : wIn % alpha);

			ostringstream tmp;
			for (int kk=k-1; kk>=0;kk--){
				tmp <<  cSize[kk] << " ";
			}	
			REPORT( DEBUG, " Chunk " << k<< "  Sizes: " << tmp.str());
			shortLatencyVersion = 3;
		} else if  (wIn > (U3+2)*alpha){
		 	//LATENCY 2+ architecture
			if ( wIn % alpha == 0 )
				k = wIn / alpha;
			else
				k = int(ceil (double(wIn) / double(alpha)) );
			
			for (int p=0; p<k-1; p++)
				cSize[p] = alpha;
			
			cSize[k-1] = ( wIn % alpha == 0 ? alpha : wIn % alpha);

			ostringstream tmp;
			for (int kk=k-1; kk>=0;kk--){
				tmp <<  cSize[kk] << " ";
			}	
			REPORT( DEBUG, " Chunk " << k<< "  Sizes: " << tmp.str());
			shortLatencyVersion = 4;
		}

		REPORT( DEBUG, "Selected Short-Latency Version is " << shortLatencyVersion);
		shortLatencyKValue = k;

//		REPORT(DETAILED, "Trying to optimize chunk splitting for short-latency architecture ...");

//		target->suggestSubaddSize(chunkSize_ ,wIn_);
//		REPORT(DEBUG, "The chunkSize for first two chunks would be: " << chunkSize_ );
//	
//		if (2*chunkSize_ >= wIn_){
//			REPORT(DETAILED, "Failed ... the input size is too small to use the short latency algorithm");
//			return false;
//		}
//			
//		cSize[0] = chunkSize_;
//		cSize[1] = chunkSize_;
//	
//		bool finished = false; /* detect when finished the first the first
//			                   phase of the chunk selection algo */
//		int width = wIn - 2*chunkSize_; /* remaining size to split into chunks */
//		int propagationSize = 0; /* carry addition size */
//		int chunkIndex = 2; /* the index of the chunk for which the size is
//			                to be determined at the current step */
//		bool invalid = false; /* the result of the first phase of the algo */
//	
//		/* FIRST PHASE */
//		while (not (finished))	 {
//			propagationSize+=2;
//			double delay = objectivePeriod - target->adderDelay(width)- target->adderDelay(propagationSize); //2*target->localWireDelay()  -
//			if ((delay > 0) || (width < 4)) {
//				cSize[chunkIndex] = width;
//				finished = true;
//			}else{
//				int cs; 
//				double slack =  target->adderDelay(propagationSize) ; //+ 2*target->localWireDelay()
//				target->suggestSlackSubaddSize( cs, width, slack);
//				width = width - cs;
//				cSize[chunkIndex] = cs;
//				if ( (cSize[chunkIndex-1]==cSize[chunkIndex]) && (cSize[chunkIndex-1]==2) && ( invalid == false) ){
//					invalid = true; /* invalidate the current splitting */
//				}
//				chunkIndex++; /* as this is not the last pair of chunks,
//					          pass to the next pair */
//			}
//		}
//			
//		REPORT(DETAILED, "Optimized splitting algorithm result: " << (invalid?"failure":"succeeded"));
//		
//		if (invalid){
//			k = -1;
//			return false;
//		}else{
//			k=chunkIndex+1;
//			REPORT(DETAILED, "Optimized splitting algorithm result: will set k " << k );
//			return true;
//		}
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


