/*
  An unsigned integer multiplier for FloPoCo
  
  Authors: Bogdan Pasca, Sebastian Banescu

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
  All rights reserved.

 */

#include <typeinfo>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "LogicIntMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;
	
	// Stuff for the small multiplier table 
	LogicIntMultiplier::SmallMultTable::SmallMultTable(Target* target, int dx, int dy) : 
		Table(target, dx+dy, dx+dy, 0, -1, true), // logic table
		dx(dx), dy(dy) {
		ostringstream name; 
		srcFileName="LogicIntMultiplier::SmallMultTable";
				name <<"SmallMultTable_" << dy << "x" << dx;
				setName(name.str());
				
				//				outDelayMap["Y"] = target->lutDelay();
	};
	
	mpz_class LogicIntMultiplier::SmallMultTable::function(int yx){
		mpz_class r;
		int y = yx>>dx;
		int x = yx -(y<<dx);
		r = x*y;
		return r;
	};



	string PP(int i, int j) {
		std::ostringstream p;
		p << "PPX" << i << "Y" << j;
		return p.str();
	};

	string PPTbl( int i, int j) {
		std::ostringstream p;
		p << "PPX" << i << "Y" << j << "_Tbl";
		return p.str();
	};

	string XY( int i, int j) {
		std::ostringstream p;
		p  << "Y" << j<< "X" << i;
		return p.str();
	};

	LogicIntMultiplier:: LogicIntMultiplier(Target* target, int wInXp, int wInYp, bool sign, map<string, double> inputDelays) :
		Operator(target, inputDelays), wInX_(wInXp), wInY_(wInYp), wOut_(wInXp + wInYp), sign_(sign){
		
		ostringstream name;
		double target_freq = target->frequency();
		
		double maxDSPFrequency = int(floor(1.0/ (target->DSPMultiplierDelay() + 1.0e-10)));
		if (target_freq > maxDSPFrequency)
			target->setFrequency(maxDSPFrequency);
		
		
		/* Name Setup procedure
		 *  The name has the format: LogicIntMultiplier_wInX_wInY
		 *  wInX = width of the X input
		 *  wInY = width of the Y input
		 */  
		name <<"LogicIntMultiplier_"<<wInX_<<"_"<<wInY_<<"_uid"<<Operator::getNewUId();;
		if (sign_)
			name << "_signed";
		srcFileName="LogicIntMultiplier";
		setName(name.str());
		
		setCopyrightString("Bogdan Pasca, Sebastian Banescu, Florent de Dinechin (2008-2012)");
	
		addInput ("X", wInX_,true);
		addInput ("Y", wInY_,true);
		addOutput("R", wOut_,1, true); /* wOut_ = wInX_ + wInY_ */

		// Halve number of cases by assuming wInY<=wInX:
		// interchange x and y in case wInY>wInX

		if(wInY_> wInX_){
			wInX=wInY_;
			wInY=wInX_;
			vhdl << tab << declare("XX", wInX, true) << " <= Y;" << endl; 
			vhdl << tab << declare("YY", wInY, true) << " <= X;" << endl; 
		}
		else{
			wInX=wInX_;
			wInY=wInY_;
			vhdl << tab << declare("XX", wInX, true) << " <= X;" << endl; 
			vhdl << tab << declare("YY", wInY, true) << " <= Y;" << endl; 
		}

		/***************************************************************************
		 ***************************************************************************/
		//LogicIntMultiplier Version -- multiplication performed in LUTs. Works for both Xilinx and Altera
		//architectures for corner cases 
		if ((wInY == 1)){
			setCriticalPath( getMaxInputDelays( inputDelays ));
			if (sign){
				manageCriticalPath( target->localWireDelay(wInX) + target->adderDelay(wInX+1) );
				vhdl << tab << "R <= (" << zg(wInX+1)  << " - (XX" << of(wInX-1) << " & XX)) when YY(0)='1' else "<< zg(wInX+1,0)<<";"<<endl;	
			}
			else {
				manageCriticalPath( target->localWireDelay(wInX) + target->lutDelay() );
				vhdl << tab << "R <= (\"0\" &  XX) when YY(0)='1' else "<<zg(wInX+1,0)<<";"<<endl;	
			}
				outDelayMap["R"] = getCriticalPath();
		}
		else if ((wInY == 2)) {
			if(target->lutInputs()-2>=wInX){
				if(sign) {
					std::ostringstream o;
					o << srcFileName << " (" << uniqueName_ << "): SORRY : ultrasmall _signed_ LogicIntMultiplier not implemented";
					throw o.str();
				} // FIXME of course
				setCriticalPath( getMaxInputDelays( inputDelays ));
				manageCriticalPath( target->localWireDelay() + target->lutDelay() );
				vhdl << tab << "R <= XX * YY;  -- fits in one look-up table"<<endl;	
				outDelayMap["R"] = getCriticalPath();
			}
			
			else { // wInY=2, one addition
				setCriticalPath( getMaxInputDelays( inputDelays ));
				// No need for the following, the mult should be absorbed in the addition
				// manageCriticalPath( target->localWireDelay() + target->lutDelay() );
				vhdl << tab << declare("R0",wInX+2) << " <= (";
				if (sign) 
					vhdl << "XX" << of(wInX-1) << " & XX" << of(wInX-1);  
				else  
					vhdl << "\"00\"";
				vhdl <<  " & XX) when YY(0)='1' else "<<zg(wInX+2,0)<<";"<<endl;	
				vhdl << tab << declare("R1i",wInX+2) << " <= ";
				if (sign) 
					vhdl << "(XX" << of(wInX-1) << "  & XX & \"0\")";
				else  
					vhdl << "(\"0\" & XX & \"0\")";
				vhdl << " when YY(1)='1' else " << zg(wInX+2,0) << ";"<<endl;	
				vhdl << tab << declare("R1",wInX+2) << " <= ";
				if (sign) 
					vhdl << "not R1i;" <<endl;
				else  
					vhdl <<  "R1i;" <<endl;
				
				IntAdder *resultAdder = new IntAdder( target, wInX+2, inDelayMap("X", target->localWireDelay() + getCriticalPath() ) );
				oplist.push_back(resultAdder);
				
				inPortMap(resultAdder, "X", "R0");
				inPortMap(resultAdder, "Y", "R1");
				inPortMapCst(resultAdder, "Cin", (sign? "'1'" : "'0'"));
				outPortMap( resultAdder, "R", "RAdder");
				vhdl << tab << instance(resultAdder, "ResultAdder") << endl;
				syncCycleFromSignal("RAdder");
				setCriticalPath( resultAdder->getOutputDelay("R"));
				
				vhdl << tab << "R <= RAdder;"<<endl;	
				outDelayMap["R"] = getCriticalPath();
			}
		} 

				
#if 1 // The Romanian version: shorter latency (when unpipelined), larger area
		else { //generic case
			int chunkSize_ = target->lutInputs()/2;
			int chunksX =  int(ceil( ( double(wInX_) / (double) chunkSize_) ));
			int chunksY =  int(ceil( ( double(wInY_) / (double) chunkSize_) ));
			
			setCriticalPath( getMaxInputDelays(inputDelays) );

			REPORT(DEBUG, "X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
			if (chunksX + chunksY > 2) { //we do more than 1 subproduct
				int widthX = wInX_;
				int widthY = wInY_;	
				bool swap = false;

				if (chunksX > chunksY){ //interchange X with Y
					int tmp = chunksX;
					chunksX = chunksY;
					chunksY = tmp;
	
					tmp = widthX;
					widthX = widthY;
					widthY = tmp;
					swap = true;
				}		

				vhdl<<tab<<declare("sX",chunkSize_*chunksX)<<" <= "<< (swap? "Y":"X")<<" & "<<zg( chunkSize_ * chunksX -widthX,0)<< ";"<<endl;
				vhdl<<tab<<declare("sY",chunkSize_*chunksY)<<" <= "<< (swap? "X":"Y")<<" & "<<zg(chunkSize_*chunksY-widthY,0)<< ";"<<endl;	
				//SPLITTINGS
				for (int k=0; k<chunksX ; k++)
					vhdl<<tab<<declare(join("x",k),chunkSize_)<<" <= sX"<<range((k+1)*chunkSize_-1,k*chunkSize_)<<";"<<endl;
				for (int k=0; k<chunksY ; k++)
					vhdl<<tab<<declare(join("y",k),chunkSize_)<<" <= sY"<<range((k+1)*chunkSize_-1,k*chunkSize_)<<";"<<endl;

				manageCriticalPath( target->lutDelay() + target->localWireDelay(target->lutInputs()) );
		
				if (sign_){
					/* ----------------------------------------------------*/
					//COMPUTE PARTIAL PRODUCTS
					for (int i=0; i<chunksY; i++)
						for (int j=0; j<chunksX; j++){
							if ((i<chunksY-1) && (j<chunksX-1)){
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_+2)<<" <= ( \"0\" & " << join("x",j) <<") * " 
								    << "( \"0\" & " << join("y",i)<< ");" << endl;
							} else if ((i==chunksY-1) && (j<chunksX-1)){
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_+1)<<" <= ( \"0\" & " << join("x",j) <<") * " 
								    << "(" << join("y",i)<< ");" << endl;
							} else if ((i<chunksY-1) && (j==chunksX-1)){	
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_+1)<<" <= ( " << join("x",j) <<") * " 
								    << "( \"0\" & " << join("y",i)<< ");" << endl;
							} else{
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_)<<" <= "<< join("x",j)<<" * " << join("y",i)<< ";" << endl;
							}								
						}		

					map<string, int> availableZeros;
					map<string, int> neededSignExtensions;

					int adderWidth = (chunksX+chunksY)*chunkSize_;
					// CONCATENATE PARTIAL PRODUCTS
					for (int i=0; i<2; i++){
						for (int j=0; j<chunksX; j++){
							int startIdx = chunksY-1-i;
							int paddWidth = adderWidth - chunkSize_*2*(startIdx/2+1);
							int endPaddWidth = chunkSize_*(j+startIdx%2);

							if (j!=chunksX-1){ //the last chunk makes us problems 
								if (i!=0)
									availableZeros[join("cp",i,j)] = paddWidth-endPaddWidth; 
							}else{
								//each of these subproducts requires sign extension
								for (int k=startIdx; k>=0; k-=2)
									if (k!=startIdx)
										neededSignExtensions[join("px",j,"y",k)]= adderWidth - (j+k+2)*chunkSize_;
							}
						}
					}

					map<string,int>::iterator it, it2;

					REPORT(DEBUG, "available ");						
					for (it = availableZeros.begin(); it!=availableZeros.end(); it++)
						REPORT(DEBUG, it->first << " -> " << it->second);

					REPORT(DEBUG, "needed");
					for (it = neededSignExtensions.begin(); it!=neededSignExtensions.end(); it++)
						REPORT(DEBUG, it->first << " -> " << it->second);
			
					map< string, pair<string,int> > reallocations;
					map< string, int> unallocated; 
			
					for (it=neededSignExtensions.begin(); it!= neededSignExtensions.end(); it++){
						bool reallocationsStatus = false;
						REPORT(DEBUG, "trying to allocate " << it->first << "," << it->second);
						for (it2=availableZeros.begin(); it2!= availableZeros.end() && (!reallocationsStatus); it2++){
							if (it2->second >= it->second){ //we have enough zeros to allocate
								reallocations[ it2->first ] = pair<string,int>(it->first, it->second); 
								reallocationsStatus = true;
								availableZeros.erase(it2->first); //delete the entry.
							}
						}
						if (!reallocationsStatus) // we could not realocate
							unallocated[it->first] = it->second;
					}
			
					map< string, pair<string,int> >::iterator it3;
			
					REPORT(DEBUG, "reallocated ");
					for (it3 = reallocations.begin(); it3!=reallocations.end(); it3++)
						REPORT(DEBUG, it3->first << " --> " << it3->second.first << ", " << it3->second.second);

					REPORT(DEBUG, "unallocated");
					for (it = unallocated.begin(); it!=unallocated.end(); it++)
						REPORT(DEBUG, it->first << " -> " << it->second);
			
					vector<string> operandNames;
			
					for (int i=0; i<2; i++){
						for (int j=0; j<chunksX; j++){
							int startIdx = chunksY-1-i;
							int paddWidth = adderWidth - chunkSize_*2*(startIdx/2+1);
							int endPaddWidth = chunkSize_*(j+startIdx%2);

							operandNames.push_back(join("cp",i,j));
							vhdl << tab << declare(join("cp",i,j),adderWidth) << " <= ";

							if (j!=chunksX-1){ //the last chunk makes us problems 
								if (i==0){ //multiply with the last chunk of y, so sign extend
									vhdl << rangeAssign(paddWidth-endPaddWidth-1, 0, join("px",j,"y",startIdx)+of(2*chunkSize_-1));
								}else{
									//we check the reallocations. Some zeros might be reallocated
									it3 = reallocations.find(join("cp",i,j));

									if (it3!=reallocations.end()){ // some reallocation happended here
										REPORT(DEBUG, "reallocated " << it3->second.second << " sign bits over " << paddWidth-endPaddWidth << " zero bits");
										vhdl << rangeAssign( it3->second.second-1, 0, it3->second.first+of(2*chunkSize_-1)) << " & " << zg(paddWidth-endPaddWidth-it3->second.second,0) ; 	
									}else{
										vhdl << zg(paddWidth-endPaddWidth,0) ; 	
									}
								}
								for (int k=startIdx; k>=0; k-=2)
									vhdl << " & " << use(join("px",j,"y",k))<<range(2*chunkSize_-1,0);
								vhdl << " & " << zg((startIdx<0?2:endPaddWidth),0) << ";" << endl;
							}else{
								//the ms one does't need any reallocations
								vhdl << rangeAssign(paddWidth-endPaddWidth-1, 0, join("px",j,"y",startIdx)+of(2*chunkSize_-1));
								for (int k=startIdx; k>=0; k-=2){
									if (k==startIdx){
										vhdl << " & " << join("px",j,"y",k)<<range(2*chunkSize_-1,0);
									}else{
										//check if the sign extension was reallocated
										bool reallocationsStatus = false;
										for(it3 = reallocations.begin(); it3!=reallocations.end() && (!reallocationsStatus); it3++){
											if (it3->second.first == join("px",j,"y",k)){
												reallocationsStatus = true;
											}
										}
								
										if (reallocationsStatus){ //if the sign extension was reallocated we simply continue
											vhdl << " & " << join("px",j,"y",k)<<range(2*chunkSize_-1,0);
										}else{
											//not reallocated.
											//we need to finish the current cp and start a new one
											vhdl << " & " << zg( (j+1 + k+1)*chunkSize_, 0) << ";" << endl;
											operandNames.push_back(join("cp",i,j,k));
											vhdl << tab << declare(join("cp",i,j,k),adderWidth) << " <= ";
											vhdl << rangeAssign(adderWidth - (j+k+2)*chunkSize_ -1 , 0, join("px",j,"y",k)+of(2*chunkSize_-1))
											     << " & " << join("px",j,"y",k) << range(2*chunkSize_-1,0);
										}											
									}
									if (k-2 <0) //we need to finish
										vhdl << " & " << zg( (j + k)*chunkSize_, 0) << ";--" << endl;
								}
							}
						}
					}

					IntMultiAdder* add =  new IntMultiAdder(target, adderWidth, operandNames.size(), inDelayMap("X0", getCriticalPath()));
					oplist.push_back(add);

					for (int i=0; i< signed(operandNames.size()); i++)
						inPortMap (add, join("X",i) , operandNames[i]);

					outPortMap(add, "R", "addRes");
					vhdl << instance(add, "adder");
					syncCycleFromSignal("addRes");
					setCriticalPath(add->getOutputDelay("R"));

					outDelayMap["R"] = getCriticalPath();	
					vhdl<<tab<<"R<=addRes" << range(adderWidth-1,adderWidth-wInX_-wInY_) << ";" << endl;		
				}else{ //unsigned version
					/* ----------------------------------------------------*/
					//COMPUTE PARTIAL PRODUCTS
					for (int i=0; i<chunksY; i++)
						for (int j=0; j<chunksX; j++)
							vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_)<<" <= "<< join("x",j)<< " * " << join("y",i) << ";" << endl;
				
					int adderWidth = (chunksX+chunksY)*chunkSize_;
					// CONCATENATE PARTIAL PRODUCTS
					for (int i=0; i<2; i++){
						for (int j=0; j<chunksX; j++){
							vhdl << tab << declare(join("cp",i,j),adderWidth) << " <= ";

							int startIdx = chunksY-1-i;
							int paddWidth = adderWidth - chunkSize_*2*(startIdx/2+1);
							int endPaddWidth = chunkSize_*(j+startIdx%2);
							vhdl << zg(paddWidth-endPaddWidth,0); 

							for (int k=startIdx; k>=0; k-=2)
								vhdl << " & " << use(join("px",j,"y",k));
			
							vhdl << " & " << zg((startIdx<0?2:endPaddWidth),0) << ";" << endl;
						}
					}

					IntMultiAdder* add =  new IntMultiAdder(target, adderWidth, 2*chunksX, inDelayMap("X0", getCriticalPath()));
					oplist.push_back(add);

					for (int i=0; i<2; i++)
						for (int j=0; j<chunksX; j++)
							inPortMap (add, join("X",i*chunksX+j) , join("cp",i,j));

					outPortMap(add, "R", "addRes");
					vhdl << instance(add, "adder");
					syncCycleFromSignal("addRes");
					setCriticalPath(add->getOutputDelay("R"));
					vhdl<<tab<<"R<=addRes" << range(adderWidth-1,adderWidth-wInX-wInY) << ";" << endl;		
					outDelayMap["R"] = getCriticalPath();				
				}
			}
		}
			
#else /////////////////// Here begins Florent Version/////////////////////
		else { //generic case
			if(sign_) {
				cout << "Not implemented" << endl;
				exit(1);
			}
			else{ 
				int dx, dy0, dy;				// Here we need to split in small sub-multiplications
				int li=target->lutInputs();
				
				dx = li>>1;
				dy0 = li - dx; 
				// TODO for Altera remove the -1
#if 1 // If the synthesis tool can merge the sub-product and the LUT: it seems the case
				dy = dy0 - 1;
#else
				dy=dy0;
#endif
				int hole =  (dx==dy?0:1);
				int hole0 = (dx==dy0?0:1);
				REPORT(DEBUG, "dx="<< dx << "  dy=" << dy << " dy0=" << dy0);

				int chunksX =  int(ceil( ( double(wInX_) / (double) dx) ));
				int chunksY =  int(ceil( ( 1+ double(wInY_-dy0) / (double) dy) ));
				int sizeXPadded=chunksX*dx; 
				int sizeYPadded=dy*(chunksY-1) + dy0; 
					// TODO remove one MSB  in case of hole
				int sizeAddendA = (chunksX+1)/2 *li;
				int sizeAddendB = chunksX/2*li;
				
				setCriticalPath( getMaxInputDelays(inputDelays) );

				REPORT(DEBUG, "X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
				REPORT(DEBUG, " sizeXPadded="<< sizeXPadded << "  sizeYPadded="<< sizeYPadded);
				REPORT(DEBUG, " sizeAddendA=" << sizeAddendA << "  sizeAddendB=" << sizeAddendB);
				if (chunksX + chunksY > 2) { //we do more than 1 subproduct
					vhdl << tab << "-- padding to the right" << endl;
					int sX=dx*chunksX;
					int sY=dy*(chunksY-1)+dy0;
					vhdl<<tab<<declare("Xp",sX)<<" <= XX & "<< zg(dx * chunksX - wInX, 0)<< ";"<<endl;
					vhdl<<tab<<declare("Yp",sY)<<" <= YY & "<< zg(dy*(chunksY-1)+dy0 - wInY, 0)<< ";"<<endl;	
					//SPLITTINGS
					for (int k=0; k<chunksX ; k++)
						vhdl<<tab<<declare(join("x",k),dx)<<" <= Xp"<<range((k+1)*dx-1,k*dx)<<";"<<endl;
					vhdl<<tab<<declare("y0",dy0)<<" <= Yp"<<range(dy0-1, 0)<<";"<<endl;
					for (int k=1; k<chunksY ; k++)
						vhdl<<tab<<declare(join("y",k),dy)<<" <= Yp"<<range(dy0+k*dy-1, dy0+(k-1)*dy)<<";"<<endl;
					
					// The first row
					
					manageCriticalPath(target->localWireDelay(chunksX) + target->lutDelay());

					SmallMultTable *t0 = new SmallMultTable( target, dx, dy0 );
					oplist.push_back(t0);
					
					for (int i=0; i<chunksX; i++) { 
						vhdl<<tab<<declare(XY(i,0), dx+dy0) << " <= y0 & x" << i << ";"<<endl;
						inPortMap(t0, "X", XY(i,0));
						outPortMap(t0, "Y", PP(i,0));
						vhdl << instance(t0, PPTbl(i,0));		
					}
					vhdl << tab << "-- initializing the accumulators" << endl;
					//			vhdl << tab << declare("R0", dx) << " <= " << PP(0,0) << range(dx-1, 0) << ";  -- go directly to output" <<endl;
					vhdl << tab << declare("SumA0", sizeAddendA) << " <= ";
					for(int j=(chunksX+1)/2-1; j>0; j--) {
						if (dy0<dx) 
							vhdl << " \"0\" & ";
						vhdl << PP(2*j,0) << " & ";
					}
					vhdl << PP(0,0) << range(dx+dy0-1, 0) << ";" <<endl; 
					vhdl << tab << declare("SumB0", sizeAddendB) << " <= ";
					for(int j=chunksX/2-1; j>=0; j--) {
						if (dy0<dx) 
							vhdl << " \"0\" & ";
						vhdl << PP(2*j+1, 0) ;
						if(j!=0)  vhdl << " & ";
					}
					vhdl << ";" <<endl; 

					// The following rows
					
					SmallMultTable *t = new SmallMultTable( target, dx, dy );
					oplist.push_back(t);

					// Estimation of the critical path
					// We don't use manageCriticalPath because the additions are cascaded

					bool pipe[1000]; // pipe[i] == pipe before row i+1

					double cp=getCriticalPath();
					double currentFanoutCP = target->localWireDelay(chunksX);
					int rowsInStage=1;
					REPORT(DEBUG, "");

					for (int iy=1; iy<chunksY; iy++){
						rowsInStage ++;
						cp -=  currentFanoutCP; 
						currentFanoutCP = target->localWireDelay((rowsInStage)*chunksX);
						cp += currentFanoutCP;
						if(iy==1) // first adder entails along carry propagation
						 cp += target->adderDelay(sizeAddendA);
						else
							cp += target->adderDelay(li); // one LUT and just li bits of carry propagation
						// following line cut from manageCriticalPath()
						if ( target->ffDelay() + cp > (1.0/target->frequency())){
							pipe[iy] = true;
							rowsInStage = 1;
							currentFanoutCP = target->localWireDelay(chunksX);
							cp = currentFanoutCP + target->adderDelay(sizeAddendA);
						}
						else pipe[iy] = false;
					}
					// Final addition
					rowsInStage ++;
					cp -=  currentFanoutCP; 
					currentFanoutCP = target->localWireDelay((rowsInStage)*chunksX);
					cp += currentFanoutCP;
					cp += target->adderDelay(sizeAddendA); 
					if ( target->ffDelay() + cp > (1.0/target->frequency()))
						pipe[chunksY] = true;
					else
						pipe[chunksY] = false;

					
					// Now actual code generation
					for (int iy=1; iy<chunksY; iy++){
						// Not easy to understand how the addition and the table are merged on virtex5
						// The largest of the two additions is always on sumA 

						if(pipe[iy]) nextCycle();

						vhdl << tab << "-- Partial product row number " << iy << endl;
						for (int ix=0; ix<chunksX; ix++) { 
							vhdl<<tab<<declare(XY(ix,iy), dx+dy) << " <= y" << iy << " & x" << ix << ";"<<endl;
							inPortMap(t, "X", XY(ix,iy));
							outPortMap(t, "Y", PP(ix,iy));
							vhdl << instance(t, PPTbl(ix,iy));		
						}
						vhdl << tab << "-- Building addend for addition << "<< iy-1 << endl;
						vhdl << tab << declare(join("addendA",iy), sizeAddendA) << " <= ";
						for(int ix=(chunksX+1)/2-1; ix>=0; ix--) {
							if (dy<dx) 
								vhdl << "\"0\" & ";
							vhdl << PP(2*ix, iy);
							if (ix>0) 
								vhdl << " & ";
						}
						vhdl << ";" <<endl; 
						vhdl << tab << declare(join("addendB",iy), sizeAddendB) << " <= ";
						for(int ix=chunksX/2-1; ix>=0; ix--) {
							if (dy<dx) 
								vhdl << "\"0\" & ";
							vhdl << PP(2*ix+1,iy);
							if(ix>0)  
								vhdl << " & ";
						}
						vhdl << ";" <<endl; 
						vhdl << tab << "-- The additions" << endl;
						int shift;
						if(iy==1) 
							shift=dy0; 
						else 
							shift=dy;
						vhdl << tab << declare(join("SumAl",iy), shift) << " <= " << join("SumA",iy-1) << range(shift-1, 0) << ";" <<endl;
						vhdl << tab << declare(join("SumA",iy), sizeAddendA) << " <= " << join("SumA",iy-1) << range(sizeAddendA-1, shift) << " + " << join("addendA",iy) << ";" <<endl;
						vhdl << tab << declare(join("SumBl",iy), shift) << " <= " << join("SumB",iy-1) << range(shift-1, 0) << ";" <<endl;
						vhdl << tab << declare(join("SumB",iy), sizeAddendB) << " <= " << join("SumB",iy-1) << range(sizeAddendB-1, shift) << " + " << join("addendB",iy) << ";" <<endl;
					}

					if(pipe[chunksY]) nextCycle();

					manageCriticalPath(target->localWireDelay() + target->adderDelay(sizeXPadded+sizeYPadded-dx));

					vhdl << tab << "-- Final sum " << endl;
					vhdl << tab << declare("SumAfinal", sizeAddendA+(chunksY-2)*dy+dy0) << " <= " << join("SumA", chunksY-1); 
					for (int iy=chunksY-1; iy>0; iy--){
						vhdl << " & " << join("SumAl",iy);
					}
					vhdl << ";" <<endl; 
					vhdl << tab << declare("SumBfinal", sizeAddendB+(chunksY-2)*dy+dy0) << " <= " << join("SumB", chunksY-1); 
					for (int iy=chunksY-1; iy>0; iy--){
						vhdl << " & " << join("SumBl",iy);
					}
					vhdl << ";" <<endl; 
					// 					vhdl << tab << declare("Ph", chunksX*dx) << " <= " << join("SumA",chunksY-1) << " + " << join("SumB",chunksY-1) << ";" <<endl;
					vhdl << tab << declare("P", sizeXPadded+sizeYPadded+dy0-dy) << " <= SumAfinal + (SumBfinal & "<< zg(dx) <<");" << endl; 
					vhdl << tab << "R <= " << "P" << range(sizeXPadded+sizeYPadded-1, sizeXPadded+sizeYPadded-wInX-wInY) << ";" << endl; 
				}
			}
			
		}	 
#endif 
		
	}

	/*
Example dx=3 dy=3

row 0
                                        xxxxxx
                                     xxxxxx  
                                  xxxxxx
                               xxxxxx
                           xxxxxx   


 ......................    xxxxxxxxxxxxxxxxxxx  sumA0
................              xxxxxxxxxxxxx     sumB


 ......................    xxxxxxxxxxxxxxxxxxx  sumA0
...................     xxxxxxxxxxxxxxxxxxx     suma1
................           xxxxxxxxxxxxx   



*/
	
			
	LogicIntMultiplier::~LogicIntMultiplier() {
	}

#if 0 // Xilinx specific
	void LogicIntMultiplier::outputVHDL(std::ostream& o, std::string name) {
		licence(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.std_logic_arith.all;" << endl;
		if (sign_)
			o << "use ieee.std_logic_signed.all;" << endl;
		else
			o << "use ieee.std_logic_unsigned.all;" << endl;
		o << "library work;" << endl;
		outputVHDLEntity(o);
		newArchitecture(o,name);
		o << "	attribute multstyle : string;"<<endl;
		o << "	attribute multstyle of arch : architecture is \"logic\";"<<endl;
		o << buildVHDLComponentDeclarations();	
		o << buildVHDLSignalDeclarations();
		beginArchitecture(o);		
		o<<buildVHDLRegisters();
		o << vhdl.str();
		endArchitecture(o);
	}
#endif


	void LogicIntMultiplier::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
		if (! sign_){

			mpz_class svR = svX * svY;

			tc->addExpectedOutput("R", svR);
		}else{
			mpz_class big1 = (mpz_class(1) << (wInX_));
			mpz_class big1P = (mpz_class(1) << (wInX_-1));
			mpz_class big2 = (mpz_class(1) << (wInY_));
			mpz_class big2P = (mpz_class(1) << (wInY_-1));

			if ( svX >= big1P)
				svX = svX-big1;

			if ( svY >= big2P)
				svY = svY -big2;
			
			mpz_class svR = svX * svY;
			if ( svR < 0){
				mpz_class tmpSUB = (mpz_class(1) << (wOut_));
				svR = tmpSUB + svR; 
			}

			tc->addExpectedOutput("R", svR);
		}

	}
}
