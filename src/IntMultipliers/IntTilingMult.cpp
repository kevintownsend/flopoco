/*
  Tilling Multiplier for FloPoCo
 
  Authors:  Sebastian Banescu, Radu Tudoran, Bogdan Pasca
 
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <string.h>
#include <limits.h>

#include <gmp.h>


#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "../IntMultiplier.hpp"
#include "IntTilingMult.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <math.h>
#include <locale>
#include <cfloat>

#include <stdio.h>
#include <mpfr.h>

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

#define DEBUGVHDL 0

	
	IntTilingMult:: IntTilingMult(Target* target, int wInX, int wInY, float ratio, int maxTimeInMinutes, bool interactive, map<string, double> inputDelays) :
		Operator(target), wInX(wInX), wInY(wInY), wOut(wInX + wInY),ratio(ratio), maxTimeInMinutes(maxTimeInMinutes-1){
 
		ostringstream name;
			srcFileName="IntTilingMultiplier";
			
			start = clock(); /* time management */
			
			name <<"IntMultiplier_"<<wInX<<"_"<<wInY<<"_uid"<<getNewUId();
			setName(name.str());
		
			setCopyrightString("Sebastian Banescu, Radu Tudoran, Bogdan Pasca 2009-2011");
	
			addInput ("X", wInX);
			addInput ("Y", wInY);
			addOutput("R", wOut); /* wOut = wInX + wInY */
		
			warningInfo();
		
			/* the board is extended horizontally and vertically. This is done
			in order to fit DSPs on the edges. */
			vn = wInX + 2*getExtraWidth(); 
			vm = wInY + 2*getExtraHeight();
			
			/* one of the important optimizations consist in removing the extra
			stripes on the top and on the bottom */	
			vnme = vn - getExtraWidth();		
			vmme = vm - getExtraHeight();
			
			
			REPORT( INFO, "Input board size is width="<<wInX<<" height="<<wInY);
			REPORT( INFO, "Extra width:="<<getExtraWidth()<<" extra height:="<<getExtraHeight());
			REPORT( INFO, "Extended board width="<<vnme<<" height="<<vmme);

			/* detailed info about the abstracted DSP block */
			int x,y;
			target->getDSPWidths(x,y);
			REPORT( INFO, "DSP block: width= "<< x <<" height=" << y);

			/* get an estimated number of DSPs needed to tile the board with */
			nrDSPs = estimateDSPs();
			REPORT( INFO, "Estimated DSPs = " <<nrDSPs);

			/* the maximum number of DSPs that can be chained on Virtex devices.
			The number is related with the number of register levels inside the
			DSPs. */
			if ( ( target_->getID() == "Virtex4") ||
			 ( target_->getID() == "Spartan3"))  // then the target is A Xilinx FPGA 
				nrOfShifts4Virtex=5;//int(sqrt(nrDSPs));			
			else
				nrOfShifts4Virtex=20;
			
			//~ float tempDist =	 (movePercentage  * getExtraWidth() * getExtraWidth()) /4.0 + (movePercentage *getExtraHeight() * getExtraHeight()) /4.0;
			/* each time a tile is placed maxDist2Move is the maximum distance
			it may be placed from the already placed blocks */
			float tempDist =	0;
			maxDist2Move = (int) ( sqrt(tempDist) );
			REPORT( INFO, "maxDist2Move = "<<maxDist2Move);
		
			float const scale=100.0;
//			costDSP = ( (1.0+scale) - scale * ratio );

			costDSP = 1.0;
			REPORT(INFO, "Cost DSP is " << costDSP);
//			costLUT = ( (1.0+scale) - scale * (1-ratio) ) / ((float) target_->getEquivalenceSliceDSP() );

			costLUT = 1.0 / double(target_->getEquivalenceSliceDSP());
			REPORT(INFO, "Equivalent slices to implement one DSP=" << target_->getEquivalenceSliceDSP());
			REPORT(INFO, "Cost LUT is " << costLUT);
			
		
			//the one
			
			if (interactive){
				cout << " DO you want to run the algorithm? (y/n)" << endl;
				string myc;
				cin >> myc;
				if ( myc.compare("y")!=0)
					exit(-1);
			}						
			runAlgorithm();		
			
		
			REPORT(INFO, "Estimated DSPs:= "<<nrDSPs);
			target->getDSPWidths(x,y);
			REPORT(INFO,"Width of DSP is := "<<x<<" Height of DSP is:="<<y);
			REPORT(INFO,"Extra width:= "<<getExtraWidth()<<" \nExtra height:="<<getExtraHeight());
			REPORT(INFO,"maxDist2Move:= "<<maxDist2Move);
	
		
			generateVHDLCode4CompleteTilling();

		
	}
	
	IntTilingMult::~IntTilingMult() {
	}



	
	void IntTilingMult::generateVHDLCode4CompleteTilling(){
		
//		bestConfig = splitLargeBlocks(bestConfig, nrDSPs);
		bindDSPs(bestConfig);      
		int nrDSPOperands = multiplicationInDSPs(bestConfig);
		int nrSliceOperands = multiplicationInSlices(bestConfig);
		REPORT(DETAILED, "nr of DSP operands="<<nrDSPOperands);
		REPORT(DETAILED, "nr of softDSP operands="<<nrSliceOperands);
		map<string, double> inMap;
		
		for (int j=0; j<nrDSPOperands; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOpDSP" << j;
				syncCycleFromSignal(concatPartialProd.str());
			}	
			
		for (int j=0; j<nrSliceOperands; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOpSlice" << j;
				syncCycleFromSignal(concatPartialProd.str());
			}
		nextCycle();/////////////////////////////////////////////////////////	

		if (nrDSPOperands+nrSliceOperands>1){
			Operator* add;
//			if ( ( target_->getID() == "Virtex4") ||
//				 ( target_->getID() == "Spartan3"))  // then the target is A Xilinx FPGA 
//				add =  new IntNAdder(getTarget(), wInX+wInY, nrDSPOperands+nrSliceOperands, inMap);
//			else
				add =  new IntCompressorTree(getTarget(), wInX+wInY, nrDSPOperands+nrSliceOperands,inMap);
			REPORT(DEBUG, "Finished generating the compressor tree");
		
			oplist.push_back(add);


			for (int j=0; j<nrDSPOperands; j++)
				{
					ostringstream concatPartialProd;
					concatPartialProd  << "addOpDSP" << j;

					inPortMap (add, join("X",j) , concatPartialProd.str());
				}	
			REPORT(DEBUG, "Finished mapping the DSP Operands");
					
			for (int j=0; j<nrSliceOperands; j++)
				{
					ostringstream concatPartialProd;
					concatPartialProd  << "addOpSlice" << j;
					REPORT(DETAILED, "@ In Port Map Current Cycle is " << getCurrentCycle());
					inPortMap (add, join("X",j+nrDSPOperands) , concatPartialProd.str());
				}	
			REPORT(DEBUG, "Finished mapping the SLICE Operands");


//			if ( ( target_->getID() == "Virtex4") ||
//				 ( target_->getID() == "Spartan3"))  // then the target is A Xilinx FPGA 
//				inPortMapCst(add, "Cin", "'0'");
			outPortMap(add, "R", "addRes");
			vhdl << instance(add, "adder");

			syncCycleFromSignal("addRes");
			setCriticalPath( add->getOutputDelay("R") );
			outDelayMap["R"] = getCriticalPath();
			vhdl << tab << "R <= addRes;" << endl;
		}else{
			if (nrDSPOperands==1){
				syncCycleFromSignal("addOpDSP0");
				vhdl << tab << "R <= addOpDSP0;" << endl;
			}else{
				syncCycleFromSignal("addOpSlice0");
				vhdl << tab << "R <= addOpSlice0;" << endl;
			}
		}
			
	}
	
	
	

	void IntTilingMult::runAlgorithm()
	{
		
		int n = vn, m = vm;
	
		/* declare and initialize a matrix of integers for printing purposes */
		mat = new int*[m];
		for(int i=0;i<m;i++){
			mat[i] = new int [n];
			for(int j=0;j<n;j++)
				mat[i][j]=0;
		}
		
		/* if the estimated number of DSPs is grather then 0 then we should 
		run the algorithm */
		if(nrDSPs>0){
			tempc = new DSP*[nrDSPs];
			for(int ii=0;ii<nrDSPs;ii++)
				tempc[ii]= new DSP();

			/* declare an array to store if multipliers are rotated. Say 
			a multiplier in Virtex5 returns x=17 y=24. If it is used with
			x=24 and y=17 the corresponding position in the rotate table 
			will be true */
			rot = new bool[nrDSPs];
			for(int i =0;i<nrDSPs;i++)
				rot[i]=false; /* initially none are rotated */
		
			//The second
			numberDSP4Overlap = nrDSPs;
			initTiling (globalConfig,nrDSPs);
			display(globalConfig);
			
			//this will initialize the bestConfig with the first configuration
			bestCost = FLT_MAX ;
			REPORT(INFO, "Max score is" << bestCost);
			bestConfig = new DSP*[nrDSPs];
			for(int i=0;i<nrDSPs;i++)
				bestConfig[i]= new DSP();
			
			compareCost();

			//the one
			numberDSP4Overlap = nrDSPs;

			tilingAlgorithm (nrDSPs-1, nrDSPs-1, false, nrDSPs-1);
			
			bindDSPs(bestConfig);
				
			display(bestConfig);
			REPORT(INFO, "Best cost is "<<bestCost);
		}else{
			n=vn;
			m=vm;
	
			mat = new int*[m];
			for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
			tempc= new DSP*[0];
			bestConfig = new DSP*[1];
			globalConfig = new DSP*[1];
			tempc[0]=bestConfig[0]=globalConfig[0]=NULL;
			
			bestCost = FLT_MAX ;
			REPORT(INFO, "Max score is "<<bestCost);
			compareCost();
			REPORT(INFO, "New best score is "<<bestCost);
			display(bestConfig);
		}
	
	}

	
	/** The movement of the DSP blocks with values belonging to their widths 
	and heights still needs to be done. Now it only runs with one type of move 
	on one of the directions, which is not ok for the cases when the DSPs 
	are not squares.
	 */
	void IntTilingMult::tilingAlgorithm(int i, int n, bool repl, int lastMovedDSP)
	{
		finish = (clock() - start)/(CLOCKS_PER_SEC*60);
		if (finish > maxTimeInMinutes){
			cout << "Time's up!"<<endl;
			return;
		}
				
		if (i == n){
			if (repl == true){ // if previous DSPs were moved this one needs to recompute all positions 
				if(replace(globalConfig,i)){ // repostioned the DSP
					compareCost();
					rot[i]=false;
					tilingAlgorithm(i,n,false,lastMovedDSP);
				}else{ 
				/* could not reposition the DSP in the bounds of the tiling board */
					rot[i]=false;
					if( lastMovedDSP>=0) // go one level up the backtracking stack
					{	
						tilingAlgorithm(lastMovedDSP,n,false, lastMovedDSP);
					}
				}
			}else{ 
				/* the last DSP is being moved on the tiling board */
				/* if successfuly moved the last block */
				if( move(globalConfig,i)){ 
					compareCost();
					finish = (clock() - start)/(CLOCKS_PER_SEC*60);

					
					tilingAlgorithm(i,n,repl,i);		//repl should be false
				} else {
					/* if we didn't rotate it and it is worth rotating */
					if (rot[i] == false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() )){ 
						globalConfig[i]->rotate();
						rot[i]=true;
						/* try repositioning the rotated block */
						if(replace(globalConfig,i)){
							compareCost();
							finish = (clock() - start)/(CLOCKS_PER_SEC*60);
							tilingAlgorithm(i, n, repl, i);		//repl should be false
						} else{
						 /* go to the previous block */
							if( i-1 >= 0){	
								finish = (clock() - start)/(CLOCKS_PER_SEC*60);
								tilingAlgorithm(i-1,n,false,i);
							}
						}
					}else {
						/* the DSP was either rotated already or is square */
						if(i-1 >= 0){	
							finish = (clock() - start)/(CLOCKS_PER_SEC*60);
							tilingAlgorithm(i-1,n,repl,i);		//repl should be false
						}
					}
				}
			}
		} else { 
			/* we are not at the last DSP */
			if(repl == true){ 
				/* the previuos DSPs were successfuly repositioned */
				if (replace(globalConfig,i)){ 
					/* and now we reposition the current DSP */
					rot[i]=false;
					tilingAlgorithm(i+1, n, repl, lastMovedDSP);
				}else{ 
					/* the current DSP could not be repositioned */
					/* go to the DSP block that was moved (not repostioned) 
					the last time */
					rot[i]=false;
					if( lastMovedDSP>=0) {
						finish = (clock() - start)/(CLOCKS_PER_SEC*60);
						tilingAlgorithm( lastMovedDSP,n,false, lastMovedDSP);
					}
				}
			}else{ 
			/* the folling DSP could not be moved or repositioned */
				if(move(globalConfig,i)){ 
				/* the current DSP was successfuly moved*/
					finish = (clock() - start)/(CLOCKS_PER_SEC*60);
					tilingAlgorithm(i+1,n,true,i);
				}else{ 
				/* the current DSP was not moved successfuly */
					if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() )){
						/* if the DSP was not rotated and is not sqare then roteate it */
						globalConfig[i]->rotate();
						if(replace(globalConfig,i)){ // the current DSP was successfuly repositioned
							rot[i]=true;
							tilingAlgorithm(i+1,n,true,i);
						}else{ 
						/* the current DSP was not successfuly repositioned */
							if(i-1>=0){
								tilingAlgorithm(i-1,n,repl,i);
							}
						}
					}else{ 
						/* the DSP is either square or has been already rotated */
						if(i-1>=0){
							tilingAlgorithm(i-1,n,repl,i);		//repl should be false
						}
					}
				}
			}
		}
	}


	bool IntTilingMult::compareOccupation(DSP** config)
	{
		int totalSize = (vnme-getExtraWidth()) *(vmme-getExtraHeight());
		int sizeBest = totalSize;
		int sizeConfig = totalSize;
		//~ cout<<"Total size is "<<totalSize<<endl;
		int c1X,c2X,c1Y,c2Y;
		int nj,ni,njj,nii;
		int ii,jj,i,j;
		int n,m;
		//~ n=wInX + 2* getExtraWidth();
		//~ m=wInY + 2* getExtraHeight();
		n=vm;
		m=vm;
	
		int nmew = vnme;
		int ew = getExtraWidth();
		int mmeh = vmme;
		int eh = getExtraHeight();
	
		for(int ti=0;ti<nrDSPs;ti++)
			if(config[ti]!=NULL)
				{
			
			
					config[ti]->getTopRightCorner(c1X,c1Y);
					config[ti]->getBottomLeftCorner(c2X,c2Y);
					c1X=n-c1X-1;
					c2X=n-c2X-1;
					//fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
					j = c2X;
					i = c1Y;
					jj = c1X;
					ii = c2Y;
			
					//cout<<" New coordinates are ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
			
			
						
					if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
						{
							nj=ni=njj=nii=0;
						}
					else
						{
							if( j < getExtraWidth() )
								nj = getExtraWidth() ;
							else
								nj = j;
							if( jj >= n - getExtraWidth() )
								njj = n -getExtraWidth() -1;
							else
								njj = jj;
							
							if( i < getExtraHeight() )
								ni = getExtraHeight() ;
							else
								ni = i;
							if( ii >= m - getExtraHeight() )
								nii = m -getExtraHeight() -1;
							else
								nii = ii;
							//costSlice +=target->getIntMultiplierCost(njj-nj+1,nii-ni+1);
							sizeConfig-=(njj-nj+1)*(nii-ni+1);
							
						}
			
				}
		
		//cout<<"Size config is"<<sizeConfig<<endl;
		//display(config);
		
		
		
		for(int ti=0;ti<nrDSPs;ti++)
			if(bestConfig[ti]!=NULL)
				{
			
			
					bestConfig[ti]->getTopRightCorner(c1X,c1Y);
					bestConfig[ti]->getBottomLeftCorner(c2X,c2Y);
					c1X=n-c1X-1;
					c2X=n-c2X-1;
					//fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
					j = c2X;
					i = c1Y;
					jj = c1X;
					ii = c2Y;
			
					//cout<<" New coordinates are ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
			
			
						
					if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
						{
							nj=ni=njj=nii=0;
						}
					else
						{
							if( j < getExtraWidth() )
								nj = getExtraWidth() ;
							else
								nj = j;
							if( jj >= n - getExtraWidth() )
								njj = n -getExtraWidth() -1;
							else
								njj = jj;
							
							if( i < getExtraHeight() )
								ni = getExtraHeight() ;
							else
								ni = i;
							if( ii >= m - getExtraHeight() )
								nii = m -getExtraHeight() -1;
							else
								nii = ii;
							//costSlice +=target->getIntMultiplierCost(njj-nj+1,nii-ni+1);
							sizeBest-=(njj-nj+1)*(nii-ni+1);
							
						}
			
				}
		
		
		//cout<<"Size best is "<<sizeBest<<endl;
		//display(bestConfig);
		
		
		if(sizeBest >= sizeConfig)		
			return true;
		else
			return false;
	
	}

	/* fill matrix with DSP index */
	void IntTilingMult::fillMatrix(int **&matrix, int lw, int lh,int topleftX,int topleftY, int botomrightX, int botomrightY, int value){
		for(int j=topleftX; j<= botomrightX; j++){
			for(int i=topleftY; i <= botomrightY; i++){
						if(j>-1 && i>-1 && i<lh && j<lw)
							matrix[i][j]=value;
					}
			}
	}


	void IntTilingMult::display(DSP** config)
	{
		ofstream fig;
		fig.open ("tiling.fig", ios::trunc);
		fig << "#FIG 3.2  Produced by xfig version 3.2.5a" << endl;
		fig << "Landscape" << endl;
		fig << "Center" << endl;
		fig << "Metric" << endl;
		fig << "A4      " << endl;
		fig << "100.00" << endl;
		fig << "Single" << endl;
		fig << "-2" << endl;
		fig << "1200 2" << endl;
	
	
		int **mat;
		int n,m;
		int count=1;
		//~ n=wInX + 2* getExtraWidth();
		//~ m=wInY + 2* getExtraHeight();
		n = vn;
		m= vm;
		REPORT(INFO, "real width"<<vn - 2* getExtraWidth()<<"real height"<<vm - 2* getExtraHeight());
		REPORT(INFO, "width "<<n<<"height "<<m);
		mat = new int*[m];
	
		int nmew = vnme;
		int ew = getExtraWidth();
		int mmeh = vmme;
		int eh = getExtraHeight();
		int nj,ni,njj,nii;
		for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				REPORT(INFO, "DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")");
				fig << " 2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5 " << endl;
				fig << "	  " << (-c2X+getExtraWidth()-1)*45 << " " << (c1Y-getExtraHeight())*45 
				         << " " << (-c1X+getExtraWidth())*45 << " " << (c1Y-getExtraHeight())*45 
				         << " " << (-c1X+getExtraWidth())*45 << " " << (c2Y-getExtraHeight()+1)*45 
				         << " " << (-c2X+getExtraWidth()-1)*45 << " " << (c2Y-getExtraHeight()+1)*45 
				         << " " << (-c2X+getExtraWidth()-1)*45 << " " << (c1Y-getExtraHeight())*45 << endl;
				
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
	
		count++;
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					{
						if(mat[i][j]==0)
							{
								int ver =0;
								int ii=i,jj=j;
								while(ver<6&&(ii<m-1||jj<n-1))
									{
							
							
										if(ver<3)
											{
								
												if(ver==0||ver==1)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=2;							
													}
							
												if(ver==0||ver==2)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=1;
													}
								
								
							
												for(int k=ii,q=jj;k>i-1&&(ver==0||ver==2);k--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=1;
															else
																ver=3;
															jj--;
														}
									
												for(int k=ii,q=jj;q>j-1&&(ver==0||ver==1);q--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=2;
															else
																ver=3;
															ii--;
														}
								
							
											}
										else
											{
												if(ver==3||ver==5)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=4;
													}
								
								
												if(ver==3||ver==4)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=5;							
													}
							
								
															
												for(int k=ii,q=jj;q>j-1&&(ver==3||ver==4);q--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=5;
															else
																ver=6;
															ii--;
														}
								
												for(int k=ii,q=jj;k>i-1&&(ver==3||ver==5);k--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=4;
															else
																ver=6;
															jj--;
														}
												if(ver==5&&jj==n-1)
													ver=6;
												if(ver==4&&ii==m-1)
													ver=6;
										
															
								
											}
									}
						
						
					
						
						
								if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
									{
										REPORT(INFO,"Partition number "<<count<<" is totally out of the real multiplication bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")");
									}
								else
									{
										if( j < getExtraWidth() )
											nj = getExtraWidth() ;
										else
											nj = j;
										if( jj >= n - getExtraWidth() )
											njj = n -getExtraWidth() -1;
										else
											njj = jj;
							
										if( i < getExtraHeight() )
											ni = getExtraHeight() ;
										else
											ni = i;
										if( ii >= m - getExtraHeight() )
											nii = m -getExtraHeight() -1;
										else
											nii = ii;
										REPORT(INFO, "Partition number "<<count<<" with bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<") has now bounds ("<<nj<<" , "<<ni<<" , "<<njj<<" , "<<nii<<")");
									}
						
								REPORT(INFO,j<<" "<<i<<" "<<jj<<" "<<ii);
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
						
							}
					}
			
			}
		
		
		
		char af;
		int afi;
		if (verbose>1){
			for(int i=0;i<m;i++)
			{
				if(i==getExtraHeight())
					cout<<endl;
				
				for(int j=0;j<n;j++)
				{
					if(j==getExtraWidth())
						cout<<" ";
					if(j==n-getExtraWidth())
						cout<<" ";
					
					if(mat[i][j]<10)
						afi=mat[i][j];
					else
						afi=mat[i][j]+7;
					af=(int)afi+48;
					cout<<af;
				}
				cout<<endl;
				if(i==m-getExtraHeight()-1)
					cout<<endl;
			}
		}
		for(int ii=0;ii<m;ii++)
		   delete [] (mat[ii]);
	
		delete[] (mat);
		
		

		fig << "		2 2 1 1 0 7 50 -1 -1 4.000 0 0 -1 0 0 5" << endl;
		fig << "	  " << (-wInX)*45 << " " << 0 
         << " " << 0 << " " << 0  
         << " " << 0 << " " << (wInY)*45 
         << " " << (-wInX)*45 << " " << (wInY)*45 
         << " " << (-wInX)*45 << " " << 0 << endl;

		
		fig.close();
	}


	/* the rest of the board not tiled by DSPs needs to be implemented using 
	soft-core multipliers */
	int IntTilingMult::partitionOfGridSlices(DSP** config,int &partitions)
	{
		cout << " repartitioning logic multipliers " << endl;
		
		int costSlice = 0, count = 1, nj,ni,njj,nii;

		int n = vn, m = vm, nmew = vnme, mmeh = vmme;
	
		int ew = getExtraWidth();
		int eh = getExtraHeight();

		/* init display matrix */
		for(int i=0; i<m; i++){
			for(int j=0;j<n;j++)
				mat[i][j]=0;
		}
				
		/* interate through the tiled DSPs and fill the matrix accordingly */
		for(int i=0;i<nrDSPs;i++){
			int trX,trY, blX,blY;
		
			config[i]->getTopRightCorner  ( trX, trY);
			config[i]->getBottomLeftCorner( blX, blY);

			trX = n - trX -1;
			blX = n - blX -1;
		
			fillMatrix( mat, n, m, blX, trY, trX, blY, count);
			count++;			
		}

		
		partitions = 0;
	
		/* for every little pixel non-tiled */
		/* ALL THIS PART: can surely be done much simpler */
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				if(mat[i][j] == 0){ 
				/* this is not yet tiled */
					int ver =0;
					int ii = i, jj = j;
					/* iterate through 6 versions of logic multipliers 
					bounded by the tiling board */
					while( ver < 6 && (ii < m-1 || jj < n-1)){
						if (ver < 3){
							if (ver == 0 || ver == 1){
							/* increase the multiplier size on the i direction (Y) */
								ii++;
							} 
							if (ii > m-1){
							/* if we have increased it too much, we set it to the maximum
							direction on Y and we set the filling direction to 2 */
								ii = m-1;
								ver = 2;							
							}
							if (ver == 0 || ver == 2){
							/* increase the dimension in the X direction */
								jj++; 
							}
							if(jj > n-1){
							/* if we have increased it too much, then we reset
							it */
								jj = n-1;
								ver = 1;
							}
										
							for(int k = ii, q = jj; k > i-1 && ( ver == 0 || ver == 2); k--){
								if ( mat[k][q] != 0){
								/* the soft-core multiplier overlaps with some other entity */
									if (ver == 0)
										ver = 1;
									else
										ver = 3;
									
									jj--; 
								}
							}
									
							for(int k = ii, q = jj; q > j-1 && ( ver == 0 || ver == 1); q--){
								if ( mat[k][q] != 0){
									if (ver == 0)
										ver = 2;
									else
										ver = 3;
									
									ii--;
								}
							}
						} else {
							if (ver == 3 || ver == 5)
								jj++;
							
							if (jj > n-1){
								jj = n-1;
								ver = 4;
							}
							
							if(ver==3||ver==4)
								ii++;
								
							if (ii > m-1){
								ii=m-1;
								ver = 5;							
							}
							
							for(int k=ii, q=jj; q>j-1 && (ver==3||ver==4);q--){
								if(mat[k][q]!=0){
									if( ver==3 )
										ver=5;
									else
										ver=6;
									
									ii--;
								}
							}
								
							for(int k=ii,q=jj;k>i-1&&(ver==3||ver==5);k--){
								if(mat[k][q]!=0){
									if(ver==3)
										ver=4;
									else
										ver=6;
									jj--;
								}
							}

							if(ver==5&&jj==n-1)
								ver=6;
							if(ver==4&&ii==m-1)
								ver=6;
						}
					}
					
					if (! (j>=nmew || jj< ew || i >= mmeh || ii < eh)){
						cerr << "j=" <<j << " jj="<<jj<<" i="<<i<<" ii="<<ii<<endl;

						/* if the logic block was not outside the board */
									
						if( j < ew )
							nj = ew ;
						else
							nj = j;
				
						if( jj >= nmew )
							njj = nmew -1;
						else
							njj = jj;
						
						if( i < eh )
							ni = eh ;
						else
							ni = i;
				
						if( ii >=mmeh)
							nii = mmeh -1;
						else
							nii = ii;
						
						partitions++;
						cerr << " Need Logic Multiplier of size = " << njj-nj+1 << " X "<< nii-ni+1 << endl;
						costSlice += target_->getIntMultiplierCost(njj-nj+1,nii-ni+1);
						
						fillMatrix(mat,n,m,j,i,jj,ii,count);
						count++;
					}
				}
			}
		}
		return costSlice;
	}	

	int IntTilingMult::bindDSPs4Virtex(DSP** &config)
	{
		int nrOfUsedDSPs = nrDSPs;
		DSP* ref;
		
		sortDSPs(config);
			
		int itx,ity, ibx,iby, jtx,jty, count;
		
		for(int i=0; i< nrDSPs; i++){
			if ((config[i]!=NULL) && (config[i]->getShiftIn()==NULL)){
				ref=config[i];
				
				count = ref->getNrOfPrimitiveDSPs(); 

				bool ver = true;
				int rw,rh;
				int sa;
				
				/* now we start binding DSPs */		
				while ( ver==true && ref->getShiftOut() == NULL && count < nrOfShifts4Virtex){
					ver=false;

					/* BINDING DSPs in crucial so this function MUST be optimized
					 TODO */
					for(int j=0; j < nrDSPs &&  ref->getShiftOut()==NULL; j++){
						
						ref->getTopRightCorner  (itx,ity);
						ref->getBottomLeftCorner(ibx,iby);

						rw = ref->getMaxMultiplierWidth();
						rh = ref->getMaxMultiplierHeight();

						/* we can potentially bind two if by combining the two
						we don't exceed the maximum number of shift combinations */
						if ( config[j]!=NULL && j!=i && ((count + config[j]->getNrOfPrimitiveDSPs()) <= nrOfShifts4Virtex)){
							config[j]->getTopRightCorner(jtx,jty);
							/* if this one is in the bounds of the board */
							if ( jtx<=vnme && jty<vmme){
								cout<<" itx ity "<<itx<<" "<<ity<<" jtx jty "<< jtx<<" "<<jty<<endl;
								sa = config[j]->getShiftAmount(); //17 for Xilinx
								/* for now this condition binds DSPs on horizontal
								we might need to change this one for optimality */
								if ( (jtx + jty == itx + ity) && config[j]->getShiftIn()==NULL){
									cout<<"DSP #"<<i<<" bind with DSP# "<<j<<"on direct line"<<endl;
									ver = true;
									ref->setShiftOut(config[j]);
									config[j]->setShiftIn(ref);
								
									nrOfUsedDSPs--;
									ref = config[j];
									count += ref->getNrOfPrimitiveDSPs();
//								}else if ( (jtx == ibx+1 && jty==ity     && rw==sa && config[j]->getShiftIn()==NULL) || 
//										   (jtx == itx   && jty==iby+1   && rh==sa && config[j]->getShiftIn()==NULL)){
								}else if (jtx + jty - itx - ity == sa &&  config[j]->getShiftIn()==NULL) {
									/* this is the binding condition */
									cout<<"DSP #"<<i<<" bind with DSP# "<<j<<" on shift line" << endl;
									ver = true;
									ref->setShiftOut(config[j]);
									config[j]->setShiftIn(ref);
								
									nrOfUsedDSPs--;
									ref = config[j];
									count += ref->getNrOfPrimitiveDSPs();
								}else{}
							}
						}
					}
					cout << "exit for ..." << endl;
				}
			}
		}
		return nrOfUsedDSPs;
	}

	/* sort the DSP blocks in the config list according to their Top
	left corner coordinates */
	void IntTilingMult::sortDSPs(DSP** &config){
		int ix,iy,jx,jy;
		DSP* temp;
		for(int i=0; i< nrDSPs-1;i++){		
			for(int j=i+1;j<nrDSPs;j++){
				config[i]->getTopRightCorner(ix,iy);
				config[j]->getTopRightCorner(jx,jy);
				if (iy + ix > jy + jx)  {
					temp      = config[i];
					config[i] = config[j];
					config[j] = temp;				
				}
			}
		}
	}


	/* binds DSPs into supertiles */
	int IntTilingMult::bindDSPs(DSP** &config){
		if  ( (target_->getID() == "Virtex4") ||
		      (target_->getID() == "Virtex5") ||
		      (target_->getID() == "Spartan3")){  // then the target is a Xilinx FPGA
			return bindDSPs4Virtex(config);
		}else{ // the target is Stratix
			return bindDSPs4Stratix(config);
		}
	}


	/* compares the cost of the global config with best config and sets 
	the best config accordingly */
	void IntTilingMult::compareCost()
	{
		for(int ii=0;ii<nrDSPs;ii++)
			memcpy(tempc[ii],globalConfig[ii],sizeof(DSP));
	
		float temp = computeCost(tempc);
		
		if(temp < bestCost){
			bestCost = temp;
			cerr << "---------------------------------------New best score is: " << bestCost << endl;

			for(int ii=0;ii<nrDSPs;ii++)
				memcpy(bestConfig[ii],globalConfig[ii],sizeof(DSP) );
		}else if (temp == bestCost ){
			if (compareOccupation(tempc)==true) {
				REPORT(INFO, "Interchange for equal cost. Now best has cost "<<temp);
				bestCost = temp;
				for(int ii=0;ii<nrDSPs;ii++)
					memcpy(bestConfig[ii],globalConfig[ii],sizeof(DSP) );
			}
		}
	}

	/* computes the cost of a given configuration */
	float IntTilingMult::computeCost(DSP** &config)
	{
		float acc=0.0;
		int nrOfUsedDSPs = nrDSPs;
		acc += ((float)nrDSPs)*costDSP;
		
		int partitions = 0; 
		float LUTs4Multiplication =  partitionOfGridSlices (config, partitions);
	
		cerr << " -----------------------------------Number of slices 4 multiplication of the rest is "<<LUTs4Multiplication<<endl;
		acc += costLUT * LUTs4Multiplication;
		
		/* some of the DSPs can be bined in order to save some additions */
		nrOfUsedDSPs = bindDSPs(config); 
		cout << "------------------------------------ number of used DSPs = "<< nrOfUsedDSPs << endl;

		float LUTs4NAdder=((float)target_->getIntNAdderCost(wInX + wInY, nrOfUsedDSPs+partitions) );
		cerr << "  ---------------------------------- AdderCost = " << LUTs4NAdder << endl;
		acc +=  LUTs4NAdder * costLUT;	
		return acc;
	}


	/* this function might need to be changed depending on the meaning of 
	ratio */	
	int IntTilingMult::estimateDSPs()
	{
		if (ratio > 1){
			REPORT(INFO, "Ratio is " << ratio << " and should be in [0,1]. Will saturate to 1");
			ratio = 1;
		}
			
		float t1,t2,t3,t4; /* meaningful vars ;) */
		int Xd, Yd;  /* the dimension of the multiplier on X and on Y */
		target_->getDSPWidths(Xd,Yd);
		bool fitMultInDSPs = true;

		int maxDSP, mDSP = target_->getNumberOfDSPs();
		int wInXt = wInX;
		int wInYt = wInY;
		if (wInY > wInX){
			wInYt = wInX;
			wInXt = wInY;
		}

		/* all trivial tiling possibilities */
		t1 = ((float) wInXt) / ((float) Xd);
		t2 = ((float) wInYt) / ((float) Yd);
		t3 = ((float) wInXt) / ((float) Yd);
		t4 = ((float) wInYt) / ((float) Xd);
		
		maxDSP = int ( max( ceil(t1)*ceil(t2), ceil(t3)*ceil(t4)) );
	
		if(maxDSP > mDSP){ /* there are not enough DSPs on this FPGA to perform 
			this multiplication using just DSPs */
			fitMultInDSPs = false;
			maxDSP = mDSP; //set the maximum number of DSPs to the multiplication size
		}
			
		return int(ceil((float)maxDSP * ratio));
	}
	
	
	
	
	int  IntTilingMult::getExtraHeight()
	{

		
		int x,y;	
		target_->getDSPWidths(x,  y);
		float temp = ratio * 0.75 * ((float) y);
		return ((int)temp);
		//return 4;

	}

	
	int  IntTilingMult::getExtraWidth()
	{

		int x,y;	
		target_->getDSPWidths(x,y);
		float temp = ratio * 0.75 * ((float) x);
		return ((int)temp);
		//return 4;

	}





	void IntTilingMult::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		mpz_class svR = svX * svY;

		tc->addExpectedOutput("R", svR);
	}

	void IntTilingMult::buildStandardTestCases(TestCaseList* tcl){
	
	}


	int IntTilingMult::checkFarness(DSP** config,int index)
	{
		int xtr1, ytr1, xbl1, ybl1, xtr2, ytr2, xbl2, ybl2;
		//if we run the algorithm with one less DSP then we should verify against the DSP nrDSP-2
	
		if(index  < 1)
			return 0;
	
		config[index]->getTopRightCorner(xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
		int dist=0;
		bool ver=false;
		int sqrDist = maxDist2Move * maxDist2Move;
	
		dist=0;	
		for (int i=0; i<index; i++)
			if ((config[i] != NULL) )
				{
					config[i]->getTopRightCorner(xtr2, ytr2);
					config[i]->getBottomLeftCorner(xbl2, ybl2);
					if(xtr1 > xbl2+1 + maxDist2Move)
						dist++;			
				}
		if(dist == index)
			return 1;
	
		//cout<<"Sqr max is "<<sqrDist<<endl;
		for (int i=0; i<index; i++)
			if ((config[i] != NULL) )
				{
					dist=0;
					config[i]->getTopRightCorner(xtr2, ytr2);
					config[i]->getBottomLeftCorner(xbl2, ybl2);
					if(xtr1 > xbl2)
						dist +=(xtr1 - xbl2-1)*(xtr1 - xbl2-1);
					//~ cout<<"xtr1= "<<xtr1<<" xbl2= "<<xbl2<<" Dist = "<<dist<<"  ";
					if(ybl2 < ytr1)
						dist +=(1+ybl2-ytr1) * (1+ybl2-ytr1) ;
					else
						if(ybl1 < ytr2)
							{
								dist +=(1+ybl1 - ytr2) *(1+ybl1 - ytr2) ;
								ver =true;
							}
			
					//~ cout<<"The distance for DSP "<<i<<" is "<<dist<<endl;
					if(dist <= sqrDist)
						{
							//cout<<"Dist  was = "<<dist<<endl;
							return 0; //OK
					
						}			
				}
	
	
	
		if(!ver)
			return 2;
	
		return 3;
	
	}

	/* we check the validity of our tiling */
	bool IntTilingMult::checkOverlap(DSP** config, int index)
	{	
		int xtr1, ytr1, xbl1, ybl1, xtr2, ytr2, xbl2, ybl2;

		config[index]->getTopRightCorner  (xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
		
		if(verbose)
			cout << tab << tab << "checkOverlap: ref is block #" << index << ". Top-right is at (" << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
	
		for (int i=0; i<index; i++){
			if (config[i] != NULL) /* it should be different, morally ... */
			{
				config[i]->getTopRightCorner  (xtr2, ytr2);
				config[i]->getBottomLeftCorner(xbl2, ybl2);
			
				if(verbose)
					cout << tab << tab << "checkOverlap: comparing with block #" << i << ". Top-right is at (" << xtr2 << ", " << ytr2 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
		
				if (((xtr2 <= xbl1) && (ytr2 <= ybl1) && (xtr2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the upper and/or right part(s) of config[i]
				((xbl2 <= xbl1) && (ybl2 <= ybl1) && (xbl2 >= xtr1) && (ybl2 >= ytr1)) || // config[index] overlaps the bottom and/or left part(s) of config[i]
				((xtr2 <= xbl1) && (ybl2 <= ybl1) && (xtr2 >= xtr1) && (ybl2 >= ytr1)) || // config[index] overlaps the upper and/or left part(s) of config[i]
				((xbl2 <= xbl1) && (ytr2 <= ybl1) && (xbl2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the bottom and/or right part(s) of config[i]
				((xbl2 >= xbl1) && (ybl2 <= ybl1) && (ytr2 >= ytr1) && (xtr2 <= xtr1)) || // config[index] overlaps the center part of config[i]
			
				((xtr1 <= xbl2) && (ytr1 <= ybl2) && (xtr1 >= xtr2) && (ytr1 >= ytr2)) || // config[i] overlaps the upper and/or right part(s) of config[index]
				((xbl1 <= xbl2) && (ybl1 <= ybl2) && (xbl1 >= xtr2) && (ybl1 >= ytr2)) || // config[i] overlaps the bottom and/or left part(s) of config[index]
				((xtr1 <= xbl2) && (ybl1 <= ybl2) && (xtr1 >= xtr2) && (ybl1 >= ytr2)) || // config[i] overlaps the upper and/or left part(s) of config[index]
				((xbl1 <= xbl2) && (ytr1 <= ybl2) && (xbl1 >= xtr2) && (ytr1 >= ytr2)) || // config[i] overlaps the bottom and/or right part(s) of config[index]
				((xbl1 >= xbl2) && (ybl1 <= ybl2) && (ytr1 >= ytr2) && (xtr1 <= xtr2)))    // config[i] overlaps the center part of config[index]
					return true;
			}
		}		
		if(verbose)
			cout << tab << tab << "checkOverlap: return false" << endl;	
		return false;
	}


	/**
		There is one case that is not resolved yet. For DSP with widths 
		different then height the algorithm should move the dsp with both values
	*/

	bool IntTilingMult::move(DSP** config, int index)
	{
		int xtr1, ytr1, xbl1, ybl1;
		int w, h;
		//target->getDSPWidths(w,h);
		w= config[index]->getMaxMultiplierWidth();
		h= config[index]->getMaxMultiplierHeight();
	
		config[index]->getTopRightCorner(xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
	
		if(verbose)
			cout << tab << "replace : DSP #" << index << " width is " << w << ", height is " << h << endl; 
		/*
		 * Initialized in the constructor:
		 * vn = wInX+2*getExtraWidth(); width of the tiling grid including extensions
		 * vm = wInY+2*getExtraHeight(); height of the tiling grid including extensions
		 * */
		//~ int exh = getExtraHeight();
		//~ int exw = getExtraWidth();
		int pos; // index for list of positions of a DSP
		
		
		if(index==0) // the first DSP block can move freely on the tiling grid
		{
			return false;
		}
		else // all DSP blocks except the first one can move only in fixed positions
		{
			do{
				// move to next position
				pos = config[index]->pop();
				if(pos >= 0){
					ytr1 = config[index]->Ypositions[pos];
					ybl1 = ytr1 + h-1;
					xtr1 = config[index]->Xpositions[pos];
					xbl1 = xtr1 + w-1;
				}
				else
					return false;
					
				config[index]->setTopRightCorner(xtr1, ytr1);
				config[index]->setBottomLeftCorner(xbl1, ybl1);
			}while (checkOverlap(config, index));
		}
		return true;
	}

	bool IntTilingMult::replace(DSP** config, int index)
	{
		int xtr1, ytr1, xbl1, ybl1;
		int w, h;
		string targetID = target_->getID();
		
		w= config[index]->getMaxMultiplierWidth();
		h= config[index]->getMaxMultiplierHeight();
		config[index]->setPosition(0);
		
		if (index > 1)
		{// take all positions from the previous DSP
			int curpos =config[index-1]->getCurrentPosition();
			int avpos = config[index-1]->getAvailablePositions();
				memcpy(config[index]->Xpositions, config[index-1]->Xpositions + curpos, sizeof(int)*(avpos - curpos));	
				memcpy(config[index]->Ypositions, config[index-1]->Ypositions + curpos, sizeof(int)*(avpos - curpos));	
				config[index]->setPosition((avpos - curpos));
		}
		
		if (index > 0){
			/* starting with the second DSP */
			w = config[index-1]->getMaxMultiplierWidth();
			h = config[index-1]->getMaxMultiplierHeight();

			/* calulate the minimum width of the DSP */
			int mind = min(w,h);
			int mindX = mind;
			int mindY = mind;
			
			if (targetID == "Virtex5"){ 
				/* align top right of current 24-17 bits left, 24 bits lower */
				mindX = -abs(w-h);
				mindY = h;
			}
			else if ((targetID == "StratixII") || (targetID == "StratixIII") || (targetID == "StratixIV")){ 
				// align on diagonal in order to use internal adders 
				mindX = -w;
				mindY = h;	
			}

			//TODO add positions from article
			int positionDisplacementX[] = {0, w, mindX};
			int positionDisplacementY[] = {h, 0, mindY};

			int x,y, x1,y1, pos;
			config[index-1]->getTopRightCorner(x,y);
			
			bool extraPosition = ( (w!=h) || ((targetID == "StratixII") || (targetID == "StratixIV") || (targetID == "StratixIV") ));
			for (int i=0; i<3; i++){
				x1 = x + positionDisplacementX[i];
				y1 = y + positionDisplacementY[i];
				if (x1 < vnme && y1 < vmme && x1>0 && y1>0) /* if inside the board */ 
					if ( (i!=2) || extraPosition) /* for Virtex5 and Altera FPGAs */
					config[index]->push(x1, y1);
			}
			
			w= config[index]->getMaxMultiplierWidth();
			h= config[index]->getMaxMultiplierHeight();

			config[index]->resetPosition();
			do{
				/* we previously added a few trial positions for this block,
				now it's time to try them out */
				/* go to next position in list */
				pos = config[index]->pop();
				if(pos >= 0){
					ytr1 = config[index]->Ypositions[pos];
					ybl1 = ytr1 + h-1;
					xtr1 = config[index]->Xpositions[pos];
					xbl1 = xtr1 + w-1;
				}else{
					/* place DSP outside board: it should never get here */
					config[index]->setTopRightCorner(vn, vm);
					config[index]->setBottomLeftCorner(vn+w, vm+h);
					return false;
				}
				
				/* we properly write-out the position of the next DSP */
				config[index]->setTopRightCorner(xtr1, ytr1);
				config[index]->setBottomLeftCorner(xbl1, ybl1);
			
			/* we validate this position if no overlaping exists with the 
			curent index-1 configuration */
			}while (checkOverlap(config, index)); 
			
			/* if we got here then we found a good position for our multiplier*/	
			return true;
		}else{	
			/* if this DSP is the firs one to be placed */			
			int exh = getExtraHeight();
			int exw = getExtraWidth();
		
			/* tiling board has a band around it. 
			
			------------------------
			' ___________________  '
			' |                  | '
			' |                  | '
			' |                  | '
			' |__________________| '
			'                      '
			'----------------------
			the next places the first multiplier in the top right corner
			*/
			 
			xtr1 = exw ;
			ytr1 = exh ;	
			ybl1 = ytr1 + h-1;
			xbl1 = xtr1 + w-1;
			
			config[index]->setTopRightCorner(xtr1, ytr1);
			config[index]->setBottomLeftCorner(xbl1, ybl1);
		
			REPORT( DETAILED, "replace : DSP width is " << w << ", height is " << h);
		}
		return true;
	}


	/* do an initial tiling of the board using the dspCount multipliers */
	void IntTilingMult::initTiling(DSP** &config, int dspCount)
	{
		int w, h; 
		target_->getDSPWidths(w, h);
		dsplimit = w*h;
		
		/* we initially allocate an array of empty pointers to our 
		tiling configuration */
		config = new DSP*[nrDSPs];
		for (int i=0; i<nrDSPs; i++)
			config[i] = NULL;
		
		for (int i=0; i<dspCount; i++){
			REPORT(DETAILED, "initTiling : iteration #" << i); 
			config[i] = target_->createDSP();						
			config[i]->setNrOfPrimitiveDSPs(1);
			
			/* each DSP offers 8 positions, actually a lot less but 8 is 
			for margin */
			config[i]->allocatePositions(3*i); 

			/* try to place the DSP ... */
			if(!replace(config, i)){
				w = config[i]->getMaxMultiplierWidth();
				h = config[i]->getMaxMultiplierHeight();
				config[i]->setTopRightCorner(vn, vm);
				config[i]->setBottomLeftCorner(vn+w, vm+h);
			}
		}
	}
	
	void IntTilingMult::initTiling2(DSP** &config, int dspCount)
	{
		//nrOfShifts4Virtex =2; 
		int w, h;
		target_->getDSPWidths(w, h);
		int min = h;
		int mw,mh;
		mw = vn - 2*getExtraWidth();
		mh = vm - 2*getExtraHeight();
		if(w < h)
		{
			min = w;
			w = h;
			h = min;
		}	
		
		// allocate the maximum number of DSP objects
		config = new DSP*[nrDSPs];
		for (int i=0; i<nrDSPs; i++)
		{	
			config[i] = NULL;
		}
		
		/* NOTE: In case the multiplication is narrower than twice the minimum dimension of a DSP 
		 * then we do not form paris of multipliers
		 */
		if ((vm < 2*min) || (vn < 2*min))
		{
			for (int i=0; i<dspCount; i++)
			{
				config[i] = target_->createDSP();
				config[i]->setNrOfPrimitiveDSPs(1);
				config[i]->allocatePositions(3*i); // each DSP offers 3 positions
				if(!replace(config, i))
				{
					w=config[i]->getMaxMultiplierWidth();
					h= config[i]->getMaxMultiplierHeight();
					config[i]->setTopRightCorner(vn, vm);
					config[i]->setBottomLeftCorner(vn+w, vm+h);
				}
			}
			dsplimit = w*h;
			return; // force exit
		}
		dsplimit = w*h*2;
		// else form paris
		// verify if we have an even or odd nr of DSPs
		bool singleDSP = false;
		if (nrDSPs % 2 == 1)
			singleDSP = true;
		
		int nrDSPsprim = nrDSPs;
			
		// compute new number of DSPs

		numberDSP4Overlap = dspCount = nrDSPs = (int) ceil((double) dspCount/2);

		
		// set shift amount according to target
		int shift = 0;
		if ((target_->getID() == "Virtex4") || (target_->getID() == "Virtex5"))
			shift = 17;
		
		int start = 0; // starting position
		
		
		//~ cout<< "mw "<<mw<<" mh "<<mh<<" w "<<w<<" h "<<h<<"estimate dsps "<<nrDSPsprim<<endl;
		//~ cout<<"h*0.9 + h*4="<<h*0.9 + h*4<<endl;
		//~ cout<<"cond "<<(mw == h*4 || mh ==h*4)<<" tout "<<((mw == h*4 || mh ==h*4) ||(  (h*0.9 + h*4 <=mw) &&  mh>=2*w  ) || (   (h*0.9 + h*4 <=mh) &&  mw>=2*w  ))<<endl;
		//~ cout<<" area= "<<mw*mh<<" area covered by DSPs= "<<h*w*nrDSPsprim<<" condition "<<( mw*mh*0.8>= h*w*nrDSPsprim )<<endl;
		
		
		
		
		if (dspCount >= 16 &&
			((mw == h*4 || mh ==h*4) ||
			(  (h*0.9 + h*4 <=mw) &&  mh>=2*w  ) ||
			(   (h*0.9 + h*4 <=mh) &&  mw>=2*w  ) ||
			( mw*mh>= h*w*nrDSPsprim )
			)) 
		{ // we have more than 16 paris of multipliers we can group 4 such pairs into a super-block
			cout<<"A super DSP was created"<<endl;
			config[start] = new DSP(0, w*2, h*4);
			config[start]->setNrOfPrimitiveDSPs(4);
			config[start]->allocatePositions(3*start); // each DSP offers 3 positions
			if(!replace(config, start))
			{
				int w=config[start]->getMaxMultiplierWidth();
				int h= config[start]->getMaxMultiplierHeight();
				config[start]->setTopRightCorner(vn, vm);
				config[start]->setBottomLeftCorner(vn+w, vm+h);
			}
			start++;
			dspCount -= 3;

			numberDSP4Overlap = nrDSPs = dspCount;

			
			int i;
			for (i=start; i<3; i++)
			{
				config[i] = new DSP(shift, w, h*2);	
				config[i]->setNrOfPrimitiveDSPs(2);
				config[i]->allocatePositions(3*i); // each DSP offers 3 positions
				if(!replace(config, i))
				{
					int w=config[i]->getMaxMultiplierWidth();
					int h= config[i]->getMaxMultiplierHeight();
					config[i]->setTopRightCorner(vn, vm);
					config[i]->setBottomLeftCorner(vn+w, vm+h*2);
				}
			}
			start = i;	
		}

		/*NOTE: if the program entered the previous if clause it will also enter the following 
		 * 	if clause. This is the effect we want to obtain */
		if (dspCount >= 10 		&&
			((mw == h*4 || mh ==h*4) ||
			(  (h*0.9 + h*4 <=mw) &&  mh>=2*w  ) ||
			(   (h*0.9 + h*4 <=mh) &&  mw>=2*w  ) ||
			( mw*mh>= h*w*nrDSPsprim )
			)	
		)
		{ // we have more than 10 paris of multipliers we can group 4 such pairs into a super-block
			cout<<"A super DSP was created"<<endl;
			config[start] = new DSP(0, w*2, h*4);
			config[start]->setNrOfPrimitiveDSPs(4);
			config[start]->allocatePositions(3*start); // each DSP offers 3 positions
			dspCount -= 3;
			numberDSP4Overlap = nrDSPs = dspCount;
			if(!replace(config, start))
			{
				int w=config[start]->getMaxMultiplierWidth();
				int h= config[start]->getMaxMultiplierHeight();
				config[start]->setTopRightCorner(vn, vm);
				config[start]->setBottomLeftCorner(vn+w, vm+h);
			}
			start++;

		}
		
		// initialize all DSP paris except first one
		dspCount--;
		for (int i=start; i<dspCount; i++)
		{
			config[i] = new DSP(shift, w, h*2);	
			config[i]->setNrOfPrimitiveDSPs(2);
			config[i]->allocatePositions(3*i); // each DSP offers 3 positions
			if(!replace(config, i))
			{
				int w=config[i]->getMaxMultiplierWidth();
				int h= config[i]->getMaxMultiplierHeight();
				config[i]->setTopRightCorner(vn, vm);
				config[i]->setBottomLeftCorner(vn+w, vm+h*2);
			}
		}
		
		//initializa last DSP or DSP pair
		if (singleDSP) // then the last position is a single DSP 
		{ // allocate single DSP
			config[dspCount] = target_->createDSP();
			config[dspCount]->setNrOfPrimitiveDSPs(1);
			config[dspCount]->allocatePositions(3*dspCount); // each DSP offers 3 positions
			if(!replace(config, dspCount))
			{
				int w=config[dspCount]->getMaxMultiplierWidth();
				int h= config[dspCount]->getMaxMultiplierHeight();
				config[dspCount]->setTopRightCorner(vn, vm);
				config[dspCount]->setBottomLeftCorner(vn+w, vm+h);
			}
			
		}
		else // then the last position is a DSP pair
		{ // allocate DSP pair
			config[dspCount] = new DSP(shift, w, h*2);
			config[dspCount]->setNrOfPrimitiveDSPs(2);
			config[dspCount]->allocatePositions(3*dspCount); // each DSP offers 3 positions
			if(!replace(config, dspCount))
			{
				int w=config[dspCount]->getMaxMultiplierWidth();
				int h= config[dspCount]->getMaxMultiplierHeight();
				config[dspCount]->setTopRightCorner(vn, vm);
				config[dspCount]->setBottomLeftCorner(vn+w, vm+h*2);
			}
			
		}	
		
	}
	
	DSP** IntTilingMult::neighbour(DSP** config)
	{
		DSP** returnConfig;
		initTiling(returnConfig, nrDSPs);
		// select random DSP to move
		srand (time(NULL));
		int index = rand() % nrDSPs;
		// get multiplier width and height
		int w, h, pos;
		w= returnConfig[index]->getMaxMultiplierWidth();
		h= returnConfig[index]->getMaxMultiplierHeight();
		
		if (index == 0) // we have a special case for the first DSP
		{
			returnConfig[index]->rotate();
			index++;
		}
		
		// for all DSPs except the first one
		
		// nr of possible positions is 3 if multiplier is sqare and 6 if it is rectangular
		int nrPositions = 6;
		
		if (w == h)
			nrPositions = 3;
		
		// starting from the selected DSP move all susequent DSPs to random positions
		for (int i=index; i<nrDSPs; i++)
		{
			// move selected DSP in a random position
			pos = rand() % nrPositions;
			if (pos >= 3) // rotate the multiplier
			{
				config[i]->rotate();
				pos -=3;
			}
			if (replace(returnConfig, i))
			{	
				pos--; // one position was used on replace
				for (int j=0; j<pos; j++) // try to move DSP a number of times
					if (!move(returnConfig, i)) 
					{ // just replace if not posible to move	
						replace(returnConfig, i);
						break;
					}
			}
			else // cannot replace the DSP
				memcpy(returnConfig[i], config[i], sizeof(DSP)); // copy back the initial position
		}
		
		return returnConfig;	
	}
	
	float IntTilingMult::temp(float f)
	{
		return f;
	}
	
	float IntTilingMult::probability(float e, float enew, float t)
	{
		if (enew < e)
			return 1.0;
		else
		{
			float dif = enew-e;
			if (dif < 1)
				return dif*t;
			else
				return dif*t/enew;
		}
	}
	
	void IntTilingMult::simulatedAnnealing()
	{
		// Initial state
		initTiling(globalConfig, nrDSPs);
		int exw = getExtraWidth();
		int exh = getExtraHeight();
		int h = globalConfig[0]->getMaxMultiplierHeight();
		int w = globalConfig[0]->getMaxMultiplierWidth();
		globalConfig[0]->setTopRightCorner(exw, exh);
		globalConfig[0]->setBottomLeftCorner(exw+w-1, exh+h-1);
		for (int i=1; i<nrDSPs; i++)
			replace(globalConfig, i);
			
		display(globalConfig);
		float e = computeCost(globalConfig);// Initial energy.
		bestConfig = new DSP*[nrDSPs];// Allocate blocks
			
		for (int i=0; i<nrDSPs; i++)// Initial "best" solution
		{
			bestConfig[i] = new DSP();
			memcpy(bestConfig[i], globalConfig[i], sizeof(DSP)); 				
		}
		bestCost = e;// Initial "best" cost
		int k = 0;// Energy evaluation count.
		int kmax = 3*nrDSPs*nrDSPs;// Arbitrary limit
		DSP **snew;
		while (k < kmax)// While time left & not good enough:
		{
			// TODO deallocate previous snew;
			snew = neighbour(globalConfig);// Pick some neighbour.
			float enew = computeCost(snew);// Compute its energy.
			
			display(snew);
			
			if (enew < bestCost)// Is this a new best?
			{// Save 'new neighbour' to 'best found'.	
				for(int i=0; i<nrDSPs; i++)
					memcpy(bestConfig[i], snew[i], sizeof(DSP));
				bestCost = enew;                  
			}
			
			if (probability(e, enew, temp(k/kmax)) > rand())// Should we move to it?
			{// Yes, change state.	
				globalConfig = snew; 
				e = enew;         
			}                 
			k = k + 1; // One more evaluation done
		}
		display(bestConfig);
		// Return the best solution found.
		for(int i=0; i<nrDSPs; i++)
			memcpy(globalConfig[i], bestConfig[i], sizeof(DSP));
	}
	
	int IntTilingMult::bindDSPs4Stratix(DSP** config)
	{
		int nrDSPs_ = nrDSPs;
		int DSPcount = nrDSPs_;
	
		for (int i=0; i<nrDSPs_; i++)
			if (config[i] == NULL)
				DSPcount--;
		
			
		for (int i=0; i<DSPcount; i++)
			if (config[i]->getNumberOfAdders() < 3) // serach for an aligned DSP block that can be added with this one
				{
					int xtri, ytri, wx, wy, nrOp;
					bool bound;
					config[i]->getTopRightCorner(xtri, ytri);
					target_->getDSPWidths(wx,wy);
				
					if (verbose)
						cout << "bindDSP4Stratix: DSP #" << i << " has less than 3 operands. Top-right is at ( " << xtri << ", " << ytri << ") width is " << wx << endl;
					DSP** operands;
					DSP** operandsj;
					DSP** opsk;
					for (int j=0; j<DSPcount; j++)
						if (i != j)
							{
								nrOp = config[i]->getNumberOfAdders();
								if (nrOp == 3)
									break;
			
								operands = config[i]->getAdditionOperands();
								bound = false;
								// check if the DSP blocks are bound
								for (int k=0; k<nrOp; k++)
									if (operands[k] == config[j]) // the DSPs are allready bound
										{
											if (verbose)
												cout << "bindDSP4Stratix: DSP #" << j << " has #" << i << " as one of its operands allready." << endl;
											bound = true;
											break;
										}
				
								if (bound)
									continue;
								// if they are not bound
								int xtrj, ytrj, nrOpj;
								config[j]->getTopRightCorner(xtrj, ytrj);
								nrOpj = config[j]->getNumberOfAdders();
				
								if (verbose)
									cout << "bindDSP4Stratix:" << tab << " Checking against DSP #" << j << " Top-right is at ( " << xtrj << ", " << ytrj << ") width is " << wx << endl;
			
								if ((((xtrj-xtri == -wx) && (ytrj-ytri == wy)) || ((xtrj-xtri == wx) && (ytrj-ytri == -wy))) && // if they have the same alignment
									 (nrOpj + nrOp < 3)) // if the DSPs we want to bind have less than 3 other DSPs to add together
									{ // copy the operands from one another and bind them 
					
										if (verbose)
											cout << "bindDSP4Stratix : DSP #" << j << " together with #" << i << " have fewer than 3 operands. We can bind them." << endl;
										operandsj = config[j]->getAdditionOperands();
					
										for (int k=0; k<nrOp; k++)
											{
												operandsj[nrOpj+k] = operands[k];
												// each operand of congif[i] also gets bounded 
												int opcntk = operands[k]->getNumberOfAdders();
												operands[k]->setNumberOfAdders(nrOp+nrOpj+1);
												opsk = operands[k]->getAdditionOperands();
												for (int l=0; l < nrOpj; l++)
													opsk[l+opcntk] = operandsj[l];
												opsk[nrOpj+opcntk] = config[j];
												operands[k]->setAdditionOperands(opsk);
											}
										operandsj[nrOp+nrOpj] = config[i];
										config[j]->setAdditionOperands(operandsj);
										config[j]->setNumberOfAdders(nrOp+nrOpj+1);
					
										for (int k=0; k<nrOpj; k++)
											{
												operands[nrOp+k] = operandsj[k];
												// each operand of congif[j] also gets bounded 
												int opcntk = operandsj[k]->getNumberOfAdders();
												operandsj[k]->setNumberOfAdders(nrOp+nrOpj+1);
												opsk = operandsj[k]->getAdditionOperands();
												for (int l=0; l < nrOp; l++)
													opsk[l+opcntk] = operands[l];
												opsk[nrOp+opcntk] = config[i];
												operandsj[k]->setAdditionOperands(opsk);
											}
										operands[nrOp+nrOpj] = config[j];
										config[i]->setAdditionOperands(operands);
										config[i]->setNumberOfAdders(nrOp+nrOpj+1);
									}
							}
				}
		
		/* We now have:
		 * 	- pairs of DSP objects which have the number of operands set to 1
		 * 	- triplets of DSP objects which have the number of operands set to 2
		 * 	- quadruplets of DSP objects which have the number of operands set to 3 
		 * We keep a counter for each possible number of operands. */
		int pair[3] = {0, 0, 0};
	
	
		for (int i=0; i<DSPcount; i++)
			{
				if (verbose)
					cout << "bindDSP4Stratix : DSP #" << i << " has " << config[i]->getNumberOfAdders() << " adders" << endl;
				pair[config[i]->getNumberOfAdders()-1]++;
			}
		if (verbose)
			{
				cout << "bindDSP4Stratix : one " << pair[0] << endl;
				cout << "bindDSP4Stratix : two " << pair[1] << endl;
				cout << "bindDSP4Stratix : three " << pair[2] << endl;
			}
		return (DSPcount - pair[0]/2 - pair[1]*2/3 - pair[2]*3/4);
	}

	DSP** IntTilingMult::splitLargeBlocks(DSP** config, int &numberOfDSPs)
	{
		int h, w, dspH, dspW, tmp, nrDSPonHeight, nrDSPonWidth, shiftAmount, newNrDSPs=0;
		getTarget()->getDSPWidths(dspW, dspH);
		
		// count total number of DSPs
		for (int i=0; i<numberOfDSPs; i++)
		{
			h = config[i]->getMaxMultiplierHeight();
			w = config[i]->getMaxMultiplierWidth();
			
			if ((h % dspH) != 0) // match width and height
			{
				tmp = dspH;
				dspH = dspW;
				dspW = tmp;
			}
			
			nrDSPonHeight = h/dspH;
			nrDSPonWidth = w/dspW;
			
			newNrDSPs += (nrDSPonHeight*nrDSPonWidth);
		}
		
		DSP** returnConfig = new DSP*[newNrDSPs];
		int index = 0;
		int xtr, xbl, ytr, ybl;
		for (int i=0; i<numberOfDSPs; i++)
		{
			h = config[i]->getMaxMultiplierHeight();
			w = config[i]->getMaxMultiplierWidth();
			shiftAmount = config[i]->getShiftAmount();
			config[i]->getTopRightCorner(xtr, ytr);
			config[i]->getBottomLeftCorner(xbl, ybl);
			
			if ((h % dspH) != 0) // match width and height
			{
				tmp = dspH;
				dspH = dspW;
				dspW = tmp;
			}
			
			nrDSPonHeight = h/dspH;
			nrDSPonWidth = w/dspW;
			// create DSP blocks
			for (int j=0; j<nrDSPonHeight; j++)
				for (int k=0; k<nrDSPonWidth; k++)
				{
					returnConfig[index] = getTarget()->createDSP();
					returnConfig[index]->setTopRightCorner(xtr+k*dspW, ytr+j*dspH);
					returnConfig[index]->setBottomLeftCorner(xtr+(k+1)*dspW-1, ytr+(j+1)*dspH-1);
					returnConfig[index]->setMultiplierHeight(dspH);
					returnConfig[index]->setMultiplierWidth(dspW);
					// take care of shiftings between DSPs
					if (shiftAmount == dspH)
					{
						if (j > 0) 
						{
							returnConfig[index]->setShiftIn(returnConfig[index-(j*nrDSPonWidth)]);
							returnConfig[index-(j*nrDSPonWidth)]->setShiftOut(returnConfig[index]);	
						}
					}
					else if (shiftAmount == dspW)
					{
						if (k > 0)
						{
							returnConfig[index]->setShiftIn(returnConfig[index-1]);
							returnConfig[index-1]->setShiftOut(returnConfig[index]);
						}
					}
					index++;
				}
		}
		numberOfDSPs = newNrDSPs;
		return returnConfig;
	}


	void IntTilingMult::convertCoordinates(int &tx, int &ty, int &bx, int &by)
	{
		tx -= getExtraWidth();
		ty -= getExtraHeight();
		bx -= getExtraWidth();
		by -= getExtraHeight();

		if (bx>=wInX) bx = wInX-1;
		if (by>=wInY) by = wInY-1;
	}

	void IntTilingMult::convertCoordinatesKeepNeg(int &tx, int &ty, int &bx, int &by)
	{
		tx -= getExtraWidth();
		ty -= getExtraHeight();
		bx -= getExtraWidth();
		by -= getExtraHeight();
	}



	int IntTilingMult::multiplicationInDSPs(DSP** config)
	{
		int nrOp = 0;		 			// number of resulting adder operands
		int trx1, try1, blx1, bly1; 	// coordinates of the two corners of a DSP block 
		int fpadX, fpadY, bpadX, bpadY;	// zero padding on both axis
		int extW, extH;					// extra width and height of the tiling grid
		int multW, multH; 				// width and height of the multiplier the DSP block is using
		ostringstream xname, yname, mname, cname, sname;
		DSP** tempc = new DSP*[nrDSPs];	// buffer that with hold a copy of the global configuration
			
		memcpy(tempc, config, sizeof(DSP*) * nrDSPs );
	
		if ( ( target_->getID() == "Virtex4") ||
			 ( target_->getID() == "Virtex5") ||
			 ( target_->getID() == "Spartan3"))  // then the target is A Xilinx FPGA 
			{	
				for (int i=0; i<nrDSPs; i++)
					if (tempc[i] != NULL)
						{
							//~ cout << "At DSP#"<< i+1 << " tempc["<<i<<"]" << endl; 
							DSP* d = tempc[i];
							int j=0;
							int connected = 0;
							bool outside = false;
							
							while (d != NULL)
								{
									connected++;
									d = d->getShiftOut();
								}
							REPORT(DETAILED, "CONNECTED ================ " <<connected);
							d = tempc[i];
							
							int startXOld=0, startYOld=0;							
							while (d != NULL)
								{
									d->getTopRightCorner(trx1, try1);
									d->getBottomLeftCorner(blx1, bly1);
									extW = getExtraWidth();
									extH = getExtraHeight();
									
									fpadX = blx1-wInX-extW+1;
									//~ cout << "fpadX = " << fpadX << endl;
									fpadX = (fpadX<0)?0:fpadX;
									
									fpadY = bly1-wInY-extH+1;
									//~ cout << "fpadY = " << fpadY << endl;
									fpadY = (fpadY<0)?0:fpadY;
									
									bpadX = extW-trx1;
									bpadX = (bpadX<0)?0:bpadX;
									
									bpadY = extH-try1;
									bpadY = (bpadY<0)?0:bpadY;
									
									multW = blx1 - trx1 + 1;
									multH = bly1 - try1 + 1;

//									multW = d->getMultiplierWidth();
//									multH = d->getMultiplierHeight();
									
									int startX = blx1-fpadX-extW;
									int endX = trx1+bpadX-extW;
									int startY = bly1-fpadY-extH;
									int endY = try1+bpadY-extH;
									
									if ((startX < endX) || (startY < endY))
									{
										outside = true;
										break;
									}
										
									setCycle(0);
									
									xname.str("");
									xname << "x" << i << "_" << j;
									vhdl << tab << declare(xname.str(), multW+2, true) << " <= \"10\" & " << zg(fpadX,0) << " & X" << range(startX, endX) << " & " << zg(bpadX,0) << ";" << endl;
									
									yname.str("");
									yname << "y" << i << "_" << j;
									vhdl << tab << declare(yname.str(), multH+2, true) << " <= \"10\" & " << zg(fpadY,0) << " & Y" << range(startY, endY) << " & " << zg(bpadY,0) << ";" << endl;
				
									if ((d->getShiftIn() != NULL) && (j>0)) // multiply accumulate
										{
											mname.str("");
											mname << "pxy" << i;
											cname.str("");
											cname << "txy" << i << j;
											setCycle(j);
											vhdl << tab << declare(cname.str(), multW+multH+2) << " <= " << use(xname.str())<<range(multW,0) << " * " << use(yname.str())<<range(multH,0) << ";" << endl;
											
											if (startX+startY == startXOld+startYOld)
												vhdl << tab << declare(join(mname.str(),j), multW+multH+2) << " <= (" << cname.str()<<range(multW+multH+1,0) << ") + (" << join(mname.str(), j-1) << ");" << endl;	
											else
												vhdl << tab << declare(join(mname.str(),j), multW+multH+2) << " <= (" << cname.str()<<range(multW+multH+1,0) << ") + (" <<zg(d->getShiftAmount(),0)<< " &" << join(mname.str(), j-1) << range(multW+multH+1, d->getShiftAmount()) << ");" << endl;	
											if (d->getShiftOut() == NULL) // concatenate the entire partial product
												{
													setCycle(connected);
													sname.seekp(ios_base::beg);
													//sname << zg(wInX-(blx1-fpadX-extW)+wInY-(bly1-fpadY-extH)-2, 0) << " & " << use(join(mname.str(),j)) << range(multW-fpadX + multH-fpadY-1, 0) << " & " << sname.str();
													int nrZeros = wInX-(blx1-extW)-fpadX + wInY-(bly1-extH)-fpadY-3;
													if (nrZeros < 0)
														sname << zg(wInX-(blx1-fpadX-extW)+wInY-(bly1-fpadY-extH)-2, 0) << " & " << use(join(mname.str(),j)) << range(multW-fpadX + multH-fpadY-1, 0) << " & " << sname.str();
//														sname << use(join(mname.str(),j)) << range(multW-fpadX + multH-fpadY-1, 0) << " & " << sname.str();

													else
														sname << zg(wInX-(blx1-extW) + wInY-(bly1-extH)-3, 0) << " & " << use(join(mname.str(),j)) << range(multW + multH, 0) << " & " << sname.str();
	
												}
											else // concatenate only the lower portion of the partial product
												{
													setCycle(connected);
													sname.seekp(ios_base::beg);
													
													/* check if the next element has as 17-bit shift wr to this one */
													int ttrx1,ttry1,bblx1,bbly1, ffpadX, ffpadY, bbpadX, bbpadY;
													d->getShiftOut()->getTopRightCorner(ttrx1, ttry1);
													d->getShiftOut()->getBottomLeftCorner(bblx1, bbly1);
									
													ffpadX = bblx1-wInX-extW+1;
													ffpadX = (ffpadX<0)?0:ffpadX;
									
													ffpadY = bbly1-wInY-extH+1;
													ffpadY = (ffpadY<0)?0:ffpadY;
									
													bbpadX = extW-ttrx1;
													bbpadX = (bbpadX<0)?0:bbpadX;
									
													bbpadY = extH-ttry1;
													bbpadY = (bbpadY<0)?0:bbpadY;
									
									
													int startXNext = bblx1-ffpadX-extW;
													int endXNext = ttrx1+bbpadX-extW;
													int startYNext = bbly1-ffpadY-extH;
													int endYNext = ttry1+bbpadY-extH;
													
													if //((startX+startY != startXOld+startYOld))// && 
													( startX+startY != startXNext+startYNext) //)
														sname << join(mname.str(),j) << range(d->getShiftAmount()-1, 0) << " & " << sname.str();
												}
										}
									else // only multiplication
										{
											mname.str("");
											mname << "pxy" << i << j;
											vhdl << tab << declare(mname.str(), multW+multH+2) << " <= " << use(xname.str())<<range(multW,0) << " * " << use(yname.str())<<range(multH,0) << ";" << endl;
											sname.str("");
											if (d->getShiftOut() == NULL) // concatenate the entire partial product
												{
//													setCycle(connected);
													//sname << zg(wInX-(blx1-fpadX-extW)+wInY-(bly1-fpadY-extH)-2, 0) << " & " << use(mname.str()) << range(multH-fpadY+multW-fpadX-1, bpadX+bpadY)<< " & " << zg(trx1-extW,0) << " & " << zg(try1-extH,0) <<  ";" << endl;
													int nrZeros = wInX-(blx1-extW)-fpadX + wInY-(bly1-extH)-fpadY-2;
													if (nrZeros < 0)
														sname << zg(wInX-(blx1-fpadX-extW)+wInY-(bly1-fpadY-extH)-2, 0) << " & " << use(mname.str()) << range(multH-fpadY+multW-fpadX-1, 0)<< " & " << zg(trx1-extW-bpadX,0) << " & " << zg(try1-extH-bpadY,0) <<  ";" << endl;
													else
														sname << zg(wInX-(blx1-extW) + wInY-(bly1-extH)-2, 0) << " & " << use(mname.str()) << range(multW + multH-1, 0) << " & " << zg(trx1-extW-bpadX,0) << " & " << zg(try1-extH-bpadY,0) <<  ";" << endl;
	
												}
											else // concatenate only the lower portion of the partial product
												{
													/* check if the next element has as 17-bit shift wr to this one */
													int ttrx1,ttry1,bblx1,bbly1, ffpadX, ffpadY, bbpadX, bbpadY;
													d->getShiftOut()->getTopRightCorner(ttrx1, ttry1);
													d->getShiftOut()->getBottomLeftCorner(bblx1, bbly1);
									
													ffpadX = bblx1-wInX-extW+1;
													ffpadX = (ffpadX<0)?0:ffpadX;
									
													ffpadY = bbly1-wInY-extH+1;
													ffpadY = (ffpadY<0)?0:ffpadY;
									
													bbpadX = extW-ttrx1;
													bbpadX = (bbpadX<0)?0:bbpadX;
									
													bbpadY = extH-ttry1;
													bbpadY = (bbpadY<0)?0:bbpadY;
									
									
													int startXNext = bblx1-ffpadX-extW;
													int endXNext = ttrx1+bbpadX-extW;
													int startYNext = bbly1-ffpadY-extH;
													int endYNext = ttry1+bbpadY-extH;
													
													if //((startX+startY != startXOld+startYOld))// && 
													( startX+startY != startXNext+startYNext) //)
														sname << use(mname.str()) << range(d->getShiftAmount()-1, bpadX+bpadY) << " & ";
														
													sname <<  zg(trx1-extW,0) << " & " << zg(try1-extH,0) << ";" << endl;
												}
										}
				
									// erase d from the tempc buffer to avoid handleing it twice
									for (int k=i+1; k<nrDSPs; k++)
										{
											if ((tempc[k] != NULL) && (tempc[k] == d))
												{
													//~ cout << "tempc[" << k << "] deleted" << endl;
													tempc[k] = NULL;
													break;
												}
										}
				
				
									startXOld = startX;
									startYOld = startY;
									
									d = d->getShiftOut();
									j++;
								}
								
							if (outside) continue;	
							sname.seekp(ios_base::beg);
							sname << tab << declare(join("addOpDSP", nrOp),wInX+wInY) << " <= " << sname.str();
							vhdl << sname.str();
							nrOp++;		
						}
		
				return nrOp;
			}
		else // the target is Stratix
		{
				int boundDSPs;  				// number of bound DSPs in a group
				DSP** addOps;					// addition operands bound to a certain DSP
	
				for (int i=0; i<nrDSPs; i++)
					if (tempc[i] != NULL){
							cout << "At DSP#"<< i+1 << " tempc["<<i<<"]" << endl; 
							tempc[i]->getTopRightCorner(trx1, try1);
							tempc[i]->getBottomLeftCorner(blx1, bly1);
							convertCoordinates(trx1, try1, blx1, bly1);

							int trx2,try2,blx2,bly2;
							tempc[i]->getTopRightCorner(trx2, try2);
							tempc[i]->getBottomLeftCorner(blx2, bly2);
							convertCoordinatesKeepNeg(trx2, try2, blx2, bly2);
							multW = blx2-trx2+1;
							multH = bly2-try2+1;
							
							int maxX = blx1;
							int maxY = bly1;
			
							setCycle(0);
							vhdl << tab << declare(join("x",i,"_0"), multW, true) << " <= " << zg(blx2-blx1) << " & X"<<range(blx1, trx1) << ";" << endl;
							vhdl << tab << declare(join("y",i,"_0"), multH, true) << " <= " << zg(bly2-bly1) << " & Y"<<range(bly1, try1) << ";" << endl;

							boundDSPs = tempc[i]->getNumberOfAdders();
							int ext = 0;        // the number of carry bits of the addtion accumulation. 
							if (boundDSPs > 0){ // need to traverse the addition operands list and perform addtion
								ext = (boundDSPs>1)?2:1;
								REPORT(INFO, "boundDSPs = " << boundDSPs);
								nextCycle();
								mname.str("");
								mname << "mult_" << i << "_0";
								vhdl << tab << declare(mname.str(), multW+multH, true) << " <= " << join("x",i,"_0") << " * " << join("y",i,"_0") << ";" << endl;
								addOps = tempc[i]->getAdditionOperands();
			
								/* At most 4 operands */
								for (int j=0; j<3; j++)
									if (addOps[j] == NULL)
										cout << "addOps["<< j << "]=NULL" << endl;
									else
										cout << "addOps["<< j << "]=not null" << endl;
			
								for (int j=0; j<boundDSPs; j++){
										cout << "j = " << j << endl;
										// erase addOps[j] from the tempc buffer to avoid handleing it twice
										for (int k=i+1; k<nrDSPs; k++){
											if ((tempc[k] != NULL) && (tempc[k] == addOps[j])){
												REPORT( DETAILED, "tempc[" << k << "] deleted");
												tempc[k] = NULL;
												break;
											}
										}
				
										addOps[j]->getTopRightCorner(trx1, try1);
										addOps[j]->getBottomLeftCorner(blx1, bly1);
										convertCoordinates(trx1, try1, blx1, bly1);
										multW = addOps[j]->getMaxMultiplierWidth();
										multH = addOps[j]->getMaxMultiplierHeight();
										
										addOps[j]->getTopRightCorner(trx2, try2);
										addOps[j]->getBottomLeftCorner(blx2, bly2);
										convertCoordinatesKeepNeg(trx2, try2, blx2, bly2);

										multW = blx2-trx2+1;
										multH = bly2-try2+1;

										if ( bly2+blx2>maxX+maxY){
											maxY = bly1;
											maxX = blx1;
										}
										
										setCycle(0); ////////////////////////////////////
										vhdl << tab << declare(join("x",i,"_",j+1), multW, true) << " <= " << zg(blx2-blx1) << " & X"<<range(blx1, trx1) << ";" << endl;
										vhdl << tab << declare(join("y",i,"_",j+1), multH, true) << " <= " << zg(bly2-bly1) << " & Y"<<range(bly1, try1) << ";" << endl;


										nextCycle(); ////////////////////////////////////
										vhdl << tab << declare(join("mult_",i,"_",j+1), multW+multH, true) << " <= " << join("x",i,"_",j+1) << " * " << join("y",i,"_",j+1) << ";" << endl;
								}
			
								nextCycle();
								vhdl << tab << declare(join("addDSP", nrOp), multW+multH+ext, true) << " <= ";
			
								for (int j=0; j<boundDSPs; j++){
									vhdl << "(" << zg(ext,0) << " & " << join("mult_",i,"_",j) << ") + "; 
								}
								vhdl << "(" << zg(ext,0) << " & " << join("mult_",i,"_",boundDSPs) << ");" << endl; 
							}else{ // multiply the two terms and you're done
								nextCycle();
								vhdl << tab << declare(join("addDSP", nrOp), multW+multH, true) << " <= " << join("x",i,"_0") << " * " << join("y",i,"_0") << ";" << endl;
							}
							vhdl << tab << declare(join("addOpDSP", nrOp), wInX+wInY) << " <= " << zg(wInX+wInY-(maxX+maxY+2+ext),0) << " & " << join("addDSP", nrOp)<<range( (maxX+maxY-(trx2+try2)+2) +ext-1,0)  << " & " << zg(trx2+try2,0) << ";" << endl;
							nrOp++;
					}
				return nrOp;
		}
	}

	int IntTilingMult::multiplicationInSlices(DSP** config)
	{
		//~ cout<<"Incepe"<<endl;
		int partitions=0;
		int **mat;
		int n,m;
		int count=1;
		vector<SoftDSP*> softConfig;
		//~ n=wInX + 2* getExtraWidth();
		//~ m=wInY + 2* getExtraHeight();
		n=vn;
		m=vm;
		//~ cout<<"width "<<n<<"height "<<m<<endl;
		mat = new int*[m];
		for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
		
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				//~ cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
		
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
		partitions=0;
		
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					{
						if(mat[i][j]==0)
							{
								int ver =0;
								int ii=i,jj=j;
								while(ver<6&&(ii<m-1||jj<n-1))
									{
										if(ver<3)
											{
												if(ver==0||ver==1)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=2;							
													}
					
												if(ver==0||ver==2)
													jj++;
					
												if(jj>n-1)
													{
														jj=n-1;
														ver=1;
													}
					
												for(int k=ii,q=jj;k>i-1&&(ver==0||ver==2);k--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=1;
															else
																ver=3;
															jj--;
														}
						
												for(int k=ii,q=jj;q>j-1&&(ver==0||ver==1);q--)
													if(mat[k][q]!=0)
														{
															if (ver==0)
																ver=2;
															else
																ver=3;
															ii--;
														}
											}
										else
											{
												if(ver==3||ver==5)
													jj++;
					
												if(jj>n-1)
													{
														jj=n-1;
														ver=4;
													}
						
												if(ver==3||ver==4)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=5;							
													}
					
												for(int k=ii,q=jj;q>j-1&&(ver==3||ver==4);q--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=5;
															else
																ver=6;
															ii--;
														}
						
												for(int k=ii,q=jj;k>i-1&&(ver==3||ver==5);k--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=4;
															else
																ver=6;
															jj--;
														}
						
												if(ver==5&&jj==n-1)
													ver=6;
												if(ver==4&&ii==m-1)
													ver=6;
											}
									}
				
								int nj,ni,njj,nii;
								int extH = getExtraHeight();
								int extW = getExtraWidth();
				
				
								if( j >= n-extW || jj < extW || i >= m-extH || ii < extH)
									{
										REPORT(DETAILED, "Partition number "<<count<<" is totally out of the real multiplication bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")");
									}
								else
									{
										if( j < extW )
											nj = extW ;
										else
											nj = j;
										if( jj >= n-extW )
											njj = n-extW -1;
										else
											njj = jj;
					
										if( i < extH )
											ni = extH;
										else
											ni = i;
										if( ii >= m-extH )
											nii = m-extH-1;
										else
											nii = ii;
										
										setCycle(0);
										SoftDSP *sdsp = new SoftDSP(wInX-nj+2*extW, nii+1, wInX-njj-2+2*extW, ni-1);
										softConfig.push_back(sdsp);
//										target_->setUseHardMultipliers(false);
										LogicIntMultiplier* mult =  new LogicIntMultiplier(target_, njj-nj+1, nii-ni+1); //unsigned
										ostringstream cname;
										cname << mult->getName() << "_" << partitions;
										mult->changeName(cname.str());
										oplist.push_back(mult);
										// TODO: compute width of x and y + corretc range for X and Y
										vhdl << tab << declare(join("x_",partitions), njj-nj+1, true) << " <= X" << range(wInX-nj-1+extW, wInX-njj-1+extW) << ";" << endl;
										inPortMap(mult, "X", join("x_",partitions));
										vhdl << tab << declare(join("y_",partitions), nii-ni+1, true) << " <= Y" << range(nii-extH, ni-extH) << ";" << endl;
										inPortMap(mult, "Y", join("y_",partitions));
					
										outPortMap(mult, "R", join("result", partitions));
					
										vhdl << instance(mult, join("Mult", partitions));

										syncCycleFromSignal(join("result", partitions));
										
										vhdl << tab << declare(join("addOpSlice", partitions), wInX+wInY) << " <= " << zg(wInX+wInY-(wInX-nj-1+extW+nii-extH)-2, 0) << " & " << join("result", partitions) << " & " << zg(wInX-njj-1+extW+ni-extH, 0) << ";" << endl;
										REPORT(DETAILED, "Partition number "<<count<<" with bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<") has now bounds ("<<nj<<" , "<<ni<<" , "<<njj<<" , "<<nii<<")");
										REPORT(DETAILED, "partitions " << partitions << " @ cycle " << getCurrentCycle());
										
										partitions++;
									}
				
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
				
							}
					}
	
			}
	
		printConfiguration(bestConfig, softConfig);
		//de verificat
		
		//cout<<"Count "<<count<<" Partitions "<<partitions<<endl;
		
		//partitions =count -partitions;
		 
		
		//~ char af;
		//~ int afi;
		//~ for(int i=0;i<m;i++)
		//~ {
		//~ for(int j=0;j<n;j++)
		//~ {
		//~ if(mat[i][j]<10)
		//~ afi=mat[i][j];
		//~ else
		//~ afi=mat[i][j]+7;
		//~ af=(int)afi+48;
		//~ cout<<af;
		//~ }
		//~ cout<<endl;
		//~ }
	
		//~ cout<<"gata"<<endl;
		
		//~ for(int i=0;i<nrDSPs;i++)
		//~ {
			//~ if(config[i]->getShiftOut()!=NULL)
			//~ {
				//~ config[i]->getShiftOut()->getTopRightCorner(n,m);
				//~ cout<<"There is a shift connection from DSP# "<<i<<" to DSP with coordinates "<<n<<","<<m<<endl;
			//~ }
		//~ }
	
		for(int ii=0;ii<m;ii++)
			delete[](mat[ii]);
	
		delete[] (mat);
		
		return partitions;
	}	

		void IntTilingMult::outputVHDL(std::ostream& o, std::string name) {
		licence(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.std_logic_arith.all;" << endl;
		if  ( (target_->getID() == "Virtex4") ||
		      (target_->getID() == "Virtex5") ||
		      (target_->getID() == "Spartan3"))  // then the target is a Xilinx FPGA
			{
		o << "use ieee.std_logic_signed.all;" << endl;
		}else{
		o << "use ieee.std_logic_unsigned.all;" << endl;
		}
		
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
	
	void IntTilingMult::printConfiguration(DSP** configuration, vector<SoftDSP*> softDSPs){
		ofstream fig;
		fig.open ("tiling.fig", ios::trunc);
		fig << "#FIG 3.2  Produced by xfig version 3.2.5a" << endl;
		fig << "Landscape" << endl;
		fig << "Center" << endl;
		fig << "Metric" << endl;
		fig << "A4      " << endl;
		fig << "100.00" << endl;
		fig << "Single" << endl;
		fig << "-2" << endl;
		fig << "1200 2" << endl;
	
		if (configuration!=NULL){
			int i=0;
			int xB,xT,yB,yT;
			for(i=0; i<nrDSPs; i++){
				configuration[i]->getTopRightCorner(xT,yT);
				configuration[i]->getBottomLeftCorner(xB,yB);
				REPORT(DETAILED, "HARD DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB);
				fig << " 2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5 " << endl;
				fig << "	  " << (-xB+getExtraWidth()-1)*45 << " " << (yT-getExtraHeight())*45 
				         << " " << (-xT+getExtraWidth())*45 << " " << (yT-getExtraHeight())*45 
				         << " " << (-xT+getExtraWidth())*45 << " " << (yB-getExtraHeight()+1)*45 
				         << " " << (-xB+getExtraWidth()-1)*45 << " " << (yB-getExtraHeight()+1)*45 
				         << " " << (-xB+getExtraWidth()-1)*45 << " " << (yT-getExtraHeight())*45 << endl;
				
				int dpx = (-(xT+xB)/2+getExtraWidth()-1)*45;
				int dpy = ((yB+yT)/2-getExtraHeight())*45;         
				REPORT(DETAILED, "x="<<dpx<<" y="<<dpy);
				fig << " 4 1 0 50 -1 0 12 0.0000 4 195 630 "<<dpx<<" "<<dpy<<" DSP"<<i<<"\\001" << endl;

				//annotations
				fig << " 4 1 0 50 -1 0 7 0.0000 4 195 630 "<<(-xB+getExtraWidth())*45<<" "<<-45<<" "<<xB-getExtraWidth()<<"\\001" << endl;
				fig << " 4 1 0 50 -1 0 7 0.0000 4 195 630 "<<(-xT+getExtraWidth()-1)*45<<" "<<-45<<" "<<xT-getExtraWidth()<<"\\001" << endl;
				fig << " 4 0 0 50 -1 0 7 0.0000 4 195 630 "<<45<<" "<<(yT-getExtraHeight()+2)*45<<" "<<yT-getExtraHeight()<<"\\001" << endl;
				fig << " 4 0 0 50 -1 0 7 0.0000 4 195 630 "<<45<<" "<<(yB-getExtraHeight()+1)*45<<" "<<yB-getExtraHeight()<<"\\001" << endl;
			}
		}

		int xB,xT,yB,yT;
		for (unsigned k=0; k < softDSPs.size(); k++){
			softDSPs[k]->trim(vnme, vmme);
			softDSPs[k]->getTopRightCorner(xT,yT);
			softDSPs[k]->getBottomLeftCorner(xB,yB);
			
			REPORT(DETAILED,  "SOFT DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB);
			fig << " 2 2 0 1 0 7 50 -1 19 0.000 0 0 -1 0 0 5 " << endl;
			fig << "	  " << (-xB+getExtraWidth()-1)*45 << " " << (yT-getExtraHeight())*45 
			<< " " << (-xT+getExtraWidth())*45 << " " << (yT-getExtraHeight())*45 
			<< " " << (-xT+getExtraWidth())*45 << " " << (yB-getExtraHeight()+1)*45 
			<< " " << (-xB+getExtraWidth()-1)*45 << " " << (yB-getExtraHeight()+1)*45 
			<< " " << (-xB+getExtraWidth()-1)*45 << " " << (yT-getExtraHeight())*45 << endl;
			int dpx = (-(xT+xB)/2+getExtraWidth()-1)*45;
			int dpy = ((yB+yT)/2-getExtraHeight())*45;         
			REPORT(DETAILED,  "x="<<dpx<<" y="<<dpy);
			fig << " 4 1 0 50 -1 0 12 0.0000 4 195 630 "<<dpx<<" "<<dpy<<" M"<<k<<"\\001" << endl;
			
			xT++;
			yT++;
			xB++;
			yB++;
			int tmp;
			tmp=xT;
			xT=xB;
			xB=tmp;
			tmp=yT;
			yT=yB;
			yB=tmp;
			
			//annotations
			fig << " 4 1 0 50 -1 0 7 0.0000 4 195 630 "<<(-xB+getExtraWidth())*45<<" "<<-45<<" "<<xB-getExtraWidth()<<"\\001" << endl;
			fig << " 4 1 0 50 -1 0 7 0.0000 4 195 630 "<<(-xT+getExtraWidth()-1)*45<<" "<<-45<<" "<<xT-getExtraWidth()<<"\\001" << endl;
			fig << " 4 0 0 50 -1 0 7 0.0000 4 195 630 "<<45<<" "<<(yT-getExtraHeight()+2)*45<<" "<<yT-getExtraHeight()<<"\\001" << endl;
			fig << " 4 0 0 50 -1 0 7 0.0000 4 195 630 "<<45<<" "<<(yB-getExtraHeight()+1)*45<<" "<<yB-getExtraHeight()<<"\\001" << endl;
		}
		
		fig << "		2 2 1 1 0 7 50 -1 -1 4.000 0 0 -1 0 0 5" << endl;
		fig << "	  " << (-wInX)*45 << " " << 0 
         << " " << 0 << " " << 0  
         << " " << 0 << " " << (wInY)*45 
         << " " << (-wInX)*45 << " " << (wInY)*45 
         << " " << (-wInX)*45 << " " << 0 << endl;

		//big X and Y
		fig << " 4 1 0 50 -1 0 16 0.0000 4 195 630 "<<-wInX/2*45<<" "<<-3*45<<" X\\001" << endl;
		fig << " 4 0 0 50 -1 0 16 0.0000 4 195 630 "<<3*45<<" "<<wInY/2*45<<" Y\\001" << endl;

		
		fig.close();
		int hh,ww,aa;
		getTarget()->getDSPWidths(hh,ww);
		aa = hh * ww;
		REPORT(INFO, " DSP area is: " << aa ); 
		
		
		if (configuration!=NULL){
			int i=0;
			int xB,xT,yB,yT;
			int underutilized = 0;
			int sunderutilized = 0;

			
			for(i=0; i<nrDSPs; i++){
				configuration[i]->getTopRightCorner(xT,yT);
				configuration[i]->getBottomLeftCorner(xB,yB);
				
				yT=yT-getExtraHeight();
				yB=yB-getExtraHeight();
				xT=xT-getExtraWidth();
				xB=xB-getExtraWidth();

				xB = min(xB, wInX-1);
				yB = min(yB, wInY-1);
				
				if ( float((xB-xT+1)*(yB-yT+1))< float(aa)) {
					if ( float((xB-xT+1)*(yB-yT+1))/float(aa) < 0.5 ){
						REPORT(INFO, "HARD DSP is SEVERELY under-utilized, just " << float((xB-xT+1)*(yB-yT+1))/float(aa) << "%");
						sunderutilized++;
					}else{
						REPORT(INFO, "HARD DSP utilized " << float((xB-xT+1)*(yB-yT+1))/float(aa) << "%");
						underutilized++;
					}
				}
			}
			REPORT(INFO, "********************************************************************************");
			REPORT(INFO, "*      underutilized = " << underutilized  + sunderutilized);
			REPORT(INFO, "*      suggested ratio = " << (float(nrDSPs)*ratio - float(sunderutilized)) / float(nrDSPs));
			REPORT(INFO, "********************************************************************************");
		}
		
		
	}

}

