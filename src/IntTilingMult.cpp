/*
 * Tilling Multiplier for FloPoCo
 *
 * Author :  Sebastian Banescu, Radu Tudoran, Bogdan Pasca
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
#include <typeinfo>
#include <math.h>
#include <string.h>
#include <limits.h>

#include <gmp.h>


#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntMultiplier.hpp"
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

	IntTilingMult:: IntTilingMult(Target* target, int wInX, int wInY, float ratio) :
		Operator(target), wInX(wInX), wInY(wInY), wOut(wInX + wInY),ratio(ratio){
 
		ostringstream name;

		name <<"IntMultiplier_"<<wInX<<"_"<<wInY;
		setName(name.str());
		
		setCopyrightString("Sebastian Banescu , Radu Tudoran, Bogdan Pasca 2009-2010");
	
		addInput ("X", wInX);
		addInput ("Y", wInY);
		addOutput("R", wOut); /* wOut = wInX + wInY */
		
		nrDSPs = estimateDSPs();
		cout<<"Estimated DSPs:= "<<nrDSPs <<endl;
		int x,y;
		target->getDSPWidths(x,y);
		
		cout<<"Width of DSP is := "<<x<<" Height of DSP is:="<<y<<endl;
		cout<<"Extra width:= "<<getExtraWidth()<<" \nExtra height:="<<getExtraHeight()<<endl;
		
		vn=wInX + 2* getExtraWidth();
		vm=wInY + 2* getExtraHeight();

		//~ float tempDist =	 (movePercentage  * getExtraWidth() * getExtraWidth()) /4.0 + (movePercentage *getExtraHeight() * getExtraHeight()) /4.0;
		float tempDist =	0;
		maxDist2Move = (int) ( sqrt(tempDist) );
		//~ if(maxDist2Move==0)
		//~ maxDist2Move=1;
		cout<<"maxDist2Move:= "<<maxDist2Move<<endl;
		
	
		
		
		
		
		
		
		
		//~ for(int i=0;i<nrDSPs;i++)
		//~ {
		//~ if(globalConfig[i]!=NULL)
		//~ {
		//~ if(globalConfig[i]->getShiftIn()!=NULL)
		//~ {
		//~ for(int j=0;j<nrDSPs;j++)
		//~ if(globalConfig[j]==globalConfig[i]->getShiftIn())
		//~ cout<<"Exists a bind of type Out from dsp "<<i+1<<" to DSPul"<<j+1<<endl;
		//~ }
		//~ if(globalConfig[i]->getShiftOut()!=NULL)
		//~ {
		//~ for(int j=0;j<nrDSPs;j++)
		//~ if(globalConfig[j]==globalConfig[i]->getShiftOut())
		//~ cout<<"Exists a bind of type Out from dsp "<<i+1<<" to DSPul"<<j+1<<endl;
		//~ }
		//~ }
		//~ }
		
		
		
		runAlgorithm();
	 
	 
	 
	 
	 
		//~ ///experiment(bogdan)  ./flopoco -target=Virtex5 IntTilingMult 58 58 0.44
		//~ tempc= new DSP*[nrDSPs];
		//~ for(int ii=0;ii<nrDSPs;ii++)
		//~ tempc[ii]= new DSP();
		//~ int n,m;
		//~ int count=1;
	
		//~ n=vn;
		//~ m=vm;
		//~ mat = new int*[m];
		//~ for(int i=0;i<m;i++)
		//~ {
		//~ mat[i] = new int [n];
		//~ for(int j=0;j<n;j++)
		//~ mat[i][j]=0;
		//~ }
	
	
	
		//~ initTiling(bestConfig,nrDSPs);
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(5,4);
		//~ globalConfig[0]->setBottomLeftCorner(21,27);
		//~ globalConfig[1]->setTopRightCorner(22,4);
		//~ globalConfig[1]->setBottomLeftCorner(38,27);
		//~ globalConfig[2]->setTopRightCorner(39,4);
		//~ globalConfig[2]->setBottomLeftCorner(55,27);
		//~ globalConfig[3]->setTopRightCorner(56,4);
		//~ globalConfig[3]->setBottomLeftCorner(72,27);
		//~ globalConfig[4]->setTopRightCorner(0,28);
		//~ globalConfig[4]->setBottomLeftCorner(16,51);
		//~ globalConfig[5]->setTopRightCorner(17,28);
		//~ globalConfig[5]->setBottomLeftCorner(33,51);
		//~ globalConfig[6]->setTopRightCorner(34,28);
		//~ globalConfig[6]->setBottomLeftCorner(50,51);
		//~ globalConfig[7]->setTopRightCorner(51,28);
		//~ globalConfig[7]->setBottomLeftCorner(67,51);
	
		
	
	
		//~ bestCost = 333222;
		//~ compareCost();
	 
		//~ display(bestConfig);
	
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(7,5);
		//~ globalConfig[0]->setBottomLeftCorner(23,28);
		//~ globalConfig[1]->setTopRightCorner(24,5);
		//~ globalConfig[1]->setBottomLeftCorner(40,28);
		//~ globalConfig[2]->setTopRightCorner(41,5);
		//~ globalConfig[2]->setBottomLeftCorner(64,21);
		//~ globalConfig[3]->setTopRightCorner(41,22);
		//~ globalConfig[3]->setBottomLeftCorner(64,38);
		//~ globalConfig[4]->setTopRightCorner(7,29);
		//~ globalConfig[4]->setBottomLeftCorner(30,45);
		//~ globalConfig[5]->setTopRightCorner(7,46);
		//~ globalConfig[5]->setBottomLeftCorner(30,62);
		//~ globalConfig[6]->setTopRightCorner(31,39);
		//~ globalConfig[6]->setBottomLeftCorner(47,62);
		//~ globalConfig[7]->setTopRightCorner(48,39);
		//~ globalConfig[7]->setBottomLeftCorner(64,62);
	
		//~ display(globalConfig);
	
	
		//~ cout<<checkFarness(globalConfig,2);
		//~ move(globalConfig,1);
		//~ display(globalConfig);
	
		//~ compareCost();
	
		//experiment pentru ./flopoco -target=Virtex5 IntTilingMult 18 34 0.35
		//~ tempc= new DSP*[nrDSPs];
		//~ for(int ii=0;ii<nrDSPs;ii++)
		//~ tempc[ii]= new DSP();
		//~ int n,m;
		//~ int count=1;
	
		//~ n=vn;
		//~ m=vm;
		//~ mat = new int*[m];
		//~ for(int i=0;i<m;i++)
		//~ {
		//~ mat[i] = new int [n];
		//~ for(int j=0;j<n;j++)
		//~ mat[i][j]=0;
		//~ }
	
	
	
		//~ initTiling(bestConfig,nrDSPs);
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(6,4);
		//~ globalConfig[0]->setBottomLeftCorner(22,27);
		//~ display(globalConfig);
	
		//~ bestCost = 333222;
		//~ compareCost();
	
	
	
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(4,4);
		//~ globalConfig[0]->setBottomLeftCorner(20,27);
		//~ display(globalConfig);
	
	
	
		//~ compareCost();
	 
		///////////////////////////////////////////
	
	
		//~ globalConfig[1]->setTopRightCorner(1,17);
		//~ globalConfig[1]->setBottomLeftCorner(17,33);
		//~ globalConfig[2]->setTopRightCorner(21,-12);
		//~ globalConfig[2]->setBottomLeftCorner(37,4);

		//~ display(globalConfig);
	
		//~ maxDist2Move=2;
	
		//~ //cout<<endl<<checkFarness(globalConfig,2) <<endl;
	
		//~ //cout<<move(globalConfig,2)<<endl;
		//~ cout<<move(globalConfig,0)<<endl;
		//~ display(globalConfig);
		//~ replace(globalConfig,1);
		//~ replace(globalConfig,2);
		//~ display(globalConfig);
		//~ cout<<move(globalConfig,2)<<endl;
	
		//~ display(globalConfig);
	
	



		//~ initTiling(bestConfig,nrDSPs);
		//~ bestCost=367.455;
	
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(1,2);
		//~ globalConfig[0]->setBottomLeftCorner(9,10);
		//~ globalConfig[1]->setTopRightCorner(9,13);
		//~ globalConfig[1]->setBottomLeftCorner(17,21);
		//~ globalConfig[2]->setTopRightCorner(18,4);
		//~ globalConfig[2]->setBottomLeftCorner(26,12);
	
	
		//~ cout<<"Initial configuration"<<endl<<endl;
		//~ display(globalConfig);
	
		//~ maxDist2Move=2;
	
		//~ cout<<endl<<checkFarness(globalConfig,1) <<endl;
	
		//~ compareCost();
	
	
		
	
		cout<<"Estimated DSPs:= "<<nrDSPs <<endl;
		target->getDSPWidths(x,y);
		cout<<"Width of DSP is := "<<x<<" Height of DSP is:="<<y<<endl;
		cout<<"Extra width:= "<<getExtraWidth()<<" \nExtra height:="<<getExtraHeight()<<endl;
		cout<<"maxDist2Move:= "<<maxDist2Move<<endl;
	
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(2,15);
		//~ globalConfig[0]->setBottomLeftCorner(18,31);
		//~ globalConfig[1]->setTopRightCorner(19,15);
		//~ globalConfig[1]->setBottomLeftCorner(35,31);
	
		//~ int t;
		//~ cout<<"cost of partitions is "<<partitionOfGridSlices(globalConfig,t);
		//~ cout<<"Cost of obtained best is "<<computeCost(globalConfig)<<endl;
		//~ display(globalConfig);
	
		//~ cout<<endl<<endl;
	
		//~ initTiling(globalConfig,nrDSPs);
	
		//~ globalConfig[0]->setTopRightCorner(2,14);
		//~ globalConfig[0]->setBottomLeftCorner(18,30);
		//~ globalConfig[1]->setTopRightCorner(19,14);
		//~ globalConfig[1]->setBottomLeftCorner(35,30);
	
		//~ cout<<"cost of partitions is "<<partitionOfGridSlices(globalConfig,t);
		//~ cout<<"Cost of obtained best is "<<computeCost(globalConfig)<<endl;
		//~ display(globalConfig);
		
		
		

		/*
	
	
		  initTiling(globalConfig, 4);
		  for (int i=0; i<19; i++)
		  move(globalConfig, 2);
		  //replace(globalConfig, 1);
		  //move(globalConfig, 0);
		  //int x, y;
	
		  for (int i=0; i<4; i++)
		  {
		  globalConfig[i]->getTopRightCorner(x,y);
		  cout << "DSP block #" << i << " is positioned at top-right: (" << x << ", " << y << ")";
		  globalConfig[i]->getBottomLeftCorner(x,y);
		  cout << "and bottom-left: (" << x << ", " << y << ")" << endl;
		  }
		*/        

	
	}
	
	IntTilingMult::~IntTilingMult() {
	}




	void IntTilingMult::runAlgorithm()
	{
		counterfirst =0 ;
		int n,m;

		//~ n=wInX + 2* getExtraWidth();
		//~ m=wInY + 2* getExtraHeight();

		n=vn;
		m=vm;
	
		mat = new int*[m];
		for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		
		tempc= new DSP*[nrDSPs];
		for(int ii=0;ii<nrDSPs;ii++)
			tempc[ii]= new DSP();

		/*we will try the algorithm with 2 values of nrDSPs	
		  One will be the estimated value(nrDSPs) and the second one will be nrDSPs-1	
		*/
		rot = new bool[nrDSPs];
		for(int i =0;i<nrDSPs;i++)
			rot[i]=false;
	
	
	
		initTiling(globalConfig,nrDSPs);	
			
		//this will initialize the bestConfig with the first configuration
		bestCost = FLT_MAX ;
		cout<<"Max score is"<<bestCost<<endl;
		//bestConfig = (DSP**)malloc(nrDSPs * sizeof(DSP*));
		bestConfig = new DSP*[nrDSPs];
		for(int i=0;i<nrDSPs;i++)
			bestConfig[i]= new DSP();
		compareCost();
		cout<<"New best score is"<<bestCost<<endl;
	
		//display(bestConfig);
	
	
		//the best configuration should be consider initially the first one. So the bestConfig parameter will be initialized with global config and hence the bestCost will be initialized with the first cost
	
	
		tilingAlgorithm(nrDSPs-1,nrDSPs-1,false,nrDSPs-1);
	
		display(bestConfig);
		cout<<"Best cost is "<<bestCost<<endl;
	
	
	
	

		//~ // After all configurations with the nrDSPs number of DSPs were evaluated then a new search is carryed with one DSP less
		//~ // After the initialization of the new configuration with nrDSPs-1, the cost must be evaluated and confrunted with the best score obtained so far.
	
		//~ if(nrDSPs-1>0)
		//~ {
		
	
		//~ for(int i =0;i<nrDSPs;i++)
		//~ rot[i]=false;
		
		//~ initTiling(globalConfig,nrDSPs -1);	
		//~ compareCost();
		//~ tilingAlgorithm(nrDSPs-2,nrDSPs-2,false);
		//~ }	
	
	
		for(int ii=0;ii<m;ii++)
			delete[](mat[ii]);
	
		delete[] (mat);
	
		for(int ii=0;ii<nrDSPs;ii++)
			free(tempc[ii]);
	
		delete[] (tempc);
	
	}

	/** The movement of the DSP blocks with values belonging to their widths and heights still needs to be done. Now it only runs with one type of move on one of the directions, which is not ok for the cases when the DSPs are not squares.
	 */
	void IntTilingMult::tilingAlgorithm(int i, int n,bool repl,int lastMovedDSP)
	{

		if(i==n)
			{	
				if(repl==true) // if previous DSPs were moved this one needs to recompute all positions 
					{
						//~ cout<<" Pas 4_1 "<<i<<endl;
						if(replace(globalConfig,i)) // repostioned the DSP
							{
								compareCost();
								rot[i]=false;
								tilingAlgorithm(i,n,false,lastMovedDSP);	
							}
						else // could not reposition the DSP in the bounds of the tiling board
							{
								rot[i]=false;
								if( lastMovedDSP>=0) // go one level up the backtracking stack
									tilingAlgorithm(lastMovedDSP,n,false, lastMovedDSP);
							}
					}
				else // the last DSP is being moved on the tiling board
					{
	
						if(move(globalConfig,i,5,5)) // successfuly moved the last block
							{
								//~ cout<<" Pas 1_1 "<<i<<endl;
								compareCost();
								tilingAlgorithm(i,n,repl,i);		//repl should be false
							}
						//~ else
						//~ if(move(globalConfig,i,DSPh,DSPw))
						//~ {
						//~ cout<<" Pas 1_1 "<<i<<endl;
						//~ compareCost();
						//~ tilingAlgorithm(i,n,repl,i);		//repl should be false
						//~ }
						else // could not find a position for the last block
							{
								if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() ))
									{ // if the DSP was not rotated and is not sqare then roteate it
										//display(globalConfig);
										//~ cout<<" Pas 2_1 "<<i<<endl;
										globalConfig[i]->rotate();
										//display(globalConfig);
										rot[i]=true;
										if(replace(globalConfig,i)) // the DSP could be repositioned
											{
												//display(globalConfig);
												compareCost();
												tilingAlgorithm(i,n,repl,i);		//repl should be false
											}
										else // go to the previous block 
											{
												if(i-1>=0)
													tilingAlgorithm(i-1,n,false,i);
											}
									}
								else // the DSP was either rotated already or is square
									{
										//~ cout<<" Pas 3_1 "<<i<<endl;
										if(i-1>=0)
											tilingAlgorithm(i-1,n,repl,i);		//repl should be false
									}
							}
					}
			}
		else // we are not at the last DSP
			{
				if(repl==true) // the previuos DSPs were successfuly repositioned
					{
						//~ cout<<" Pas 4_2 "<<i<<endl;
						if(replace(globalConfig,i)) // the current DSP was successfuly repositioned
							{
								rot[i]=false;
								tilingAlgorithm(i+1,n,repl, lastMovedDSP);
							}
						else // the current DSP could not be repositioned
							{// go to the DSP block that was moved (not repostioned) the last time
								rot[i]=false;
								if( lastMovedDSP>=0) 
									tilingAlgorithm( lastMovedDSP,n,false, lastMovedDSP);
							}
			
		
					}
				else // the folling DSP could not be moved or repositioned 
					{	
						//~ if(i==0)
						//~ display(globalConfig);
						if(move(globalConfig,i,17,17)) // the current DSP was successfuly moved
							{
						
								if(i==0){
									counterfirst++;
				
				
				
				
									if(counterfirst%100==0)
										cout<<counterfirst<<"DSP #1 has made 100 steps!"<<endl;
									//~ display(globalConfig);
									//~ cout<<endl<<endl<<endl;
				
								}
								tilingAlgorithm(i+1,n,true,i);
							}
						//~ else
						//~ if(move(globalConfig,i,DSPh,DSPw))
						//~ {
						
						//~ if(i==0){
						//~ counterfirst++;
				
				
				
				
						//~ if(counterfirst%100==0)
						//~ cout<<counterfirst<<"DSP #1 has made 100 steps!"<<endl;
						//~ display(globalConfig);
						//~ cout<<endl<<endl<<endl;
				
						//~ }
						//~ tilingAlgorithm(i+1,n,true,i);
						//~ }
						else // the current DSP was not moved successfuly
							{
								if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() ))
									{// if the DSP was not rotated and is not sqare then roteate it
										//~ cout<<" Pas 2_2 "<<i<<endl;
										globalConfig[i]->rotate();
										if(replace(globalConfig,i)) // the current DSP was successfuly repositioned
											{
												rot[i]=true;
												tilingAlgorithm(i+1,n,true,i);
											}
										else // the current DSP was not successfuly repositioned
											{
												if(i-1>=0)
													tilingAlgorithm(i-1,n,repl,i);
											}
									}
								else // the DSP is either square or has been already rotated
									{
										//~ cout<<" Pas 3_2 "<<i<<endl;
										if(i-1>=0)
											tilingAlgorithm(i-1,n,repl,i);		//repl should be false
									}
							}
					}
			}

	
	}


	bool IntTilingMult::compareOccupation(DSP** config)
	{
		int totalSize = wInX * wInY;
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
	
		int nmew = n-getExtraWidth();
		int ew = getExtraWidth();
		int mmeh = m - getExtraHeight() ;
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

	void IntTilingMult::fillMatrix(int **&matrix,int lw,int lh,int topleftX,int topleftY,int botomrightX,int botomrightY,int value)
	{
		for(int j=topleftX;j<=botomrightX;j++)
			{
				for(int i=topleftY;i<=botomrightY;i++)
					{
						if(j>-1&&i>-1&&i<lh&&j<lw)
							matrix[i][j]=value;
					}
			}
	
	}


	void IntTilingMult::display(DSP** config)
	{
	
	
		int **mat;
		int n,m;
		int count=1;
		n=wInX + 2* getExtraWidth();
		m=wInY + 2* getExtraHeight();
		cout<<"real width"<<wInX<<"real height"<<wInY<<endl;
		cout<<"width "<<n<<"height "<<m<<endl;
		mat = new int*[m];
	
		int nmew = n-getExtraWidth();
		int ew = getExtraWidth();
		int mmeh = m - getExtraHeight() ;
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
				cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
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
										cout<<"Partition number "<<count<<" is totally out of the real multiplication bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
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
										cout<<"Partition number "<<count<<" with bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<" has now bounds ("<<nj<<" , "<<ni<<" , "<<njj<<" , "<<nii<<")"<<endl;
									}
						
								cout<<j<<" "<<i<<" "<<jj<<" "<<ii<<endl;
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
						
							}
					}
			
			}
		
		
		
		char af;
		int afi;
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
		
		for(int ii=0;ii<m;ii++)
		   delete [] (mat[ii]);
	
		delete[] (mat);
		
		
	
	}



	int IntTilingMult::partitionOfGridSlices(DSP** config,int &partitions)
	{
		//~ cout<<"Incepe"<<endl;
		int costSlice=0;
	
		int n,m;
		int count=1;
		n=vn;
		m=vm;
	
		int nmew = n-getExtraWidth();
		int ew = getExtraWidth();
		int mmeh = m - getExtraHeight() ;
		int eh = getExtraHeight();
		int nj,ni,njj,nii;
	
		
		//~ cout<<"width "<<n<<"height "<<m<<endl;
	
		for(int i=0;i<m;i++)
			{
			
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
		//partitions = count;
		partitions = 0;
	
		//~ cout<<"Partea 2"<<endl;
		
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
												//~ if(count==3)
												//~ cout<<"P0  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
								
								
							
												for(int k=ii,q=jj;k>i-1&&(ver==0||ver==2);k--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=1;
															else
																ver=3;
															jj--;
														}
												//~ if(count==3)
												//~ {
												//~ cout<<"P1   ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl; 
												//~ }
									
												for(int k=ii,q=jj;q>j-1&&(ver==0||ver==1);q--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=2;
															else
																ver=3;
															ii--;
														}
												//~ if(count==3)
												//~ {cout<<"P2  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
												//~ cout<<endl;
												//~ }
							
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
								
												//~ if(count==3)
												//~ cout<<"P3  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
								
												if(ver==3||ver==4)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=5;							
													}
							
								
												//~ if(count==3)
												//~ cout<<"P3  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;

								
												for(int k=ii,q=jj;q>j-1&&(ver==3||ver==4);q--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=5;
															else
																ver=6;
															ii--;
														}
												//~ if(count==3)
												//~ {
												//~ cout<<"P4   ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl; 
												//~ }
								
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
										
												//~ if(count==3)
												//~ {cout<<"P5  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
												//~ cout<<endl;
												//~ }
							
								
											}
									}
						
								//~ cout<<count<<endl;
						
					
												
						
								if( j>=nmew || jj< ew || i >= mmeh || ii < eh)
									{
							
									}
								else
									{
										if( j < ew )
											nj = ew ;
										else
											nj = j;
										if( jj >=nmew )
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
							
										costSlice += getTarget()->getIntMultiplierCost(njj-nj+1,nii-ni+1);
							
							
														
									}
						
						
						
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
						
							}
					}
			
			}
		
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
	
		//~ for(int ii=0;ii<m;ii++)
		//~ delete[](mat[ii]);
	
		//~ delete[] (mat);
	
		return costSlice;
	}	



	int IntTilingMult::bindDSPs4Virtex(DSP** &config)
	{
		int nrOfUsedDSPs=0;
		for(int i=0;i<nrDSPs;i++)
			if(config[i]!=NULL)
				{
					nrOfUsedDSPs++;
				}
		DSP* ref;
	
		sortDSPs(config);
		
			
		int itx,ity,jtx,jty,ibx,iby;//,jbx,jby;
		
		for(int i=0;i<nrDSPs;i++)
			{
				if(config[i]!=NULL)
					{
						ref=config[i];
						bool ver=true;
						while(ver==true&&ref->getShiftOut()==NULL)
							{
								ver=false;
								ref->getTopRightCorner(itx,ity);
								ref->getBottomLeftCorner(ibx,iby);
					
					
								for(int j=0;j<nrDSPs&&ver==false;j++)
									{
										if(config[j]!=NULL &&j!=i)
											{
												config[j]->getTopRightCorner(jtx,jty);
												//cout<<"Now considering taking(in left) dsp nr. "<<i<<" with tx:="<<itx<<" ty:="<<ity<<" bx:="<<ibx<<"by:="<<iby<<" with dsp nr. "<<j<<" with tx:="<<jtx<<" ty:="<<jty<<endl;
												//config[j]->getBottomLeftCorner(jbx,jby);
												if(jtx==ibx+1&&jty==ity&&ref->getMaxMultiplierWidth()==config[j]->getShiftAmount()&&config[j]->getShiftIn()==NULL)
													{
														ver=true;
														ref->setShiftOut(config[j]);
														config[j]->setShiftIn(ref);
														nrOfUsedDSPs--;
														ref=config[j];
													}
											}
									}
					
								for(int j=0;j<nrDSPs&&ver==false;j++)
									{
										if(config[j]!=NULL &&j!=i)
											{
												config[j]->getTopRightCorner(jtx,jty);
												//cout<<"Now considering taking(down) dsp nr. "<<i<<" with tx:="<<itx<<" ty:="<<ity<<" bx:="<<ibx<<"by:="<<iby<<" with dsp nr. "<<j<<" with tx:="<<jtx<<" ty:="<<jty<<endl;
												//config[j]->getBottomLeftCorner(jbx,jby);
												if(iby+1==jty&&itx==jtx&&ref->getMaxMultiplierHeight()==config[j]->getShiftAmount()&&config[j]->getShiftIn()==NULL)
													{
														ver=true;
														ref->setShiftOut(config[j]);
														config[j]->setShiftIn(ref);
														nrOfUsedDSPs--;
														ref=config[j];								
													}
							
											}						
									}					
							}
				
					}
			
			}
	
		return nrOfUsedDSPs;
	
	}

	void IntTilingMult::sortDSPs(DSP** &config)
	{
		int ix,iy,jx,jy;
		DSP* temp;
		for(int i=0;i<nrDSPs-1;i++)
			{		
				for(int j=i+1;j<nrDSPs;j++)
					{
						config[i]->getTopRightCorner(ix,iy);
						config[j]->getTopRightCorner(jx,jy);
						if(iy>jy)
							{
								temp=config[i];
								config[i]=config[j];
								config[j]=temp;				
							}
						else
							if(iy==jy)
								{
									if(ix>jx)
										{
											temp=config[i];
											config[i]=config[j];
											config[j]=temp;						
										}
								}
					}
			}
	
	}

	int IntTilingMult::bindDSPs(DSP** &config)
	{
		if  ( (getTarget()->getID() == "Virtex4") ||
		      (getTarget()->getID() == "Virtex5") ||
		      (getTarget()->getID() == "Spartan3"))  // then the target is a Xilinx FPGA
			{
				return bindDSPs4Virtex(config);
			}
		else // the target is Stratix
			{
				return bindDSPs4Stratix(config);
			}
	}


	void IntTilingMult::compareCost()
	{
		//~ cout<<"Inta la cost"<<endl;
	
		//~ DSP** tempc;
		
		//~ tempc= new DSP*[nrDSPs];
		//~ for(int ii=0;ii<nrDSPs;ii++)
		//~ tempc[ii]= new DSP();
		
		//display(globalConfig);
	
		//memcpy(tempc,globalConfig,sizeof(DSP*) *nrDSPs );
		for(int ii=0;ii<nrDSPs;ii++)
			memcpy(tempc[ii],globalConfig[ii],sizeof(DSP) );
	
		//display(tempc);
		//~ cout<<"intra la display cost"<<endl;
	
		float temp = computeCost(tempc);
	
		//~ cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
	
		if(temp < bestCost)
			{
				//~ cout<<"Costul e mai bun la cel curent!Schimba"<<endl;
				cout<<"Interchange! Score for current is"<<temp<<" and current best is"<<bestCost<<endl;
				//~ int c1X,c2X,c1Y,c2Y;
				//~ int n=wInX + 2* getExtraWidth();
				//~ tempc[0]->getTopRightCorner(c1X,c1Y);
				//~ tempc[0]->getBottomLeftCorner(c2X,c2Y);
				//~ cout<<"DSP #"<<1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				//~ c1X=n-c1X-1;
				//~ c2X=n-c2X-1;
				//~ //cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
				//~ cout<<"matrix DSP #"<<1<<"has topleft ("<<c2X<<","<<c1Y<<") and botomright ("<<c1X<<","<<c2Y<<")"<<endl;
		
		
		
				bestCost=temp;
				//memcpy(bestConfig,tempc,sizeof(DSP*) *nrDSPs );	
				for(int ii=0;ii<nrDSPs;ii++)
					memcpy(bestConfig[ii],tempc[ii],sizeof(DSP) );
				//display(bestConfig);
			}
		else
			if(temp == bestCost )
				{
					//cout<<"Cost egal!!!"<<endl;
					//cout<<"Rezult compare is"<<compareOccupation(tempc)<<endl;
					//cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
					//display(bestConfig);
					if(compareOccupation(tempc)==true)
						{
							cout<<"Interchange for equal cost. Now best has cost "<<temp<<endl;
							bestCost=temp;
			
							for(int ii=0;ii<nrDSPs;ii++)
								memcpy(bestConfig[ii],tempc[ii],sizeof(DSP) );
							//	display(bestConfig);
						}
				}
	
	
		//~ for(int ii=0;ii<nrDSPs;ii++)
		//~ free(tempc[ii]);
	
		//~ delete[] (tempc);
	
	}



	float IntTilingMult::computeCost(DSP** &config)
	{
	
		float const scale=100.0;
		float acc=0.0;
		float costDSP,costLUT;
		costDSP = ( (1.0+scale) - scale * ratio );
		costLUT = ( (1.0+scale) - scale * (1-ratio) ) /  ((float) getTarget()->getEquivalenceSliceDSP() );
	
		//~ cout<<"Cost of a DSP is "<<costDSP<<endl<<"Cost of a Slice is "<<costLUT<<endl;
	
		int nrOfUsedDSPs=0;
	
		for(int i=0;i<nrDSPs;i++)
			if(config[i]!=NULL)
				{
					acc+=costDSP;
					nrOfUsedDSPs++;
				}
	
		
		//~ cout<<"Number of used DSP blocks is "<<nrOfUsedDSPs<<endl;
		
		int partitions;
		float LUTs4Multiplication =  partitionOfGridSlices(config,partitions);
	
		//~ cout<<"Number of slices 4 multiplication of the rest is "<<LUTs4Multiplication<<endl;
		
		acc =((float)nrOfUsedDSPs)*costDSP + costLUT * LUTs4Multiplication;
		
		//~ cout<<"Number of partitions for LUTs is "<<partitions<<endl;
		nrOfUsedDSPs = bindDSPs(config);
		//~ cout<<"Number of operands coming from DSPs is "<<nrOfUsedDSPs<<endl;
		
	
		float LUTs4NAdder=((float)getTarget()->getIntNAdderCost(wInX + wInY,nrOfUsedDSPs+partitions) );
	
		//~ cout<<"LUTs used for last "<<nrOfUsedDSPs+partitions<<" adder are"<<LUTs4NAdder<<endl;
		
		acc +=  LUTs4NAdder* costLUT;	
		
		//~ Substracting the cost of different additions that can be done in DSPs(Virtex) or the cost of a DSP if they can be combined(Altera)
		
		
		return acc;
			
	}

	int IntTilingMult::estimateDSPs()
	{
		float t1,t2;
		int Xd, Yd; //the dimension of the multiplier on X and on Y
		int multX, multY;
		bool fitMultiplicaication = false;

		getTarget()->getDSPWidths(Xd,Yd);
		int maxDSP= getTarget()->getNumberOfDSPs();
	
		t1= ((float) wInX) / ((float) Xd);
		t2= ((float) wInY) / ((float) Yd);
		
		multX = int (ceil(t1));
		multY = int (ceil(t2));	
		
		if(maxDSP >= (multX * multY) ){
			fitMultiplicaication = true;
			maxDSP = ((int)t1) * ((int)t2); //set the maximum number of DSPs to the multiplication size
		}
			
		if (ratio == 1){
			if (! fitMultiplicaication){
				REPORT(INFO, "Warning!!! The number of existing DSPs on this FPGA is not enough to cover the whole multiplication!");
			}else{
				REPORT(INFO, "Warning!!! The minimum number of DSP that is neccessary to cover the whole addition will be used!");
			}
			
			return maxDSP;
		}else{	
			float temp = ( float(getTarget()->getIntMultiplierCost(wInX,wInY)) * ratio)  /   ((1.-ratio)*float(getTarget()->getEquivalenceSliceDSP())) ;
			int i_tmp = int(ceil(temp));
	
			if(i_tmp > maxDSP){
				if (fitMultiplicaication){
					REPORT(INFO, "Warning!!! The number of estimated DSPs with respect with this ratio of preference is grather then the needed number of DSPs to perform this multiplication!");
				}else{
					REPORT(INFO, "Warning!!! The number of estimated DSPs with respect with this ratio of preference is grather then the total number of DSPs that exist on this board!");
				}
				i_tmp = maxDSP;
			}
	
			return i_tmp ;
		}
	}

	int  IntTilingMult::getExtraHeight()
	{
		int x,y;	
		getTarget()->getDSPWidths(x,  y);
		float temp = ratio * 0.75 * ((float) y);
		return ((int)temp);
	}

	
	int  IntTilingMult::getExtraWidth()
	{
		int x,y;	
		getTarget()->getDSPWidths(x,y);
		float temp = ratio * 0.75 * ((float) x);
		return ((int)temp);
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

	bool IntTilingMult::checkOverlap(DSP** config, int index)
	{	
		int xtr1, ytr1, xbl1, ybl1, xtr2, ytr2, xbl2, ybl2;
		//int nrdsp = sizeof(config)/sizeof(DSP*);
		config[index]->getTopRightCorner(xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
	

	
		bool a1 ;
		bool a2 ;
		bool a3 ;
		bool a4 ;
		bool a5 ;
		bool a6 ;
				
		bool b1 ;
		bool b2 ;
		bool b3 ;
		bool b4 ;
		bool b5 ;
		bool b6 ;
	
	
		if(verbose)
			cout << tab << tab << "checkOverlap: ref is block #" << index << ". Top-right is at (" << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
	
		for (int i=0; i<index; i++)
			if (config[i] != NULL)
				{
					config[i]->getTopRightCorner(xtr2, ytr2);
					config[i]->getBottomLeftCorner(xbl2, ybl2);
			
					if ( 
						 (((xtr2 >= xtr1) && (ytr2 > ybl1)) || 	// config[index] is above and to the right of config[i]
						  (xbl1 < xtr2) 							// config[index] is to the right of config[i]
						  )
						  )
						return true;
			
				
					if(verbose)
						cout << tab << tab << "checkOverlap: comparing with block #" << i << ". Top-right is at (" << xtr2 << ", " << ytr2 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
			
					//~ if (((xtr2 <= xbl1) && (ytr2 <= ybl1) && (xtr2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the upper and/or right part(s) of config[i]
					//~ ((xbl2 <= xbl1) && (ybl2 <= ybl1) && (xbl2 >= xtr1) && (ybl2 >= ytr1)) || // config[index] overlaps the bottom and/or left part(s) of config[i]
					//~ ((xtr2 <= xbl1) && (ybl2 <= ybl1) && (xtr2 >= xtr1) && (ybl2 >= ytr1)) || // config[index] overlaps the upper and/or left part(s) of config[i]
					//~ ((xbl2 <= xbl1) && (ytr2 <= ybl1) && (xbl2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the bottom and/or right part(s) of config[i]
					//~ ((xbl2 >= xbl1) && (ybl2 <= ybl1) && (ytr2 >= ytr1) && (xtr2 <= xtr1)) || // config[index] overlaps the center part of config[i]
				
					//~ ((xtr1 <= xbl2) && (ytr1 <= ybl2) && (xtr1 >= xtr2) && (ytr1 >= ytr2)) || // config[i] overlaps the upper and/or right part(s) of config[index]
					//~ ((xbl1 <= xbl2) && (ybl1 <= ybl2) && (xbl1 >= xtr2) && (ybl1 >= ytr2)) || // config[i] overlaps the bottom and/or left part(s) of config[index]
					//~ ((xtr1 <= xbl2) && (ybl1 <= ybl2) && (xtr1 >= xtr2) && (ybl1 >= ytr2)) || // config[i] overlaps the upper and/or left part(s) of config[index]
					//~ ((xbl1 <= xbl2) && (ytr1 <= ybl2) && (xbl1 >= xtr2) && (ytr1 >= ytr2)) || // config[i] overlaps the bottom and/or right part(s) of config[index]
					//~ ((xbl1 >= xbl2) && (ybl1 <= ybl2) && (ytr1 >= ytr2) && (xtr1 <= xtr2))    // config[i] overlaps the center part of config[index]
					//~ )
					//~ return true;
					
					// the optimisation of the above if
					a1 = (xtr2 <= xbl1);
					a2 = (xtr2 >= xtr1);
					a3 = (xbl2 <= xbl1);
					a4 = (xbl2 >= xtr1);
					a5 = (xbl2 >= xbl1);
					a6 = (xtr1 >= xtr2);
				
					b1 = (ytr2 <= ybl1);
					b2 = (ytr2 >= ytr1);
					b3 = (ybl2 <= ybl1);
					b4 = (ybl2 >= ytr1);
					b5 = (ytr1 >= ytr2);
					b6 = (ybl1 <= ybl2);
				
					if ((((a1 && a2)||(a3 && a4)) && ((b1 && b2)||(b3 && b4))) || 
					    (((a4 && a6)||(a5 && a1)) && ((b6 && b1)||(b4 && b5))) || 
					    (((a5 && b3) && ( b2 && a6)) || ((a3 && b6) && (b5 && a2))))
						return true;
				
				
			
			
				}	
		if(verbose)
			cout << tab << tab << "checkOverlap: return false" << endl;	
		return false;
	}


	/**
		There is one case that is not resolved w=yet. For DSP with widths different then height the algorithm should move the dsp with both values
	*/

	bool IntTilingMult::move(DSP** config, int index,int w,int h)
	{
		int xtr1, ytr1, xbl1, ybl1;
		//~ int w, h;
		//target->getDSPWidths(w,h);
		//~ w= config[index]->getMaxMultiplierWidth();
		//~ h= config[index]->getMaxMultiplierHeight();
	
	
	
	
		config[index]->getTopRightCorner(xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
	
		if(verbose)
			cout << tab << "replace : DSP #" << index << " width is " << w << ", height is " << h << endl; 
		/*
		 * Initialized in the constructor:
		 * vn = wInX+2*getExtraWidth(); width of the tiling grid including extensions
		 * vm = wInY+2*getExtraHeight(); height of the tiling grid including extensions
		 * */
		int gw = wInX+getExtraWidth()-1; // the left limit of the tiling grid
		int gh = wInY+getExtraHeight()-1; // the bottom limit of the tiling grid
		int exh = getExtraHeight();
		int exw = getExtraWidth();
	
		if ((xtr1 < 0) || (ytr1 < 0) || (xbl1 > vn-1) || (ybl1 > vm-1))
			{// then the DSP block is placed outside the bounds 		
	
				if(index==0)
					{
						do{
							// move down one unit
							ytr1++;
							ybl1++;
			
							if (ytr1 > gh) // the DSP block has reached the bottom limit of the tiling grid
								{
									// move to top of grid and one unit to the left 
									xtr1++;
									xbl1++;
					
									if (!( (getTarget()->getID() == "Virtex4") ||
										   (getTarget()->getID() == "Virtex5") ||
										   (getTarget()->getID() == "Spartan3")) )  // then the target is NOT Xilinx FPGA 
										{
											ytr1 = exh ;
											ybl1 = ytr1 + h-1;	
										}
									else
										{
											ytr1 = exh - h+1;
											ybl1 = ytr1 + h-1;
										}
					
									if (xtr1 > gw) // the DSP block has reached the left limit of the tiling grid
										return false;
								}						
							config[index]->setTopRightCorner(xtr1, ytr1);
							config[index]->setBottomLeftCorner(xbl1, ybl1);
						}while (checkOverlap(config, index));
					}
				else
					{
						do{
							// move down one unit
							if(ytr1 < exh)
								{
									ytr1++;
									ybl1++;
								}
							else
								{
									ytr1 += h;
									ybl1 += h;
								}
							if (ytr1 > gh) // the DSP block has reached the bottom limit of the tiling grid
								{
									// move to top of grid and one unit to the left 
									if(xtr1< exw)
										{
											xtr1+=1;
											xbl1+=1;
										}
									else
										{
											xtr1+=w;
											xbl1+=w;
										}
					
									if (!( (getTarget()->getID() == "Virtex4") ||
										  ( getTarget()->getID() == "Virtex5") ||
										  ( getTarget()->getID() == "Spartan3")))  // then the target is NOT Xilinx FPGA 
										{
											ytr1 = exh ;
											ybl1 = ytr1 + h-1;	
										}
									else
										{
											ytr1 = exh - h+1;
											ybl1 = ytr1 + h-1;
										}
					
				
									if (xtr1 > gw) // the DSP block has reached the left limit of the tiling grid
										return false;
								}						
							config[index]->setTopRightCorner(xtr1, ytr1);
							config[index]->setBottomLeftCorner(xbl1, ybl1);
						}while (checkOverlap(config, index));
					}
		
			}
		else
			{
				if(index==0)
					{
						do{
							// move down one unit
				
							ytr1++;
							ybl1++;
			
							if (ybl1 > vm-1) // the DSP block has reached the bottom limit of the tiling grid
								{
									// move to top of grid and one unit to the left 
									xtr1++;
									xbl1++;
										
					
					
									if (!( (getTarget()->getID() == "Virtex4") ||
										  ( getTarget()->getID() == "Virtex5") ||
										  ( getTarget()->getID() == "Spartan3")))  // then the target is NOT Xilinx FPGA 
										{
											ytr1 = exh ;
											ybl1 = ytr1 + h-1;	
										}
									else
										{
											ytr1 = 0;
											ybl1 = h-1;
										}
				
									if (xbl1 > vn-1) // the DSP block has reached the left limit of the tiling grid
										return false;
								}			
			
							config[index]->setTopRightCorner(xtr1, ytr1);
							config[index]->setBottomLeftCorner(xbl1, ybl1);
						}while (checkOverlap(config, index));
					}
				else
					{
						do{
							// move down one unit
							if(ytr1 < exh)
								{
									ytr1++;
									ybl1++;
								}
							else
								{
									ytr1+=h;
									ybl1+=h;
								}
							if (ybl1 > vm-1) // the DSP block has reached the bottom limit of the tiling grid
								{
									// move to top of grid and one unit to the left 
									if(xtr1<exw)
										{
											xtr1+=1;
											xbl1+=1;
										}
									else
										{
											xtr1+=w;
											xbl1+=w;
										}
									
					
									if (!( (getTarget()->getID() == "Virtex4") ||
										  ( getTarget()->getID() == "Virtex5") ||
										  ( getTarget()->getID() == "Spartan3")))  // then the target is NOT Xilinx FPGA 
										{
											ytr1 = exh ;
											ybl1 = ytr1 + h-1;	
										}
									else
										{
											ytr1 = 0;
											ybl1 = h-1;
										}
				
									if (xbl1 > vn-1) // the DSP block has reached the left limit of the tiling grid
										return false;
								}			
			
							config[index]->setTopRightCorner(xtr1, ytr1);
							config[index]->setBottomLeftCorner(xbl1, ybl1);
						}while (checkOverlap(config, index));
					}
		
		
			}
		// set the current position of the DSP block within the tiling grid
		config[index]->setTopRightCorner(xtr1, ytr1);
		config[index]->setBottomLeftCorner(xbl1, ybl1);
	
		int f = checkFarness(config,index);	
		if (f == 0)
			return true;
		else if (f == 1)
			return false;
		else if (f  == 2)
			{
				// move to top of grid and one unit to the left 
				xtr1++;
				xbl1++;
				ytr1 = exh - h+1;
				ybl1 = ytr1 + h-1;
		
				if (xbl1 > vn-1) // the DSP block has reached the left limit of the tiling grid
					return false;
				config[index]->setTopRightCorner(xtr1, ytr1);
				config[index]->setBottomLeftCorner(xbl1, ybl1);
				move(config, index,w,h);
			}
		else if (f == 3)
			{
				move(config, index,w,h);	
			}
		return false;
	}

	bool IntTilingMult::replace(DSP** config, int index)
	{
		int xtr1, ytr1, xbl1, ybl1;
		int w, h;
		//target->getDSPWidths(w,h);
		w= config[index]->getMaxMultiplierWidth();
		h= config[index]->getMaxMultiplierHeight();
	
		xtr1 = ytr1 = 0;
		xbl1 = w-1;
		ybl1 = h-1;
	
		config[index]->setTopRightCorner(xtr1, ytr1);
		config[index]->setBottomLeftCorner(xbl1, ybl1);
	
		if(verbose)
			cout << tab << "replace : DSP width is " << w << ", height is " << h << endl; 
		// try yo place the DSP block inside the extended region of the tiling grid
		while (checkOverlap(config, index))
			{
				// move down one unit
				ytr1++;
				ybl1++;
		
				if(verbose)
					cout << tab << "replace : moved DSP one unit down." << endl;
				if (ybl1 > vm) // the DSP block has reached the bottom limit of the tiling grid
					{
						// move to top of grid and one unit to the left 
						xtr1++;
						xbl1++;
						ytr1 = 0;
						ybl1 = h-1;
						if(verbose)
							cout << tab << "replace : moved DSP up and one unit left." << endl;
					}			
		
				config[index]->setTopRightCorner(xtr1, ytr1);
				config[index]->setBottomLeftCorner(xbl1, ybl1);
				if(verbose)
					cout << tab << "replace : Top-right is at ( " << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
			}
	
		if (xbl1 > vn)
			return false;
	 
		return true;
	
		//~ if (xbl1 > vn) // it was not possible to place the DSP inside the extended grid
		//~ { // then try to place it anywhere around the tiling grid such that it covers at leat one 1x1 square
		//~ xtr1 = getExtraWidth() - w+1;
		//~ ytr1 = getExtraHeight() - h+1;
		//~ xbl1 = xtr1 + w-1;
		//~ ybl1 = ytr1 + h-1;
		
		//~ config[index]->setTopRightCorner(xtr1, ytr1);
		//~ config[index]->setBottomLeftCorner(xbl1, ybl1);
		
		//~ if(verbose)
		//~ cout << tab << "replace : DSP width is " << w << ", height is " << h << endl; 
			
		//~ int bLimit = wInY+getExtraHeight();
		//~ int exh = getExtraHeight();
		
		//~ while (checkOverlap(config, index))
		//~ {
		//~ // move down one unit
		//~ ytr1++;
		//~ ybl1++;
			
		//~ if(verbose)
		//~ cout << tab << "replace : moved DSP one unit down." << endl;
		//~ if (ytr1 > bLimit) // the DSP block has reached the bottom limit of the tiling grid
		//~ {
		//~ // move to top of grid and one unit to the left 
		//~ xtr1++;
		//~ //if (xtr1 > wInX+getExtraWidth()) // the DSP block has reached the left limit of the tiling grid
		//~ xbl1++;
		//~ ytr1 = exh - h+1;
		//~ ybl1 = ytr1 + h-1;
				
		//~ if(verbose)
		//~ cout << tab << "replace : moved DSP up and one unit left." << endl;
		//~ }			
			
		//~ config[index]->setTopRightCorner(xtr1, ytr1);
		//~ config[index]->setBottomLeftCorner(xbl1, ybl1);
		//~ if(verbose)
		//~ cout << tab << "replace : Top-right is at ( " << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
		//~ }
	
		//~ }
		//~ // set the current position of the DSP block within the tiling grid
		//~ config[index]->setTopRightCorner(xtr1, ytr1);
		//~ config[index]->setBottomLeftCorner(xbl1, ybl1);
	}

	void IntTilingMult::initTiling(DSP** &config, int dspCount)
	{
		config = new DSP*[nrDSPs];
	
	
	
		for (int i=0; i<nrDSPs; i++)
			{
				config[i] = NULL;
			}
		if (!( ( getTarget()->getID() == "Virtex4") ||
			   ( getTarget()->getID() == "Virtex5") ||
			   ( getTarget()->getID() == "Spartan3")))  // then the target is NOT Xilinx FPGA 
			{
				for (int i=0; i<dspCount; i++)
					{
						if(verbose)
							cout << "initTiling : iteration #" << i << endl; 
						config[i] = getTarget()->createDSP();						
						replace(config, i);
					}
			}
		else
			{
				int w,h;	
				getTarget()->getDSPWidths(w,h);
			
				//~ config[0] = target->createDSP();
				//~ config[0]->setTopRightCorner(getExtraWidth(), getExtraHeight());
				//~ config[0]->setBottomLeftCorner(getExtraWidth() +w-1 , getExtraHeight() + h-1);
		
				int xtr,ytr,xbl,ybl;
	
				xtr = 	getExtraWidth();
				ytr = 	getExtraHeight();
				xbl =	getExtraWidth() +w-1;
				ybl = 	getExtraHeight() + h-1;
		
				for(int i=0;i<dspCount;i++)
					{
						config[i] = getTarget()->createDSP();
						
						if( ytr<vm&& ybl<vm)
							{
								config[i]->setTopRightCorner(xtr, ytr);
								config[i]->setBottomLeftCorner(xbl , ybl);
								ytr = ybl+1;
								ybl = ybl + h;
							}
						else
							{
								xtr = xbl+1;
								xbl = xbl + w;
								ytr = 	getExtraHeight();
								ybl = 	getExtraHeight() + h-1;
								config[i]->setTopRightCorner(xtr, ytr);
								config[i]->setBottomLeftCorner(xbl , ybl);
								ytr = ybl+1;
								ybl = ybl + h;
			
							}
		
					}		
		
				//~ for (int i=0; i<dspCount; i++)
				//~ {
				//~ if(verbose)
				//~ cout << "initTiling : iteration #" << i << endl; 
				//~ config[i] = target->createDSP();
				//~ replace(config, i);
				//~ }	
			}
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
					getTarget()->getDSPWidths(wx,wy);
				
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

	int IntTilingMult::multiplicationInDSPs(DSP** config)
	{
		int nrOp = 1;		 			// number of resulting adder operands
		int trx1, try1, blx1, bly1; 	// coordinates of the two corners of a DSP block 
		int fpadX, fpadY, bpadX, bpadY;	// zero padding on both axis
		int extW, extH;					// extra width and height of the tiling grid
		int multW, multH; 				// width and height of the multiplier the DSP block is using
		ostringstream xname, yname, mname, cname, sname;
		DSP** tempc = new DSP*[nrDSPs];	// buffer that with hold a copy of the global configuration
			
		memcpy(tempc, config, sizeof(DSP*) * nrDSPs );
	
		if ( ( getTarget()->getID() == "Virtex4") ||
			 ( getTarget()->getID() == "Virtex5") ||
			 ( getTarget()->getID() == "Spartan3"))  // then the target is A Xilinx FPGA 
			{
				for (int i=0; i<nrDSPs; i++)
					if (tempc[i] != NULL)
						{
							cout << "At DSP#"<< i+1 << " tempc["<<i<<"]" << endl; 
							DSP* d = tempc[i];
							int j=0;
							int connected = 0;
			
							while (d != NULL)
								{
									connected++;
									d = d->getShiftOut();
								}
			
							d = tempc[i];
			
							while (d != NULL)
								{
									d->getTopRightCorner(trx1, try1);
									d->getBottomLeftCorner(blx1, bly1);
									extW = getExtraWidth();
									extH = getExtraHeight();
									fpadX = blx1-wInX-extW+1;
									fpadX = (fpadX<0)?0:fpadX;
									cout << "fpadX = " << fpadX << endl;
									fpadY = bly1-wInY-extH+1;
									cout << "fpadY = " << fpadY << endl;
									fpadY = (fpadY<0)?0:fpadY;
									bpadX = extW-trx1;
									bpadX = (bpadX<0)?0:bpadX;
									bpadY = extH-try1;
									bpadY = (bpadY<0)?0:bpadY;
									multW = d->getMaxMultiplierWidth();
									multH = d->getMaxMultiplierHeight();
			
									setCycle(0);
									xname.str("");
									xname << "x" << i << "_" << j;
									vhdl << tab << declare(xname.str(), multW) << " <= " << zg(fpadX,0) << " & " << "X" << range(blx1-fpadX, trx1+bpadX) << " & " << zg(bpadX,0) << ";" << endl;
									yname.str("");
									yname << "y" << i << "_" << j;
									vhdl << tab << declare(yname.str(), multH) << " <= " << zg(fpadY,0) << " & " << "Y" << range(bly1-fpadY, try1+bpadY) << " & " << zg(bpadY,0) << ";" << endl;
				
									if (d->getShiftIn() != NULL) // multiply accumulate
										{
											mname.str("");
											mname << "pxy" << i;
											cname.str("");
											cname << "txy" << i << j;
											setCycle(j);
											vhdl << tab << declare(cname.str(), multW+multH) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
											vhdl << tab << declare(join(mname.str(),j), multW+multH+1) << " <= (\"0\" & " << use(cname.str()) << ") + " << use(join(mname.str(), j-1)) << range(multW+multH-1, d->getShiftAmount()) << ";" << endl;	
											if (d->getShiftOut() == NULL) // concatenate the entire partial product
												{
													setCycle(connected);
													sname.seekp(ios_base::beg);
													sname << zg(wInX+wInY+extW+extH-blx1-bly1-3, 0) << " & " << use(join(mname.str(),j)) << " & " << sname.str();
												}
											else // concatenate only the lower portion of the partial product
												{
													setCycle(connected);
													sname.seekp(ios_base::beg);
													sname << use(join(mname.str(),j)) << range(d->getShiftAmount()-1, 0) << " & " << sname.str();
												}
										}
									else // only multiplication
										{
											mname.str("");
											mname << "pxy" << i << j;
											vhdl << tab << declare(mname.str(), multW+multH) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
											sname.str("");
											if (d->getShiftOut() == NULL) // concatenate the entire partial product
												{
													setCycle(connected);
													sname << zg(wInX+wInY+extW+extH-blx1-bly1-2, 0) << " & " << use(mname.str()) << range(multH+multW-1, bpadX+bpadY)<< " & " << zg(trx1-extW,0) << " & " << zg(try1-extH,0) <<  ";" << endl;
												}
											else // concatenate only the lower portion of the partial product
												{
													setCycle(connected);
													sname << use(mname.str()) << range(d->getShiftAmount()-1, bpadX+bpadY) << " & " << zg(trx1-extW,0) << " & " << zg(try1-extH,0) << ";" << endl;
												}
										}
				
									// erase d from the tempc buffer to avoid handleing it twice
									for (int k=i+1; k<nrDSPs; k++)
										{
											if ((tempc[k] != NULL) && (tempc[k] == d))
												{
													cout << "tempc[" << k << "] deleted" << endl;
													tempc[k] = NULL;
													break;
												}
										}
				
				
									d = d->getShiftOut();
									j++;
								}	
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
				int mPadX, mPadY, minPad;		// minimum padding of addition operands
	
				for (int i=0; i<nrDSPs; i++)
					if (tempc[i] != NULL)
						{
							mPadX = INT_MIN;
							mPadY = INT_MIN;
							cout << "At DSP#"<< i+1 << " tempc["<<i<<"]" << endl; 
							tempc[i]->getTopRightCorner(trx1, try1);
							tempc[i]->getBottomLeftCorner(blx1, bly1);
							fpadX = blx1-wInX-getExtraWidth()+1;
							fpadX = (fpadX<0)?0:fpadX;
							fpadY = bly1-wInY-getExtraHeight()+1;
							fpadY = (fpadY<0)?0:fpadY;
							bpadX = getExtraWidth()-trx1;
							bpadX = (bpadX<0)?0:bpadX;
							mPadX = (bpadX>mPadX)?bpadX:mPadX;
							bpadY = getExtraHeight()-try1;
							bpadY = (bpadY<0)?0:bpadY;
							mPadY = (bpadY>mPadY)?bpadY:mPadY;
							minPad = bpadY+bpadX;
							multW = tempc[i]->getMaxMultiplierWidth();
							multH = tempc[i]->getMaxMultiplierHeight();
			
							setCycle(0);
							xname.str("");
							xname << "x" << i << "_0";
							vhdl << tab << declare(xname.str(), multW, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadX,0) << " & " << "X" << range(blx1-fpadX, trx1+bpadX) << " & " << zg(bpadX,0) << ";" << endl;
							yname.str("");
							yname << "y" << i << "_0";
							vhdl << tab << declare(yname.str(), multH, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadY,0) << " & " << "Y" << range(bly1-fpadY, try1+bpadY) << " & " << zg(bpadY,0) << ";" << endl;
			
							boundDSPs = tempc[i]->getNumberOfAdders();
							int ext = 0; // the number of carry bits of the addtion
							if (boundDSPs > 0) // need to traverse the addition operands list and perform addtion
								{
									ext = (boundDSPs>1)?2:1;
									cout << "boundDSPs = " << boundDSPs << endl;
									nextCycle();
									mname.str("");
									mname << "mult_" << i << "_0";
									vhdl << tab << declare(mname.str(), multW+multH, true, Signal::registeredWithAsyncReset) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
				
									addOps = tempc[i]->getAdditionOperands();
				
									//
									for (int j=0; j<3; j++)
										if (addOps[j] == NULL)
											cout << "addOps["<< j << "]=NULL" << endl;
										else
											cout << "addOps["<< j << "]=not null" << endl;
				
				
					
									//
				
				
									for (int j=0; j<boundDSPs; j++)
										{
											cout << "j = " << j << endl;
											// erase addOps[j] from the tempc buffer to avoid handleing it twice
											for (int k=i+1; k<nrDSPs; k++)
												{
													if ((tempc[k] != NULL) && (tempc[k] == addOps[j]))
														{
															cout << "tempc[" << k << "] deleted" << endl;
															tempc[k] = NULL;
															break;
														}
												}
					
											addOps[j]->getTopRightCorner(trx1, try1);
											addOps[j]->getBottomLeftCorner(blx1, bly1);
											fpadX = blx1-wInX-getExtraWidth()+1;
											fpadX = (fpadX<0)?0:fpadX;
											fpadY = bly1-wInY-getExtraHeight()+1;
											fpadY = (fpadY<0)?0:fpadY;
											bpadX = getExtraWidth()-trx1;
											bpadX = (bpadX<0)?0:bpadX;
											mPadX = (bpadX>mPadX)?bpadX:mPadX;
											bpadY = getExtraHeight()-try1;
											bpadY = (bpadY<0)?0:bpadY;
											mPadY = (bpadY>mPadY)?bpadY:mPadY;
											minPad = (minPad<(bpadY+bpadX))?minPad:(bpadY+bpadX);
											multW = addOps[j]->getMaxMultiplierWidth();
											multH = addOps[j]->getMaxMultiplierHeight();
					
											setCycle(0);
											xname.str("");
											xname << "x" << i << "_" << j+1;
											vhdl << tab << declare(xname.str(), multW, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadX,0) << " & " << "X" << range(blx1-fpadX, trx1+bpadX) << " & " << zg(bpadX,0) << ";" << endl;
											yname.str("");
											yname << "y" << i << "_" << j+1;
											vhdl << tab << declare(yname.str(), multH, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadY,0) << " & " << "Y" << range(bly1-fpadY, try1+bpadY) << " & " << zg(bpadY,0) << ";" << endl;
					
											nextCycle();
											mname.str("");
											mname << "mult_" << i << "_" << j+1;
											vhdl << tab << declare(mname.str(), multW+multH, true, Signal::registeredWithAsyncReset) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
										}
				
									nextCycle();
									vhdl << tab << declare(join("addDSP", nrOp), multW+multH+ext, true, Signal::registeredWithAsyncReset) << " <= ";
				
									int j=0;
				
									for (j=0; j<boundDSPs; j++)
										{
											mname.str("");
											mname << "mult_" << i << "_" << j;
											vhdl << "(" << zg(ext,0) << " & " << use(mname.str()) << ") + "; 
										}
				
									mname.str("");
									mname << "mult_" << i << "_" << j;
									vhdl << "(" << zg(ext,0) << " & " << use(mname.str()) << ");" << endl; 
				
								} 
							else // multiply the two terms and you're done
								{
									nextCycle();
									vhdl << tab << declare(join("addDSP", nrOp), multW+multH, true, Signal::registeredWithAsyncReset) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
								}
							vhdl << tab << declare(join("addOpDSP", nrOp), wInX+wInY) << " <= " << zg(wInX+minPad-blx1-ext,0) << " & " << zg(wInY-bly1,0) << " & " << use(join("addDSP", nrOp)) << range(multW+multH+ext-1, minPad) << " & " << zg(trx1+try1+minPad-mPadX-mPadY,0) << ";" << endl;
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
										cout<<"Partition number "<<count<<" is totally out of the real multiplication bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
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
										getTarget()->setUseHardMultipliers(false);
										IntMultiplier* mult =  new IntMultiplier(getTarget(), njj-nj+1, nii-ni+1);
										ostringstream cname;
										cname << mult->getName() << "_" << partitions;
										mult->changeName(cname.str());
										oplist.push_back(mult);
					
										vhdl << tab << declare(join("x_",partitions)) << " <= " << use("X") << range(wInX-nj-1+extW, wInX-njj-1+extW) << ";" << endl;
										inPortMap(mult, "X", join("x_",partitions));
										vhdl << tab << declare(join("y_",partitions)) << " <= " << use("Y") << range(nii, ni) << ";" << endl;
										inPortMap(mult, "Y", join("y_",partitions));
					
										outPortMap(mult, "R", join("addOpSlice", partitions));
					
										vhdl << instance(mult, join("Mult", partitions));

										syncCycleFromSignal(join("addOpSlice", partitions));
										cout<<"Partition number "<<count<<" with bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<" has now bounds ("<<nj<<" , "<<ni<<" , "<<njj<<" , "<<nii<<")"<<endl;
										partitions++;
									}
				
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
				
							}
					}
	
			}
	
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
	
		for(int ii=0;ii<m;ii++)
			delete[](mat[ii]);
	
		delete[] (mat);
	
		return partitions;
	}	


}

