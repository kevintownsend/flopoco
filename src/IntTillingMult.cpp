/*
 * Floating Point Adder for FloPoCo
 *
 * Author :  Radu Tudoran, Bogdan Pasca
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

#include <gmp.h>


#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "IntTillingMult.hpp"

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
extern vector<Operator*> oplist;

#define DEBUGVHDL 0

IntTilingMult:: IntTilingMult(Target* target, int wInX, int wInY,float ratio) :
	Operator(target), target(target),wInX(wInX), wInY(wInY), wOut(wInX + wInY),ratio(ratio){
 
	ostringstream name;

	name <<"IntMultiplier_"<<wInX<<"_"<<wInY;
	setName(name.str());
		
	setCopyrightString("Bogdan Pasca, Sebastian Banescu , Radu Tudoran 2009");
	
	addInput ("X", wInX);
	addInput ("Y", wInY);
	addOutput("R", wOut); /* wOut = wInX + wInY */
		
	nrDSPs = estimateDSPs();
	cout<<"Estimated DSPs:= "<<nrDSPs <<endl;
	int x,y;
	target->getDSPWidths(x,y);
	cout<<"Width of DSP is := "<<x<<" Height of DSP is:="<<y<<endl;
	cout<<"Extra width:= "<<getExtraWidth()<<" \nExtra height:="<<getExtraHeight()<<endl;
		

	/*we will try the algorithm with 2 values of nrDSPs	
	One will be the estimated value(nrDSPs) and the second one will be nrDSPs-1	
	*/
		
	
		
		//~ cout<<"Number of slices for multiplication is "<<partitionOfGridSlices(globalConfig,t)<<endl;
		//~ sortDSPs(globalConfig);
		
		//~ cout<<"Number of slices for multiplication is "<<partitionOfGridSlices(globalConfig,t)<<endl;
		
		
		
		//cout<<"Cost of configuration is "<<computeCost(globalConfig)<<endl;
		
		//cout<<"Initial configuration"<<endl<<endl;
		//display(globalConfig);
		
		
		
		
		//~ cout<<"Partitions of grid before binding "<<partitionOfGridSlices(globalConfig,t)<<endl;
		//~ cout<<"Operands from DSPs, from a maximum of "<<nrDSPs<<" we use "<<bindDSPs(globalConfig)<<endl;
		//~ cout<<"Partitions of grid after binding "<<partitionOfGridSlices(globalConfig,t)<<endl;
		
		
		
		
		
		//~ for(int i=0;i<nrDSPs;i++)
		//~ {
			//~ if(globalConfig[i]!=NULL)
			//~ {
				//~ if(globalConfig[i]->getShiftIn()!=NULL)
				//~ {
					//~ for(int j=0;j<nrDSPs;j++)
						//~ if(globalConfig[j]==globalConfig[i]->getShiftIn())
							//~ cout<<"Exista legatura In la dspul "<<i+1<<" de la DSPul"<<j+1<<endl;
				//~ }
				//~ if(globalConfig[i]->getShiftOut()!=NULL)
				//~ {
					//~ for(int j=0;j<nrDSPs;j++)
						//~ if(globalConfig[j]==globalConfig[i]->getShiftOut())
							//~ cout<<"Exista legatura Out la dspul "<<i+1<<" catre DSPul"<<j+1<<endl;
				//~ }
			//~ }
		//~ }
		
		
		
	runAlgorithm();
	
	//initTiling(bestConfig,nrDSPs);
	
	//~ bestCost=350.455;
	/*
	initTiling(globalConfig,nrDSPs);
	
	globalConfig[0]->setTopRightCorner(19,0);
	globalConfig[0]->setBottomLeftCorner(35,16);
	globalConfig[1]->setTopRightCorner(10,19);
	globalConfig[1]->setBottomLeftCorner(26,35);
	replace(globalConfig, 1);
	
	display(globalConfig);
	*/
	
	
	//~ compareCost();
	
	//~ cout<<"etapa 2"<<endl<<endl;

	
	//~ initTiling(globalConfig,nrDSPs);
	
	//~ globalConfig[0]->setTopRightCorner(2,5);
	//~ globalConfig[0]->setBottomLeftCorner(18,21);
	//~ globalConfig[1]->setTopRightCorner(19,5);
	//~ globalConfig[1]->setBottomLeftCorner(35,21);
	
	
		
	
	//~ cout<<endl<<endl;
	//~ display(globalConfig);
	//~ cout<<endl<<endl;	
	//~ //display(bestConfig);
	
	//~ compareCost();
	
	
	cout<<"Estimated DSPs:= "<<nrDSPs <<endl;
	target->getDSPWidths(x,y);
	cout<<"Width of DSP is := "<<x<<" Height of DSP is:="<<y<<endl;
	cout<<"Extra width:= "<<getExtraWidth()<<" \nExtra height:="<<getExtraHeight()<<endl;
	
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
		
		
	
		
	//~ /*we will try the algorithm with 2 values of nrDSPs	
	//~ One will be the estimated value(nrDSPs) and the second one will be nrDSPs-1	
	//~ */
	//~ rot = new bool[nrDSPs];
	//~ for(int i =0;i<nrDSPs;i++)
		//~ rot[i]=false;
	
	
	
	 //~ initTiling(globalConfig,nrDSPs);	
			
	//~ //this will initialize the bestConfig with the first configuration
	//~ bestCost = FLT_MAX ;
	//~ cout<<"Max score is"<<bestCost<<endl;
	//~ //bestConfig = (DSP**)malloc(nrDSPs * sizeof(DSP*));
	//~ bestConfig = new DSP*[nrDSPs];
	//~ for(int i=0;i<nrDSPs;i++)
			//~ bestConfig[i]= new DSP();
	//~ compareCost();
	//~ cout<<"New best score is"<<bestCost<<endl;
	
	//~ display(bestConfig);
	
	
	//~ //the best configuration should be consider initially the first one. So the bestConfig parameter will be initialized with global config and hence the bestCost will be initialized with the first cost
	
	
	//~ tilingAlgorithm(nrDSPs-1,nrDSPs-1,false);
	
	//~ display(bestConfig);
	//~ cout<<"Best cost is "<<bestCost<<endl;
	
	
	
	
	
	


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
	
	display(bestConfig);
	
	
	//the best configuration should be consider initially the first one. So the bestConfig parameter will be initialized with global config and hence the bestCost will be initialized with the first cost
	
	
	tilingAlgorithm(nrDSPs-1,nrDSPs-1,false);
	
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
	
}

void IntTilingMult::tilingAlgorithm(int i, int n,bool repl)
{

if(i==n)
{	
	if(repl==true)
	{
		//~ cout<<" Pas 4_1 "<<i<<endl;
		replace(globalConfig,i);
		compareCost();
		rot[i]=false;
		tilingAlgorithm(i,n,false);	
	}
	else
	{
	
		if(move(globalConfig,i))
		{
			//~ cout<<" Pas 1_1 "<<i<<endl;
			compareCost();
			tilingAlgorithm(i,n,repl);		//repl should be false
		}
		else
		{
			if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() ))
			{
				//~ cout<<" Pas 2_1 "<<i<<endl;
				globalConfig[i]->rotate();
				rot[i]=true;
				replace(globalConfig,i);
				compareCost();
				tilingAlgorithm(i,n,repl);		//repl should be false
			}
			else
			{
				//~ cout<<" Pas 3_1 "<<i<<endl;
				if(i-1>=0)
				tilingAlgorithm(i-1,n,repl);		//repl should be false
			}
		}
	}
}
else
{
	if(repl==true)
	{
		//~ cout<<" Pas 4_2 "<<i<<endl;
		replace(globalConfig,i);
		rot[i]=false;
		tilingAlgorithm(i+1,n,repl);
		
	}
	else
	{	
		if(move(globalConfig,i))
		{
			//~ cout<<" Pas 1_2 "<<i<<endl;
			if(i==0)
				cout<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Primul a mai facut un pas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl<<endl;
			tilingAlgorithm(i+1,n,true);
		}
		else
		{
			if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() ))
			{
				//~ cout<<" Pas 2_2 "<<i<<endl;
				globalConfig[i]->rotate();
				replace(globalConfig,i);
				rot[i]=true;
				tilingAlgorithm(i+1,n,true);
			}
			else
			{
				//~ cout<<" Pas 3_2 "<<i<<endl;
				if(i-1>=0)
				tilingAlgorithm(i-1,n,repl);		//repl should be false
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
	//cout<<"Total size is "<<totalSize<<endl;
	int c1X,c2X,c1Y,c2Y;
	int nj,ni,njj,nii;
	int ii,jj,i,j;
	int n,m;
	n=wInX + 2* getExtraWidth();
	m=wInY + 2* getExtraHeight();
	
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
			
			
						
						if( j>= n-getExtraWidth() || jj< getExtraWidth() || i >= m - getExtraHeight() || ii < getExtraHeight())
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
			
			
						
						if( j>= n-getExtraWidth() || jj< getExtraWidth() || i >= m - getExtraHeight() || ii < getExtraHeight())
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
											if(ver=0)
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
											if(ver=3)
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
						
						if( j>= n-getExtraWidth() || jj< getExtraWidth() || i >= m - getExtraHeight() || ii < getExtraHeight())
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
	int **mat;
	int n,m;
	int count=1;
	n=wInX + 2* getExtraWidth();
	m=wInY + 2* getExtraHeight();
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
	//partitions = count;
	partitions = 0;
		
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
											if(ver=0)
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
											if(ver=3)
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
						
						int nj,ni,njj,nii;
						
						
						if( j>=n-getExtraWidth() || jj< getExtraWidth() || i >= m - getExtraHeight() || ii < getExtraHeight())
						{
							
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
							
							partitions++;
							
							costSlice +=target->getIntMultiplierCost(njj-nj+1,nii-ni+1);
							
							
														
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
	
	return costSlice;
}	

//~ void IntTilingMult::resetConnections(DSP** &config)
//~ {
//~ for(int i=0;i<nrDSPs;i++)
	//~ {
		//~ if(config[i]!=NULL)
		//~ {
			
		//~ }
		
	//~ }

	
//~ }

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
	if (strncmp(typeid(*target).name(), "7Virtex", 7) == 0) // then the target is Virtex
	{
		//cout<<"Virtex!"<<endl;
		return bindDSPs4Virtex(config);
	}
	else // the target is Stratix
	{
		//cout<<"Stratix"<<endl;
		return bindDSPs4Stratix(config);
	}
		
}


void IntTilingMult::compareCost()
{
	DSP** tempc;
		
	tempc= new DSP*[nrDSPs];
	for(int ii=0;ii<nrDSPs;ii++)
		tempc[ii]= new DSP();
		
	//display(globalConfig);
	
	//memcpy(tempc,globalConfig,sizeof(DSP*) *nrDSPs );
	for(int ii=0;ii<nrDSPs;ii++)
		memcpy(tempc[ii],globalConfig[ii],sizeof(DSP) );
	
	//display(tempc);
	
	float temp = computeCost(tempc);
	
	//cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
	
	if(temp < bestCost)
	{
		cout<<"Costul e mai bun la cel curent!Schimba"<<endl;
		cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
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
		display(bestConfig);
	}
	else
		if(temp == bestCost )
		{
			cout<<"Cost egal!!!"<<endl;
			cout<<"Rezult compare is"<<compareOccupation(tempc)<<endl;
			cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
			//display(bestConfig);
			if(compareOccupation(tempc)==true)
			{
			
				bestCost=temp;
			
				for(int ii=0;ii<nrDSPs;ii++)
					memcpy(bestConfig[ii],tempc[ii],sizeof(DSP) );
				display(bestConfig);
			}
		}
	
	
	for(int ii=0;ii<nrDSPs;ii++)
		free(tempc[ii]);
	
	delete[] (tempc);
	
}


// de revazut costul sa fie calculat nu si pentru extensii!!!!!!!!!
float IntTilingMult::computeCost(DSP** &config)
{
	
	float const scale=100.0;
	float acc=0.0;
	float costDSP,costLUT;
	costDSP = ( (1.0+scale) - scale * ratio );
	costLUT = ( (1.0+scale) - scale * (1-ratio) ) /  ((float) target->getEquivalenceSliceDSP() );
	
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
		acc =((float)nrOfUsedDSPs)*costDSP + costLUT * partitionOfGridSlices(config,partitions);
		
	//~ cout<<"Number of partitions for LUTs is "<<partitions<<endl;
	nrOfUsedDSPs = bindDSPs(config);
	//~ cout<<"Number of operands coming from DSPs is "<<nrOfUsedDSPs<<endl;
		
	
	 float LUTs4NAdder=((float)target->getIntNAdderCost(wInX + wInY,nrOfUsedDSPs+partitions) );
	
	//~ cout<<"LUTs used for last "<<nrOfUsedDSPs+partitions<<" adder are"<<LUTs4NAdder<<endl;
		
	acc +=  LUTs4NAdder* costLUT;	
		
			//~ Substracting the cost of different additions that can be done in DSPs(Virtex) or the cost of a DSP if they can be combined(Altera)
		
		
	return acc;
			
}

int IntTilingMult::estimateDSPs()
{
	float t1,t2;
	int Xd,Yd;
	target->getDSPWidths(Xd,Yd);
	bool ver=true;
	int maxDSP= target->getNumberOfDSPs();
	
	//~ cout<<"Existing DSP on the board "<< target->getNumberOfDSPs()<<endl;
	
	t1= ((float) wInX) / ((float) Xd);
		t2= ((float) wInY) / ((float) Yd);
		
		if(t1 - ((int)t1) >0)
		{
			t1++;
		}
		if(t2 - ((int)t2) >0)
		{
			t2++;
		}
		
		if(maxDSP> ((int)t1) * ((int)t2) )
		{
			ver=false;
			maxDSP = ((int)t1) * ((int)t2);
			
		}
	
	if(1==ratio) 
	{
		if(ver==true)
			cout<<"Warning!!! The number of existing DSPs on this board is not enough to cover the whole multiplication!"<<endl;
		else
			cout<<"Warning!!! The minimum number of DSP that is neccessary to cover the whole addition will be used!"<<endl;
		return maxDSP;
	}
	else
	{	
		//cout<<target->multiplierLUTCost(wInX,wInY) <<"    "<<target->getEquivalenceSliceDSP()<<endl;
		float temp = (target->getIntMultiplierCost(wInX,wInY) * ratio)  /   ((1- ratio)  *  target->getEquivalenceSliceDSP() ) ;
		//cout<<temp<<endl;
		if(temp - ((int)temp)  > 0 ) // if the estimated number of dsps is a number with nonzero real part than we will return the integer number grather then this computed value. eg. 2.3 -> 3
		{
		//cout<<"Intra"<<endl;	
		temp++;	
		}
		
		if(temp>maxDSP)
		{
			if(ver==false)
				cout<<"Warning!!! The number of estimated DSPs with respect with this ratio of preference is grather then the needed number of DSPs to perform this multiplication!"<<endl;
			else
				cout<<"Warning!!! The number of estimated DSPs with respect with this ratio of preference is grather then the total number of DSPs that exist on this board!"<<endl;
			
			temp=maxDSP;
		}
		
		return ((int)  temp) ;
	}
	
	return 0;
}


int  IntTilingMult::getExtraHeight()
{
int x,y;	
target->getDSPWidths(x,  y);
float temp = ratio * 0.75 * ((float) y);
return ((int)temp);
}

	
int  IntTilingMult::getExtraWidth()
{
int x,y;	
target->getDSPWidths(x,y);
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


bool IntTilingMult::checkOverlap(DSP** config, int index)
{	
	int xtr1, ytr1, xbl1, ybl1, xtr2, ytr2, xbl2, ybl2;
	//int nrdsp = sizeof(config)/sizeof(DSP*);
	config[index]->getTopRightCorner(xtr1, ytr1);
	config[index]->getBottomLeftCorner(xbl1, ybl1);
	
	if(verbose)
		cout << tab << tab << "checkOverlap: ref is block #" << index << ". Top-right is at (" << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
	
	for (int i=0; i<nrDSPs; i++)
		if ((config[i] != NULL) && (i != index))
		{
			config[i]->getTopRightCorner(xtr2, ytr2);
			config[i]->getBottomLeftCorner(xbl2, ybl2);
			
			if ((i < index) && 
				(((xtr2 >= xtr1) && (ytr2 > ybl1)) || 	// config[index] is above and to the right of config[i]
				 (xbl1 < xtr2) 							// config[index] is to the right of config[i]
				)
			   )
			   return true;
			
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
				((xbl1 >= xbl2) && (ybl1 <= ybl2) && (ytr1 >= ytr2) && (xtr1 <= xtr2))    // config[i] overlaps the center part of config[index]
				)
				return true;
		}	
	if(verbose)
		cout << tab << tab << "checkOverlap: return false" << endl;	
	return false;
}

bool IntTilingMult::move(DSP** config, int index)
{
	int xtr1, ytr1, xbl1, ybl1;
	int w, h;
	target->getDSPWidths(w,h);
	
	
	config[index]->getTopRightCorner(xtr1, ytr1);
	config[index]->getBottomLeftCorner(xbl1, ybl1);
	
	if(verbose)
		cout << tab << "replace : DSP #" << index << " width is " << w << ", height is " << h << endl; 
	
	if ((xtr1 < 0) || (ytr1 < 0) || (xbl1 > wInX+2*getExtraWidth()-1) || (ybl1 > wInY+2*getExtraHeight()-1))
	{// then the DSP block is placed outside the bounds 
		do{
			// move down one unit
			ytr1++;
			ybl1++;
			
			if (ytr1 > wInY+getExtraHeight()-1) // the DSP block has reached the bottom limit of the tiling grid
			{
				// move to top of grid and one unit to the left 
				xtr1++;
				xbl1++;
				ytr1 = getExtraHeight() - h+1;
				ybl1 = ytr1 + h-1;
				
				if (xtr1 > wInX+getExtraWidth()-1) // the DSP block has reached the left limit of the tiling grid
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
			ytr1++;
			ybl1++;
			
			if (ybl1 > wInY+2*getExtraHeight()-1) // the DSP block has reached the bottom limit of the tiling grid
			{
				// move to top of grid and one unit to the left 
				xtr1++;
				xbl1++;
				ytr1 = 0;
				ybl1 = h-1;
				
				if (xbl1 > wInX+2*getExtraWidth()-1) // the DSP block has reached the left limit of the tiling grid
					return false;
			}			
			
			config[index]->setTopRightCorner(xtr1, ytr1);
			config[index]->setBottomLeftCorner(xbl1, ybl1);
		}while (checkOverlap(config, index));
	}
	// set the current position of the DSP block within the tiling grid
	config[index]->setTopRightCorner(xtr1, ytr1);
	config[index]->setBottomLeftCorner(xbl1, ybl1);
	return true;
}

void IntTilingMult::replace(DSP** config, int index)
{
	int xtr1, ytr1, xbl1, ybl1;
	int w, h;
	target->getDSPWidths(w,h);
	
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
 		if (ybl1 > wInY+2*getExtraHeight()) // the DSP block has reached the bottom limit of the tiling grid
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
	
	if (xbl1 > wInX+2*getExtraWidth()) // it was not possible to place the DSP inside the extended grid
	{ // then try to place it anywhere around the tiling grid such that it covers at leat one 1x1 square
		xtr1 = getExtraWidth() - w+1;
		ytr1 = getExtraHeight() - h+1;
		xbl1 = xtr1 + w-1;
		ybl1 = ytr1 + h-1;
		
		config[index]->setTopRightCorner(xtr1, ytr1);
		config[index]->setBottomLeftCorner(xbl1, ybl1);
		
		if(verbose)
			cout << tab << "replace : DSP width is " << w << ", height is " << h << endl; 
		while (checkOverlap(config, index))
		{
			// move down one unit
			ytr1++;
			ybl1++;
			
			if(verbose)
				cout << tab << "replace : moved DSP one unit down." << endl;
			if (ytr1 > wInY+getExtraHeight()) // the DSP block has reached the bottom limit of the tiling grid
			{
				// move to top of grid and one unit to the left 
				xtr1++;
				//if (xtr1 > wInX+getExtraWidth()) // the DSP block has reached the left limit of the tiling grid
				xbl1++;
				ytr1 = getExtraHeight() - h+1;
				ybl1 = ytr1 + h-1;
				
				if(verbose)
					cout << tab << "replace : moved DSP up and one unit left." << endl;
			}			
			
			config[index]->setTopRightCorner(xtr1, ytr1);
			config[index]->setBottomLeftCorner(xbl1, ybl1);
			if(verbose)
				cout << tab << "replace : Top-right is at ( " << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
		}
	
	}
	// set the current position of the DSP block within the tiling grid
	config[index]->setTopRightCorner(xtr1, ytr1);
	config[index]->setBottomLeftCorner(xbl1, ybl1);
}

void IntTilingMult::initTiling(DSP** &config, int dspCount)
{
	config = new DSP*[nrDSPs];
	
	for (int i=0; i<nrDSPs; i++)
	{
		config[i] = NULL;
	}
	
	for (int i=0; i<dspCount; i++)
	{
		if(verbose)
			cout << "initTiling : iteration #" << i << endl; 
		config[i] = target->createDSP();
		replace(config, i);
	}
}
	
int IntTilingMult::bindDSPs4Stratix(DSP** config)
{
	int nrDSPs_ = estimateDSPs();
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
			target->getDSPWidths(wx,wy);
				
			if (verbose)
				cout << "bindDSP4Stratix: DSP #" << i << " has less than 3 operands. Top-right is at ( " << xtri << ", " << ytri << ") width is " << wx << endl;
			
			for (int j=0; j<DSPcount; j++)
			if (i != j)
			{
				nrOp = config[i]->getNumberOfAdders();
				if (nrOp == 3)
					break;
			
				DSP** operands = config[i]->getAdditionOperands();
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
					DSP** operandsj = config[j]->getAdditionOperands();
					
					for (int k=0; k<nrOp; k++)
					{
						operandsj[nrOpj+k] = operands[k];
						operands[k]->setNumberOfAdders(nrOp+nrOpj+1);
					}
					operandsj[nrOp+nrOpj] = config[i];
					config[j]->setAdditionOperands(operandsj);
					config[j]->setNumberOfAdders(nrOp+nrOpj+1);
					
					for (int k=0; k<nrOpj; k++)
					{
						operands[nrOp+k] = operandsj[k];
						operandsj[k]->setNumberOfAdders(nrOp+nrOpj+1);
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
