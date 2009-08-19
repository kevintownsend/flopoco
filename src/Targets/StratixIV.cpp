/*
 * A model of Stratix IV FPGA optimized for (EP4S40G2F40C2ES1 speed grade 2)
 *
 * Author : Florent de Dinechin, Sebastian Banescu
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
#include "StratixIV.hpp"
#include <iostream>
#include <sstream> 
#include "../utils.hpp"

double StratixIV::adderDelay(int size) {
  return (distantWireDelay(10) + fdCtoQ_ + lutDelay() + 
			((size-3) * fastcarryDelay_) + 
			((size/almsPerLab_) * (innerLABcarryDelay_- fastcarryDelay_)) + 
			((size/(almsPerLab_*2)) * (interLABcarryDelay_ - innerLABcarryDelay_)) + 
			carryInToSumOut_ + ffDelay_); 
};

double StratixIV::carryPropagateDelay() {
  return  fastcarryDelay_; 
};

double StratixIV::localWireDelay(){
  return lut2lutDelay_;
};

double StratixIV::distantWireDelay(int n){
  return n*elemWireDelay_;
};

double StratixIV::lutDelay(){
  return lutDelay_;
};

double StratixIV::ffDelay(){
  return ffDelay_;
};

long StratixIV::sizeOfMemoryBlock()
{
return sizeOfBlock;	
};

 DSP* StratixIV::createDSP() 
{
	int multW, multH;
	getDSPWidths(multW, multH);
	
	// create a DSP block having a shifting amount of 0
	DSP* dsp_ = new DSP(0, multW, multH);
	
	return dsp_;
};

bool StratixIV::suggestSubmultSize(int &x, int &y, int wInX, int wInY){
// TODO This is the VirtexIV function. Stratix II is more interesting
// (DSP blocks are 36x36 and my be split as 9x9 or 18x18)
	if (getUseHardMultipliers()){
		x = y = 0;//max(wInX, wInY);
		int padX[3], padY[3]; // nr of zero padding for a specific width multiplier
		double maxF; // will hold the maximum possible freqeuncy for each multiplier width
		for (int i=0; i<3; i++){ // for each multiplier width available
			maxF = 1/(inputRegDelay_[i] + multiplierDelay_[i]); // maximum possible freqeuncy 
			padX[i] = multiplierWidth_[i] - (wInX % multiplierWidth_[i]);
			padY[i] = multiplierWidth_[i] - (wInY % multiplierWidth_[i]);
			
			if ((wInY < multiplierWidth_[i]) && (y != 0)){
				if ((i > 0) && (padY[i] > 2*padY[i-1]))
					y = multiplierWidth_[i-1];
				else
					y = wInY;
			}
				
			if ((wInX < multiplierWidth_[i]) && (x != 0)){
				if ((i > 0) && (padX[i] > 2*padY[i-1]))
					x = multiplierWidth_[i-1];
				else
					x = wInX;
			}
			
			if ((x != 0) && (y != 0))
			{
				if (frequency() < maxF)
					return true;
				else
					return false;
			}
		}
		
		// we have the maximum frequency for 36x36 in maxF
		if (maxF > frequency()){ // for lower freqency we prefer 36x36
			if (x == 0)
				x = 36;
			if (y == 0)
				y = 36;
			return true;
		}else{	// to obtain the highest freqency with lower logic utilization we need 18x18
			if (x == 0)
				x = 9;
			if (y == 0)
				y = 36;
			return true;
		}
	}else{
		// TODO functional approximation of multiplier size based on frequency
		x = y = 18;
		return true;
	}
		
	// control should never get here
	return false;
}	 
	 	 
bool StratixIV::suggestSubaddSize(int &x, int wIn){

	//int chunkSize = (int)floor( (1./frequency() - (fdCtoQ_ + lutDelay() + carryInToSumOut_ + ffDelay_ + interLABcarryDelay_)) / carryPropagateDelay()); // 1 if no need for pipeline
	suggestSlackSubaddSize(x, wIn, 0);
	
	//x = chunkSize;		
	if (x>0) return true;
	else {
		x=1;		
		return false;
	} 
}

bool  StratixIV::suggestSlackSubaddSize(int &x, int wIn, double slack){

	float time = 1./frequency() - slack - (distantWireDelay(10) + fdCtoQ_ + lutDelay() + carryInToSumOut_ + ffDelay_);
	cout << "suggestSubaddSize : time = " << time << endl;
	int carryFlag = 0;
	int chunkSize = 0;
	
	while (time > 0)
	{
		chunkSize++;
		cout << "suggestSubaddSize : chunkSize = " << chunkSize << " time = " << time << endl;
		
		if (carryFlag == 0) 
		{
			time -= fastcarryDelay_;
			cout << "suggestSubaddSize : subtracted fastcarryDelay"<< endl;
			if (chunkSize % (almsPerLab_*2) == 0)
				carryFlag = 2;
			else if (chunkSize % almsPerLab_ == 0)
				carryFlag = 1;
			
		}
		else if (carryFlag == 1)
		{
			time -= innerLABcarryDelay_;
			cout << "suggestSubaddSize : subtracted innerLABCarryDelay"<< endl;
			carryFlag = 0;
		}	
		else if (carryFlag == 2)
		{
			time -= interLABcarryDelay_;
			cout << "suggestSubaddSize : subtracted interLABCarryDelay"<< endl;
			carryFlag = 0;
		}
	}
	chunkSize--; // decremented because of the loop condition (time > 0). When exiting the loop the time is negative
	
	x = chunkSize;		
	if (x>0) return true;
	else {
		x=1;		
		return false;
	} 


}
  
int StratixIV::getIntMultiplierCost(int wInX, int wInY){
	
	int lutCost = 0;
	//int chunkSize_ = this->lutInputs()/2; SYNTHESIS NOT EFFICIENT WITH 6-LUTs
	int chunkSize_ = 2; 
	
	int chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
	int chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));

	for (int i=0; i<chunksX; i++)
		for (int j=0; j<chunksY; j++)
			lutCost += 4; // one lut for each bit of a partial product
						
	return lutCost + this->getIntNAdderCost(chunksX*chunkSize_, chunksY*chunkSize_);

}

void StratixIV::getDSPWidths(int &x, int &y){
	// set the multiplier width acording to the desired frequency
	for (int i=0; i<3; i++)
		if (this->frequency() < 1/multiplierDelay_[i])
			x = y = multiplierWidth_[i];
}

int StratixIV::getEquivalenceSliceDSP(){
	int lutCost = 0;
	int x, y;
	getDSPWidths(x,y);
	// add multiplier cost
	lutCost += getIntMultiplierCost(x, y);
	// add accumulator cost
	//lutCost += accumulatorLUTCost(x, y);
	// add partial cost of final adder
	//lutCost += adderLUTCost(x,y);
	return lutCost;
}
	
int StratixIV::getNumberOfDSPs() 
{
	int dsps = 1288; // number of 9-bit elements
	int x, y;
	getDSPWidths(x, y);
	
	switch (x)
	{
		case 9: y = dsps;
			break;
		case 12: y = dsps*3/4;
			break;
		case 18: y = dsps/2;
			break;
		case 36: y = dsps/4;
			break;
	}
	return y;		
};

int StratixIV::getIntNAdderCost(int wIn, int n)
{
	int chunkSize, lastChunkSize, nr, a, b, cost;
	
	suggestSubaddSize(chunkSize, wIn);
	lastChunkSize = wIn%chunkSize;
	nr = ceil((double) wIn/chunkSize);
	a = (nr-1)*wIn + (nr-1)*nr/2 + wIn*((n-1)*n/2-1)+nr*(n-2);
	b = nr*lastChunkSize + nr*(nr-1)*(chunkSize+1)/2 + nr*(n-1) + (n-2)*wIn;
	cost = (a+b)/2;
	return cost;
}

