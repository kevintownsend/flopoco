/*
 * A model of Stratix II FPGA 
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
#include "StratixII.hpp"
#include <iostream>
#include <sstream> 
#include "../utils.hpp"

double StratixII::adderDelay(int size) {
  return (fdCtoQ_ + lut2_ + muxStoO_ + shareOutToCarryOut_ + 
			((size-3) * fastcarryDelay_) + 
			((size/8) * (innerLABcarryDelay_- fastcarryDelay_)) + 
			((size/16) * (interLABcarryDelay_ - innerLABcarryDelay_)) + 
			carryInToSumOut_ + ffDelay_); 
};

double StratixII::carryPropagateDelay() {
  return  fastcarryDelay_; 
};

double StratixII::localWireDelay(){
  return lut2lutDelay_;
};

double StratixII::distantWireDelay(int n){
  return n*elemWireDelay_;
};

double StratixII::lutDelay(){
  return lutDelay_;
};

double StratixII::ffDelay(){
  return ffDelay_;
};

long StratixII::sizeOfMemoryBlock()
{
return sizeOfBlock;	
};

 DSP* StratixII::createDSP() 
{
	int multW, multH;
	getDSPWidths(multW, multH);
	
	// create a DSP block having a shifting amount of 0
	DSP* dsp_ = new DSP(0, multW, multH);
	
	return dsp_;
};

bool StratixII::suggestSubmultSize(int &x, int &y, int wInX, int wInY){
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
	 	 
bool StratixII::suggestSubaddSize(int &x, int wIn){

	int chunkSize = (int)floor( (1./frequency() - (fdCtoQ_ + lutDelay() + muxStoO_ + shareOutToCarryOut_ + carryInToSumOut_ + ffDelay_ + interLABcarryDelay_)) / carryPropagateDelay()); // 1 if no need for pipeline
	x = chunkSize;		
	if (x>0) return true;
	else {
		x=1;		
		return false;
	} 
}

bool  StratixII::suggestSlackSubaddSize(int &x, int wIn, double slack){

	int chunkSize = (int)floor( (1./frequency() - slack - (fdCtoQ_ + lutDelay() + muxStoO_ + shareOutToCarryOut_ + carryInToSumOut_ + ffDelay_ + interLABcarryDelay_)) / carryPropagateDelay()); // 1 if no need for pipeline
	
	//int chunkSize = 2 + (int)floor( (1./frequency() - slack - (fdCtoQ_ + slice2sliceDelay_ + lut2_ + muxcyStoO_ + xorcyCintoO_ + ffd_)) / muxcyCINtoO_ );
	
	x = chunkSize;		
	if (x>0) return true;
	else {
		x=1;		
		return false;
	} 


}
  
int StratixII::multiplierLUTCost(int wInX, int wInY){
	
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

void StratixII::getDSPWidths(int &x, int &y){
	// set the multiplier width acording to the desired frequency
	for (int i=0; i<3; i++)
		if (this->frequency() < 1/multiplierDelay_[i])
			x = y = multiplierWidth_[i];
}

int StratixII::getEquivalenceSliceDSP(){
	int lutCost = 0;
	int x, y;
	getDSPWidths(x,y);
	// add multiplier cost
	lutCost += multiplierLUTCost(x, y);
	// add accumulator cost
	//lutCost += accumulatorLUTCost(x, y);
	// add partial cost of final adder
	//lutCost += adderLUTCost(x,y);
	return lutCost;
}
	
int StratixII::getNumberOfDSPs() 
{
	int dsps = 96; // number of 9-bit elements
	int x, y;
	getDSPWidths(x, y);
	
	switch (x)
	{
		case 9: y = dsps;
			break;
		case 18: y = dsps/2;
			break;
		case 36: y = dsps/8;
			break;
	}
	return y;		
};

int StratixII::getIntNAdderCost(int wIn, int n)
{
	int chunkSize, lastChunkSize, nr, a, b, cost;
	
	suggestSubaddSize(chunkSize, wIn);
	lastChunkSize = wIn%chunkSize;
	nr = ceil((double) wIn/chunkSize);
	a = (nr-1)*wIn + (nr-1)*nr/2 + wIn*((n-1)*n/2-1);
	b = nr*lastChunkSize + (nr-2)*(nr-1)*chunkSize/2 + nr*(nr-1)/2 + (n-1)*chunkSize;
	cost = max(a,b)/2;
	return cost;
}

