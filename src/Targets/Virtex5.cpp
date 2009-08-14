/*
 * A model of FPGA that works well enough for Virtex-5 chips  (LX, speedgrade -3, e.g. xc5vlx30-3) 
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
#include "Virtex5.hpp"
#include <iostream>
#include <sstream>
#include "../utils.hpp"

double Virtex5::adderDelay(int size) {
  return lut2_ + muxcyStoO_ + double(size-1)*muxcyCINtoO_ + xorcyCintoO_ ; 
};


double Virtex5::ffDelay() {
  return fdCtoQ_ + ffd_; 
};

double Virtex5::carryPropagateDelay() {
  return  fastcarryDelay_; 
};

double Virtex5::localWireDelay(){
  return  elemWireDelay_ ;
};

double Virtex5::distantWireDelay(int n){
  return n*elemWireDelay_;
};

double Virtex5::lutDelay(){
  return lutDelay_;
};

long Virtex5::sizeOfMemoryBlock()
{
return sizeOfBlock;	
};


 DSP* Virtex5::createDSP() 
{
	int x, y;
	getDSPWidths(x, y);
	
	/* create DSP block with constant shift of 17
	 * and a maxium unsigned multiplier width (17x17) for this target
	 */
	DSP* dsp_ = new DSP(dspFixedShift_, x, y);
	
	return dsp_;
};

bool Virtex5::suggestSubmultSize(int &x, int &y, int wInX, int wInY){
	
	getDSPWidths(x, y);
	
	if (getUseHardMultipliers())
	{
		if (wInX <= x)
			x = wInX;
			
		if (wInY <= y)
			y = wInY;
	}
	return true;
};	 
	 
bool Virtex5::suggestSubaddSize(int &x, int wIn){

//	int chunkSize = (int)floor( (1./frequency() - lutDelay()) / carryPropagateDelay()); // 1 if no need for pipeline


	int chunkSize = 2 + (int)floor( (1./frequency() - (fdCtoQ_ + slice2sliceDelay_ + lut2_ + muxcyStoO_ + xorcyCintoO_ + ffd_)) / muxcyCINtoO_ );
	x = chunkSize;		
	if (x > 0) 
		return true;
	else {
		x = 2;		
		return false;
	} 
};

bool Virtex5::suggestSlackSubaddSize(int &x, int wIn, double slack){

//	int chunkSize = (int)floor( (1./frequency() - lutDelay()) / carryPropagateDelay()); // 1 if no need for pipeline


	int chunkSize = 2 + (int)floor( (1./frequency() - slack - (fdCtoQ_ + slice2sliceDelay_ + lut2_ + muxcyStoO_ + xorcyCintoO_ + ffd_)) / muxcyCINtoO_ );
	x = chunkSize;		
	if (x > 0) 
		return true;
	else {
		x = 2;		
		return false;
	} 
};

int Virtex5::getIntMultiplierCost(int wInX, int wInY){
	
	int lutCost = 0;
	int chunkSize_ = this->lutInputs()/2;
	
	int chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
	int chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
	
	for (int i=0; i<chunksX; i++)
		for (int j=0; j<chunksY; j++)
			lutCost += 3; // one LUT for each: CY<4>, lut<1>, product of 1 pair of underlying LSBs
	// NOTE: each 6 input LUT function has 2 outputs
			
	return lutCost + this->getIntNAdderCost(chunksX*chunkSize_, chunksY*chunkSize_);
};
  
void Virtex5::getDSPWidths(int &x, int &y){
	// unsigned multiplication
	x = multXInputs_-1;
	y = multYInputs_-1;
}

int Virtex5::getEquivalenceSliceDSP(){
	int lutCost = 0;
	int x, y;
	getDSPWidths(x,y);
	// add multiplier cost
	lutCost += getIntMultiplierCost(x, y);
	// add shifter and accumulator cost
	//lutCost += accumulatorLUTCost(x, y);
	return lutCost;
}

int Virtex5::getNumberOfDSPs() 
{
	return nrDSPs_; 		
};

int Virtex5::getIntNAdderCost(int wIn, int n)
{
	int chunkSize, lastChunkSize, nr, a, b, cost;
	
	suggestSubaddSize(chunkSize, wIn);
	lastChunkSize = wIn%chunkSize;
	nr = ceil((double) wIn/chunkSize);
	a = (nr-2)*lastChunkSize/2 + nr*wIn/2 + 2*(nr-1 + wIn) + chunkSize + (2*wIn+(nr-1))*(n-3);
	b = nr*lastChunkSize/2 + nr*wIn/2 + (nr-2)*chunkSize + wIn + 2*wIn*(n-3);
	cost = a+b*0.25;
	return cost;
}
