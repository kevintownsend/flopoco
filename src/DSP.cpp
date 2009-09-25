/*
 * The class that models different digital signal processors. 
 * Classes for real chips instantiate this one giving it apropriate parameter values.
 *
 * Author : Sebastian Banescu, Radu Tudoran
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
#include "DSP.hpp"
using namespace std;

extern int verbose;

DSP::DSP(int Shift,int maxMultWidth,int maxMultHeight):
maxMultiplierWidth_(maxMultWidth) , maxMultiplierHeight_(maxMultHeight), fixedShift_(Shift)
{
			multiplierWidth_ = 0;
			multiplierHeight_ = 0;
			nrAdders_ = 0;
			multAccumulate_ = false;	
			
			DSP** ops = new DSP*[3];
	
			for (int j=0; j<3; j++)
			{
				ops[j] = NULL;
			}
			
			addOperands_ = ops;
			shiftIn_ = NULL;
			shiftOut_ = NULL;
}

int DSP::getMultiplierWidth(){
	return multiplierWidth_;
}

int DSP::getMaxMultiplierWidth(){
	return maxMultiplierWidth_;
}

int DSP::getMultiplierHeight(){
	return multiplierHeight_;
}

int DSP::getMaxMultiplierHeight(){
	return maxMultiplierHeight_;
}


void DSP::setMultiplierWidth(int w){
	multiplierWidth_ = w;
}

void DSP::setMultiplierHeight(int h){
	multiplierHeight_ = h;
}

int DSP::getShiftAmount(){
	return fixedShift_;
}

void DSP::setShiftAmount(int s){
	fixedShift_ = s;
}

int DSP::getNumberOfAdders(){
	return nrAdders_;
}

void DSP::setNumberOfAdders(int o){
	nrAdders_ = o;
}

bool DSP::isMultiplyAccumulate(){
	return multAccumulate_;
}

void DSP::setMultiplyAccumulate(bool m){
	multAccumulate_ = m;
}

void DSP::getTopRightCorner(int &x, int &y){
	x = positionX_[0];
	y = positionY_[0];
}

void DSP::setTopRightCorner(int x, int y){
	positionX_[0] = x;
	positionY_[0] = y;
}

void DSP::getBottomLeftCorner(int &x, int &y){
	x = positionX_[1];
	y = positionY_[1];
}

void DSP::setBottomLeftCorner(int x, int y){
	positionX_[1] = x;
	positionY_[1] = y;
}

DSP* DSP::getShiftIn(){
	return shiftIn_;
}

void DSP::setShiftIn(DSP* d){
	shiftIn_ = d;
}

DSP* DSP::getShiftOut(){
	return shiftOut_;
}

void DSP::setShiftOut(DSP* d){
	shiftOut_ = d;
}

DSP** DSP::getAdditionOperands(){
	return addOperands_;
}

void DSP::setAdditionOperands(DSP** o){
	addOperands_ = o;
}

void DSP::rotate(){
	int tx,ty;
	getTopRightCorner(tx,ty);
	//cout<<"Old dimensions were "<<maxMultiplierHeight_<<" "<<maxMultiplierWidth_<<endl;
	int tmp = maxMultiplierHeight_;
	maxMultiplierHeight_ = maxMultiplierWidth_;
	maxMultiplierWidth_ = tmp;
	//cout<<"New dimensions are "<<maxMultiplierHeight_<<" "<<maxMultiplierWidth_<<endl;
	
	//getBottomLeftCorner(tx,ty);
	//cout<<"Old were "<<tx<<" "<<ty<<endl;
	setBottomLeftCorner(tx + maxMultiplierWidth_-1 , ty + maxMultiplierHeight_-1);
	//getBottomLeftCorner(tx,ty);
	//cout<<"New are "<<tx<<" "<<ty<<endl;
	
}
