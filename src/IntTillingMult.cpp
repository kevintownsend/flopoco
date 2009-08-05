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

#include <stdio.h>
#include <mpfr.h>

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0

IntTillingMult:: IntTillingMult(Target* target, int wInX, int wInY,float ratio) :
	Operator(target), target(target),wInX(wInX), wInY(wInY), wOut(wInX + wInY),ratio(ratio){
 
	ostringstream name;

	name <<"IntMultiplier_"<<wInX<<"_"<<wInY;
	setName(name.str());
		
	setCopyrightString("Bogdan Pasca, Sebastian Banescu , Radu Tudoran 2009");
	
	addInput ("X", wInX);
	addInput ("Y", wInY);
	addOutput("R", wOut); /* wOut = wInX + wInY */
		
	nrDSPs = estimateDSPs();
	cout<< nrDSPs <<endl;
	cout<<"Extra width:= "<<getExtraWidth()<<" \n Extra height"<<getExtraHeight()<<endl;
		
	
	initTiling(globalConfig, 4);
	int x, y;
	
		for (int i=0; i<4; i++)
		{
			globalConfig[i]->getTopRightCorner(x,y);
			cout << "DSP block #" << i << " is positioned at top-right: (" << x << ", " << y << ")";
			globalConfig[i]->getBottomLeftCorner(x,y);
			cout << "and bottom-left: (" << x << ", " << y << ")" << endl;
		}
	}
	
IntTillingMult::~IntTillingMult() {
}


int IntTillingMult::estimateDSPs()
{
	if(1==ratio) 
	{
		return target->getNumberOfDSPs();
	}
	else
	{	
		float temp = (target->multiplierLUTCost(wInX,wInY) * ratio)  /   ((1- ratio)  *  target->getEquivalenceSliceDSP() ) ;
		if(temp - ((int)temp)  > 0 ) // if the estimated number of dsps is a number with nonzero real part than we will return the integer number grather then this computed value. eg. 2.3 -> 3
		{
		temp++;	
		}
		
		return ((int)  temp) ;
	}
	
	return 0;
}


int  IntTillingMult::getExtraHeight()
{
int x,y;	
target->getDSPWidths(x,  y);
float temp = ratio * 0.75 * ((float) y);
return ((int)temp);
}

	
int  IntTillingMult::getExtraWidth()
{
int x,y;	
target->getDSPWidths(x,y);
float temp = ratio * 0.75 * ((float) x);
return ((int)temp);
}

void IntTillingMult::tillingAlgorithm()
{
	
}



void IntTillingMult::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");

	mpz_class svR = svX * svY;

	tc->addExpectedOutput("R", svR);
}

void IntTillingMult::buildStandardTestCases(TestCaseList* tcl){
	
}


bool IntTillingMult::checkOverlap(DSP** config, int index)
{	
	int xtr1, ytr1, xbl1, ybl1, xtr2, ytr2, xbl2, ybl2;
	//int nrdsp = sizeof(config)/sizeof(DSP*);
	config[index]->getTopRightCorner(xtr1, ytr1);
	config[index]->getBottomLeftCorner(xbl1, ybl1);
	
	cout << tab << tab << "checkOverlap: ref is block #" << index << ". Top-right is at (" << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
	
	for (int i=0; i<nrDSPs; i++)
		if ((config[i] != NULL) && (i != index))
		{
			config[i]->getTopRightCorner(xtr2, ytr2);
			config[i]->getBottomLeftCorner(xbl2, ybl2);
			
			cout << tab << tab << "checkOverlap: comparing with block #" << i << ". Top-right is at (" << xtr2 << ", " << ytr2 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
			
			if (((xtr2 < xbl1) && (ytr2 < ybl1) && (xtr2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the upper and/or right part(s) of config[i]
				((xbl2 < xbl1) && (ybl2 < ybl1) && (xbl2 > xtr1) && (ybl2 > ytr1)) || // config[index] overlaps the bottom and/or left part(s) of config[i]
				((xtr2 < xbl1) && (ybl2 < ybl1) && (xtr2 > xtr1) && (ybl2 > ytr1)) || // config[index] overlaps the upper and/or left part(s) of config[i]
				((xbl2 < xbl1) && (ytr2 < ybl1) && (xbl2 > xtr1) && (ytr2 > ytr1)) || // config[index] overlaps the bottom and/or right part(s) of config[i]
				((xbl2 > xbl1) && (ybl2 < ybl1) && (ytr2 > ytr1) && (xtr2 < xtr1)) || // config[index] overlaps the center part of config[i]
				
				((xtr1 < xbl2) && (ytr1 < ybl2) && (xtr1 > xtr2) && (ytr1 > ytr2)) || // config[i] overlaps the upper and/or right part(s) of config[index]
				((xbl1 < xbl2) && (ybl1 < ybl2) && (xbl1 > xtr2) && (ybl1 > ytr2)) || // config[i] overlaps the bottom and/or left part(s) of config[index]
				((xtr1 < xbl2) && (ybl1 < ybl2) && (xtr1 > xtr2) && (ybl1 > ytr2)) || // config[i] overlaps the upper and/or left part(s) of config[index]
				((xbl1 < xbl2) && (ytr1 < ybl2) && (xbl1 > xtr2) && (ytr1 > ytr2)) || // config[i] overlaps the bottom and/or right part(s) of config[index]
				((xbl1 > xbl2) && (ybl1 < ybl2) && (ytr1 > ytr2) && (xtr1 < xtr2))    // config[i] overlaps the center part of config[index]
				)
				return true;
		}		
	cout << tab << tab << "checkOverlap: return false" << endl;	
	return false;
}

bool IntTillingMult::move(DSP** config, int index)
{
	int xtr1, ytr1, xbl1, ybl1;
	int w, h;
	target->getDSPWidths(w,h);
	
	
	config[index]->getTopRightCorner(xtr1, ytr1);
	config[index]->getBottomLeftCorner(xbl1, ybl1);
	
	do{
		// move down one unit
		ytr1++;
		ybl1++;
		
		if (ybl1 > h+2*getExtraHeight()-1) // the DSP block has reached the bottom limit of the tiling grid
		{
			// move to top of grid and one unit to the left 
			xtr1++;
			xbl1++;
			ytr1 = 0;
			ybl1 = h-1;
			
			if (xbl1 > w+2*getExtraWidth()-1) // the DSP block has reached the left limit of the tiling grid
				return false;
		}			
		
		config[index]->setTopRightCorner(xtr1, ytr1);
		config[index]->setBottomLeftCorner(xbl1, ybl1);
	}while (checkOverlap(config, index));
	
	// set the current position of the DSP block within the tiling grid
	config[index]->setTopRightCorner(xtr1, ytr1);
	config[index]->setBottomLeftCorner(xbl1, ybl1);
	return true;
}

void IntTillingMult::replace(DSP** config, int index)
{
	int xtr1, ytr1, xbl1, ybl1;
	int w, h;
	target->getDSPWidths(w,h);
	
	xtr1 = ytr1 = 0;
	xbl1 = w-1;
	ybl1 = h-1;
	
	config[index]->setTopRightCorner(xtr1, ytr1);
	config[index]->setBottomLeftCorner(xbl1, ybl1);
	
	cout << tab << "replace : DSP width is " << w << ", height is " << h << endl; 
	while (checkOverlap(config, index))
	{
		// move down one unit
		ytr1++;
		ybl1++;
		
		cout << tab << "replace : moved DSP one unit down." << endl;
 		if (ybl1 > h+2*getExtraHeight()) // the DSP block has reached the bottom limit of the tiling grid
		{
			// move to top of grid and one unit to the left 
			xtr1++;
			xbl1++;
			ytr1 = 0;
			ybl1 = h-1;
			cout << tab << "replace : moved DSP up and one unit left." << endl;
		}			
		
		config[index]->setTopRightCorner(xtr1, ytr1);
		config[index]->setBottomLeftCorner(xbl1, ybl1);
		cout << tab << "replace : Top-right is at ( " << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
	}
	
	// set the current position of the DSP block within the tiling grid
	config[index]->setTopRightCorner(xtr1, ytr1);
	config[index]->setBottomLeftCorner(xbl1, ybl1);
}

void IntTillingMult::initTiling(DSP** config, int dspCount)
{
	config = new DSP*[nrDSPs];
	
	for (int i=0; i<dspCount; i++)
	{
		cout << "initTiling : iteration #" << i << endl; 
		config[i] = target->createDSP();
		replace(config, i);
	}
}
	
