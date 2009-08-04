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
		
	
		
	}
	
IntTillingMult::~IntTillingMult() {
}


int IntTillingMult::estimateDSPs()
{
	if(1==ratio) //all dsps can be used. A function from target must be created in order to get the maximum number of dsp
	{
		return 50;
	}
	else
	{	
		if(0==ratio)
		{
		ratio = ratio + 0.000000001;
		}
		return ((int) (target->multiplierLUTCost(wInX,wInY) * ratio  / (1- ratio) ) );
	}
	
	return 0;
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



	