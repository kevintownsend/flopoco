/*
 * An integer multiplier for FloPoCo
 *
 * Authors : Bogdan Pasca
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

#include <typeinfo>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntTruncMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	IntTruncMultiplier::IntTruncMultiplier(Target* target, DSP** configuration, vector<SoftDSP*> softDSPs, int wX, int wY, int k):
		Operator(target), wX(wX), wY(wY){
 
		ostringstream name;
		name <<"IntTruncMultiplier_"<<wX<<"_"<<wY;
		setName(name.str());
	
		setCopyrightString("Bogdan Pasca 2010");
	
		addInput ("X", wX);
		addInput ("Y", wY);
		addOutput("R", wX + wY- k); 
	
		printConfiguration(configuration, softDSPs);
		
	}

	IntTruncMultiplier::~IntTruncMultiplier() {
	}
	
	
	void IntTruncMultiplier::printConfiguration(DSP** configuration, vector<SoftDSP*> softDSPs){
		if (configuration!=NULL){
			int i=0;
			int xB,xT,yB,yT;
			while(configuration[i]!=NULL){
				configuration[i]->getTopRightCorner(xT,yT);
				configuration[i]->getBottomLeftCorner(xB,yB);
				cout << "HARD DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
				i++;
			}
		}

		int xB,xT,yB,yT;
		for (int k=0; k < softDSPs.size(); k++){
			softDSPs[k]->getTopRightCorner(xT,yT);
			softDSPs[k]->getBottomLeftCorner(xB,yB);
			cout << "SOFT DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
		} 
	}

//	void IntTruncMultiplier::emulate(TestCase* tc)
//	{
//		mpz_class svX = tc->getInputValue("X");
//		mpz_class svY = tc->getInputValue("Y");

//		mpz_class svR = svX * svY;

//		tc->addExpectedOutput("R", svR);
//	}
}
