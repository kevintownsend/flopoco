/*
 * An integer comparator unit
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

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntComparator.hpp"

using namespace std;



namespace flopoco{


	IntComparator::IntComparator(Target* target, int wIn, int criteria, map<string, double> inputDelays) :
		Operator(target, inputDelays), wIn_(wIn), criteria_(criteria) {
	
		// -------- Parameter set up -----------------
		srcFileName = "IntComaprator";
		if (criteria<-2 && criteria>2)
			criteria_=0;
		
		ostringstream name;
		switch(criteria){
			case -2: name << "IntComparator_"<< wIn<<"_"<<"less";   break;
			case -1: name << "IntComparator_"<< wIn<<"_"<<"leq";    break;
			case  0: name << "IntComparator_"<< wIn<<"_"<<"eq";     break;
			case  1: name << "IntComparator_"<< wIn<<"_"<<"geq";    break;
			case  2: name << "IntComparator_"<< wIn<<"_"<<"greater";break;
			default: name << name << "IntComparator_"<< wIn<<"_"<<"eq";
		}
		setName(name.str());
		setCopyrightString("Bogdan Pasca (2010)");

		addInput ("X", wIn_);
		addInput ("Y", wIn_);
		addOutput("R"); 

		switch(criteria_){
			case -2: vhdl << tab << "R <= '1' when (X <  Y) else '0';"<<endl; break;
			case -1: vhdl << tab << "R <= '1' when (X <= Y) else '0';"<<endl; break;
			case  0: vhdl << tab << "R <= '1' when (X =  Y) else '0';"<<endl; break;
			case  1: vhdl << tab << "R <= '1' when (X >  Y) else '0';"<<endl; break;
			case  2: vhdl << tab << "R <= '1' when (X >= Y) else '0';"<<endl; break;
			default:;
		}
		
		
		}

	IntComparator::~IntComparator() {
	}
	
	void IntComparator::emulate(TestCase* tc){
		mpz_class svX = tc->getInputValue ( "X" );
		mpz_class svY = tc->getInputValue ( "Y" );
		
		mpz_class svR;
		
		switch(criteria_){
			case -2: if ( svX <  svY ) svR = 1; else svR = 0; break;
			case -1: if ( svX <= svY ) svR = 1; else svR = 0; break;
			case  0: if ( svX == svY ) svR = 1; else svR = 0; break;
			case  1: if ( svX >= svY ) svR = 1; else svR = 0; break;
			case  2: if ( svX > svY  ) svR = 1; else svR = 0; break;
			default:; 
		}
		
		tc->addExpectedOutput ( "R", svR );
	}


}
