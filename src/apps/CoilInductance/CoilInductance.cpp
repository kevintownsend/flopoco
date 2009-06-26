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
#include "../../utils.hpp"
#include "../../Operator.hpp"

#include "CoilInductance.hpp"

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


CoilInductance::CoilInductance(Target* target, int LSBI, int MSBI, int LSBO, int MSBO, char *filepath) :
	Operator(target), MSBI(MSBI), LSBI(LSBI), LSBO(LSBO), MSBO(MSBO) ,filepath(filepath){
	
	if ((MSBI < LSBI)){
		cerr << 
			" CoilInductance: Input constraint LSBI <= MSBI not met."<<endl;
		exit (EXIT_FAILURE);
	}
	
	if ((MSBO < LSBO)){
		cerr << 
			" CoilInductance: Input constraint LSBO <= MSBO not met."<<endl;
		exit (EXIT_FAILURE);
	}
		
	ostringstream name; 
	name <<"CoilInductance_"<<abs(LSBI)<<"_"<<abs(MSBI)<<"_"<<abs(LSBO)<<"_"<<abs(MSBO);
	setName(name.str());
	
		
	}

CoilInductance::~CoilInductance() {
}

void CoilInductance::emulate(TestCase * tc)
{
}

void CoilInductance::buildStandardTestCases(TestCaseList* tcl){
	
}