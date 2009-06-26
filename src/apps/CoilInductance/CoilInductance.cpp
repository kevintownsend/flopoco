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
	
	setCopyrightString("Bogdan Pasca, Radu Tudoran (2009)");
	
	addOutput("O",MSBO-LSBO);
	
	//Counters for addressing the memories and for frequency division
	
	//Memories instantiation
	
	//Computing the segments x1-x0 x3-x2 y1-y0 y3-y2 z1-z0 z3-z2
	
	//syncCycleFromSignal("????"); sincronization with memories
	
	inputWidth= MSBI-LSBI;
	
	vhdl<<tab<<declare("signal_x0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_x1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_x2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_x3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	
	vhdl<<tab<<declare("signal_y0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_y1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_y2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_y3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	vhdl<<tab<<declare("signal_z0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_z1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_z2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_z3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
		
	
	//performing x1-x0
	
	vhdl<<tab<<declare("not_x0",inputWidth)<<"<= not( "<<use("signal_x0")<<");"<<endl;
	
	segment1X = new IntAdder(target,inputWidth);
	segment1X->changeName(getName()+"segment1X");
	oplist.push_back(segment1X);
	inPortMap  (segment1X, "X", use("signal_x1"));
	inPortMap  (segment1X, "Y", use("not_x0"));
	inPortMapCst(segment1X, "Cin", "'1'");
	outPortMap (segment1X, "R","segmentX1mX0");
	vhdl << instance(segment1X, "segment1X");
	
	syncCycleFromSignal("segmentX1mX0");
	
	
	
	
	}

CoilInductance::~CoilInductance() {
}

void CoilInductance::emulate(TestCase * tc)
{
}

void CoilInductance::buildStandardTestCases(TestCaseList* tcl){
	
}

