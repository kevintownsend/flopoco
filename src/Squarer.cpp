/*
 * An squarer for FloPoCo
 *
 * It may be pipelined to arbitrary frequency.
 * Also useful to derive the carry-propagate delays for the subclasses of Target
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
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Squarer.hpp"

using namespace std;
extern vector<Operator *> oplist;

Squarer::Squarer(Target* target, int wIn, map<string, double> inputDelays):
Operator(target), wIn_(wIn), inputDelays_(inputDelays)
{
	ostringstream name;
	name << "Squarer_" << wIn_;
	setName(name.str());

	// Set up the IO signals
	addInput ("X"  , wIn_);
	addOutput("R"  , 2*wIn_);

	// --------- Sub-components ------------------

	intadder = new IntAdder(target, 84);
	oplist.push_back(intadder);

	if (verbose){
		cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
	}

	if (isSequential()){
		if (wIn==51) {
			vhdl << declare("sigX",51) << "<= X;"<<endl;
			vhdl << declare("x0_16_sqr",34) << "<= " << use("sigX")<<range(16,0) << " * " << use("sigX")<<range(16,0)<<";"<<endl;
			vhdl << declare("x17_33_sqr",34) << "<= " << use("sigX")<<range(33,17) << " * " << use("sigX")<<range(33,17)<<";"<<endl;
			vhdl << declare("x34_50_sqr",34) << "<= " << use("sigX")<<range(50,34) << " * " << use("sigX")<<range(50,34)<<";"<<endl;
			vhdl << declare("x0_16_x17_33",34) << "<= "<< use("sigX")<<range(16,0) << " * " << use("sigX")<<range(33,17)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x0_16_x34_50",34) << "<= "<< use("x0_16_x17_33")<<range(33,17) << " + "<< use("sigX")<<range(16,0)   << " * " << use("sigX")<<range(50,34)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x17_33_x34_50",34) << "<= "<< use("x0_16_x34_50")<<range(33,17) << " + "<< use("sigX")<<range(33,17) << " * " << use("sigX")<<range(50,34)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("op1",84) << "<= "<< use("x34_50_sqr") << " & " << use("x17_33_sqr") << " & "<< use("x0_16_sqr")<<range(33,18)<<";"<<endl;
			vhdl << declare("op2",84) << "<= \"0000000000000000\" & " << use("x17_33_x34_50") << " & " << use("x0_16_x34_50")<<range(16,0)<<" & " << use("x0_16_x17_33")<<range(16,0)<<";"<<endl;
			
			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << instance(intadder, "ADDER1");
			
			syncCycleFromSignal("adderOutput", false);
			
			vhdl << "R <= " << use("adderOutput") << " & " << use("x0_16_sqr")<<range(17,0)<<";"<<endl;
		}
	}
}

Squarer::~Squarer() {
}


//void Squarer::outputVHDL(std::ostream& o, std::string name) {
//	ostringstream signame;
//	licence(o,"Bogdan Pasca (2009)");
//	Operator::stdLibs(o);
//	outputVHDLEntity(o);
//	newArchitecture(o,name);
//	o << buildVHDLSignalDeclarations();
//	o << output
//	beginArchitecture(o);
//	o << buildVHDLRegisters();
//	o << vhdl.str();
//	endArchitecture(o);
//}



void Squarer::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svR = svX * svX ;
	tc->addExpectedOutput("R", svR);
}


