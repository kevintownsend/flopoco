/*
 * An IntSquarer for FloPoCo
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
#include <cstdlib>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntSquarer.hpp"

using namespace std;
extern vector<Operator *> oplist;

IntSquarer::IntSquarer(Target* target, int wIn, map<string, double> inputDelays):
Operator(target), wIn_(wIn), inputDelays_(inputDelays)
{
	ostringstream name;
	name << "IntSquarer_" << wIn_;
	setName(name.str());
	setCopyrightString("Bogdan Pasca (2009)");

	// Set up the IO signals
	addInput ("X"  , wIn_);
	addOutput("R"  , 2*wIn_);


	if (verbose){
		cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
	}

	if (isSequential()){
		if (wIn < 17 ) {
			vhdl << tab << "R <= X * X;" << endl; 
		}
		else if ((wIn>17) && (wIn<=34)) {

			vhdl << declare("x0_16",17) << " <= "<<use("X") << range(16,0) << ";" << endl;
			if (wIn<34)
				vhdl << declare("x17_33",17) << " <= "<<zeroGenerator(34-wIn,0) << " & " <<  use("X") << range(wIn-1,17) << ";" << endl;
			else
				vhdl << declare("x17_33",17) << " <= "<<use("X") << range(33,17) << ";" << endl;
			
			nextCycle();//////////////////////////////////////////////
			
			vhdl << declare("p0",34) << " <= " << use("x0_16") << " * " << use("x0_16") << ";" <<endl;

			nextCycle();//////////////////////////////////////////////
			
			vhdl << declare("p1",34) << " <= " << use("x17_33") << " * " << use("x0_16") << " + " << "( \"0\" & "<< use("p0")<<range(33,18) <<")" << ";" <<endl;
			vhdl << declare("f0",18) << " <= " << use("p0") << range(17,0) << ";" << endl;
			
			nextCycle();//////////////////////////////////////////////

			vhdl << declare("p2",34) << " <= " << use("x17_33") << " * " << use("x17_33") << " + " << "( \"0\" & "<< use("p1")<<range(33,16) <<")" << ";" <<endl;
			vhdl << declare("f1",16) << " <= " << use("p1") << range(15,0) << ";" << endl;

			nextCycle();//////////////////////////////////////////////

			if (wIn<34)
				vhdl << "R <= " << use("p2")<<range(2*wIn-34-1,0) << " & " << use("f1") << " & " << use("f0") << ";" << endl;
			else
				vhdl << "R <= " << use("p2") << " & " << use("f1") << " & " << use("f0") << ";" << endl;

		}
		else if ((wIn>34) && (wIn<=51)) {
			// --------- Sub-components ------------------
			intadder = new IntAdder(target, 84);
			oplist.push_back(intadder);
			if (wIn<51)  
				vhdl << declare("sigX",51) << "<= "<<zeroGenerator(51-wIn,0) <<" & X;"<<endl;
			else
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
			
			vhdl << "R <= " << use("adderOutput")<<range(2*wIn-19,0) << " & " << use("x0_16_sqr")<<range(17,0)<<";"<<endl;
		}
		else if ((wIn>51) && (wIn<=68) && (wIn!=53) ) {
			// --------- Sub-components ------------------
			intadder = new IntAdder(target, 101);
			oplist.push_back(intadder);

			if (wIn<68)  
			vhdl << declare("sigX",68) << "<= "<<zeroGenerator(68-wIn,0) <<" & X;"<<endl;
			else
			vhdl << declare("sigX",68) << "<= X;"<<endl;

			vhdl << declare("x0_16_x17_33",34) << "<= "<< use("sigX")<<range(16,0) << " * " << use("sigX")<<range(33,17)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x0_16_sqr",36) << "<= " << "(\"00\" & " << use("x0_16_x17_33")<<range(15,0)<<" & \"000000000000000000\") + " 
			                                << "( \"0\" & "<< use("sigX")<<range(16,0) << ") * (\"0\" & " << use("sigX")<<range(16,0)<<");"<<endl;
			vhdl << declare("x17_33_x34_50",34) << "<= "<< use("sigX")<<range(33,17) << " * " << use("sigX")<<range(50,34)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x17_33_sqr",36) << "<= " << "(\"00\" & " << use("x17_33_x34_50")<<range(15,0)<<" & "<<use("x0_16_x17_33")<<range(33,16) <<") + " << use("x0_16_sqr")<<"(34) + "
			                                << "( \"0\" & "<< use("sigX")<<range(33,17) << ") * (\"0\" & " << use("sigX")<<range(33,17)<<");"<<endl;
			vhdl << declare("x51_67_x34_50",34) << "<= "<< use("sigX")<<range(67,51) << " * " << use("sigX")<<range(50,34)<<";"<<endl;
			vhdl << declare("x_0_16_34_50",34) << " <= " << "( "<< use("sigX")<<range(16,0) << ") * (" << use("sigX")<<range(50,34)<<");"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x34_50_sqr",36) << "<= " << "(\"00\" & " << use("x51_67_x34_50")<<range(15,0)<<" & "<<use("x17_33_x34_50")<<range(33,16) <<") + " << use("x17_33_sqr")<<"(34) + "
			                                << "( \"0\" & "<< use("sigX")<<range(50,34) << ") * (\"0\" & " << use("sigX")<<range(50,34)<<");"<<endl;
			vhdl << declare("x_0_16_51_67_pshift", 34) << " <= " << use("x_0_16_34_50")<<range(33,17) << " + "
			                                << "( "<< use("sigX")<<range(16,0) << ") * (" << use("sigX")<<range(67,51)<<");"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x51_67_sqr",34) << "<= " << "( \"00000000000000\" & "<<use("x51_67_x34_50")<<range(33,16) <<") + " << use("x34_50_sqr")<<"(34) + "
			                                << "( "<< use("sigX")<<range(67,51) << ") * (" << use("sigX")<<range(67,51)<<");"<<endl;
			vhdl << declare("x_17_33_51_67_pshift", 34) << " <= " << use("x_0_16_51_67_pshift")<<range(33,17) << " + "
			                                << "( "<< use("sigX")<<range(33,17) << ") * (" << use("sigX")<<range(67,51)<<");"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("op1",101) << "<= "<< use("x51_67_sqr")<<" & " <<  use("x34_50_sqr")<<range(33,0) << " & " << use("x17_33_sqr")<<range(33,1) <<  ";"<<endl;
			vhdl << declare("op2",101) << "<="<< zeroGenerator(101-68,0)<<" & " << use("x_17_33_51_67_pshift") << " & " << use("x_0_16_51_67_pshift")<<range(16,0)<<" & " << use("x_0_16_34_50")<<range(16,0)<<";"<<endl;
			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << instance(intadder, "ADDER1");
			
			syncCycleFromSignal("adderOutput", false);
			
			vhdl << "R <= " << use("adderOutput")<<range(2*wIn-36,0) << " & " <<use("x17_33_sqr")<<range(0,0) << " & " << use("x0_16_sqr")<<range(33,0)<<";"<<endl;
		}
		else if (wIn==53){
			//instantiate a 51bit squarer
			intsquarer = new IntSquarer(target, 51);
			oplist.push_back(intsquarer);
			
			bool tempPipelineStatus = target->isPipelined();
			bool tempDSPStatus = target->getUseHardMultipliers();
			target->setNotPipelined();
			target->setUseHardMultipliers(false);
			
			if (tempPipelineStatus) 
				target->setPipelined();
			if (tempDSPStatus)
				target->setUseHardMultipliers(true);

			intadder = new IntAdder(target, 54);
			oplist.push_back(intadder);
			intadd2 = new IntAdder(target, 53);
			oplist.push_back(intadd2);

			
			vhdl << declare("sigX",53) << "<= X;"<<endl;
			
			inPortMapCst(intsquarer, "X", "sigX(50 downto 0)");
			outPortMap(intsquarer, "R", "out_Squarer_51");
			vhdl << instance(intsquarer, "SQUARER51");


			vhdl << declare("op1mul2",53) << "<= "<< "(\"00\" & "<< use("sigX")<<range(50,0) <<") when "<<use("sigX")<<"(51)='1' else "<<zeroGenerator(53,0) << ";"<<endl;
			vhdl << declare("op2mul2",53) << "<= "<< "(\"0\" & "<< use("sigX")<<range(50,0) <<" & \"0\") when "<<use("sigX")<<"(52)='1' else "<<zeroGenerator(53,0) << ";"<<endl;

			nextCycle(); ////////////////////////////////////////////////

			inPortMap(intadd2, "X", "op1mul2");
			inPortMap(intadd2, "Y", "op2mul2");
			inPortMapCst(intadd2, "Cin", "'0'");
			outPortMap(intadd2, "R", "x51_52_times_x_0_50");
			vhdl << instance(intadd2, "MULT2");

			syncCycleFromSignal("out_Squarer_51",true);
			
			vhdl << declare("x51_52_sqr",4) << " <= " << use("sigX")<<range(52,51) << " * " << use("sigX")<<range(52,51) <<";"<<endl;
			
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("op1",54) << "<= "<< use("x51_52_sqr")<<" & " <<  use("out_Squarer_51")<<range(101,52) <<  ";"<<endl;
			vhdl << declare("op2",54) << "<="<< zeroGenerator(1,0)<<" & " <<  use("x51_52_times_x_0_50") <<";"<<endl;

			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << instance(intadder, "ADDER54");
			
			syncCycleFromSignal("adderOutput", false);
			
			vhdl << "R <= " << use("adderOutput") << " & " <<use("out_Squarer_51")<<range(51,0)<<";"<<endl;
		} else {
			cerr << " For the moment IntSquarer does not support inputs larger than 68 bits. " << endl;
			exit (EXIT_FAILURE);
		}
	}
}

IntSquarer::~IntSquarer() {
}


//void IntSquarer::outputVHDL(std::ostream& o, std::string name) {
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



void IntSquarer::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svR = svX * svX ;
	tc->addExpectedOutput("R", svR);
}


