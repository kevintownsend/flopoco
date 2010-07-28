/*
  An integer squarer for FloPoCo
 
  Author: Bogdan Pasca

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.
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

namespace flopoco{

	extern vector<Operator *> oplist;

	IntSquarer::IntSquarer(Target* target, int wIn, map<string, double> inputDelays):
		Operator(target), wIn_(wIn), inputDelays_(inputDelays)
	{
		ostringstream name;
		name << "IntSquarer_" << wIn_;
		setName(name.str());
		setCopyrightString("Bogdan Pasca (2009)");

		srcFileName = "IntSquarer";

		// Set up the IO signals
		addInput ("X"  , wIn_);
		addOutput("R"  , 2*wIn_);

		setCriticalPath( getMaxInputDelays(inputDelays) );

		if (wIn <= 17 ) {
			vhdl << tab << declare( "sX", wIn) << " <= X;" << endl;
			vhdl << tab << declare( "sY", wIn) << " <= X;" << endl;
			
			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );
			vhdl << tab << "R <= sX * sY;" << endl; 
			outDelayMap["R"] = getCriticalPath();
		}
		else if ((wIn>17) && (wIn<=34)) {

			vhdl << declare("x0_16",17) << " <= X" << range(16,0) << ";" << endl;
			if (wIn<34)
				vhdl << declare("x17_33",17) << " <= "<<zg(34-wIn,0) << " & X" << range(wIn-1,17) << ";" << endl;
			else
				vhdl << declare("x17_33",17) << " <= X" << range(33,17) << ";" << endl;
						
			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );
			vhdl << declare("p0",34) << " <= x0_16 * x0_16;" <<endl;
			vhdl << declare("f0",18) << " <= p0" << range(17,0) << ";" << endl;
			vhdl << declare("p1",34) << " <= x17_33 * x0_16;" <<endl;

			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << declare("s1",34) << " <= p1 + ( \"0\" & p0" << range(33,18) <<");" <<endl;
			vhdl << declare("f1",16) << " <= s1" << range(15,0) << ";" << endl;
			vhdl << declare("p2",34) << " <= x17_33 * x17_33;" <<endl;

			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << declare("s2",34) << " <= p2 + ( \"0\" & s1" << range(33,16) <<");" <<endl;

			if (wIn<34)
				vhdl << "R <= s2" << range(2*wIn-34-1,0) << " & f1 & f0;" << endl;
			else
				vhdl << "R <= s2 & f1 & f0;" << endl;

			outDelayMap["R"] = getCriticalPath();
		}
		else if ((wIn>34) && (wIn<=51)) {
			// --------- Sub-components ------------------
			intadder = new IntAdder(target, 84);
			oplist.push_back(intadder);
			if (wIn<51)  
				vhdl << declare("sigX",51) << "<= "<<zg(51-wIn,0) <<" & X;"<<endl;
			else
				vhdl << declare("sigX",51) << "<= X;"<<endl;
			vhdl << declare("x0_16_sqr",34) << "<= sigX" << range(16,0) << " * sigX" << range(16,0)<<";"<<endl;
			vhdl << declare("x17_33_sqr",34) << "<= sigX" << range(33,17) << " * sigX" << range(33,17)<<";"<<endl;
			vhdl << declare("x34_50_sqr",34) << "<= sigX" << range(50,34) << " * sigX" << range(50,34)<<";"<<endl;
			vhdl << declare("x0_16_x17_33",34) << "<= sigX" << range(16,0) << " * sigX" << range(33,17)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x0_16_x34_50",34) << "<= x0_16_x17_33" << range(33,17) << " + sigX" << range(16,0)   << " * sigX" << range(50,34)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("x17_33_x34_50",34) << "<= x0_16_x34_50" << range(33,17) << " + sigX" << range(33,17) << " * sigX" << range(50,34)<<";"<<endl;
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("op1",84) << "<= x34_50_sqr & x17_33_sqr & x0_16_sqr" << range(33,18)<<";"<<endl;
			vhdl << declare("op2",84) << "<= \"0000000000000000\" & x17_33_x34_50 & x0_16_x34_50" << range(16,0)<<" & x0_16_x17_33" << range(16,0)<<";"<<endl;
			
			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << instance(intadder, "ADDER1");
			
			syncCycleFromSignal("adderOutput", false);
			
			vhdl << "R <= adderOutput" << range(2*wIn-19,0) << " & x0_16_sqr" << range(17,0)<<";"<<endl;
		}
		else if ((wIn>51) && (wIn<=68) && (wIn!=53) ) {
			// --------- Sub-components ------------------

			if (wIn<68)  
				vhdl << declare("sigX",68) << "<= "<<zg(68-wIn,0) <<" & X;"<<endl;
			else
				vhdl << declare("sigX",68) << "<= X;"<<endl;
				
			setCriticalPath( getMaxInputDelays(inputDelays) );	

			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() ); 
			vhdl << declare("x0_16_x17_33",34) << "<= sigX" << range(16,0) << " * sigX" << range(33,17)<<";"<<endl;
			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << declare("x0_16_sqr",36) << "<= (\"00\" & x0_16_x17_33" << range(15,0)<<" & \"000000000000000000\") + " 
				  << "( \"0\" & sigX" << range(16,0) << ") * (\"0\" & sigX" << range(16,0)<<");"<<endl;
			vhdl << declare("x17_33_x34_50",34) << "<= sigX" << range(33,17) << " * sigX" << range(50,34)<<";"<<endl;
			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << declare("x17_33_sqr",36) << "<= (\"00\" & x17_33_x34_50" << range(15,0)<<" & x0_16_x17_33" << range(33,16) <<") + x0_16_sqr(34) + "
				  << "( \"0\" & sigX" << range(33,17) << ") * (\"0\" & sigX" << range(33,17)<<");"<<endl;
			vhdl << declare("x51_67_x34_50",34) << "<= sigX" << range(67,51) << " * sigX" << range(50,34)<<";"<<endl;
			vhdl << declare("x_0_16_34_50",34) << " <= ( sigX" << range(16,0) << ") * (sigX" << range(50,34)<<");"<<endl;
			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << declare("x34_50_sqr",36) << "<= (\"00\" & x51_67_x34_50" << range(15,0)<<" & x17_33_x34_50" << range(33,16) <<") + x17_33_sqr(34) + "
				  << "( \"0\" & sigX" << range(50,34) << ") * (\"0\" & sigX" << range(50,34)<<");"<<endl;
			vhdl << declare("x_0_16_51_67_pshift", 34) << " <= x_0_16_34_50" << range(33,17) << " + "
				  << "( sigX" << range(16,0) << ") * (sigX" << range(67,51)<<");"<<endl;
			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << declare("x51_67_sqr",34) << "<= ( \"00000000000000\" & x51_67_x34_50" << range(33,16) <<") + x34_50_sqr(34) + "
				  << "( sigX" << range(67,51) << ") * (sigX" << range(67,51)<<");"<<endl;
			vhdl << declare("x_17_33_51_67_pshift", 34) << " <= x_0_16_51_67_pshift" << range(33,17) << " + "
				  << "( sigX" << range(33,17) << ") * (sigX" << range(67,51)<<");"<<endl;

//			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("op1",101) << "<= x51_67_sqr & x34_50_sqr" << range(33,0) << " & x17_33_sqr" << range(33,1) <<  ";"<<endl;
			vhdl << declare("op2",101) << "<="<< zg(101-68,0)<<" & x_17_33_51_67_pshift & x_0_16_51_67_pshift" << range(16,0)<<" & x_0_16_34_50" << range(16,0)<<";"<<endl;
			intadder = new IntAdder(target, 101, inDelayMap("X", target->DSPToLogicWireDelay() + getCriticalPath() ));
			oplist.push_back(intadder);

			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << instance(intadder, "ADDER1");
			
			syncCycleFromSignal("adderOutput", false);
			setCriticalPath( intadder->getOutputDelay("R"));
						
			vhdl << "R <= adderOutput" << range(2*wIn-36,0) << " & x17_33_sqr" << range(0,0) << " & x0_16_sqr" << range(33,0)<<";"<<endl;
			outDelayMap["R"] = getCriticalPath();
		}
		else if (wIn==53){
			//TODO -> port to new pipeline framework
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


			vhdl << declare("op1mul2",53) << "<= (\"00\" & sigX" << range(50,0) <<") when sigX(51)='1' else "<<zg(53,0) << ";"<<endl;
			vhdl << declare("op2mul2",53) << "<= (\"0\" & sigX" << range(50,0) <<" & \"0\") when sigX(52)='1' else "<<zg(53,0) << ";"<<endl;

			nextCycle(); ////////////////////////////////////////////////

			inPortMap(intadd2, "X", "op1mul2");
			inPortMap(intadd2, "Y", "op2mul2");
			inPortMapCst(intadd2, "Cin", "'0'");
			outPortMap(intadd2, "R", "x51_52_times_x_0_50");
			vhdl << instance(intadd2, "MULT2");

			syncCycleFromSignal("out_Squarer_51",true);
			
			vhdl << declare("x51_52_sqr",4) << " <= sigX" << range(52,51) << " * sigX" << range(52,51) <<";"<<endl;
			
			nextCycle(); ////////////////////////////////////////////////
			vhdl << declare("op1",54) << "<= x51_52_sqr & out_Squarer_51" << range(101,52) <<  ";"<<endl;
			vhdl << declare("op2",54) << "<="<< zg(1,0)<<" & x51_52_times_x_0_50;"<<endl;

			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << instance(intadder, "ADDER54");
			
			syncCycleFromSignal("adderOutput", false);
			
			vhdl << "R <= adderOutput & out_Squarer_51" << range(51,0)<<";"<<endl;
		} else {
			cerr << " For the moment IntSquarer does not support inputs larger than 68 bits. " << endl;
			exit (EXIT_FAILURE);
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


}
