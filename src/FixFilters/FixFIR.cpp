#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixFIR.hpp"

#include "ShiftReg.hpp"
#include "FixSOPC.hpp"

using namespace std;

namespace flopoco {

	FixFIR::FixFIR(Target* target, int p_, vector<string> coeff_, bool useBitheap_, map<string, double> inputDelays) : 
		Operator(target), p(p_), coeff(coeff_), useBitheap(useBitheap_)
	{
		srcFileName="FixFIR";
		setCopyrightString ( "Louis Beseme, Florent de Dinechin (2014)" );
		useNumericStd_Unsigned();

		ostringstream name;
		name << "FixFIR_"<< p << "_uid" << getNewUId();
		setNameWithFreq( name.str() );

		n=coeff.size();
		
		//manage the critical path
		setCriticalPath(getMaxInputDelays(inputDelays));

		addInput("X", 1+p, true);

		ShiftReg *shiftReg = new ShiftReg(target, 1+p, n, inputDelays);

		addSubComponent(shiftReg);
		inPortMap(shiftReg, "X", "X");

		for(int i = 0; i<n; i++) {
			outPortMap(shiftReg, join("Xd", i), join("Y", i));
		}

		vhdl << instance(shiftReg, "shiftReg");

		// REPORT(INFO,"pouet");
		// REPORT(INFO,getCycleFromSignal("Y0", false));
		// REPORT(INFO,getCurrentCycle());
		syncCycleFromSignal("Y0");

		FixSOPC *fixSOPC = new FixSOPC(target, p, coeff, useBitheap, inputDelays);
		addSubComponent(fixSOPC);

		for(int i=0; i<n; i++) {
			inPortMap(fixSOPC, join("X",i), join("Y", i));
		}

		outPortMap(fixSOPC, "R", "Rtmp");

		vhdl << instance(fixSOPC, "fixSOPC");
		syncCycleFromSignal("Rtmp");

		addOutput("R", fixSOPC->wO, true);
		vhdl << "R <= Rtmp;" << endl;

	};

	FixFIR::~FixFIR(){

	};

	void FixFIR::emulate(TestCase * tc){

	};

	void FixFIR::buildStandardTestCases(TestCaseList* tcl){

	};

}
	
