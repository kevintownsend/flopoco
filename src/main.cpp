/*
  the FloPoCo command-line interface

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Authors : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr
			Bogdan Pasca, Bogdan.Pasca@ens-lyon.org

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA-Lyon
  2008-2014.
  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <mpfr.h>
#include <sollya.h>
#include <cstdlib>


#include "FloPoCo.hpp"
#include "utils.hpp"
#include "main.hpp"

using namespace std;
using namespace flopoco;

int main(int argc, char* argv[] )
{
	try {

		Fix2FP::registerFactory();
		FP2Fix::registerFactory();
		InputIEEE::registerFactory();
		OutputIEEE::registerFactory();

		Shifter::registerFactory();
		LZOC::registerFactory();
		LZOCShifterSticky::registerFactory();

		IntAdder::registerFactory();
#if 0 // Plug them for debug purpose only
		IntAdderClassical::registerFactory();
		IntAdderAlternative::registerFactory();
		IntAdderShortLatency::registerFactory();
		IntAdderSpecific::registerFactory();
		LongIntAdderAddAddMuxGen1::registerFactory();
		LongIntAdderAddAddMuxGen2::registerFactory();
		LongIntAdderCmpAddIncGen1::registerFactory();
		LongIntAdderCmpAddIncGen2::registerFactory();
		LongIntAdderCmpCmpAddGen1::registerFactory();
		LongIntAdderCmpCmpAddGen2::registerFactory();
		LongIntAdderMuxNetwork::registerFactory();
		IntComparatorSpecific::registerFactory();
#endif
		IntComparator::registerFactory();
		IntDualSub::registerFactory();
		IntMultiplier::registerFactory();
		IntSquarer::registerFactory();

		FPConstMult::registerFactory();
		FPRealKCM::registerFactory();
		IntConstDiv::registerFactory();
		FPConstDiv::registerFactory();
		FixFunctionByTable::registerFactory();
		FixFunctionBySimplePoly::registerFactory();
		FixFunctionByPiecewisePoly::registerFactory();
		FixFunctionByMultipartiteTable::registerFactory();
		BasicPolyApprox::registerFactory();
		PiecewisePolyApprox::registerFactory();
		FixRealKCM::registerFactory();
		TestBench::registerFactory();
		Wrapper::registerFactory();
		FPAdd::registerFactory();
		FPAddSub::registerFactory();
		FPAddDualPath::registerFactory();
		FPAdd3Input::registerFactory();
		FPAddSinglePath::registerFactory();
		FPMult::registerFactory();
		//FPMultKaratsuba::registerFactory();
		FPSquare::registerFactory();
		FPDiv::registerFactory();
		FPDiv::NbBitsMinRegisterFactory();
		FPSqrt::registerFactory();

		FPLargeAcc::registerFactory();
		LargeAccToFP::registerFactory();
		FPDotProduct::registerFactory();
		
		FPExp::registerFactory();
		IterativeLog::registerFactory();
		FPPow::registerFactory();
		FixSinCos::registerFactory();
		CordicSinCos::registerFactory();
		FixAtan2::registerFactory();
		//FixedComplexAdder::registerFactory();
		FixFIR::registerFactory();
		FixSOPC::registerFactory();
		FixIIR::registerFactory();

		// hidden for now
		// Fix2DNorm::registerFactory();

		TargetModel::registerFactory();
		// Uncomment me to play within FloPoCo operator development
	  UserDefinedOperator::registerFactory();
	}
	catch (const std::exception &e) {
		cerr << "Error while registering factories: " << e.what() <<endl;
		exit(EXIT_FAILURE);
	}
	// cout << UserInterface::getFactoryCount() << " factories registered " << endl ;

	try {
		UserInterface::main(argc, argv);
	}
	catch (const std::exception &e) {
		cerr << "Error in main: " << e.what() <<endl;
		exit(EXIT_FAILURE);
	}

	return 0;
}



#if 0

	//------------ Resource Estimation --------------------------------
	int reLevel;
	bool resourceEstimationDebug = false;
	//-----------------------------------------------------------------

	//------------------ Floorplanning --------------------------------
	bool floorplanning = false;
	bool floorplanningDebug = false;
	ostringstream floorplanMessages;
	//-----------------------------------------------------------------


	//------------------------ Resource Estimation ---------------------
	for (vector<Operator*>::iterator it = oplist->begin(); it!=oplist->end(); ++it) {
		Operator* op = *it;

		if(reLevel!=0){
			if(op->reActive)
				cerr << op->generateStatistics(reLevel);
			else{
				cerr << "Resource estimation option active for an operator that has NO estimations in place." << endl;
			}
		}
	}
	//------------------------------------------------------------------


	//------------------------------------------------------------------
#endif




