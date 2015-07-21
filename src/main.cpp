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
#include "FPDivSqrt/Tools/NbBitsMin.cpp"

using namespace std;
using namespace flopoco;

int main(int argc, char* argv[] )
{
	try {
		Shifter::registerFactory();
		FPAdd::registerFactory();
		FPExp::registerFactory();
		BasicPolyApprox::registerFactory();
		FixFunctionBySimplePoly::registerFactory();
		FixFunctionByPiecewisePoly::registerFactory();
		FixFunctionByTable::registerFactory();
		TestBench::registerFactory();
		NbBitsMinRegisterFactory();
		FPDiv::registerFactory();
		FPSqrt::registerFactory();
		FPAddSub::registerFactory();
		FPAddDualPath::registerFactory();
		FPAdd3Input::registerFactory();
		FPAddSinglePath::registerFactory();
		FPMult::registerFactory();
		//FPMultKaratsuba::registerFactory();
		FPSquare::registerFactory();
		IntAdder::registerFactory();
	}
	catch (std::string s) {
		cerr << "Error while registering factories: " << s <<endl;
		exit(EXIT_FAILURE);
	}
	// cout << UserInterface::getFactoryCount() << " factories registered " << endl ;

	UserInterface::main(argc, argv);

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




