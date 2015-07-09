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
#include "src/FPDivSqrt/Tools/NbBitsMin.cpp"

using namespace std;
using namespace flopoco;


	//------------ Resource Estimation --------------------------------
	int reLevel;
	bool resourceEstimationDebug = false;
	//-----------------------------------------------------------------

	//------------------ Floorplanning --------------------------------
	bool floorplanning = false;
	bool floorplanningDebug = false;
	ostringstream floorplanMessages;
	//-----------------------------------------------------------------




int main(int argc, char* argv[] )
{
	sollya_lib_init();



	Shifter::registerFactory();
	FPExp::registerFactory();
	FPDiv::registerFactory();
	FPSqrt::registerFactory();
	FPAddSub::registerFactory();	
	FPAddDualPath::registerFactory();
	FPAdd3Input::registerFactory();
	FPAddSinglePath::registerFactory();
	FPMult::registerFactory();
	NbBitsMinRegisterFactory();

	//	cout << UserInterface::getFactoryCount() << " factories registered " << endl ;

	UserInterface::parseAll(argc, argv);

	//cout << "Successfuly built "<< UserInterface::globalOpList.size() << " operator(s)" << endl;



	UserInterface::finalReport(cerr);






#if 0

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

	sollya_lib_close();

	return 0;
}



