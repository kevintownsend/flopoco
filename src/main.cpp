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

// Global variables, useful in this main to avoid parameter passing


	string filename="flopoco.vhdl";
	string cl_name=""; // used for the -name option
	Target* target;
	
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

	
	Target* target = new Virtex5(); // this also creates a global operator list TODO move it to the factory, or somewhere.

	verbose=3;
	Shifter::registerFactory();
	FPExp::registerFactory();

	//	cout << UserInterface::getFactoryCount() << " factories registered " << endl ;

	cout << UserInterface::getFullDoc();

	UserInterface::parseAll(target, argc, argv);

	//cout << "Successfuly built "<< UserInterface::globalOpList.size() << " operator(s)" << endl;

	
	ofstream file;
	file.open(filename.c_str(), ios::out);
	UserInterface::outputVHDLToFile(file); 
	file.close();

	UserInterface::finalReport(cerr); 
	cerr << "Output file: " << filename <<endl;




	
		
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



