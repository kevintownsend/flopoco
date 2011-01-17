/*
An integer adder for FloPoCo

It may be pipelined to arbitrary frequency.
Also useful to derive the carry-propagate delays for the subclasses of Target

Authors:  Bogdan Pasca, Florent de Dinechin

This file is part of the FloPoCo project
developed by the Arenaire team at Ecole Normale Superieure de Lyon

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
CeCILL license, 2008-2010.
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntAdder.hpp"
#include "IntAdderClassical.hpp"
#include "IntAdderAlternative.hpp"
#include "IntAdderShortLatency.hpp"
// 
using namespace std;

namespace flopoco {
	extern vector<Operator*> oplist;
	
	IntAdder::IntAdder ( Target* target, int wIn, map<string, double> inputDelays, int optimizeType, bool srl, int implementation):
	Operator ( target, inputDelays, false ), wIn_ ( wIn )  {
		ostringstream name;
		srcFileName="IntAdder";
		setCopyrightString ( "Bogdan Pasca, Florent de Dinechin (2008-2010)" );
		
		name << "IntAdder_" << wIn_<<"_f"<<target->frequencyMHz()<<"_uid"<<getNewUId();
		setName ( name.str() );
		
		// Set up the IO signals
		addInput ( "X"  , wIn_, true );
		addInput ( "Y"  , wIn_, true );
		addInput ( "Cin", 1 );
		addOutput ( "R"  , wIn_, 1 , true );
		
		REPORT(INFO, "Implementing IntAdder " << wIn);
		
		Operator* intAdderInstantiation;
		intAdderInstantiation = new IntAdderClassical(target, wIn, name.str() , inputDelays, optimizeType, srl);
		addImplementationList.push_back(intAdderInstantiation);
		
		intAdderInstantiation = new IntAdderAlternative(target, wIn, name.str() , inputDelays, optimizeType, srl);
		addImplementationList.push_back(intAdderInstantiation);

		intAdderInstantiation = new IntAdderShortLatency(target, wIn, name.str() , inputDelays, optimizeType, srl);
		addImplementationList.push_back(intAdderInstantiation);
				
		int currentCost = 16384;
		for (unsigned j=0; j< addImplementationList.size(); j++)
			if (currentCost > addImplementationList[j]->getOperatorCost()){
				currentCost = addImplementationList[j]->getOperatorCost();
				selectedVersion = j;
			}
		
		REPORT(INFO, "Selected implementation is "<<selectedVersion<<" with cost="<<currentCost);
		
		addImplementationList[selectedVersion]->setuid(getuid()); //the selected implemetation becomes this operator 
		
		oplist.push_back(addImplementationList[selectedVersion]); //the code of the selected implementation 
		outDelayMap["R"] = addImplementationList[selectedVersion]->getOutputDelay("R"); //populate output delays
		setCycle(addImplementationList[selectedVersion]->getPipelineDepth());
		
		//cleanup; clear the oplist of the components that will be unused, and the components used therein 
		for (unsigned j=0; j< addImplementationList.size(); j++)
			if ( int(j)!= selectedVersion){
				REPORT(INFO, "deleting version "<<int(j));
				cleanup(&oplist, addImplementationList[j]);
			}
		REPORT(INFO, "Finished implementing the adder");
	}
	
	/**************************************************************************/
	IntAdder::~IntAdder() {
	}
	
	void IntAdder::outputVHDL(std::ostream& o, std::string name) {

	
	}
	
	/******************************************************************************/
	void IntAdder::emulate ( TestCase* tc ) {
		mpz_class svX = tc->getInputValue ( "X" );
		mpz_class svY = tc->getInputValue ( "Y" );
		mpz_class svC = tc->getInputValue ( "Cin" );
		
		mpz_class svR = svX + svY + svC;
		// Don't allow overflow
		mpz_clrbit ( svR.get_mpz_t(),wIn_ );
		
		tc->addExpectedOutput ( "R", svR );
	}
	
        void IntAdder::changeName(std::string operatorName){
                Operator::changeName(operatorName);
		addImplementationList[selectedVersion]->changeName(operatorName);
        }
}


