/*
An integer adder for FloPoCo

It may be pipelined to arbitrary frequency.
Also useful to derive the carry-propagate delays for the subclasses of Target

Authors:  Bogdan Pasca, Florent de Dinechin

This file is part of the FloPoCo project
developed by the Arenaire team at Ecole Normale Superieure de Lyon

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
2008-2010.
  All rights reserved.
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
#include "IntAdder.hpp"

#include "IntAdderClassical.hpp"
#include "IntAdderAlternative.hpp"
#include "IntAdderShortLatency.hpp"

using namespace std;
namespace flopoco {

	IntAdder::IntAdder ( Target* target, int wIn, map<string, double> inputDelays, int optimizeType, bool srl, int implementation):
	Operator ( target, inputDelays), wIn_ ( wIn )  {
		ostringstream name;
		srcFileName="IntAdder";
		setCopyrightString ( "Bogdan Pasca, Florent de Dinechin (2008-2010)" );

		name << "IntAdder_" << wIn_<<"_f"<<target->frequencyMHz()<<"_uid"<<getNewUId();

		// Set up the IO signals
		addInput ( "X"  , wIn_, true );
		addInput ( "Y"  , wIn_, true );
		addInput( "Cin");
		addOutput ( "R"  , wIn_, 1 , true );

		REPORT(DETAILED, "Implementing IntAdder " << wIn << " implementation="<<implementation);

		Operator* intAdderInstantiation;

		if (implementation == -1){ // we must explore
			intAdderInstantiation = new IntAdderClassical(target, wIn, inputDelays, optimizeType, srl);
			addImplementationList.push_back(intAdderInstantiation);

			intAdderInstantiation = new IntAdderAlternative(target, wIn, inputDelays, optimizeType, srl);
			addImplementationList.push_back(intAdderInstantiation);

//			intAdderInstantiation = new IntAdderShortLatency(target, wIn, inputDelays, optimizeType, srl);
//			addImplementationList.push_back(intAdderInstantiation);
		}else{
			switch (implementation){
				case 0:
					intAdderInstantiation = new IntAdderClassical(target, wIn, inputDelays, optimizeType, srl);
					addImplementationList.push_back(intAdderInstantiation);
					break;
				case 1:
					intAdderInstantiation = new IntAdderAlternative(target, wIn, inputDelays, optimizeType, srl);
					addImplementationList.push_back(intAdderInstantiation);
					break;
				case 2:
					intAdderInstantiation = new IntAdderShortLatency(target, wIn, inputDelays, optimizeType, srl);
					addImplementationList.push_back(intAdderInstantiation);
					break;
				default:
					intAdderInstantiation = new IntAdderClassical(target, wIn, inputDelays, optimizeType, srl);
					addImplementationList.push_back(intAdderInstantiation);
			}
		}

		int currentCost = 16384;
		selectedVersion = 0;
		for (unsigned j=0; j< addImplementationList.size(); j++)
			if (currentCost > addImplementationList[j]->getOperatorCost()){
				currentCost = addImplementationList[j]->getOperatorCost();
				selectedVersion = j;
			}

		cloneOperator(addImplementationList[selectedVersion]);
		changeName ( name.str() );

		REPORT(DETAILED, "Selected implementation for IntAdder"<< wIn << " is "<<selectedVersion<<" with cost="<<currentCost);

//		//cleanup; clear the oplist of the components that will be unused, and the components used therein
//		for (unsigned j=0; j< addImplementationList.size(); j++){
//			REPORT(DEBUG, "deleting version "<<int(j));
//			cleanup(&oplist, addImplementationList[j]);
//		}
//		REPORT(DEBUG, "Finished implementing the adder");
	}

	/**************************************************************************/
	IntAdder::~IntAdder() {
	}

	/******************************************************************************/
	void IntAdder::emulate ( TestCase* tc ) {
		// get the inputs from the TestCase
		mpz_class svX = tc->getInputValue ( "X" );
		mpz_class svY = tc->getInputValue ( "Y" );
		mpz_class svC = tc->getInputValue ( "Cin" );

		// compute the multiple-precision output
		mpz_class svR = svX + svY + svC;
		// Don't allow overflow: the output is modulo 2^wIn
		svR = svR & ((mpz_class(1)<<wIn_)-1);

		// complete the TestCase with this expected output
		tc->addExpectedOutput ( "R", svR );
	}

	OperatorPtr IntAdder::parseArguments(Target *target, const vector<string> &args) {
		int wIn = UserInterface::checkStrictlyPositiveInt(args, "wIn");
		return new IntAdder(target, wIn);
	}

	void IntAdder::registerFactory(){
		UserInterface::add("IntAdder", // name
											 "Integer adder. In modern VHDL, integer addition is expressed by a + and one usually needn't define an entity for it. However, this operator will be pipelined if the addition is too large to be performed at the target frequency.",
											 "operator; floating point; floating-point adders", // categories
											 "wE(int): exponent size in bits; ",
											 "",
											 IntAdder::parseArguments
											 ) ;

	}

//    void IntAdder::changeName(std::string operatorName){
//		Operator::changeName(operatorName);
//		addImplementationList[selectedVersion]->changeName(operatorName);
//		//		cout << "changin IntAdder name to " << operatorName << endl;
//    }

}


