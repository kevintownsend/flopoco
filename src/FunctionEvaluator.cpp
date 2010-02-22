/*
 * Function Evaluator for FloPoCo
 *
  * Authors: Bogdan Pasca, Mioara Joldes
  * Copyright ENS-Lyon, INRIA, CNRS, UCBL
  *
  * This file is part of the FloPoCo project developed by the Arenaire team at Ecole Normale Superieure de Lyon
  */


#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "FunctionEvaluator.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

#define DEBUGVHDL 0


	FunctionEvaluator::FunctionEvaluator(Target* target, string func, int wInX, int wOutX, int n,double xmin, double xmax, double scale):
		Operator(target){

		ostringstream name;
		srcFileName="FunctionEvaluator";
		
		name<<"FunctionEvaluator"; 
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca, Mioara Joldes (2010)");		

		tg = new TableGenerator(target, func, wInX, wOutX, n, xmin, xmax, scale);
		oplist.push_back(tg);
		
		YVar* y = new YVar(wInX - tg->wIn, -tg->wIn);
		
		pe = new PolynomialEvaluator(target, tg->getCoeffParamVector(), y, wOutX );
		oplist.push_back(pe);

		

		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addInput ("X", wInX);
		addOutput("R", pe->getOutputSize());
		
		vhdl << tab << declare("addr", tg->wIn) << " <= X"<<range(wInX-1, wInX-tg->wIn)<<";"<<endl;
		
		inPortMap ( tg, "X", "addr");
		outPortMap ( tg, "Y", "Coef");
		
		vhdl << instance ( tg, "GeneratedTable" );
		
		syncCycleFromSignal("Coef");
		
		vhdl << tab << declare ("y",y->getSize()) << " <= X"<<range(y->getSize()-1 ,0) << ";" << endl;
		
		/* get the coefficients */
		int lsb = 0, sizeS = 0;
		for (uint32_t i=0; i< pe->getCoeffParamVector().size(); i++){
			lsb += sizeS;
			sizeS = pe->getCoeffParamVector()[i]->getSize()+1;
			vhdl << tab << declare(join("a",i), sizeS ) << "<= Coef"<< range (lsb+sizeS-1, lsb) << ";" << endl;
		}
		
		inPortMap( pe, "Y", "y");
		for (uint32_t i=0; i< pe->getCoeffParamVector().size(); i++){
			inPortMap(pe,  join("a",i), join("a",i));
		}
		
		outPortMap( pe, "R", "Rpe");
		vhdl << instance( pe, "PolynomialEvaluator");
		
		syncCycleFromSignal("Rpe");
		vhdl << tab << "R <= Rpe;"<< endl;  
	}

	FunctionEvaluator::~FunctionEvaluator() {
	}
}

