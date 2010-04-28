/*
  Function Evaluator for FloPoCo

  Authors: Bogdan Pasca, Mioara Joldes

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

  */

#ifdef HAVE_SOLLYA

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


	FunctionEvaluator::FunctionEvaluator(Target* target, string func, int wInX, int wOutX, int n):
		Operator(target), wInX_(wInX), wOutX_(wOutX){

		ostringstream name;
		srcFileName="FunctionEvaluator";
		
		name<<"FunctionEvaluator"; 
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca, Mioara Joldes (2010)");		

    pf=new  PiecewiseFunction(func);
		tg = new TableGenerator(target, pf, wInX, wOutX+1, n);
		oplist.push_back(tg);
		combinatorialOperator = false;
		
		YVar* y = new YVar(wInX - tg->wIn, -tg->wIn+1);
		
		
		pe = new PolynomialEvaluator(target, tg->getCoeffParamVector(), y, wOutX+1, tg->getMaxApproxError() );
		
		oplist.push_back(pe);
		wR = pe->getRWidth();
		weightR = pe->getRWeight();


		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addInput ("X", wInX);

		
		vhdl << tab << declare("addr", tg->wIn) << " <= X"<<range(wInX-1, wInX-tg->wIn)<<";"<<endl;

		nextCycle();/////////////////////////////////// The Coefficent ROM has a registered iunput
		
		inPortMap ( tg, "X", "addr");
		outPortMap ( tg, "Y", "Coef");
		
		vhdl << instance ( tg, "GeneratedTable" );
		
		syncCycleFromSignal("Coef");
		nextCycle();/////////////////////////////////// The Coefficent ROM has a registered output
		
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
//		nextCycle();/////////////////
		//now we round
//		vhdl << tab << declare( "op", wOutX+1 + (pe->getOutputWeight())) << " <= Rpe"<<range(pe->getOutputSize()-1, pe->getOutputSize() - (wOutX+pe->getOutputWeight()+ 1)) << ";" << endl;
		
//		IntAdder *ia = new IntAdder(target, wOutX+1 + (pe->getOutputWeight()+1));
//		oplist.push_back(ia);
//		
//		inPortMap    (ia, "X", "op1_rnd");
//		inPortMapCst (ia, "Y", zg(wOutX+1 + (pe->getOutputWeight()+1), 0));
//		inPortMapCst (ia, "Cin", "'1'");
//		outPortMap   (ia, "R", "postRoundR");   

//		vhdl << instance( ia, "Final_Round");
//		syncCycleFromSignal( "postRoundR" );  
		
		addOutput("R", pe->getRWidth());
		vhdl << tab << "R <= Rpe;" << endl;  
	}

	FunctionEvaluator::~FunctionEvaluator() {
	}
	
	void FunctionEvaluator::emulate(TestCase* tc)
	{
		
    int  nrFunctions=(pf->getPiecewiseFunctionArray()).size();

      Function *f=pf->getPiecewiseFunctionArray(0);
    
    if (nrFunctions !=1) {
      cout<<"Warning: we are dealing with a piecewise function; only the first branch will be evaluated"<<endl;    
    }
    mpz_class rd; // rounded down
		mpz_class ru;
		 
    mpz_class svX = tc->getInputValue("X");
    mpfr_t mpX, mpR;
  	mpfr_init2(mpX,165);
    mpfr_init2(mpR,165);
    int outSign = 0;
    /* Convert a random signal to an mpfr_t in [0,1[ */
		mpfr_set_z(mpX, svX.get_mpz_t(), GMP_RNDN);
		mpfr_div_2si(mpX, mpX, wInX_, GMP_RNDN);

		/* Compute the function */
		f->eval(mpR, mpX);

		/* Compute the signal value */
		if (mpfr_signbit(mpR))
			{
				outSign = 1;
				mpfr_abs(mpR, mpR, GMP_RNDN);
			}
		mpfr_mul_2si(mpR, mpR, wOutX_, GMP_RNDN);

		/* NOT A TYPO. Function Evaluator only guarantees faithful
		 * rounding, so we will round down here,
		 * add both the upper and lower neighbor.
		 */
		mpfr_get_z(rd.get_mpz_t(), mpR, GMP_RNDD);
		ru = rd + 1;
  
    tc->addExpectedOutput("RD", rd);
    tc->addExpectedOutput("RU", ru);
    mpfr_clear(mpX);
    mpfr_clear(mpR);
	}

}
	
#endif //HAVE_SOLLYA	
