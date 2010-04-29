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


	FunctionEvaluator::FunctionEvaluator(Target* target, string func, int wInX, int wOutX, int n, bool finalRounding):
		Operator(target), wInX_(wInX), wOutX_(wOutX), finalRounding_(finalRounding){

		ostringstream name;
		srcFileName="FunctionEvaluator";
		
		name<<"FunctionEvaluator"; 
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca, Mioara Joldes (2010)");
		
		REPORT(INFO, "Will try to implement the function: " << func );
		REPORT(INFO, "The input precision of x is: " << wInX << " with x in [0,1[");
		REPORT(INFO, "The required number fractional bits of the output is: "<<wOutX);
		REPORT(INFO, "The degree of the polynomial used to aproximate this function is: " << n );
		
		pf = new PiecewiseFunction(func);
		tg = new TableGenerator(target, pf, wInX, wOutX+1, n);
		oplist.push_back(tg);
		
		REPORT(INFO, "The number of intervals of the function: "<< intpow2(tg->wIn));
		REPORT(INFO, "The number of bits used for y: "<< wInX-tg->wIn << " and its weight is: " << -tg->wIn);
		
		
		/* the nmber of bits of the address */
		int addrSize = tg->wIn;
		/* the number of bits of y is */
		int ySize = wInX - tg->wIn;
		/* the weight of y considering that y is in 1.XXXX */
		int yWeight = -tg->wIn;
		
		YVar* y = new YVar(ySize, yWeight);
		pe = new PolynomialEvaluator(target, tg->getCoeffParamVector(), y, wOutX+1, tg->getMaxApproxError() );
		oplist.push_back(pe);

		wR = pe->getRWidth();
		weightR = pe->getRWeight()+1;
		
		REPORT(INFO, "The width of the polynomial evaluator output is: "<< wR);
		/* number of digits to the left of the . */
		REPORT(INFO, "The weight of the PE output is: "<<weightR);


		addInput ("X", wInX);
		vhdl << tab << declare("addr", addrSize) << " <= X"<<range(wInX-1, wInX-addrSize)<<";"<<endl;

		/* Instantiate the table generator */
		nextCycle();/////////////////////////////////// The Coefficent ROM has a registered iunput
		inPortMap ( tg, "X", "addr");
		outPortMap ( tg, "Y", "Coef");
		vhdl << instance ( tg, "GeneratedTable" );
		syncCycleFromSignal("Coef");
		nextCycle();/////////////////////////////////// The Coefficent ROM has a registered output
		
		vhdl << tab << declare ("y",ySize) << " <= X"<<range(ySize-1 ,0) << ";" << endl;
		
		/* get the coefficients */
		int lsb = 0, sizeS = 0;
		for (uint32_t i=0; i< pe->getCoeffParamVector().size(); i++){
			lsb += sizeS;
			sizeS = pe->getCoeffParamVector()[i]->getSize()+1;
			vhdl << tab << declare(join("a",i), sizeS ) << "<= Coef"<< range (lsb+sizeS-1, lsb) << ";" << endl;
		}
		
		/* Instantiate the polynomial evaluator generator*/
		inPortMap( pe, "Y", "y");
		for (uint32_t i=0; i< pe->getCoeffParamVector().size(); i++){
			inPortMap(pe,  join("a",i), join("a",i));
		}
		
		outPortMap( pe, "R", "Rpe");
		vhdl << instance( pe, "PolynomialEvaluator");

		syncCycleFromSignal("Rpe");
		if (finalRounding_){
			/* number of bits to recover is */
			int recover = weightR + wOutX;
			vhdl << tab << declare("sticky",1) << " <= '0' when Rpe"<<range(pe->getRWidth()-(pe->getRWeight()+recover)-3,0) <<" = " 
		                                                        <<   zg(pe->getRWidth()-(pe->getRWeight()+recover)-2,0) << " else '1';"<<endl;
			vhdl << tab << declare("extentedf", 1 + recover + 2) << " <= Rpe"<<range(pe->getRWidth()-pe->getRWeight(), pe->getRWidth()-(pe->getRWeight()+recover)-2)<<";"<<endl; 
		                  
			nextCycle();                              
			IntAdder *a = new IntAdder(target, 1 + recover + 2);
			oplist.push_back(a);
	
			inPortMap(a, "X", "extentedf");
			inPortMapCst(a, "Y", zg(1 + recover + 2,0) );
			inPortMapCst(a, "Cin", "sticky");
			outPortMap(a, "R", "fPostRound");
			vhdl << instance(a, "Rounding_Adder");

			syncCycleFromSignal("fPostRound");
		
			addOutput("R", 1 + recover + 1);
			vhdl << tab << " R <= fPostRound"<<range(recover+2, 1)<<";"<<endl;
		}else{
			addOutput("R", pe->getRWidth());
			vhdl << tab << "R <= Rpe;" << endl;  
		}
		
	}

	FunctionEvaluator::~FunctionEvaluator() {
	}
	
	
	void FunctionEvaluator::emulate(TestCase* tc){
		setToolPrecision(165);
		int  nrFunctions=(pf->getPiecewiseFunctionArray()).size();
		
		Function *f = pf->getPiecewiseFunctionArray(0);

		if (nrFunctions !=1) {
			cerr<<"Warning: we are dealing with a piecewise function; only the first branch will be evaluated"<<endl;    
		}

		mpz_class rd, ru; // rounded down

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
		if (verbose){
		if (verbose==3){
		cout<<"Input is:"<<sPrintBinary(mpX)<<endl;
		cout<<"Output is:"<<sPrintBinary(mpR)<<endl;
	}
		}
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

		tc->addExpectedOutput("R", rd);
		tc->addExpectedOutput("R", ru);
		mpfr_clear(mpX);
		mpfr_clear(mpR);
	}

}
	
#endif //HAVE_SOLLYA	
