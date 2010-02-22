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
		
		
		pe = new PolynomialEvaluator(target, tg->getCoeffParamVector(), y, wOutX, tg->getMaxApproxError() );
		
		oplist.push_back(pe);

		

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
		vhdl << tab << declare( "op", wOutX+1 + (pe->getOutputWeight())) << " <= Rpe"<<range(pe->getOutputSize()-1, pe->getOutputSize() - (wOutX+pe->getOutputWeight()+ 1)) << ";" << endl;
		
//		IntAdder *ia = new IntAdder(target, wOutX+1 + (pe->getOutputWeight()+1));
//		oplist.push_back(ia);
//		
//		inPortMap    (ia, "X", "op1_rnd");
//		inPortMapCst (ia, "Y", zg(wOutX+1 + (pe->getOutputWeight()+1), 0));
//		inPortMapCst (ia, "Cin", "'1'");
//		outPortMap   (ia, "R", "postRoundR");   

//		vhdl << instance( ia, "Final_Round");
//		syncCycleFromSignal( "postRoundR" );  
		
		addOutput("R", wOutX+1 + pe->getOutputWeight());
		vhdl << tab << "R <= op;" << endl;  
	}

	FunctionEvaluator::~FunctionEvaluator() {
	}
	
	void FunctionEvaluator::emulate(TestCase* tc)
	{
#if 0	
		mpz_class svX = tc->getInputValue("X");

		/* Get inputs / outputs */
		mpz_class &x  = a[0];
		mpz_class &rd = a[1]; // rounded down
		mpz_class &ru = a[2]; // rounded up

//		int outSign = 0;

//		mpfr_t mpX, mpR;
//		mpfr_inits(mpX, mpR, 0, NULL);

//		/* Convert a random signal to an mpfr_t in [0,1[ */
//		mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
//		mpfr_div_2si(mpX, mpX, wI, GMP_RNDN);

//		/* Compute the function */
//		f.eval(mpR, mpX);

//		/* Compute the signal value */
//		if (mpfr_signbit(mpR))
//			{
//				outSign = 1;
//				mpfr_abs(mpR, mpR, GMP_RNDN);
//			}
//		mpfr_mul_2si(mpR, mpR, wO, GMP_RNDN);

//		/* NOT A TYPO. HOTBM only guarantees faithful
//		 * rounding, so we will round down here,
//		 * add both the upper and lower neighbor.
//		 */
//		mpfr_get_z(rd.get_mpz_t(), mpR, GMP_RNDD);
//		ru = rd + 1;
#endif

	}

}
	
	
