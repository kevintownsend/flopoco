/*
  Floating Point Squarer for FloPoCo
  
  Author : Bogdan Pasca
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

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
#include "FPSquarer.hpp"
#include "FPNumber.hpp"
#include "IntSquarer.hpp"

using namespace std;

namespace flopoco{


	extern vector<Operator*> oplist;


	FPSquarer::FPSquarer(Target* target, int wE, int wFX, int wFR) :
		Operator(target), wE_(wE), wFX_(wFX), wFR_(wFR){

		ostringstream name;
		name<<"FPSquarer_"<<wE_<<"_"<<wFX_<<"_"<<wFR_; 
		setName(name.str());
		setCopyrightString("Bogdan Pasca (2009)");
	
		/* Set up the IO signals */
		/* Inputs: 2b(Exception) + 1b(Sign) + wE_ bits (Exponent) + wFX_ bits(Fraction) */
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addFPInput  ("X"   , wE_, wFX_);
		addFPOutput ("R"   , wE_, wFR_);  

		vhdl << tab << declare("exc",2) << " <= " << use("X") << range(wFX+wE+3-1,wFX+1+wE) << ";" << endl;
		vhdl << tab << declare("exp" ,wE_)   << " <= " << use("X") << range(wFX+wE-1,wFX) << ";" << endl;
		vhdl << tab << declare("frac",wFX_+1) << " <= " << "\"1\" & " << use("X") << range(wFX-1,0) << ";" << endl;
	
		//process the exponent 
		vhdl << tab << declare("extExponent", wE+2) << "<=" << "\"0\" & " << use("exp") << " & \"0\";"<<endl;
		vhdl << tab << declare("negBias",     wE+2) << "<=" << "CONV_STD_LOGIC_VECTOR("<< ((1<<(wE+2))-1)-((1<<(wE-1))-1)<<","<<wE+2<<");"<<endl;
	
		vhdl << tab << declare("extExpPostBiasSub",wE+2) << " <= " << use("extExponent") << " + " << use("negBias") << " + '1';"<<endl; 
	
		//instantiate an IntSquarer
		IntSquarer *sqr = new IntSquarer(target,wFX_+1);
		oplist.push_back(sqr);
	
		inPortMap(sqr, "X", use("frac"));
		outPortMap(sqr, "R", "sqrFrac");
		vhdl << tab << instance(sqr, "FractionSquarer");
		syncCycleFromSignal("sqrFrac");
		nextCycle();////
	
		//we have 2 cases
		if (wFR > 2*(wFX_+1)-1){
			//no rounding needed
			vhdl << tab << declare("extExp",wE+2) << " <= " << use("extExpPostBiasSub") << " + " << use("sqrFrac")<<"("<<2*(wFX_+1)-1<<");" << endl;
			vhdl << tab << declare("finalFrac",wFR) << "<= " << use("sqrFrac")<<range(2*(wFX_+1)-2,0) << " & " << zg(wFR - (2*(wFX_+1)-1),0) << " when " <<use("sqrFrac")<<"("<<2*(wFX_+1)-1<<")='1' else "  << endl;
			vhdl << tab << tab << use("sqrFrac")<<range(2*(wFX_+1)-3,0) << " & " << zg(wFR - (2*(wFX_+1)-2),0) << ";" << endl;
	    
			vhdl << tab << declare ("excConcat",4) << " <= " << use("exc") << " & " <<  use("extExp")<<range(wE+1,wE) << ";" <<endl;
			//exception bits
			vhdl << tab << "with " << use("excConcat") << " select " << endl;
			vhdl << tab << declare("excR",2) << "<=""\"00\" when \"0000\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0001\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0010\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0011\"," << endl;
			vhdl << tab << tab << "\"01\" when \"0100\"," << endl;
			vhdl << tab << tab << "\"10\" when \"0101\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0110\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0111\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1000\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1001\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1010\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1011\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1100\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1101\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1110\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1111\"," << endl;
			vhdl << tab << tab << "\"11\" when others;" << endl;
	    
			//compose result
			vhdl << tab << "R" << " <= " << use("excR") << " & " << " \"0\"  & " << use("extExp") <<range(wE-1, 0) << " & " << use("finalFrac") << ";" << endl;
		}else{
			//rounding will be needed
			vhdl << tab << declare("sticky",1) << "<=" << "'0' when "<<use("sqrFrac")<<range( 2*(wFX+1)-(wFR+3)-1,0)<<"="<<zg(2*(wFX+1)-(wFR+3),0) << "else '1';"<<endl;
			vhdl << tab << declare("guard",1) << " <= " << use("sqrFrac")<<of(2*(wFX+1)-(wFR+3))<<" when "<<use("sqrFrac")<<of(2*(wFX+1)-1)<<"='0' else " 
				  << use("sqrFrac")<<of(2*(wFX+1)-(wFR+3)+1) << ";"<< endl;
			vhdl << tab << declare("fracULP",1) << "<=" << use("sqrFrac")<<of(2*(wFX+1)-(wFR+3)+1)<<" when "<<use("sqrFrac")<<of(2*(wFX+1)-1)<<"='0' else " 
				  << use("sqrFrac")<<of(2*(wFX+1)-(wFR+3)+2) << ";"<< endl;

			vhdl << tab << declare("extExp",wE+2) << " <= " << use("extExpPostBiasSub") << " + " << use("sqrFrac")<<"("<<2*(wFX_+1)-1<<");" << endl; //the normalization
			//not really final
			vhdl << tab << declare("finalFrac",wFR) << "<= " << use("sqrFrac")<<range(2*(wFX_+1)-2,2*(wFX_+1)-1-wFR) << " when " <<use("sqrFrac")<<"("<<2*(wFX_+1)-1<<")='1' else "  << endl; 
			vhdl << tab << tab << use("sqrFrac")<<range(2*(wFX_+1)-3,2*(wFX_+1)-3+1-wFR) << ";" << endl;
	    
			//the rounding phase
			IntAdder* add = new IntAdder(target, wE+2 + wFR);
			oplist.push_back(add); 
		
			vhdl << tab << declare("concatExpFrac", wE+2 + wFR) << " <= " << use("extExp") << " & " << use("finalFrac") << ";" << endl;
			vhdl << tab << declare("addCin",1) << " <= " << "(" << use("guard") << " and " << use("sticky") << ") or (" <<use("fracULP") << " and " << use("guard") << " and " << "not(" << use("sticky")<<"));"<<endl;
	
			inPortMap(add,"X", use("concatExpFrac"));
			inPortMapCst(add,"Y",zg(wE+2 + wFR,0));
			inPortMap(add, "Cin", use("addCin"));
			outPortMap(add, "R", "postRound");
			vhdl << instance (add, "Rounding_Instance");

			syncCycleFromSignal("postRound");	


		
			vhdl << tab << declare ("excConcat",4) << " <= " << use("exc") << " & " <<  use("postRound")<<range(wE+2 + 
	    
	    
																																				 wFR-1,wE+2 + wFR-2) << ";" <<endl;
			//exception bits
			vhdl << tab << "with " << use("excConcat") << " select " << endl;
			vhdl << tab << declare("excR",2) << "<=""\"00\" when \"0000\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0001\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0010\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0011\"," << endl;
			vhdl << tab << tab << "\"01\" when \"0100\"," << endl;
			vhdl << tab << tab << "\"10\" when \"0101\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0110\"," << endl;
			vhdl << tab << tab << "\"00\" when \"0111\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1000\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1001\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1010\"," << endl;
			vhdl << tab << tab << "\"10\" when \"1011\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1100\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1101\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1110\"," << endl;
			vhdl << tab << tab << "\"11\" when \"1111\"," << endl;
			vhdl << tab << tab << "\"11\" when others;" << endl;

			vhdl << tab << "R" << " <= " << use("excR") << " & " << " \"0\"  & " << use("postRound") <<range(wE + wFR -1, wFR) << " & " << use("postRound")<<range(wFR-1,0) << ";" << endl;
		
		}
	
	}

	FPSquarer::~FPSquarer() {
	}


	// FIXME: the following is only functional for a correctly rounded multiplier.
	// The old code for the non-normalized case is commented out below, just in case.
	void FPSquarer::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");

		/* Compute correct value */
		FPNumber fpx(wE_, wFX_);
		fpx = svX;
		mpfr_t x, r;
		mpfr_init2(x, 1+wFX_);
		mpfr_init2(r, 1+wFR_); 
		fpx.getMPFR(x);
		mpfr_mul(r, x, x, GMP_RNDN);

		// Set outputs 
		FPNumber  fpr(wE_, wFR_, r);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		// clean up
		mpfr_clears(x,r, NULL);
	}


}
