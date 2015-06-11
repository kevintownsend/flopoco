// general c++ library for manipulating streams
#include <iostream>
#include <sstream>


// include the header of the Operator
#include "FixComplexKCM.hpp"

using namespace std;
namespace flopoco {

	FixComplexKCM::FixComplexKCM(
			Target* target,
			int msb_in, 
			int lsb_in, 
			int lsb_out,
			string constant
		): 	
			Operator(target),
			signedInput(signedInput),
			msb_in(msb_in),
			lsb_in(lsb_in),
			lsb_out(lsb_out),
			constant(constant)
	{
		if(lsb_in>msb_in) 
		{
			throw string("FixComplexKCM: Error, lsbIn>msbIn");
		}

		// definition of the source file name, used for info and error reporting
		// using REPORT 
		srcFileName="FixComplexKCM";

		// definition of the name of the operator
		ostringstream name;
		name << "FixComplexKCM_" << msb_in <<"_" << lsb_in << "_" << 
			lsb_out << "_" << constant ;

		setName(name.str());
		
		// Copyright 
		setCopyrightString("3IF 2015 dev team (2015)");

		int inputWidth = 1 + msb_in - lsb_in;
		
		// declaring inputs
		addInput ("ReIN" , inputWidth);
		addInput ("ImIN" , inputWidth);

		//Since it's a sum of two products
		msb_out =  * msb_in + 1;
		if(msb_out < lsb_out)
		{
			throw string(
				"FixComplexKCM: Error, the result computed "
				"would always be zero (msb_out < lsb_out)");
		}

		int outputWidth = msb_out - lsb_out + 1;

		// declaring output
		addOutput("ReOut", outputWidth);
		addOutput("ImOut", outputWidth);

		// basic message
		REPORT(INFO,"Declaration of FixComplexKCM\n");

		//
		//VHDL outputs goes here
		//
		
	};

	
	void FixComplexKCM::emulate(TestCase * tc) {
		int inputWidth = msb_in - lsb_in + 1;		
		int outputWidth = msb_out - lsb_out + 1;

		/* first we are going to format the entries */
		mpz_class reIn = tc->getInputValue("ReIN");
		mpz_class imIn = tc->getInputValue("ImIN");

		//Cast to floating point number
		mpfr_t reIn_mpfr, imIn_mpfr;
		mpfr_init2(reIn_mpfr, inputWidth + 1);
		mpfr_init2(imIn_mpfr, inputWidth + 1);

		//Exact
		mpft_set_z(reIn_mpfr, reIn.get_mpz_t(), GMP_RNDN); 
		mpft_set_z(imIn_mpfr, imIn.get_mpz_t(), GMP_RNDN);

		//Exact
		mpfr_mul_2si(reIn_mpfr, reIn_mpfr, lsb_in, GMP_RNDN);
		mpfr_mul_2si(imIn_mpfr, imIn_mpfr, lsb_in, GMP_RNDN);

		mpfr_t re_prod, im_prod, crexim_prod, xrecim_prod;
		mpfr_t reOut, imOut;

		mpfr_inits2(inputWidth + 1, re_prod, reIn_mpfr, crexim_prod, xrecim_prod);
		mpfr_inits2(outputWidth + 1, reOut, imOut);

		//Products : Error <= 1/4 * 2^lsb_in 
		mpfr_mult(re_prod, reIn_mpfr, constant_re, GMP_RNDN);
		mpfr_mult(im_prod, imIn_mpfr, constant_im, GMP_RNDN);
		mpfr_mult(crexim_prod, constant_re, imIn_mpfr, GMP_RNDN);
		mpfr_mult(xrecim_prod, reIn_mpfr, constant_im, GMP_RNDN);

		//Sums : Error <= 2^(lsb_in -1)	
		mpfr_sub(reOut, re_prod, im_prod, GMP_RNDN);
		mpfr_add(imOut, crexim_prod, xrecim_prod, GMP_RNDN);

		//Scale back (Exact)
		mpfr_mul_2si(reOut, reOut, -lsb_out);
		mpfr_mul_2si(imOut, imOut, -lsb_out);

		//Get bits vector
		mpz_class reUp, reDown, imUp, imDown;

		mprfr_get_z(reUp.get_mpz_t(), reOut, GMP_RNDU);
		mprfr_get_z(reDown.get_mpz_t(), reOut, GMP_RNDD);
		mprfr_get_z(imDown.get_mpz_t(), imOut, GMP_RNDD);
		mprfr_get_z(imUp.get_mpz_t(), imOut, GMP_RNDU);
		
		//Add expected results to corresponding outputs
		tc->addExpectedOutput("ReOut", reUp);	
		tc->addExpectedOutput("ReOut", reDown);	
		tc->addExpectedOutput("ImOut", imUp);	
		tc->addExpectedOutput("ImOut", imDown);	
	}


	void FixComplexKCM::buildStandardTestCases(TestCaseList * tcl) {
		// please fill me with regression tests or corner case tests!
	}
}//namespace
