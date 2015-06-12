#include <iostream>
#include <sstream>
#include <sollya.h>


// include the header of the Operator
#include "FixComplexKCM.hpp"

using namespace std;
namespace flopoco {

	FixComplexKCM::FixComplexKCM(
			Target* target,
			bool signedInput,
			int msb_in, 
			int lsb_in, 
			int lsb_out,
			string constant_re,
			string constant_im
		): 	
			Operator(target),
			signedInput(signedInput),
			msb_in(msb_in),
			lsb_in(lsb_in),
			lsb_out(lsb_out),
			constant_re(constant_re),
			constant_im(constant_im)
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
		name << "FixComplexKCM_" << vhdlize(msb_in) <<"_" << vhdlize(lsb_in) 
			<< "_" << vhdlize(lsb_out) << "_" << vhdlize(constant_re) << "_" <<
			vhdlize(constant_im) << "_" << ((signedInput) ? "" : "un") <<
			"signed" ;

		setName(name.str());
		
		// Copyright 
		setCopyrightString("3IF 2015 dev team (2015)");

		int inputWidth = 1 + msb_in - lsb_in;
		
		// declaring inputs
		addInput ("ReIN" , inputWidth);
		addInput ("ImIN" , inputWidth);

		//Since it's a sum of two products
		msb_out =  2 * msb_in + 1;
		if(msb_out < lsb_out)
		{
			throw string(
				"FixComplexKCM: Error, the result computed "
				"would always be zero (msb_out < lsb_out)");
		}

		int outputWidth = msb_out - lsb_out + 1;
		if(!signedInput)
		{
			outputWidth ++;
		}

		// declaring output
		addOutput("ReOut", outputWidth);
		addOutput("ImOut", outputWidth);

		// basic message
		REPORT(INFO,"Declaration of FixComplexKCM\n");

		//Computing constants for testBench
		sollya_obj_t nodeIm, nodeRe;	
		nodeRe = sollya_lib_parse_string(constant_re.c_str());
		nodeIm = sollya_lib_parse_string(constant_im.c_str());
		
		if(sollya_lib_obj_is_error(nodeRe))
		{
			ostringstream error;
			error << srcFileName <<" : Unable to parse string \""  <<
				constant_re << "\" as a numeric constant" << endl;
		}
		
		if(sollya_lib_obj_is_error(nodeIm))
		{
			ostringstream error;
			error << srcFileName <<" : Unable to parse string \""  <<
				constant_im << "\" as a numeric constant" << endl;
		}

		mpfr_inits2(10000, mpfr_constant_re, mpfr_constant_im, NULL);

		sollya_lib_get_constant(mpfr_constant_re, nodeRe);
		sollya_lib_get_constant(mpfr_constant_im, nodeIm);

		//
		//VHDL outputs goes here
		//
	};

	
	void FixComplexKCM::emulate(TestCase * tc) {
		int inputWidth = msb_in - lsb_in + 1;		
		int outputWidth = msb_out - lsb_out + 1;
		if(!signedInput)
		{
			outputWidth++;
		}

		/* first we are going to format the entries */
		mpz_class reIn = tc->getInputValue("ReIN");
		mpz_class imIn = tc->getInputValue("ImIN");

		/* Sign handling */
		// Msb index counting from one
		bool reInNeg = (
				signedInput &&
				(mpz_tstbit(reIn.get_mpz_t(), inputWidth - 1) == 1)
			);

		bool imInNeg = (
				signedInput && 
				(mpz_tstbit(imIn.get_mpz_t(), inputWidth - 1) == 1)
			);

		// 2's complement -> absolute value unsigned representation
		if(reInNeg)
		{
			reIn = (mpz_class(1) << inputWidth) - reIn;
		}

		if(imInNeg)
		{
			imIn = (mpz_class(1) << inputWidth) - imIn;
		}

		//Cast to mp floating point number
		mpfr_t reIn_mpfr, imIn_mpfr;
		mpfr_init2(reIn_mpfr, inputWidth + 1);
		mpfr_init2(imIn_mpfr, inputWidth + 1);

		//Exact
		mpfr_set_z(reIn_mpfr, reIn.get_mpz_t(), GMP_RNDN); 
		mpfr_set_z(imIn_mpfr, imIn.get_mpz_t(), GMP_RNDN);

		//Exact
		mpfr_mul_2si(reIn_mpfr, reIn_mpfr, lsb_in, GMP_RNDN);
		mpfr_mul_2si(imIn_mpfr, imIn_mpfr, lsb_in, GMP_RNDN);

		mpfr_t re_prod, im_prod, crexim_prod, xrecim_prod;
		mpfr_t reOut, imOut;

		mpfr_inits2(
				2 * inputWidth + 1, 
				re_prod, 
				im_prod, 
				crexim_prod, 
				xrecim_prod, 
				NULL
			);

		mpfr_inits2(5 * outputWidth + 1, reOut, imOut, NULL);

		// c_r * x_r -> re_prod
		mpfr_mul(re_prod, reIn_mpfr, mpfr_constant_re, GMP_RNDN);

		// c_i * x_i -> im_prod
		mpfr_mul(im_prod, imIn_mpfr, mpfr_constant_im, GMP_RNDN);

		// c_r * x_i -> crexim_prod
		mpfr_mul(crexim_prod, mpfr_constant_re, imIn_mpfr, GMP_RNDN);

		// x_r * c_im -> xrecim_prod
		mpfr_mul(xrecim_prod, reIn_mpfr, mpfr_constant_im, GMP_RNDN);

		/* Input sign correction */
		if(reInNeg)
		{
			//Exact
			mpfr_neg(re_prod, re_prod, GMP_RNDN);
			mpfr_neg(xrecim_prod, xrecim_prod, GMP_RNDN);
		}

		if(imInNeg)
		{
			//Exact
			mpfr_neg(im_prod, im_prod, GMP_RNDN);
			mpfr_neg(crexim_prod, crexim_prod, GMP_RNDN);
		}

		mpfr_sub(reOut, re_prod, im_prod, GMP_RNDN);
		mpfr_add(imOut, crexim_prod, xrecim_prod, GMP_RNDN);

		bool reOutNeg = (mpfr_sgn(reOut) < 0);
		bool imOutNeg = (mpfr_sgn(imOut) < 0);

		if(reOutNeg)
		{
			//Exact
			mpfr_abs(reOut, reOut, GMP_RNDN);
		}

		if(imOutNeg)
		{
			//Exact
			mpfr_abs(imOut, imOut, GMP_RNDN);
		}

		//Scale back (Exact)
		mpfr_mul_2si(reOut, reOut, -lsb_out, GMP_RNDN);
		mpfr_mul_2si(imOut, imOut, -lsb_out, GMP_RNDN);

		//Get bits vector
		mpz_class reUp, reDown, imUp, imDown;

		mpfr_get_z(reUp.get_mpz_t(), reOut, GMP_RNDU);
		mpfr_get_z(reDown.get_mpz_t(), reOut, GMP_RNDD);
		mpfr_get_z(imDown.get_mpz_t(), imOut, GMP_RNDD);
		mpfr_get_z(imUp.get_mpz_t(), imOut, GMP_RNDU);

		cout << "OutputWidth : " << outputWidth << endl;
		cout << "re :" << reUp << endl << "im :" << imUp <<endl;

		//If result was negative, compute 2's complement
		if(reOutNeg)
		{
			reUp = (mpz_class(1) << outputWidth) - reUp;
			reDown = (mpz_class(1) << outputWidth) - reDown;
		}

		if(imOutNeg)
		{
			imUp = (mpz_class(1) << outputWidth) - imUp;
			imDown = (mpz_class(1) << outputWidth) - imDown;
		}

		//Add expected results to corresponding outputs
		tc->addExpectedOutput("ReOut", reUp);	
		tc->addExpectedOutput("ReOut", reDown);	
		tc->addExpectedOutput("ImOut", imUp);	
		tc->addExpectedOutput("ImOut", imDown);	
	}


	void FixComplexKCM::buildStandardTestCases(TestCaseList * tcl) {
		TestCase* tc;

		tc = new TestCase(this);		
		tc->addInput("ReIN", 0);
		tc->addInput("ImIN", 0);
		emulate(tc);
		tcl->add(tc);
		
		tc = new TestCase(this);		
		tc->addInput("ReIN", 1);
		tc->addInput("ImIN", 0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);		
		tc->addInput("ReIN", 0);
		tc->addInput("ImIN", 1);
		emulate(tc);
		tcl->add(tc);
		
		tc = new TestCase(this);		
		tc->addInput("ReIN", 1);
		tc->addInput("ImIN", 1);
		emulate(tc);
		tcl->add(tc);

		if(signedInput)
		{
			tc = new TestCase(this);		
			tc->addInput("ReIN", -1);
			tc->addInput("ImIN", -1);
			emulate(tc);
			tcl->add(tc);
		}

		tc = new TestCase(this);		
		tc->addInput("ReIN", 2);
		tc->addInput("ImIN", 0);
		emulate(tc);
		tcl->add(tc);

	}
}//namespace

//TODO : Free memory (mpfr_t)
