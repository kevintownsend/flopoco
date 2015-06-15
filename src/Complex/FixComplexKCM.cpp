#include <iostream>
#include <sstream>
#include <sollya.h>


#include "ConstMult/FixRealKCM.hpp"
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

		input_width = 1 + msb_in - lsb_in;
		
		// declaring inputs
		addInput ("ReIN" , input_width);
		addInput ("ImIN" , input_width);

		//Computing constants for testBench and in order to know constant width
		sollya_obj_t nodeIm, nodeRe;	
		nodeRe = sollya_lib_parse_string(constant_re.c_str());
		
		if(sollya_lib_obj_is_error(nodeRe))
		{
			ostringstream error;
			error << srcFileName <<" : Unable to parse string \""  <<
				constant_re << "\" as a numeric constant" << endl;
			throw error.str();
		}
		
		nodeIm = sollya_lib_parse_string(constant_im.c_str());
		if(sollya_lib_obj_is_error(nodeIm))
		{
			ostringstream error;
			error << srcFileName <<" : Unable to parse string \""  <<
				constant_im << "\" as a numeric constant" << endl;
			throw error.str();
		}

		mpfr_inits2(10000, mpfr_constant_re, mpfr_constant_im, NULL);

		sollya_lib_get_constant(mpfr_constant_re, nodeRe);
		sollya_lib_get_constant(mpfr_constant_im, nodeIm);

		mpfr_t log2C;
		mpfr_init2(log2C, 100); 
		
		//Constant real part width
		mpfr_log2(log2C, mpfr_constant_re, GMP_RNDN);
		constantReMsb = mpfr_get_si(log2C, GMP_RNDU);

		//Constant imaginary part width
		mpfr_log2(log2C, mpfr_constant_im, GMP_RNDN);
		constantImMsb = mpfr_get_si(log2C, GMP_RNDU);

		//Free
		mpfr_clear(log2C);

		int constantMaxMsb = (constantImMsb > constantReMsb) ? 	constantImMsb :
																constantReMsb;

		//Since it's a sum of two products
		msb_out =  msb_in + constantMaxMsb + 1;

		if(msb_out < lsb_out)
		{
			throw string(
					"FixComplexKCM: Error, the result computed "
					"would always be zero (msb_out < lsb_out)"
					);
		}

		output_width = msb_out - lsb_out + 1;
		if(!signedInput)
		{
			output_width ++;
		}

		// declaring output
		addOutput("ReOut", output_width);
		addOutput("ImOut", output_width);

		// basic message
		REPORT(INFO,"Declaration of FixComplexKCM\n");
	};

	FixComplexKCM::~FixComplexKCM()
	{
		mpfr_clears(mpfr_constant_im, mpfr_constant_re, NULL);
	}
	
	void FixComplexKCM::emulate(TestCase * tc) {

		if(!signedInput)
		{
			output_width++;
		}

		/* first we are going to format the entries */
		mpz_class reIn = tc->getInputValue("ReIN");
		mpz_class imIn = tc->getInputValue("ImIN");

		/* Sign handling */
		// Msb index counting from one
		bool reInNeg = (
				signedInput &&
				(mpz_tstbit(reIn.get_mpz_t(), input_width - 1) == 1)
			);

		bool imInNeg = (
				signedInput && 
				(mpz_tstbit(imIn.get_mpz_t(), input_width - 1) == 1)
			);

		// 2's complement -> absolute value unsigned representation
		if(reInNeg)
		{
			reIn = (mpz_class(1) << input_width) - reIn;
		}

		if(imInNeg)
		{
			imIn = (mpz_class(1) << input_width) - imIn;
		}

		//Cast to mp floating point number
		mpfr_t reIn_mpfr, imIn_mpfr;
		mpfr_init2(reIn_mpfr, input_width + 1);
		mpfr_init2(imIn_mpfr, input_width + 1);

		//Exact
		mpfr_set_z(reIn_mpfr, reIn.get_mpz_t(), GMP_RNDN); 
		mpfr_set_z(imIn_mpfr, imIn.get_mpz_t(), GMP_RNDN);

		//Exact
		mpfr_mul_2si(reIn_mpfr, reIn_mpfr, lsb_in, GMP_RNDN);
		mpfr_mul_2si(imIn_mpfr, imIn_mpfr, lsb_in, GMP_RNDN);

		mpfr_t re_prod, im_prod, crexim_prod, xrecim_prod;
		mpfr_t reOut, imOut;

		mpfr_inits2(
				2 * input_width + 1, 
				re_prod, 
				im_prod, 
				crexim_prod, 
				xrecim_prod, 
				NULL
			);

		mpfr_inits2(5 * output_width + 1, reOut, imOut, NULL);

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
		mpz_class reUp, reDown, imUp, imDown, carry;

		mpfr_get_z(reUp.get_mpz_t(), reOut, GMP_RNDU);
		mpfr_get_z(reDown.get_mpz_t(), reOut, GMP_RNDD);
		mpfr_get_z(imDown.get_mpz_t(), imOut, GMP_RNDD);
		mpfr_get_z(imUp.get_mpz_t(), imOut, GMP_RNDU);
		carry = 0;

		//If result was negative, compute 2's complement
		if(reOutNeg)
		{
			reUp = (mpz_class(1) << output_width) - reUp;
			reDown = (mpz_class(1) << output_width) - reDown;
		}

		if(imOutNeg)
		{
			imUp = (mpz_class(1) << output_width) - imUp;
			imDown = (mpz_class(1) << output_width) - imDown;
		}

		//Handle border cases
		if(imUp > (mpz_class(1) << output_width) - 1 )
		{
			cout << "OhOhOhh" ;
			imUp = 0;
		}

		if(reUp > (mpz_class(1) << output_width) - 1)
		{
			reUp = 0;
		}

		if(imDown > (mpz_class(1) << output_width) - 1 )
		{
			cout << "OhOhOhh" ;
			imDown = 0;
		}

		if(reDown > (mpz_class(1) << output_width) - 1)
		{
			reDown = 0;
		}

		//Add expected results to corresponding outputs
		tc->addExpectedOutput("ReOut", reUp);	
		tc->addExpectedOutput("ReOut", reDown);	
		tc->addExpectedOutput("ImOut", imUp);	
		tc->addExpectedOutput("ImOut", imDown);	

		mpfr_clears(
				reOut,
				imOut, 
				re_prod, 
				im_prod, 
				crexim_prod, 
				xrecim_prod, 
				reIn_mpfr, 
				imIn_mpfr,
				NULL
			);
	}


	void FixComplexKCM::buildStandardTestCases(TestCaseList * tcl) {
		TestCase* tc;

		int one = 1;
		if(lsb_in < 0 && msb_in >= 0)
		{
			//Real one
			one = one << -lsb_in;
		}

		tc = new TestCase(this);		
		tc->addInput("ReIN", 0);
		tc->addInput("ImIN", 0);
		emulate(tc);
		tcl->add(tc);
		
		tc = new TestCase(this);		
		tc->addInput("ReIN", one);
		tc->addInput("ImIN", 0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);		
		tc->addInput("ReIN", 0);
		tc->addInput("ImIN", one);
		emulate(tc);
		tcl->add(tc);
		
		tc = new TestCase(this);		
		tc->addInput("ReIN", one);
		tc->addInput("ImIN", one);
		emulate(tc);
		tcl->add(tc);

		if(signedInput)
		{
			tc = new TestCase(this);		
			tc->addInput("ReIN", -1 * one);
			tc->addInput("ImIN", -1 * one);
			emulate(tc);
			tcl->add(tc);
		}

		tc = new TestCase(this);		
		tc->addInput("ReIN", 2 * one);
		tc->addInput("ImIN", 0);
		emulate(tc);
		tcl->add(tc);

	}
}//namespace
