/*
  FixFunction object for FloPoCo

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  

  All rights reserved.

*/

#include "FixFunction.hpp"
#include <sstream>

namespace flopoco{


	FixFunction::FixFunction(string sollyaString_)
	{
 		ostringstream completeDescription;
		completeDescription << sollyaString_  ; // << " on [" << printMPFR(xmin) << ", " << printMPFR(xmax) << "]"; 
		description = completeDescription.str();
		srcFileName = "FixFunction";  // useful only for reporting
		uniqueName_=description;  // useful only for reporting

		// Now do the parsing in Sollya
		fS= sollya_lib_parse_string(sollyaString_.c_str());

		/* If  parse error throw an exception */
		if (sollya_lib_obj_is_error(fS))
			THROWERROR("Unable to parse input function.");

		REPORT(DEBUG, "Function: " << description << " successfully parsed");
	}



		
	FixFunction::~FixFunction()
	{
	  sollya_lib_clear_obj(fS);
	}

	string FixFunction::getDescription() const
	{
		return description;
	}

	void FixFunction::eval(mpfr_t r, mpfr_t x) const
	{
		sollya_lib_evaluate_function_at_point(r, fS, x, NULL);
	}

	double FixFunction::eval(double x) const
	{
		mpfr_t mpX, mpR;
		double r;

		mpfr_inits(mpX, mpR, NULL);
		mpfr_set_d(mpX, x, GMP_RNDN);
		sollya_lib_evaluate_function_at_point(mpR, fS, mpX, NULL);
		r = mpfr_get_d(mpR, GMP_RNDN);
		mpfr_clears(mpX, mpR, NULL);

		return r;
	}

	sollya_obj_t FixFunction::getSollyaObj() const
	{
		return fS;
	}

} //namespace
