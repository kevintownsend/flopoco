#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixedPointSinOrCos.hpp"

using namespace std;

namespace flopoco{

	// TODO free scale in the destructor

	FixedPointSinOrCos::FixedPointSinOrCos(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
		int degree = 4;

		srcFileName="FixedPointSinOrCos";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "FixedPointSinOrCos_" << 1+w <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "FixedPointSinOrCos_" << 1+w << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		addInput("X",1+w,true);
		addInput("SelectFunction");		// 0 for cosine 1 for sine

		// declaring output
		addOutput("SorC",1+w,2);
		
		//reduce the argument X to [0, 1/2)
		vhdl << tab << declare("absX", 1+w ) << "<= (X xor (" << w << " downto 0 => X(" << w << ")))"
											<< " + " 
											<< "(" << zg(w, 0) << " & X(" << w << "));"<< endl;
		vhdl << tab << declare("reducedX", 1+w) 
						<< "<= (absX(" << w << " downto " << w-1 << ") - (\'0\' & (absX(" << w << ") xor absX(" << w-1 << "))))" 
						<< " & absX(" << w-2 << " downto 0);" << endl;
		vhdl << tab << declare("quadrantX", 2) << " <= X(" << w << " downto " << w-1 << ");" << endl;
		
		vhdl << tab << declare("XinSin", 1+w) << " <= reducedX;" << endl; 
		vhdl << tab << declare("XinCos", 1+w) << " <= (\"01\" & " << zg(w-1) << ") - reducedX;" << endl; 
		
		//remove the sign bit (msb)
		vhdl << tab << declare("shortXinSin", w) << " <= XinSin(" << w-1 << " downto 0);" << endl; 
		vhdl << tab << declare("shortXinCos", w) << " <= XinCos(" << w-1 << " downto 0);" << endl; 
		
		//compute the sine and the cosine
		FunctionEvaluator* sinCosEvaluator = new FunctionEvaluator(target, "sin(x*Pi),0,1,1", w, w, degree);
		oplist.push_back(sinCosEvaluator);
		
		inPortMap(sinCosEvaluator, "X", "shortXinSin");
		outPortMap(sinCosEvaluator, "R", "intSin");
		vhdl << instance(sinCosEvaluator, "sinEvaluator") << endl;
		
		inPortMap(sinCosEvaluator, "X", "shortXinCos");
		outPortMap(sinCosEvaluator, "R", "intCos");
		vhdl << instance(sinCosEvaluator, "cosEvaluator") << endl;
		
		//extract the needed bits from the function output
		vhdl << tab << declare("shortIntSin", 1+w) << " <= intSin(" << w << " downto 0);" << endl; 
		vhdl << tab << declare("shortIntCos", 1+w) << " <= intCos(" << w << " downto 0);" << endl; 
		
		//assign the correct value to the output
		vhdl << tab << declare("reducedC", 1+w) << "<= shortIntCos;" << endl;
		vhdl << tab << declare("reducedS", 1+w) << "<= shortIntSin;" << endl;
		
		vhdl << tab << declare("negReducedC", 1+w) << "<= (reducedC xor (" << w << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(w, 0) << " & \'1\');"<< endl;
		vhdl << tab << declare("negReducedS", 1+w) << "<= (reducedS xor (" << w << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(w, 0) << " & \'1\');"<< endl;
												   
		declare("effectiveC", 1+w);
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "effectiveC <= reducedC when \"00\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"10\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedC when others;" << endl;	//edit to signal error
		
		declare("effectiveS", 1+w);
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "effectiveS <= reducedS when \"00\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedC when \"10\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedS when others;" << endl;	//edit to signal error
		
		vhdl << tab << "SorC <= effectiveS when SelectFunction=\'1\' else effectiveC;" << endl;
		
	};

	void FixedPointSinOrCos::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpz_class svSelect = tc->getInputValue("SelectFunction");
		mpfr_t z, select, constPi, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		
		/* Compute correct value */
		mpfr_init2(z, 10*w);
		mpfr_init(select);
		
		mpfr_init2(constPi, 10*w);
		
		mpfr_set_z (select, svSelect.get_mpz_t(), GMP_RNDN); 
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, w, GMP_RNDN); // this rounding is acually exact
		
		mpfr_const_pi( constPi, GMP_RNDN);
		mpfr_mul(z, z, constPi, GMP_RNDN);
		
		mpfr_init2(rsin, 10*w); 
		mpfr_init2(rcos, 10*w); 


		mpz_init2 (rsin_z, 2*w);
		mpz_init2 (rcos_z, 2*w);
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);

		mpfr_add_d(rsin, rsin, 6.0, GMP_RNDN);
		mpfr_add_d(rcos, rcos, 6.0, GMP_RNDN);
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDN); // exact rnd here

		// Rounding down
		mpfr_get_z (rsin_z, rsin, GMP_RNDD); // there can be a real rounding here
		mpfr_get_z (rcos_z, rcos, GMP_RNDD); // there can be a real rounding here
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		sin_zc -= mpz_class(6)<<w;
		cos_zc -= mpz_class(6)<<w;
		
		if(mpfr_cmp_d (select, 0.0) == 0)
			tc->addExpectedOutput ("SorC", cos_zc);
		else
			tc->addExpectedOutput ("SorC", sin_zc);
		
		// Rounding up
		mpfr_get_z (rsin_z, rsin, GMP_RNDU); // there can be a real rounding here
		mpfr_get_z (rcos_z, rcos, GMP_RNDU); // there can be a real rounding here
		mpz_class sin_zcu (rsin_z), cos_zcu (rcos_z);
		sin_zcu -= mpz_class(6)<<w;
		cos_zcu -= mpz_class(6)<<w;

		if(mpfr_cmp_d (select, 0.0) == 0)
			tc->addExpectedOutput ("SorC", cos_zcu);
		else
			tc->addExpectedOutput ("SorC", sin_zcu);
		


		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
		mpfr_free_cache();
	}


	void FixedPointSinOrCos::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpf_t zinit;
		mpfr_t z;
		mpz_t z_z;
		
		//mpf_set_default_prec (1+wI+wF+guard);
		
		mpfr_init2(z, 1+w+ceil(log2(1 + w)));
		mpz_init2 (z_z, 1+w+ceil(log2(1 + w)));
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		tc -> addInput ("SelectFunction",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "1.5707963267949e0", 10);
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("SelectFunction",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "0.5235987755983e0", 10);
		mpf_set_str (zinit, "0.16666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("SelectFunction",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "0.78539816339745e0", 10);
		mpf_set_str (zinit, "0.25e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("SelectFunction",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "1.0471975511966e0", 10);
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD);
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("SelectFunction",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		
		mpfr_clears (z, NULL);
	}

}
