#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCos.hpp"

using namespace std;

namespace flopoco{

	CordicSinCos::CordicSinCos(Target* target, int w_, int reducedIterations_, map<string, double> inputDelays) 
		: Operator(target), w(w_), reducedIterations(reducedIterations_)
	{
		srcFileName="CordicSinCos";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCos_" << 1+w <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCos_" << 1+w << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		addInput  ( "X"  , wcs, true );

		// declaring output
		addOutput  ( "C"  , wcs, 2 );
		addOutput  ( "S"  , wcs, 2 );
		
		Operator* CordicInstantiation;

		if(reducedIterations == 0)
			CordicInstantiation = new CordicSinCosClassic(target, w, inputDelays);
		else
			CordicInstantiation = new CordicSinCosRedIter(target, w, inputDelays);

		//setIndirectOperator(CordicInstantiation);
		cloneOperator(CordicInstantiation);
		oplist.push_back(CordicInstantiation);
	};


	void CordicSinCos::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpfr_t z, constPi, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		int g = ceil(log2(1 + w));
		
		/* Compute correct value */
		mpfr_init2(z, 1+w+g);
		
		mpfr_init2(constPi, 1+w+g);
		
		mpfr_init2(rsin, 1+w+g); 
		mpfr_init2(rcos, 1+w+g); 
		mpz_init2 (rsin_z, 1+w+g);
		mpz_init2 (rcos_z, 1+w+g);
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDD); // this rounding is exact
		mpfr_div_2si (z, z, w, GMP_RNDD); // this rounding is acually exact
		
		mpfr_const_pi( constPi, GMP_RNDD);
		mpfr_mul(z, z, constPi, GMP_RNDD);
		
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDD); // exact rnd here
		mpfr_get_z (rsin_z, rsin, GMP_RNDN); // there can be a real rounding here
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDD); // exact rnd here
		mpfr_get_z (rcos_z, rcos, GMP_RNDN); // there can be a real rounding here

		// Set outputs 
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		tc->addExpectedOutput ("C", cos_zc);
		tc->addExpectedOutput ("S", sin_zc);

		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
		mpfr_free_cache();
	}


	void CordicSinCos::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpf_t zinit;
		mpfr_t z;
		mpz_t z_z;
		
		mpfr_init2(z, 1+w+ceil(log2(1 + w)));
		mpz_init2 (z_z, 1+w+ceil(log2(1 + w)));
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		mpf_set_str (zinit, "0.16666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		mpf_set_str (zinit, "0.25e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD);
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		mpfr_clears (z, NULL);
	}
	
	void CordiSinCos::changeName(std::string operatorName){
	        Operator::changeName(operatorName);
	        if(getIndirectOperator())  getIndirectOperator()->changeName(operatorName);
        }

}











