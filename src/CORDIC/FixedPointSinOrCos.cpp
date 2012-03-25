#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixedPointSinOrCos.hpp"

using namespace std;

namespace flopoco{

	//choice of sine/cosine is made based on the msb of the input
	FixedPointSinOrCos::FixedPointSinOrCos(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
		int degree = 4;
		
		//initialize testing random number
		gmp_randinit_mt (state);
		
		srcFileName="FixedPointSinOrCos";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "FixedPointSinOrCos_" << 1+w << "_f" << target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "FixedPointSinOrCos_" << 1+w << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		// one extra bit in the beginning to switch between sine and cosine
		// 0 = cosine, 1 = sine
		addInput  ( "X"  , 1+w, true );
		addInput  ( "Control"		 );

		// declaring output
		addOutput  ( "C"  , 1+w, 2 );
		addOutput  ( "S"  , 1+w, 2 );
		
		//prepare inputs 
		vhdl << tab << declare("Xin", 1+w) << "<= X;" << endl;
		vhdl << tab << declare("selectFunction") << "<= Control;" << endl;
		
		//get absolute value of the input
		vhdl << tab << declare("preNewXin", 1+w) << "<= Xin xor (" << w << " downto 0 => Xin(" << w << "));" << endl;
		vhdl << tab << declare("cInNewXin") << "<= X(" << w << ");" << endl;
		
		IntAdder *posAdder = new IntAdder(target, 1+w, inDelayMap("X",getCriticalPath()));
		oplist.push_back(posAdder);
		
		inPortMap(posAdder, "X", "preNewXin");
		inPortMapCst(posAdder, "Y", zg(1+w, 0));
		inPortMap(posAdder, "Cin", "cInNewXin");
		outPortMap (posAdder, "R", "newXin");
		vhdl << instance(posAdder, "posAdder") << endl;
		
		vhdl << tab << declare("newXinshort", w) << "<= X(" << w-1 << " downto 0);" << endl;
		
		//create the function evaluators
		FunctionEvaluator* sinEvaluator = new FunctionEvaluator(target, "sin(x*Pi/2),0,1,1", w, w, degree);
		oplist.push_back(sinEvaluator);
		
		inPortMap(sinEvaluator, "X", "newXinshort");
		outPortMap (sinEvaluator, "R", "intS");
		vhdl << instance(sinEvaluator, "sinEvaluator") << endl;
		
		FunctionEvaluator* cosEvaluator = new FunctionEvaluator(target, "cos(x*Pi/2),0,1,1", w, w, degree);
		oplist.push_back(cosEvaluator);
		
		inPortMap(cosEvaluator, "X", "newXinshort");
		outPortMap (cosEvaluator, "R", "intC");
		vhdl << instance(cosEvaluator, "cosEvaluator") << endl;
		
		//get the true value for sin and cos
		vhdl << tab << declare("preSMask", 1+w) << "<= (" << w << " downto 0 => Xin(" << w << "));" << endl;
		vhdl << tab << declare("preS", 1+w) << "<= intS(intS'length-2 downto intS'length-" << w << "-2) xor preSMask;" << endl;
		vhdl << tab << declare("cInS") << "<= X(" << w << ");" << endl;
		
		inPortMap(posAdder, "X", "preS");
		inPortMapCst(posAdder, "Y", zg(1+w, 0));
		inPortMap(posAdder, "Cin", "cInS");
		outPortMap (posAdder, "R", "correctSignS");
		vhdl << instance(posAdder, "sinAdder") << endl;
		
		vhdl << tab << "C <= intC(intC'length-3 downto intC'length-" << w << "-3) and (" << w << " downto 0 => (not selectFunction));" << endl;
		vhdl << tab << "S <= correctSignS and (" << w << " downto 0 => selectFunction);" << endl;
	};


	void FixedPointSinOrCos::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpz_class svControl = tc->getInputValue("Control");
		mpfr_t z, control, constPi, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		
		/* Compute correct value */
		mpfr_init2(z, 1+w+2);
		mpfr_init(control);
		
		mpfr_init2(constPi, 10*(1+w));
		
		mpfr_init2(rsin, 10*(1+w)); 
		mpfr_init2(rcos, 10*(1+w)); 
		mpz_init2 (rsin_z, 10*(1+w));
		mpz_init2 (rcos_z, 10*(1+w));
		
		mpfr_set_z (control, svControl.get_mpz_t(), GMP_RNDN); // this rounding is exact
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, w, GMP_RNDN); // this rounding is acually exact
		
		mpfr_mul_2si(z, z, -1, GMP_RNDN);
		mpfr_const_pi( constPi, GMP_RNDN);
		mpfr_mul(z, z, constPi, GMP_RNDN);
		
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDN); // exact rnd here
		mpfr_get_z (rsin_z, rsin, GMP_RNDN); // there can be a real rounding here
		
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDN); // exact rnd here
		mpfr_get_z (rcos_z, rcos, GMP_RNDN); // there can be a real rounding here

		// Set outputs 
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		if(mpfr_cmp_si(control, 0) == 0)
			tc->addExpectedOutput ("C", cos_zc);
		else
			tc->addExpectedOutput ("C", mpz_class(0));
		if(mpfr_cmp_si(control, 0) == 0)
			tc->addExpectedOutput ("S", mpz_class(0));
		else
			tc->addExpectedOutput ("S", sin_zc);

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
		
		//mpf_set_default_prec (1+wI+wF+guardxy);
		
		mpfr_init2(z, 1+w);
		mpz_init2 (z_z, 1+w);
		
		//z=0 ---------------------------
		//z=0 cos
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		tc -> addInput ("Control",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=0 sin
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		tc -> addInput ("Control",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		//z=0 ---------------------------
		
		//z=pi/2 ---------------------------
		//z=pi/2 cos
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "1.5707963267949e0", 10);
		mpf_set_str (zinit, "1.0e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2 sin
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "1.5707963267949e0", 10);
		mpf_set_str (zinit, "1.0e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		//z=pi/2 ---------------------------
		
		//z=pi/6 ---------------------------
		//z=pi/6 cos
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "0.5235987755983e0", 10);
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6 sin
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "0.5235987755983e0", 10);
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		//z=pi/6 ---------------------------
		
		//z=pi/4 ---------------------------
		//z=pi/4 cos
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "0.78539816339745e0", 10);
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4 sin
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "0.78539816339745e0", 10);
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		//z=pi/4 ---------------------------
		
		//z=pi/3 ---------------------------
		//z=pi/3 cos
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "1.0471975511966e0", 10);
		mpf_set_str (zinit, "0.66666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN);
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3 sin
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w);
		//mpf_set_str (zinit, "1.0471975511966e0", 10);
		mpf_set_str (zinit, "0.66666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDN);
		
		mpfr_mul_2si (z, z, w, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		
		tc -> addInput ("X",mpz_class(z_z));
		tc -> addInput ("Control",mpz_class(1));
		emulate(tc);
		tcl->add(tc);
		//z=pi/3 ---------------------------
		
		mpfr_clears (z, NULL);
	}

	//still testing
	TestCase* FixedPointSinOrCos::buildRandomTestCase(int i) 
	{
		TestCase* tc = new TestCase(this);
		mpz_class h;
		mpfr_t randomNumber;
		
		mpfr_init2 (randomNumber, 10*(1+w));
		mpfr_urandomb (randomNumber, state);
		mpfr_mul_2si(randomNumber, randomNumber, w, GMP_RNDN);
		mpfr_get_z(h.get_mpz_t(), randomNumber,  GMP_RNDN);
		
		tc = new TestCase (this);
		tc -> addInput ("X", h);
		tc -> addInput ("Control", mpz_class(0));
		emulate(tc);
		
		return tc;
	}

	std::string FixedPointSinOrCos::generateFixPointNumber(float x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		float xcopy = x;
		mpfr_t mx;
		mpz_class h;
		
		mpfr_init2 (mx, 1+wI+wF);
		
		if(xcopy<0){
			xcopy = xcopy * (-1);
		}
		
		mpfr_set_d(mx, xcopy, GMP_RNDD);
		mpfr_mul_2si(mx, mx, wF, GMP_RNDD);
		
		mpfr_get_z(h.get_mpz_t(), mx,  GMP_RNDD); 
        
        result = unsignedBinary(h, size);
        
        return result;
	}
	
	std::string FixedPointSinOrCos::generateFixPointNumber(mpf_t x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		mpfr_t mx;
		mpz_class h;
		
		mpfr_init2 (mx, 1+wI+wF);
		
		if(x<0){
			mpf_neg (x, x);
		}
		
		mpfr_set_f(mx, x, GMP_RNDD);
		mpfr_mul_2si(mx, mx, wF, GMP_RNDD);
		
		mpfr_get_z(h.get_mpz_t(), mx,  GMP_RNDD);         
        result = unsignedBinary(h, size);
        
        return result;
	}
	
	std::string FixedPointSinOrCos::getParamName(std::string s, int number)
	{
		std::stringstream aux;
		
		aux << s << number;
		
		return aux.str();
		
	}
	
	mpz_class FixedPointSinOrCos::fp2fix(mpfr_t x, int wI, int wF){
		mpz_class h;
		
		mpfr_mul_2si(x, x, wF, GMP_RNDD);
        mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDD);  
		
		return h;
	}

}











