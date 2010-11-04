/*
   A multiplier by a floating-point constant for FloPoCo

  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License,  2008-2010.
*/


// TODO Case mantissa=1, wE_out>wE_in
// TODO standard test vectors: 1, 0, various exn, xcut borders


#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "FPConstMult.hpp"
#include "../FPNumber.hpp"

using namespace std;


namespace flopoco{

extern vector<Operator*> oplist;


	FPConstMult::FPConstMult(Target* target, int wE_in_, int wF_in_, int wE_out_, int wF_out_, int cst_sgn_, int cst_exp_, mpz_class cst_sig_):
		Operator(target), 
		wE_in(wE_in_), wF_in(wF_in_), wE_out(wE_out_), wF_out(wF_out_), 
		cst_sgn(cst_sgn_), cst_exp_when_mantissa_int(cst_exp_), cst_sig(cst_sig_)
	{
		srcFileName="FPConstMult";
		ostringstream name;
		name <<"FPConstMult_"<<(cst_sgn==0?"":"M") <<mpz2string(cst_sig)<<"b"<<(cst_exp_when_mantissa_int<0?"M":"")<<abs(cst_exp_when_mantissa_int)<<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
		uniqueName_=name.str();

		if(cst_sig==0) {
			REPORT(INFO, "building a multiplier by 0, it will be easy");
			vhdl  << tab << "r <= " << rangeAssign(wE_out+wF_out+1, 0, "'0'") << ";"<<endl;
		}



		else {
			// Constant normalization
			while ((cst_sig % 2) ==0) {
				REPORT(INFO, "Significand is even, normalising");
				cst_sig = cst_sig >>1;
				cst_exp_when_mantissa_int+=1;
			}
			mantissa_is_one = false;
			if(cst_sig==1) {
				REPORT(INFO, "Constant mantissa is 1, multiplying by it will be easy"); 
				mantissa_is_one = true;
			}

			int cst_width = intlog2(cst_sig);
			cst_exp_when_mantissa_1_2 = cst_exp_when_mantissa_int + cst_width - 1; 

			// initialize mpfr constant
			mpfr_init2(mpfr_cst_sig, max(cst_width, 2));
			mpfr_set_z(mpfr_cst_sig, cst_sig.get_mpz_t(), GMP_RNDN); // exact op
			mpfr_mul_2si(mpfr_cst_sig, mpfr_cst_sig, -(cst_width-1), GMP_RNDN);  // exact op, gets mpfr_cst_sig in 1..2
			
			mpfr_init(mpfr_cst);
			mpfr_set(mpfr_cst, mpfr_cst_sig, GMP_RNDN);
			mpfr_mul_2si(mpfr_cst, mpfr_cst, cst_exp_when_mantissa_1_2, GMP_RNDN);
			
			if(cst_sgn==1)
				mpfr_neg(mpfr_cst,  mpfr_cst, GMP_RNDN);
			REPORT(INFO, "mpfr_cst = " << mpfr_get_d(mpfr_cst, GMP_RNDN));


			if(!mantissa_is_one) {			
				// initialize mpfr_xcut_sig = 2/cst_sig, will be between 1 and 2
				mpfr_init2(mpfr_xcut_sig, 4*(cst_width+wE_in+wE_out));
				mpfr_set_d(mpfr_xcut_sig, 2.0, GMP_RNDN);               // exaxt op
				mpfr_div(mpfr_xcut_sig, mpfr_xcut_sig, mpfr_cst_sig, GMP_RNDD);

				// now  round it up to wF_in+1 bits 
				mpfr_t xcut_wF;
				mpfr_init2(xcut_wF, wF_in+1);
				mpfr_set(xcut_wF, mpfr_xcut_sig, GMP_RNDD);
				mpfr_mul_2si(xcut_wF, xcut_wF, wF_in, GMP_RNDN);
				// It should now be an int; cast it into a mpz, then a mpz_class 
				mpz_t zz;
				mpz_init2(zz, wF_in+1);
				mpfr_get_z(zz, xcut_wF, GMP_RNDN);
				xcut_sig_rd= mpz_class(zz);

				REPORT(DETAILED, "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN));
				REPORT(DETAILED, "mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN)); 
				if(verbose) { // can't use REPORT because of GMP's C++ wrapper problem under windoze
					cerr << "xcut_sig_rd   = " << mpz2string(xcut_sig_rd) << "   ";
					printBinNumGMP(cerr, xcut_sig_rd, wF_in+1);  cerr << endl;
				}
				mpfr_clear(xcut_wF);
				mpz_clear(zz);
				// sub component
				icm = new IntConstMult(target, wF_in+1, cst_sig);
				oplist.push_back(icm);
			}


			// do all the declarations. Pushed into a method so that CRFPConstMult can inherit it
			setup();

		}

	}



	void FPConstMult::setup() {


		// Set up the IO signals
		addFPInput("X", wE_in, wF_in);
		addFPOutput("R", wE_out, wF_out);


		setCopyrightString("Florent de Dinechin (2007)");

		// bit width of constant exponent
		int wE_cst=intlog2(abs(cst_exp_when_mantissa_1_2));
		REPORT(DETAILED, "wE_cst = "<<wE_cst<<endl);
	
		// We have to compute Er = E_X - bias(wE_in) + E_C + bias(wE_R)
		// Let us pack all the constants together
		mpz_class expAddend = -bias(wE_in) + cst_exp_when_mantissa_1_2  + bias(wE_out);
		int expAddendSign=0;
		if (expAddend < 0) {
			expAddend = -expAddend;
			expAddendSign=1;
		}
		int wE_sum; // will be the max size of all the considered  exponents
		wE_sum = intlog2(expAddend);
		if(wE_in > wE_sum)
			wE_sum = wE_in;
		if(wE_out > wE_sum) 
			wE_sum = wE_out;


		vhdl << tab << declare("x_exn",2) << " <=  X("<<wE_in<<"+"<<wF_in<<"+2 downto "<<wE_in<<"+"<<wF_in<<"+1);"<<endl;
		vhdl << tab << declare("x_sgn") << " <=  X("<<wE_in<<"+"<<wF_in<<");"<<endl;
		vhdl << tab << declare("x_exp", wE_in) << " <=  X("<<wE_in<<"+"<<wF_in<<"-1 downto "<<wF_in<<");"<<endl;


		if(mantissa_is_one) {			
			vhdl << tab << "-- The mantissa of the constant is  1" << endl;
			if(wF_out == wF_in) {
				vhdl << tab << declare("r_frac", wF_out) << " <= X("<<wF_in-1 <<" downto 0);"<<endl;
			}
			else if(wF_out > wF_in){
				vhdl << tab << tab << declare("r_frac", wF_out) << " <= X("<<wF_in-1 <<" downto 0)  &  " << rangeAssign(wF_out-wF_in-1, 0, "'0'") << ";"<<endl;
				}
			else{ // wF_out < wF_in, this is a rounding of the mantissa TODO
				throw string("FPConstMult: multiplication by a power of two when  wF_out < wF_in not yet implemented, please complain to the FloPoCo team if you need it");
			}
			vhdl << tab << declare("gt_than_xcut") << " <= '0';"<<endl;
	
		}
		else{ // normal case, mantissa is not one
			vhdl << tab << declare("x_sig", wF_in+1) << " <= '1' & X("<<wF_in-1 <<" downto 0);"<<endl;
			vhdl << tab << declare("xcut_rd", wF_in+1) << " <= \""
				  << unsignedBinary(xcut_sig_rd, wF_in+1) << "\";"<<endl;
			vhdl << tab << declare("gt_than_xcut") << " <= '1' when ( x_sig("<<wF_in-1<<" downto 0) > xcut_rd("<<wF_in-1<<" downto 0) ) else '0';"<<endl;
			inPortMap  (icm, "inX", "x_sig");
			outPortMap (icm, "R","sig_prod");
			vhdl << instance(icm, "sig_mult");
			setCycleFromSignal("sig_prod"); 
			nextCycle();
			// Possibly shift the significand one bit left, and remove implicit 1 
			vhdl << tab << declare("shifted_frac",    wF_out+1) << " <= sig_prod("<<icm->rsize -2<<" downto "<<icm->rsize - wF_out-2 <<")  when gt_than_xcut = '1'"<<endl
				  << tab << "           else sig_prod("<<icm->rsize -3<<" downto "<<icm->rsize - wF_out - 3<<");"<<endl;  
			// add the rounding bit
			vhdl << tab << tab << declare("rounded_frac",   wF_out+1) << " <= (("<<wF_out <<" downto 1 => '0') & '1') + shifted_frac;"<<endl;
			vhdl << tab << tab << declare("r_frac", wF_out) << " <= rounded_frac("<<wF_out <<" downto  1);"<<endl;
		}

		// Handling signs is trivial
		if(cst_sgn==0)
			vhdl << tab << declare("r_sgn") << " <= x_sgn; -- positive constant"<<endl;
		else
			vhdl << tab << declare("r_sgn") << " <= not x_sgn; -- negative constant"<<endl;

		// exponent handling
		vhdl << tab << declare("abs_unbiased_cst_exp",wE_sum+1) << " <= \""
			  << unsignedBinary(expAddend, wE_sum+1) << "\";" << endl;
		vhdl << tab << declare("r_exp_nopb",    wE_out+1) << " <= "
			  << "((" << wE_sum << " downto " << wE_in << " => '0')  & x_exp)  "
			  << (expAddendSign==0 ? "+" : "-" ) << "  abs_unbiased_cst_exp"
			  << "  +  (("<<wE_sum<<" downto 1 => '0') & gt_than_xcut);"<<endl;

		// overflow handling
		vhdl << tab << declare("overflow") << " <= " ;
		if (maxExp(wE_in) + cst_exp_when_mantissa_1_2 + 1 < maxExp(wE_out)) // no overflow can ever happen
			vhdl << "'0'; --  overflow never happens for this constant and these (wE_in, wE_out)" << endl;
		else 
			vhdl <<  "'0' when r_exp_nopb(" << wE_sum << " downto " << wE_out << ") = (" << wE_sum << " downto " << wE_out <<" => '0')     else '1';" << endl;

		// underflow handling
		vhdl << tab << declare("underflow") << " <= " ;
		if (minExp(wE_in) + cst_exp_when_mantissa_1_2 > minExp(wE_out)) // no underflow can ever happen
			vhdl << "'0'; --  underflow never happens for this constant and these (wE_in, wE_out)" << endl;
		else 
			vhdl <<  "r_exp_nopb(" << wE_sum << ");" << endl;
			 
	
		vhdl << tab << declare("r_exp", wE_out) << " <= r_exp_nopb("<<wE_out-1<<" downto 0) ;"<<endl;

		vhdl << tab << declare("r_exn", 2) << " <=      \"00\" when ((x_exn = \"00\") or (x_exn = \"01\" and underflow='1'))  -- zero"<<endl 
			  << tab << "         else \"10\" when ((x_exn = \"10\") or (x_exn = \"01\" and overflow='1'))   -- infinity" << endl
			  << tab << "         else \"11\" when  (x_exn = \"11\")                      -- NaN" << endl
			  << tab << "         else \"01\";                                          -- normal number" << endl;

		vhdl  << tab << "r <= r_exn & r_sgn & r_exp & r_frac;"<<endl;


	}




	FPConstMult::FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out):
		Operator(target),
		wE_in(wE_in), wF_in(wF_in), wE_out(wE_out), wF_out(wF_out) 
	{
	}


	FPConstMult::~FPConstMult() {
		// TODO but who cares really
		// clean up memory -- with ifs, because in some cases they have not been init'd
		if(mpfr_cst) mpfr_clear(mpfr_cst);
		if(mpfr_cst_sig) mpfr_clear(mpfr_cst_sig);
		if(mpfr_xcut_sig) mpfr_clear(mpfr_xcut_sig);
	}




	void FPConstMult::emulate(TestCase *tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");

		/* Compute correct value */
		FPNumber fpx(wE_in, wF_in);
		fpx = svX;
		mpfr_t x, r;
		mpfr_init2(x, 1+wF_in);
		mpfr_init2(r, 1+wF_out); 
		fpx.getMPFR(x);
		mpfr_mul(r, x, mpfr_cst, GMP_RNDN);

		// Set outputs 
		FPNumber  fpr(wE_out, wF_out, r);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		// clean up
		mpfr_clears(x, r, NULL);

	}

}
