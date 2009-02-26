/*
 * A multiplier by a floating-point constant for FloPoCo
 *
 * Author : Florent de Dinechin
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/


// TODO  we discard the lower bits, try not to compute them in IntConstMult, or at least not to register them.

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

extern vector<Operator*> oplist;






FPConstMult::FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, int cst_sgn, int cst_exp, mpz_class cst_sig):
	Operator(target), 
	wE_in(wE_in), wF_in(wF_in), wE_out(wE_out), wF_out(wF_out), 
	cst_sgn(cst_sgn), cst_exp_when_mantissa_int(cst_exp), cst_sig(cst_sig)
{
	ostringstream name;

#ifdef _WIN32
	name <<"FPConstMult_"<<(cst_sgn==0?"":"M") <<mpz2string(cst_sig)<<"b"<<(cst_exp<0?"M":"")<<abs(cst_exp)<<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
#else
	name <<"FPConstMult_"<<(cst_sgn==0?"":"M") <<cst_sig<<"b"<<(cst_exp<0?"M":"")<<abs(cst_exp)<<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
#endif
	uniqueName_=name.str();

	int cst_width = intlog2(cst_sig);
	cst_exp_when_mantissa_1_2 = cst_exp_when_mantissa_int + cst_width - 1; 

	// TODO if the constant is zero or a power of two
	// TODO normalize the constant significand

	// initialize mpfr constant
	// all this is ugly because no mpfr equivalent of mpz_class
	mpfr_init2(mpfr_cst_sig, cst_width);
	mpfr_set_z(mpfr_cst_sig, cst_sig.get_mpz_t(), GMP_RNDN); // exact op
	mpfr_mul_2si(mpfr_cst_sig, mpfr_cst_sig, -(cst_width-1), GMP_RNDN);  // exact op, gets mpfr_cst_sig in 1..2
	
 	mpfr_init(mpfr_cst);
 	mpfr_set(mpfr_cst, mpfr_cst_sig, GMP_RNDN);
 	mpfr_mul_2si(mpfr_cst, mpfr_cst, cst_exp_when_mantissa_1_2, GMP_RNDN);

	if(cst_sgn==1)
		mpfr_neg( mpfr_cst,  mpfr_cst, GMP_RNDN);
	
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

#ifdef _WIN32
	if(verbose) {
		cout << "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN) <<endl;
		cout << "mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN) <<endl;
		//GMP's C++ wrapper problem
		cout << "xcut_sig_rd   = " << mpz2string(xcut_sig_rd) << "   ";
		printBinNumGMP(cout, xcut_sig_rd, wF_in+1);  cout << endl;
	}
#else
	if(verbose) {
		cout << "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN) <<endl;
		cout << "mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN) <<endl;
		cout << "xcut_sig_rd   = " << xcut_sig_rd << "   ";
		printBinNumGMP(cout, xcut_sig_rd, wF_in+1);  cout << endl;
	}
#endif


	// do all the declarations. Pushed into a method so that CRFPConstMult can inherit it
	icm = new IntConstMult(target, wF_in+1, cst_sig);
	oplist.push_back(icm);

	setup();



}



void FPConstMult::setup() {


 	// Set up the IO signals
 	addFPInput("X", wE_in, wF_in);
 	addFPOutput("R", wE_out, wF_out);


	setCopyrightString("Florent de Dinechin (2007)");

	// bit width of constant exponent
	int wE_cst=intlog2(abs(cst_exp_when_mantissa_1_2));
	if(verbose)
		cout << "  wE_cst="<<wE_cst<<endl;
	
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


 	vhdl << tab << declare("xcut_rd", wF_in+1) << " <= \""
		  << unsignedBinary(xcut_sig_rd, wF_in+1) << "\";"<<endl;
	
	vhdl << tab << declare("x_exn",2) << " <=  x("<<wE_in<<"+"<<wF_in<<"+2 downto "<<wE_in<<"+"<<wF_in<<"+1);"<<endl;
	vhdl << tab << declare("x_sgn") << " <=  x("<<wE_in<<"+"<<wF_in<<");"<<endl;
	vhdl << tab << declare("x_exp", wE_in) << " <=  x("<<wE_in<<"+"<<wF_in<<"-1 downto "<<wF_in<<");"<<endl;
	vhdl << tab << declare("x_sig", wF_in+1) << " <= '1' & x("<<wF_in-1 <<" downto 0);"<<endl;

	vhdl << tab << declare("gt_than_xcut") << " <= '1' when ( x_sig("<<wF_in-1<<" downto 0) > xcut_rd("<<wF_in-1<<" downto 0) ) else '0';"<<endl;



	inPortMap  (icm, "inX", "x_sig");
	outPortMap (icm, "R","sig_prod");
	vhdl << instance(icm, "sig_mult");
	setCycleFromSignal("sig_prod"); 
	nextCycle();

	// Possibly shift the significand one bit left, and remove implicit 1 
	vhdl << tab << declare("shifted_frac",    wF_out+1) << " <= " << use("sig_prod") << "("<<icm->rsize -2<<" downto "<<icm->rsize - wF_out-2 <<")  when " << use("gt_than_xcut") << " = '1'"<<endl
		  << tab << "           else " << use("sig_prod") << "("<<icm->rsize -3<<" downto "<<icm->rsize - wF_out - 3<<");"<<endl;  
	// add the rounding bit
	vhdl << tab << tab << declare("rounded_frac",   wF_out+1) << " <= (("<<wF_out <<" downto 1 => '0') & '1') + " << use("shifted_frac") << ";"<<endl;
	vhdl << tab << tab << declare("r_frac", wF_out) << " <= rounded_frac("<<wF_out <<" downto  1);"<<endl;
	// Handling signs is trivial
	if(cst_sgn==0)
		vhdl << tab << declare("r_sgn") << " <= " << use("x_sgn") << "; -- positive constant"<<endl;
	else
		vhdl << tab << "r_sgn <= not " << use("x_sgn") << "; -- negative constant"<<endl;

	// exponent handling
	vhdl << tab << declare("abs_unbiased_cst_exp",wE_sum+1) << " <= \""
		  << unsignedBinary(expAddend, wE_sum+1) << "\";" << endl;
	vhdl << tab << declare("r_exp_nopb",    wE_out+1) << " <= "
		  << "((" << wE_sum << " downto " << wE_in << " => '0')  & " << use("x_exp") << ")  "
		  << (expAddendSign==0 ? "+" : "-" ) << "  abs_unbiased_cst_exp"
		  << "  +  (("<<wE_sum<<" downto 1 => '0') & " <<  use("gt_than_xcut") << ");"<<endl;

	// overflow handling
	vhdl << tab << declare("overflow") << " <= " ;
	if (maxExp(wE_in) + cst_exp_when_mantissa_1_2 + 1 < maxExp(wE_out)) // no overflow can ever happen
		vhdl << "'0'; --  overflow never happens for this constant and these (wE_in, wE_out)" << endl;
	else 
		vhdl <<  "'0' when r_exp_nopb(" << wE_sum << " downto " << wE_out << ") = (" << wE_sum << " downto " << wE_out <<" => '0')     else '1';" << endl;

	// underflow handling
	vhdl << tab << declare("underflow") << " <= " ;
	if (minExp(wE_in) + cst_exp_when_mantissa_1_2 > minExp(wE_out)) // no overflow can ever happen
		vhdl << "'0'; --  underflow never happens for this constant and these (wE_in, wE_out)" << endl;
	else 
		vhdl <<  "r_exp_nopb(" << wE_sum << ");" << endl;
			 
	
	vhdl << tab << declare("r_exp", wE_out) << " <= r_exp_nopb("<<wE_out-1<<" downto 0) ;"<<endl;

	vhdl << tab << declare("r_exn", 2) << " <=      \"00\" when (("<< use("x_exn")<<" = \"00\") or ("<<use("x_exn")<<" = \"01\" and underflow='1'))  -- zero"<<endl 
		  << tab << "         else \"10\" when (("<<use("x_exn")<<" = \"10\") or ("<<use("x_exn")<<" = \"01\" and overflow='1'))   -- infinity"<<endl
		  << tab << "         else \"11\" when  ("<<use("x_exn")<<" = \"11\")                      -- NaN"<<endl
		  << tab << "         else \"01\";                                          -- normal number"<<endl;

	vhdl  << tab << "r <= " << use("r_exn") << " & " << use("r_sgn") << " & " << use("r_exp") << " & r_frac;"<<endl;


}




FPConstMult::FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out):
	Operator(target),
	wE_in(wE_in), wF_in(wF_in), wE_out(wE_out), wF_out(wF_out) 
{
}


FPConstMult::~FPConstMult() {
	// TODO but who cares really
	// delete icm; Better not: it has been added to oplist
	mpfr_clears(mpfr_cst, mpfr_cst_sig, mpfr_xcut_sig, NULL);
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

