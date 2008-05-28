/*
 * A correctly-rounded multiplier by an arbitrary real constant for FloPoCo
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

#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../HOTBM/sollya.h" // TODO : fix upstream Sollya, or fix in FloPoCo
#include "../utils.hpp"
#include "../Operator.hpp"
#include "FPConstMult.hpp"
#include "CRFPConstMult.hpp"
#include "../FloFP.hpp"

using namespace std;

extern vector<Operator*> oplist;



CRFPConstMult::CRFPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, string _constant):
	FPConstMult(target, wE_in, wF_in, wE_out, wF_out), 
	constant (_constant) 
{
	sollya_node_t node;
   mpfr_t mpR;
	mpz_t zz;

	/* Convert the input string into a sollya evaluation tree */
	node = parseString(constant.c_str());	/* If conversion did not succeed (i.e. parse error) */
	if (node == 0) {
		cerr << "Unable to parse string "<< constant << " as a numeric constant" <<endl;
		exit(EXIT_FAILURE);
	}

	mpfr_inits(mpR, NULL);
	evaluateConstantExpression(mpR, node,  getToolPrecision());


	if(verbose){
		double r;
		r = mpfr_get_d(mpR, GMP_RNDN);
		cout << "  Constant evaluates to " <<r <<endl; 
	}
	// compute the precision -- TODO with NWB

	cst_width =  2*wF_in+4;
	cerr << "***** WARNING Taking constant with 2*wF_in+4 bits. Correct rounding is not yet guaranteed. This is being implemented." <<endl;


	if(verbose)
		cout << "  Required significand precision to reach correct rounding is " << cst_width<<endl; 

	evaluateConstantExpression(mpR, node,  cst_width);

	// Now convert mpR into exponent + integral significand


	// sign
	cst_sgn = mpfr_sgn(mpR);
	if(cst_sgn<0)
		mpfr_neg(mpR, mpR, GMP_RNDN);

	// compute exponent and mantissa
	cst_exp_when_mantissa_1_2 = mpfr_get_exp(mpR) - 1; //mpfr_get_exp() assumes significand in [1/2,1)  
	cst_exp_when_mantissa_int = cst_exp_when_mantissa_1_2 - cst_width + 1; 
	mpfr_init2(mpfr_cst_sig, cst_width);
	mpfr_div_2si(mpfr_cst_sig, mpR, cst_exp_when_mantissa_1_2, GMP_RNDN);

	if(verbose) {
		cout << "  mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN) <<endl;
	}

	// Build the corresponding FPConstMult.

	
	// initialize mpfr_xcut_sig = 2/cst_sig, will be between 1 and 2
	mpfr_init2(mpfr_xcut_sig, 32*(cst_width+wE_in+wE_out)); // should be accurate enough
	mpfr_set_d(mpfr_xcut_sig, 2.0, GMP_RNDN);               // exaxt op
	mpfr_div(mpfr_xcut_sig, mpfr_xcut_sig, mpfr_cst_sig, GMP_RNDD);

	// now  round it down to wF_in+1 bits 
	mpfr_t xcut_wF;
	mpfr_init2(xcut_wF, wF_in+1);
	mpfr_set(xcut_wF, mpfr_xcut_sig, GMP_RNDD);
	mpfr_mul_2si(xcut_wF, xcut_wF, wF_in, GMP_RNDN);

	// It should now be an int; cast it into a mpz, then a mpz_class 
	mpz_init2(zz, wF_in+1);
	mpfr_get_z(zz, xcut_wF, GMP_RNDN);
	xcut_sig_rd = mpz_class(zz);
   mpz_clear(zz);

	if(verbose) {
		cout << "  mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN) <<endl;
	}
	// Now build the mpz significand
	mpfr_mul_2si(mpfr_cst_sig,  mpfr_cst_sig, cst_width, GMP_RNDN);
   
	// It should now be an int; cast it into a mpz, then a mpz_class 
	mpz_init2(zz, cst_width);
	mpfr_get_z(zz, mpfr_cst_sig, GMP_RNDN);
	cst_sig = mpz_class(zz);
   mpz_clear(zz);


	if(verbose) {
		cout << "  mpzclass cst_sig  = " << cst_sig;
	}

	icm = new IntConstMult(target, wF_in+1, cst_sig);
	oplist.push_back(icm);



	// build the name
	ostringstream name; 
	name <<"FPConstMult_"<<(cst_sgn==0?"":"M") <<cst_sig<<"b"
		  <<(cst_exp_when_mantissa_int<0?"M":"")<<abs(cst_exp_when_mantissa_int)
		  <<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
	unique_name=name.str();

	// Set up the IO signals
	add_FP_input("X", wE_in, wF_in);
	add_FP_output("R", wE_out, wF_out);

	// cleaning up
	mpfr_clears(mpR, mpfr_xcut_sig, xcut_wF, mpfr_cst_sig, NULL);

}



CRFPConstMult::~CRFPConstMult() {
	// TODO but who cares really
}



#endif //HAVE_SOLLYA
