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
#include "../FloFP.hpp"

using namespace std;

extern vector<Operator*> oplist;






FPConstMult::FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, int cst_sgn, int cst_exp, mpz_class cst_sig):
	Operator(target), 
	wE_in(wE_in), wF_in(wF_in), wE_out(wE_out), wF_out(wF_out), 
	cst_sgn(cst_sgn), cst_exp_when_mantissa_int(cst_exp), cst_sig(cst_sig)
{
	ostringstream name; 
	name <<"FPConstMult_"<<(cst_sgn==0?"":"M") <<cst_sig<<"b"<<(cst_exp<0?"M":"")<<abs(cst_exp)<<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
	unique_name=name.str();

	int cst_width = intlog2(cst_sig);
	cst_exp_when_mantissa_1_2 = cst_exp_when_mantissa_int + cst_width - 1; 

	// TODO if the constant is zero or a power of two
	// TODO normalize the constant significand

	// initialize mpfr constant
	// all this is ugly because no mpfr equivalent of mpz_class
	mpfr_init2(mpfr_cst_sig, cst_width);
	mpfr_set_z(mpfr_cst_sig, cst_sig.get_mpz_t(), GMP_RNDN); // exact op
	mpfr_mul_2si(mpfr_cst_sig, mpfr_cst_sig, -(cst_width-1), GMP_RNDN);  // exact op
	
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

	if(verbose) {
		cout << "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN) <<endl;
		cout << "mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN) <<endl;
		cout << "xcut_sig_rd   = " << xcut_sig_rd << "   ";
		printBinNumGMP(cout, xcut_sig_rd, wF_in+1);  cout << endl;
	}

 	/* Initialize second operand */
 	mpfr_init(mpY);
 	mpfr_set(mpY, mpfr_cst_sig, GMP_RNDN);
 	mpfr_mul_2si(mpY, mpY, cst_exp_when_mantissa_1_2, GMP_RNDN);

	// do all the declarations. Pushed into a method so that CRFPConstMult can inherit it
	icm = new IntConstMult(target, wF_in+1, cst_sig);
	oplist.push_back(icm);

	if (target->is_pipelined())
		set_sequential();
	else 
		set_combinatorial();
	setup();
}



void FPConstMult::setup() {


 	// Set up the IO signals
 	add_FP_input("X", wE_in, wF_in);
 	add_FP_output("R", wE_out, wF_out);


	// Set up non-registered signals 
	add_signal_bus      ("x_sig",          wF_in+1);
	add_signal_bus      ("sig_prod",      icm->rsize);   // register the output of the int const mult
	add_signal_bus      ("rounded_frac",   wF_out+1);
	add_signal          ("overflow");
	add_signal          ("underflow");
	add_signal_bus      ("r_frac", wF_out);
 
	// Set up other signals and pipeline-related stuff
 	if (is_sequential()) {
 		icm_depth = icm->pipeline_depth();
 		set_pipeline_depth(icm_depth+1);

		add_delay_signal_bus_no_reset("x_exn",          2,         1);
		add_delay_signal_no_reset    ("x_sgn",          1,         1);
		add_delay_signal_bus_no_reset("x_exp",          wE_in,     1);
		add_delay_signal_bus_no_reset("shifted_frac",    wF_out+1,  1);

		add_delay_signal_no_reset    ("gt_than_xcut",   1,         icm_depth);

		add_delay_signal_bus_no_reset("r_exp_nopb",    wE_out+1,   1);

		add_delay_signal_bus_no_reset("r_exn",         2,          icm_depth);
		add_delay_signal_no_reset    ("r_sgn",         1,          icm_depth);
		add_delay_signal_bus_no_reset("r_exp",         wE_out,     icm_depth);


 	}
 	else {
		add_signal_bus("x_exn",          2);
		add_signal    ("x_sgn"           );
		add_signal_bus("x_exp",          wE_in);

		add_signal_bus      ("shifted_frac",   wF_out+1);
		add_signal    ("gt_than_xcut");

		add_signal_bus("r_exp_nopb",    wE_out+1);

		add_signal_bus      ("r_exn", 2);
		add_signal          ("r_sgn");
		add_signal_bus      ("r_exp", wE_out);
 	}
}




FPConstMult::FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out):
	Operator(target),
	wE_in(wE_in), wF_in(wF_in), wE_out(wE_out), wF_out(wF_out) 
{
}


FPConstMult::~FPConstMult() {
	// TODO but who cares really
	// delete icm; Better not: it has been added to oplist
	// mpfr_clear(mpfr_xcut_sig);
	// mpfr_clear(mpfr_cst_sig);
	mpfr_clear(mpY);
}


void FPConstMult::output_vhdl(ostream& o, string name) {
	Licence(o,"Florent de Dinechin (2007)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	o << "architecture arch of " << name  << " is" << endl;

	icm->Operator::output_vhdl_component(o);

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
	output_vhdl_signal_declarations(o);


	o << tab << "signal abs_unbiased_cst_exp  : std_logic_vector("<<wE_sum<<" downto 0) := \"";
	printBinNumGMP(o,  expAddend, wE_sum+1);    o << "\";"<<endl;
 	o << tab << "signal xcut_rd      : std_logic_vector("<<wF_in<<" downto 0) := \"";
 	printBinPosNumGMP(o, xcut_sig_rd, wF_in+1);     o<<"\";"<<endl;

	o << tab << "begin"<<endl;
	o << tab << tab << "x_exn <=  x("<<wE_in<<"+"<<wF_in<<"+2 downto "<<wE_in<<"+"<<wF_in<<"+1);"<<endl;
	o << tab << tab << "x_sgn <=  x("<<wE_in<<"+"<<wF_in<<");"<<endl;
	o << tab << tab << "x_exp <=  x("<<wE_in<<"+"<<wF_in<<"-1 downto "<<wF_in<<");"<<endl;
	o << tab << tab << "x_sig <= '1' & x("<<wF_in-1 <<" downto 0);"<<endl;
	o << tab << tab << "sig_mult : "
	  << icm->unique_name<<endl
	  << tab<<tab<<tab << "port map(inx => x_sig, r => sig_prod";
	if(is_sequential())
		o << ", clk => clk, rst => rst";
	o << ");"<<endl;
	// Possibly shift the significand one bit left, and remove implicit 1 
	o << tab << tab << "gt_than_xcut <= '1' when ( x_sig("<<wF_in-1<<" downto 0) > xcut_rd("<<wF_in-1<<" downto 0) ) else '0';"<<endl;
	o << tab << tab << "shifted_frac <=  sig_prod("<<icm->rsize -2<<" downto "<<icm->rsize - wF_out - 2<<")  when " << get_delay_signal_name("gt_than_xcut", icm_depth) << " = '1'"<<endl
		<< tab << tab << "           else sig_prod("<<icm->rsize -3<<" downto "<<icm->rsize - wF_out - 3<<");"<<endl;  
	// add the rounding bit
	o << tab << tab << "rounded_frac <= (("<<wF_out <<" downto 1 => '0') & '1') + " << get_delay_signal_name("shifted_frac", 1) << ";"<<endl;
	o << tab << tab << "r_frac <= rounded_frac("<<wF_out <<" downto  1);"<<endl;
	// Handling signs is trivial
	if(cst_sgn==0)
		o << tab << tab << "r_sgn <= " << get_delay_signal_name("x_sgn",1) << "; -- positive constant"<<endl;
	else
		o << tab << tab << "r_sgn <= not " << get_delay_signal_name("x_sgn",1) << "; -- negative constant"<<endl;

	// exponent handling
	o << tab << tab << "r_exp_nopb <= ";
	o << "((" << wE_sum << " downto " << wE_in << " => '0')  & " << get_delay_signal_name("x_exp", 1) << ")  ";
	o << (expAddendSign==0 ? "+" : "-" ) << "  abs_unbiased_cst_exp";
	o << "  +  (("<<wE_sum<<" downto 1 => '0') & " <<  get_delay_signal_name("gt_than_xcut", 1) << ");"<<endl;

	// overflow handling
	if (maxExp(wE_in) + cst_exp_when_mantissa_1_2 + 1 < maxExp(wE_out)) // no overflow can ever happen
		o << tab << tab << "overflow <= '0'; --  overflow never happens for this constant and these (wE_in, wE_out)" << endl;
	else 
		o << tab << tab << "overflow <= '0' when r_exp_nopb(" << wE_sum << " downto " << wE_out << ") = (" << wE_sum << " downto " << wE_out <<" => '0')     else '1';" << endl;

	// underflow handling
	if (minExp(wE_in) + cst_exp_when_mantissa_1_2 > minExp(wE_out)) // no overflow can ever happen
		o << tab << tab << "underflow <= '0'; --  underflow never happens for this constant and these (wE_in, wE_out)" << endl;
	else 
		o << tab << tab << "underflow <= '0' when r_exp_nopb(" << wE_sum << ") = '1'     else '1';" << endl;
			 
	
	o << tab << tab << "r_exp <= r_exp_nopb("<<wE_out-1<<" downto 0) ;"<<endl;

	o << tab << tab << "r_exn <=      \"00\" when ((x_exn = \"00\") or (underflow='1'))  -- zero"<<endl 
		<< tab << tab << "         else \"10\" when ((x_exn = \"10\") or (overflow='1'))   -- infinity"<<endl
		<< tab << tab << "         else \"11\" when  (x_exn = \"11\")                      -- NaN"<<endl
		<< tab << tab << "         else \"01\";                                          -- normal number"<<endl;

	output_vhdl_registers(o);

	o << tab << tab << "r <= " << get_delay_signal_name("r_exn", icm_depth) << " & " << get_delay_signal_name("r_sgn", icm_depth) << " & r_exp & r_frac;"<<endl;
	o << "end architecture;" << endl << endl;		
}

TestIOMap FPConstMult::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("X"));
	tim.add(*get_signal_by_name("R"));
	return tim;
}

void FPConstMult::fillTestCase(mpz_class a[])
{
	mpz_class& svX = a[0];
	mpz_class& svR = a[1];
	FloFP x(wE_in, wF_in), r(wE_out, wF_out);
	x = svX;
	r = x * mpY;
	svR = r.getSignalValue();
}

