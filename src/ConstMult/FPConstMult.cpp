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


#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "FPConstMult.hpp"

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

  icm = new IntConstMult(target, wF_in+1, cst_sig);
  oplist.push_back(icm);

  // Set up the IO signals
  add_FP_input("X", wE_in, wF_in);
  add_FP_output("R", wE_out, wF_out);
}



FPConstMult::~FPConstMult() {
  // TODO but who cares really
  // delete icm; Better not: it has been added to oplist
  // mpfr_free(mpfr_xcut_sig);
  // mpfr_free(mpfr_cst_sig);
}



// TODO this one is nicer because it shows the FP ins as FP. Upgrade Signal this way.
#if 0
void FPConstMult::output_vhdl_component(ostream& o, string name)
{
  o  << tab << "component " << name << " is" << endl
     << tab << "   port ( x : in  std_logic_vector(" << wE_in << "+"<< wF_in << "+2 downto 0);" << endl
     << tab << "          r : out std_logic_vector(" << wE_out << "+"<< wF_out << "+2 downto 0) );" << endl
     << tab << "end component;" << endl << endl;

}
#endif

void FPConstMult::output_vhdl(ostream& o, string name) {
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);
#if 0
  o << "entity  " << name << " is" << endl
    << "    port ( x : in  std_logic_vector(" << wE_in << "+"<< wF_in << "+2 downto 0);" << endl
    << "           r : out std_logic_vector(" << wE_in << "+"<< wF_in << "+2 downto 0) );" << endl
    << "end entity;" << endl << endl;
#else
  output_vhdl_entity(o);
#endif
  o << "architecture arch of " << name  << " is" << endl;

  icm->Operator::output_vhdl_component(o);

  // bit width of constant exponent
  int wE_cst=intlog2(abs(cst_exp_when_mantissa_1_2));
  cout << "wE_cst="<<wE_cst<<endl;

  int wE_sum;
  // bias_in is 2^(wE_in-1) -1  (011...11)
  // define wE_sum  such that the sum of x_exp and cst_exp is always representable on wE_sum bits
  if (wE_cst < wE_in)
    wE_sum = wE_in+1;
  else {
    wE_sum = wE_cst+1;
    cerr<<"TODO: exponent datapath wrong when cst_exp too large"<<endl;
  }

  o << tab << "signal x_exn        : std_logic_vector(1 downto 0);"<<endl;
  o << tab << "signal x_sgn        : std_logic;"<<endl;
  o << tab << "signal x_exp        : std_logic_vector("<<wE_in-1<<" downto 0);"<<endl;  
  o << tab << "signal x_sig        : std_logic_vector("<<wF_in<<" downto 0);"<<endl;

  o << tab << "signal xcut_rd      : std_logic_vector("<<wF_in<<" downto 0) := \"";
  printBinPosNumGMP(o, xcut_sig_rd, wF_in+1);     o<<"\";"<<endl;

  if(wE_in==wE_out){
    // TODO add check that the exponent fits the representation
    o << tab << "signal abs_cst_exp  : std_logic_vector("<<wE_in<<" downto 0) := \"";
    printBinNum(o,  abs(cst_exp_when_mantissa_1_2), wE_in+1);    o << "\";"<<endl;
  }
  else {
    cerr << "TODO: exponent datapath not implemented when wE_in != wE_out";
    o << tab << "signal exp_bias_in  : std_logic_vector("<<wE_cst-1<<" downto 0) := 2^(wE_in-1) -1";
  }

  o << tab << "signal gt_than_xcut : std_logic;"<<endl;
  o << tab << "signal sig_prod     : std_logic_vector("<<icm->rsize -1<<" downto 0);"<<endl;
  o << tab << "signal shifted_sig  : std_logic_vector("<<wF_out<<" downto 0);"<<endl;  
  o << tab << "signal rounded_sig  : std_logic_vector("<<wF_out<<" downto 0);"<<endl;  

  o << tab << "signal r_exp_nopb   : std_logic_vector("<<wE_out<<" downto 0);  -- no overflow or underflow "<<endl;  
  o << tab << "signal overflow     : std_logic;"<<endl;
  o << tab << "signal underflow    : std_logic;"<<endl;

  o << tab << "signal r_exn        : std_logic_vector(1 downto 0);"<<endl;
  o << tab << "signal r_sgn        : std_logic;"<<endl;
  o << tab << "signal r_exp        : std_logic_vector("<<wE_out-1<<" downto 0);"<<endl;  
  o << tab << "signal r_sig        : std_logic_vector("<<wF_out-1<<" downto 0);"<<endl;  
  o << tab << "begin"<<endl;
  o << tab << tab << "x_exn <=  x("<<wE_in<<"+"<<wF_in<<"+2 downto "<<wE_in<<"+"<<wF_in<<"+1);"<<endl;
  o << tab << tab << "x_sgn <=  x("<<wE_in<<"+"<<wF_in<<");"<<endl;
  o << tab << tab << "x_exp <=  x("<<wE_in<<"+"<<wF_in<<"-1 downto "<<wF_in<<");"<<endl;
  o << tab << tab << "x_sig <= '1' & x("<<wF_in-1 <<" downto 0);"<<endl;
  o << tab << tab << "sig_mult : "<<icm->unique_name<<endl<<tab<<tab<<tab<<"port map(x => x_sig, r => sig_prod);"<<endl;
  // Possibly shift the significand one bit left, and remove implicit 1 
  o << tab << tab << "gt_than_xcut <= '1' when ( x_sig("<<wF_in-1<<" downto 0) > xcut_rd("<<wF_in-1<<" downto 0) ) else '0';"<<endl;
  o << tab << tab << "shifted_sig <=  sig_prod("<<icm->rsize -2<<" downto "<<icm->rsize - wF_out - 2<<")  when gt_than_xcut = '1'"<<endl
    << tab << tab << "           else sig_prod("<<icm->rsize -3<<" downto "<<icm->rsize - wF_out - 3<<");"<<endl;  
  // add the rounding bit
  o << tab << tab << "rounded_sig <= (("<<wF_out-1 <<" downto 0 => '0') & '1') + shifted_sig;"<<endl;
  o << tab << tab << "r_sig <= rounded_sig("<<wF_out <<" downto  1);"<<endl;
  // Handling signs is trivial
  if(cst_sgn==0)
    o << tab << tab << "r_sgn <= x_sgn; -- positive constant"<<endl;
  else
    o << tab << tab << "r_sgn <= not x_sgn; -- negative constant"<<endl;
  // Handling exponents is slightly more complex
  if(wE_in==wE_out){
    o << tab << tab << "r_exp_nopb <= ('0' & x_exp)  "
      <<(cst_exp_when_mantissa_1_2 >= 0 ? "+" : "-" )<<"  abs_cst_exp"
      <<"  +  (("<<wE_in<<" downto 1 => '0') & gt_than_xcut);"<<endl;
    if (cst_exp_when_mantissa_1_2 >= 0) {     // Is the constant greater than 1?
      o << tab << tab << "overflow  <= r_exp_nopb("<<wE_out<<") ;"<<endl;
      o << tab << tab << "underflow <= '0'; --  constant greater than 1, underflow never happens"<<endl;
    }
    else { // constant strictly smaller than 1
      o << tab << tab << "overflow  <= '0'; --  constant strictly smaller than 1, overflow never happens"<<endl;  
      o << tab << tab << "underflow <= r_exp_nopb("<<wE_out<<") ;"<<endl;
    }
  }
  else {
    cerr << "TODO: exponent datapath not implemented when wE_in != wE_out";
  }
  
  o << tab << tab << "r_exp <= r_exp_nopb("<<wE_out-1<<" downto 0) ;"<<endl;

  o << tab << tab << "r_exn <=      \"00\" when ((x_exn = \"00\") or (underflow='1'))  -- zero"<<endl 
    << tab << tab << "         else \"10\" when ((x_exn = \"10\") or (overflow='1'))   -- infinity"<<endl
    << tab << tab << "         else \"11\" when  (x_exn = \"11\")                      -- NaN"<<endl
    << tab << tab << "         else \"01\";                                          -- normal number"<<endl;

  o << tab << tab << "r <= r_exn & r_sgn & r_exp & r_sig;"<<endl;
  o << "end architecture;" << endl << endl;
    
}
