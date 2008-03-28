/*
 * A long accumulator for FloPoCo
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
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "LongAcc.hpp"

using namespace std;

extern vector<Operator*> oplist;

LongAcc::LongAcc(Target* target, int wEX, int wFX, int MaxMSBX, int LSBA, int MSBA): 
  Operator(target), 
  wEX(wEX), wFX(wFX), MaxMSBX(MaxMSBX), LSBA(LSBA), MSBA(MSBA)
{
  ostringstream name; 
  name <<"LongAcc_"<<wEX<<"_"<<wFX<<"_"
       <<(MaxMSBX>=0?"":"M")<<abs(MaxMSBX)<<"_"
       <<(LSBA>=0?"":"M")<<abs(LSBA)<<"_"
       <<(MSBA>=0?"":"M")<<abs(MSBA) ;
  unique_name=name.str();

  // This operator is a sequential one
  set_sequential();

  // Set up various architectural parameters
  sizeAcc = MSBA-LSBA+1;
    
  maxShift = MaxMSBX-LSBA; // shift is 0 when the implicit 1 is at LSBA
  sizeShift = intlog2(maxShift);
  sizeSummand = MaxMSBX-LSBA+1; 
  sizeShiftedFrac = maxShift + wFX+1;
  E0X = (1<<(wEX-1)) -1; // exponent bias
  // The left shift value is MaxMSBX - (EX-E0X) if it is positive
  // If it is negative, it means the summand is shifted out completely: it will be set to 0.
  // So we have to compute (MaxMSBX+E0X) - EX
  // The constant should fit on WEX bits, because MaxMSBX is one possible exponent of X
  // so we have a subtraction of two wEX-bit numbers.
  int biasedMaxMSBX = MaxMSBX + E0X;
  if(biasedMaxMSBX < 0 || intlog2(biasedMaxMSBX)>wEX) {
    cerr << "ERROR in LongAcc: maxMSBX="<<MaxMSBX<<" is not a valid exponent of X (range " << (-E0X) << " to " << ((1<<wEX)-1)-E0X <<endl;
    exit (EXIT_FAILURE);
  }

  // Create an instance of the required input shifter. 
  shifter = new Shifter(target, wFX+1, maxShift, Left);
  oplist.push_back(shifter);

  add_FP_input ("X", wEX,wFX);
  add_output ("A", sizeAcc);  
  //  add_output ("XOverflow");  
  //  add_output ("XUnderflow");  
  //  add_output ("AccOverflow");  

  // Unregistered signal
  add_signal("expX", wEX);
  add_signal("fracX", wFX+1);
  add_signal("shifted_frac", sizeShiftedFrac);
  add_signal("shiftval", wEX+1); //, "includes a sign bit");

  // setup the pipeline 

  if(target->is_pipelined()) {
    // TODO here code to pipeline the 2's complement if the target frequency is high, see IntAdder 
    // meanwhile we just do it in 1 level, assuming there is a register at the out of the shifter.
    c2_pipeline_depth = 1;
  } 
  else {
    c2_pipeline_depth = 0;
  }

  // on one side, add delays for the non-complemented signal
  add_delay_signal("summand", sizeSummand, c2_pipeline_depth); 
  add_delay_signal("flushedToZero", 1, shifter->pipeline_depth());

  //on the other side, add a delay for 2'complement of summand
  add_registered_signal("summand2c", sizeSummand);

  // final pipeline depth is the sum of shifter pipeline and 2's complement pipeline
  set_pipeline_depth(shifter->pipeline_depth() + c2_pipeline_depth);

  add_delay_signal("exnX", 2, pipeline_depth());
  add_delay_signal("signX", 1, pipeline_depth());

  add_registered_signal_with_reset("ext_summand2c", sizeAcc);
  add_registered_signal_with_reset("acc", sizeAcc); //,  "includes overflow bit");

  if(verbose)
    cout << tab <<unique_name<< " pipeline depth is " << pipeline_depth() << " cycles" <<endl;
}

LongAcc::~LongAcc() {
}


void LongAcc::output_vhdl(ostream& o, string name) {
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);
  output_vhdl_entity(o);
  o << "architecture arch of " << name  << " is" << endl;
  shifter->output_vhdl_component(o);
  output_vhdl_signal_declarations(o);
  o << "begin" << endl;
  o << tab << "fracX <= \"1\" & X("<<wFX-1<<" downto 0);" << endl;
  o << tab << "expX <= X("<<wEX+wFX-1<<" downto "<<wFX<<");" << endl;
  o << tab << "signX <= X("<<wEX+wFX<<");" << endl;
  o << tab << "exnX <= X("<<wEX+wFX+2<<" downto "<<wEX+wFX+1<<");" << endl;
  int64_t exp_offset = E0X+LSBA;
  if(exp_offset >=0) {
    o << tab << "shiftval <=  (\"0\"&expX) - \"";  printBinNum(o, exp_offset,  wEX+1); o << "\";" << endl;
  }
  else {
    o << tab << "shiftval <=  (\"0\"&expX) + \"";  printBinNum(o, -exp_offset,  wEX+1); o  << "\";" << endl;
  }
  o << endl;
  o << tab << "-- shift of the input into the proper place " << endl;
  o << tab << "input_shifter: " << shifter->unique_name << endl;
  o << tab << "    port map ( X => fracX, " << endl;
  o << tab << "               S => shiftval("<< shifter->wShiftIn - 1 <<" downto 0), " << endl;
  o << tab << "               R => shifted_frac";
  if (shifter->is_sequential()) {
    o<<"," << endl;
    o << tab << "             clk => clk," << endl;
    o << tab << "             rst => rst " << endl;
  }
  o << tab << "             );" << endl; 
  o << endl;
 
  o << tab << "flushedToZero <=     '1' when (shiftval("<<wEX<<")='1' -- negative left shift " << endl;
  o << tab << "                               or exnX=\"00\")" << endl;
  o << tab << "                 else'0';" << endl;
  o << tab << "summand <= ("<<sizeSummand-1<<" downto 0 => '0')  when "<< get_delay_signal_name("flushedToZero", shifter->pipeline_depth()) << "='1'  else shifted_frac("<<sizeShiftedFrac-1<<" downto "<<wFX<<");" << endl;
  o << endl;
  o << tab << "-- 2's complement of the summand" << endl;
  // This is the line that should be pipelined
  o << tab << "summand2c <= summand when "<< get_delay_signal_name("signX", shifter->pipeline_depth()) <<"='0' else ("<<sizeSummand-1<<" downto 0 => '0') - summand; "<< endl;
  o << endl;
  o << tab << "-- extension of the summand to accumulator size" << endl;
  o << tab << "ext_summand2c <= ("<<sizeAcc-1<<" downto "<<sizeSummand<<"=>'0') & summand2c   when  " << get_delay_signal_name("signX", shifter->pipeline_depth()) <<"='0'" << endl;
  o << tab << "            else ("<<sizeAcc-1<<" downto "<<sizeSummand<<"=> not "<< get_delay_signal_name("flushedToZero", shifter->pipeline_depth()) << ") & summand2c;" << endl;
  o << endl;
  o << tab << "-- accumulation itself" << endl;
  o << tab << "acc <= ext_summand2c_d   +   acc_d;" << endl;
  o << endl;

  output_vhdl_registers(o);

  o << tab << "  A <=   acc_d;" << endl;
  
  //TODO sticky overflow for the accumulator.
  //TODO sticky overflow for the input.
  //TODO sticky underflow for the input.

  o << "end architecture;" << endl << endl;
    
}



void LongAcc::test_precision(int n) {
  mpfr_t ref_acc, long_acc, fp_acc, r, d, one, two, msb;
  double sum, error;

  // initialisation
  mpfr_init2(ref_acc, 10000);
  mpfr_init2(long_acc, sizeAcc+1);
  mpfr_init2(fp_acc, wFX+1);
  mpfr_init2(r, wFX+1);
  mpfr_init2(d, 100);
  mpfr_init2(one, 100);
  mpfr_init2(two, 100);
  mpfr_init2(msb, 100);
  mpfr_set_d(one, 1.0, GMP_RNDN);
  mpfr_set_d(two, 2.0, GMP_RNDN);
  mpfr_set_d(msb, (double)(1<<(MSBA+1)), GMP_RNDN); // works for MSBA<32

  //cout<<"%-------Acc. of positive numbers--------------- "<<endl;
  mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
  mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
  //put a one in the MSBA+1 bit will turn all the subsequent additions
  // into fixed-point ones
  mpfr_set_d(long_acc, (double)(1<<(MSBA+1)), GMP_RNDN);

  for(int i=0; i<n; i++){
    mpfr_random(r); // deprecated function; r is [0,1]
    //    mpfr_add(r, one, r, GMP_RNDN); // r in [1 2[
    
    mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
    mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);
    mpfr_add(long_acc, long_acc, r, GMP_RNDN);
    if(mpfr_greaterequal_p(long_acc, msb)) 
      mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);
      
  }

  // remove the leading one from long acc
  if(mpfr_greaterequal_p(long_acc, msb)) 
    mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);
    
  sum=mpfr_get_d(ref_acc, GMP_RNDN);
  cout  << "% unif[0 1] :sum="<< sum;
  sum=mpfr_get_d(fp_acc, GMP_RNDN);
  cout << "   FPAcc="<< sum;
  sum=mpfr_get_d(long_acc, GMP_RNDN);
  cout << "   LongAcc="<< sum;

  cout <<endl << n << " & ";
  // compute the error for the FP adder
  mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
  mpfr_div(d, d, ref_acc, GMP_RNDN);
  error=mpfr_get_d(d, GMP_RNDN);
  // cout << " Relative error between fp_acc and ref_acc is "<< error << endl;
  cout << scientific << setprecision(2)  << error << " & ";
  // compute the error for the long acc
  mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
  mpfr_div(d, d, ref_acc, GMP_RNDN);
  error=mpfr_get_d(d, GMP_RNDN);
  //  cout << "Relative error between long_acc and ref_acc is "<< error << endl;
  cout << scientific << setprecision(2)  << error << " & ";


    //cout<<"%-------Acc. of positive/negative numbers--------------- "<<endl;

  mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
  mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
  //put a one in the MSBA+1 bit will turn all the subsequent additions
  // into fixed-point ones
  mpfr_set_d(long_acc, (double)(1<<(MSBA+1)), GMP_RNDN);

  for(int i=0; i<n; i++){
    mpfr_random(r); // deprecated function; r is [0,1]
    mpfr_mul(r, r, two, GMP_RNDN); 
    mpfr_sub(r, r, one, GMP_RNDN); 
    
    mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
    mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);
    mpfr_add(long_acc, long_acc, r, GMP_RNDN);
  }

  // remove the leading one from long acc
  mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);


  // compute the error for the FP adder
  mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
  mpfr_div(d, d, ref_acc, GMP_RNDN);
  error=mpfr_get_d(d, GMP_RNDN);
  // cout << "Relative error between fp_acc and ref_acc is "<< error << endl;
  cout << scientific << setprecision(2) << error << " & ";

  // compute the error for the long acc
  mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
  mpfr_div(d, d, ref_acc, GMP_RNDN);
  error=mpfr_get_d(d, GMP_RNDN);
  //  cout << "Relative error between long_acc and ref_acc is "<< error << endl;
  cout << scientific << setprecision(2)  << error << " \\\\ \n     \\hline \n";

  sum=mpfr_get_d(ref_acc, GMP_RNDN);
  cout << "% unif[-1 1] : sum="<< sum;
  sum=mpfr_get_d(fp_acc, GMP_RNDN);
   cout << "   FPAcc="<< sum;
  sum=mpfr_get_d(long_acc, GMP_RNDN);
   cout << "   LongAcc="<< sum;
    cout <<endl;

}



// read the values from a file and accumulate them
void LongAcc::test_precision2() {

  mpfr_t ref_acc, long_acc, fp_acc, r, d, one, two, msb;
  double dr, dacc, sum, error;

  // initialisation
#define DPFPAcc 1
  mpfr_init2(ref_acc, 10000);
  mpfr_init2(long_acc, sizeAcc+1);
#if DPFPAcc
  mpfr_init2(fp_acc,52 +1);
#else
  mpfr_init2(fp_acc, wFX+1);
#endif
  mpfr_init2(r, wFX+1);
  mpfr_init2(d, 100);
  mpfr_init2(one, 100);
  mpfr_init2(two, 100);
  mpfr_init2(msb, 100);
  mpfr_set_d(one, 1.0, GMP_RNDN);
  mpfr_set_d(two, 2.0, GMP_RNDN);
  mpfr_set_d(msb, (double)(1<<(MSBA+1)), GMP_RNDN); // works for MSBA<32

  mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
  mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
  //put a one in the MSBA+1 bit will turn all the subsequent additions
  // into fixed-point ones
  mpfr_set_d(long_acc, (double)(1<<(MSBA+1)), GMP_RNDN);


  ifstream myfile ("/home/fdedinec/accbobines.txt");
  int i=0;
  if (myfile.is_open())
  {
    while (! myfile.eof() )
    {
      myfile >> dr;
      mpfr_set_d(r, dr, GMP_RNDN);
      mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
      mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);
      mpfr_add(long_acc, long_acc, r, GMP_RNDN);
      i++;
      if (i%100000==0) {
	cout<<"i="<<i<<" "<<endl;

	// compute the error for the FP adder
	mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
	mpfr_div(d, d, ref_acc, GMP_RNDN);
	error=mpfr_get_d(d, GMP_RNDN);
	cout << "Relative error between fp_acc and ref_acc is "<< error << endl;
	//cout << scientific << setprecision(2) << error << " & ";

	// remove the leading one from long acc
	mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);
	// compute the error for the long acc
	mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
	mpfr_div(d, d, ref_acc, GMP_RNDN);
	error=mpfr_get_d(d, GMP_RNDN);
	cout << "Relative error between long_acc and ref_acc is "<< error << endl;
	//cout << scientific << setprecision(2)  << error << " \\\\ \n     \\hline \n";

	sum=mpfr_get_d(ref_acc, GMP_RNDN);
	cout << " exact sum="<< sum;
	sum=mpfr_get_d(fp_acc, GMP_RNDN);
	cout << "   FPAcc="<< sum;
	sum=mpfr_get_d(long_acc, GMP_RNDN);
	cout << "   LongAcc="<< sum;
	cout <<endl;
	// add the leading one back
	mpfr_add(long_acc, long_acc, msb, GMP_RNDN);
      }
    }
    myfile.close();
  }
  // remove the leading one from long acc
  mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);


  // compute the error for the FP adder
  mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
  mpfr_div(d, d, ref_acc, GMP_RNDN);
  error=mpfr_get_d(d, GMP_RNDN);
  cout << "Relative error between fp_acc and ref_acc is "<< error << endl;
   //cout << scientific << setprecision(2) << error << " & ";

  // compute the error for the long acc
  mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
  mpfr_div(d, d, ref_acc, GMP_RNDN);
  error=mpfr_get_d(d, GMP_RNDN);
    cout << "Relative error between long_acc and ref_acc is "<< error << endl;
    //cout << scientific << setprecision(2)  << error << " \\\\ \n     \\hline \n";

  sum=mpfr_get_d(ref_acc, GMP_RNDN);
  cout << " exact sum="<< sum;
  sum=mpfr_get_d(fp_acc, GMP_RNDN);
   cout << "   FPAcc="<< sum;
  sum=mpfr_get_d(long_acc, GMP_RNDN);
   cout << "   LongAcc="<< sum;
    cout <<endl;


  exit(0);


}
