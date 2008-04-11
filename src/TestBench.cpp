/*
 * A generic wrapper generator for FloPoCo. 
 *
 * A wrapper is a VHDL entity that places registers before and after
 * an operator, so that you can synthesize it and get delay and area,
 * without the synthesis tools optimizing out your design because it
 * is connected to nothing.
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
#include "utils.hpp"
#include "Operator.hpp"
#include "TestBench.hpp"



TestBench::TestBench(Target* target, Operator *op, int n):
  Operator(target), op(op)
{
  unique_name="TestBench_" + op->unique_name;

  add_output ("done");    // true when test is finished
  add_output ("OK");      // true when one test is successful
  add_output ("passed");  // true when test is finished and all tests successful

  // declare internal registered signals
  for(int i=0; i<op->ioList.size(); i++){
    string idext = op->ioList[i]->id() ;
    add_signal(idext, op->ioList[i]->width());
  }

  // Build the test case list
  op->add_standard_test_cases(test_case_list);
  op->add_random_test_cases(test_case_list, n);
  
}



TestBench::~TestBench() {
}





/** builds a ticking clock */

void output_clock_def(std::ostream& o) {
}





void TestBench::output_vhdl(ostream& o, string name) {
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);
  set_sequential(); // a test bench does have a clock and rst. Actually that's all that it has.

  output_vhdl_entity(o);
  o << "architecture behavorial of " << name  << " is" << endl;

  // the operator to wrap
  op->Operator::output_vhdl_component(o);
  // The local signals
  output_vhdl_signal_declarations(o);

  o << "begin\n";

  // the instance
  o << tab << "test:" << op->unique_name << "\n"
    << tab << tab << "port map ( ";
  for(int i=0; i<op->ioList.size(); i++) {
    Operator::Signal s = *op->ioList[i];
    if(i>0) 
      o << tab << tab << "           ";
    string idext =  op->ioList[i]->id() ;
    if(op->ioList[i]->type() == Signal::in)
      o << op->ioList[i]->id()  << " =>  " << idext;
    else
      o << op->ioList[i]->id()  << " =>  " << idext;
    if (i < op->ioList.size()-1) 
      o << "," << endl;
  }
  if(op->is_sequential()) { // add the clock and reset
    o << "," << endl << tab << tab << "           clk =>  clk,"<<endl
      << tab << tab << "           rst =>  rst );"<<endl;
  }
  else 
    o << ");" <<endl;

  o <<endl;
  o << tab << "-- Ticking clock signal" <<endl;
  o << tab << "process" <<endl;
  o << tab << "begin" <<endl;
  o << tab << tab << "if not done then" <<endl;
  o << tab << tab << tab << "clk <= '1';" <<endl;
  o << tab << tab << tab << "wait for 5 ns;" <<endl;
  o << tab << tab << tab << "clk <= '0';" <<endl;
  o << tab << tab << tab << "wait for 5 ns;" <<endl;
  o << tab << tab << "else" <<endl;
  o << tab << tab << tab << "clk <= 'U';" <<endl;
  o << tab << tab << tab << "wait;" <<endl;
  o << tab << tab << "end if;" <<endl;
  o << tab << "end process;" <<endl;

  o <<endl;
  o << tab << "-- Setting the inputs" <<endl;
  o << tab << "process" <<endl;
  o << tab << "begin" <<endl;

  for (int i=0; i<test_case_list.size(); i++) {
    TestCaseInput in = test_case_list[i].input;
    TestCaseOutput out = test_case_list[i].expected_output;
    for (TestCaseInput::iterator it=in.begin(); it!=in.end(); it++) {
      string name = it->first;
      mpz_class val = it->second;
      Signal *s = get_signal_by_name(name);
      // assign this input
      o  << tab << tab << name << " <= \"" << unsigned_binary(val, s->width()) << "\"; -- " ;
      // add a comment with the value in human-readable form
      if(s->isFP()) 
	o << "FP number" << endl; //TODO convert
      else 
	o << val <<endl;
    }
    o  << tab << tab << "wait 10ns;" <<endl;
    
  }
  o << tab << "end process;" <<endl;
  o <<endl;


  o << tab << "-- Checking the outputs" <<endl;
  o << tab << "process" <<endl;
  o << tab << "begin" <<endl;
  o  << tab << tab << "wait "<< op->pipeline_depth()*10 <<"ns; -- wait for pipeline to flush" <<endl;

  for (int i=0; i<test_case_list.size(); i++) {
    TestCaseOutput out = test_case_list[i].expected_output;
    for (TestCaseOutput::iterator it=out.begin(); it!=out.end(); it++) {
      o << tab << tab << "Assert "<< it->first << "=" << it->second << "; -- " << test_case_list[i].comment << endl; 
      
    }
    o  << tab << tab << "wait 10ns;" <<endl;
    
  }
  o << tab << "end process;" <<endl;


  
  o << "end architecture;" << endl << endl;
    
}
