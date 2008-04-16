/* vim: set tabstop=8 softtabstop=2 shiftwidth=2: */
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

TestBench::TestBench(Target* target, Operator* op, int n):
  Operator(target), op(op), n(n)
{
  unique_name = "TestBench_" + op->unique_name;

  // declare internal registered signals
  for(int i=0; i<op->ioList.size(); i++){
    string idext = op->ioList[i]->id() ;
    add_signal(idext, op->ioList[i]->width());
  }

  /* add bogus clk and rst signal if the UUT does not have them */
  try {
    add_signal("clk", 1);
    add_signal("rst", 1);
  } catch (std::string) {
    /* silently ignore */
  }
}

TestBench::~TestBench() { }

void TestBench::output_vhdl(ostream& o, string name) {
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);

  output_vhdl_entity(o);
  o << "architecture behavorial of " << name  << " is" << endl;

  // the operator to wrap
  op->output_vhdl_component(o);
  // The local signals
  output_vhdl_signal_declarations(o);

  o << "begin\n";

  // the instance
  // XXX: Ugly! Will write something to encapsulate this.
  o << tab << "uut:" << op->unique_name << "\n"
    << tab << tab << "port map ( ";
  for(int i=0; i<op->ioList.size(); i++) {
    Signal s = *op->ioList[i];
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
  o << ");" <<endl;

  o <<endl;

  o << tab << "-- Ticking clock signal" <<endl;
  o << tab << "process" <<endl;
  o << tab << "begin" <<endl;
  o << tab << tab << "clk <= '1';" <<endl;
  o << tab << tab << "wait for 5 ns;" <<endl;
  o << tab << tab << "clk <= '0';" <<endl;
  o << tab << tab << "wait for 5 ns;" <<endl;
  o << tab << "end process;" <<endl;
  o <<endl;

  /* Generate some TestCases */
  TestCaseList tcl = op->generateStandardTestCases(n) + op->generateRandomTestCases(n);

  o << tab << "-- Setting the inputs" <<endl;
  o << tab << "process" <<endl;
  o << tab << "begin" <<endl;
  o << tab << tab << "-- Send reset" <<endl;
  o << tab << tab << "rst <= '1';" << endl;
  o << tab << tab << "wait for 10 ns;" << endl;
  o << tab << tab << "rst <= '0';" << endl;
  for (int i = 0; i < tcl.getNumberOfTestCases(); i++)
  {
    o << tcl.getTestCase(i).getInputVHDL(tab + tab);
    o << tab << tab << "wait for 10 ns;" <<endl;
  } 
  o << tab << tab << "wait for 100000 ns; -- allow simulation to finish" << endl;
  o << tab << "end process;" <<endl;
  o <<endl;

  int currentOutputTime = 0;
  o << tab << "-- Checking the outputs" <<endl;
  o << tab << "process" <<endl;
  o << tab << "begin" <<endl;
  o << tab << tab << "wait for 10 ns; -- wait for reset to complete" <<endl;
  currentOutputTime += 10;
  o << tab << tab << "wait for "<< op->pipeline_depth()*10 <<" ns; -- wait for pipeline to flush" <<endl;
  currentOutputTime += op->pipeline_depth()*10;
  for (int i = 0; i < tcl.getNumberOfTestCases(); i++)
  {
    o << tab << tab << "wait for 5 ns;" <<endl;
    currentOutputTime += 5;
    o << tab << tab << "-- " << "current time: " << currentOutputTime <<endl;
    o << tcl.getTestCase(i).getInputVHDL(tab + tab + "-- input: ");
    o << tcl.getTestCase(i).getExpectedOutputVHDL(tab + tab);
    o << tab << tab << "wait for 5 ns;" <<endl;
    currentOutputTime += 5;
  } 
  o << tab << tab << "assert false report \"End of simulation\" severity failure;" <<endl;
  o << tab << "end process;" <<endl;
  
  o << "end architecture;" << endl << endl;

  cerr << "To run the simulation, type the following in 'rlwrap vsim -c':" <<endl;
  cerr << tab << "vdel -all -lib work" <<endl;
  cerr << tab << "vlib work" <<endl;
  cerr << tab << "vcom flopoco.vhdl" <<endl;
  cerr << tab << "vsim " << name <<endl;
  cerr << tab << "add wave -r *" <<endl;
  cerr << tab << "run " << currentOutputTime << endl;
}
