/*
 * A test bench generator for FloPoCo. 
 *
 * Author : Cristian Klein, Florent de Dinechin
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
#include <set>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "TestBench.hpp"

extern int LongAccN;
TestBench::TestBench(Target* target, Operator* op, int n):
	Operator(target), op_(op), n_(n)
{
	LongAccN = n;
	setOperatorName();
	setPipelineDepth(42);	// could be any number

	// declare internal registered signals
	for(int i=0; i < op_->getIOListSize(); i++){
		string idext = op_->getIOListSignal(i)->getName() ;
		addSignal(idext, op_->getIOListSignal(i)->width());
	}

	/* add bogus clk and rst signal if the UUT does not have them */
	try {
		addSignal("clk", 1);
		addSignal("rst", 1);
	} catch (std::string) {
		/* silently ignore */
	}

	// Generate the standard and random test cases for this operator
	op-> buildStandardTestCases(&tcl_);
	op-> buildRandomTestCases(&tcl_, n);
}

TestBench::~TestBench() { 
}

void TestBench::setOperatorName(){
	uniqueName_ = "TestBench_" + op_->getOperatorName();
}

void TestBench::outputVHDL(ostream& o, string name) {
	licence(o,"Florent de Dinechin, Cristian Klein (2007)");
	Operator::stdLibs(o);

	outputVHDLEntity(o);
	o << "architecture behavorial of " << name  << " is" << endl;

	// the operator to wrap
	op_->outputVHDLComponent(o);
	// The local signals
	outputVHDLSignalDeclarations(o);

	o << endl <<  // Fixed by Bodgan
		tab << "-- FP compare function (found vs. real)\n" <<
		tab << "function fp_equal(a : std_logic_vector; b : std_logic_vector) return boolean is\n" <<
		tab << "begin\n" <<
		tab << tab << "if b(b'high downto b'high-1) = \"01\" then\n" <<
		tab << tab << tab << "return a = b;\n" <<
		tab << tab << "elsif b(b'high downto b'high-1) = \"11\" then\n" <<
		tab << tab << tab << "return (a(a'high downto a'high-1)=b(b'high downto b'high-1));\n" <<
		tab << tab << "else\n" <<
		tab << tab << tab << "return a(a'high downto a'high-2) = b(b'high downto b'high-2);\n" <<
		tab << tab << "end if;\n" <<
		tab << "end;\n";

	/* In VHDL, literals may be incorrectly converted to „std_logic_vector(... to ...)” instead
	 * of „downto”. So, for each FP output width, create a subtype used for casting.
	 */
	{
		std::set<int> widths;
		for (int i=0; i<op_->getIOListSize(); i++){
			Signal* s = op_->getIOListSignal(i);
			
			if (s->type() != Signal::out) continue;
			if (s->isFP() != true) continue;
			widths.insert(s->width());
		}

		if (widths.size() > 0)
			o << endl << tab << "-- FP subtypes for casting\n";
		for (std::set<int>::iterator it = widths.begin(); it != widths.end(); it++)
		{
			int w = *it;
			o << tab << "subtype fp" << w << " is std_logic_vector(" << w-1 << " downto 0);\n";
		}
	}
	
	o << "begin\n";

	// the instance
	o << tab << "uut:" << op_->getOperatorName() << "\n"
	  << tab << tab << "port map ( clk => clk, rst => rst," << endl;
	for(int i=0; i<op_->getIOListSize(); i++) {
		o << tab << tab << "           ";
		string idext =  op_->getIOListSignal(i)->getName() ;
		if(op_->getIOListSignal(i)->type() == Signal::in)
			o << op_->getIOListSignal(i)->getName()  << " =>  " << idext;
		else
			o << op_->getIOListSignal(i)->getName()  << " =>  " << idext;
		if (i < op_->getIOListSize()-1) 
			o << "," << endl;
	}
	o << ");" <<endl;

	o <<endl;

	o << tab << "-- Ticking clock signal" <<endl;
	o << tab << "process" <<endl;
	o << tab << "begin" <<endl;
	o << tab << tab << "clk <= '0';" <<endl;
	o << tab << tab << "wait for 5 ns;" <<endl;
	o << tab << tab << "clk <= '1';" <<endl;
	o << tab << tab << "wait for 5 ns;" <<endl;
	o << tab << "end process;" <<endl;
	o <<endl;

	o << tab << "-- Setting the inputs" <<endl;
	o << tab << "process" <<endl;
	o << tab << "begin" <<endl;
	o << tab << tab << "-- Send reset" <<endl;
	o << tab << tab << "rst <= '1';" << endl;
	o << tab << tab << "wait for 10 ns;" << endl;
	o << tab << tab << "rst <= '0';" << endl;
	for (int i = 0; i < tcl_.getNumberOfTestCases(); i++)
	{
		o << tcl_.getTestCase(i)->getInputVHDL(tab + tab);
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
	o << tab << tab << "wait for "<< op_->getPipelineDepth()*10 <<" ns; -- wait for pipeline to flush" <<endl;
	currentOutputTime += op_->getPipelineDepth()*10;
	for (int i = 0; i < tcl_.getNumberOfTestCases(); i++)
	{
		o << tab << tab << "wait for 5 ns;" <<endl;
		currentOutputTime += 5;
		o << tab << tab << "-- " << "current time: " << currentOutputTime <<endl;
		o << tcl_.getTestCase(i)->getInputVHDL(tab + tab + "-- input: ");
		o << tcl_.getTestCase(i)->getExpectedOutputVHDL(tab + tab);
		o << tab << tab << "wait for 5 ns;" <<endl;
		currentOutputTime += 5;
	} 
	o << tab << tab << "assert false report \"End of simulation\" severity failure;" <<endl;
	o << tab << "end process;" <<endl;
	
	o << "end architecture;" << endl << endl;

	cerr << "To run the simulation, type the following in 'vsim -c':" <<endl;
	cerr << tab << "vdel -all -lib work" <<endl;
	cerr << tab << "vlib work" <<endl;
	cerr << tab << "vcom flopoco.vhdl" <<endl;
	cerr << tab << "vsim " << name <<endl;
	cerr << tab << "add wave -r *" <<endl;
	cerr << tab << "run " << currentOutputTime << endl;
}
