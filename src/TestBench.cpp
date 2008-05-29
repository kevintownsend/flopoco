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
#include <set>
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
	set_pipeline_depth(42);	// could be any number

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
	/*
	 * Generate some TestCases
	 * We do this as early as possible, so that no VHDL is generated
	 * if test cases are not implemented for the given operator
	 */
	TestIOMap::TestIOMapContainer tim = op->getTestIOMap().get();
	typedef TestIOMap::TestIOMapContainer::iterator TIMit;

	TestCaseList tcl;

	/* Get the size of the I/O Map */
	int size = 0;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		int maxNumValues = it->second;
		size += maxNumValues;
	}

	/* Allocate buffer which store I/O values metadata */
	bool *a_isIn = new bool[size];
	bool *a_isFP = new bool[size];
	int *a_w = new int[size];
	const Signal **a_s = new const Signal*[size];
	
	/* Get the metadata */
	size = 0;
	int max_bits = 0;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		const Signal& s = it->first;
		int maxNumValues = it->second;

		for (int j = 0; j < maxNumValues; j++)
		{
			a_s   [size] = &s;
			a_isIn[size] = (s.type() == Signal::in);
			a_isFP[size] = s.isFP();
			a_w   [size] = s.width();
			max_bits = max(max_bits, a_w[size]);
			size++;
		}
	}

	/* Generate test cases */
	mpz_class *a = new mpz_class[size];
	for (int i = 0; i < n; i++)
	{
		/* Fill inputs, erase outputs */
		for (int j = 0; j < size; j++)
		{
			if (a_isIn[j])
			{
				if (a_isFP[j])
					a[j] = (mpz_class(1) << (a_w[j]-2)) + getLargeRandom(a_w[j]-2);
				else
					a[j] = getLargeRandom(a_w[j]);
			}
			else
				a[j] = -1;
		}

		/* Get correct outputs */
		op->fillTestCase(a);

		/* Store test case */
		TestCase tc;
		for (int j = 0; j < size; j++)
		{
			if (a_isIn[j])
				tc.addInput(*a_s[j], a[j]);
			else
			{
				if (a[j] >= 0)
					tc.addExpectedOutput(*a_s[j], a[j]);
			}
		}
		tcl.add(tc);
	}
	delete[] a;
	delete[] a_s;
	delete[] a_isIn;
	delete[] a_isFP;
	delete[] a_w;

	Licence(o,"Cristian KLEIN (2007)");
	Operator::StdLibs(o);

	output_vhdl_entity(o);
	o << "architecture behavorial of " << name  << " is" << endl;

	// the operator to wrap
	op->output_vhdl_component(o);
	// The local signals
	output_vhdl_signal_declarations(o);

	o << endl <<
		tab << "-- FP compare function (found vs. real)\n" <<
		tab << "function fp_equal(a : std_logic_vector; b : std_logic_vector) return boolean is\n" <<
		tab << "begin\n" <<
		tab << tab << "if b(b'high downto b'high-1) = \"01\" then\n" <<
		tab << tab << tab << "return a = b;\n" <<
		tab << tab << "else\n" <<
		tab << tab << tab << "return a(a'high downto a'high-2) = b(b'high downto b'high-2);\n" <<
		tab << tab << "end if;\n" <<
		tab << "end;\n";

	/* In VHDL, literals may be incorrectly converted to „std_logic_vector(... to ...)” instead
	 * of „downto”. So, for each FP output width, create a subtype used for casting.
	 */
	{
		std::set<int> widths;
		for (TIMit it = tim.begin(); it != tim.end(); it++)
		{
			const Signal& s = it->first;
			if (s.type() != Signal::out) continue;
			if (s.isFP() != true) continue;
			widths.insert(s.width());
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
