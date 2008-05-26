/*
 * Author: Cristian KLEIN 
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
#include <fstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "BigTestBench.hpp"

/** Rounds a number to the next multiple of 4 */
int inline round4(int x)
{
	return (x % 4) ? x + (4 - x % 4) : x;
}

BigTestBench::BigTestBench(Target* target, Operator* op, int n):
	Operator(target), op(op), n(n)
{
	unique_name = "BigTestBench_" + op->unique_name;
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

BigTestBench::~BigTestBench() { }

/**
 * Gets the name of the variable which is somehow related to a signal.
 * Also does escaping of special characters by replacing them with 'C'.
 * @param type prefix to add to the name of the signal
 * @param s the signal from which to take the name
 * @param var the zero-based expected output variant number, or -1 for input signal
 * @return a string which can be used as a VHDL identifier
 */
inline std::string s2vg(std::string type, Signal s, int var = -1)
{
	std::stringstream o;

	o << type << "_";
	char *x = strdup(s.id().c_str());
	for (char *p = x; *p; p++)
		switch (*p)
		{
			case '(':
			case ')':
			case ' ':
				*p = 'C';
				break;
		}
	o << x;
	free(x);

	if (var >= 0)
		o << "_var" << var;
	return o.str();
}

/**
 * Gets the name of the variable which should temporarly
 * store a value retrieved from the file, for a certain signal.
 * @param s the signal for which this value is retrieved.
 * @param var the zero-based expected output variant number, or -1 for input signal
 * @return a string which can be used as a VHDL identifier
 */
inline std::string s2v(Signal s, int var = -1)
{
	return s2vg("v", s, var);
}

/**
 * Like s2v, but gets the „good” signal.
 * @see s2v.
 */
inline std::string s2g(Signal s, int var = -1)
{
	return s2vg("g", s, var);
}

void BigTestBench::output_vhdl(ostream& o, string name) {
	/* Generate some TestCases
	 * We do this as early as possible, so that no VHDL is generated
	 * if test cases are not implemented for the given operator
	 */
	cerr << "Generating BIG test bench ... ";

	TestCaseList tclStd = op->generateStandardTestCases(1) +
		op->generateRandomTestCases(1000);

	/* Get the number of possible output values */
	TestCaseList::Inputs&  inputs  = tclStd.getInputMap();
	TestCaseList::Outputs& outputs = tclStd.getOutputMap();

	Licence(o,"Cristian KLEIN (2008)");
	o << "library ieee;" <<endl;
	o << "use ieee.std_logic_textio.all;" <<endl;
	o << "use ieee.std_logic_1164.all;" <<endl;
	o << "use std.textio.all;" <<endl;
	o << endl;
	o << "library work;" <<endl;

	output_vhdl_entity(o);
	o << "architecture behavorial of " << name  << " is" << endl;

	// the operator to wrap
	op->output_vhdl_component(o);
	// The local signals
	output_vhdl_signal_declarations(o);

	o << tab << "file fTest : text open read_mode is \""
		<< op->unique_name << ".test\";" <<endl;

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

	o << tab << "-- Input Process" << endl;
	o << tab << "process" <<endl;
	o << tab << tab << "variable buf : line;" <<endl;

	for (TestCaseList::Inputs::iterator it = inputs.begin(); it != inputs.end(); it++)
	{
		Signal s = it->first;
		o << tab << tab << "variable v_" << s.id() << ": std_ulogic"; 
		if (s.width() > 1)
			o << "_vector(" << round4(s.width())-1 << " downto 0)";
		o << ";" << endl;
	}

	for (TestCaseList::Outputs::iterator it = outputs.begin(); it != outputs.end(); it++)
	{
		Signal s = it->first;
		int maxValNum = it->second;
		for (int i = 0; i < maxValNum; i++)
		{
			o << tab << tab << "variable " << s2v(s,i) << ": std_ulogic"; 
			if (s.width() > 1)
				o << "_vector(" << round4(s.width())-1 << " downto 0)";
			o << ";" << endl;
			o << tab << tab << "variable " << s2g(s,i) << ": boolean;\n";
		}
	}

	o << tab << "begin" <<endl;
	o << tab << tab << "rst <= '1';" <<endl;
	o << tab << tab << "wait until falling_edge(clk);" <<endl;
	o << tab << tab << "rst <= '0';" <<endl;
	o <<endl;
	o << tab << tab << "while not (endfile(fTest)) loop" <<endl;
	o << tab << tab << tab << "-- wait for clk" << endl;
	o << tab << tab << tab << "wait until falling_edge(clk);" << endl;
	o << tab << tab << tab << "-- read inputs" << endl;
	for (TestCaseList::Inputs::iterator it = inputs.begin(); it != inputs.end(); it++)
	{
		Signal s = it->first;
		o << tab << tab << tab << "readline(fTest,buf); ";
		if (s.width() == 1)
			o << "read(buf," << s2v(s) << "); " << s.id() << " <= " << s2v(s) << ";" <<endl; 
		else
			o << "hread(buf," << s2v(s) << "); " << s.id() << " <= to_stdlogicvector(" << s2v(s) << "(" << s.width()-1 << " downto 0));" <<endl; 
	}
	o << tab << tab << tab << "-- check results" << endl;
	o << tab << tab << tab << "wait for 1 ns;" <<endl;
	for (TestCaseList::Outputs::iterator it = outputs.begin(); it != outputs.end(); it++)
	{
		Signal s = it->first;
		int maxValNum = it->second;

		for (int i = 0; i < maxValNum; i++)
		{
			o << tab << tab << tab
				<< "readline(fTest,buf); ";
			if (s.width() == 1)
				o << "read";
			else
				o << "hread";
			o << "(buf," << s2v(s,i) << "," << s2g(s,i) << ");\n";
		}

		o << tab << tab << tab << "assert not " << s2g(s,0) << "\n";
		for (int i = 0; i < maxValNum; i++)
		{	
			o << tab << tab << tab << tab << "or "
				<< "(" << s2g(s,i) << " and " << s.id() << "=";
			if (s.width() == 1)
				o << s2v(s,i);
			else 
				o << "to_stdlogicvector(" << s2v(s,i) << "(" << s.width()-1 << " downto 0))";
			o << ")" << endl;
		}
		o << tab << tab << tab << tab << "report \"Incorrect value for " << s.id() << "\" severity error;\n";
	}
 
	o << tab << tab << "end loop;" <<endl;
	o << tab << tab << "assert false report \"End of simulation\" severity failure;" << endl;
	o << tab << "end process;" <<endl;
	o << "end architecture;" <<endl;

	/* Generate the test file */
	TestCaseList tclIn, tclOut;
	ofstream f;
	f.open((op->unique_name + ".test").c_str(), ios::out);
	for (int i = 0; i < n; i++)
	{
		/* Generate a new test case list, if we depleted the one before */
		if (i % 10000 == 0)
		{
			tclIn = op->generateRandomTestCases(10000);
			cerr << ".";
		}
		if ((i - op->pipeline_depth()) % 10000 == 0)
			tclOut = tclIn;
		
		/* Input vector */
		TestCase tc = tclIn.getTestCase(i % 10000);
		for (TestCaseList::Inputs::iterator it = inputs.begin(); it != inputs.end(); it++)
		{
			Signal s = it->first;
			f << tc.signalValueToVHDLHex(s, tc.getInput(s), false) << endl;
		}

		/* Output vector */
		if (i >= op->pipeline_depth())
			tc = tclOut.getTestCase((i - op->pipeline_depth()) % 10000);
		for (TestCaseList::Outputs::iterator it = outputs.begin(); it != outputs.end(); it++)
		{
			Signal s = it->first;
			int maxValNum = it->second;

			for (int k = 0; k < maxValNum; k++)
			{
				try {
					if (i < op->pipeline_depth()) throw 0;
					f << tc.signalValueToVHDLHex(s, tc.getExpectedOutput(s, k), false) << endl;
				} catch (int) {
					f << "(pipeline fill)" << endl;
				} catch (std::string) {
					f << "N/A" << endl;
				}
			}
		}
	}
	f.close();
	cerr << endl << endl;

	cerr << "To run the simulation, type the following in 'rlwrap vsim -c':" <<endl;
	cerr << tab << "vdel -all -lib work" <<endl;
	cerr << tab << "vlib work" <<endl;
	cerr << tab << "vcom flopoco.vhdl" <<endl;
	cerr << tab << "vsim " << name <<endl;
	cerr << tab << "run -all" << endl;
}
