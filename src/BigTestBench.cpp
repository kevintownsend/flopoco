/* vim: set tabstop=8 softtabstop=2 shiftwidth=2: */
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

void BigTestBench::output_vhdl(ostream& o, string name) {
	// XXX: Not actually used
	/* Generate some TestCases
	 * We do this as early as possible, so that no VHDL is generated
	 * if test cases are not implemented for the given operator
	 */
	op->generateStandardTestCases(n);

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
	o << tab << tab << "variable good : boolean;" <<endl <<endl;
	for (int i = 0; i < op->ioList.size(); i++)
	{
		Signal *s = op->ioList[i];
		o << tab << tab << "variable v_" << s->id() << ": std_ulogic"; 
		if (s->width() > 1)
			o << "_vector(" << s->width() - 1 << " downto 0)";
		o << ";" << endl;
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
	for (int i = 0; i < op->ioList.size(); i++)
	{
		Signal *s = op->ioList[i];
		if (s->type() != Signal::in) continue;
		if ((s->id() == "clk") || (s->id() == "rst")) continue;
		o << tab << tab << tab << "readline(fTest,buf); ";
		if (s->width() == 1)
			o << "read(buf,v_" << s->id() << "); " << s->id() << " <= v_" << s->id() << ";" <<endl; 
		else
			o << "hread(buf,v_" << s->id() << "); " << s->id() << " <= to_stdlogicvector(v_" << s->id() << ");" <<endl; 
	}
	o << tab << tab << tab << "-- check results" << endl;
	o << tab << tab << tab << "wait for 1 ns;" <<endl;
	for (int i = 0; i < op->ioList.size(); i++)
	{
		Signal *s = op->ioList[i];
		if (s->type() != Signal::out) continue;
		o << tab << tab << tab
			<< "readline(fTest,buf); ";
		if (s->width() == 1)
			o << "read";
		else
			o << "hread";

		o << "(buf,v_" << s->id() << ", good); ";
		o << "assert not good or " << s->id() << " = ";
		if (s->width() == 1)
			o << "v_" << s->id() << " ";
		else 
			o << "to_stdlogicvector(v_" << s->id() << ") ";
		o << "report \"Incorrect value for " << s->id() << "\" severity error;"
			<<endl; 
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
			tclIn = op->generateRandomTestCases(10000);
		if ((i - op->pipeline_depth()) % 10000 == 0)
			tclOut = tclIn;
		
		/* Input vector */
		TestCase tc = tclIn.getTestCase(i % 10000);
		for (int j = 0; j < op->ioList.size(); j++)
		{
			Signal& s = *op->ioList[j];
			if (s.type() != Signal::in) continue;
			if ((s.id() == "clk") || (s.id() == "rst")) continue;
			f << tc.signalValueToVHDLHex(s, tc.getInput(s), false) << endl;
		}

		/* Output vector */
		if (i >= op->pipeline_depth())
			tc = tclOut.getTestCase((i - op->pipeline_depth()) % 10000);
		for (int j = 0; j < op->ioList.size(); j++)
		{
			Signal& s = *op->ioList[j];
			if (s.type() != Signal::out) continue;
			try {
				if (i < op->pipeline_depth()) throw std::string();
				// Happens when no value is expected from this output signal
				f << tc.signalValueToVHDLHex(s, tc.getOneExpectedOutput(s), false) << endl;
			} catch (std::string) {
				f << "(none expected)" << endl;
			}
		}
	}
	f.close();

	cerr << "To run the simulation, type the following in 'rlwrap vsim -c':" <<endl;
	cerr << tab << "vdel -all -lib work" <<endl;
	cerr << tab << "vlib work" <<endl;
	cerr << tab << "vcom flopoco.vhdl" <<endl;
	cerr << tab << "vsim " << name <<endl;
	cerr << tab << "run -all" << endl;
}
