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

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <gmpxx.h>
#include <cstdlib>

#include "Operator.hpp"
#include "TestIOMap.hpp"
#include "BigTestBench.hpp"

mpz_class getLargeRandom(int w);

/** Rounds a number to the next multiple of 4 */
int inline round4(int x)
{
	return (x % 4) ? x + (4 - x % 4) : x;
}

BigTestBench::BigTestBench(Target* target, Operator* op, int n):
	Operator(target), op_(op), n_(n)
{
	setOperatorName();
	setPipelineDepth(42);	// could be any number

	// declare internal registered signals
	for(int i=0; i<op_->getIOListSize(); i++){
		string idext = op_->getIOListSignal(i)->getSignalName() ;
		addSignal(idext, op_->getIOListSignal(i)->width());
	}

	/* add bogus clk and rst signal if the UUT does not have them */
	try {
		addSignal("clk", 1);
		addSignal("rst", 1);
	} catch (std::string) {
		/* silently ignore */
	}
}

BigTestBench::~BigTestBench() { }

void BigTestBench::setOperatorName(){
	uniqueName_ = "BigTestBench_" + op_->getOperatorName();
}

/**
 * Gets the name of the variable which is somehow related to a signal.
 * Also does escaping of special characters by replacing them with 'C'.
 * @param type prefix to add to the name of the signal
 * @param s the signal from which to take the name
 * @param var the zero-based expected output variant number, or -1 for input signal
 * @return a string which can be used as a VHDL identifier
 */
inline std::string s2vg(std::string type, const Signal& s, int var = -1)
{
	std::stringstream o;

	o << type << "_";
	char *x = strdup(s.getSignalName().c_str());
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
inline std::string s2v(const Signal& s, int var = -1)
{
	return s2vg("v", s, var);
}

/**
 * Like s2v, but gets the „good” signal.
 * @see s2v.
 */
inline std::string s2g(const Signal& s, int var = -1)
{
	return s2vg("g", s, var);
}

void BigTestBench::outputVHDL(ostream& o, string name) {
	size_t res;
	cerr << "Generating BIG test bench ... ";

	/* Get IO Map */
	TestIOMap::TestIOMapContainer tim = op_->getTestIOMap().get();
	typedef TestIOMap::TestIOMapContainer::iterator TIMit;

	licence(o,"Cristian KLEIN (2008)");
	o << "library ieee;" <<endl;
	o << "use ieee.std_logic_textio.all;" <<endl;
	o << "use ieee.std_logic_1164.all;" <<endl;
	o << "use std.textio.all;" <<endl;
	o << endl;
	o << "library work;" <<endl;

	outputVHDLEntity(o);
	o << "architecture behavorial of " << name  << " is" << endl;

	// the operator to wrap
	op_->outputVHDLComponent(o);
	// The local signals
	outputVHDLSignalDeclarations(o);

	o << tab << "file fTest : text open read_mode is \""
		<< op_->getOperatorName() << ".test\";" <<endl;
	
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

	o << "begin\n";

	// the instance
	// XXX: Ugly! Will write something to encapsulate this.
	o << tab << "uut:" << op_->getOperatorName() << "\n"
		<< tab << tab << "port map ( ";
	for(int i=0; i<op_->getIOListSize(); i++) {
		//Signal& s = *op_->getIOListSignal(i);
		if(i>0) 
			o << tab << tab << "           ";
		string idext =  op_->getIOListSignal(i)->getSignalName() ;
		if(op_->getIOListSignal(i)->type() == Signal::in)
			o << op_->getIOListSignal(i)->getSignalName()  << " =>  " << idext;
		else
			o << op_->getIOListSignal(i)->getSignalName()  << " =>  " << idext;
		if (i < op_->getIOListSize()-1) 
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

	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		const Signal& s = it->first;
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
	o << tab << tab << tab << "-- read everything from file" << endl;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		const Signal& s = it->first;
		int numMaxValues = it->second;

		for (int i = 0; i < numMaxValues; i++)
		{
			o << tab << tab << tab << "readline(fTest,buf); ";
			if (s.width() == 1)
				o << "read(buf," << s2v(s,i) << "," << s2g(s,i) << ");\n";
			else
				o << "hread(buf," << s2v(s,i) << "," << s2g(s,i) << ");\n"; 
		}
	}
	
	o << tab << tab << tab << "-- assign inputs" << endl;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		const Signal& s = it->first;
		int numMaxValues = it->second;
		if (s.type() != Signal::in) continue;

		/* XXX: useless, but uniform */
		for (int i = 0; i < numMaxValues; i++)
		{
			o << tab << tab << tab;
			if (s.width() == 1)
				o << s.getSignalName() << " <= " << s2v(s,i) << "; ";
			else
				o << s.getSignalName() << " <= to_stdlogicvector(" << s2v(s,i) << "(" << s.width()-1 << " downto 0)); "; 
			o << "assert " << s2g(s,i) << " report \"Invalid input for " << s.getSignalName() << "\" severity failure;\n";
		}
	}

	o << tab << tab << tab << "-- wait a bit\n";
	o << tab << tab << tab << "wait for 1 ns;\n";
 	o << tab << tab << tab << "-- check outputs" << endl;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		const Signal& s = it->first;
		int numMaxValues = it->second;
		if (s.type() != Signal::out) continue;

		o << tab << tab << tab << "assert (not " << s2g(s,0) << ") -- we don't care at all about this signal or\n";
		for (int i = 0; i < numMaxValues; i++)
		{
			o << tab << tab << tab << tab << "or (" << s2g(s,i) << " and "; 
			if (s.isFP())
			{
				o << "fp_equal(" << s.getSignalName() << ", to_stdlogicvector(" << s2v(s,i) << "(" << s.width()-1 << " downto 0)))";
			}
			else
			{
				if (s.width() == 1)
					o << s.getSignalName() << " = " << s2v(s,i);
				else
					o << s.getSignalName() << " = to_stdlogicvector(" << s2v(s,i) << "(" << s.width()-1 << " downto 0))"; 
			}
			o << ")\n";
		}
		o << tab << tab << tab << tab << "report \"Incorrect value for " << s.getSignalName() << "\";\n";
	}
 
	o << tab << tab << "end loop;" <<endl;
	o << tab << tab << "assert false report \"End of simulation\" severity failure;" << endl;
	o << tab << "end process;" <<endl;
	o << "end architecture;" <<endl;

	/* Generate the test file */
	FILE *f = fopen((op_->getOperatorName() + ".test").c_str(), "w");

	/* Get number of TestCase I/O values */
	int size = 0;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
		size += it->second;

	/* WARNING: Highly optimized (i.e. hard to read, hard to not make mistakes)
	 * code ahead */

	/* Allocate buffer which store I/O values metadata */
	bool *a_isIn = new bool[size];
	bool *a_isFP = new bool[size];
	int *a_w = new int[size];
	
	/* Get the metadata */
	size = 0;
	int max_bits = 0;
	for (TIMit it = tim.begin(); it != tim.end(); it++)
	{
		const Signal& s = it->first;
		int maxNumValues = it->second;

		for (int j = 0; j < maxNumValues; j++)
		{
			a_isIn[size] = (s.type() == Signal::in);
			a_isFP[size] = s.isFP();
			a_w   [size] = s.width();
			max_bits = max(max_bits, a_w[size]);
			size++;
		}
	}

	/* Buffer for zero-copy padding output values
	 * size is max_bits * 2 * 2 + 1 + 1
	 *              |     |   |   |   +-- end-of string zero
	 *              |     |   |   +------ newline
	 *              |     |   +---------- number of characters per bytes
	 *              |     |               (two because of hex representation)
	 *              |     +-------------- once for number, once for possible padding
	 *              +-------------------- maximum number of bits in inputs and outputs
	 */
	char *buf = new char[4*max_bits+2];
	/* Pre-fill the left half, used only for padding */
	memset(buf, '0', 2*max_bits);

	/* TestVector list w/ pipeline
	 * implemented as a circular buffer
	 * new value go to head */
	int head = 0;
	int depth = op_->getPipelineDepth();

	mpz_class** a;
	a = new mpz_class*[depth+2];
	for (int i = 0; i < depth + 1; i++) {
		a[i] = new mpz_class[size];
		for (int j = 0; j < size; j++)
			a[i][j] = -1;
	}
	a[depth+1] = a[0];

	/* Output test I/O values */
	for (int i = 0; i < n_; i++)
	{
		/* Fill inputs, erase outputs */
		for (int j = 0; j < size; j++)
		{
			if (a_isIn[j])
			{
				if (a_isFP[j])
					a[head][j] = (mpz_class(1) << (a_w[j]-2)) + getLargeRandom(a_w[j]-2);
				else
					a[head][j] = getLargeRandom(a_w[j]);
			}
			else
				a[head][j] = -1;
		}

		/* Get correct outputs */
		op_->fillTestCase(a[head]);

		/* Write it to the file */
		for (int j = 0; j < size; j++)
		{
			/* Load the value as the current input
			 * or the output delayed by the pipeline depth */
			mpz_class& z = a_isIn[j] ? a[head][j] : a[head+1][j];

			/* Should this value be outputed */
			if (z < 0)
			{
				res=fwrite("N/A\n", 4, 1, f);
				if(res!=1) {
					cerr << "BigTestBench: Error in fwrite (disk full?)\n";
					exit(1);
				}
				continue;
			}

			/* Output the hex number from the middle of the buffer */
			char *p = &buf[2*max_bits];
			/* Get hex number and size*/
			mpz_get_str(p, 16, z.get_mpz_t());
			int l = strlen(p);
			/* Compute number of zeros to pad and
			 * adjust pointer and length */
			int toadd = (a_w[j]-1)/4+1-l;
			p -= toadd; l += toadd;
			/* Add newline */
			p[l] = '\n'; p[l+1] = 0; l++;
			/* Write to file */
			res=fwrite(p, l, 1, f);
			if(res!=1) {
				cerr << "BigTestBench: Error in fwrite (disk full?)\n";
				exit(1);
			}
		}

		/* Advance pipeline circular buffer */
		head++;
		if (head == depth + 1)
			head = 0;
	}
	
	/* Cleanup */
	delete buf;
	delete a_isIn;
	delete a_isFP;
	delete a_w;
	for (int i = 0; i < depth + 1; i++)
		delete[] a[i];
	delete a;
	fclose(f);

	cerr << endl << endl;

	cerr << "To run the simulation, type the following in 'vsim -c':" <<endl;
	cerr << tab << "vdel -all -lib work" <<endl;
	cerr << tab << "vlib work" <<endl;
	cerr << tab << "vcom flopoco.vhdl" <<endl;
	cerr << tab << "vsim " << name <<endl;
	cerr << tab << "run -all" << endl;
}

