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

namespace flopoco{


	extern int LongAccN;
	TestBench::TestBench(Target* target, Operator* op, int n):
		Operator(target), op_(op), n_(n)
	{
		LongAccN = n;
		setName("TestBench_" + op_->getName());
		setPipelineDepth(42);	// could be any number


		// Generate the standard and random test cases for this operator
		op-> buildStandardTestCases(&tcl_);
		op-> buildRandomTestCases(&tcl_, n);
	

		// The instance
		//  portmap inputs and outputs
		string idext;
		for(int i=0; i < op->getIOListSize(); i++){
			Signal* s = op->getIOListSignal(i);
			if(s->type() == Signal::out) 
				outPortMap (op, s->getName(), s->getName());
			if(s->type() == Signal::in) {
				declare(s->getName(), s->width());
				inPortMap (op, s->getName(), s->getName());
			}
		}
		// add clk and rst
		declare("clk");
		declare("rst");


		// The VHDL for the instance
		vhdl << instance(op, "test");

		vhdl << tab << "-- Ticking clock signal" <<endl;
		vhdl << tab << "process" <<endl;
		vhdl << tab << "begin" <<endl;
		vhdl << tab << tab << "clk <= '0';" <<endl;
		vhdl << tab << tab << "wait for 5 ns;" <<endl;
		vhdl << tab << tab << "clk <= '1';" <<endl;
		vhdl << tab << tab << "wait for 5 ns;" <<endl;
		vhdl << tab << "end process;" <<endl;
		vhdl <<endl;

		vhdl << tab << "-- Setting the inputs" <<endl;
		vhdl << tab << "process" <<endl;
		vhdl << tab << "begin" <<endl;
		vhdl << tab << tab << "-- Send reset" <<endl;
		vhdl << tab << tab << "rst <= '1';" << endl;
		vhdl << tab << tab << "wait for 10 ns;" << endl;
		vhdl << tab << tab << "rst <= '0';" << endl;
		for (int i = 0; i < tcl_.getNumberOfTestCases(); i++)
			{
				vhdl << tcl_.getTestCase(i)->getInputVHDL(tab + tab);
				vhdl << tab << tab << "wait for 10 ns;" <<endl;
			} 
		vhdl << tab << tab << "wait for 100000 ns; -- allow simulation to finish" << endl;
		vhdl << tab << "end process;" <<endl;
		vhdl <<endl;

		int currentOutputTime = 0;
		vhdl << tab << "-- Checking the outputs" <<endl;
		vhdl << tab << "process" <<endl;
		vhdl << tab << "begin" <<endl;
		vhdl << tab << tab << "wait for 10 ns; -- wait for reset to complete" <<endl;
		currentOutputTime += 10;
		if (op_->getPipelineDepth() > 0){
			vhdl << tab << tab << "wait for "<< op_->getPipelineDepth()*10 <<" ns; -- wait for pipeline to flush" <<endl;
			currentOutputTime += op_->getPipelineDepth()*10;
		}
		else{
			vhdl << tab << tab << "wait for "<< 2 <<" ns; -- wait for pipeline to flush" <<endl;
			currentOutputTime += 2;
		}
	
		for (int i = 0; i < tcl_.getNumberOfTestCases(); i++)
			{
				//		vhdl << tab << tab << "wait for 5 ns;" <<endl;
				//		currentOutputTime += 5;
				vhdl << tab << tab << "-- " << "current time: " << currentOutputTime <<endl;
				TestCase* tc = tcl_.getTestCase(i);
				if (tc->getComment() != "")
					vhdl << tab <<  "-- " << tc->getComment() << endl;
				vhdl << tc->getInputVHDL(tab + tab + "-- input: ");
				vhdl << tc->getExpectedOutputVHDL(tab + tab);
				vhdl << tab << tab << "wait for 10 ns;" <<endl;
				currentOutputTime += 10;
			} 
		vhdl << tab << tab << "assert false report \"End of simulation\" severity failure;" <<endl;
		vhdl << tab << "end process;" <<endl;
	
		simulationTime=currentOutputTime;
	}

	TestBench::~TestBench() { 
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
		o << vhdl.str() << endl;
		o << "end architecture;" << endl << endl;


	}
}
