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
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "TestBench.hpp"

using namespace std;

namespace flopoco{


	extern int LongAccN;
	TestBench::TestBench(Target* target, Operator* op, int n, bool fromFile):
		Operator(target), op_(op), n_(n)
	{
		LongAccN = n;
		setName("TestBench_" + op_->getName());
		
		setCombinatorial(); // this is a combinatorial operator
		setPipelineDepth(0);	// could be any number
		setCycle(0);


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
                // declaration and assignation of the recirculation signal
                if (op_->isRecirculatory()) {
                  declare("stall_s");
                  vhdl << tab << "stall_s <= '0';" << endl;
                }
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
		vhdl << endl;


                if (fromFile) generateTestFromFile();
                else generateTestInVhdl();
        }








        /* Generating the tests using a file to store the IO, allow to have a lot of IOs without
         * increasing the VHDL compilation time
         */ 
        void TestBench::generateTestFromFile() {
                // we reordonate the Signal in order to put all the output 
                // TODO :could be clean by using two list, directly retrieved from the operator
                vector<Signal*> inputSignalVector;
                vector<Signal*> outputSignalVector;
		for(int i=0; i < op_->getIOListSize(); i++){
			Signal* s = op_->getIOListSignal(i);
                        if (s->type() == Signal::out) outputSignalVector.push_back(s);
                        else if (s->type() == Signal::in) inputSignalVector.push_back(s);
                };

                // decleration of test time
		int currentOutputTime = 0;

                // In order to generate the file containing inputs and expected output in a correct order
                // we will store the use order for file decompression
                list<string> IOorderInput;
                list<string> IOorderOutput;


                vhdl << tab << "-- Reading the input from a file " << endl;
		vhdl << tab << "process" <<endl;


                /* Variable declaration */
                vhdl << tab << tab << "variable inline : line; " << endl;                    // variable to read a line
                vhdl << tab << tab << "variable counter : integer := 1;" << endl;
                vhdl << tab << tab << "variable tmpChar : character;" << endl;                        // variable to store a character (escape between inputs)
                vhdl << tab << tab << "file inputsFile : text is \"essai.input\"; " << endl; // declaration of the input file


                /* Variable to store value for inputs and expected outputs*/
		for(int i=0; i < op_->getIOListSize(); i++){
			Signal* s = op_->getIOListSignal(i);
			vhdl << tab << tab << "variable V_" << s->getName(); 
                                      vhdl << " : bit_vector("<< s->width() - 1 << " downto 0);" << endl;
		}

                /* Process Beginning */
                vhdl << tab << "begin" << endl;

                /* Reset Sending */
		vhdl << tab << tab << "-- Send reset" <<endl;
		vhdl << tab << tab << "rst <= '1';" << endl;
		vhdl << tab << tab << "wait for 10 ns;" << endl;
		vhdl << tab << tab << "rst <= '0';" << endl;
                currentOutputTime += 10;

                /* File Reading */
                vhdl << tab << tab << "while not endfile(inputsFile) loop" << endl;
                vhdl << tab << tab << tab << " -- positionning inputs" << endl;

                /* All inputs and the corresponding expected outputs will be on the same line
                 * so we begin by reading this line, once and for all (once by test) */
                vhdl << tab << tab << tab << "readline(inputsFile,inline);" << endl;

                // input reading and forwarding to the operator
		for(unsigned int i=0; i < inputSignalVector.size(); i++){
			Signal* s = inputSignalVector[i];
                        vhdl << tab << tab << tab << "read(inline ,V_"<< s->getName() << ");" << endl;
                        vhdl << tab << tab << tab << "read(inline,tmpChar);" << endl; // we consume the character between each inputs
			if (s->width() > 1) vhdl << tab << tab << tab << s->getName() << " <= " << "to_stdlogicvector(V_" << s->getName() << ");" << endl;
			else vhdl << tab << tab << tab << s->getName() << " <= " << "to_stdlogicvector(V_" << s->getName() << ")(0);" << endl;
                        // adding the IO to IOorder
                        IOorderInput.push_back(s->getName());
		}

                /**
                 * Declaration of a waiting time between sending the input and comparing
                 * the result with the output
                 * in case of a pipelined operator we have to wait the complete latency of all the operator
                 * that means all the pipeline stages each step
                 * TODO : entrelaced the inputs / outputs in order to avoid this wait
                 */
                vhdl << tab << tab << tab << " -- verifying the corresponding output" << endl;
                vhdl << tab << tab << tab << " wait for 5 ns;" << endl;
                currentOutputTime += 5 * tcl_.getNumberOfTestCases();
		vhdl << tab << tab << tab << "wait for "<< op_->getPipelineDepth()*10 <<" ns; -- wait for pipeline to flush" <<endl;
		currentOutputTime += op_->getPipelineDepth()*10* tcl_.getNumberOfTestCases();

                // consume the "? "
                vhdl << tab << tab << tab << "read(inline,tmpChar);" << endl; // we consume the character between each inputs
                vhdl << tab << tab << tab << "read(inline,tmpChar);" << endl; // we consume the character between each inputs



		for(unsigned int i=0; i < outputSignalVector.size(); i++){
			Signal* s = outputSignalVector[i];
                        vhdl << tab << tab << tab << "read(inline ,V_"<< s->getName() << ");" << endl;
                        vhdl << tab << tab << tab << "read(inline,tmpChar);" << endl; // we consume the character between each inputs
                        if (s->isFP()) //cerr << "managing fp equality is not implemented yet, see line ~ 197 TestBench.cpp  . "<< endl;
                        vhdl << tab << tab << tab << "assert false or fp_equal(fp"<< s->width() << "'(" << s->getName() << ") ,to_stdlogicvector(V_" <<  s->getName() << ")) report(\"Incorrect output for R, expected \" & str(to_stdlogicvector(V_R)) & \" and it outputs \" & str(R)) &  \"|| line : \" & integer'image(counter) & \" of input file \" ;"<< endl;
			else if (s->isIEEE()) vhdl << tab << tab << tab << "assert false or fp_equal_ieee(" << s->getName() << " ,to_stdlogicvector(V_" <<  s->getName() << "),"<<s->wE()<<" , "<<s->wF()<<") report(\"Incorrect output for R, expected \" & str(to_stdlogicvector(V_R)) & \" and it outputs \" & str(R)) &  \"|| line : \" & integer'image(counter) & \" of input file \" ;"<< endl;
                        else vhdl << tab << tab << tab << "assert false or (" << s->getName() << "= to_stdlogicvector(V_" << s->getName() << "))" << "report(\"Incorrect output for R, expected \" & str(to_stdlogicvector(V_" << s->getName() << ")) & \" and it outputs \" & str(" << s->getName() <<")) &  \"|| line : \" & integer'image(counter) & \" of input file \" ;"<< endl;
                        //else cerr << " Le test à partir d'un fichier n'est pas encore implémenté pour les vecteurs non IEEE" << endl;

                        /* adding the IO to the IOorder list */
                        IOorderOutput.push_back(s->getName());
                };
                vhdl << tab << tab << tab << " wait for 5 ns; -- wait for pipeline to flush" << endl;
                vhdl << tab << tab << tab << "counter := counter + 1;" << endl;
                currentOutputTime += 5 * tcl_.getNumberOfTestCases();
                vhdl << tab << tab << "end loop;" << endl;

		vhdl << tab << tab << "assert false report \"End of simulation\" severity failure;" <<endl;
		vhdl << tab << "end process;" <<endl;

                /* Setting the computed simulation Time */	
		simulationTime=currentOutputTime;

	        /* Generating a file of inputs */ 
                // opening a file to write down the output (for text-file based test)
                // if n < 0 we do not generate a file
                if (n_ >= 0) {
                  string inputFileName = "essai.input";
                  ofstream fileOut(inputFileName.c_str(),ios::out);
                  // if error at opening, let's mention it !
                  if (!fileOut) cerr << "FloPoCo was not abe to open " << inputFileName << " in order to write down inputs. " << endl;
	        	for (int i = 0; i < tcl_.getNumberOfTestCases(); i++)	{
		        	TestCase* tc = tcl_.getTestCase(i);
                                if (fileOut) fileOut << tc->generateInputString(IOorderInput,IOorderOutput);
                      } 

                      // closing input file
                      fileOut.close();
                };
          }













          void TestBench::generateTestInVhdl() {
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
				vhdl << tab << tab << "-- current time: " << currentOutputTime <<endl;
				TestCase* tc = tcl_.getTestCase(i);
				if (tc->getComment() != "")
					vhdl << tab <<  "-- " << tc->getComment() << endl;
				vhdl << tc->getInputVHDL(tab + tab + "-- input: ");
                                // if a file is open for output we generated them
                                //if (fileOut) fileOut << tc->generateInputString();
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


                o << endl << endl << endl;
                /* Generation of Vhdl function to parse file into std_logic_vector */


                o << " -- converts std_logic into a character" << endl;
                o << tab << "function chr(sl: std_logic) return character is" << endl
                  << tab << tab << "variable c: character;" << endl
                  << tab << "begin" << endl 
                  << tab << tab << "case sl is" << endl 
                  << tab << tab << tab << "when 'U' => c:= 'U';" << endl
                  << tab << tab << tab << "when 'X' => c:= 'X';" << endl  
                  << tab << tab << tab << "when '0' => c:= '0';" << endl  
                  << tab << tab << tab << "when '1' => c:= '1';" << endl
                  << tab << tab << tab << "when 'Z' => c:= 'Z';" << endl  
                  << tab << tab << tab << "when 'W' => c:= 'W';" << endl  
                  << tab << tab << tab << "when 'L' => c:= 'L';" << endl  
                  << tab << tab << tab << "when 'H' => c:= 'H';" << endl  
                  << tab << tab << tab << "when '-' => c:= '-';" << endl  
                  << tab << tab << "end case;" << endl 
                  << tab << tab << "return c;" << endl 
                  << tab <<  "end chr;" << endl; 


                o << tab << "-- converts std_logic into a string (1 to 1)" << endl  
                  << tab << "function str(sl: std_logic) return string is" << endl  
                  << tab << " variable s: string(1 to 1);" << endl  
                  << tab << " begin" << endl  
                  << tab << tab << "s(1) := chr(sl);" << endl  
                  << tab << tab << "return s;" << endl 
                  << tab << "end str;" << endl;  
     


                o << tab << "-- converts std_logic_vector into a string (binary base)" << endl  
                  << tab << "-- (this also takes care of the fact that the range of" << endl  
                  << tab << "--  a string is natural while a std_logic_vector may" << endl  
                  << tab << "--  have an integer range)" << endl  
                  << tab << "function str(slv: std_logic_vector) return string is" << endl  
                  << tab << tab << "variable result : string (1 to slv'length);" << endl  
                  << tab << tab << "variable r : integer;" << endl  
                  << tab << "begin" << endl  
                  << tab << tab << "r := 1;" << endl  
                  << tab << tab << "for i in slv'range loop" << endl  
                  << tab << tab << tab << "result(r) := chr(slv(i));" << endl  
                  << tab << tab << tab << "r := r + 1;" << endl  
                  << tab << tab << "end loop;" << endl  
                  << tab << tab << "return result;" << endl  
                  << tab << "end str;" << endl; 

              o << endl << endl << endl;

                /* If op_ is an IEEE operator (IEEE input and output, we define) the function
                 * fp_equal for the considered precision in the ieee case
                 */
		o << endl <<  // Fixed by Nicolas
			tab << "-- test isZero\n" <<
			tab << "function iszero(a : std_logic_vector) return boolean is\n" <<
			tab << "begin\n" <<
			tab << tab << "return  a = (a'high downto 0 => '0');\n" <<        // test if exponent = "0000---000"
			tab << "end;\n\n\n" <<
                        
			tab << "-- FP IEEE compare function (found vs. real)\n" <<
			tab << "function fp_equal_ieee" << "(a : std_logic_vector;" 
                                                        << " b : std_logic_vector;" 
                                                        << " we : integer;"  
                                                        << " wf : integer) return boolean is\n" 
		     <<	tab << "begin\n" <<
			tab << tab << "if a(wf+we downto wf) = b(wf+we downto 0) and b(we+wf-1 downto wf) = (we downto 1 => '1') then\n" <<        // test if exponent = "1111---111"
			tab << tab << tab << "if iszero(b(wf-1 downto 0)) then return  iszero(a(wf-1 downto 0));\n" <<               // +/- infinity cases
                        tab << tab << tab << "else return not iszero(a(wf - 1 downto 0));\n" <<         
                        tab << tab << tab << "end if;\n" <<         
			tab << tab << "else\n" <<
			tab << tab << tab << "return a = b;\n" <<
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
