// works only with sollya
#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "FixSinCos.hpp"

#include "FixedPointFunctions/FunctionEvaluator.hpp"

using namespace std;


// personalized parameter
//string FixSinCos::operatorInfo = "FixSinCos w <options>";



FixSinCos::FixSinCos(Target * target, int w_):Operator(target), w(w_)
{
/* constructor of the FixSinCos
  Target is the targeted FPGA : Stratix, Virtex ... (see Target.hpp for more informations)
  w and param1 are some parameters declared by this Operator developpers, 
  any number can be declared, you will have to modify 
      -> this function,  
      -> the prototype of this function (available in FixSinCos.hpp)
      -> the lines in main.cpp where the command line arguments are parsed in order to generate this FixSinCos
*/
	/* In this constructor we are going to generate an operator that takes as input three bit vectors X,Y,Z of lenght w, treats them as unsigned integers, sums them and then output the last param1 bit of the sum adding the first bit of the sum (most significant) in front of this output, all the vhdl code needed by this operator has to be generated in this function */

	// definition of the name of the operator
	ostringstream name;
	name << "FixSinCos_" << w	//<< "_" << param1
	    ;
	setName(name.str());

	setCopyrightString("Guillaume Sergent 2012");

	/* SET UP THE IO SIGNALS
	   Each IO signal is declared by addInput(name,n) or addOutput(name,n) 
	   where name is a string that stands for the name of the variable and 
	   n is an integer (int)   that stands for the length of the corresponding 
	   input/output */

	// declaring inputs
	addInput("X", w);
	//addFullComment(" addFullComment for a large comment ");
	//addComment("addComment for small left-aligned comment");

	// declaring output
	addOutput("S", w);
	addOutput("C", w);


	/* Some peace of informations can be delivered to the flopoco user if  the -verbose option is set
	   [eg: flopoco -verbose=0 FixSinCos 10 5 ]
	   , by using the REPORT function.
	   There is three level of details
	   -> INFO for basic information ( -verbose=1 )
	   -> DETAILED for complete information, includes the INFO level ( -verbose=2 )
	   -> DEBUG for complete and debug information, to be used for getting information 
	   during the debug step of operator development, includes the INFO and DETAILED levels ( -verbose=3 )
	 */
	// basic message
	//REPORT(INFO,"Declaration of FixSinCos \n");

	// more detailed message
	//ostringstream detailedMsg;
	//detailedMsg  << "this operator has received two parameters " << w << " and " << param1;
	//REPORT(DETAILED,detailedMsg.str());

	// debug message for developper
	//REPORT(DEBUG,"debug of FixSinCos");



	/* vhdl is the stream which receives all the vhdl code, some special functions are
	   available to smooth variable declaration and use ... 
	   -> when first using a variable (Eg: T), declare("T",64) will create a vhdl
	   definition of the variable : signal T and includes it it the header of the architecture definition of the operator

	   Each code transmited to vhdl will be parsed and the variables previously declared in a previous cycle will be delayed automatically by a pipelined register.
	 */
	//vhdl << declare("T",
	//		w + 1) << " <= ('0' & X) + ('O' & Y);" << endl;

	// the argument is reduced into (0,1/4) because one sin/cos
	// computation in this range can always compute the right sin/cos
	vhdl << declare ("A",1) << " <= X" << of (w-1) << ";" << endl
	     << declare ("B",1) << " <= X" << of (w-2) << ";" << endl
	     << declare ("Y",w-2) << " <= X " << range (w-3,0) << ";"
	     << endl;
	// now X -> A*.5 + B*.25 + Y where A,B \in {0,1} and Y \in {0,.25}

	// Y_prime = .25 - y
	// the VHDL is probably invalid
	vhdl << declare ("Y_prime",w-2) << " <= "
	     << "2**" << (w-2) << " - Y;" << endl;

	// Y_in = / Y_prime if B=1
	//        \ Y if B=0
	vhdl << declare ("Y_in",w-2) << " <= Y_prime when B='1' else Y;"
	     << endl;

	// Exch = A ^ B
	vhdl << declare ("Exch",1) << " <= A xor B;" << endl;

	FunctionEvaluator *f_sin, *f_cos;
	// 5 is hardcoded: when addition formulae are added it will
	// be replaced by dichotomy
	f_sin = new FunctionEvaluator (target,"sin(x*Pi*4)",w-2,w,5);
	f_cos = new FunctionEvaluator (target,"cos(x*Pi*4)",w-2,w,5);

	oplist.push_back (f_sin);
	oplist.push_back (f_cos);
	inPortMap (f_sin, "X", "Y_prime");
	outPortMap (f_sin, "R", "S_out");
	inPortMap (f_cos, "X", "Y_prime");
	outPortMap (f_cos, "R", "C_out");

	// now just add the signs to the signals
	/* the output format is consequently
	 * struct output_format {
	 *	unsigned int sgn: 1; // sign is 1 for neg, as in IEEE754
	 *	unsigned int mantissa: w;
	 * };
	 */
	vhdl << "S <= '0' & (C_out when Exch else S_out);" << endl;
	vhdl << "C <= A & (S_out when Exch else C_out);" << endl;


	/* declaring a new cycle, each variable used after this line will be delayed 
	   with respect to his state in the precedent cycle
	 */
	//nextCycle();
	//vhdl << declare("R",
	//		w + 2) << " <=  ('0' & T) + (\"00\" & Z);" << endl;

	/* the use(variable) is a deprecated function, that can still be encoutered 
	   in old Flopoco's operators, simply use vhdl << "S" (flopoco will generate the correct delayed value as soon as a previous declare("S") exists */
	// we first put the most significant bit of the result into R
	//vhdl << "S <= (R" << of(w + 1) << " & ";
	// and then we place the last param1 bits
	//vhdl << "R" << range(param1 - 1, 0) << ");" << endl;
};


void FixSinCos::emulate(TestCase * tc)
{
	/* This function will be used when the TestBench command is used in the command line
	   we have to provide a complete and correct emulation of the operator, in order to compare correct output generated by this function with the test input generated by the vhdl code */
	/* first we are going to format the entries */
	mpz_class sx = tc->getInputValue("X");
	//appends 010 on left of the binary representation (yeak)
	mpz_class fsx ("010" + sx.get_str(2),2);
	FPNumber fpx(0, w);
	fpx = fsx;

	mpfr_t x, sind, cosd; //round-to-nearest as of now
	mpfr_init2 (x, 1+w);
	mpfr_init2 (sind, 1+w);
	mpfr_init2 (cosd, 1+w);
	fpx.getMPFR (x);
	mpfr_sin (sind, x, GMP_RNDD); //round to lower
	mpfr_cos (cosd, x, GMP_RNDD);
	FPNumber fpsin (0, w, sind);
	FPNumber fpcos (0, w, cosd);
	mpz_class ssin = fpsin.getSignalValue();
	mpz_class scos = fpcos.getSignalValue();
	// now we have to remove the 2 first bits of S/C
	string ssin_s = ssin.get_str(2);
	string scos_s = scos.get_str(2);
	mpz_class ssin_fix_d (ssin_s.substr (2, string::npos), 2);
	mpz_class scos_fix_d (scos_s.substr (2, string::npos), 2);
	// and since FunctionEvaluator does only faithful rounding
	// we add also as expected results the upper roundings
	tc->addExpectedOutput ("S", ssin_fix_d);
	tc->addExpectedOutput ("C", scos_fix_d);
	tc->addExpectedOutput ("S", ssin_fix_d+1);
	tc->addExpectedOutput ("C", scos_fix_d+1);
	mpfr_clears (x, sind, cosd, NULL);
}


void FixSinCos::buildStandardTestCases(TestCaseList * tcl)
{
}

void FixSinCos::buildRandomTestCases(TestCaseList * tcl, int n)
{
}

TestCase *FixSinCos::buildRandomTestCases(int i)
{
	TestCase *tc = new TestCase(this);
	return tc;
}

#endif // SOLLYA

