/*
  the FloPoCo command-line interface
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <mpfr.h>
#include <cstdlib>


#include "Operator.hpp"
#include "FlopocoStream.hpp"

#include "Target.hpp"
#include "Targets/Spartan3.hpp"
#include "Targets/Virtex4.hpp"
#include "Targets/Virtex5.hpp"
#include "Targets/Virtex6.hpp"
#include "Targets/StratixII.hpp"
#include "Targets/StratixIV.hpp"
#include "Shifters.hpp"
#include "LZOC.hpp"
#include "LZOCShifterSticky.hpp"
#include "IntAdder.hpp"
#include "IntComparator.hpp"
#include "IntNAdder.hpp"
#include "IntCompressorTree.hpp"
#include "LongIntAdder.hpp"
#include "LongIntAdderCmpCmpAdd.hpp"
#include "LongIntAdderCmpAddAdd.hpp"

#include "IntDualSub.hpp"

#include "IntMultiplier.hpp"
#include "IntTilingMult.hpp"
#include "Targets/DSP.hpp"

#include "SignedIntMultiplier.hpp"
#include "IntTruncMultiplier.hpp"
#include "IntKaratsuba.hpp"

#include "FPMultiplier.hpp"
#include "FPMultiplierKaratsuba.hpp"
#include "FPMultiplierTiling.hpp"

#include "FPSquarer.hpp"
#include "FPAdder.hpp"
#include "FPAdderSinglePath.hpp"
#include "FPAdder3Input.hpp"

#include "Fix2FP.hpp"
//~ #include "apps/CoilInductance/CoordinatesTableX.hpp"
//~ #include "apps/CoilInductance/CoordinatesTableZ.hpp"
//~ #include "apps/CoilInductance/CoordinatesTableY.hpp"
//#include "apps/CoilInductance/CoilInductance.hpp"
#include "FPDiv.hpp"
#include "FPSqrt.hpp"
#include "FPSqrtPoly.hpp"

#include "LongAcc.hpp"
#include "LongAcc2FP.hpp"
#include "DotProduct.hpp"


#include "PolynomialEvaluator.hpp"


#include "Wrapper.hpp"
#include "TestBench.hpp"

#include "ConstMult/IntConstMult.hpp"
#include "ConstMult/FPConstMult.hpp"
#include "ConstMult/IntIntKCM.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "ConstMult/FPRealKCM.hpp"


#include "FPExp.hpp" 
#include "FPLog.hpp"
#include "FPPow.hpp"

#include "InputIEEE.hpp"
#include "OutputIEEE.hpp"
 
#include "apps/Collision.hpp"
#include "apps/FPFMAcc.hpp"
#include "apps/FPJacobi.hpp"


#include "IntSquarer.hpp"

#ifndef _WIN32
#include "ConstMult/CRFPConstMult.hpp"
#include "ConstMult/FPConstMultParser.hpp"


#ifdef HAVE_HOTBM
#include "HOTBM.hpp"
#endif

#ifdef HAVE_LNS
#include "LNS/LNSAddSub.hpp"
#include "LNS/LNSAdd.hpp"
#include "LNS/CotranTables.hpp"
#include "LNS/Cotran.hpp"
#include "LNS/CotranHybrid.hpp"
#include "LNS/LNSMul.hpp"
#include "LNS/LNSDiv.hpp"
#include "LNS/LNSSqrt.hpp"
#include "LNS/AtanPow.hpp"

#endif
#endif

#include "PolyTableGenerator.hpp"
#include "FunctionEvaluator.hpp"
#include "FPPipeline.hpp"

#include "UserDefinedOperator.hpp"

#define BRIGHT 1
#define RED 31
#define OPER 32
#define PARAM 34
#define OP(op,paramList)             {cerr << "    "; printf("%c[%d;%dm",27,1,OPER); cerr <<  op; printf("%c[%dm",27,0); cerr<< " "; printf("%c[%d;%dm",27,1,PARAM); cerr << paramList; printf("%c[%dm\n",27,0); } 



using namespace std;
using namespace flopoco;

// Global variables, useful through most of FloPoCo, to be encapuslated in something someday

namespace flopoco{

	vector<Operator*> oplist;

	string filename="flopoco.vhdl";
	string cl_name=""; // used for the -name option

	int verbose=0;
	Target* target;
	int LongAccN;

	/* flex vars */
	int yyTheCycle;
	vector<pair<string, int> > theUseTable;
	
	map<string, double> emptyDelayMap;
	bool combinatorialOperator;

	int Operator::uid = 0; //init of the uid static member of Operator

}

static void usage(char *name){
	cerr << "\nUsage: "<<name<<" <operator specification list>\n" ;
	cerr << "Each operator specification is one of: \n";
    OP( "UserDefinedOperator","param0 param1");
	cerr << "    ____________ SHIFTERS/LZOC _________________________________________________\n";
	OP("LeftShifter","wIn MaxShift");
	OP("RightShifter","wIn MaxShift");
	OP("LZOC","wIn");
	OP("LZOCShifter","wIn wOut");
	OP("LZCShifter","wIn wOut");
	OP("LOCShifter","wIn wOut");
	OP("LZOCShifterSticky","wIn wOut");
	OP("LZCShifterSticky","wIn wOut");
	OP("LOCShifterSticky","wIn wOut");
	cerr << "    ____________ ADDERS/SUBTRACTERS ____________________________________________\n";
	OP("IntAdder","wIn");
	cerr << "      Integer adder, possibly pipelined\n";
	OP("MyIntAdder","wIn optimizeType srl implementation bufferedInputs");
	cerr << "      Integer adder, multple parameters, possibly pipelined\n";
	cerr << "      optimizeType=<0,1,2,3> 0=LUT 1=REG 2=SLICE 3=LATENCY\n";
	cerr << "      srl=<0,1> Allow SRLs\n";
	cerr << "      implementation=<-1,0,1,2> -1=optimizeType dependent,\n";  
	cerr << "                                 0=Classical, 1=Alternative, 2=Short-Latency\n";
	cerr << "      bufferedInputs=<0,1>\n";
	OP("IntDualSub","wIn opType");
	cerr << "      Integer adder/subtracter or dual subtracter, possibly pipelined\n";
	cerr << "      opType: if 1, compute X-Y and X+Y; if 0, compute X-Y and Y-X \n";
	OP("IntNAdder","wIn N");
	cerr << "      Multi-operand addition, possibly pipelined\n";
	OP("IntCompressorTree","wIn N");
	cerr << "      Multi-operand addition using compressor trees, possibly pipelined\n";
	cerr << "    ____________ INTEGER MULTIPLIERS/SQUARER/KARATSUBA _________________________\n";
	OP("IntMultiplier","wInX wInY signed");
	cerr << "      Integer multiplier of two integers X and Y of sizes wInX and wInY \n";
	cerr << "      signed is one of {0,1} \n";	
	OP("SignedIntMultiplier","wInX wInY");
	cerr << "      Signed integer multiplier. wInX and wInY include the sign\n";	
	OP("IntTilingMultiplier","wInX wInY ratio maxTimeInMinutes");
	cerr << "      Integer multiplier of two integers X and Y of sizes wInX and wInY\n";
	cerr << "      0 <= ratio <= 1; larger ratio => DSP dominant architectures\n";
	cerr << "      maxTimeInMinutes 0..; 0=find optimal solution (no time limit)\n"; 	
	OP("IntTruncMultiplier","wInX wInY ratio error useLimits maxTimeInMinutes");
	cerr << "      Integer multiplier of two integers X and Y of sizes wInX and wInY \n"; 
	cerr << "      with a given error order.\n";	
	cerr << "      0 <= ratio <= 1; larger ratio => DSP dominant architectures\n";
	cerr << "      wInX+wInY<error<=0. The order of the error.\n";
	cerr << "      useLimits. Soft-core multipliers are size-limited\n";
	cerr << "      maxTimeInMinutes 0..; 0=find optimal solution (no time limit)\n"; 	
	OP ("IntKaratsuba","wIn");
	cerr << "      integer multiplier of two integers X and Y of sizes wIn. 17 < wIn <= 68\n";	
	OP ("IntSquarer","wIn");
	cerr << "      integer squarer. For now wIn <=68 \n";		
	OP( "IntConstMult","w c");
	cerr << "      Integer constant multiplier using shift-and-add: w - input size, c - the constant\n";
	OP( "IntIntKCM","w c signedInput");
	cerr << "      Integer constant multiplier using KCM: w - input size, c - the constant\n";
#ifdef HAVE_SOLLYA
	OP( "FixRealKCM","lsbIn msbIn signedInput lsbOut constant");
	cerr << "      Faithful multiplier of a fixed-point input by a real constant\n";
	cerr << "      The constant is provided as a Sollya expression, e.g \"log(2)\"\n";
#endif // HAVE_SOLLYA
	cerr << "    ____________ FLOATING-POINT OPERATORS ______________________________________\n";
	OP("Fix2FP","LSB MSB Signed wE wF");
	cerr << "      Convert a 2's compliment fixed-point number in the bit range MSB...LSB \n";
	cerr << "      into floating-point\n";
	OP( "FPAdder","wE wF");
	cerr << "      Floating-point adder \n";
	OP( "FPMultiplier","wE wF_in wF_out");
	cerr << "      Floating-point multiplier, supporting different in/out precision  \n";
	OP( "FPMultiplierKaratsuba","wE wF_in wF_out");
	cerr << "      Floating-point multiplier, supporting different in/out precision. \n";
	cerr << "      Mantissa multiplier uses Karatsuba\n";
	OP( "FPMultiplierTiling","wE wF_in wF_out ratio timeInMinutes");
	cerr << "      Floating-point multiplier, supporting different in/out precision. \n";
	cerr << "      Mantissa multiplier uses Tiling Algorithm  \n";
	OP( "FPSquarer","wE wFin wFout");
	cerr << "      Floating-point squarer \n";
	OP( "FPDiv","wE wF");
	cerr << "      Floating-point divider \n";
	OP("FPSqrt","wE wF");
	cerr << "      Floating-point square root, implemented using digit recurrence\n";
	cerr << "      (no DSP, long latency)\n";
#ifdef HAVE_SOLLYA
//	cerr << "    FPSqrtPoly wE wF correctlyRounded degree\n";
	OP( "FPSqrtPoly","wE wF degree");
	cerr << "      Floating-point square root, using polynomial approximation \n";
	cerr << "      (DSP-based, shorter latency and higher frequency for large wF)\n";
	cerr << "      faithful rounding\n";
//	cerr << "      correctlyRounded (0 or 1) selects between faithful and correct rounding (NYImplemented)\n";
//	cerr << "      correctlyRounded (0) selects faithful and correct rounding (NYImplemented)\n";
	cerr << "      degree (1,...k) polynomial degree. Higher degree => more DSP less BRAM\n";
#endif // HAVE_SOLLYA
	OP("FPConstMult","wE_in wF_in wE_out wF_out cst_sgn cst_exp cst_int_sig");
	cerr << "      Floating-point constant multiplier\n";
	cerr << "      The constant is provided as integral significand and integral exponent.\n";
#ifdef HAVE_SOLLYA
	OP( "FPConstMultParser","wE_in wF_in wE_out wF_out wF_C constant_expr");
	cerr << "      Floating-point constant multiplier with a parser for the constant:\n";
	cerr << "      last argument is a Sollya expression between double quotes,e.g.\"exp(pi/2)\".\n";
	OP( "CRFPConstMult","wE_in wF_in wE_out wF_out constant_expr");
	cerr << "      Correctly-rounded floating-point constant multiplier\n";
	cerr << "      The constant is provided as a Sollya expression, between double quotes.\n";
#endif // HAVE_SOLLYA
	OP( "LongAcc","wE_in wF_in MaxMSB_in LSB_acc MSB_acc");
	cerr << "      Long fixed-point accumulator\n";
	OP( "LongAcc2FP","LSB_acc MSB_acc wE_out wF_out");
	cerr << "      Post-normalisation unit for LongAcc \n";
	OP( "DotProduct","wE wFX wFY MaxMSB_in LSB_acc MSB_acc");
	cerr << "      Floating-point dot product unit \n";
	OP( "FPExp","wE wF");
	cerr << "      Floating-point exponential function\n";
	OP( "FPLog","wE wF InTableSize");
	cerr << "      Floating-point logarithm function;\n";
	cerr << "      InTableSize is the numbers of bits to input to the tables. \n";
	cerr << "      O defaults to something sensible\n";
	OP( "FPPowr","wE wF");// LogTableSize ExpTableSize ExpDegree");
	cerr << "      Floating-point powr function from IEEE-754-2008 (experimental);\n";
//	cerr << "      For the parameters, try 8 23 10 10 2 3 3 (simple), 11 52 12 12 2 33  (double)\n";
	OP( "OutputIEEE","wEI wFI wEO wFO");
	cerr << "      Conversion from FloPoCo to IEEE-754-like floating-point formats\n";
	OP( "InputIEEE","wEI wFI wEO wFO");
	cerr << "      Conversion from IEEE-754-like to FloPoCo floating-point formats\n";
#ifdef HAVE_HOTBM
	cerr << "    ____________ GENERIC FUNCTION EVALUATORS ____________________________________\n";
	cerr << "      We provide two methods to evaluate a function on [0,1]\n";
	OP( "FunctionEvaluator","function wI wO degree");
	cerr << "      Polynomial based method for fixed-point functions\n";
	cerr << "      wI - input width, wO - the number of bits at the right of the dot, \n";
	cerr << "      degree - degree of polynomial approx\n";
	cerr << "      function - sollya-syntaxed function to implement, between double quotes\n";
	OP( "HOTBM","function wI wO degree");
	cerr << "      High-Order Table-Based Method for fixed-point functions (NPY)\n";
	cerr << "      wI - input width, wO - output width, degree - degree of polynomial approx\n";
	cerr << "      function - sollya-syntaxed function to implement, between double quotes\n";
	OP( "HOTBMFX","function wE_in wF_in wE_out wF_out degree");
	cerr << "      Same as HOTBM, with explicit fixed-point formats (NPY)\n";
	cerr << "      Note: input is unsigned, output is signed.\n";
	OP( "HOTBMRange","function wI wO degree xmin xmax scale");
	cerr << "      Same as HOTBM, with explicit range and scale (NPY)\n";
	cerr << "      xmin xmax - bounds of the input range, mapped to [0,1[\n";
	cerr << "      scale - scaling factor to apply to the function output\n";
#endif // HAVE_HOTBM
#ifdef HAVE_LNS
	cerr << "    ____________ LNS OPERATORS _________________________________________________\n";
	cerr << "    Common parameters:\n";
	cerr << "      wE - width of integral part of exponent. Typically from 4 to 8.\n";
	cerr << "           Negative values allowed at your own risk.\n";
	cerr << "      wF - width of fractional part of exponent. Typically from 8 to 20.\n";
	OP( "LNSAddSub","wE wF");
	cerr << "      Addition in Logarithmic Number System.\n";
	OP( "LNSMul","wE wF");
	cerr << "      LNS multiplication.\n";
	OP( "LNSDiv","wE wF");
	cerr << "      LNS division.\n";
	OP( "LNSSqrt","wE wF");
	cerr << "      LNS square root.\n";
#if 0 // Currently broken, most of its code was commented out for using deprecated methods. TODO
	cerr << "    AtanPow wE wF o\n";
	cerr << "      (4/pi)*atan(2^x) function.\n";
#endif
#endif // HAVE_LNS
	cerr << "    ____________ APPLICATIONS __________________________________________________\n";
	OP("Collision","wE wF opt");;
	cerr << "       A collision detection operator, computes the predicate X²+Y²+Z²<R2\n";
	cerr << "       opt: assemble FP operators if 0, optimized architecture if 1 \n";

	// cerr << "  Applications: \n";
	// cerr << "    CoilInductance LSBI MSBI wEIn wFIn MaxMSBO LSBO MSBO FilePath\n";
	// cerr << "       \n";
	//To be removed from command line interface
	// cerr << "    CoordinatesTableX wIn LSB MSB FilePath\n";
	// cerr << "    CoordinatesTableY wIn LSB MSB FilePath\n";
	// cerr << "    CoordinatesTableZ wIn LSB MSB FilePath\n";
	//=====================================================

	cerr << "    ____________ TEST-BENCH ____________________________________________________\n";
	OP ("TestBench","n");
	cerr << "       Behavorial test bench for the preceding operator\n";
	cerr << "       This test bench will include standard tests, plus n random tests.\n";
	OP( "TestBenchFile","n");
	cerr << "       Behavorial test bench for the preceding operator\n";
	cerr << "       This test bench will include standard tests, plus n random tests.\n";
	cerr << "       Inputs and outputs are stored in a file to reduce VHDL compilation time.\n";
	cerr << "       if n=-2, an exhaustive test is generated (use only for small operators).\n";
	cerr << "    ____________ WRAPPER _______________________________________________________\n";
	OP ("Wrapper","");
	cerr << "       Wraps the preceding operator between registers\n";
	cerr << "(NPY) Not pipelined yet\n";
	cerr << "________________________________________________________________________________\n";
	cerr << "________________ OPTIONS________________________________________________________\n";
	cerr << "General options, affecting the operators that follow them:\n";
	cerr << "   -outputfile=<output file name>           (default=flopoco.vhdl)\n";
	cerr << "   -verbose=<1|2|3>                         (default=0)\n";
	cerr << "   -pipeline=<yes|no>                       (default=yes)\n";
	cerr << "   -frequency=<target frequency in MHz>     (default=400)\n";
	cerr << "   -target=<Spartan3|Virtex4|Virtex5|StratixII|StratixIV>      (default=Virtex4)\n";
	cerr << "   -DSP_blocks=<yes|no>\n";
	cerr << "       optimize for the use of DSP blocks   (default=yes)\n";
	cerr << "   -name=<entity name>\n";
	cerr << "       defines the name of the VHDL entity of the next operator\n";
	exit (EXIT_FAILURE);
}

int checkStrictyPositive(char* s, char* cmd) {
	int n=atoi(s);
	if (n<=0){
		cerr<<"ERROR: got "<<s<<", expected strictly positive number."<<endl;
		usage(cmd);
	}
	return n;
}

int checkPositiveOrNull(char* s, char* cmd) {
	int n=atoi(s);
	if (n<0){
		cerr<<"ERROR: got "<<s<<", expected positive-or-null number."<<endl;
		usage(cmd);
	}
	return n;
}

bool checkBoolean(char* s, char* cmd) {
	int n=atoi(s);
	if (n!=0 && n!=1) {
		cerr<<"ERROR: got "<<s<<", expected a boolean (0 or 1)."<<endl;
		usage(cmd);
	}
	return (n==1);
}


int checkSign(char* s, char* cmd) {
	int n=atoi(s);
	if (n!=0 && n!=1) {
		cerr<<"ERROR: got "<<s<<", expected a sign bit (0 or 1)."<<endl;
		usage(cmd);
	}
	return n;
}


void addOperator(Operator *op) {
	if(cl_name!="")	{
		op->changeName(cl_name);
		cl_name="";
	}
	// TODO check name not already in list...
	// TODO this procedure should be a static method of operator, so that other operators can call it
	oplist.push_back(op);
}



bool parseCommandLine(int argc, char* argv[]){
	if (argc<2) {
		usage(argv[0]); 
		return false;
	}

	Operator* op;
	int i=1;
	cl_name="";
	do {
		string opname(argv[i++]);
		if(opname[0]=='-'){
			string::size_type p = opname.find('=');
			if (p == string::npos) {
				cerr << "ERROR: Option missing an = : "<<opname<<endl; 
				return false;
			} else {
				string o = opname.substr(1, p - 1), v = opname.substr(p + 1);
				if (o == "outputfile") {
					if(oplist.size()!=0)
						cerr << "WARNING: global option "<<o<<" should come before any operator specification" << endl; 
					filename=v;
				}
				else if (o == "verbose") {
					verbose = atoi(v.c_str()); // there must be a more direct method of string
					if (verbose<0 || verbose>4) {
						cerr<<"ERROR: verbose should be 1, 2 or 3,    got "<<v<<"."<<endl;
						usage(argv[0]);
					}
				}
				else if (o == "target") {
					if(v=="Virtex4") target=new Virtex4();
					else if (v=="Virtex6") target=new Virtex6();
					else if (v=="Virtex5") target=new Virtex5();
					else if (v=="Spartan3") target=new Spartan3();
					else if (v=="StratixII") target=new StratixII();
					else if (v=="StratixIV") target=new StratixIV();
					else {
						cerr<<"ERROR: unknown target: "<<v<<endl;
						usage(argv[0]);
					}
				}
				else if (o == "pipeline") {
					if(v=="yes") target->setPipelined();
					else if(v=="no")  target->setNotPipelined();
					else {
						cerr<<"ERROR: pipeline option should be yes or no,    got "<<v<<"."<<endl; 
						usage(argv[0]);
					}
				}
				else if (o == "frequency") {
					int freq = atoi(v.c_str());
					if (freq>1 && freq<10000) {
						target->setFrequency(1e6*(double)freq);
						if(verbose) 
							cerr << "Frequency set to "<<target->frequency()<< " Hz" <<endl; 
					}
					else {
						cerr<<"WARNING: frequency out of reasonible range, ignoring it."<<endl; 
					}
				}
				else if (o == "DSP_blocks") {
					if(v=="yes") target->setUseHardMultipliers(true);
					else if(v=="no")  target->setUseHardMultipliers(false);
					else {
						cerr<<"ERROR: DSP_blocks option should be yes or no,    got "<<v<<"."<<endl; 
						usage(argv[0]);
					}
				}
				else if (o == "name") {
					cl_name=v; // TODO?  check it is a valid VHDL entity name 
				}
				else { 	
					cerr << "Unknown option "<<o<<endl; 
					return false;
				}
			}
		}
		else if(opname=="IntIntKCM"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int w = atoi(argv[i++]);
				mpz_class mpc(argv[i++]);
				int signedInput = checkBoolean(argv[i++], argv[0]);

				cerr << "> IntIntKCM , w="<<w<<", c="<<mpz2string(mpc) << (signedInput==0?" unsigned":" signed") << "\n";
				op = new IntIntKCM(target, w, mpc, signedInput);
				addOperator(op);
			}        
		}
		else if(opname=="IntConstMult"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int w = atoi(argv[i++]);
				mpz_class mpc(argv[i++]);
				cerr << "> IntConstMult , w="<<w<<", c="<<mpz2string(mpc)<<"\n";
				op = new IntConstMult(target, w, mpc);
				addOperator(op);
			}        
		}

		else if(opname=="FPConstMult"){
			int nargs = 7;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE_in = checkStrictyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictyPositive(argv[i++], argv[0]);
				int cst_sgn  = checkSign(argv[i++], argv[0]); 
				int cst_exp  = atoi(argv[i++]); // TODO no check on this arg
				mpz_class cst_sig(argv[i++]);
				cerr << "> FPConstMult, wE_in="<<wE_in<<", wF_in="<<wF_in
						 <<", wE_out="<<wE_out<<", wF_out="<<wF_out
						 <<", cst_sgn="<<cst_sgn<<", cst_exp="<<cst_exp<< ", cst_sig="<<mpz2string(cst_sig)<<endl;
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, cst_sgn, cst_exp, cst_sig);
				addOperator(op);
			}        
		} 	
#ifdef HAVE_SOLLYA
		else if(opname=="FixRealKCM"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int lsbIn = atoi(argv[i++]);
				int msbIn = atoi(argv[i++]);
				int signedInput = checkBoolean(argv[i++], argv[0]);
				int lsbOut = atoi(argv[i++]);
				string constant = argv[i++];
				cerr << "> FixRealKCM, msbIn="<<msbIn<<", lsbIn="<<lsbIn 
					  <<", lsbOut="<<lsbOut
					  << ", constant="<<constant <<endl;
				op = new FixRealKCM(target, lsbIn, msbIn, signedInput, lsbOut, constant);
				addOperator(op);
			}        
		}
		else if(opname=="FPRealKCM"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = atoi(argv[i++]);
				int wF = atoi(argv[i++]);
				string constant = argv[i++];
				cerr << "> FPRealKCM, wE="<<wE<<", wF="<<wF << ", constant="<<constant <<endl;
				op = new FPRealKCM(target, wE, wF, constant);
				addOperator(op);
			}        
		}
		else if(opname=="CRFPConstMult"){ 
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]);
			else { 
				int wE_in = checkStrictyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictyPositive(argv[i++], argv[0]);
				string constant = argv[i++];
				cerr << "> CRFPConstMult, wE_in="<<wE_in<<", wF_in="<<wF_in 
					  <<", wE_out="<<wE_out<<", wF_out="<<wF_out
					  << ", constant="<<constant <<endl;
				op = new CRFPConstMult(target, wE_in, wF_in, wE_out, wF_out, constant);
				addOperator(op);
			}        
		} 	
		else if(opname=="FPConstMultParser"){ 
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0]);
			else { 
				int wE_in = checkStrictyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictyPositive(argv[i++], argv[0]);
				int wF_C = checkStrictyPositive(argv[i++], argv[0]);
				string constant = argv[i++];
				cerr << "> FPConstMultParser, wE_in="<<wE_in<<", wF_in="<<wF_in 
					  <<", wE_out="<<wE_out<<", wF_out="<<wF_out 
					  << ", wF_C=" << wF_C
					  << ", constant="<<constant <<endl;
				op = new FPConstMultParser(target, wE_in, wF_in, wE_out, wF_out, wF_C, constant);
				addOperator(op);
			}        
		} 	
#endif // HAVE_SOLLYA
		else if(opname=="LeftShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int maxShift = checkStrictyPositive(argv[i++], argv[0]);
				map<string, double> inputDelays;
				inputDelays["X"]=0;
				inputDelays["S"]=0;
				cerr << "> LeftShifter, wIn="<<wIn<<", maxShift="<<maxShift<<"\n";
				op = new Shifter(target, wIn, maxShift, Shifter::Left);
				addOperator(op);
			}    
		}
		else if(opname=="RightShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int maxShift = checkStrictyPositive(argv[i++], argv[0]);
				map<string, double> inputDelays;
				inputDelays["X"]=0;
				inputDelays["S"]=0;
				cerr << "> RightShifter, wIn="<<wIn<<", maxShift="<<maxShift<<"\n";
				op = new Shifter(target, wIn, maxShift, Shifter::Right, inputDelays);
				addOperator(op);
			}
		}
		else if(opname=="LZOC"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> LZOC, wIn="<<wIn<<"\n";
				op = new LZOC(target, wIn);
				addOperator(op);
			}
		}
		else if(opname=="LZOCShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn  = checkStrictyPositive(argv[i++], argv[0]);
				int wOut = checkStrictyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					cerr << "> LZOCShifter, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, -1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LZCShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn  = checkStrictyPositive(argv[i++], argv[0]);
				int wOut = checkStrictyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					cerr << "> LZCShifter, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, 0);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LOCShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn  = checkStrictyPositive(argv[i++], argv[0]);
				int wOut = checkStrictyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					cerr << "> LOCShifter, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, 1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LZOCShifterSticky"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn  = checkStrictyPositive(argv[i++], argv[0]);
				int wOut = checkStrictyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					cerr << "> LZOCShifterSticky, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, -1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LZCShifterSticky"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn  = checkStrictyPositive(argv[i++], argv[0]);
				int wOut = checkStrictyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					cerr << "> LZCShifterSticky, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, 0);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LOCShifterSticky"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn  = checkStrictyPositive(argv[i++], argv[0]);
				int wOut = checkStrictyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					cerr << "> LOCShifterSticky, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, 1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="UserDefinedOperator"){
                        // the UserDefinedOperator expects 2 parameters
			int nargs = 2;
			if (i+nargs > argc)
                            /* if there is less than 2 parameters, we output 
                              the help information for Flopoco */
				usage(argv[0]);
			else {
				int param0 = checkStrictyPositive(argv[i++], argv[0]);
				int param1 = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> UserDefinedOperator, param0="<<param0<<", param1=" << param1 << endl  ;
                                op = new UserDefinedOperator(target,param0,param1);
                                addOperator(op);
			}    
		}
		else if(opname=="IntAdder"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> IntAdder, wIn="<<wIn<<endl  ;
				op = new IntAdder(target,wIn);
				addOperator(op);
			}    
		}
		else if(opname=="IntComparator"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int criteria = atoi(argv[i++]);
				
				cerr << "> IntComparator, wIn="<<wIn<<" criteria="<<criteria<<endl  ;
				op = new IntComparator(target,wIn,criteria,false,0);
				addOperator(op);
			}    
		}
		else if(opname=="IntConstComparator"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int criteria = atoi(argv[i++]);
				int constant = atoi(argv[i++]);
				
				cerr << "> IntComparator, wIn="<<wIn<<" criteria="<<criteria<<endl  ;
				op = new IntComparator(target,wIn,criteria, true, constant);
				addOperator(op);
			}    
	}
	
		else if(opname=="MyIntAdder"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int type = atoi(argv[i++]);
				int srl = atoi(argv[i++]);
				int implementation = atoi(argv[i++]);
				int bufferedIn = atoi(argv[i++]);
				cerr << "> IntAdder, wIn="<<wIn<<", frequency="<<target->frequency()<< endl  ;

				map <string, double> delayMap;

				if (!bufferedIn){
					delayMap["X"] = 1e-25;
				}
				
				switch (type) {
					case 0: op = new IntAdder(target, wIn, delayMap, 0, srl, implementation); break; //lut optimized
					case 1: op = new IntAdder(target, wIn, delayMap, 1, srl, implementation); break; //reg
					case 2: op = new IntAdder(target, wIn, delayMap, 2, srl, implementation); break; //slice
					case 3: op = new IntAdder(target, wIn, delayMap, 3, srl, implementation); break; //latency
					default: op = new IntAdder(target,wIn, delayMap, 2, srl, implementation); break;
				}
				addOperator(op);
			}    
		}

		else if(opname=="IntNAdder"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int N   = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> IntNAdder, wIn="<<wIn<<" N="<<N<<endl  ;
				op = new IntNAdder(target,wIn,N);
				addOperator(op);
			}    
		}

		else if(opname=="IntCompressorTree"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int N   = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> IntCompressorTree, wIn="<<wIn<<" N="<<N<<endl  ;
				op = new IntCompressorTree(target,wIn,N);
				addOperator(op);
			}    
		}

		else if(opname=="LongIntAdder"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdder, wIn="<<wIn<<endl  ;
				op = new LongIntAdder(target,wIn);
				addOperator(op);
			}    
		}
		
		else if(opname=="LongIntAdderCmpCmpAdd"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdderCmpCmpAdd, wIn="<<wIn<<endl  ;
				op = new LongIntAdderCmpCmpAdd(target,wIn);
				addOperator(op);
			}    
		}
		
		else if(opname=="LongIntAdderCmpAddAdd"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdderCmpAddAdd, wIn="<<wIn<<endl  ;
				op = new LongIntAdderCmpAddAdd(target,wIn);
				addOperator(op);
			}    
		}

		//HIDDEN
		else if(opname=="IntDualSub"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				int opType = checkBoolean(argv[i++], argv[0]);
				cerr << "> IntDualSub, wIn="<<wIn<<" OpType="<<opType<<endl  ;
				op = new IntDualSub(target,wIn,opType);
				addOperator(op);
			}    
		}
		else if(opname=="IntMultiplier"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wInX    = checkStrictyPositive(argv[i++], argv[0]);
				int wInY    = checkStrictyPositive(argv[i++], argv[0]);
				int sign    =  checkBoolean(argv[i++], argv[0]);
				bool sign_;
				if (sign > 0)
					sign_ = true;
				else
					sign_ = false;
				cerr << "> IntMultiplier , wInX="<<wInX<<", wInY="<<wInY<<" signed="<<sign_<<"\n";
				op = new IntMultiplier(target, wInX, wInY, emptyDelayMap, sign_);
				addOperator(op);
			}
		}

//		else if(opname=="IntTruncMultiplier"){
//			int nargs = 4;
//			if (i+nargs > argc)
//				usage(argv[0]);
//			else {
//				int wInX = checkStrictyPositive(argv[i++], argv[0]);
//				int wInY = checkStrictyPositive(argv[i++], argv[0]);
//				float ratio = atof(argv[i++]);
//				int k    = atoi(argv[i++]);
//				
//				cerr << "> IntTruncMultiplier , wInX="<<wInX<<", wInY="<<wInY<<"\n";
//				
//				op = new IntTruncMultiplier(target, wInX, wInY, ratio, k);
//				addOperator(op);
//			}
//		}
		else if(opname=="SignedIntMultiplier"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wInX = checkStrictyPositive(argv[i++], argv[0]);
				int wInY = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> SignedIntMultiplier , wInX="<<wInX<<", wInY="<<wInY<<"\n";
				op = new SignedIntMultiplier(target, wInX, wInY);
				addOperator(op);
			}
		}
		else if(opname=="IntKaratsuba"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wIn = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> IntKaratsuba , wIn="<<wIn<<"\n";
				op = new IntKaratsuba(target, wIn);
				addOperator(op);
			}    
		}   
		else if(opname=="IntTilingMultiplier"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wInX = checkStrictyPositive(argv[i++], argv[0]);
				int wInY = checkStrictyPositive(argv[i++], argv[0]);
				float r = atof(argv[i++]);
				if (r<0)
					throw "Ratio > 0";
				int maxTimeInMinutes = atoi(argv[i++]);
				
				cerr << "> IntTilingMultiplier , wInX="<<wInX<<", wInY="<<wInY<<" ratio=" << r << " maxTimeInMinutes= "<< maxTimeInMinutes << "\n";
				op = new IntTilingMult(target, wInX, wInY, r, maxTimeInMinutes);
				addOperator(op);
			}
		}
		else if(opname=="IntTruncMultiplier"){
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wInX = checkStrictyPositive(argv[i++], argv[0]);
				int wInY = checkStrictyPositive(argv[i++], argv[0]);
				float r = atof(argv[i++]);
				if (r<0)
					throw "Ratio > 0";
				int ordError = atoi(argv[i++]);
				if (ordError<0)
					throw "error order is >= 0";
				int useLimits = atoi(argv[i++]);
				if (!((useLimits==0)||(useLimits==1)))
					throw "useLimits is 0 or 1";				
				int maxTimeInMinutes = atoi(argv[i++]);
				int sign = atoi(argv[i++]);

				cerr << "> IntTruncMult , wInX="<<wInX<<", wInY="<<wInY<<" ratio=" << r <<" error order= "<<ordError<< " useLimits="<<useLimits<<" maxTimeInMinutes="<<maxTimeInMinutes<<" signed="<<sign<<"\n";
				op = new IntTruncMultiplier(target, wInX, wInY, r,ordError,useLimits, maxTimeInMinutes, true, sign);
				addOperator(op);
			}
		}		
		else if(opname=="FPAdder"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wF = checkStrictyPositive(argv[i++], argv[0]);
				
				cerr << "> FPAdder , wE="<<wE<<", wF="<<wF<<" \n";
				op = new FPAdder(target, wE, wF, wE, wF, wE, wF);
				addOperator(op);
			}
		}
		else if(opname=="FPAdderSinglePath"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wF = checkStrictyPositive(argv[i++], argv[0]);
				
				cerr << "> FPAdder , wE="<<wE<<", wF="<<wF<<" \n";
				op = new FPAdderSinglePath(target, wE, wF, wE, wF, wE, wF);
				addOperator(op);
			}
		}		
		else if(opname=="FPAdder3Input"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wF = checkStrictyPositive(argv[i++], argv[0]);
				
				cerr << "> FPAdder3Input , wE="<<wE<<", wF="<<wF<<" \n";
				op = new FPAdder3Input(target, wE, wF);
				addOperator(op);
			}
	}	
		else if(opname=="FPFMAcc"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wF = checkStrictyPositive(argv[i++], argv[0]);
				int adderLatency = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> FPFMAcc , wE="<<wE<<", wF="<<wF<<" adderLatency="<<adderLatency<< " \n";
				op = new FPFMAcc(target, wE, wF, adderLatency);
				addOperator(op);
			}
		}		
		else if(opname=="FPJacobi"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wF = checkStrictyPositive(argv[i++], argv[0]);
				int l0 = checkStrictyPositive(argv[i++], argv[0]);
				int l1 = checkStrictyPositive(argv[i++], argv[0]);
				int l2 = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> FPJacobi , wE="<<wE<<", wF="<<wF<<" l0="<<l0<<" l1="<<l1<<" l2="<<l2<< " \n";
				op = new FPJacobi(target, wE, wF, l0, l1, l2);
				addOperator(op);
			}
	}	
		else if(opname=="Fix2FP"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int LSB = atoi(argv[i++]);//checkStrictyPositive(argv[i++], argv[0]);
				int MSB = atoi(argv[i++]);//checkStrictyPositive(argv[i++], argv[0]);
				int sign = atoi(argv[i++]);
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wF = checkStrictyPositive(argv[i++], argv[0]);
				
				cerr << "> Fix2FP, LSB="<<LSB<<", MSB="<<MSB<<", wE="<<wE<<", wF="<<wF<<" \n";
				op = new Fix2FP(target, LSB, MSB, sign,wE, wF);
				addOperator(op);
			}
		}
		// else if(opname=="CoilInductance"){
		// 	int nargs = 7;
		// 	if (i+nargs > argc)
		// 		usage(argv[0]);
		// 	else {
		// 		int LSBI = atoi(argv[i++]);
		// 		int MSBI = atoi(argv[i++]);
		// 		int wE = checkStrictyPositive(argv[i++], argv[0]);
		// 		int wF = checkStrictyPositive(argv[i++], argv[0]);				
		// 		int MaxMSBO= atoi(argv[i++]);
		// 		int LSBO = atoi(argv[i++]);
		// 		int MSBO = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoilInductance  LSBI="<<LSBI<<", MSBI="<<MSBI<<",wEIn="<<wE<<",wFIn"<<wF<<", MaxMSBO="<<MaxMSBO<<", LSBO="<<LSBO<<", MSBO="<<MSBO<<" \n";
		// 		op = new CoilInductance(target, LSBI, MSBI,wE,wF,MaxMSBO,LSBO,MSBO,pa);
		// 		addOperator(op);
		// 	}
		// }
		// else if(opname=="CoordinatesTableX"){
		// 	int nargs = 4;
		// 	if (i+nargs > argc)
		// 		usage(argv[0]);
		// 	else {
		// 		int wIn = checkStrictyPositive(argv[i++], argv[0]);
		// 		int LSB = atoi(argv[i++]);
		// 		int MSB = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoordinatesTableX, wIn="<<wIn<<", LSB="<<LSB<<", MSB="<<MSB<<" \n";
		// 		op = new CoordinatesTableX(target, wIn,LSB, MSB,pa);
		// 		addOperator(op);
		// 	}
		// }
		// else if(opname=="CoordinatesTableZ"){
		// 	int nargs = 4;
		// 	if (i+nargs > argc)
		// 		usage(argv[0]);
		// 	else {
		// 		int wIn = checkStrictyPositive(argv[i++], argv[0]);
		// 		int LSB = atoi(argv[i++]);
		// 		int MSB = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoordinatesTableZ, wIn="<<wIn<<", LSB="<<LSB<<", MSB="<<MSB<<" \n";
		// 		op = new CoordinatesTableZ(target, wIn,LSB, MSB,pa);
		// 		addOperator(op);
		// 	}
		// }
		// else if(opname=="CoordinatesTableY"){
		// 	int nargs = 4;
		// 	if (i+nargs > argc)
		// 		usage(argv[0]);
		// 	else {
		// 		int wIn = checkStrictyPositive(argv[i++], argv[0]);
		// 		int LSB = atoi(argv[i++]);
		// 		int MSB = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoordinatesTableY, wIn="<<wIn<<", LSB="<<LSB<<", MSB="<<MSB<<" \n";
		// 		op = new CoordinatesTableY(target, wIn,LSB, MSB,pa);
		// 		addOperator(op);
		// 	}
		// }
		else if(opname=="FPMultiplier"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wFIn = checkStrictyPositive(argv[i++], argv[0]);
				int wFOut = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> FPMultiplier , wE="<<wE<<", wFIn="<<wFIn<<", wFOut="<<wFOut<<" \n";
				op = new FPMultiplier(target, wE, wFIn, wE, wFIn, wE, wFOut, 1);
				addOperator(op);
			}
		}
		else if(opname=="FPMultiplierKaratsuba"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wFIn = checkStrictyPositive(argv[i++], argv[0]);
				int wFOut = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> FPMultiplierKaratsuba , wE="<<wE<<", wFIn="<<wFIn<<", wFOut="<<wFOut<<" \n";
				op = new FPMultiplierKaratsuba(target, wE, wFIn, wE, wFIn, wE, wFOut, 1);
				addOperator(op);
			}
		}  
		else if(opname=="FPMultiplierTiling"){
			int nargs = 5; 
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wFIn = checkStrictyPositive(argv[i++], argv[0]);
				int wFOut = checkStrictyPositive(argv[i++], argv[0]);
				float r = atof(argv[i++]);
				int maxTimeInMinutes = atoi(argv[i++]);

				cerr << "> FPMultiplierTiling , wE="<<wE<<", wFIn="<<wFIn<<", wFOut="<<wFOut<<" ratio="<<r<<" maxTime="<<maxTimeInMinutes<<" \n";
				op = new FPMultiplierTiling(target, wE, wFIn, wE, wFIn, wE, wFOut, 1, r, maxTimeInMinutes);
				addOperator(op);
			}
		}  
		else if(opname=="FPSquarer"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wFX = checkStrictyPositive(argv[i++], argv[0]);
				int wFR = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> FPSquarer , wE="<<wE<<", wFX="<<wFX<<" wFR="<<wFR<< " \n";
				op = new FPSquarer(target, wE, wFX, wFR);
				addOperator(op);
			}
		} 
		else if (opname == "FPDiv")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> FPDiv: wE=" << wE << " wF=" << wF << endl;
			op = new FPDiv(target, wE, wF);
			addOperator(op);
		}
		else if (opname == "FPSqrt")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> FPSqrt: wE=" << wE << " wF=" << wF  << endl;
			op = new FPSqrt(target, wE, wF);
			addOperator(op);
		}
#ifdef HAVE_SOLLYA
		else if (opname == "FPSqrtPoly")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
//			int correctlyRounded = checkBoolean(argv[i++], argv[0]);
			int degree = checkStrictyPositive(argv[i++], argv[0]);
//			cerr << "> FPSqrtPoly: wE=" << wE << " wF=" << wF << " correctlyRounded="<< correctlyRounded << " degree =" << degree << endl;
			cerr << "> FPSqrtPoly: wE=" << wE << " wF=" << wF << " degree =" << degree << endl;
			op = new FPSqrtPoly(target, wE, wF, 0, degree);
			addOperator(op);
		}
#endif // HAVE_SOLLYA

		else if(opname=="LongAcc"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wEX = checkStrictyPositive(argv[i++], argv[0]);
				int wFX = checkStrictyPositive(argv[i++], argv[0]);
				int MaxMSBX = atoi(argv[i++]); // may be negative
				int LSBA = atoi(argv[i++]); // may be negative
				int MSBA = atoi(argv[i++]); // may be negative
				cerr << "> Long accumulator , wEX="<<wEX<<", wFX="<<wFX<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<"\n";
				op = new LongAcc(target, wEX, wFX, MaxMSBX, LSBA, MSBA);
				addOperator(op);
			}
		}

		// hidden and undocumented
		else if(opname=="LongAccPrecTest"){
			int nargs = 6; // same as LongAcc, plus an iteration count
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wEX = atoi(argv[i++]);
				int wFX = atoi(argv[i++]);
				int MaxMSBX = atoi(argv[i++]);
				int LSBA = atoi(argv[i++]);
				int MSBA = atoi(argv[i++]);
				int n = atoi(argv[i++]);
				cerr << "> Test of long accumulator accuracy, wEX="<<wEX<<", wFX="<<wFX<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<", "<< n << " tests\n";
				LongAcc * op = new LongAcc(target, wEX, wFX, MaxMSBX, LSBA, MSBA);
				// op->test_precision(n);
				op->test_precision2();
			}    
		}

		else if(opname=="LongAcc2FP"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int LSBA = atoi(argv[i++]); // may be negative
				int MSBA = atoi(argv[i++]); // may be negative
				int wE_out = checkStrictyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictyPositive(argv[i++], argv[0]);
				cerr << "> Post-Normalization unit for Long accumulator, LSBA="<<LSBA<<", MSBA="<<MSBA<<" wE_out="<<wE_out<<", wF_out="<<wF_out<<"\n";
				op = new LongAcc2FP(target, LSBA, MSBA, wE_out, wF_out);
				addOperator(op);
			}
		}
#ifndef _WIN32
		// hidden and undocumented
		else if(opname=="DotProdPrecTest"){
			int nargs = 7; // same as LongAcc, plus an iteration count
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wFX = checkStrictyPositive(argv[i++], argv[0]);
				int wFY = checkStrictyPositive(argv[i++], argv[0]);
				int MaxMSBX = atoi(argv[i++]);
				int LSBA = atoi(argv[i++]);
				int MSBA = atoi(argv[i++]);
				int n = atoi(argv[i++]);
				cerr << "> Test of DotProduct accuracy, wEX="<<wE<<", wFX="<<wFX<<", wFY="<<wFY<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<", "<< n << " tests\n";
				DotProduct * op = new DotProduct(target, wE, wFX, wFY, MaxMSBX, LSBA, MSBA);
				op->test_precision(n);
			}    
		}
#endif
		else if(opname=="DotProduct"){
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				int wE = checkStrictyPositive(argv[i++], argv[0]);
				int wFX = checkStrictyPositive(argv[i++], argv[0]);
				int wFY = checkStrictyPositive(argv[i++], argv[0]);
				int MaxMSBX = atoi(argv[i++]); // may be negative
				int LSBA = atoi(argv[i++]); // may be negative
				int MSBA = atoi(argv[i++]); // may be negative
				cerr << "> DotProduct , wE="<<wE<<", wFX="<<wFX<<", wFY="<<wFY<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<"\n";
				op = new DotProduct(target, wE, wFX, wFY, MaxMSBX, LSBA, MSBA);
				addOperator(op);
			}
		}
//		else if(opname=="PolynomialEvaluator"){
//			int nargs = 1;
//			if (i+nargs > argc)
//				usage(argv[0]);
//			else {
//				int prec = atoi(argv[i++]); // may be negative
//				FixedPointCoefficient* f0 = new FixedPointCoefficient( 27, 0);
//				FixedPointCoefficient* f1 = new FixedPointCoefficient( 17,-1);
//				FixedPointCoefficient* f2 = new FixedPointCoefficient( 9, -3);

//				YVar* y = new YVar(16, -6);
//				
//				vector<FixedPointCoefficient*> coef;
//				coef.push_back(f0);
//				coef.push_back(f1);
//				coef.push_back(f2);

//				op = new PolynomialEvaluator(target, coef, y, prec);
//				addOperator(op);
//			}
//		}
				
#ifndef _WIN32
#ifdef HAVE_HOTBM
		else if (opname == "HOTBM") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			string func = argv[i++];
			int wI = checkStrictyPositive(argv[i++], argv[0]);
			int wO = checkStrictyPositive(argv[i++], argv[0]);
			int n  = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> HOTBM func='" << func << "', wI=" << wI << ", wO=" << wO <<endl;
			op = new HOTBM(target, func, "", wI, wO, n);
			addOperator(op);
		}

		else if (opname == "HOTBMFX") {
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			string func = argv[i++];
			int wE_in = atoi(argv[i++]); // may be negative
			int wF_in = atoi(argv[i++]); // may be negative
			int wE_out = atoi(argv[i++]); // may be negative
			int wF_out = atoi(argv[i++]); // may be negative
			int n  = checkStrictyPositive(argv[i++], argv[0]);
			
			// Input of HOTBM is unsigned, but output is in 2's-complement!
			// Use 1 more bit for the output sign
			int wI = wE_in + wF_in;
			int wO = wE_out + wF_out + 1;
			double xmin = 0;
			double xmax = ldexp(1.0, wE_in);
			double scale = ldexp(1.0, -wE_out);
			cerr << "> HOTBM func='" << func << "', wE_in=" << wE_in << ", wF_in=" << wF_in
			     << ", wE_out=" << wE_out << ", wF_out=" << wF_out << endl;
			op = new HOTBM(target, func, "", wI, wO, n, xmin, xmax, scale);
			addOperator(op);
		}
		
		else if (opname == "HOTBMRange") {
			int nargs = 7;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			string func = argv[i++];
			int wI = checkStrictyPositive(argv[i++], argv[0]);
			int wO = checkStrictyPositive(argv[i++], argv[0]);
			int n  = checkStrictyPositive(argv[i++], argv[0]);
			double xmin = atof(argv[i++]);
			double xmax = atof(argv[i++]);

			// xmax < xmin is a valid use case...
			double scale = atof(argv[i++]);
			cerr << "> HOTBM func='" << func << "', wI=" << wI << ", wO=" << wO
			     << ", xmin=" << xmin << ", xmax=" << xmax << ", scale=" << scale <<endl;
			op = new HOTBM(target, func, "", wI, wO, n, xmin, xmax, scale);
			addOperator(op);
		}
		
#endif // HAVE_HOTBM

#endif

		else if (opname == "FPExp")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> FPExp: wE=" << wE << " wF=" << wF << endl;
			op = new FPExp(target, wE, wF, 0, 0);
			addOperator(op);
		}


		else if (opname == "FPExpExpert")
		{
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int k=atoi(argv[i++]);
			int d=atoi(argv[i++]);
			int g=atoi(argv[i++]);
			int fullInput=checkBoolean(argv[i++],  argv[0]);
			cerr << "> FPExp: wE=" << wE << " wF=" << wF << endl;
			op = new FPExp(target, wE, wF, k, d, g, fullInput);
			addOperator(op);
		}

		else if (opname == "FPLog")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int inTableSize=atoi(argv[i++]);
			cerr << "> FPLog: wE=" << wE << " wF=" << wF << endl;
			op = new FPLog(target, wE, wF, inTableSize);
			addOperator(op);
		}

		else if (opname == "FPPowrExpert")
		{
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0]); // and exit 

			//int logTableSize, int expTableSize, int expDegree, int expG, int logG
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int logTableSize=atoi(argv[i++]);
			int expTableSize=atoi(argv[i++]);
			int expDegree=atoi(argv[i++]);
			cerr << "> FPPowr: wE=" << wE << " wF=" << wF << endl;
			op = new FPPow(target, wE, wF,logTableSize, expTableSize, expDegree);
			addOperator(op);
		}

		else if (opname == "FPPowr")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit 

			//int logTableSize, int expTableSize, int expDegree, int expG, int logG
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> FPPowr: wE=" << wE << " wF=" << wF << endl;
			op = new FPPow(target, wE, wF,0, 0, 0);
			addOperator(op);
		}

		else if (opname == "InputIEEE")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wEI = checkStrictyPositive(argv[i++], argv[0]);
			int wFI = checkStrictyPositive(argv[i++], argv[0]);
			int wEO = checkStrictyPositive(argv[i++], argv[0]);
			int wFO = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> InputIEEE: wEI=" << wEI << " wFI=" << wFI << " wEO=" << wEO << " wFO=" << wFO << endl;
			op = new InputIEEE(target, wEI, wFI, wEO, wFO);
			addOperator(op);
		}

		else if (opname == "OutputIEEE")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wEI = checkStrictyPositive(argv[i++], argv[0]);
			int wFI = checkStrictyPositive(argv[i++], argv[0]);
			int wEO = checkStrictyPositive(argv[i++], argv[0]);
			int wFO = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> OutputIEEE: wEI=" << wEI << " wFI=" << wFI << " wEO=" << wEO << " wFO=" << wFO << endl;
			op = new OutputIEEE(target, wEI, wFI, wEO, wFO);
			addOperator(op);
		}

		else if (opname == "Collision")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int optimize = checkBoolean(argv[i++], argv[0]);
			cerr << "> Collision: wE=" << wE << " wF=" << wF << (optimize==0? " using FP operators" : " optimized version") << endl;
			op = new Collision(target, wE, wF, optimize);
			addOperator(op);
		}
		else if (opname == "IntSquarer")
		{
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wIn = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> IntSquarer: wIn=" << wIn << endl;
			op = new IntSquarer(target, wIn);
			addOperator(op);
		}
#ifndef _WIN32
#ifdef HAVE_LNS
		else if (opname == "LNSAddSub")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> LNSAddSub: wE=" << wE << " wF=" << wF << endl;
			op = new LNSAddSub(target, wE, wF);
			if(cl_name!="")	op->setName(cl_name);
			addOperator(op);
		}
		else if (opname == "LNSMul")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> LNSMul: wE=" << wE << " wF=" << wF << endl;
			op = new LNSMul(target, wE, wF);
			if(cl_name!="")	op->setName(cl_name);
			addOperator(op);
		}
		else if (opname == "LNSDiv")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> LNSDiv: wE=" << wE << " wF=" << wF << endl;
			op = new LNSDiv(target, wE, wF);
			addOperator(op);
		}
		else if (opname == "LNSSqrt")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> LNSSqrt: wE=" << wE << " wF=" << wF << endl;
			op = new LNSSqrt(target, wE, wF);
			if(cl_name!="")	op->setName(cl_name);
			addOperator(op);
		}
		else if (opname == "LNSAdd")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int o = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> LNSAdd: wE=" << wE << " wF=" << wF << " o=" << o << endl;
			op = new LNSAdd(target, wE, wF, o);
			addOperator(op);
		}
		// Undocumented LNS operators, for debugging purposes
#if 0
		else if (opname == "CotranF1")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int j = checkStrictyPositive(argv[i++], argv[0]);
			int wE = atoi(argv[i++]);
			cerr << "> CotranF1: wF=" << wF << " j=" << j << " wE=" << wE << endl;
			op = new TableOp(target, new CotranF1Table(wF, j, wE));
			addOperator(op);
		}
		else if (opname == "CotranF2")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int j = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> CotranF2: wF=" << wF << " j=" << j << endl;
			op = new TableOp(target, new CotranF2Table(wF, j));
			addOperator(op);
		}
		else if (opname == "CotranF3")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int j = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> CotranF3: wF=" << wF << " j=" << j << endl;
			op = new TableOp(target, new CotranF3Table(wF, j));
			addOperator(op);
		}
#endif // #if 0
		else if (opname == "Cotran")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int j = checkStrictyPositive(argv[i++], argv[0]);
			int wECotran = atoi(argv[i++]);
			int o = atoi(argv[i++]);
			cerr << "> Cotran: wE=" << wE << " wF=" << wF << " j=" << j
				<< " wECotran=" << wECotran << " o=" << o << endl;
			op = new Cotran(target, wE, wF, j, wECotran);
			addOperator(op);
		}
		else if (opname == "CotranHybrid")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int j = checkStrictyPositive(argv[i++], argv[0]);
			int wECotran = atoi(argv[i++]);
			int o = atoi(argv[i++]);
			cerr << "> Cotran: wE=" << wE << " wF=" << wF << " j=" << j 
				<< " wECotran=" << wECotran << " o=" << o << endl;
			op = new CotranHybrid(target, wE, wF, j, wECotran);
			addOperator(op);
		}
		else if (opname == "AtanPow")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			int o = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> AtanPow: wE=" << wE << " wF=" << wF << " o=" << o << endl;
			op = new AtanPow(target, wE, wF, o);
			addOperator(op);
		}
#endif

#endif //ifndef _WIN32
		else if (opname == "Wrapper") {
			int nargs = 0;
			if (i+nargs > argc)
				usage(argv[0]);
			else {
				if(oplist.empty()){
					cerr<<"ERROR: Wrapper has no operator to wrap (it should come after the operator it wraps)"<<endl;
					usage(argv[0]);
				}
				Operator* toWrap = oplist.back();
				cerr << "> Wrapper for " << toWrap->getName()<<endl;
				op =new Wrapper(target, toWrap);
				addOperator(op);
			}
		}
#ifdef HAVE_HOTBM
		else if (opname == "PolyTableGenerator") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			string func = argv[i++];
			int wI = checkStrictyPositive(argv[i++], argv[0]);
			int wO = atoi(argv[i++]);
			int n  = checkStrictyPositive(argv[i++], argv[0]);
			
						
			cerr << "> PolyTableGenerator func='" << func << "', wI=" << wI << ", wO=" << wO <<endl;	
			
			Operator* tg = new PolyTableGenerator(target, func, wI, wO, n);
				addOperator(tg);
			
		}
		

		else if (opname == "FunctionEvaluator") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			string func = argv[i++];
			int wI = checkStrictyPositive(argv[i++], argv[0]);
			int wO = atoi(argv[i++]);
			int n  = checkStrictyPositive(argv[i++], argv[0]);
			
				
			cerr << "> FunctionEvaluator func='" << func << "', wI=" << wI << ", wO=" << wO << ", degree=" << n << endl;	
			string arg=func+",0,1,1"; // we are not sure it works for other values
			Operator* tg = new FunctionEvaluator(target, arg, wI, wO, n);
			addOperator(tg);
		}

		else if (opname == "FPPipeline") {
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			string expr = argv[i++];
			int wE = checkStrictyPositive(argv[i++], argv[0]);
			int wF = checkStrictyPositive(argv[i++], argv[0]);
			cerr << "> FPPipeline expr='" << expr << "', wE=" << wE << ", wF=" << wF << endl;	
			Operator* tg = new FPPipeline(target, expr, wE, wF);
			addOperator(tg);
		}

#endif

		else if (opname == "TestBench") {
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			if(oplist.empty()){
				cerr<<"ERROR: TestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0]); // and exit
			}
			int n = checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = oplist.back();
			Operator* op = new TestBench(target, toWrap, n);
			cerr << "> TestBench for " << toWrap->getName()<<endl;
			addOperator(op);
			cerr << "To run the simulation using ModelSim, type the following in 'vsim -c':" <<endl;
			cerr << tab << "vdel -all -lib work" <<endl;
			cerr << tab << "vlib work" <<endl;
			cerr << tab << "vcom " << filename <<endl;
			cerr << tab << "vsim " << op->getName() <<endl;
			cerr << tab << "add wave -r *" <<endl;
			cerr << tab << "run " << ((TestBench*)op)->getSimulationTime() <<"ns" << endl;
			cerr << "To run the simulation using gHDL, type the following in a shell prompt:" <<endl;
			cerr << tab << "ghdl -a --ieee=synopsys -fexplicit "<< filename <<endl;
			cerr << tab << "ghdl -e --ieee=synopsys -fexplicit " << op->getName() <<endl;
			cerr << tab << "ghdl -r --ieee=synopsys " << op->getName() << " --vcd=" << op->getName() << ".vcd" <<endl;
			cerr << tab << "gtkwave " << op->getName() << ".vcd" << endl;
		}
		
		else if (opname == "TestBenchFile") {
			/* Using a file to store IO */
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			if(oplist.empty()){
				cerr<<"ERROR: TestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0]); // and exit
			}
			int n = atoi(argv[i++]);//checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = oplist.back();
			Operator* op = new TestBench(target, toWrap, n, true);
			cerr << "> TestBench for " << toWrap->getName()<<endl;
			addOperator(op);
			cerr << "To run the simulation using ModelSim, type the following in 'vsim -c':" <<endl;
			cerr << tab << "vdel -all -lib work" <<endl;
			cerr << tab << "vlib work" <<endl;
			cerr << tab << "vcom " << filename <<endl;
			cerr << tab << "vsim " << op->getName() <<endl;
			cerr << tab << "add wave -r *" <<endl;
			cerr << tab << "run " << ((TestBench*)op)->getSimulationTime() << "ns" << endl;
			cerr << "To run the simulation using gHDL, type the following in a shell prompt:" <<endl;
			cerr << tab << "ghdl -a --ieee=synopsys -fexplicit "<< filename <<endl;
			cerr << tab << "ghdl -e --ieee=synopsys -fexplicit " << op->getName() <<endl;
			cerr << tab << "ghdl -r --ieee=synopsys " << op->getName() << " --vcd=" << op->getName() << ".vcd" <<endl;
			cerr << tab << "gtkwave " << op->getName() << ".vcd" << endl;
		}
#if 0  //We should resurrect it some day
		else if (opname == "BigTestBench") {
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0]); // and exit
			if(oplist.empty()){
				cerr<<"ERROR: BigTestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0]); // and exit
			}
			int n = checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = oplist.back();
			cerr << "> BigTestBench for " << toWrap->getName()<<endl;
			Operator* op = new BigTestBench(target, toWrap, n);
			addOperator(op);
		}
 #endif
		else  {
			cerr << "ERROR: Problem parsing input line, exiting";
			usage(argv[0]);
		}
	} while (i<argc);
	return true;
}

int main(int argc, char* argv[] )
{
	uint32_t i;
	string srcFileName = "main.cpp";

	

	target = new Virtex4();

	try {
		parseCommandLine(argc, argv);
	} catch (char const * s) {
		cerr << "Exception while parsing command line: " << s << endl;
		return 1;
	} catch (std::string s){
		cerr << "Exception while parsing command line: " << s << endl;
		return 1;
	}

	if(oplist.empty()) {
		cerr << "Nothing to generate, exiting\n";
		exit(EXIT_FAILURE);
	}

	ofstream file;
	file.open(filename.c_str(), ios::out);
	
	for(i=0; i<oplist.size(); i++) {
		try {
			REPORT(FULL, "OPERATOR:"<<oplist[i]->getName());
			REPORT(FULL, "--DECLARE LIST---------------------------------------------------");
			REPORT(FULL, printMapContent(oplist[i]->getDeclareTable()) );
			REPORT(FULL, "--USE LIST-------------------------------------------------------");

//			combinatorialOperator = not(oplist[i]->isSequential());
			oplist[i]->getFlopocoVHDLStream()->flush();
			
			REPORT(FULL, printVectorContent(  (oplist[i]->getFlopocoVHDLStream())->getUseTable()) );

			/* second parse is only for sequential operators */
			if (oplist[i]->isSequential()){
				REPORT (FULL, "--2nd PASS-------------------------------------------------------");
				oplist[i]->parse2();
			}

			oplist[i]->outputVHDL(file);			
		} catch (std::string s) {
			cerr << "Exception while generating '" << oplist[i]->getName() << "': " << s <<endl;
		}
	}
	file.close();
	
	for(int k=oplist.size()-1; k>=0; k--) {
		if (unsigned(k)== (oplist.size()-1))
			oplist[k]->level = 0;
		else{
			string currentName = oplist[k]->getName();
			oplist[k]->level = 0;
			for(int j=oplist.size()-1; j>k; j--) {
				if (oplist[j]->hasComponent(currentName) ){
					oplist[k]->level = oplist[j]->level + 1;
				}	
			}
			
		}
	}	
	
	
	cerr << endl<<"Final report:"<<endl;
	for(i=0; i<oplist.size(); i++) {
		oplist[i]->outputFinalReport();
	}
	cerr<< "Output file: " << filename <<endl;
	return 0;
}
