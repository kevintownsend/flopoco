/*
  the FloPoCo command-line interface
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Authors : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr
            Bogdan Pasca, Bogdan.Pasca@ens-lyon.org

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA-Lyon  
  2008-2014.
  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <mpfr.h>
#include <sollya.h>
#include <cstdlib>



#include "FloPoCo.hpp"


#define BRIGHT 1
#define RED 31
#define OPER 32
#define NEWOPER 32
#define PARAM 34
#define OP(op,paramList)             {cerr << "  "; printf("%c[%d;%dm",27,1,OPER); cerr <<  op; printf("%c[%dm",27,0); cerr<< " "; printf("%c[%d;%dm",27,1,PARAM); cerr << paramList; printf("%c[%dm\n",27,0); } 
#define NEWOP(op,paramList)          {cerr << "  "; printf("%c[%d;%dm",27,1,NEWOPER); cerr <<  op; printf("%c[%dm",27,0); cerr<< " "; printf("%c[%d;%dm",27,1,PARAM); cerr << paramList; printf("%c[%dm\n",27,0); } 


using namespace std;
using namespace flopoco;

// Global variables, useful in this main to avoid parameter passing


	string filename="flopoco.vhdl";
	string cl_name=""; // used for the -name option
	Target* target;
	
	//------------ Resource Estimation --------------------------------
	int reLevel;
	bool resourceEstimationDebug = false;
	//-----------------------------------------------------------------
	
	//------------------ Floorplanning --------------------------------
	bool floorplanning = false;
	bool floorplanningDebug = false;
	ostringstream floorplanMessages;
	//-----------------------------------------------------------------


extern void random_usage(char *name, string opName);
	
extern bool random_parseCommandLine(
	int argc, char *argv[], Target *target,
	std::string opname, int &i
);



void usage(char *name, string opName = ""){
	bool full = (opName=="");

	if ( full ){
		cerr << "\nUsage: "<<name<<" <operator specification list>\n" ;
		cerr << "Each operator specification is one of: \n";
	}
	if ( full )
		cerr << "    ____________ SHIFTERS/LZOC _________________________________________________\n";

	if ( full || opName == "LeftShifter") 
		OP("LeftShifter","wIn MaxShift");
	if ( full || opName == "RightShifter") 
		OP("RightShifter","wIn MaxShift");
	if ( full || opName == "LZOC") 	
		OP("LZOC","wIn");
	if ( full || opName == "LZOCShifter") 	
		OP("LZOCShifter","wIn wOut");
	if ( full || opName == "LZCShifter") 	
		OP("LZCShifter","wIn wOut");
	if ( full || opName == "LOCShifter") 	
		OP("LOCShifter","wIn wOut");
	if ( full || opName == "LZOCShifterSticky") 		
		OP("LZOCShifterSticky","wIn wOut");
	if ( full || opName == "LZCShifterSticky") 		
		OP("LZCShifterSticky","wIn wOut");
	if ( full || opName == "LOCShifterSticky") 			
		OP("LOCShifterSticky","wIn wOut");

	if ( full )
		cerr << "    ____________ ADDERS/SUBTRACTERS ____________________________________________\n";

	if ( full || opName == "IntAdder"){ 	
		OP("IntAdder","wIn");
		cerr << "      Integer adder, possibly pipelined\n";
	}
	if ( full || opName == "IntAdderExpert" || opName == "IntAdder"){ 		
		OP("IntAdderExpert","wIn optimizeType srl implementation bufferedInputs inputDelay");
		cerr << "      Integer adder, multple parameters, possibly pipelined\n";
		cerr << "      optimizeType=<0,1,2,3> 0=LUT 1=REG 2=SLICE 3=LATENCY\n";
		cerr << "      srl=<0,1> Allow SRLs\n";
		cerr << "      implementation=<-1,0,1,2> -1=optimizeType dependent,\n";  
		cerr << "                                 0=Classical, 1=Alternative, 2=Short-Latency\n";
		cerr << "      bufferedInputs=<0,1>\n";
		cerr << "      inputDelay\n";
	}
	if ( full || opName == "IntAdder" || opName == "LongIntAdderAddAddMux")
		OP("LongIntAdderAddAddMux","wIn generation");
	if ( full || opName == "IntAdder" || opName == "LongIntAdderCmpAddInc")
		OP("LongIntAdderCmpAddInc","wIn generation");
	if ( full || opName == "IntAdder" || opName == "LongIntAdderCmpCmpAdd")	
		OP("LongIntAdderCmpCmpAdd","wIn generation");
	if ( full || opName == "IntAdder" || opName == "LongIntAdderAddAddMux" || opName == "LongIntAdderCmpAddInc" || opName == "LongIntAdderCmpCmpAdd")
		cerr << "      generation 1 are the 2010 ver., 2 are the 2011 ver.\n";

	if ( full || opName == "IntAdder" || opName == "IntDualSub"){
		OP("IntDualSub","wIn opType");
		cerr << "      Integer adder/subtracter or dual subtracter, possibly pipelined\n";
		cerr << "      opType: if 1, compute X-Y and X+Y; if 0, compute X-Y and Y-X \n";
	}

	if ( full )
		cerr << "    ____________ MULTIPLICATION AND DIVISION BY CONSTANTS ______________________\n";

	if ( full || opName == "IntMultiplier" || opName == "IntConstMult"){				
		OP( "IntConstMult","w c");
		cerr << "      Integer constant multiplier using shift-and-add: w - input size, c - the constant\n";
	}

	if ( full || opName == "IntMultiplier" || opName == "IntIntKCM"){					
		OP( "IntIntKCM","w c signedInput");
		cerr << "      Integer constant multiplier using KCM: w - input size, c - the constant\n";
	}

	if ( full || opName == "FixRealKCM"){					
		OP( "FixRealKCM","lsbIn msbIn signedInput lsbOut constant useBitheap");
		cerr << "      Faithful multiplier of a fixed-point input by a real constant\n";
		cerr << "      The constant is provided as a Sollya expression, e.g \"log(2)\"\n";
		cerr << "      The multiplier might or might not us bit heaps, based on the value of the useBitheap parameter\n";
	}

	if ( full || opName == "IntConstDiv"){					
		OP( "IntConstDiv","n d alpha");
		cerr << "      Euclidean division of input of size n by d (returning q and r)\n";
		cerr << "      Algorithm uses radix 2^alpha,   alpha=-1 means a sensible default.\n";
	}

	if ( full || opName == "IntConstRem"){					
		OP( "IntConstDiv","n d alpha");
		cerr << "      Remainder of Euclidean division of input of size n by d\n";
		cerr << "      Algorithm uses radix 2^alpha,   alpha=-1 means a sensible default.\n";
	}

	if ( full || opName == "FPConstDiv"){					
		OP( "FPConstDiv","wE wF d");
		cerr << "      Floating-point division by the (small) integer d\n";
		OP( "FPConstDivExpert","wE wF d e alpha");
		cerr << "      Floating-point division by  d.2^e, where d is a small integer\n";
		cerr << "      Algorithm uses radix 2^alpha,   alpha=-1 means a sensible default. \n";  // 
	}
	if ( full || opName == "FPConstMult"){					
		OP( "FPConstMult","wE_in wF_in wE_out wF_out wC constant_expr");
		cerr << "      Faithful floating-point constant multiplier\n";
		cerr << "      last argument is a Sollya expression between double quotes,e.g.\"exp(pi/2)\".\n";
		cerr << "      If wC>1, it is the size in bits on which the constant must be evaluated.\n";
		cerr << "      If wC=0 the size is computed for a faithful result.\n";
	}
	if ( full || opName == "FPConstMult" || opName == "CRFPConstMult"){					
		OP( "CRFPConstMult","wE_in wF_in wE_out wF_out constant_expr");
		cerr << "      Correctly-rounded floating-point constant multiplier\n";
		cerr << "      The constant is provided as a Sollya expression, between double quotes.\n";
	}	
	if ( full || opName == "FPConstMult" || opName == "FPConstMultRational"){					
		OP( "FPConstMultRational","wE_in wF_in wE_out wF_out a b");
		cerr << "      Floating-point constant multiplier by a rational a/b\n";
		cerr << "      Useful for multiplications by simple rational constants such as 2/3 or 1/9\n";
	}
	if ( full || opName == "FPConstMult" || opName == "FPConstMultExpert"){					
		OP("FPConstMultExpert","wE_in wF_in wE_out wF_out cst_sgn cst_exp cst_int_sig");
		cerr << "      Floating-point constant multiplier\n";
		cerr << "      The constant is provided as integral significand and integral exponent.\n";
	}
	if ( full || opName == "FPConstMult" || opName == "FPRealKCM"){					
		NEWOP("FPRealKCM","wE wF constantExpression");
		cerr << "      Floating-point constant multiplier using the KCM algorithm\n";
		cerr << "      last argument is a Sollya expression between double quotes,e.g.\"exp(pi/2)\".\n";
	}


	if ( full )
		cerr << "    ____________ FUNCTION EVALUATORS ___________________________________________\n";

	if ( full || opName == "FixFunctionEval"){					
		OP( "FixFunctionEval", "function x");	
		cerr << "  ** Helper/debug feature, does not generate VHDL **\n";
		cerr << "  Evaluates a function in arbitrary precision\n";
		cerr << "      function - sollya-syntaxed function to implement, between double quotes\n";
		cerr << "      x - arbitrary precision value\n";
	}

	if ( full || opName == "BasicPolyApprox"){					
		OP( "BasicPolyApprox", "function targetAcc addGuardBitsToConstant");	
		cerr << "  ** Helper/debug feature, does not generate VHDL **\n";
		cerr << "  Builds a polynomial approximation of a function, accurate to targetAcc on the interval [0,1) \n";
		cerr << "    function - sollya-syntaxed function to implement, between double quotes\n";
		cerr << "    targetAcc - real number, e.g. 0.00001\n";
		cerr << "    addGuardBitsToConstant: if -1, will do something sensible (relax the size of the constant by the number of guard bits that will be added to ensure faithful evaluation)\n";
		cerr << "                          : if >=0, this number of extra LSB bits will be added to the minimal size of the constant\n";
	}

	if ( full || opName == "PiecewisePolyApprox"){					
		OP( "PiecewisePolyApprox", "function targetAcc degree");	
		cerr << "  ** Helper/debug feature, does not generate VHDL **\n";
		cerr << "  Builds a piecewise polynomial approximation of a function, accurate to targetAcc on the interval [0,1) \n";
		cerr << "    function - sollya-syntaxed function to implement, between double quotes\n";
		cerr << "    targetAcc - real number, e.g. 0.00001\n";
		cerr << "    degree: degree of the polynomial approximation\n";
	}

	if ( full || opName == "FixFunctionByTable" || opName == "FixFunction"){					
		OP( "FixFunctionByTable","function lsbI  msbO lsbO");
		cerr << "  Simple tabulation of a function defined on [0,1)\n";
		cerr << "    lsbI: weight of input LSB, for instance -8 for an 8-bit input (bit weights from -1 to -8)\n";
		cerr << "    msbO and lsbO: weights of output MSB and LSB,\n";
		cerr << "    function: sollya-syntaxed function to implement, e.g. \"sin(x*Pi/2)\" \n";
	}

	if ( full || opName == "FixFunctionBySimplePoly" || opName == "FixFunction"){					
		OP( "FixFunctionBySimplePoly","function lsbI msbO lsbO ");
		cerr << "  A function evaluator using a single polynomial on the interval [0,1), evaluated with Horner scheme \n";
		cerr << "    function: sollya-syntaxed function to implement, e.g. \"sin(x*Pi/2)\" \n";
		cerr << "    lsbI: weight of input LSB, for instance -8 for an 8-bit input (bit weights from -1 to -8)\n";
		cerr << "    msbO and lsbO: weights of output MSB and LSB,\n";
	}

	if ( full || opName=="options"){
		cerr << "________________________________________________________________________________\n";
		cerr << "________________ OPTIONS________________________________________________________\n";
		cerr << "General options that should only appear before any operator specification:\n";
		cerr << "   -outputfile=<output file name>           (default=flopoco.vhdl)\n";
		cerr << "   -target=<Spartan3|Virtex4|Virtex5|Virtex6|StratixII|StratixIII|StratixIV|StratixV|CycloneII|CycloneIII|CycloneIV|CycloneV>      (default=Virtex5)\n";
		cerr << "Options affecting the operators that follow them:\n";
		cerr << "   -pipeline=<yes|no>                       (default=yes)\n";
		cerr << "   -frequency=<target frequency in MHz>     (default=400)\n";
		cerr << "   -clockenable=<yes|no>                    (default is no)\n";
		cerr << "   -DSP_blocks=<yes|no>\n";
		cerr << "       optimize for the use of DSP blocks   (default=yes)\n";
		cerr << "   -name=<entity name>\n";
		cerr << "       defines the name of the VHDL entity of the next operator\n";
		cerr << "   -resourceEstimation=level\n";
		cerr << "       level=0 disables resource estimation (default)\n";
		cerr << "       level=1..3 larger number means more details\n";
		cerr << "   -floorplanning=<yes|no>\n";
		cerr << "       generate a floorplan (experimental, Xilinx only)\n";
		cerr << "Debugging options, affecting the operators that follow them:\n";
		cerr << "   -verbose=<0|1|2|3>     verbosity level. 0: no output (default), 1: basic info, 3: full debug \n";
		cerr << "   -reDebugging=<yes|no>  debug output for resource estimation (default=no)\n";
		cerr << "   -flpDebugging=<yes|no> debug output for floorplanning (default=no)\n";
	}
	exit (EXIT_FAILURE);
}

int checkStrictlyPositive(char* s, char* cmd) {
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
		cerr << "Updating entity name to: " << cl_name << endl;
		op->changeName(cl_name);
		cl_name="";
	}
	op->addToGlobalOpList();
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
					if(!target->getGlobalOpListRef()->empty())
						cerr << "WARNING: global option "<<o<<" should come before any operator specification" << endl; 
					filename=v;
				}
				else if (o == "verbose") {
					verbose = atoi(v.c_str()); // there must be a more direct method of string
					if (verbose<0 || verbose>4) {
						cerr<<"ERROR: verbose should be 1, 2 or 3,    got "<<v<<"."<<endl;
						usage(argv[0], "options");
					}
				}
				else if (o == "target") {
					Target* oldTarget=target;
					if(!target->getGlobalOpListRef()->empty()){
								cerr<<"ERROR: target should be changed before any component is defined"<<endl; 
								usage(argv[0],"options");
					}
					if(v=="Virtex4") target=new Virtex4();
					else if (v=="Virtex5") target=new Virtex5();
					else if (v=="Virtex6") target=new Virtex6();
					else if (v=="Spartan3") target=new Spartan3();
					else if (v=="StratixII") target=new StratixII();
					else if (v=="StratixIII") target=new StratixIII();
					else if (v=="StratixIV") target=new StratixIV();
					else if (v=="StratixV") target=new StratixV();
					else if (v=="CycloneII") target=new CycloneII();
					else if (v=="CycloneIII") target=new CycloneIII();
					else if (v=="CycloneIV") target=new CycloneIV();
					else if (v=="CycloneV") target=new CycloneV();
					else {
						cerr<<"ERROR: unknown target: "<<v<<endl;
						usage(argv[0],"options");
					}
					// if previous options had changed it
					target->setFrequency(oldTarget->frequency());
					target->setUseHardMultipliers(oldTarget->hasHardMultipliers());
					if (oldTarget->isPipelined()) 
						target->setPipelined();
					else 
						target->setNotPipelined();
					target->setClockEnable(oldTarget->useClockEnable());
					delete(oldTarget);
				}

				else if (o == "pipeline") {
					if(v=="yes") target->setPipelined();
					else if(v=="no")  target->setNotPipelined();
					else {
						cerr<<"ERROR: pipeline option should be yes or no,    got "<<v<<"."<<endl; 
						usage(argv[0],"options");
					}
				}
				else if (o == "clockenable") {
					if(v=="yes") target->setClockEnable(true);
					else if(v=="no")  target->setClockEnable(false);
					else {
						cerr<<"ERROR: clockenable option should be yes or no,    got "<<v<<"."<<endl; 
						usage(argv[0],"options");
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
						usage(argv[0],"options");
					}
				}
				else if (o == "name") {
					cl_name=v; // TODO?  check it is a valid VHDL entity name 
				}
			}
		}

		//-------------------  SHIFTERS AND LZOCS ----------------------

		else if(opname=="LeftShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int maxShift = checkStrictlyPositive(argv[i++], argv[0]);
				map<string, double> inputDelays;
				inputDelays["X"]=0;
				inputDelays["S"]=0;
				op = new Shifter(target, wIn, maxShift, Shifter::Left);
				addOperator(op);
			}    
		}
		else if(opname=="RightShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int maxShift = checkStrictlyPositive(argv[i++], argv[0]);
				map<string, double> inputDelays;
				inputDelays["X"]=0;
				inputDelays["S"]=0;
				op = new Shifter(target, wIn, maxShift, Shifter::Right, inputDelays);
				addOperator(op);
			}
		}
		else if(opname=="LZOC"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				op = new LZOC(target, wIn);
				addOperator(op);
			}
		}
		else if(opname=="LZOCShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn  = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, -1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LZCShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn  = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, 0);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LOCShifter"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn  = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, 1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LZOCShifterSticky"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn  = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, -1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LZCShifterSticky"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn  = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, 0);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}
		else if(opname=="LOCShifterSticky"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn  = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				if (wIn > 1){
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, 1);
					addOperator(op);
				}else
					cerr << "wIn must be > 1"<<endl;
			}
		}

		else if(opname=="IntAdder"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0], opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				op = new IntAdder(target,wIn, inDelayMap("X",target->ffDelay() + target->localWireDelay()) );
				addOperator(op);
			}    
		}
		else if(opname=="IntAdderExpert"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int type = atoi(argv[i++]);
				int srl = atoi(argv[i++]);
				int implementation = atoi(argv[i++]);
				int bufferedIn = atoi(argv[i++]);
				double inputDelay = atof(argv[i++]);
				map <string, double> delayMap;

				if (bufferedIn){
					delayMap["X"] = target->ffDelay() + target->localWireDelay() + 1.0e-25;
				}else{
					delayMap["X"] = inputDelay;
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
		else if(opname=="IntComparator"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int criteria = atoi(argv[i++]);
				
				op = new IntComparator(target,wIn,criteria,false,0);
				addOperator(op);
			}    
		}
		else if(opname=="IntConstComparator"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int criteria = atoi(argv[i++]);
				int constant = atoi(argv[i++]);
				
				op = new IntComparator(target,wIn,criteria, true, constant);
				addOperator(op);
			}    
		}

		//-------------------CONSTANT MULTIPLICATION AND DIVISION ----------------------
		
		
		else if(opname=="IntIntKCM"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int w = atoi(argv[i++]);
				mpz_class mpc(argv[i++]);
				int signedInput = checkBoolean(argv[i++], argv[0]);

				op = new IntIntKCM(target, w, mpc, signedInput);
				addOperator(op);
			}        
		}
		
		else if(opname=="IntConstMult"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int w = atoi(argv[i++]);
				mpz_class mpc(argv[i++]);
				op = new IntConstMult(target, w, mpc);
				addOperator(op);
			}        
		}

		else if(opname=="IntConstDiv"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int n = checkStrictlyPositive(argv[i++], argv[0]);
				int d = checkStrictlyPositive(argv[i++], argv[0]);
				int alpha = atoi(argv[i++]);
				op = new IntConstDiv(target, n, d, alpha);
				addOperator(op);
			} 
		}

		else if(opname=="IntConstRem"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int n = checkStrictlyPositive(argv[i++], argv[0]);
				int d = checkStrictlyPositive(argv[i++], argv[0]);
				int alpha = atoi(argv[i++]);
				op = new IntConstDiv(target, n, d, alpha, true);
				addOperator(op);
			} 
		}


		else if(opname=="FPConstMultRational"){
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictlyPositive(argv[i++], argv[0]);
				int a  = atoi(argv[i++]); 
				int b  = checkStrictlyPositive(argv[i++], argv[0]); 
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, a, b);
				addOperator(op);
			}        
		} 	

		else if(opname=="FPConstMultExpert"){
			int nargs = 7;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictlyPositive(argv[i++], argv[0]);
				int cst_sgn  = checkSign(argv[i++], argv[0]); 
				int cst_exp  = atoi(argv[i++]); // TODO no check on this arg
				mpz_class cst_sig(argv[i++]);
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, cst_sgn, cst_exp, cst_sig);
				addOperator(op);
			}        
		} 	
		else if(opname=="FPConstDiv"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				int d = checkStrictlyPositive(argv[i++], argv[0]);
				op = new FPConstDiv(target, wE, wF, wE, wF, d, 0, -1); // exponent 0, alpha = default
				addOperator(op);
			}        
		}
		else if(opname=="FPConstDivExpert"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				int d = checkStrictlyPositive(argv[i++], argv[0]);
				int e = atoi(argv[i++]);
				int alpha = atoi(argv[i++]);
				op = new FPConstDiv(target, wE, wF, wE, wF, d, e, alpha);
				addOperator(op);
			}        
		} 	

		else if(opname=="FixRealKCM"){
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int lsbIn = atoi(argv[i++]);
				int msbIn = atoi(argv[i++]);
				int signedInput = checkBoolean(argv[i++], argv[0]);
				int lsbOut = atoi(argv[i++]);
				string constant = argv[i++];
				int useBitheap = checkBoolean(argv[i++], argv[0]);
				op = new FixRealKCM(target, lsbIn, msbIn, signedInput, lsbOut, constant, 1.0, emptyDelayMap, useBitheap);
				addOperator(op);
			}
		}
		
		else if(opname=="FixRealKCMExpert"){ // hidden, for debug
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int lsbIn = atoi(argv[i++]);
				int msbIn = atoi(argv[i++]);
				int signedInput = checkBoolean(argv[i++], argv[0]);
				int lsbOut = atoi(argv[i++]);
				string constant = argv[i++];
				float targetUlpError = atof(argv[i++]);
				op = new FixRealKCM(target, lsbIn, msbIn, signedInput, lsbOut, constant, targetUlpError);
				addOperator(op);
			}        
		}
		else if(opname=="FPRealKCM"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = atoi(argv[i++]);
				int wF = atoi(argv[i++]);
				string constant = argv[i++];
				op = new FPRealKCM(target, wE, wF, constant);
				addOperator(op);
			}        
		}
		else if(opname=="CRFPConstMult"){ 
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else { 
				int wE_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictlyPositive(argv[i++], argv[0]);
				string constant = argv[i++];
				op = new CRFPConstMult(target, wE_in, wF_in, wE_out, wF_out, constant);
				addOperator(op);
			}        
		} 	
		else if(opname=="FPConstMult"){ 
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else { 
				int wE_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_C = atoi(argv[i++]);
				string constant = argv[i++];
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, wF_C, constant);
				addOperator(op);
			}        
		} 	



		//-------------------FUNCTION EVALUATORS----------------------
		else if (opname == "BasicPolyApprox") {
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			double targetAcc = atof(argv[i++]);
			int g =  atof(argv[i++]);
			BasicPolyApprox *toto = new BasicPolyApprox(func, targetAcc, g);
			cout << "Computed degree is " << toto->degree;
			cout << "Accuracy is " << toto->approxErrorBound << " ("<< log2(toto->approxErrorBound) << " bits)";
		}
		else if (opname == "PiecewisePolyApprox") {
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			double targetAcc = atof(argv[i++]);
			int degree =  atof(argv[i++]);
			PiecewisePolyApprox *toto = new PiecewisePolyApprox(func, targetAcc, degree);
		}


		else if (opname == "FixFunctionEval") {
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			FixFunction *f = new FixFunction(argv[i++]);
			double x = strtod(argv[i++], NULL);
			double r = f -> eval(x);
			cout << "\n r=" << r << endl;  
		}

		else if (opname == "FixFunctionByTable") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int lsbI = atoi(argv[i++]);
			int msbO = atoi(argv[i++]);
			int lsbO = atoi(argv[i++]);
			Operator* tg = new FixFunctionByTable(target, func, lsbI, msbO, lsbO);
			addOperator(tg);
		}

		else if (opname == "FixFunctionBySimplePoly") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int lsbI = atoi(argv[i++]);
			int msbO = atoi(argv[i++]);
			int lsbO = atoi(argv[i++]);
			Operator* tg = new FixFunctionBySimplePoly(target, func, lsbI, msbO, lsbO);
			addOperator(tg);
		}

		else if (opname == "TestBench") {
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			if(target->getGlobalOpListRef()->empty()){
				cerr<<"ERROR: TestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0],opname); // and exit
			}
			int n = checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = target->getGlobalOpListRef()->back();
			Operator* op = new TestBench(target, toWrap, n);
			addOperator(op);
			cerr << "To run the simulation using ModelSim, type the following in 'vsim -c':" <<endl;
			cerr << tab << "vdel -all -lib work" <<endl;
			cerr << tab << "vlib work" <<endl;
			cerr << tab << "vcom " << filename <<endl;
			cerr << tab << "vsim " << op->getName() <<endl;
			cerr << tab << "add wave -r *" <<endl;
			cerr << tab << "run " << ((TestBench*)op)->getSimulationTime() <<"ns" << endl;
			cerr << "To run the simulation using gHDL, type the following in a shell prompt:" <<endl;
			string simlibs;
			if(op->getStdLibType()==0 || op->getStdLibType()==-1)
				simlibs="--ieee=synopsys ";
			if(op->getStdLibType()==1)
				simlibs="--ieee=standard ";
			cerr <<  "ghdl -a " << simlibs << "-fexplicit "<< filename <<endl;
			cerr <<  "ghdl -e " << simlibs << "-fexplicit " << op->getName() <<endl;
			cerr <<  "ghdl -r " << simlibs << op->getName() << " --vcd=" << op->getName() << ".vcd" <<endl;
			cerr <<  "gtkwave " << op->getName() << ".vcd" << endl;
		}
		
		else if (opname == "TestBenchFile") {
			/* Using a file to store IO */
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			if(target->getGlobalOpListRef()->empty()){
				cerr<<"ERROR: TestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0],opname); // and exit
			}
			int n = atoi(argv[i++]);//checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = target->getGlobalOpListRef()->back();
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
			string simlibs;
			if(op->getStdLibType()==0 || op->getStdLibType()==-1)
				simlibs="--ieee=synopsys ";
			if(op->getStdLibType()==1)
				simlibs="--ieee=standard ";
			cerr <<  "ghdl -a " << simlibs << "-fexplicit "<< filename <<endl;
			cerr <<  "ghdl -e " << simlibs << "-fexplicit " << op->getName() <<endl;
			cerr <<  "ghdl -r " << simlibs << op->getName() << " --vcd=" << op->getName() << ".vcd" <<endl;
			cerr <<  "gtkwave " << op->getName() << ".vcd" << endl;
		}

		else  {
			cerr << "ERROR: Problem parsing input line, exiting";
			usage(argv[0]);
		}



	} while (i<argc);
	return true;
}

















int main(int argc, char* argv[] )
{
#ifdef HAVE_SOLLYA
	sollya_lib_init();
#endif

	uint32_t i;
	

	target = new Virtex5(); // this also creates a global operator list

	// for historical reasons, to get rid of some day

	try {
		parseCommandLine(argc, argv);
	} catch (char const * s) {
		cerr << "Exception while parsing command line: " << s << endl;
		return 1;
	} catch (std::string s){
		cerr << "Exception while parsing command line: " << s << endl;
		return 1;
	}

	 vector<Operator*>* oplist=target->getGlobalOpListRef();




	ofstream file;
	file.open(filename.c_str(), ios::out);
	Operator::outputVHDLToFile(*oplist, file); 
	file.close();
	
	cerr << endl<<"Final report:"<<endl;
	for(i=0; i<oplist->size(); i++) {
		(*oplist)[i]->outputFinalReport(0);
	}
	
	cerr<< "Output file: " << filename <<endl;
	

	//------------------------ Resource Estimation ---------------------
	for (vector<Operator*>::iterator it = oplist->begin(); it!=oplist->end(); ++it) {
		Operator* op = *it;
		
		if(reLevel!=0){
			if(op->reActive)
				cerr << op->generateStatistics(reLevel);
			else{
				cerr << "Resource estimation option active for an operator that has NO estimations in place." << endl;
			}
		}
	}
	//------------------------------------------------------------------
	
	//------------------ Resource Estimation Debugging -----------------
	if(resourceEstimationDebug){
		ofstream file;
		file.open("flopoco.debug");
		for (vector<Operator*>::iterator it = oplist->begin(); it!=oplist->end(); ++it) {
			Operator* op = *it;
			
			if(op->reActive)
				file << op->resourceEstimate.str();
			else{
				cerr << "Resource estimation debugging option active for an operator that has NO estimations in place." << endl;
			}
		}
		file.close();
		cerr << "Resource estimation log written to the \'flopoco.debug\' file" << endl;
	}
	//------------------------------------------------------------------
	
	//------------------------ Floorplanning ---------------------------
	if(floorplanning){
		for (vector<Operator*>::iterator it = oplist->begin(); it!=oplist->end(); ++it) {
			Operator* op = *it;
			
			if(op->reActive == false){
				cerr << "Floorplanning option can be used only when the resource estimations have been performed.\n" 
						<< " Please reconsider your strategy.\n" << endl;
				exit(1);
			}
			if(floorplanning)
				floorplanMessages << op->createFloorplan();
		}
	}
	//------------------------------------------------------------------
	
	//------------------ Floorplanning Debugging -----------------------
	if(floorplanningDebug){
		if(!(reLevel>0 && floorplanning)){
			cerr << "Debugging Floorplanning option can be used only when there are resource estimations and floorplanning is enabled." 
						<< " Please rerun the program with the appropriate options" << endl;
				exit(1);
		}else{
			ofstream file;
			file.open("flopoco.floorplan.debug");
			for (vector<Operator*>::iterator it = oplist->begin(); it!=oplist->end(); ++it) {
				Operator* op = *it;
				
				file << op->floorplan.str();
				file << floorplanMessages.str();
			}
			file.close();
			cerr << "Floorplanning log (for debugging purposes) written to the \'flopoco.floorplanning.debug\' file" << endl;
		}
	}


	//------------------------------------------------------------------

#ifdef HAVE_SOLLYA
	sollya_lib_close();
#endif

	return 0;
}



