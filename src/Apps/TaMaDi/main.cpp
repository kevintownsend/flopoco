/*
  the FloPoCo command-line interface
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Authors : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr
            Bogdan Pasca, Bogdan.Pasca@ens-lyon.org

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <mpfr.h>
#include <cstdlib>
#include "TaMaDi.hpp"
#include "main.hpp"

//#define BRIGHT 1
//#define RED 31
//#define OPER 32
//#define NEWOPER 32
//#define PARAM 34
//#define OP(op,paramList)             {cerr << "    "; printf("%c[%d;%dm",27,1,OPER); cerr <<  op; printf("%c[%dm",27,0); cerr<< " "; printf("%c[%d;%dm",27,1,PARAM); cerr << paramList; printf("%c[%dm\n",27,0); } 
//#define NEWOP(op,paramList)          {cerr << "    "; printf("%c[%d;%dm",27,1,NEWOPER); cerr <<  op; printf("%c[%dm",27,0); cerr<< " "; printf("%c[%d;%dm",27,1,PARAM); cerr << paramList; printf("%c[%dm\n",27,0); } 


using namespace std;
using namespace flopoco;

// Global variables, useful in this main to avoid parameter passing


	string filename="flopoco.vhdl";
	string cl_name=""; // used for the -name option
	Target* target;



static void usage(char *name, string opName = ""){
	bool full = (opName=="");

	if ( full ){
		cerr << "\nUsage: "<<name<<" <operator specification list>\n" ;
		cerr << "Each operator specification is one of: \n";
	}

	if ( full || opName == "TaMaDiModule" || opName == "TaMaDi"){
		OP("TaMaDiModule","wP degree nbOfIterations IntervalIDWidth widthComp nbOfPE inFiFoDepth PEFiFoDepth outFiFoDepth");;
		cerr << "       Assembles a TaMaDiMultimodule formed of multiple TaMaDi cores\n";
		cerr << "       wP             : internal fixed-point working precision \n";
		cerr << "       degree         : polyonomial degree used in the TaMaDiCore \n";
		cerr << "       nbOfIterations : the number of allowed iterations per TaMaDiCore\n";
		cerr << "       IntervalIDWidth: the number of bits required for storing the interval ID\n";
		cerr << "       widthComp      : the width of the pattern detector \n";
		cerr << "       nbOfPE         : the number of TaMaDiCores instantiated in the multi-module \n";
		cerr << "       inFiFoDepth    : the multi-module's input FIFO depth \n";
		cerr << "       PEFiFoDepth    : each TaMaDiCore is associated with its own output FIFO; \n";
		cerr << "                        this argument sets its depth; \n";
		cerr << "       outFiFoDepth   : the multi-module's output FIFO depth \n";
	}	
	if ( full || opName == "TaMaDiModuleDummyWrapper"|| opName == "TaMaDi"){
		OP("TaMaDiModuleDummyWrapper","wP degree nbOfIterations IntervalIDWidth widthComp nbOfPE inFiFoDepth PEFiFoDepth outFiFoDepth");;
		cerr << "       Assembles a TaMaDiMultimodule and adds a serial interface\n";
	}	
	if ( full || opName == "TaMaDiDeserializer"|| opName == "TaMaDi"){
		OP("TaMaDiDeserializer","wP degree nbOfIterations IntervalIDWidth widthComp nbOfPE inFiFoDepth PEFiFoDepth outFiFoDepth");;
		cerr << "       Module used to deserialize the data from the DMA FIFO into initialization packets\n";
	}	
	if ( full || opName == "TaMaDiModuleWrapperInterface"|| opName == "TaMaDi"){
		OP("TaMaDiModuleWrapperInterface","wP degree nbOfIterations IntervalIDWidth widthComp nbOfPE inFiFoDepth PEFiFoDepth outFiFoDepth");;
		cerr << "       Adds the counter-based interface for the credit-based dispatch scheme\n";
	}	
	if ( full || opName == "TaMaDiModuleDispatcherInterface"|| opName == "TaMaDi"){
		OP("TaMaDiModuleDispatcherInterface","wP degree nbOfIterations IntervalIDWidth widthComp nbOfPE inFiFoDepth PEFiFoDepth outFiFoDepth InterfaceInFiFoDepth InterfaceOutFiFoDepth");;
		cerr << "       The Dispatcher module which dispatches initialization packets to the interface-wrapped TaMaDiMultiModules\n";
		cerr << "       InterfaceInFiFoDepth    : the Dispatcher's input FIFO depth (will be filled-in by the Deserializer) \n";
		cerr << "       InterfaceOutFiFoDepth   : the Dispatcher's output FIFO depth (will be filled by the TaMaDiMultimodules with HR cases\n";
	}	
	if ( full || opName == "TaMaDiSystem"|| opName == "TaMaDi"){
		OP("TaMaDiSystem","wP degree nbOfIterations IntervalIDWidth widthComp nbOfPE inFiFoDepth PEFiFoDepth outFiFoDepth InterfaceInFiFoDepth InterfaceOutFiFoDepth nbOfModules");;
		cerr << "       Assembles a full system with a Deserializer, Dispatcher and several multi-modules\n";
		cerr << "       wP             : internal fixed-point working precision \n";
		cerr << "       degree         : polyonomial degree used in the TaMaDiCore \n";
		cerr << "       nbOfIterations : the number of allowed iterations per TaMaDiCore\n";
		cerr << "       IntervalIDWidth: the number of bits required for storing the interval ID\n";
		cerr << "       widthComp      : the width of the pattern detector \n";
		cerr << "       nbOfPE         : the number of TaMaDiCores instantiated in the multi-module \n";
		cerr << "       inFiFoDepth    : the multi-module's input FIFO depth \n";
		cerr << "       PEFiFoDepth    : each TaMaDiCore is associated with its own output FIFO; \n";
		cerr << "                        this argument sets its depth; \n";
		cerr << "       outFiFoDepth   : the multi-module's output FIFO depth \n";
		cerr << "       InterfaceInFiFoDepth    : the Dispatcher's input FIFO depth (will be filled-in by the Deserializer) \n";
		cerr << "       InterfaceOutFiFoDepth   : the Dispatcher's output FIFO depth (will be filled by the TaMaDiMultimodules with HR cases\n";
		cerr << "       nbOfModules             : the number of multi-modules instantiated in the TaMaDiSystem"<<endl;
	}	
	if ( full || opName == "TaMaDiCore" || opName == "TaMaDi"){
		OP("TaMaDiCore","wP degree nbOfIterations IntervalIDWidth widthComp");;
		cerr << "       The TaMaDiCore performing the polynomial evaluation based on the tabulated differences\n";
	}	

	if ( full || opName=="TestBench" ){
		cerr << "    ____________ TEST-BENCH ____________________________________________________\n";
		OP ("TestBench","n");
		cerr << "       Behavorial test bench for the preceding operator\n";
		cerr << "       This test bench will include standard tests, plus n random tests.\n";
		OP( "TestBenchFile","n");
		cerr << "       Behavorial test bench for the preceding operator\n";
		cerr << "       This test bench will include standard tests, plus n random tests.\n";
		cerr << "       Inputs and outputs are stored in a file to reduce VHDL compilation time.\n";
		cerr << "       if n=-2, an exhaustive test is generated (use only for small operators).\n";
	}
	
	if ( full || opName=="extra"){
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
		cerr << "   -target=<Spartan3|Virtex4|Virtex5|StratixII|StratixIII|StratixIV>      (default=Virtex5)\n";
		cerr << "   -DSP_blocks=<yes|no>\n";
		cerr << "       optimize for the use of DSP blocks   (default=yes)\n";
		cerr << "   -name=<entity name>\n";
		cerr << "       defines the name of the VHDL entity of the next operator\n";
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


void addOperator(vector<Operator*> &oplist, Operator *op) {
	if(cl_name!="")	{
		cerr << "Updating entity name to: " << cl_name << endl;
		op->changeName(cl_name);
		cl_name="";
	}
	// TODO check name not already in list...
	// TODO this procedure should be a static method of operator, so that other operators can call it
	oplist.push_back(op);
}



bool parseCommandLine(int argc, char* argv[], vector<Operator*> &oplist){
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
						usage(argv[0], "extra");
					}
				}
				else if (o == "target") {
					Target* oldTarget=target;
					if(v=="Virtex4") target=new Virtex4();
					else if (v=="Virtex6") target=new Virtex6();
					else if (v=="Virtex5") target=new Virtex5();
					else if (v=="Spartan3") target=new Spartan3();
					else if (v=="StratixII") target=new StratixII();
					else if (v=="StratixIII") target=new StratixIII();
					else if (v=="StratixIV") target=new StratixIV();
					else {
						cerr<<"ERROR: unknown target: "<<v<<endl;
						usage(argv[0],"extra");
					}
					// if previous options had changed it
					target->setFrequency(oldTarget->frequency());
					target->setUseHardMultipliers(oldTarget->getUseHardMultipliers());
					if (oldTarget->isPipelined()) 
						target->setPipelined();
					else 
						target->setNotPipelined();
				}

				else if (o == "pipeline") {
					if(v=="yes") target->setPipelined();
					else if(v=="no")  target->setNotPipelined();
					else {
						cerr<<"ERROR: pipeline option should be yes or no,    got "<<v<<"."<<endl; 
						usage(argv[0],"extra");
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
						usage(argv[0],"extra");
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
				usage(argv[0],opname);
			else {
				int w = atoi(argv[i++]);
				mpz_class mpc(argv[i++]);
				int signedInput = checkBoolean(argv[i++], argv[0]);

				cerr << "> IntIntKCM , w="<<w<<", c="<<mpz2string(mpc) << (signedInput==0?" unsigned":" signed") << "\n";
				op = new IntIntKCM(target, w, mpc, signedInput);
				addOperator(oplist, op);
			}        
		}
		else if(opname=="IntConstMult"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int w = atoi(argv[i++]);
				mpz_class mpc(argv[i++]);
				cerr << "> IntConstMult , w="<<w<<", c="<<mpz2string(mpc)<<"\n";
				op = new IntConstMult(target, w, mpc);
				addOperator(oplist, op);
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
				addOperator(oplist, op);
			} 
		}


		else if(opname=="FPConstMultRational"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_in = checkStrictlyPositive(argv[i++], argv[0]);
				int wE_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictlyPositive(argv[i++], argv[0]);
				int a  = atoi(argv[i++]); 
				int b  = checkStrictlyPositive(argv[i++], argv[0]); 
				cerr << "> FPConstMultRational, wE_in="<<wE_in<<", wF_in="<<wF_in
						 <<", wE_out="<<wE_out<<", wF_out="<<wF_out
						 <<", constant="<<a<<"/"<<b<<endl;
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, a, b);
				addOperator(oplist, op);
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
				cerr << "> FPConstMultExpert, wE_in="<<wE_in<<", wF_in="<<wF_in
						 <<", wE_out="<<wE_out<<", wF_out="<<wF_out
						 <<", cst_sgn="<<cst_sgn<<", cst_exp="<<cst_exp<< ", cst_sig="<<mpz2string(cst_sig)<<endl;
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, cst_sgn, cst_exp, cst_sig);
				addOperator(oplist, op);
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
				addOperator(oplist, op);
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
				addOperator(oplist, op);
			}        
		} 	


#ifdef HAVE_SOLLYA
		else if(opname=="FixRealKCM"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
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
				addOperator(oplist, op);
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
				cerr << "> FixRealKCM, msbIn="<<msbIn<<", lsbIn="<<lsbIn 
					  <<", lsbOut="<<lsbOut
					  << ", constant="<<constant <<endl;
				op = new FixRealKCM(target, lsbIn, msbIn, signedInput, lsbOut, constant, targetUlpError);
				addOperator(oplist, op);
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
				cerr << "> FPRealKCM, wE="<<wE<<", wF="<<wF << ", constant="<<constant <<endl;
				op = new FPRealKCM(target, wE, wF, constant);
				addOperator(oplist, op);
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
				cerr << "> CRFPConstMult, wE_in="<<wE_in<<", wF_in="<<wF_in 
					  <<", wE_out="<<wE_out<<", wF_out="<<wF_out
					  << ", constant="<<constant <<endl;
				op = new CRFPConstMult(target, wE_in, wF_in, wE_out, wF_out, constant);
				addOperator(oplist, op);
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
				cerr << "> FPConstMult, wE_in="<<wE_in<<", wF_in="<<wF_in 
					  <<", wE_out="<<wE_out<<", wF_out="<<wF_out 
					  << ", wF_C=" << wF_C
					  << ", constant="<<constant <<endl;
				op = new FPConstMult(target, wE_in, wF_in, wE_out, wF_out, wF_C, constant);
				addOperator(oplist, op);
			}        
		} 	
#endif // HAVE_SOLLYA
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
				cerr << "> LeftShifter, wIn="<<wIn<<", maxShift="<<maxShift<<"\n";
				op = new Shifter(target, wIn, maxShift, Shifter::Left);
				addOperator(oplist, op);
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
				cerr << "> RightShifter, wIn="<<wIn<<", maxShift="<<maxShift<<"\n";
				op = new Shifter(target, wIn, maxShift, Shifter::Right, inputDelays);
				addOperator(oplist, op);
			}
		}
		else if(opname=="LZOC"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> LZOC, wIn="<<wIn<<"\n";
				op = new LZOC(target, wIn);
				addOperator(oplist, op);
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
					cerr << "> LZOCShifter, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, -1);
					addOperator(oplist, op);
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
					cerr << "> LZCShifter, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, 0);
					addOperator(oplist, op);
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
					cerr << "> LOCShifter, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 0, 1);
					addOperator(oplist, op);
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
					cerr << "> LZOCShifterSticky, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, -1);
					addOperator(oplist, op);
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
					cerr << "> LZCShifterSticky, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, 0);
					addOperator(oplist, op);
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
					cerr << "> LOCShifterSticky, wIn="<<wIn<<", wOut="<<wOut<< endl;
					op = new LZOCShifterSticky(target, wIn, wOut, intlog2(wIn), 1, 1);
					addOperator(oplist, op);
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
				usage(argv[0],opname);
			else {
				int param0 = checkStrictlyPositive(argv[i++], argv[0]);
				int param1 = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> UserDefinedOperator, param0="<<param0<<", param1=" << param1 << endl  ;
                                op = new UserDefinedOperator(target,param0,param1);
                                addOperator(oplist, op);
			}    
		}
		else if(opname=="IntAdder"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0], opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> IntAdder, wIn="<<wIn<<endl  ;
				op = new IntAdder(target,wIn, inDelayMap("X",target->ffDelay() + target->localWireDelay()) );
				addOperator(oplist, op);
			}    
		}
		else if(opname=="IntComparator"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int criteria = atoi(argv[i++]);
				
				cerr << "> IntComparator, wIn="<<wIn<<" criteria="<<criteria<<endl  ;
				op = new IntComparator(target,wIn,criteria,false,0);
				addOperator(oplist, op);
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
				
				cerr << "> IntComparator, wIn="<<wIn<<" criteria="<<criteria<<endl  ;
				op = new IntComparator(target,wIn,criteria, true, constant);
				addOperator(oplist, op);
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

				cerr << "> IntAdder, wIn="<<wIn<<", frequency="<<target->frequency()<< " inputDelay="<<inputDelay<< endl  ;

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
				addOperator(oplist, op);
			}    
		}

		/* Multioperand adders */

		else if(opname=="IntMultiAdder"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int N   = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> IntMultiAdder, wIn="<<wIn<<" N="<<N<<endl  ;
				op = new IntMultiAdder(target,wIn,N);
				addOperator(oplist, op);
			}    
		}

		else if(opname=="IntNAdder"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int N   = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> IntNAdder, wIn="<<wIn<<" N="<<N<<endl  ;
				op = new IntNAdder(target, wIn, N, inDelayMap("X0",target->localWireDelay()));
				addOperator(oplist, op);
			}    
		}

		else if(opname=="IntCompressorTree"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int N   = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> IntCompressorTree, wIn="<<wIn<<" N="<<N<<endl  ;
				op = new IntCompressorTree(target,wIn,N);
				addOperator(oplist, op);
			}    
		}

		/* Exploration of other fast adders */
		else if(opname=="IntAdderSpecific"){ //Hidden
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> IntAdderSpecific, wIn="<<wIn<<endl  ;
				op = new IntAdderSpecific(target,wIn);
				addOperator(oplist, op);
			}    
		}
		else if(opname=="IntComparatorSpecific"){ //Hidden
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int type = atoi(argv[i++]);
				cerr << "> IntComparatorSpecific, wIn="<<wIn<<" type="<<type<< endl  ;
				op = new IntComparatorSpecific(target,wIn,type);
				addOperator(oplist, op);
			}    
		}
		else if(opname=="CarryGenerationCircuit"){ //Hidden
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> CarryGenerationCircuit, wIn="<<wIn<<endl  ;
				op = new CarryGenerationCircuit(target,wIn);
				addOperator(oplist, op);
			}    
		}
		/* interface */
		else if(opname=="LongIntAdderAddAddMux"){ //AAM
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int g   = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdderSpecific, wIn="<<wIn<<endl  ;
				if (g==1)
					op = new LongIntAdderAddAddMuxGen1(target,wIn);
				else if (g==2)
					op = new LongIntAdderAddAddMuxGen2(target,wIn, inDelayMap("X",target->ffDelay() + target->localWireDelay() ));
				else 
					throw "Generation parameter is either 1 or 2";
				addOperator(oplist, op);
			}    
		}
		else if(opname=="LongIntAdderMuxNetwork"){ //Mux network
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdderMuxNetwork, wIn="<<wIn<<endl  ;
				op = new LongIntAdderMuxNetwork(target,wIn);
				addOperator(oplist, op);
			}    
		}
		else if(opname=="LongIntAdderCmpCmpAdd"){ //CCA
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int g   = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdderCmpCmpAddSpecific, wIn="<<wIn<<endl  ;
				if (g==1)
					op = new LongIntAdderCmpCmpAddGen1(target,wIn);
				else if (g==2)
					op = new LongIntAdderCmpCmpAddGen2(target,wIn, inDelayMap("X",target->ffDelay() + target->localWireDelay() ));
				else 
					throw "Generation parameter is either 1 or 2";
				addOperator(oplist, op);
			}    
		}
		else if(opname=="LongIntAdderCmpAddInc"){ //CAI
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int g   = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> LongIntAdderCmpAddInc, wIn="<<wIn<<endl  ;
				if (g==1)
					op = new LongIntAdderCmpAddIncGen1(target,wIn);
				else if (g==2)
					op = new LongIntAdderCmpAddIncGen2(target,wIn, inDelayMap("X",target->ffDelay() + target->localWireDelay() ));
				else 
					throw "Generation parameter is either 1 or 2";
				addOperator(oplist, op);
			}    
		}
		/*---------------------------------------------------------*/
		
		//HIDDEN
		else if(opname=="IntDualSub"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				int opType = checkBoolean(argv[i++], argv[0]);
				cerr << "> IntDualSub, wIn="<<wIn<<" OpType="<<opType<<endl  ;
				op = new IntDualSub(target,wIn,opType);
				addOperator(oplist, op);
			}    
		}
		else if(opname=="IntMultiplier"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wInX    = checkStrictlyPositive(argv[i++], argv[0]);
				int wInY    = checkStrictlyPositive(argv[i++], argv[0]);
				int sign    =  checkBoolean(argv[i++], argv[0]);
				float ratio = atof(argv[i++]);
				bool sign_;
				if (sign > 0)
					sign_ = true;
				else
					sign_ = false;
				cerr << "> IntMultiplier , wInX="<<wInX<<", wInY="<<wInY<<" signed="<<sign_<<"\n";
				op = new IntMultiplier(target, wInX, wInY, emptyDelayMap, sign_, ratio);
				addOperator(oplist, op);
			}
		}
		else if(opname=="UnsignedIntMultiplier"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wInX = checkStrictlyPositive(argv[i++], argv[0]);
				int wInY = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> UnsignedIntMultiplier , wInX="<<wInX<<", wInY="<<wInY<<"\n";
				op = new UnsignedIntMultiplier(target, wInX, wInY);
				addOperator(oplist, op);
			}
		}
		else if(opname=="SignedIntMultiplier"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wInX = checkStrictlyPositive(argv[i++], argv[0]);
				int wInY = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> SignedIntMultiplier , wInX="<<wInX<<", wInY="<<wInY<<"\n";
				op = new SignedIntMultiplier(target, wInX, wInY);
				addOperator(oplist, op);
			}
		}
		else if(opname=="IntKaratsuba"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wIn = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> IntKaratsuba , wIn="<<wIn<<"\n";
				op = new IntKaratsuba(target, wIn);
				addOperator(oplist, op);
			}    
		}   
		else if(opname=="IntTilingMultiplier"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wInX = checkStrictlyPositive(argv[i++], argv[0]);
				int wInY = checkStrictlyPositive(argv[i++], argv[0]);
				float r = atof(argv[i++]);
				if (r<0)
					throw "Ratio > 0";
				int maxTimeInMinutes = atoi(argv[i++]);
				
				cerr << "> IntTilingMultiplier , wInX="<<wInX<<", wInY="<<wInY<<" ratio=" << r << " maxTimeInMinutes= "<< maxTimeInMinutes << "\n";
				op = new IntTilingMult(target, wInX, wInY, r, maxTimeInMinutes, true);
				addOperator(oplist, op);
			}
		}
		else if(opname=="IntTruncMultiplier"){
			int nargs = 7;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wInX = checkStrictlyPositive(argv[i++], argv[0]);
				int wInY = checkStrictlyPositive(argv[i++], argv[0]);
				int wOut = checkStrictlyPositive(argv[i++], argv[0]);
				float r = atof(argv[i++]);
				// TODO I think it is more consistent to move this kind of tests insinde the operator
				if (r<0)
					throw "Ratio must be > 0";
				int useLimits = atoi(argv[i++]);
				if (!((useLimits==0)||(useLimits==1)))
					throw "useLimits is 0 or 1";				
				int maxTimeInMinutes = atoi(argv[i++]);
				int sign = atoi(argv[i++]);

				cerr << "> IntTruncMult , wInX="<<wInX<<", wInY=" <<wInY <<", wOut=" <<wOut << " ratio=" << r << " useLimits="<<useLimits<<" maxTimeInMinutes="<<maxTimeInMinutes<<" signed="<<sign<<"\n";
				op = new IntTruncMultiplier(target, wInX, wInY, wOut, r,useLimits, maxTimeInMinutes, false, sign);
				addOperator(oplist, op);
			}
		}		
		// For the FPAdd the default is the single-path design
		else if(opname=="FPAdd"){ 
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> FPAdd , wE="<<wE<<", wF="<<wF<<" \n";
				op = new FPAddSinglePath(target, wE, wF, wE, wF, wE, wF);
				addOperator(oplist, op);
			}
		}		
		else if(opname=="FPAddDualPath"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> FPAddDualPath , wE="<<wE<<", wF="<<wF<<" \n";
				op = new FPAddDualPath(target, wE, wF, wE, wF, wE, wF);
				addOperator(oplist, op);
			}
		}
		else if(opname=="FPAdd3Input"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> FPAdd3Input , wE="<<wE<<", wF="<<wF<<" \n";
				op = new FPAdd3Input(target, wE, wF);
				addOperator(oplist, op);
			}
	}	
		else if(opname=="FPFMAcc"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				int adderLatency = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> FPFMAcc , wE="<<wE<<", wF="<<wF<<" adderLatency="<<adderLatency<< " \n";
				op = new FPFMAcc(target, wE, wF, adderLatency);
				addOperator(oplist, op);
			}
		}		
#ifdef HAVE_SOLLYA
		else if(opname=="FPJacobi"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				int l0 = checkStrictlyPositive(argv[i++], argv[0]);
				int l1 = checkStrictlyPositive(argv[i++], argv[0]);
				int l2 = checkStrictlyPositive(argv[i++], argv[0]);
				int ver = atoi(argv[i++]);
				cerr << "> FPJacobi , wE="<<wE<<", wF="<<wF<<" l0="<<l0<<" l1="<<l1<<" l2="<<l2<<" ver="<<ver<< " \n";
				op = new FPJacobi(target, wE, wF, l0, l1, l2,ver);
				addOperator(oplist, op);
			}
	}	
#endif // HAVE_SOLLYA
		else if(opname=="TaMaDiModule"){
			int nargs = 9;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID  =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]); //number of pe
				int inFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int peFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int outFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> TaMaDiModule , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID<< " width comparisson="<<widthComp<<" n="<<n<<" inFifo="<<inFifoDepth<<" peFifo="<<peFifoDepth<<" outFifo"<<outFifoDepth <<" \n";
				op = new TaMaDiModule(target, wP, degree, numberOfIterations, widthOfIntervalID, widthComp, n, inFifoDepth, peFifoDepth, outFifoDepth);
				addOperator(oplist, op);
			}
		}	
		else if(opname=="TaMaDiModuleDummyWrapper"){
			int nargs = 9;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]); //number of pe
				int inFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int peFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int outFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> TaMaDiModuleDummyWrapper , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID<<" width comparisson="<<widthComp<<" n="<<n<<" inFifo="<<inFifoDepth<<" peFifo="<<peFifoDepth<<" outFifo"<<outFifoDepth <<" \n";
				op = new TaMaDiModuleDummyWrapper(target, wP, degree, numberOfIterations, widthOfIntervalID, widthComp, n, inFifoDepth, peFifoDepth, outFifoDepth);
				addOperator(oplist, op);
			}
		}
		else if(opname=="TaMaDiDeserializer"){
			int nargs = 9;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]); //number of pe
				int inFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int peFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int outFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> TaMaDiDeserializer , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID<<" width comparisson="<<widthComp<<" n="<<n<<" inFifo="<<inFifoDepth<<" peFifo="<<peFifoDepth<<" outFifo"<<outFifoDepth <<" \n";
				op = new TaMaDiDeserializer(target, wP, degree, numberOfIterations, widthOfIntervalID, n, inFifoDepth, peFifoDepth, outFifoDepth);
				addOperator(oplist, op);
			}
		}		
		else if(opname=="TaMaDiModuleWrapperInterface"){
			int nargs = 9;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]); //number of pe
				int inFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int peFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int outFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> TaMaDiModuleWrapperInterface , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID<<" width comparisson="<<widthComp<<" n="<<n<<" inFifo="<<inFifoDepth<<" peFifo="<<peFifoDepth<<" outFifo"<<outFifoDepth <<" \n";
				op = new TaMaDiModuleWrapperInterface(target, wP, degree, numberOfIterations, widthOfIntervalID, widthComp, n, inFifoDepth, peFifoDepth, outFifoDepth);
				addOperator(oplist, op);
			}
		}
		else if(opname=="TaMaDiDispatcherInterface"){
			int nargs = 11;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]); //number of pe
				int inFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int peFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int outFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int InterfaceInFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int InterfaceOutFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> TaMaDiModuleWrapperInterface , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID<<" width comparisson="<<widthComp<<" n="<<n<<" inFifo="<<inFifoDepth<<" peFifo="<<peFifoDepth<<" outFifo"<<outFifoDepth << " interfaceInFifo="<<InterfaceInFifoDepth<<" interfaceOutFifoDepth="<<InterfaceOutFifoDepth<<"\n";
				op = new TaMaDiDispatcherInterface(target, wP, degree, numberOfIterations, widthOfIntervalID, widthComp,  n, inFifoDepth, peFifoDepth, outFifoDepth, InterfaceInFifoDepth, InterfaceOutFifoDepth);
				addOperator(oplist, op);
			}
		}
		else if(opname=="TaMaDiSystem"){
			int nargs = 12;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]); //number of pe
				int inFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int peFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int outFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int InterfaceInFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int InterfaceOutFifoDepth = checkStrictlyPositive(argv[i++], argv[0]);
				int moduleCount = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> TaMaDiSystem , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID<<" width comparisson="<<widthComp<<" n="<<n<<" inFifo="<<inFifoDepth<<" peFifo="<<peFifoDepth<<" outFifo"<<outFifoDepth << " interfaceInFifo="<<InterfaceInFifoDepth<<" interfaceOutFifoDepth="<<InterfaceOutFifoDepth<<" moduleCount="<<moduleCount<<"\n";
				op = new TaMaDiSystem(target, wP, degree, numberOfIterations, widthOfIntervalID, widthComp, n, inFifoDepth, peFifoDepth, outFifoDepth, InterfaceInFifoDepth, InterfaceOutFifoDepth, moduleCount );
				addOperator(oplist, op);
			}
		}
		else if(opname=="TaMaDiCore"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wP = checkStrictlyPositive(argv[i++], argv[0]);
				int degree = checkStrictlyPositive(argv[i++], argv[0]);
				int numberOfIterations = checkStrictlyPositive(argv[i++], argv[0]);
				int widthOfIntervalID =  checkStrictlyPositive(argv[i++], argv[0]);
				int widthComp          =  checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> TaMaDiCore , wP="<<wP<<", degree="<<degree<<" numberOfIterations="<<numberOfIterations<<" widthOfIntervalID="<<widthOfIntervalID <<" width comparisson="<<widthComp<<" \n";
				op = new TaMaDiCore(target, wP, degree, numberOfIterations, widthOfIntervalID, widthComp);
				addOperator(oplist, op);
			}
		}	
		else if(opname=="TaMaDiPriorityEncoder"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int n = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> TaMaDiPriorityEncoder , n="<<n<<" \n";
				op = new TaMaDiPriorityEncoder(target, n);
				addOperator(oplist, op);
			}
		}	
		else if(opname=="TaMaDiDecoder"){
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int n = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> TaMaDiDecoder ,output n="<<n<<" \n";
				op = new TaMaDiDecoder(target, n);
				addOperator(oplist, op);
			}
		}	
		else if(opname=="TaMaDiFIFO"){
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int w = checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]);
				int limit = atoi(argv[i++]);
				cerr << "> TaMaDiFIFO ,output w="<<n<<" depth="<<n<<" with full report limit at "<<limit<<" \n";
				op = new TaMaDiFIFO(target, w, n, limit);
				addOperator(oplist, op);
			}
		}
		
		else if(opname=="TaMaDiShiftRegister"){
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int w = checkStrictlyPositive(argv[i++], argv[0]);
				int n = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> TaMaDiShiftRegister ,w="<<w<<" levels="<<n<<" \n";
				op = new TaMaDiShiftRegister(target, w, n);
				addOperator(oplist, op);
			}
		}		

		else if(opname=="Fix2FP"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int LSB = atoi(argv[i++]);//checkStrictlyPositive(argv[i++], argv[0]);
				int MSB = atoi(argv[i++]);//checkStrictlyPositive(argv[i++], argv[0]);
				int sign = atoi(argv[i++]);
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wF = checkStrictlyPositive(argv[i++], argv[0]);
				
				cerr << "> Fix2FP, LSB="<<LSB<<", MSB="<<MSB<<", wE="<<wE<<", wF="<<wF<<" \n";
				op = new Fix2FP(target, LSB, MSB, sign,wE, wF);
				addOperator(oplist, op);
			}
		}
		// else if(opname=="CoilInductance"){
		// 	int nargs = 7;
		// 	if (i+nargs > argc)
		// 		usage(argv[0],opname);
		// 	else {
		// 		int LSBI = atoi(argv[i++]);
		// 		int MSBI = atoi(argv[i++]);
		// 		int wE = checkStrictlyPositive(argv[i++], argv[0]);
		// 		int wF = checkStrictlyPositive(argv[i++], argv[0]);				
		// 		int MaxMSBO= atoi(argv[i++]);
		// 		int LSBO = atoi(argv[i++]);
		// 		int MSBO = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoilInductance  LSBI="<<LSBI<<", MSBI="<<MSBI<<",wEIn="<<wE<<",wFIn"<<wF<<", MaxMSBO="<<MaxMSBO<<", LSBO="<<LSBO<<", MSBO="<<MSBO<<" \n";
		// 		op = new CoilInductance(target, LSBI, MSBI,wE,wF,MaxMSBO,LSBO,MSBO,pa);
		// 		addOperator(oplist, op);
		// 	}
		// }
		// else if(opname=="CoordinatesTableX"){
		// 	int nargs = 4;
		// 	if (i+nargs > argc)
		// 		usage(argv[0],opname);
		// 	else {
		// 		int wIn = checkStrictlyPositive(argv[i++], argv[0]);
		// 		int LSB = atoi(argv[i++]);
		// 		int MSB = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoordinatesTableX, wIn="<<wIn<<", LSB="<<LSB<<", MSB="<<MSB<<" \n";
		// 		op = new CoordinatesTableX(target, wIn,LSB, MSB,pa);
		// 		addOperator(oplist, op);
		// 	}
		// }
		// else if(opname=="CoordinatesTableZ"){
		// 	int nargs = 4;
		// 	if (i+nargs > argc)
		// 		usage(argv[0],opname);
		// 	else {
		// 		int wIn = checkStrictlyPositive(argv[i++], argv[0]);
		// 		int LSB = atoi(argv[i++]);
		// 		int MSB = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoordinatesTableZ, wIn="<<wIn<<", LSB="<<LSB<<", MSB="<<MSB<<" \n";
		// 		op = new CoordinatesTableZ(target, wIn,LSB, MSB,pa);
		// 		addOperator(oplist, op);
		// 	}
		// }
		// else if(opname=="CoordinatesTableY"){
		// 	int nargs = 4;
		// 	if (i+nargs > argc)
		// 		usage(argv[0],opname);
		// 	else {
		// 		int wIn = checkStrictlyPositive(argv[i++], argv[0]);
		// 		int LSB = atoi(argv[i++]);
		// 		int MSB = atoi(argv[i++]);
		// 		char *pa=argv[i++];
		// 		cerr << "> CoordinatesTableY, wIn="<<wIn<<", LSB="<<LSB<<", MSB="<<MSB<<" \n";
		// 		op = new CoordinatesTableY(target, wIn,LSB, MSB,pa);
		// 		addOperator(oplist, op);
		// 	}
		// }
		else if(opname=="FPMult"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0],opname);
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wFIn = checkStrictlyPositive(argv[i++], argv[0]);
			int wFOut = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPMult , wE="<<wE<<", wFIn="<<wFIn<<", wFOut="<<wFOut<<" \n";
			op = new FPMult(target, wE, wFIn, wE, wFIn, wE, wFOut, true /*normd*/, true /*CR*/);
			addOperator(oplist, op);
		} 
		else if(opname=="FPMultFaithful"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0],opname);
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wFIn = checkStrictlyPositive(argv[i++], argv[0]);
			int wFOut = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPMultFaithful , wE="<<wE<<", wFIn="<<wFIn<<", wFOut="<<wFOut<<" \n";
			op = new FPMult(target, wE, wFIn, wE, wFIn, wE, wFOut, true, false);
			addOperator(oplist, op);
		}
		else if(opname=="FPMultKaratsuba"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wFIn = checkStrictlyPositive(argv[i++], argv[0]);
				int wFOut = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> FPMultKaratsuba , wE="<<wE<<", wFIn="<<wFIn<<", wFOut="<<wFOut<<" \n";
				op = new FPMultKaratsuba(target, wE, wFIn, wE, wFIn, wE, wFOut, 1);
				addOperator(oplist, op);
			}
		}  
		else if(opname=="FPMultExpert"){
			int nargs = 7; 
			if (i+nargs > argc)
				usage(argv[0],opname);
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wFXIn = checkStrictlyPositive(argv[i++], argv[0]);
			int wFYIn = checkStrictlyPositive(argv[i++], argv[0]);
			int wFOut = checkStrictlyPositive(argv[i++], argv[0]);
			int correctRounding = checkBoolean(argv[i++], argv[0]);
			float r = atof(argv[i++]);
			int maxTimeInMinutes = atoi(argv[i++]);

			op = new FPMult(target, wE, wFXIn, wE, wFYIn, wE, wFOut, true, correctRounding, r, maxTimeInMinutes);
			addOperator(oplist, op);
		}  
		else if(opname=="FPSquare"){
			int nargs = 3; 
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wFX = checkStrictlyPositive(argv[i++], argv[0]);
				int wFR = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> FPSquare , wE="<<wE<<", wFX="<<wFX<<" wFR="<<wFR<< " \n";
				op = new FPSquare(target, wE, wFX, wFR);
				addOperator(oplist, op);
			}
		} 
		else if (opname == "FPDiv")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPDiv: wE=" << wE << " wF=" << wF << endl;
			op = new FPDiv(target, wE, wF);
			addOperator(oplist, op);
		}
		else if (opname == "FPSqrt")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPSqrt: wE=" << wE << " wF=" << wF  << endl;
			op = new FPSqrt(target, wE, wF);
			addOperator(oplist, op);
		}
#ifdef HAVE_SOLLYA
		else if (opname == "FPSqrtPoly")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int degree = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPSqrtPoly: wE=" << wE << " wF=" << wF << " degree =" << degree << endl;
			op = new FPSqrtPoly(target, wE, wF, false, degree);
			addOperator(oplist, op);
		}
#endif // HAVE_SOLLYA

		else if(opname=="FPLargeAcc"){
			int nargs = 5;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wEX = checkStrictlyPositive(argv[i++], argv[0]);
				int wFX = checkStrictlyPositive(argv[i++], argv[0]);
				int MaxMSBX = atoi(argv[i++]); // may be negative
				int LSBA = atoi(argv[i++]); // may be negative
				int MSBA = atoi(argv[i++]); // may be negative
				cerr << "> Long accumulator , wEX="<<wEX<<", wFX="<<wFX<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<"\n";
				op = new FPLargeAcc(target, wEX, wFX, MaxMSBX, LSBA, MSBA);
				addOperator(oplist, op);
			}
		}

		// hidden and undocumented
		else if(opname=="FPLargeAccPrecTest"){
			int nargs = 6; // same as FPLargeAcc, plus an iteration count
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wEX = atoi(argv[i++]);
				int wFX = atoi(argv[i++]);
				int MaxMSBX = atoi(argv[i++]);
				int LSBA = atoi(argv[i++]);
				int MSBA = atoi(argv[i++]);
				int n = atoi(argv[i++]);
				cerr << "> Test of long accumulator accuracy, wEX="<<wEX<<", wFX="<<wFX<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<", "<< n << " tests\n";
				FPLargeAcc * op = new FPLargeAcc(target, wEX, wFX, MaxMSBX, LSBA, MSBA);
				// op->test_precision(n);
				op->test_precision2();
			}    
		}

		else if(opname=="LargeAccToFP"){
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int LSBA = atoi(argv[i++]); // may be negative
				int MSBA = atoi(argv[i++]); // may be negative
				int wE_out = checkStrictlyPositive(argv[i++], argv[0]);
				int wF_out = checkStrictlyPositive(argv[i++], argv[0]);
				cerr << "> Post-Normalization unit for Long accumulator, LSBA="<<LSBA<<", MSBA="<<MSBA<<" wE_out="<<wE_out<<", wF_out="<<wF_out<<"\n";
				op = new LargeAccToFP(target, LSBA, MSBA, wE_out, wF_out);
				addOperator(oplist, op);
			}
		}
#ifndef _WIN32
		// hidden and undocumented
		else if(opname=="DotProdPrecTest"){
			int nargs = 7; // same as FPLargeAcc, plus an iteration count
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wFX = checkStrictlyPositive(argv[i++], argv[0]);
				int wFY = checkStrictlyPositive(argv[i++], argv[0]);
				int MaxMSBX = atoi(argv[i++]);
				int LSBA = atoi(argv[i++]);
				int MSBA = atoi(argv[i++]);
				int n = atoi(argv[i++]);
				cerr << "> Test of FPDotProduct accuracy, wEX="<<wE<<", wFX="<<wFX<<", wFY="<<wFY<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<", "<< n << " tests\n";
				FPDotProduct * op = new FPDotProduct(target, wE, wFX, wFY, MaxMSBX, LSBA, MSBA);
				op->test_precision(n);
			}    
		}
#endif
		else if(opname=="FPDotProduct"){
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				int wE = checkStrictlyPositive(argv[i++], argv[0]);
				int wFX = checkStrictlyPositive(argv[i++], argv[0]);
				int wFY = checkStrictlyPositive(argv[i++], argv[0]);
				int MaxMSBX = atoi(argv[i++]); // may be negative
				int LSBA = atoi(argv[i++]); // may be negative
				int MSBA = atoi(argv[i++]); // may be negative
				double ratio = atof(argv[i++]); // may be negative
				cerr << "> FPDotProduct , wE="<<wE<<", wFX="<<wFX<<", wFY="<<wFY<<", MaxMSBX="<<MaxMSBX<<", LSBA="<<LSBA<<", MSBA="<<MSBA<<", ratio="<<ratio<<"\n";
				op = new FPDotProduct(target, wE, wFX, wFY, MaxMSBX, LSBA, MSBA, ratio);
				addOperator(oplist, op);
			}
		}
//		else if(opname=="PolynomialEvaluator"){
//			int nargs = 1;
//			if (i+nargs > argc)
//				usage(argv[0],opname);
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
//				addOperator(oplist, op);
//			}
//		}
				
#ifndef _WIN32
#ifdef HAVE_SOLLYA
		else if (opname == "HOTBM") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int wI = checkStrictlyPositive(argv[i++], argv[0]);
			int wO = checkStrictlyPositive(argv[i++], argv[0]);
			int n  = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> HOTBM func='" << func << "', wI=" << wI << ", wO=" << wO <<endl;
			op = new HOTBM(target, func, "", wI, wO, n);
			addOperator(oplist, op);
		}

		else if (opname == "HOTBMFX") {
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int wE_in = atoi(argv[i++]); // may be negative
			int wF_in = atoi(argv[i++]); // may be negative
			int wE_out = atoi(argv[i++]); // may be negative
			int wF_out = atoi(argv[i++]); // may be negative
			int n  = checkStrictlyPositive(argv[i++], argv[0]);
			
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
			addOperator(oplist, op);
		}
		
		else if (opname == "HOTBMRange") {
			int nargs = 7;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int wI = checkStrictlyPositive(argv[i++], argv[0]);
			int wO = checkStrictlyPositive(argv[i++], argv[0]);
			int n  = checkStrictlyPositive(argv[i++], argv[0]);
			double xmin = atof(argv[i++]);
			double xmax = atof(argv[i++]);

			// xmax < xmin is a valid use case...
			double scale = atof(argv[i++]);
			cerr << "> HOTBM func='" << func << "', wI=" << wI << ", wO=" << wO
			     << ", xmin=" << xmin << ", xmax=" << xmax << ", scale=" << scale <<endl;
			op = new HOTBM(target, func, "", wI, wO, n, xmin, xmax, scale);
			addOperator(oplist, op);
		}
		
#endif // HAVE_SOLLYA

#endif

		else if (opname == "FPExp")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPExp: wE=" << wE << " wF=" << wF << endl;
			op = new FPExp(target, wE, wF, 0, 0);
			addOperator(oplist, op);
		}


		else if (opname == "FPExpExpert")
		{
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int k=atoi(argv[i++]);
			int d=atoi(argv[i++]);
			int g=atoi(argv[i++]);
			int fullInput=checkBoolean(argv[i++],  argv[0]);
			cerr << "> FPExp: wE=" << wE << " wF=" << wF << endl;
			op = new FPExp(target, wE, wF, k, d, g, fullInput);
			addOperator(oplist, op);
		}

		else if (opname == "FPLog")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int inTableSize=atoi(argv[i++]);
			cerr << "> FPLog: wE=" << wE << " wF=" << wF << endl;
			op = new FPLog(target, wE, wF, inTableSize);
			addOperator(oplist, op);
		}

		else if (opname == "FPPowerExpert")
		{
			int nargs = 6;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit 

			//int logTableSize, int expTableSize, int expDegree, int expG, int logG int type
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int type=atoi(argv[i++]);
			int logTableSize=atoi(argv[i++]);
			int expTableSize=atoi(argv[i++]);
			int expDegree=atoi(argv[i++]);
			cerr << "> "<< (type==0?"FPPow":"FPPowr") << ": wE=" << wE << " wF=" << wF 
			     << "logTableSize" << logTableSize << "expTableSize" << expTableSize << "expDegree" << expDegree << endl;
			op = new FPPow(target, wE, wF, type, logTableSize, expTableSize, expDegree);
			addOperator(oplist, op);
		}

		else if (opname == "FPPow")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit 

			//int logTableSize, int expTableSize, int expDegree, int expG, int logG
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPPowr: wE=" << wE << " wF=" << wF << endl;
			op = new FPPow(target, wE, wF, 0);
			addOperator(oplist, op);
		}

		else if (opname == "FPPowr")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit 

			//int logTableSize, int expTableSize, int expDegree, int expG, int logG
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPPowr: wE=" << wE << " wF=" << wF << endl;
			op = new FPPow(target, wE, wF, 1);
			addOperator(oplist, op);
		}

		else if (opname == "InputIEEE")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wEI = checkStrictlyPositive(argv[i++], argv[0]);
			int wFI = checkStrictlyPositive(argv[i++], argv[0]);
			int wEO = checkStrictlyPositive(argv[i++], argv[0]);
			int wFO = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> InputIEEE: wEI=" << wEI << " wFI=" << wFI << " wEO=" << wEO << " wFO=" << wFO << endl;
			op = new InputIEEE(target, wEI, wFI, wEO, wFO);
			addOperator(oplist, op);
		}

		else if (opname == "OutputIEEE")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wEI = checkStrictlyPositive(argv[i++], argv[0]);
			int wFI = checkStrictlyPositive(argv[i++], argv[0]);
			int wEO = checkStrictlyPositive(argv[i++], argv[0]);
			int wFO = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> OutputIEEE: wEI=" << wEI << " wFI=" << wFI << " wEO=" << wEO << " wFO=" << wFO << endl;
			op = new OutputIEEE(target, wEI, wFI, wEO, wFO);
			addOperator(oplist, op);
		}

		else if (opname == "Collision")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int optimize = checkBoolean(argv[i++], argv[0]);
			cerr << "> Collision: wE=" << wE << " wF=" << wF << (optimize==0? " using FP operators" : " optimized version") << endl;
			op = new Collision(target, wE, wF, optimize);
			addOperator(oplist, op);
		}
		else if (opname == "FPSumOf3Squares")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int optimize = checkBoolean(argv[i++], argv[0]);
			cerr << "> FPSumOf3Squares: wE=" << wE << " wF=" << wF << (optimize==0? " using FP operators" : " optimized version") << endl;
			op = new FPSumOf3Squares(target, wE, wF, optimize);
			addOperator(oplist, op);
		}
		else if (opname == "IntSquarer")
		{
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wIn = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> IntSquarer: wIn=" << wIn << endl;
			op = new IntSquarer(target, wIn);
			addOperator(oplist, op);
		}
#ifndef _WIN32
#ifdef HAVE_LNS
		else if (opname == "LNSAddSub")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> LNSAddSub: wE=" << wE << " wF=" << wF << endl;
			op = new LNSAddSub(target, wE, wF);
			if(cl_name!="")	op->setName(cl_name);
			addOperator(oplist, op);
		}
		else if (opname == "LNSMul")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> LNSMul: wE=" << wE << " wF=" << wF << endl;
			op = new LNSMul(target, wE, wF);
			if(cl_name!="")	op->setName(cl_name);
			addOperator(oplist, op);
		}
		else if (opname == "LNSDiv")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> LNSDiv: wE=" << wE << " wF=" << wF << endl;
			op = new LNSDiv(target, wE, wF);
			addOperator(oplist, op);
		}
		else if (opname == "LNSSqrt")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> LNSSqrt: wE=" << wE << " wF=" << wF << endl;
			op = new LNSSqrt(target, wE, wF);
			if(cl_name!="")	op->setName(cl_name);
			addOperator(oplist, op);
		}
		else if (opname == "LNSAdd")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int o = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> LNSAdd: wE=" << wE << " wF=" << wF << " o=" << o << endl;
			op = new LNSAdd(target, wE, wF, o);
			addOperator(oplist, op);
		}
		// Undocumented LNS operators, for debugging purposes
#if 0
		else if (opname == "CotranF1")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int j = checkStrictlyPositive(argv[i++], argv[0]);
			int wE = atoi(argv[i++]);
			cerr << "> CotranF1: wF=" << wF << " j=" << j << " wE=" << wE << endl;
			op = new TableOp(target, new CotranF1Table(wF, j, wE));
			addOperator(oplist, op);
		}
		else if (opname == "CotranF2")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int j = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> CotranF2: wF=" << wF << " j=" << j << endl;
			op = new TableOp(target, new CotranF2Table(wF, j));
			addOperator(oplist, op);
		}
		else if (opname == "CotranF3")
		{
			int nargs = 2;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int j = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> CotranF3: wF=" << wF << " j=" << j << endl;
			op = new TableOp(target, new CotranF3Table(wF, j));
			addOperator(oplist, op);
		}
#endif // #if 0
		else if (opname == "Cotran")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int j = checkStrictlyPositive(argv[i++], argv[0]);
			int wECotran = atoi(argv[i++]);
			int o = atoi(argv[i++]);
			cerr << "> Cotran: wE=" << wE << " wF=" << wF << " j=" << j
				<< " wECotran=" << wECotran << " o=" << o << endl;
			op = new Cotran(target, wE, wF, j, wECotran);
			addOperator(oplist, op);
		}
		else if (opname == "CotranHybrid")
		{
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = atoi(argv[i++]);	// can be null or negative
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int j = checkStrictlyPositive(argv[i++], argv[0]);
			int wECotran = atoi(argv[i++]);
			int o = atoi(argv[i++]);
			cerr << "> Cotran: wE=" << wE << " wF=" << wF << " j=" << j 
				<< " wECotran=" << wECotran << " o=" << o << endl;
			op = new CotranHybrid(target, wE, wF, j, wECotran);
			addOperator(oplist, op);
		}
		else if (opname == "AtanPow")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			int o = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> AtanPow: wE=" << wE << " wF=" << wF << " o=" << o << endl;
			op = new AtanPow(target, wE, wF, o);
			addOperator(oplist, op);
		}
		else if (opname == "LogSinCos")
		{
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			int fL = checkStrictlyPositive(argv[i++], argv[0]);
			int fTheta = checkStrictlyPositive(argv[i++], argv[0]);
			int o = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> LogSinCos: fL=" << fL
				<< " fTheta=" << fTheta << " o=" << o << endl;
			op = new LogSinCos(target, fL, fTheta, o);
			addOperator(oplist, op);
		}
#endif

#endif //ifndef _WIN32
		else if (opname == "Wrapper") {
			int nargs = 0;
			if (i+nargs > argc)
				usage(argv[0],opname);
			else {
				if(oplist.empty()){
					cerr<<"ERROR: Wrapper has no operator to wrap (it should come after the operator it wraps)"<<endl;
					usage(argv[0],opname);
				}
				Operator* toWrap = oplist.back();
				cerr << "> Wrapper for " << toWrap->getName()<<endl;
				op =new Wrapper(target, toWrap);
				addOperator(oplist, op);
			}
		}
#ifdef HAVE_SOLLYA
#if 0
		else if (opname == "PolyTableGenerator") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int wO = atoi(argv[i++]);
			int n  = checkStrictlyPositive(argv[i++], argv[0]);
			
						
			cerr << "> PolyTableGenerator func='" << func << "', wO=" << wO <<endl;	
			
			Operator* tg = new PolyTableGenerator(target, func,  wO, n);
				addOperator(oplist, tg);
			
		}
#endif		

		else if (opname == "FunctionTable") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int wI = checkStrictlyPositive(argv[i++], argv[0]);
			int lsbO = atoi(argv[i++]);
			int msbO = atoi(argv[i++]);
			
				
			cerr << "> FunctionTable func='" << func << "', wI=" << wI << ", lsbO=" << lsbO << ", msbO=" << msbO << endl;	
			Operator* tg = new FunctionTable(target, func, wI, lsbO, msbO);
			addOperator(oplist, tg);
		}

		else if (opname == "FunctionEvaluator") {
			int nargs = 4;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string func = argv[i++];
			int wI = checkStrictlyPositive(argv[i++], argv[0]);
			int wO = atoi(argv[i++]);
			int n  = checkStrictlyPositive(argv[i++], argv[0]);
			
				
			cerr << "> FunctionEvaluator func='" << func << "', wI=" << wI << ", wO=" << wO << ", degree=" << n << endl;	
			string arg=func+",0,1,1"; // we are not sure it works for other values
			Operator* tg = new FunctionEvaluator(target, arg, wI, wO, n);
			addOperator(oplist, tg);
		}

		else if (opname == "FPPipeline") {
			int nargs = 3;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			string filename = argv[i++];
			int wE = checkStrictlyPositive(argv[i++], argv[0]);
			int wF = checkStrictlyPositive(argv[i++], argv[0]);
			cerr << "> FPPipeline filename='" << filename << "', wE=" << wE << ", wF=" << wF << endl;	
			Operator* tg = new FPPipeline(target, filename, wE, wF);
			addOperator(oplist, tg);
		}

#endif

		else if (opname == "TestBench") {
			int nargs = 1;
			if (i+nargs > argc)
				usage(argv[0],opname); // and exit
			if(oplist.empty()){
				cerr<<"ERROR: TestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0],opname); // and exit
			}
			int n = checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = oplist.back();
			Operator* op = new TestBench(target, toWrap, n);
			cerr << "> TestBench for " << toWrap->getName()<<endl;
			addOperator(oplist, op);
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
				usage(argv[0],opname); // and exit
			if(oplist.empty()){
				cerr<<"ERROR: TestBench has no operator to wrap (it should come after the operator it wraps)"<<endl;
				usage(argv[0],opname); // and exit
			}
			int n = atoi(argv[i++]);//checkPositiveOrNull(argv[i++], argv[0]);
			Operator* toWrap = oplist.back();
			Operator* op = new TestBench(target, toWrap, n, true);
			cerr << "> TestBench for " << toWrap->getName()<<endl;
			addOperator(oplist, op);
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
	

	target = new Virtex5();

	vector<Operator*> oplist;

	try {
		parseCommandLine(argc, argv, oplist);
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
	Operator::outputVHDLToFile(oplist, file);
	file.close();
	

	cerr << endl<<"Final report:"<<endl;
	for(i=0; i<oplist.size(); i++) {
		oplist[i]->outputFinalReport(0);
	}
	
	cerr<< "Output file: " << filename <<endl;
	return 0;
}



