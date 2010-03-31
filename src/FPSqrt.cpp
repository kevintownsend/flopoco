/*
  Floating Point Square Root for FloPoCo
 
  Authors : 
  Jeremie Detrey, Florent de Dinechin (digit-recurrence version)
  Mioara Joldes, Guillaume Revy, Bogdan Pasca (polynomial version)

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

 */
 


#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntAdder.hpp"
#include "IntMultiplier.hpp"
#include "IntSquarer.hpp"
#include "FPSqrt.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

#define DEBUGVHDL 0
	//#define LESS_DSPS

	FPSqrt::FPSqrt(Target* target, int wE, int wF, bool useDSP, bool correctlyRounded) :
		Operator(target), wE(wE), wF(wF), useDSP(useDSP), correctRounding(correctlyRounded) {

		ostringstream name;

		name<<"FPSqrt_"<<wE<<"_"<<wF;

		uniqueName_ = name.str(); 

		// -------- Parameter set up -----------------

		addFPInput ("X", wE, wF);
		addFPOutput("R", wE, wF);



		if(useDSP) {
			////////////////////////////////////////////////////////////////////////////////////
			//      Original polynomial version 

			/*These are the amount of shifts with respect to 0 that the coefficients a2 a1 and a0 are shifted */
			if(wF!=23)
				throw "Only wF=23 at the moment";
			int coeff_msb[3];
			coeff_msb[0] = 0; //a0
			coeff_msb[1] =-1; //a1
			coeff_msb[2] =-3; //a2
	
			int coeffStorageSizes[3];
			coeffStorageSizes[0] = 27; 
			coeffStorageSizes[1] = 17;
			coeffStorageSizes[2] = 9;
	
			int coeffTableWidth = coeffStorageSizes[0] + coeffStorageSizes[1] + coeffStorageSizes[2];
	
	
			int msb_x = -6; //number of zeros after the dot
	
			int tableAddressWidth = -msb_x + 2; //the +1 comes from the LSB of the exponent

			int sizeOfX = 16;
			//make the truncature before the first addition
			int keepBits;
	
			//how many extra bits I keep from the first addition result for correct rounding
			int keepBits2;
	
			if (correctRounding){
				keepBits = 17;
				keepBits2 = 2; //was 2

			}
			else{
				keepBits = 0;
				keepBits2 = 0;
			}
	
			vhdl << "--Split the Floating-Point input"<<endl;
			vhdl << tab << declare("fracX",wF) << " <= X" << range(wF-1, 0) << ";"  << endl; 
			vhdl << tab << declare("expX", wE)<< "  <= X" << range(wE+wF-1, wF) << ";"  << endl; 
			vhdl << "--A concatenation of the exception bits and the sign bit"<<endl;
			vhdl << tab << declare("excsX", 3) << " <= X" << range(wE+wF+2, wE+wF) << ";"  << endl; 
			vhdl << "--If the real exponent is odd"<<endl;
			vhdl << tab << declare("OddExp")   << " <= not(expX(0));"  << endl;  
	
			//first estimation of the exponent
			vhdl << tab << declare("expBiasPostDecrement", wE+1) << " <= " << "CONV_STD_LOGIC_VECTOR("<< (1<<(wE-1))-2 <<","<<wE+1<<")"<<";"<<endl;
			vhdl << tab << declare("expPostBiasAddition", wE+1) << " <= " << "( \"0\" & expX) + expBiasPostDecrement + not(OddExp);"<<endl;
	
			//the addres bits for the coefficient ROM
			vhdl << tab << declare("address", tableAddressWidth) << " <= OddExp & X" << range(wF-1, wF-tableAddressWidth+1) << ";"  << endl; //the MSB of address is the LSB of the exponent

			//get the correct size of x for the multiplication
			vhdl << tab << declare("lowX", sizeOfX + 1) << " <= " << "(\"0\" & X"<<range(sizeOfX-1,0) << ") when OddExp='0' else "
				  << "(X"<<range(sizeOfX-1,0) << " & \"0\");"<<endl<<endl;
			//instantiate the coefficient table
			Table* t;
			if (correctRounding)
				t = new PolynomialTableCorrectRounded(target, tableAddressWidth, coeffTableWidth);
			else
				t = new PolynomialTable(target, tableAddressWidth, coeffTableWidth);
			oplist.push_back(t);
			combinatorialOperator = false;
	
			nextCycle();//// this pipeline level is needed in order to infer a BRAM here
	
			inPortMap    (t, "X", "address");
			outPortMap   (t, "Y", "data");
			vhdl << instance(t, "SQRT_Coeffs_Table");

			syncCycleFromSignal("data");

			nextCycle();/////////////////////////////////// The Coefficent ROM has a registered output
                                      
			//get a2 from memory
			vhdl <<endl << tab << declare("a2", coeffStorageSizes[2]) << "<= data"<<range(coeffTableWidth-1, coeffTableWidth-coeffStorageSizes[2]) <<";"<<endl;
			//perform (-a2)*x, each term is <= than 17 bits so * will be mapped into one DSP block
			vhdl << tab << declare("prodA2X",coeffStorageSizes[2] + sizeOfX + 1) << " <= lowX  * a2 ;" << endl;

			nextCycle();/////////////////////////////////// Will be absorbed by the DSP macro
	
			//get a1 from memory
			vhdl <<endl << tab << declare("a1",coeffStorageSizes[1]) << " <= data"<<range(coeffStorageSizes[0] + coeffStorageSizes[1] - 1, coeffStorageSizes[0]) << ";" <<endl;
			//sign-extend and pad a1 for a1 - (-a2*x)
			vhdl << tab << declare("signExtA1ZeroPad", 1 + coeffStorageSizes[1] + keepBits) << " <= " << "\"0\" & a1" << " & " << zg(keepBits, 0) << ";" << endl;
			//sign-extend and align the -a2*x product
			vhdl << tab << declare("signExtAlignedProdA2X", 1 + coeffStorageSizes[1] + keepBits) << " <= " << "\"0\" & " << zg(coeff_msb[1]-(coeff_msb[2]+msb_x),0) << " & " 
				  << "prodA2X"<< range(coeffStorageSizes[2] + sizeOfX, coeffStorageSizes[2] + sizeOfX - (1 + coeffStorageSizes[1] + keepBits + coeff_msb[2]+msb_x )+1) << ";"<<endl;
			vhdl << tab << declare("negSignExtAlignedProdA2X", 1 + coeffStorageSizes[1] + keepBits) << " <= not(signExtAlignedProdA2X);"<<endl;
			//subtract (-a2*x) from a1, => instantiate IntAdder
			IntAdder* add1 = new IntAdder(target, 1 + coeffStorageSizes[1] + keepBits);
			oplist.push_back(add1);
	
			inPortMap(add1,"X", "signExtA1ZeroPad");
			inPortMap(add1,"Y", "negSignExtAlignedProdA2X");
			inPortMapCst(add1,"Cin", "'1'");
			outPortMap(add1,"R", "add1Res");
			vhdl << instance(add1, "Adder_a1_prod_x_a2");
	
			syncCycleFromSignal("add1Res"); 
			nextCycle();//////////////////
	
			//perform the multiplication between x and ( a1 + a2x )  //TODO pipeline if keepbits2 > 0  
			//FIXME for now we instantiate an int Multiplier, but we can do better

#ifdef LESS_DSPS		
			if ((correctRounding)&&(false)){
#else
				if (correctRounding){
#endif
					IntMultiplier * my_mul = new IntMultiplier(target, (sizeOfX+1), 34);
					oplist.push_back(my_mul);
		
					vhdl << tab << declare("justASignal",34) << " <= " << zg(15,0) << " & add1Res"<<range(coeffStorageSizes[1] + keepBits-1, keepBits - keepBits2) << ";" << endl;
	
					inPortMap(my_mul,"X", "lowX");
					//		inPortMapCst(my_mul,"Y", use("add1Res")+range(coeffStorageSizes[1] + keepBits-1, keepBits - keepBits2) );
					inPortMap(my_mul,"Y", "justASignal");
					//		outPortMap(my_mul, "R", "prodXA1sumA2X");
					outPortMap(my_mul, "R", "prodXA1sumA2X_large");
					vhdl << tab << instance(my_mul,"SecondMultiplier");
	
					syncCycleFromSignal("prodXA1sumA2X_large"); 
		
					vhdl << tab << declare("prodXA1sumA2X",36) << " <= prodXA1sumA2X_large"<<range(35,0) << ";" << endl;
				}else{
					IntMultiplier * my_mul = new IntMultiplier(target, (sizeOfX+1), coeffStorageSizes[1] + keepBits2);
					oplist.push_back(my_mul);
					inPortMap(my_mul,"X", "lowX");
					inPortMapCst(my_mul,"Y", "add1Res"+range(coeffStorageSizes[1] + keepBits-1, keepBits - keepBits2) );
					outPortMap(my_mul, "R", "prodXA1sumA2X");
					vhdl << tab << instance(my_mul,"SecondMultiplier");
	
					syncCycleFromSignal("prodXA1sumA2X"); 
				}
				//vhdl << tab << declare("prodXA1sumA2X", (sizeOfX+1)+ coeffStorageSizes[1] + keepBits2 ) << " <= " << use("lowX") << " * " 
				//                                                                                    << use("add1Res")<<range( coeffStorageSizes[1] + keepBits-1, keepBits - keepBits2)<<";" <<endl;
#ifndef LESS_DSPS    
				nextCycle();/////////////////                                                                               
#endif
				//compose the operands for the addition a0 + [ prev_computation ]
				//fetch a0 from memory
				vhdl << endl << tab << declare("a0",coeffStorageSizes[0]) << " <= data" << range(coeffStorageSizes[0]-1,0) << ";" << endl;
				vhdl << tab << declare ("ovfGuardA0", 1 + coeffStorageSizes[0]) << " <= " << " \"0\" & a0"<< ";" << endl;
				vhdl << tab << declare ("ovfGuardAlignProdXA1sumA2X", 1 + coeffStorageSizes[0]) << " <= " << " \"0\" & " << zg((1+coeff_msb[0])+1-msb_x , 0)  << " & "
					  << "prodXA1sumA2X"<<range((sizeOfX+1)+ coeffStorageSizes[1] -1 + keepBits2, keepBits2 + (sizeOfX+1)+ coeffStorageSizes[1] - (coeffStorageSizes[0]- (-msb_x+1+(1+ coeff_msb[0])))) << ";" <<endl; 
	
				IntAdder * add2 =  new IntAdder(target, 1 + coeffStorageSizes[0]);
				oplist.push_back(add2);
				inPortMap(add2,"X", "ovfGuardA0");
				inPortMap(add2,"Y", "ovfGuardAlignProdXA1sumA2X");
				inPortMapCst(add2,"Cin", "'0'");
				outPortMap(add2,"R", "sumA0ProdXA1sumA2X");
				vhdl << instance(add2, "FinalAdder");
	
				syncCycleFromSignal("sumA0ProdXA1sumA2X"); 
		
				if (!correctlyRounded){
	
					vhdl << tab << declare("finalFrac", wF) << " <= sumA0ProdXA1sumA2X" << range(coeffStorageSizes[0]-2, coeffStorageSizes[0]-wF-1) << ";" << endl;
					vhdl << tab << declare("finalExp", wE) << " <= expPostBiasAddition" << range(wE,1) <<";"<<endl;

					vhdl << tab << "-- sign/exception handling" << endl;
					vhdl << tab << "with excsX select" <<endl
						  << tab << tab <<  declare("exnR", 2) << " <= " << "\"01\" when \"010\", -- positive, normal number" << endl
						  << tab << tab << "excsX" << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
						  << tab << tab << "\"11\" when others;"  << endl;

					vhdl << tab << "R <= exnR & excsX(0) & finalExp & finalFrac;" << endl; 
				}else{
					//		vhdl << tab << declare("normalizeBit",1) << " <= " << use("sumA0ProdXA1sumA2X") << "(" << coeffStorageSizes[0] << ")"<<";"<<endl;
					//		nextCycle();/////////////////////////
					//		vhdl << tab << declare("preSquareFrac", wF+2) << " <= " << use("sumA0ProdXA1sumA2X") << range(coeffStorageSizes[0], coeffStorageSizes[0]-wF-1) << " when " << use("normalizeBit") <<"='1' else "
					//			                                           << use("sumA0ProdXA1sumA2X") << range(coeffStorageSizes[0]-1, coeffStorageSizes[0]-wF-2) << ";" << endl;
					//		vhdl << tab << declare("preSquareExp", wE) << " <= " << use("expPostBiasAddition") << range(wE,1) << " + " << use("normalizeBit")<<";"<<endl;

					vhdl << tab << declare("preSquareFrac", wF+2) << " <= sumA0ProdXA1sumA2X" << range(coeffStorageSizes[0]-1, coeffStorageSizes[0]-wF-2) << ";" << endl;
					vhdl << tab << declare("preSquareExp", wE) << " <= expPostBiasAddition" << range(wE,1) <<";"<<endl;


					vhdl << tab << declare("preSquareConcat", 1 + wE + wF+1) << " <= " << zg(1,0) << " & preSquareExp & preSquareFrac"<<range(wF,0)<<";"<<endl;;

#ifndef LESS_DSPS
					nextCycle();//////////////////	    
#endif
					IntAdder *predictAdder = new IntAdder(target, 1 + wE + wF+1);
					oplist.push_back(predictAdder);

					inPortMap(predictAdder,"X","preSquareConcat");	
					inPortMapCst(predictAdder,"Y",zg(1 + wE + wF+1,0));
					inPortMapCst(predictAdder,"Cin","'1'");
					outPortMap(predictAdder,"R","my_predictor");
					vhdl << tab << instance(predictAdder, "Predictor");
	    	    
					IntSquarer *iSquarer = new IntSquarer(target,max(wF+2,34));
					oplist.push_back(iSquarer);
	    
					vhdl << tab << declare("op1",max(wF+2,34)) << " <= " << zg(34-(wF+2),0) << " & preSquareFrac;" << endl;
	    
					inPortMap(iSquarer, "X", "op1");
					outPortMap(iSquarer,"R", "sqrResult");
					vhdl << instance(iSquarer,"FractionSquarer");
	    
					syncCycleFromSignal("sqrResult");
	    
					vhdl << tab << declare("approxSqrtXSqr", 2*(wF+2) + 1) << " <= " << zg(1,0) << " & sqrResult"<<range(2*(wF+2)-1,0)<<";"<<endl;
					vhdl << tab << declare("realXFrac", 2*(wF+2) + 1) << " <= " << "( \"001\" & fracX & " << zg(2*(wF+2) + 1-2-wF-1 ,0)<<") when OddExp='0' else "<<endl;
					vhdl << tab << tab << "( \"01\" & fracX & " << zg(2*(wF+2) + 1-2-wF,0)<<");"<<endl;
	    
					vhdl << tab << declare("negRealXFrac", 2*(wF+2) + 1) << " <= " << "not(realXFrac);"<<endl;
	    
					IntAdder *myIntAdd = new IntAdder(target, 2*(wF+2) + 1);
					oplist.push_back(myIntAdd);
	    
					inPortMap(myIntAdd,"X","approxSqrtXSqr");	
					inPortMap(myIntAdd,"Y","negRealXFrac");
					inPortMapCst(myIntAdd,"Cin","'1'");
					outPortMap(myIntAdd,"R","my_add_result");
					vhdl << tab << instance(myIntAdd, "Comparator");

					syncCycleFromSignal("my_add_result");
					vhdl << tab << declare("greater",1) << " <= my_add_result"<<of(2*(wF+2))<<";"<<endl;

					syncCycleFromSignal("my_predictor");    

					vhdl << tab << "-- sign/exception handling" << endl;
					vhdl << tab << "with excsX select" <<endl
						  << tab << tab <<  declare("exnR", 2) << " <= " << "\"01\" when \"010\", -- positive, normal number" << endl
						  << tab << tab << "excsX" << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
						  << tab << tab << "\"11\" when others;"  << endl;

					vhdl << tab << "R <= (exnR & excsX(0) & preSquareConcat"<<range(wE + wF,1) << ") when greater='0' else "<<endl;;
					vhdl << tab << tab << "(exnR & excsX(0) & my_predictor"<<range(wE + wF,1) << ");"<<endl;
	    
				}
			}
			////////////////////////////////////////////////////////////////////////////////////
			else {
				// Digit-recurrence implementation recycled from FPLibrary
				//cout << "   DDDD" <<  target->adderDelay(10) << "  " <<  target->localWireDelay() << "  " << target->lutDelay();
				vhdl << tab << declare("fracX", wF) << " <= X" << range(wF-1, 0) << "; -- fraction"  << endl; 
				vhdl << tab << declare("eRn0", wE) << " <= \"0\" & X" << range(wE+wF-1, wF+1) << "; -- exponent" << endl;
				vhdl << tab << declare("xsX", 3) << " <= X"<< range(wE+wF+2, wE+wF) << "; -- exception and sign" << endl;

				vhdl << tab << declare("eRn1", wE) << " <= eRn0 + (\"00\" & " << rangeAssign(wE-3, 0, "'1'") << ") + X(" << wF << ");" << endl;

				vhdl << tab << declare(join("w",wF+3), wF+4) << " <= \"111\" & fracX & \"0\" when X(" << wF << ") = '0' else" << endl
					  << tab << "       \"1101\" & fracX;" << endl;
				//		vhdl << tab << declare(join("d",wF+3)) << " <= '0';" << endl;
				//		vhdl << tab << declare(join("s",wF+3),1) << " <= '1';" << endl;

				double delay= target->lutDelay() + target->localWireDelay() + target->ffDelay(); // estimated delay so far (one mux)
				for(int step=1; step<=wF+2; step++) {
					int i = wF+3-step; // to have the same indices as FPLibrary 
					vhdl << tab << "-- Step " << i << endl;
					string di = join("d", i);
					string xi = join("x", i);
					string wi = join("w", i);
					string wip = join("w", i+1);
					string si = join("s", i);
					string sip = join("s", i+1);
					//			string zs = join("zs", i);
					string ds = join("ds", i);
					string xh = join("xh", i);
					string wh = join("wh", i);
					vhdl << tab << declare(di) << " <= "<< use(wip) << "("<< wF+3<<");" << endl;
					vhdl << tab << declare(xi,wF+5) << " <= "<< use(wip) << " & \"0\";" << endl;
					//			vhdl << tab << declare(zs,step,true) << " <= \"0\" & " << use(sip) << ";" << endl;
					vhdl << tab << declare(ds,step+3) << " <=  \"0\" & ";
					if (step>1)
						vhdl 	<< use(sip) << " & ";
					vhdl << " (not " << di << ") & " << di << " & \"1\";" << endl;
					vhdl << tab << declare(xh,step+3) << " <= " << xi << range(wF+4, wF+2-step) << ";" << endl;
					vhdl << tab << "with " << di << " select" << endl
						  << tab << tab <<  declare(wh, step+3) << " <= " << xh << " - " << ds << " when '0'," << endl
						  << tab << tab << "      " << xh << " + " << ds << " when others;" << endl;
					vhdl << tab << declare(wi, wF+4) << " <= " << wh << range(step+1,0);
					if(step <= wF+1) 
						vhdl << " & " << xi << range(wF+1-step, 0) << ";" << endl;  
					else
						vhdl << ";" << endl; 
					vhdl << tab << declare(si, step) << " <= ";
					if(step==1)
						vhdl << "not " << di << " ;"<< endl; 
					else
						vhdl << use(sip) /*<< range(step-1,1)*/ << " & not " << di << ";"<< endl; 
				
					// Pipeline management
					double stageDelay= target->adderDelay(step) + target->localWireDelay() + 2*target->lutDelay();
					delay += stageDelay;
					if (verbose>=2) {
						cout << "estimated delay for stage "<< step << " is " << stageDelay << "s" << endl;
						cout << "   cumulated delay would be " << delay << "s,   target is " << 1/target->frequency()<< endl;
					}
					if(delay > 1/target->frequency()) {
						// insert a pipeline register and reset the cumulated delay
						nextCycle();
						delay= target->ffDelay() + stageDelay;
						if (verbose>=2) 
							cout << "----inserted a register level" << endl;
					}
				}
				vhdl << tab << declare("d0") << " <= " << use("w1") << "(" << wF+3 << ") ;" << endl;
				vhdl << tab << declare("fR", wF+4) << " <= " << use("s1") << " & not d0 & '1';" << endl;

				// end of component FPSqrt_Sqrt in fplibrary
				vhdl << tab << "-- normalisation of the result, removing leading 1" << endl;
				vhdl << tab <<  "with fR(" << wF+3 << ") select" << endl
					  << tab << tab << declare("fRn1", wF+2) << " <= fR" << range(wF+2, 2) << " & (fR(1) or fR(0)) when '1'," << endl
					  << tab << tab << "        fR" <<range(wF+1, 0) << "                    when others;" << endl;
				vhdl << tab << declare("round") << " <= fRn1(1) and (fRn1(2) or fRn1(0)) ; -- round  and (lsb or sticky) : that's RN, tie to even" << endl;

				nextCycle();
		
				vhdl << tab << declare("fRn2", wF) << " <= " <<  use("fRn1") << range(wF+1, 2) <<" + (" << rangeAssign(wF-1, 1, "'0'") << " & " << use("round")  << "); -- rounding sqrt never changes exponents " << endl;
				vhdl << tab << declare("Rn2", wE+wF) << " <= " << use("eRn1") << " & " << use("fRn2") << ";" << endl;
		
				vhdl << tab << "-- sign and exception processing" << endl;
				vhdl << tab <<  "with " << use("xsX") << " select" << endl
					  << tab << tab << declare("xsR", 3) << " <= \"010\"  when \"010\",  -- normal case" << endl
					  << tab << tab <<  "       \"100\"  when \"100\",  -- +infty" << endl
					  << tab << tab <<  "       \"000\"  when \"000\",  -- +0" << endl
					  << tab << tab <<  "       \"001\"  when \"001\",  -- the infamous sqrt(-0)=-0" << endl
					  << tab << tab <<  "       \"110\"  when others; -- return NaN" << endl;
			
				vhdl << tab << "R <= xsR & Rn2; " << endl; 
			}


		}

		FPSqrt::~FPSqrt() {
		}






		void FPSqrt::emulate(TestCase * tc)
		{
			/* Get I/O values */
			mpz_class svX = tc->getInputValue("X");

			/* Compute correct value */
			FPNumber fpx(wE, wF);
			fpx = svX;
			mpfr_t x, r;
			mpfr_init2(x, 1+wF);
			mpfr_init2(r, 1+wF); 
			fpx.getMPFR(x);

			if(correctRounding) {
				mpfr_sqrt(r, x, GMP_RNDN);
				FPNumber  fpr(wE, wF, r);
				/* Set outputs */
				mpz_class svr= fpr.getSignalValue();
				tc->addExpectedOutput("R", svr);
			}
			else { // faithful rounding 
				mpfr_sqrt(r, x, GMP_RNDU);
				FPNumber  fpru(wE, wF, r);
				mpz_class svru = fpru.getSignalValue();
				tc->addExpectedOutput("R", svru);

				mpfr_sqrt(r, x, GMP_RNDD);
				FPNumber  fprd(wE, wF, r);
				mpz_class svrd = fprd.getSignalValue();
				/* Set outputs */
				tc->addExpectedOutput("R", svrd);
			}

			mpfr_clears(x, r, NULL);
		}




		// One test out of 4 fully random (tests NaNs etc)
		// All the remaining ones test positive numbers.

		void FPSqrt::buildRandomTestCases(TestCaseList* tcl, int n){

			TestCase *tc;
			mpz_class a;

			for (int i = 0; i < n; i++) {
				tc = new TestCase(this); 
				/* Fill inputs */
				if ((i & 3) == 0)
					a = getLargeRandom(wE+wF+3);
				else
					a  = getLargeRandom(wE+wF) + (mpz_class(1)<<(wE+wF+1)); // 010xxxxxx
				tc->addInput("X", a);

				/* Get correct outputs */
				emulate(tc);

				tcl->add(tc);
			}
		}
	}
