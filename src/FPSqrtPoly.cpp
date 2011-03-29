/*
  Floating Point Square Root, polynomial version, for FloPoCo
 
  Authors : Mioara Joldes, Bogdan Pasca (2010)

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

#ifdef HAVE_SOLLYA
#include "IntAdder.hpp"
#include "IntMultiplier.hpp"
#include "IntSquarer.hpp"
#include "FPSqrtPoly.hpp"

#define KEEP_HANDCRAFTED_VERSION 1

using namespace std;

namespace flopoco{
	extern vector<Operator*> oplist;

	FPSqrtPoly::FPSqrtPoly(Target* target, int wE, int wF, bool correctlyRounded, int degree):
		Operator(target), wE(wE), wF(wF), correctRounding(correctlyRounded) {

		ostringstream name;
		name<<"FPSqrtPoly_"<<wE<<"_"<<wF;
		uniqueName_ = name.str(); 

		// -------- Parameter set up -----------------
		addFPInput ("X", wE, wF);
		addFPOutput("R", wE, wF);



#if KEEP_HANDCRAFTED_VERSION
		if(wF==23) {
			////////////////////////////////////////////////////////////////////////////////////
			//      Original hand-crafted polynomial version, beats the automatically generated one 

			/*These are the amount of shifts with respect to 0 that the coefficients a2 a1 and a0 are shifted */
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
			vhdl << tab << declare("expBiasPostDecrement", wE+1) << " <= CONV_STD_LOGIC_VECTOR("<< (1<<(wE-1))-2 <<","<<wE+1<<");"<<endl;
			vhdl << tab << declare("expPostBiasAddition", wE+1) << " <= ( \"0\" & expX) + expBiasPostDecrement + not(OddExp);"<<endl;
	
			//the addres bits for the coefficient ROM
			vhdl << tab << declare("address", tableAddressWidth) << " <= OddExp & X" << range(wF-1, wF-tableAddressWidth+1) << ";"  << endl; //the MSB of address is the LSB of the exponent

			//get the correct size of x for the multiplication
			vhdl << tab << declare("lowX", sizeOfX + 1) << " <= (\"0\" & X"<<range(sizeOfX-1,0) << ") when OddExp='0' else "
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
			vhdl << tab << declare("signExtA1ZeroPad", 1 + coeffStorageSizes[1] + keepBits) << " <= \"0\" & a1 & " << zg(keepBits, 0) << ";" << endl;
			//sign-extend and align the -a2*x product
			vhdl << tab << declare("signExtAlignedProdA2X", 1 + coeffStorageSizes[1] + keepBits) << " <= \"0\" & " << zg(coeff_msb[1]-(coeff_msb[2]+msb_x),0) << " & " 
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
				//vhdl << tab << declare("prodXA1sumA2X", (sizeOfX+1)+ coeffStorageSizes[1] + keepBits2 ) << " <= lowX * " 
				// << "add1Res" << range( coeffStorageSizes[1] + keepBits-1, keepBits - keepBits2)<<";" <<endl;
#ifndef LESS_DSPS    
				nextCycle();/////////////////                                                                               
#endif
				//compose the operands for the addition a0 + [ prev_computation ]
				//fetch a0 from memory
				vhdl << endl << tab << declare("a0",coeffStorageSizes[0]) << " <= data" << range(coeffStorageSizes[0]-1,0) << ";" << endl;
				vhdl << tab << declare ("ovfGuardA0", 1 + coeffStorageSizes[0]) << " <=  \"0\" & a0;" << endl;
				vhdl << tab << declare ("ovfGuardAlignProdXA1sumA2X", 1 + coeffStorageSizes[0]) << " <=  \"0\" & " << zg((1+coeff_msb[0])+1-msb_x , 0)  << " & "
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
						  << tab << tab <<  declare("exnR", 2) << " <= \"01\" when \"010\", -- positive, normal number" << endl
						  << tab << tab << "excsX" << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
						  << tab << tab << "\"11\" when others;"  << endl;

					vhdl << tab << "R <= exnR & excsX(0) & finalExp & finalFrac;" << endl; 
				}

		}else{
#endif // KEEP_HANDCRAFTER_VERSION



		vhdl << tab << declare("excsX",3) << " <= X"<<range(wE+wF+2,wE+wF)<<";"<<endl;
		vhdl << tab << declare("sX",1) << "  <= X"<<of(wE+wF)<<";"<<endl;
		vhdl << tab << declare("expX",wE) << " <= X"<<range(wE+wF-1,wF)<<";"<<endl;
		vhdl << tab << declare("fX",wF+1) << " <= \"1\" & X"<<range(wF-1,0 )<<";"<<endl;
		
		vhdl << "--If the real exponent is odd"<<endl;
		vhdl << tab << declare("OddExp")    << " <= not(expX(0));"  << endl;  

		//first estimation of the exponent
		vhdl << tab << declare("expBiasPostDecrement", wE+1) << " <= CONV_STD_LOGIC_VECTOR("<< (1<<(wE-1))-2 <<","<<wE+1<<");"<<endl;
		vhdl << tab << declare("expPostBiasAddition", wE+1) << " <= ( \"0\" & expX) + expBiasPostDecrement + not(OddExp);"<<endl;
	
		vhdl << tab << "-- sign/exception handling" << endl;
		vhdl << tab << "with excsX select" <<endl
			  << tab << tab <<  declare("exnR", 2) << " <= \"01\" when \"010\", -- positive, normal number" << endl
			  << tab << tab << "excsX" << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
			  << tab << tab << "\"11\" when others;"  << endl;

//		FunctionEvaluator *fixpsqrt = new FunctionEvaluator(target, "sqrt(2+2*x),0,1,1;sqrt(1+x),0,1,1", wF+1, wF, degree); //TODO
//		oplist.push_back(fixpsqrt);

//***************************************************************************
		PiecewiseFunction *pf= new  PiecewiseFunction("sqrt(2+2*x),0,1,1;sqrt(1+x),0,1,1");
		PolyTableGenerator *tg = new PolyTableGenerator(target, pf, wF+1, degree);
		oplist.push_back(tg);
		combinatorialOperator = false;
		
		int k1, k2;
		k1 = (tg->getNrIntArray())[0];
		k2 = (tg->getNrIntArray())[1];
		
		int aa = tg->wIn;
	
		int dk1k2 = k1-k2;
		
		vhdl << tab << declare("theAddr",aa) << " <= (expX(0) & X"<<range(wF-1,wF-k1)<<") when expX(0)='0' else "<<endl
		                                 << "(expX(0) & " << zg(dk1k2,0) << " & X"<<range(wF-1,wF-k2)<<");"<<endl;
		vhdl << tab << declare("theFrac",wF-k2) << " <= ("<< zg(dk1k2,0) <<" & X"<<range(wF-k1-1,0)<<") when expX(0)='0' else "<<endl
		                                 << "X"<<range(wF-k2-1,0)<<";"<<endl;
		
		REPORT(DEBUG, "k1="<<k1<<" k2="<<k2);
		
		YVar* y = new YVar(wF-k2, -k2);
		
		PolynomialEvaluator *pe = new PolynomialEvaluator(target, tg->getCoeffParamVector(), y, wF+1, tg->getMaxApproxError() );
		oplist.push_back(pe);
		
//		wR = pe->getRWidth();
//		weightR = pe->getRWeight();

//		vhdl << tab << declare("addr", tg->wIn) << " <= Xaddr;"<<endl;
		nextCycle();/////////////////////////////////// The Coefficent ROM has a registered iunput
		
		inPortMap ( tg, "X", "theAddr");
		outPortMap ( tg, "Y", "Coef");
		vhdl << instance ( tg, "GeneratedTable" );
		
		syncCycleFromSignal("Coef");
		nextCycle();/////////////////////////////////// The Coefficent ROM has a registered output
		
		vhdl << tab << declare ("y",y->getSize()) << " <= theFrac;" << endl;
		
		/* get the coefficients */
		int lsb = 0, sizeS = 0;
		for (uint32_t i=0; i< pe->getCoeffParamVector().size(); i++){
			lsb += sizeS;
			sizeS = pe->getCoeffParamVector()[i]->getSize()+1;
			vhdl << tab << declare(join("a",i), sizeS ) << "<= Coef"<< range (lsb+sizeS-1, lsb) << ";" << endl;
		}
		
		inPortMap( pe, "Y", "y");
		for (uint32_t i=0; i< pe->getCoeffParamVector().size(); i++){
			inPortMap(pe,  join("a",i), join("a",i));
		}
		outPortMap( pe, "R", "rfx");
		vhdl << instance( pe, "PolynomialEvaluator");
		
		syncCycleFromSignal("rfx");


//***************************************************************************


//		inPortMap(fixpsqrt, "Xaddr", "theAddr");
//		inPortMap(fixpsqrt, "Xfrac", "theFrac");
//		outPortMap(fixpsqrt, "R", "rfx");
//		vhdl << instance(fixpsqrt, "FixPointSQRT");
//		
//		syncCycleFromSignal("rfx");
		
		
//		vhdl << tab << declare("sticky",1) << " <=  rfx"<<of(pe->getRWidth()-(pe->getRWeight()+wF)-1)<<";"<<endl;
		
		
		vhdl << tab << declare("sticky",1) << " <= '0' when rfx"<<range(pe->getRWidth()-(pe->getRWeight()+wF)-3,0) <<" = " 
		                                                        <<   zg(pe->getRWidth()-(pe->getRWeight()+wF)-2,0) << " else '1';"<<endl;
		
		vhdl << tab << declare("extentedf", 1 + wF + 2) << " <= rfx"<<range(pe->getRWidth()-pe->getRWeight(), pe->getRWidth()-(pe->getRWeight()+wF)-2)<<";"<<endl; 
//		                                                << " & sticky;"<<endl;
		                  
		nextCycle();                              
		IntAdder *a = new IntAdder(target, 1 + wF + 2);
		oplist.push_back(a);
		
		inPortMap(a, "X", "extentedf");
		inPortMapCst(a, "Y", zg(1 + wF + 2,0) );
		inPortMapCst(a, "Cin", "sticky");
		outPortMap(a, "R", "fPostRound");
		vhdl << instance(a, "Rounding_Adder");

		syncCycleFromSignal("fPostRound");
		

		addOutput("RFull", pe->getRWidth());
		vhdl << tab << " RFull <= rfx;" << endl; 
		
		vhdl << tab << " R <= exnR & sX & expPostBiasAddition"<<range(wE,1)<<" & fPostRound"<<range(wF, 1)<<";"<<endl;
		
		
//		cout << "The result number of bits is " << 	pe->getRWidth();
//		cout << "The result weight is " << 	        pe->getRWeight();
			
#if KEEP_HANDCRAFTED_VERSION
		}
#endif			
		
	}



	FPSqrtPoly::~FPSqrtPoly() {
	}

	void FPSqrtPoly::emulate(TestCase * tc)
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

	void FPSqrtPoly::buildRandomTestCases(TestCaseList* tcl, int n){

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
	TestCase* FPSqrtPoly::buildRandomTestCases(int i){

		TestCase *tc;
		mpz_class a;

		tc = new TestCase(this); 
		/* Fill inputs */
		if ((i & 3) == 0)
			a = getLargeRandom(wE+wF+3);
		else
			a  = getLargeRandom(wE+wF) + (mpz_class(1)<<(wE+wF+1)); // 010xxxxxx
		tc->addInput("X", a);

		/* Get correct outputs */
		emulate(tc);

		return tc;
	}
}

#endif //HAVE_SOLLYA

