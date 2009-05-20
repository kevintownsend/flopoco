/*
 * Floating Point Square Root for FloPoCo
 *
 * Authors : 
 * Jeremie Detrey, Florent de Dinechin (digit-recurrence version)
 * Mioara Joldes, Bogdan Pasca (polynomial version)
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
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntAdder.hpp"
#include "IntSquarer.hpp"
#include "FPSqrt.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPSqrt::FPSqrt(Target* target, int wE, int wF, bool useDSP, bool correctlyRounded) :
	Operator(target), wE(wE), wF(wF), useDSP(useDSP), correctRounding(correctlyRounded) {

	ostringstream name;

	name<<"FPSqrt_"<<wE<<"_"<<wF;

	uniqueName_ = name.str(); 

	// -------- Parameter set up -----------------

	addFPInput ("X", wE, wF);
	addFPOutput("R", wE, wF);

	if(useDSP) {
	
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
	int keepBits = 0;
	
	vhdl << tab << declare("fX",   wF) << " <= X" << range(wF-1, 0) << ";"  << endl; 
	vhdl << tab << declare("expX", wE) << " <= X" << range(wE+wF-1, wF) << ";"  << endl; 
	vhdl << tab << declare("exnsX", 3) << " <= X" << range(wE+wF+2, wE+wF) << ";"  << endl; 
	vhdl << tab << declare("OddExp")   << " <= not(expX(0));"  << endl;  
	
	//first estimation of the exponent
	vhdl << tab << declare("biasm1", wE+1) << " <= " << "CONV_STD_LOGIC_VECTOR("<< (1<<(wE-1))-2 <<","<<wE+1<<")"<<";"<<endl;
	vhdl << tab << declare("exp_addition", wE+1) << " <= " << "( \"0\" & " << use("expX") << ") + "<<use("biasm1") << " + not(" << use("OddExp") << ");"<<endl;
	
	vhdl << tab << declare("address", tableAddressWidth) << " <= "<< use("OddExp") << " & X" << range(wF-1, wF-tableAddressWidth+1) << ";"  << endl; //the MSB of address is the LSB of the exponent

	//get the correct size of x for the multiplicat	ion
	vhdl << tab << declare("lowX", sizeOfX + 1) << " <= " << "(\"0\" & " << use("X")<<range(sizeOfX-1,0) << ") when " << use("OddExp")<<"='0' else "
	                                            << "(" << use("X")<<range(sizeOfX-1,0) << " & \"0\") ;"<<endl;
	//instantiate the coefficient table
	PolynomialTable* t = new PolynomialTable(target, tableAddressWidth, coeffTableWidth);
	oplist.push_back(t);
	
	nextCycle();////
	
	inPortMapCst (t, "X", use("address"));
	outPortMap(t, "Y", "data");
	vhdl << instance(t, "SQRT_Coeffs_Table");

	syncCycleFromSignal("data");

	nextCycle();///////////////////////////////////

	                                      
	//get a2 from memory
	vhdl << tab << declare("a2", coeffStorageSizes[2]) << "<=" << use("data")<<range(coeffTableWidth-1, coeffTableWidth-coeffStorageSizes[2]) <<";"<<endl;
	//perform (-a2)*x
	vhdl << tab << declare("prod_a2_x",coeffStorageSizes[2] + sizeOfX + 1) << " <= " << use("lowX") << " * " << use("a2") << ";" << endl;
	nextCycle();////
	
	//sign-extend and pad a1 for a1 - (-a2*x)
	vhdl << tab << declare("ext_a1_pad", 1 + coeffStorageSizes[1] + keepBits) << " <= " << "\"0\" & " 
	                                                                                    << use("data")<<range(coeffStorageSizes[0] + coeffStorageSizes[1] - 1, coeffStorageSizes[0]) << " & "
	                                                                                    << zeroGenerator(keepBits, 0) << ";" << endl;
	vhdl << tab << declare("ext_prod_a2_x_pad", 1 + coeffStorageSizes[1] + keepBits) << " <= " << zeroGenerator(9,0) << " & " 
	                                                                                      << use("prod_a2_x")<< range(coeffStorageSizes[2] + sizeOfX, coeffStorageSizes[2] + sizeOfX - (1 + coeffStorageSizes[1] + keepBits - 9 )+1)
	                                                                                      << ";"<<endl;
	vhdl << tab << declare("neg_ext_prod_a2_x_pad", 1 + coeffStorageSizes[1] + keepBits) << " <= not(" << use("ext_prod_a2_x_pad") << ");"<<endl;
	//subtract (-a2*x) from a1, => instantiate IntAdder
	IntAdder* add1 = new IntAdder(target, 1 + coeffStorageSizes[1] + keepBits);
	oplist.push_back(add1);
	
	inPortMap(add1,"X", "ext_a1_pad");
	inPortMap(add1,"Y", "neg_ext_prod_a2_x_pad");
	inPortMapCst(add1,"Cin", "'1'");
	outPortMap(add1,"R", "add1Res");
	vhdl << instance(add1, "Adder_a1_prod_x_a2");
	
	syncCycleFromSignal("add1Res"); 
	
	nextCycle();//////////////////
	//perform the multiplication between x and ( a1 + a2x )    
	vhdl << tab << declare("prod_x_prev_prod", (sizeOfX+1)+ coeffStorageSizes[1] + keepBits ) << " <= " << use("lowX") << " * " 
	                                                                                         << use("add1Res")<<range( coeffStorageSizes[1] + keepBits-1, 0)<<";" <<endl;
	
	nextCycle();/////////////////                                                                               
	//compose the operands for the addition a0 + [ prev_computation ]
	vhdl << tab << declare ("right_term", 1 + coeffStorageSizes[0]) << " <= " << zeroGenerator(1 + (1+coeff_msb[0])+1-msb_x , 0)  << " & "
         << use("prod_x_prev_prod")<<range((sizeOfX+1)+ coeffStorageSizes[1] -1 + keepBits, keepBits + (sizeOfX+1)+ coeffStorageSizes[1] - (coeffStorageSizes[0]- (-msb_x+1+(1+ coeff_msb[0])))) << ";" <<endl; 
	vhdl << tab << declare ("left_term", 1+	coeffStorageSizes[0]) << " <= " << " \"0\" & " << use("data") << range(coeffStorageSizes[0]-1,0) << ";" << endl;
	
	IntAdder * add2 =  new IntAdder(target, 1+	coeffStorageSizes[0]);
	oplist.push_back(add2);
	inPortMap(add2,"X", "left_term");
	inPortMap(add2,"Y", "right_term");
	inPortMapCst(add2,"Cin", "'0'");
	outPortMap(add2,"R", "add2Res");
	vhdl << instance(add2, "FinalAdder");
	
	syncCycleFromSignal("add2Res"); 
	
	if (!correctlyRounded){
		vhdl << declare("norm_bit",1) << " <= " << use("add2Res") << "(" << coeffStorageSizes[0] << ")"<<";"<<endl;
	
		vhdl << declare("finalFrac", wF) << " <= " << use("add2Res") << range(coeffStorageSizes[0]-1, coeffStorageSizes[0]-wF) << " when " << use("norm_bit") <<"='1' else "
			                                       << use("add2Res") << range(coeffStorageSizes[0]-2, coeffStorageSizes[0]-wF-1) << ";" << endl;
		vhdl << declare("finalExp", wE) << " <= " << use("exp_addition")<<range(wE,1) << " + " << use("norm_bit")<<";"<<endl;

		vhdl << tab << "-- sign/exception handling" << endl;
		vhdl << tab << "with " << use("exnsX") << " select" <<endl
			  << tab << tab <<  declare("exnR", 2) << " <= " << "\"01\" when \"010\", -- positive, normal number" << endl
			  << tab << tab << use("exnsX") << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
			  << tab << tab << "\"11\" when others;"  << endl;

		vhdl << tab << "R <= "<<use("exnR")<<" & "<< use("exnsX") <<"(0) & " << use("finalExp") << " & " << use("finalFrac")<< ";" << endl; 
	}else{
		vhdl << declare("norm_bit",1) << " <= " << use("add2Res") << "(" << coeffStorageSizes[0] << ")"<<";"<<endl;
		vhdl << declare("preSquareFrac", wF+2) << " <= " << use("add2Res") << range(coeffStorageSizes[0], coeffStorageSizes[0]-wF-1) << " when " << use("norm_bit") <<"='1' else "
			                                           << use("add2Res") << range(coeffStorageSizes[0]-1, coeffStorageSizes[0]-wF-2) << ";" << endl;
		vhdl << declare("preSquareExp", wE) << " <= " << use("exp_addition") << range(wE,1) << " + " << use("norm_bit")<<";"<<endl;
	    vhdl << declare("preSquareConcat", 1 + wE + wF+1) << " <= " << zeroGenerator(1,0) << " & " << use("preSquareExp") << " & " << use("preSquareFrac")<<range(wF,0)<<";"<<endl;;
	    
	    IntAdder *predictAdder = new IntAdder(target, 1 + wE + wF+1);
	    oplist.push_back(predictAdder);

	    inPortMap(predictAdder,"X",use("preSquareConcat"));	
	    inPortMapCst(predictAdder,"Y",zeroGenerator(1 + wE + wF+1,0));
	    inPortMapCst(predictAdder,"Cin","'1'");
	    outPortMap(predictAdder,"R","my_predictor");
	    vhdl << tab << instance(predictAdder, "Predictor");
	    	    
	    IntSquarer *iSquarer = new IntSquarer(target,wF+2);
	    oplist.push_back(iSquarer);
	    
	    inPortMap(iSquarer, "X", use("preSquareFrac"));
	    outPortMap(iSquarer,"R", "sqrResult");
	    vhdl << instance(iSquarer,"FractionSquarer");
	    
	    syncCycleFromSignal("sqrResult");
	    
	    vhdl << tab << declare("approxSqrtXSqr", 2*(wF+2) + 1) << " <= " << zeroGenerator(1,0) << " & " << use("sqrResult")<<";"<<endl;
	    vhdl << tab << declare("realXFrac", 2*(wF+2) + 1) << " <= " << "( \"10\" & not(" <<  use("fX") << ") & not(" << zeroGenerator(2*(wF+2) + 1-2-wF ,0)<<")) when " << use("OddExp") <<"='0' else "<<endl;
	    vhdl << tab << tab << "( \"110\" & not(" <<  use("fX") << ") & not(" << zeroGenerator(2*(wF+2) + 1-2-wF-1,0)<<"));"<<endl;
	    
	    IntAdder *myIntAdd = new IntAdder(target, 2*(wF+2) + 1);
	    oplist.push_back(myIntAdd);
	    
	    inPortMap(myIntAdd,"X",use("approxSqrtXSqr"));	
	    inPortMap(myIntAdd,"Y",use("realXFrac"));
	    inPortMapCst(myIntAdd,"Cin","'1'");
	    outPortMap(myIntAdd,"R","my_add_result");
	    vhdl << tab << instance(myIntAdd, "Comparator");

	    syncCycleFromSignal("my_add_result");
	    vhdl << tab << declare("greater",1) << " <= " << use("my_add_result")<<of(2*(wF+2))<<";"<<endl;

		syncCycleFromSignal("my_predictor");    

		vhdl << tab << "-- sign/exception handling" << endl;
		vhdl << tab << "with " << use("exnsX") << " select" <<endl
			  << tab << tab <<  declare("exnR", 2) << " <= " << "\"01\" when \"010\", -- positive, normal number" << endl
			  << tab << tab << use("exnsX") << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
			  << tab << tab << "\"11\" when others;"  << endl;

		vhdl << tab << "R <= (" << use("exnR") << " & " << use("exnsX") <<"(0) & " << use("preSquareConcat")<<range(wE + wF,1) << ") when " << use("greater")<<"='1' else "<<endl;;
		vhdl << tab << tab << "(" << use("exnR") << " & " << use("exnsX") <<"(0) & " << use("my_predictor")<<range(wE + wF,1) << ");"<<endl;
	    
	}
	
	}
	////////////////////////////////////////////////////////////////////////////////////
	else {
		// Digit-recurrence implementation recycled from FPLibrary

		vhdl << tab << declare("fX", wF) << " <= X" << range(wF-1, 0) << "; -- fraction"  << endl; 
		vhdl << tab << declare("eRn0", wE) << " <= \"0\" & X" << range(wE+wF-1, wF+1) << "; -- exponent" << endl;
		vhdl << tab << declare("xsX", 3) << " <= X"<< range(wE+wF+2, wE+wF) << "; -- exception and sign" << endl;

		vhdl << tab << declare("eRn1", wE) << " <= eRn0 + (\"00\" & " << rangeAssign(wE-3, 0, "'1'") << ") + X(" << wF << ");" << endl;

		vhdl << tab << declare(join("w",wF+3), wF+4) << " <= \"111\" & fX & \"0\" when X(" << wF << ") = '0' else" << endl
			  << tab << "       \"1101\" & fX;" << endl;
		//		vhdl << tab << declare(join("d",wF+3)) << " <= '0';" << endl;
		vhdl << tab << declare(join("s",wF+3),1) << " <= '1';" << endl;

		for(int step=1; step<=wF+2; step++) {
			int i = wF+3-step; // to have the same indices as FPLibrary 
			vhdl << tab << "-- Step " << i << endl;
			string di = join("d", i);
			string xi = join("x", i);
			string wi = join("w", i);
			string wip = join("w", i+1);
			string si = join("s", i);
			string sip = join("s", i+1);
			string zs = join("zs", i);
			string ds = join("ds", i);
			string xh = join("xh", i);
			string wh = join("wh", i);
			vhdl << tab << declare(di) << " <= "<< use(wip) << "("<< wF+3<<");" << endl;
			vhdl << tab << declare(xi,wF+5) << " <= "<< use(wip) << " & \"0\";" << endl;
			vhdl << tab << declare(zs,step+1) << " <= \"0\" & " << use(sip) << ";" << endl;
			vhdl << tab << declare(ds,step+3) << " <= " << zs << range(step,1) << "& (not " << di << ") & " << di << " & \"1\";" << endl;
			vhdl << tab << declare(xh,step+3) << " <= " << xi << range(wF+4, wF+2-step) << ";" << endl;
			vhdl << tab << "with " << di << " select" << endl
				  << tab << tab <<  declare(wh, step+3) << " <= " << xh << " - " << ds << " when '0'," << endl
				  << tab << tab << "      " << xh << " + " << ds << " when others;" << endl;
			vhdl << tab << declare(wi, wF+4) << " <= " << wh << range(step+1,0);
			if(step <= wF+1) 
				vhdl << " & " << xi << range(wF+1-step, 0) << ";" << endl;  
			else
				vhdl << ";" << endl; 
			vhdl << tab << declare(si, step+1) << " <= ";
			if(step==1)
				vhdl << "not " << di << " & '1';"<< endl; 
			else
				vhdl << use(sip) << range(step-1,1) << " & not " << di << " & '1';"<< endl; 
				
			nextCycle();
		}
		vhdl << tab << declare("d0") << " <= " << use("w1") << "(" << wF+3 << ") ;" << endl;
		vhdl << tab << declare("fR", wF+4) << " <= " << use("s1") << range(wF+2, 1) << " & not d0 & '1';" << endl;

		// end of component FPSqrt_Sqrt in fplibrary
		vhdl << tab << "-- normalisation of the result, removing leading 1" << endl;
		vhdl << tab <<  "with fR(" << wF+3 << ") select" << endl
			  << tab << tab << declare("fRn1", wF+2) << " <= fR" << range(wF+2, 2) << " & (fR(1) or fR(0)) when '1'," << endl
			  << tab << tab << "        fR" <<range(wF+1, 0) << "                    when others;" << endl;
		vhdl << tab << declare("round") << " <= fRn1(1) and (fRn1(2) or fRn1(0)) ; -- round  and (lsb or sticky) : that's RN, tie to even" << endl;

		nextCycle();
		
		vhdl << tab << declare("Rn1", wE+wF) << " <= " << use("eRn1") << " & " << use("fRn1") << range(wF+1, 2) << ";" << endl;
		vhdl << tab << declare("Rn2", wE+wF) << " <= Rn1 + (" << rangeAssign(wE+wF-2, 0, "'0'") << " & " << use("round")  << "); -- never overflows for sqrt" << endl;
		
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



