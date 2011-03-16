/*
   A multiplier by a floating-point constant for FloPoCo

  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License,  2008-2010.
*/


// TODO Case mantissa=1, wE_out>wE_in
// TODO standard test vectors: 1, 0, various exn, xcut borders


#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../sollya.h"
#include "../Operator.hpp"
#include "FPConstMult.hpp"
#include "../FPNumber.hpp"

using namespace std;


namespace flopoco{

extern vector<Operator*> oplist;


	// The expert version

	FPConstMult::FPConstMult(Target* target, int wE_in_, int wF_in_, int wE_out_, int wF_out_, int cstSgn_, int cst_exp_, mpz_class cst_sig_):
		Operator(target), 
		wE_in(wE_in_), wF_in(wF_in_), wE_out(wE_out_), wF_out(wF_out_), 
		cstSgn(cstSgn_), cst_exp_when_mantissa_int(cst_exp_), cst_sig(cst_sig_)
	{
		srcFileName="FPConstMult";
		ostringstream name;
		name <<"FPConstMult_"<<(cstSgn==0?"":"M") <<mpz2string(cst_sig)<<"b"<<(cst_exp_when_mantissa_int<0?"M":"")<<abs(cst_exp_when_mantissa_int)<<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
		uniqueName_=name.str();

		if(cst_sig==0) {
			REPORT(INFO, "building a multiplier by 0, it will be easy");
			vhdl  << tab << "r <= " << rangeAssign(wE_out+wF_out+1, 0, "'0'") << ";"<<endl;
		}
		else {
			// Constant normalization
			while ((cst_sig % 2) ==0) {
				REPORT(INFO, "Significand is even, normalising");
				cst_sig = cst_sig >>1;
				cst_exp_when_mantissa_int+=1;
			}
			mantissa_is_one = false;
			if(cst_sig==1) {
				REPORT(INFO, "Constant mantissa is 1, multiplying by it will be easy"); 
				mantissa_is_one = true;
			}

			cstWidth = intlog2(cst_sig);
			cst_exp_when_mantissa_1_2 = cst_exp_when_mantissa_int + cstWidth - 1; 

			// initialize mpfr constant
			mpfr_init2(mpfr_cst_sig, max(cstWidth, 2));
			mpfr_set_z(mpfr_cst_sig, cst_sig.get_mpz_t(), GMP_RNDN); // exact op
			mpfr_mul_2si(mpfr_cst_sig, mpfr_cst_sig, -(cstWidth-1), GMP_RNDN);  // exact op, gets mpfr_cst_sig in 1..2
			
			mpfr_init(mpfrC);
			mpfr_set(mpfrC, mpfr_cst_sig, GMP_RNDN);
			mpfr_mul_2si(mpfrC, mpfrC, cst_exp_when_mantissa_1_2, GMP_RNDN);
			
			if(cstSgn==1)
				mpfr_neg(mpfrC,  mpfrC, GMP_RNDN);
			REPORT(INFO, "mpfrC = " << mpfr_get_d(mpfrC, GMP_RNDN));


			if(!mantissa_is_one) {			
				// initialize mpfr_xcut_sig = 2/cst_sig, will be between 1 and 2
				mpfr_init2(mpfr_xcut_sig, 4*(cstWidth+wE_in+wE_out));
				mpfr_set_d(mpfr_xcut_sig, 2.0, GMP_RNDN);               // exaxt op
				mpfr_div(mpfr_xcut_sig, mpfr_xcut_sig, mpfr_cst_sig, GMP_RNDD);

				// now  round it up to wF_in+1 bits 
				mpfr_t xcut_wF;
				mpfr_init2(xcut_wF, wF_in+1);
				mpfr_set(xcut_wF, mpfr_xcut_sig, GMP_RNDD);
				mpfr_mul_2si(xcut_wF, xcut_wF, wF_in, GMP_RNDN);
				// It should now be an int; cast it into a mpz, then a mpz_class 
				mpz_t zz;
				mpz_init2(zz, wF_in+1);
				mpfr_get_z(zz, xcut_wF, GMP_RNDN);
				xcut_sig_rd= mpz_class(zz);

				REPORT(DETAILED, "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN));
				REPORT(DETAILED, "mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN)); 
				if(verbose) { // can't use REPORT because of GMP's C++ wrapper problem under windoze
					cerr << "xcut_sig_rd   = " << mpz2string(xcut_sig_rd) << "   ";
					printBinNumGMP(cerr, xcut_sig_rd, wF_in+1);  cerr << endl;
				}
				mpfr_clear(xcut_wF);
				mpz_clear(zz);
				// sub component
				icm = new IntConstMult(target, wF_in+1, cst_sig);
				oplist.push_back(icm);
			}


			// do all the declarations. Pushed into a method so that CRFPConstMult can inherit it
			buildVHDL();

		}

	}

















#ifdef HAVE_SOLLYA




	// The parser version
	FPConstMult::FPConstMult(Target* target, int wE_in_, int wF_in_, int wE_out_, int wF_out_, int wF_C, string constant):
		Operator(target), 
		wE_in(wE_in_), wF_in(wF_in_), wE_out(wE_out_), wF_out(wF_out_), cstWidth(wF_C)
	{
		sollya_node_t node;


		bool periodic_constant;

		// Ugly interface hack !
		if(wF_C==-1)
			periodic_constant= true;

		srcFileName="FPConstMult";
		/* Convert the input string into a sollya evaluation tree */
		node = parseString(constant.c_str());	/* If conversion did not succeed (i.e. parse error) */
		if (node == 0) {
			ostringstream error;
			error << srcFileName << ": Unable to parse string "<< constant << " as a numeric constant" <<endl;
			throw error.str();
		}
		//     flopoco -verbose=1 FPConstMultParser 8 23 8 23 30 "1/7"


		mpfr_inits(mpfrC, NULL);
		evaluateConstantExpression(mpfrC, node,  getToolPrecision());
		REPORT(DEBUG, "Constant evaluates to " << mpfr_get_d(mpfrC, GMP_RNDN));

		int expWhenVeryLargeInt;


		if(periodic_constant) {
			int periodSize, headerSize; 
			mpz_class periodicPattern, header; // The two values to look for

			// Evaluate to a very large number of bits, then look for a period
			// evaluate to 2000 bits
			#define EVAL_PRECISION 2000
			int zSize=2*EVAL_PRECISION;
			mpfr_set_prec(mpfrC, zSize);

			evaluateConstantExpression(mpfrC, node, zSize);
			int maxPeriodSize=128;
			mpz_t z_mpz;
			mpz_class z, x0, x1;
			mpz_init(z_mpz);
			expWhenVeryLargeInt = mpfr_get_z_exp(z_mpz, mpfrC);
			z = mpz_class(z_mpz);
			mpz_clear(z_mpz);
			REPORT(DETAILED, "Looking for a period in a constant which looks like " << z.get_str(2)); 
			periodSize=2;

			bool found=false;
			while (!found &&  periodSize < maxPeriodSize) {
				// first go to the middle of these 2000 bits
				mpz_class t = z >> (EVAL_PRECISION);
				// First filter
				x0 = t-((t >> periodSize)<<periodSize);
				t = t >> periodSize;
				x1 = t-((t >> periodSize)<<periodSize);
				int i=2;
 				int max_number_of_periods = (EVAL_PRECISION/maxPeriodSize) >> 1;
				REPORT(DEBUG, "Trying periodSize=" << periodSize << " max_number_of_periods=" << max_number_of_periods);
				while(x0==x1 && i<max_number_of_periods) {
					// REPORT(DEBUG, "i=" <<i << " x0=" << x0.get_str(2) << " x1=" << x1.get_str(2) );
					t = t >> periodSize;
					x1 = t-((t >> periodSize)<<periodSize);
					i++;
				}
				if(i==max_number_of_periods)
					found=true;
				else
					periodSize ++;
			}


			if(!found){
				ostringstream error;
				error << srcFileName << ": period not found" <<endl;
				throw error.str();
			}
			
			REPORT(DEBUG, "Found periodic pattern " << x0 << " (" << x0.get_str(2) << ") " << " of size "<< periodSize);
			// Now look for header bits
			// first build a table of rotations of this constant
			vector<mpz_class> periods;
			for (int i=0; i<periodSize; i++) {
				periods.push_back(x0); 
				x1 = x0>>(periodSize-1); // MSB
				x0 = ((x0-(x1 << (periodSize-1)))<<1) + x1;
				// cerr << x0.get_str(2) << endl;
			}
			// Now compare to the first bits of the constant
			headerSize=0;
			header=0;
			x1 = z >> (zSize-periodSize);
			// cerr << x1.get_str(2) << endl;
			bool header_found=false;
			while (!header_found && headerSize<maxPeriodSize) {
				for (int i=0; i<periodSize; i++){
					if (x1==periods[i]) {
						header_found=true;
						periodicPattern=periods[i];
					}
				}
				if(!header_found) {
					headerSize++;
					header = z >> (zSize-headerSize); 
					x1 = (z - (header << (zSize-headerSize)))  >> (zSize-headerSize - periodSize) ;
				}
			} // end while

			// TODO The previous is wrong

			REPORT(DETAILED, "Found header " << header.get_str(2) << " of size "<< headerSize 
			       << " and period " << periodicPattern.get_str(2) << " of size " << periodSize);

			// Now go on
			int wC = wF_out + 3; // for faithful rounding, but could come from the interface: TODO
			int r = ceil(   ((double)(wC-headerSize)) / ((double)periodSize)   ); // Needed repetitions
			REPORT(DETAILED, "wC=" << wC << ", need to repeat the period " << r << " times");
			int i = intlog2(r) -1; // 2^i < r < 2^{i+1}
			int rr = r - (1<<i);
			int j;
			if (rr==0)
				j=-1;
			else
				j= intlog2(rr-1);

			// now round up the number of needed repetitions
			if(j==-1) {
				r=(1<<i);
				REPORT(DETAILED, "... Will repeat 2^i with i=" << i);
			}
			else{
				r=(1<<i)+(1<<j);
				REPORT(DETAILED, "... Will repeat 2^i+2^j with i=" << i << " and j=" << j);
			}


			// Now rebuild the mpfrC constant (for emulate() etc)
			// First, as an integer mantissa
			cst_sig = header;
			for(int k=0; k<r; k++)
				cst_sig = (cst_sig<<periodSize) + periodicPattern;
			REPORT(DEBUG, "Constant mantissa rebuilt as " << cst_sig << " ==  " << cst_sig.get_str(2) );

			// now as an MPFR
			cstWidth = headerSize  + r*periodSize;
			mpfr_set_prec(mpfrC, cstWidth);
			// get the exponent when the mantissa is on this size
			evaluateConstantExpression(mpfrC, node, zSize); // mpfrC is a dummy here
			mpz_class dummy;
			cst_exp_when_mantissa_int = mpfr_get_z_exp(dummy.get_mpz_t(), mpfrC);
			mpfr_set_z(mpfrC, cst_sig.get_mpz_t(), GMP_RNDN);
			mpfr_mul_2si(mpfrC, mpfrC, cst_exp_when_mantissa_int, GMP_RNDN); // exact
			REPORT(DEBUG, "Constant rebuilt as " << mpfr_get_d( mpfrC, GMP_RNDN) );

			icm = new IntConstMult(target, wF_in+1, cst_sig, periodicPattern, periodSize, header, headerSize, i, j);
			oplist.push_back(icm);


		}
				// else {
				// 	ostringstream error;
				// 	error << srcFileName << ": Found no header for periodic pattern " << periodicPattern.get_str(2) << " of size " << periodSize ;
				// 	throw error.str();
				// }






		else{  // if (periodic)

			// Nonperiodic version

			if(wF_C==0) //  means: please compute wF_C for faithful rounding
				cstWidth=wF_out+3;
			
			mpfr_set_prec(mpfrC, wF_out+3);
			evaluateConstantExpression(mpfrC, node,  cstWidth);

			setup();

			if(!constant_is_zero && !mantissa_is_one) {
				icm = new IntConstMult(target, wF_in+1, cst_sig);
				oplist.push_back(icm);
			}
			
		}
		
		
		// build the name
		ostringstream name; 
		name <<"FPConstMult_"<<(cstSgn==0?"":"M") <<cst_sig<<"b"
				 <<(cst_exp_when_mantissa_int<0?"M":"")<<abs(cst_exp_when_mantissa_int)
				 <<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
		uniqueName_=name.str();
		
		// cleaning up : TODO
		//			mpfr_clears(mpfrC, mpfr_xcut_sig, xcut_wF, mpfr_cst_sig, NULL);
		
		
		buildVHDL();
	}










	// Set up the various constants out of an MPFR constant
	// The constant precision must be set up properly in 
	void FPConstMult::setup()
	{

		mpz_t zz;
		cstWidth=mpfr_get_prec(mpfrC);

		REPORT(DEBUG, "cstWidth=" << cstWidth);

		if(mpfr_zero_p(mpfrC)) {
			REPORT(INFO, "building a multiplier by 0, it will be easy");
			constant_is_zero=true;
			return;
		}

		// sign
		if(mpfr_sgn(mpfrC)<0) {
			mpfr_neg(mpfrC, mpfrC, GMP_RNDN);
			cstSgn=1;
		} 
		else {
			cstSgn=0;
		}

		// compute exponent and mantissa
		cst_exp_when_mantissa_1_2 = mpfr_get_exp(mpfrC) - 1; //mpfr_get_exp() assumes significand in [1/2,1)  
		cst_exp_when_mantissa_int = cst_exp_when_mantissa_1_2 - cstWidth + 1; 
		mpfr_init2(mpfr_cst_sig, cstWidth);
		mpfr_div_2si(mpfr_cst_sig, mpfrC, cst_exp_when_mantissa_1_2, GMP_RNDN);
		REPORT(INFO, "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN));
		
		// Build the corresponding FPConstMult.

		// initialize mpfr_xcut_sig = 2/cst_sig, will be between 1 and 2
		mpfr_init2(mpfr_xcut_sig, 32*(cstWidth+wE_in+wE_out)); // should be accurate enough
		mpfr_set_d(mpfr_xcut_sig, 2.0, GMP_RNDN);               // exaxt op
		mpfr_div(mpfr_xcut_sig, mpfr_xcut_sig, mpfr_cst_sig, GMP_RNDD);

		// now  round it down to wF_in+1 bits 
		mpfr_t xcut_wF;
		mpfr_init2(xcut_wF, wF_in+1);
		mpfr_set(xcut_wF, mpfr_xcut_sig, GMP_RNDD);
		mpfr_mul_2si(xcut_wF, xcut_wF, wF_in, GMP_RNDN);
		
		// It should now be an int; cast it into a mpz, then a mpz_class 
		mpz_init2(zz, wF_in+1);
		mpfr_get_z(zz, xcut_wF, GMP_RNDN);
		xcut_sig_rd = mpz_class(zz);
		mpz_clear(zz);
		REPORT(DETAILED, "mpfr_xcut_sig = " << mpfr_get_d(mpfr_xcut_sig, GMP_RNDN) );

		// Now build the mpz significand
		mpfr_mul_2si(mpfr_cst_sig,  mpfr_cst_sig, cstWidth, GMP_RNDN);
		
		// It should now be an int; cast it into a mpz, then a mpz_class 
		mpz_init2(zz, cstWidth);
		mpfr_get_z(zz, mpfr_cst_sig, GMP_RNDN);
		cst_sig = mpz_class(zz);
		mpz_clear(zz);
		REPORT(DETAILED, "mpzclass cst_sig = " << cst_sig);
		
		// Constant normalization
		while ((cst_sig % 2) ==0) {
			REPORT(DETAILED, "Significand is even, normalising");
			cst_sig = cst_sig >>1;
			cst_exp_when_mantissa_int += 1;
			cstWidth -= 1;
		}

		if(cst_sig==1) {
			REPORT(INFO, "Constant mantissa is 1, multiplying by it will be easy"); 
			mantissa_is_one = true;
			return;
		}
		return; // normal constant
	}
	




	FPConstMult::FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out):
		Operator(target),
		wE_in(wE_in), wF_in(wF_in), wE_out(wE_out), wF_out(wF_out) 
	{
	}
#endif //HAVE_SOLLYA


	FPConstMult::~FPConstMult() {
		// TODO but who cares really
		// clean up memory -- with ifs, because in some cases they have not been init'd
		if(mpfrC) mpfr_clear(mpfrC);
		if(mpfr_cst_sig) mpfr_clear(mpfr_cst_sig);
		if(mpfr_xcut_sig) mpfr_clear(mpfr_xcut_sig);
	}





	void FPConstMult::buildVHDL() {

		// Set up the IO signals
		addFPInput("X", wE_in, wF_in);
		addFPOutput("R", wE_out, wF_out);

		setCopyrightString("Florent de Dinechin (2007)");

		if(constant_is_zero) {
			vhdl << tab << declare("x_exn",2) << " <=  X("<<wE_in<<"+"<<wF_in<<"+2 downto "<<wE_in<<"+"<<wF_in<<"+1);"<<endl;
			vhdl << tab << declare("x_sgn") << " <=  X("<<wE_in<<"+"<<wF_in<<");"<<endl;
			
			vhdl << tab << declare("r_exn", 2) << " <=      \"00\" when ((x_exn = \"00\") or (x_exn = \"01\"))  -- zero"<<endl 
					 << tab << "         else \"11\" ;-- 0*inf = 0*NaN = NaN" << endl;

			vhdl  << tab << "R <= r_exn & x_sgn & " << rangeAssign(wE_out+wF_out-1, 0, "'0'") << ";"<<endl;
			
			return;
		}

		// bit width of constant exponent
		int wE_cst=intlog2(abs(cst_exp_when_mantissa_1_2));
		REPORT(DETAILED, "wE_cst = "<<wE_cst<<endl);
	
		// We have to compute Er = E_X - bias(wE_in) + E_C + bias(wE_R)
		// Let us pack all the constants together
		mpz_class expAddend = -bias(wE_in) + cst_exp_when_mantissa_1_2  + bias(wE_out);
		int expAddendSign=0;
		if (expAddend < 0) {
			expAddend = -expAddend;
			expAddendSign=1;
		}
		int wE_sum; // will be the max size of all the considered  exponents
		wE_sum = intlog2(expAddend);
		if(wE_in > wE_sum)
			wE_sum = wE_in;
		if(wE_out > wE_sum) 
			wE_sum = wE_out;


		vhdl << tab << declare("x_exn",2) << " <=  X("<<wE_in<<"+"<<wF_in<<"+2 downto "<<wE_in<<"+"<<wF_in<<"+1);"<<endl;
		vhdl << tab << declare("x_sgn") << " <=  X("<<wE_in<<"+"<<wF_in<<");"<<endl;
		vhdl << tab << declare("x_exp", wE_in) << " <=  X("<<wE_in<<"+"<<wF_in<<"-1 downto "<<wF_in<<");"<<endl;


		if(mantissa_is_one) {			
			vhdl << tab << "-- The mantissa of the constant is  1" << endl;
			if(wF_out == wF_in) {
				vhdl << tab << declare("r_frac", wF_out) << " <= X("<<wF_in-1 <<" downto 0);"<<endl;
			}
			else if(wF_out > wF_in){
				vhdl << tab << tab << declare("r_frac", wF_out) << " <= X("<<wF_in-1 <<" downto 0)  &  " << rangeAssign(wF_out-wF_in-1, 0, "'0'") << ";"<<endl;
				}
			else{ // wF_out < wF_in, this is a rounding of the mantissa TODO
				throw string("FPConstMult: multiplication by a power of two when  wF_out < wF_in not yet implemented, please complain to the FloPoCo team if you need it");
			}
			vhdl << tab << declare("gt_than_xcut") << " <= '0';"<<endl;
	
		}
		else{ // normal case, mantissa is not one
			vhdl << tab << declare("x_sig", wF_in+1) << " <= '1' & X("<<wF_in-1 <<" downto 0);"<<endl;
			vhdl << tab << declare("xcut_rd", wF_in+1) << " <= \""
				  << unsignedBinary(xcut_sig_rd, wF_in+1) << "\";"<<endl;
			vhdl << tab << declare("gt_than_xcut") << " <= '1' when ( x_sig("<<wF_in-1<<" downto 0) > xcut_rd("<<wF_in-1<<" downto 0) ) else '0';"<<endl;
			inPortMap  (icm, "inX", "x_sig");
			outPortMap (icm, "R","sig_prod");
			vhdl << instance(icm, "sig_mult");
			setCycleFromSignal("sig_prod"); 
			nextCycle();
			// Possibly shift the significand one bit left, and remove implicit 1 
			vhdl << tab << declare("shifted_frac",    wF_out+1) << " <= sig_prod("<<icm->rsize -2<<" downto "<<icm->rsize - wF_out-2 <<")  when gt_than_xcut = '1'"<<endl
				  << tab << "           else sig_prod("<<icm->rsize -3<<" downto "<<icm->rsize - wF_out - 3<<");"<<endl;  
			// add the rounding bit
			vhdl << tab << tab << declare("rounded_frac",   wF_out+1) << " <= (("<<wF_out <<" downto 1 => '0') & '1') + shifted_frac;"<<endl;
			vhdl << tab << tab << declare("r_frac", wF_out) << " <= rounded_frac("<<wF_out <<" downto  1);"<<endl;
		}

		// Handling signs is trivial
		if(cstSgn==0)
			vhdl << tab << declare("r_sgn") << " <= x_sgn; -- positive constant"<<endl;
		else
			vhdl << tab << declare("r_sgn") << " <= not x_sgn; -- negative constant"<<endl;

		// exponent handling
		vhdl << tab << declare("abs_unbiased_cst_exp",wE_sum+1) << " <= \""
			  << unsignedBinary(expAddend, wE_sum+1) << "\";" << endl;
		vhdl << tab << declare("r_exp_nopb",    wE_out+1) << " <= "
			  << "((" << wE_sum << " downto " << wE_in << " => '0')  & x_exp)  "
			  << (expAddendSign==0 ? "+" : "-" ) << "  abs_unbiased_cst_exp"
			  << "  +  (("<<wE_sum<<" downto 1 => '0') & gt_than_xcut);"<<endl;

		// overflow handling
		vhdl << tab << declare("overflow") << " <= " ;
		if (maxExp(wE_in) + cst_exp_when_mantissa_1_2 + 1 < maxExp(wE_out)) // no overflow can ever happen
			vhdl << "'0'; --  overflow never happens for this constant and these (wE_in, wE_out)" << endl;
		else 
			vhdl <<  "'0' when r_exp_nopb(" << wE_sum << " downto " << wE_out << ") = (" << wE_sum << " downto " << wE_out <<" => '0')     else '1';" << endl;

		// underflow handling
		vhdl << tab << declare("underflow") << " <= " ;
		if (minExp(wE_in) + cst_exp_when_mantissa_1_2 > minExp(wE_out)) // no underflow can ever happen
			vhdl << "'0'; --  underflow never happens for this constant and these (wE_in, wE_out)" << endl;
		else 
			vhdl <<  "r_exp_nopb(" << wE_sum << ");" << endl;
			 
	
		vhdl << tab << declare("r_exp", wE_out) << " <= r_exp_nopb("<<wE_out-1<<" downto 0) ;"<<endl;

		vhdl << tab << declare("r_exn", 2) << " <=      \"00\" when ((x_exn = \"00\") or (x_exn = \"01\" and underflow='1'))  -- zero"<<endl 
			  << tab << "         else \"10\" when ((x_exn = \"10\") or (x_exn = \"01\" and overflow='1'))   -- infinity" << endl
			  << tab << "         else \"11\" when  (x_exn = \"11\")                      -- NaN" << endl
			  << tab << "         else \"01\";                                          -- normal number" << endl;

		vhdl  << tab << "r <= r_exn & r_sgn & r_exp & r_frac;"<<endl;


	}






	void FPConstMult::emulate(TestCase *tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");

		/* Compute correct value */
		FPNumber fpx(wE_in, wF_in);
		fpx = svX;
		mpfr_t x, r;
		mpfr_init2(x, 1+wF_in);
		mpfr_init2(r, 1+wF_out); 
		fpx.getMPFR(x);
		mpfr_mul(r, x, mpfrC, GMP_RNDN);

		// Set outputs 
		FPNumber  fpr(wE_out, wF_out, r);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		// clean up
		mpfr_clears(x, r, NULL);

	}

}
