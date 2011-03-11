/*
 * A correctly-rounded multiplier by an arbitrary real constant for FloPoCo

  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License, 2008-2010.
*/

#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../FixedPointFunctions/HOTBM/sollya.h" // TODO : fix upstream Sollya, or fix in FloPoCo
#include "../utils.hpp"
#include "../Operator.hpp"
#include "FPConstMult.hpp"
#include "FPConstMultParser.hpp"
#include "../FPNumber.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;



	FPConstMultParser::FPConstMultParser(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, int wF_C, string _constant):
		FPConstMult(target, wE_in, wF_in, wE_out, wF_out), 
		cst_width(wF_C+1),
		constant (_constant) 
	{
		sollya_node_t node;
		mpfr_t mpR;
		mpz_t zz;

		bool periodic_constant;

		// Ugly interface hack !
		if(wF_C==1)
			periodic_constant= true;

		srcFileName="FPConstMultParser";
		/* Convert the input string into a sollya evaluation tree */
		node = parseString(constant.c_str());	/* If conversion did not succeed (i.e. parse error) */
		if (node == 0) {
			ostringstream error;
			error << srcFileName << ": Unable to parse string "<< constant << " as a numeric constant" <<endl;
			throw error.str();
		}
		//     flopoco -verbose=1 FPConstMultParser 8 23 8 23 30 "1/7"


		mpfr_inits(mpR, NULL);
		evaluateConstantExpression(mpR, node,  getToolPrecision());
		
		REPORT(DEBUG, "Constant evaluates to " << mpfr_get_d(mpR, GMP_RNDN));

		if(periodic_constant) {
			int periodSize, headerSize; 
			mpz_class periodicPattern, header; // The two values to look for
			// Evaluate to a very large number of bits, then look for a period
			// evaluate to 2000 bits
			#define EVAL_PRECISION 10000
			int zSize=2*EVAL_PRECISION;
			mpfr_set_prec(mpR, zSize);
			evaluateConstantExpression(mpR, node, zSize);
			int exponent;
			int maxPeriodSize=128;
			mpz_t z_mpz;
			mpz_class z, x0, x1;
			mpz_init(z_mpz);
			exponent = mpfr_get_z_exp(z_mpz, mpR);
			z = mpz_class(z_mpz);
			mpz_clear(z_mpz);
			REPORT(DEBUG, "Constant is " << z.get_str(2)); 
			REPORT(DETAILED, "Looking for a period");
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
					REPORT(DEBUG, "i=" <<i << " x0=" << x0.get_str(2) << " x1=" << x1.get_str(2) );
					t = t >> periodSize;
					x1 = t-((t >> periodSize)<<periodSize);
					i++;
				}
				if(i==max_number_of_periods)
					found=true;
				else
					periodSize ++;
			}

			if(found){
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
					//cerr << "DDDDD " << headerSize << " " <<   header.get_str(2) << " " << x1.get_str(2) << endl;
					for (int i=0; i<periodSize; i++){
						//cerr << "zzzzzzz " <<  i << "  "<< x1.get_str(2) << " " << periods[i].get_str(2) << endl;
						if (x1==periods[i]) {
							header_found=true;
							periodicPattern=periods[i];
						}
					}
					if(!header_found) {
						headerSize++;
						header = z >> (zSize-headerSize); 
						mpz_class a, b, c;
						// a = (z - (header << (zSize-headerSize)));
						// b = a >> ((zSize-headerSize - periodSize));
						// cerr << "KKKKKKKKKKKKKK " << ((zSize-headerSize - periodSize))  << "  "<< a.get_str(2) << endl << b.get_str(2) << endl;
						x1 = (z - (header << (zSize-headerSize)))  >> (zSize-headerSize - periodSize) ;
					}
				}
				if(header_found) {
					REPORT(INFO, "Found header " << header.get_str(2) << " and periodic pattern " << periodicPattern.get_str(2) << " of size " << periodSize);
				}
				else {
					REPORT(INFO, "Found no header for periodic pattern " << periodicPattern.get_str(2) << " of size " << periodSize);
				}
			}			
			else
				REPORT(INFO, "Periodic pattern not found");
			exit(0);
		}
		else{
			evaluateConstantExpression(mpR, node,  cst_width);
		}
		// Now convert mpR into exponent + integral significand

		if(0==mpfr_get_d(mpR, GMP_RNDN)) {
			REPORT(INFO, "building a multiplier by 0, it will be easy");
			vhdl  << tab << "r <= " << rangeAssign(wE_out+wF_out+1, 0, "'0'") << ";"<<endl;
		}
		else {
			// sign
			if(mpfr_sgn(mpR)<0) {
				mpfr_neg(mpR, mpR, GMP_RNDN);
				cst_sgn=1;
			} 
			else {
				cst_sgn=0;
			}

			// compute exponent and mantissa
			cst_exp_when_mantissa_1_2 = mpfr_get_exp(mpR) - 1; //mpfr_get_exp() assumes significand in [1/2,1)  
			cst_exp_when_mantissa_int = cst_exp_when_mantissa_1_2 - cst_width + 1; 
			mpfr_init2(mpfr_cst_sig, cst_width);
			mpfr_div_2si(mpfr_cst_sig, mpR, cst_exp_when_mantissa_1_2, GMP_RNDN);
			REPORT(INFO, "mpfr_cst_sig  = " << mpfr_get_d(mpfr_cst_sig, GMP_RNDN));

			// Build the corresponding FPConstMult.

			// initialize mpfr_xcut_sig = 2/cst_sig, will be between 1 and 2
			mpfr_init2(mpfr_xcut_sig, 32*(cst_width+wE_in+wE_out)); // should be accurate enough
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
			mpfr_mul_2si(mpfr_cst_sig,  mpfr_cst_sig, cst_width, GMP_RNDN);
   
			// It should now be an int; cast it into a mpz, then a mpz_class 
			mpz_init2(zz, cst_width);
			mpfr_get_z(zz, mpfr_cst_sig, GMP_RNDN);
			cst_sig = mpz_class(zz);
			mpz_clear(zz);
			REPORT(DETAILED, "mpzclass cst_sig = " << cst_sig);

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
			
			// build the name
			ostringstream name; 
			name <<"FPConstMult_"<<(cst_sgn==0?"":"M") <<cst_sig<<"b"
				  <<(cst_exp_when_mantissa_int<0?"M":"")<<abs(cst_exp_when_mantissa_int)
				  <<"_"<<wE_in<<"_"<<wF_in<<"_"<<wE_out<<"_"<<wF_out;
			uniqueName_=name.str();

			// cleaning up
			mpfr_clears(mpR, mpfr_xcut_sig, xcut_wF, mpfr_cst_sig, NULL);

			icm = new IntConstMult(target, wF_in+1, cst_sig);
			oplist.push_back(icm);
			
			setup();
		}
	}



	FPConstMultParser::~FPConstMultParser() {
		// TODO but who cares really
	}


}
#endif //HAVE_SOLLYA
