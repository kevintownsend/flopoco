#include <iostream>
#include <sstream>
#include "gmp.h"
#include "mpfr.h"

#include "sollya.h"

#include "FixedPointFIRBH.hpp"

#include "ConstMult/FixRealKCMBH.hpp"

using namespace std;
namespace flopoco{


	FixedPointFIRBH::FixedPointFIRBH(Target* target, int p_, vector<string> coeff_) : Operator(target), p(p_), coeff(coeff_) {
		srcFileName="FixedPointFIRBH";
					
		ostringstream name;
		name<<"FixedPointFIRBH_"<<p<<"_uid"<<getNewUId(); 
		setName(name.str()); 
	
		setCopyrightString("Florent de Dinechin, Matei Istoan (2013)");		

		n=coeff.size();

		for (int i=0; i< n; i++)
			addInput(join("X",i), 1+p); // sign + p bits, from weights -1 to -p
		

		ostringstream clist;
		for (int i=0; i< n; i++)
			clist << "    " << coeff[i] << ", ";
		REPORT(INFO, "Building a "<< n << "-tap FIR of precision "<< p << " for coefficients " << clist.str());

		// guard bits for a faithful result
		int g= 1+ intlog2(n-1); 
		REPORT(INFO, "g=" << g);

		mpfr_t absCoeff, sumAbsCoeff;
		mpfr_init2 (absCoeff, 10*(1+p));
		mpfr_init2 (sumAbsCoeff, 10*(1+p));
		mpfr_set_d (sumAbsCoeff, 0.0, GMP_RNDN);
		
		for(int i=0; i< n; i++)
		{
			// parse the coeffs from the string, with Sollya parsing
			sollya_node_t node;
			mpfr_t mpC;
			
			node = parseString(coeff[i].c_str());
			// If conversion did not succeed (i.e. parse error)
			if(node == 0)
			{
				ostringstream error;
				error << srcFileName << ": Unable to parse string " << coeff[i] << " as a numeric constant" << endl;
				throw error.str();
			}

			mpfr_init2(mpC, 10000);
			setToolPrecision(10000);
			evaluateConstantExpression(mpC, node,  getToolPrecision());
			
			mpfr_init_set(mpcoeff[i], mpC, GMP_RNDN);
			if(mpfr_get_d(mpcoeff[i], GMP_RNDN) < 0)
				coeffsign[i] = 1;
			else
				coeffsign[i] = 0;
			
			mpfr_abs(mpcoeff[i], mpcoeff[i], GMP_RNDN);
				
			// Accumulate the absolute values
			mpfr_add(sumAbsCoeff, sumAbsCoeff, mpcoeff[i], GMP_RNDU);
		}
		
		// now sumAbsCoeff is the max value that the filter can take.
		double sumAbs = mpfr_get_d(sumAbsCoeff, GMP_RNDU); // just to make the following loop easier
		int leadingBit=0;
		while(sumAbs >= 2.0)
		{
			sumAbs *= 0.5;
			leadingBit++;
		}
		while(sumAbs < 1.0)
		{
			sumAbs *= 2.0;
			leadingBit--;
		}
		REPORT(INFO, "Worst-case weight of MSB of the result is "<< leadingBit);

		wO = 1+ (leadingBit - (-p)) + 1; //1 + sign  ; 

		addOutput("R", wO, 2); // sign + 


#if 0		
		ostringstream dlist;
		for (int i=0; i< n; i++)
		{
			char *ptr;
			mp_exp_t exp;
			ptr = mpfr_get_str (0,  &exp, 10, 0, mpcoeff[i], GMP_RNDN);
			dlist << "  ." << ptr  << "e"<< exp <<  ", ";
			if(exp>0)
				THROWERROR("Coefficient #" << i << " ("<< coeff[i]<<") larger than one in absolute value")
			mpfr_free_str(ptr);
		}
		REPORT(DEBUG, "After conversion to MPFR: " << dlist.str());
#endif		


		int size = 1 + (leadingBit - (-p) +1) + g; // sign + overflow  bits on the left, guard bits on the right
		REPORT(INFO, "Sum size is: "<< size);
		
		//compute the guard bits from the KCM mulipliers
		int guardBitsKCM = 0;
		
		for(int i=0; i<n; i++)
		{
			int wIn = p + 1;	//p bits + 1 sign bit
			int lsbOut = -p-g;
			double targetUlpError = 1.0;
			int temp = FixRealKCMBH::neededGuardBits(target, wIn, lsbOut, targetUlpError);
			
			if(temp > guardBitsKCM)
				guardBitsKCM = temp;
		}
		
		size += guardBitsKCM; // sign + overflow  bits on the left, guard bits + guard bits from KCMs on the right
		REPORT(INFO, "Sum size with KCM guard bits is: "<< size);
		
		//create the bitheap that computes the sum
		bitHeap = new BitHeap(this, size);
		
		for (int i=0; i<n; i++) 
		{
			// Multiplication: instantiating a KCM object.
			FixRealKCMBH* mult = new FixRealKCMBH(this,			// the envelopping operator
													target, 	// the target FPGA
													getSignalByName(join("X",i)),
													-p, 		// input LSB weight
													-1, 		// input MSB, but one sign bit will be added
													true, 		// signed
													-p-g, 		// output LSB weight -- the output MSB is computed out of the constant
													coeff[i], 	// pass the string unmodified
													bitHeap		// pass the reference to the bitheap that will accumulate the intermediary products
												);
		}
		
		//rounding - add 1/2 ulps
		bitHeap->addConstantOneBit(g+guardBitsKCM-1);
		
		//compress the bitheap
		bitHeap -> generateCompressorVHDL();
		
		vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << range(size-1, g+guardBitsKCM) << ";" << endl;
			
	};
	
	
	void FixedPointFIRBH::emulate(TestCase * tc)
	{
		// Not completely safe: we compute everything on 10 times the required precision, and hope that rounding this result is equivalent to rounding the exact result

		mpfr_t x, t, s, rd, ru;
		mpfr_init2 (x, 1+p);
		mpfr_init2 (t, 10*(1+p));
		mpfr_init2 (s, 10*(1+p));
		mpfr_init2 (rd, 1+p);
		mpfr_init2 (ru, 1+p);
	
		mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0

		for (int i=0; i<n; i++)
		{
			mpz_class sx = tc->getInputValue(join("X", i));			// get the input bit vector as an integer
			sx = bitVectorToSigned(sx, 1+p);						// convert it to a signed mpz_class
			mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD);				// convert this integer to an MPFR; this rounding is exact
			mpfr_div_2si (x, x, p, GMP_RNDD);						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact

			mpfr_mul(t, x, mpcoeff[i], GMP_RNDN);					// Here rounding possible, but precision used is ridiculously high so it won't matter
					
			if(coeffsign[i]==1)
				mpfr_neg(t, t, GMP_RNDN);

			mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
		}

		// now we should have in s the (exact in most cases) sum
		// round it up and down

		// make s an integer -- no rounding here
		mpfr_mul_2si (s, s, p, GMP_RNDN);

		mpz_class rdz, ruz;

		mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD);					// there can be a real rounding here
		rdz=signedToBitVector(rdz, wO);
		tc->addExpectedOutput ("R", rdz);

		mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU);					// there can be a real rounding here	
		ruz=signedToBitVector(ruz, wO);
		tc->addExpectedOutput ("R", ruz);
		
		mpfr_clears (x, t, s, rd, ru, NULL);
	}






	// please fill me with regression tests or corner case tests
	void FixedPointFIRBH::buildStandardTestCases(TestCaseList * tcl) {
		TestCase *tc;

		// first few cases to check emulate()
		// All zeroes
		tc = new TestCase(this);
		for(int i=0; i<n; i++)
			tc->addInput(join("X",i), mpz_class(0) );
		emulate(tc);
		tcl->add(tc);

		// All ones (0.11111)
		tc = new TestCase(this);
		for(int i=0; i<n; i++)
			tc->addInput(join("X",i), (mpz_class(1)<<p)-1 );
		emulate(tc);
		tcl->add(tc);

		// n cases with one 0.5 and all the other 0s
		for(int i=0; i<n; i++){
			tc = new TestCase(this);
			for(int j=0; j<n; j++){
				if(i==j)
					tc->addInput(join("X",j), (mpz_class(1)<<(p-1)) );
				else
					tc->addInput(join("X",j), mpz_class(0) );					
			}
			emulate(tc);
			tcl->add(tc);
		}
	}
}
