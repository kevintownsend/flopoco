#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixFIR.hpp"

#include "ShiftReg.hpp"
#include "FixSOPC.hpp"

using namespace std;

namespace flopoco {

	FixFIR::FixFIR(Target* target, int lsb_, vector<string> coeff_, bool useBitheap_, map<string, double> inputDelays) : 
		Operator(target), p(-lsb_), coeff(coeff_), useBitheap(useBitheap_)
	{
		srcFileName="FixFIR";
		setCopyrightString ( "Louis Beseme, Florent de Dinechin (2014)" );
		useNumericStd_Unsigned();

		if(p<1) {
			THROWERROR("Can't build an architecture for this value of LSB")
		}
		ostringstream name;
		name << "FixFIR_"<< p << "_uid" << getNewUId();
		setNameWithFreq( name.str() );

		n=coeff.size();

		// initialize stuff for emulate
		for(int i=0; i<=n; i++) {
			xHistory[i]=0;
		}
		currentIndex=0;

		//manage the critical path
		setCriticalPath(getMaxInputDelays(inputDelays));

		addInput("X", 1+p, true);

		ShiftReg *shiftReg = new ShiftReg(target, 1+p, n, inputDelays);

		addSubComponent(shiftReg);
		inPortMap(shiftReg, "X", "X");

		for(int i = 0; i<n; i++) {
			outPortMap(shiftReg, join("Xd", i), join("Y", i));
		}

		vhdl << instance(shiftReg, "shiftReg");


		// REPORT(INFO,getCycleFromSignal("Y0", false));
		// REPORT(INFO,getCurrentCycle());
#if 0
		syncCycleFromSignal("Y0");
#else 
		setCycle(0);
#endif
		mpfr_set_default_prec(10000);

		FixSOPC *fixSOPC = new FixSOPC(target, -p, coeff, useBitheap, inputDelays);
		for (int i=0; i<n; i++) {
			mpfr_init_set(mpcoeff[i], fixSOPC->mpcoeff[i], GMP_RNDN);
			coeffsign[i] = fixSOPC->coeffsign[i];
		}
		wO = fixSOPC->wO;

		addSubComponent(fixSOPC);

		for(int i=0; i<n; i++) {
			inPortMap(fixSOPC, join("X",i), join("Y", i));
		}

		outPortMap(fixSOPC, "R", "Rtmp");

		vhdl << instance(fixSOPC, "fixSOPC");
		syncCycleFromSignal("Rtmp");

		addOutput("R", fixSOPC->wO, true);
		vhdl << "R <= Rtmp;" << endl;

	};

	FixFIR::~FixFIR(){

	};

	void FixFIR::emulate(TestCase * tc){

#if 1
		mpz_class sx;

		sx = tc->getInputValue("X"); 		// get the input bit vector as an integer
		xHistory[currentIndex] = sx;

		mpfr_t x, t, s, rd, ru;
		mpfr_init2 (x, 1+p);
		mpfr_init2 (t, 10*(1+p));
		mpfr_init2 (s, 10*(1+p));
		mpfr_init2 (rd, 1+p);
		mpfr_init2 (ru, 1+p);		
		mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0
		
		for (int i=0; i< n; i++)	{
			sx = xHistory[(currentIndex+n-i)%n];
			sx = bitVectorToSigned(sx, 1+p); 						// convert it to a signed mpz_class
			mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 	 // convert this integer to an MPFR; this rounding is exact
			mpfr_div_2si (x, x, p, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact
			
			mpfr_mul(t, x, mpcoeff[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter

			if(coeffsign[i]==1)
				mpfr_neg(t, t, GMP_RNDN); 
			
			mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
			}

			// now we should have in s the (exact in most cases) sum
			// round it up and down

			// make s an integer -- no rounding here
			mpfr_mul_2si (s, s, p, GMP_RNDN);

			mpz_class rdz, ruz;

			mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
			rdz=signedToBitVector(rdz, wO);
			tc->addExpectedOutput ("R", rdz);
			// tc->addExpectedOutput ("R", rdz);

			mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
			ruz=signedToBitVector(ruz, wO);
			tc->addExpectedOutput ("R", ruz);

			mpfr_clears (x, t, s, rd, ru, NULL);
			
			currentIndex=(currentIndex+1)%n; //  circular buffer to store the inputs

	
#else
		static int idx = 0;
		static bool full = false; 							// set to true when the fir start to output valid data (after n input) 
		static TestCase * listTC [10000]; // should be enough for everybody


		listTC[idx] = tc;

		if(n == 1)					// if the fir part has only one tap we don't wait to get the output
			full = true; 

		// We are waiting until the first meaningful value comes out of the FIR
		if (full) {
			mpfr_t x, t, s, rd, ru;
			mpfr_init2 (x, 1+p);
			mpfr_init2 (t, 10*(1+p));
			mpfr_init2 (s, 10*(1+p));
			mpfr_init2 (rd, 1+p);
			mpfr_init2 (ru, 1+p);		

			mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0


			int k = idx; // We start to sum from the last input

			for (int i=0; i< n; i++)
			{

				mpz_class sx = listTC[k]->getInputValue("X"); 		// get the input bit vector as an integer
				sx = bitVectorToSigned(sx, 1+p); 						// convert it to a signed mpz_class
				mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 				// convert this integer to an MPFR; this rounding is exact
				mpfr_div_2si (x, x, p, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact

				mpfr_mul(t, x, mpcoeff[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter

				if(coeffsign[i]==1)
					mpfr_neg(t, t, GMP_RNDN); 

				mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
			
				k = (k+1)%n;	
			}

			k = (k-1+n)%n; //to get the corresponding testCase to the outputed value

			// now we should have in s the (exact in most cases) sum
			// round it up and down

			// make s an integer -- no rounding here
			mpfr_mul_2si (s, s, p, GMP_RNDN);

			mpz_class rdz, ruz;

			mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
			rdz=signedToBitVector(rdz, wO);
			listTC[k]->addExpectedOutput ("R", rdz);
			// tc->addExpectedOutput ("R", rdz);

			mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
			ruz=signedToBitVector(ruz, wO);
			listTC[k]->addExpectedOutput ("R", ruz);

			mpfr_clears (x, t, s, rd, ru, NULL);
		}
		
		idx = (idx-1+n)%n; // We use a circular buffer to store the inputs

		if (idx ==  1) {
			full = true;
		}

#endif
	};




}
	
