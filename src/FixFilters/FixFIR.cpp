#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixFIR.hpp"

#include "ShiftReg.hpp"
#include "FixSOPC.hpp"

using namespace std;

namespace flopoco {

	FixFIR::FixFIR(Target* target, int lsbInOut_) : 
		Operator(target), lsbInOut(lsbInOut_)	{	}

	FixFIR::FixFIR(Target* target, int lsbInOut_, vector<string> coeff_, map<string, double> inputDelays) : 
		Operator(target), lsbInOut(lsbInOut_), coeff(coeff_)
	{
		srcFileName="FixFIR";
		setCopyrightString ( "Louis Beseme, Florent de Dinechin (2014)" );

		ostringstream name;
		name << "FixFIR_uid" << getNewUId();
		setNameWithFreq( name.str() );

		buildVHDL();
	};

	

	// The method that does the work once coeff[] is known
	void FixFIR::buildVHDL(){
#if 0
		n=coeff.size();

		useNumericStd_Unsigned();
		if(p<1) {
			THROWERROR("Can't build an architecture for this value of LSB")
		}
		addInput("X", 1+p, true);

		ShiftReg *shiftReg = new ShiftReg(getTarget(), 1+p, n);

		addSubComponent(shiftReg);
		inPortMap(shiftReg, "X", "X");

		for(int i = 0; i<n; i++) {
			outPortMap(shiftReg, join("Xd", i), join("Y", i));
		}

		vhdl << instance(shiftReg, "shiftReg");

		setCycle(0);
		mpfr_set_default_prec(10000);

		FixSOPC *fixSOPC = new FixSOPC(getTarget(), -p, coeff);

		wO = fixSOPC->wO;

		// Here we should set the critical path of x[0] only by some setCriticalPath();
		addSubComponent(fixSOPC);

		for(int i=0; i<n; i++) {
			inPortMap(fixSOPC, join("X",i), join("Y", i));
		}

		outPortMap(fixSOPC, "R", "Rtmp");

		vhdl << instance(fixSOPC, "fixSOPC");
		syncCycleFromSignal("Rtmp");

		addOutput("R", fixSOPC->wO, true);
		vhdl << "R <= Rtmp;" << endl;

		// initialize stuff for emulate
		for(int i=0; i<=n; i++) {
			xHistory[i]=0;
		}
		currentIndex=0;
		for (int i=0; i<n; i++) {
			mpfr_init_set(mpcoeff[i], fixSOPC->mpcoeff[i], GMP_RNDN);
			coeffsign[i] = fixSOPC->coeffsign[i];
		}
#endif
	};

	FixFIR::~FixFIR(){

	};

	void FixFIR::emulate(TestCase * tc){
		// TODO: delegate the computation to a method of FixSOPC
#if 0
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
#endif
	};




}
	
