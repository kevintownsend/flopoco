#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicAtan2.hpp"

using namespace std;

namespace flopoco{


	//The wIn+1 below is for consistency with FixSinCos and FixSinOrCos interfaces.
	// TODO possibly fix all the code instead... This would enable sharing emulate() etc.
 
	CordicAtan2::CordicAtan2(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{

		int stage;
		srcFileName="CordicAtan2";
		setCopyrightString ( "Matei Istoan, Florent de Dinechin (2012-...)" );

		ostringstream name;
		name << "CordicAtan2_"<< w_;
		if(target->isPipelined())
			name  <<"_f" << target->frequencyMHz();
		else 
			name << "_comb";
		name << "_uid" << getNewUId();
		setName( name.str() );
		


#define ROUNDED_ROTATION 1 // 0:trunc 

#if ROUNDED_ROTATION
		REPORT(DEBUG, "Using rounded rotation trick");
#endif



#if 0
		//error analysis
		double eps;  //error in ulp
		eps=0.5; //initial rounding of kfactor
		double shift=0.5;
		for(stage=1; stage<=maxIterations; stage++){
#if ROUNDED_ROTATION
			eps = eps + eps*shift + 0.5; // 0.5 assume rounding in the rotation.
#else
			eps = eps + eps*shift + 1.0; // 1.0 assume truncation in the rotation.
#endif
			shift *=0.5;
		}

		if (reducedIterations == 1) {
			eps+=2; // two multiplications in sequence, each truncated
		}

		eps+=1; // the final neg-by-not
		REPORT(DEBUG, "Error analysis computes eps=" << eps << " ulps (before final rounding)");

		// guard bits depend only on the number of iterations
		g = 1+(int) ceil(log2(eps)); // +1 for the final rounding 
#endif
		
		// *********    internal precision and fixed-point alignment **************

 		// The input is as follows:
		// s has weight 2^0
		// q has weight 2^-1
		// o has weight 2^-2
		// Purpose: have the LSB of z, cosine, sine at weight -w
		// This means that the Cos and Sin datapath will have w+1 bits
		//   (sign at weight zero)
		// while the Z datapath starts on w-1 bits (sign bit at weight  -2,
		//  and decreasing, so the invariant is: sign bit at weight -stage-1)

		// everybody needs many digits of Pi
		mpfr_init2(constPi, 10*w);
		mpfr_const_pi( constPi, GMP_RNDN);

		//compute the scale factor		
		mpfr_init2(scale, w+2);
		mpfr_set_d(scale, -1.0, GMP_RNDN);           // exact
		mpfr_mul_2si(scale, scale, -w+1, GMP_RNDN); // exact
		mpfr_add_d(scale, scale, 1.0, GMP_RNDN);     // exact
		REPORT(DEBUG, "scale=" << printMPFR(scale, 15));
		

		// declaring inputs
		addInput  ( "X"  , w, true );
		addInput  ( "Y"  , w, true );

		// declaring output
		addOutput  ( "A"  , w, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		manageCriticalPath( target->lutDelay());
		
	};


	CordicAtan2::~CordicAtan2(){
		mpfr_clears (kfactor, constPi, NULL);		
	 };


	void CordicAtan2::emulate(TestCase * tc) 
	{
		mpfr_t x,y,a;
		mpfr_init2(x, 10*w);
		mpfr_init2(y, 10*w);
		mpfr_init2(a, 10*w);

		mpz_class az;

		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		// interpret as signed two'ss complement
		if (1==(svX >> (w-1))) // sign bit
			svX -= (1<<w);
		if (1==(svY >> (w-1))) // sign bit
			svY -= (1<<w);
		/* Compute correct value */
		
		mpfr_set_z (x, svX.get_mpz_t(), GMP_RNDN); //  exact
		mpfr_set_z (y, svY.get_mpz_t(), GMP_RNDN); //  exact
	
		mpfr_atan2(a, x, y, GMP_RNDN); // a between -pi and pi
		mpfr_div(a, a, constPi, GMP_RNDN); // a between -1 and 1

		// Now convert a to fix point
		// Align to fix point by adding 6 -- if we just add 4 there is a 1-bit shift in case a<0
		mpfr_add_d(a, a, 6.0, GMP_RNDN);
		mpfr_mul_2si (a, a, w, GMP_RNDN); // exact scaling 

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDD); // there can be a real rounding here
		az -= mpz_class(6)<<w;
 		tc->addExpectedOutput ("A", az);

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDU); // there can be a real rounding here
		az -= mpz_class(6)<<w;
 		tc->addExpectedOutput ("A", az);
		
		// clean up
		mpfr_clears (x,y,a, NULL);		
	}






	void CordicAtan2::buildStandardTestCases(TestCaseList * tcl) 
	{
	}


	

}

