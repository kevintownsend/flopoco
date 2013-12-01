#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicAtan2.hpp"

using namespace std;

namespace flopoco{

	// TODO Options:
	// An option to output the angle in the (atan(2^-i)) basis
	// an an option for outputting the norm of the vector as well (scaled or not)

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
		


#define ROUNDED_ROTATION 0 // 0:trunc 

#if ROUNDED_ROTATION
		REPORT(DEBUG, "Using rounded rotation trick");
#endif

		// ulp = weight of the LSB of the result is 2^(-w+1)
		// half-ulp is 2^-w
		// atan(2^-w) < 2^-w therefore after w interations,  method error will be bounded by 2^-w 
		maxIterations = w;

		//error analysis
		double eps;  //error in ulp
		eps=2; // initial neg-by-not

		double shift=0.5;
		for(stage=1; stage<=maxIterations; stage++){
#if ROUNDED_ROTATION
			eps = eps + eps*shift + 0.5; // 0.5 assume rounding in the rotation.
#else
			eps = eps + eps*shift + 1.0; // 1.0 assume truncation in the rotation.
#endif
			shift *=0.5;
		}

		eps+=1; // the final neg-by-not

		// guard bits depend only on the number of iterations
		g = 1+(int) ceil(log2(eps)); // +1 for the final rounding 
		REPORT(DEBUG, "Error analysis computes eps=" << eps << " ulps (before final rounding), hence  g=" << g);
		
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
		

		// declaring inputs. 
		// man atan2 says "atan2(y,x) is atan(y/x) so we will provide the inputs in the same order...
		addInput  ( "Y"  , w, true );
		addInput  ( "X"  , w, true );

		// declaring output
		addOutput  ( "A"  , w, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		//		manageCriticalPath( target->lutDelay());


		mpfr_t  zatan;
		mpfr_init2(zatan, 10*w);
		int sizeZ=w+g; 
		int zMSB=0;
		int zLSB=-w-g+1;
		///////////// VHDL GENERATION

		// manage symmetries: we need to reduce (X,Y) to a quadrant (-pi/4, pi/4)

		vhdl << tab << declare("sgnX") << " <= X" << of(w-1) << ";" << endl;
		vhdl << tab << declare("sgnY") << " <= X" << of(w-1) << ";" << endl;

		// TODO: replace the following with LUT-based comparators
		// and first maybe experiment with synthesis tools
		vhdl << tab << declare("XmY", w+1) << " <= (sgnX & X)-(sgnY & Y);" << endl;
		vhdl << tab << declare("XpY", w+1) << " <= (sgnX & X)+(sgnY & Y);" << endl;
 		vhdl << tab << declare("XltY") << " <= XmY" << of(w) <<";" << endl;
		vhdl << tab << declare("mYltX") << " <= XpY" << of(w) <<";" << endl;

		vhdl << tab << "-- quadrant will also be the angle to add at the end" <<endl;
		vhdl << tab << declare("quadrant", 2) << " <= " << endl;
		vhdl << tab << tab << "\"00\"  when (not sgnX and not XltY and     mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"01\"  when (not sgnY and     XltY and     mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"10\"  when (    sgnX and     XltY and not mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"11\";"    << endl;

		vhdl << tab << declare("mX", w+g) << " <= (not X) & " << og(g)<< ";  -- negation by not, implies one ulp error." << endl;
		vhdl << tab << declare("mY", w+g) << " <= (not Y) & " << og(g)<< ";  -- negation by not, implies one ulp error. " << endl;

		vhdl << tab << declare("X1", w+g) << " <= " << endl;
		vhdl << tab << tab << " X when quadrant=\"00\"   else " << endl;
		vhdl << tab << tab << " Y when quadrant=\"01\"   else " << endl;
		vhdl << tab << tab << "mX when quadrant=\"10\"   else " << endl;
		vhdl << tab << tab << "mY;"    << endl;

		vhdl << tab << declare("Y1", w+g) << " <= " << endl;
		vhdl << tab << tab << " Y when quadrant=\"00\"   else " << endl;
		vhdl << tab << tab << "mX when quadrant=\"01\"   else " << endl;
		vhdl << tab << tab << "mY when quadrant=\"10\"   else " << endl;
		vhdl << tab << tab << " X;"    << endl;


		vhdl << tab << declare("Z1", sizeZ) << " <= " << zg(sizeZ)<< ";   " << endl;



		//CORDIC iterations

		for(stage=1; stage<=maxIterations; stage++){
			//shift Xin and Yin with 2^n positions to the right
			// From there on X is always positive, but Y may be negative and thus need sign extend
			vhdl << tab << declare(join("sgnY", stage))  << " <= " <<  join("Y", stage)  <<  of(w+g-1) << ";" << endl; 
			vhdl << tab << declare(join("XShift", stage), w+g) << " <= " << zg(stage) << " & X" << stage << range(w+g-1, stage) << ";" <<endl;			
			vhdl << tab << declare(join("YShift", stage), w+g) 
			     << " <= " << rangeAssign(w+g-1, w+g-stage, join("sgnY", stage))   
			     << " & Y" << stage << range(w+g, stage) << ";" << endl;
			
			// Critical path delay for one stage:
			// The data dependency is from one Z to the next 
			// We may assume that the rotations themselves overlap once the DI are known
			//manageCriticalPath(target->localWireDelay(w) + target->adderDelay(w) + target->lutDelay()));
			// 			manageCriticalPath(target->localWireDelay(sizeZ) + target->adderDelay(sizeZ));
			
			
			
#if ROUNDED_ROTATION  // rounding of the shifted operand, should save 1 bit in each addition
			vhdl << tab << declare(join("XShiftRoundBit", stage)) << " <= " << join("X", stage)  << of(stage-1) << ";" << endl;
			vhdl << tab << declare(join("YShiftRoundBit", stage)) << " <= " << join("Y", stage) << of(stage-1) << ";" <<endl;
			vhdl << tab << declare(join("XShiftNeg", stage), w+1) << " <= " << rangeAssign(w, 0, join("sgnY", stage)) << " xor " << join("XShift", stage)   << " ;" << endl;
			vhdl << tab << declare(join("YShiftNeg", stage), w+1) << " <= (not " << rangeAssign(w, 0, join("sgnY", stage)) << ") xor " << join("YShift", stage)   << " ;" << endl;

			vhdl << tab << declare(join("X", stage+1), w+1) << " <= " 
			     << join("X", stage) << " + " << join("YShiftNeg", stage) << " +  not (" << join("sgnY", stage) << " xor " << join("YShiftRoundBit", stage) << ") ;" << endl;

			vhdl << tab << declare(join("Y", stage+1), w+1) << " <= " 
			     << join("Y", stage) << " + " << join("XShiftNeg", stage) << " + (" << join("sgnY", stage) << " xor " << join("XShiftRoundBit", stage) << ") ;" << endl;

#else
// truncation of the shifted operand
			vhdl << tab << declare(join("X", stage+1), w+1) << " <= " 
			     << join("X", stage) << " - " << join("YShift", stage) << " when " << join("sgnY", stage) << "=\'0\' else "
			     << join("X", stage) << " + " << join("YShift", stage) << " ;" << endl;
			
			vhdl << tab << declare(join("Y", stage+1), w+1) << " <= " 
			     << join("Y", stage) << " + " << join("XShift", stage) << " when " << join("sgnY", stage) << "=\'0\' else "
			     << join("Y", stage) << " - " << join("XShift", stage) << " ;" << endl;
#endif			
			
			//create the constant signal for the arctan
			mpfr_set_d(zatan, 1.0, GMP_RNDN);
			mpfr_div_2si(zatan, zatan, stage, GMP_RNDN);
			mpfr_atan(zatan, zatan, GMP_RNDN);
			mpfr_div(zatan, zatan, constPi, GMP_RNDN);
			REPORT(DEBUG, "stage=" << stage << "  atancst=" << printMPFR(zatan, 15));		
			//create the arctangent factor to be added to Zin
									
			vhdl << tab << declare(join("atan2PowStage", stage), sizeZ) << " <= " << unsignedFixPointNumber(zatan, zMSB, zLSB) << ";" <<endl;
			
			vhdl << tab << declare(join("Z", stage+1), sizeZ) << " <= " 
					 << join("Z", stage) << " + " << join("atan2PowStage", stage) << " when " << join("sgnY", stage) << "=\'1\' else "
					 << join("Z", stage) << " - " << join("atan2PowStage", stage) << " ;" << endl;
		}
		
		// Give the time to finish the last rotation
		// manageCriticalPath( target->localWireDelay(w+1) + target->adderDelay(w+1) // actual CP delay
		//                     - (target->localWireDelay(sizeZ+1) + target->adderDelay(sizeZ+1))); // CP delay that was already added

		// reconstruction

		vhdl << tab << "A <= (quadrant & " << zg(w-2) << ")   + ("
				 << join("Z", maxIterations) << of(sizeZ-1) << " & " << join("Z", maxIterations) << of(sizeZ-1) // sign extension 
				 << " & "<<join("Z", maxIterations) << range(sizeZ-1, sizeZ-w+2) << ");" << endl;
		
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
	
		mpfr_atan2(a, y, x, GMP_RNDN); // a between -pi and pi
		mpfr_div(a, a, constPi, GMP_RNDN); // a between -1 and 1

		// Now convert a to fix point
		// Align to fix point by adding 6 -- if we just add 4 there is a 1-bit shift in case a<0
		mpfr_add_d(a, a, 6.0, GMP_RNDN);
		mpfr_mul_2si (a, a, w-1, GMP_RNDN); // exact scaling 

		mpz_class mask = (mpz_class(1)<<w) -1; 

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDD); // there can be a real rounding here
		az -= mpz_class(6)<<(w-1);
		az &= mask; 
 		tc->addExpectedOutput ("A", az);

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDU); // there can be a real rounding here
		az -= mpz_class(6)<<(w-1);
		az &= mask; 
 		tc->addExpectedOutput ("A", az);
		
		// clean up
		mpfr_clears (x,y,a, NULL);		
	}






	void CordicAtan2::buildStandardTestCases(TestCaseList * tcl) 
	{
	}


	

}

