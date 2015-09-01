#include <iostream>
#include <iomanip>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"
#include "sollya.h"

#include "FixFIR.hpp"

#include "ShiftReg.hpp"

using namespace std;

namespace flopoco {

	const int veryLargePrec = 6400;  /*6400 bits should be enough for anybody */

	FixFIR::FixFIR(Target* target, int lsbInOut_, bool rescale_) :
		Operator(target), lsbInOut(lsbInOut_), rescale(rescale_)	{	}



	FixFIR::FixFIR(Target* target, int lsbInOut_, vector<string> coeff_, bool rescale_, map<string, double> inputDelays) :
		Operator(target), lsbInOut(lsbInOut_), coeff(coeff_), rescale(rescale_)
	{
		srcFileName="FixFIR";
		setCopyrightString ( "Louis Besème, Florent de Dinechin (2014)" );

		ostringstream name;
		name << "FixFIR_uid" << getNewUId();
		setNameWithFreq( name.str() );

		buildVHDL();
	};



	// The method that does the work once coeff[] is known
	void FixFIR::buildVHDL(){
		n=coeff.size();

		useNumericStd_Unsigned();
		if(-lsbInOut<1) {
			THROWERROR("Can't build an architecture for this value of lsbInOut: " << lsbInOut)
		}
		addInput("X", 1-lsbInOut, true);

		ShiftReg *shiftReg = new ShiftReg(getTarget(), 1-lsbInOut, n);

		addSubComponent(shiftReg);
		inPortMap(shiftReg, "X", "X");

		for(int i = 0; i<n; i++) {
			outPortMap(shiftReg, join("Xd", i), join("Y", i));
		}

		vhdl << instance(shiftReg, "shiftReg");

		setCycle(0);

		if (rescale) {
			// Most of this code is copypasted from SOPC.
			// parse the coeffs from the string, with Sollya parsing
			mpfr_t sumAbsCoeff;
			mpfr_init2 (sumAbsCoeff, 1-lsbInOut);
			mpfr_set_d (sumAbsCoeff, 0.0, GMP_RNDN);

			for (int i=0; i< n; i++)	{
				mpfr_t absCoeff;
				mpfr_init2 (absCoeff, veryLargePrec);
				sollya_obj_t node;
				node = sollya_lib_parse_string(coeff[i].c_str());
				// If conversion did not succeed (i.e. parse error)
				if(node == 0)	{
					ostringstream error;
					error << srcFileName << ": Unable to parse string " << coeff[i] << " as a numeric constant" << endl;
					throw error.str();
				}
				mpfr_init2(absCoeff, veryLargePrec);
				sollya_lib_get_constant(absCoeff, node);
				sollya_lib_clear_obj(node);
				// Accumulate the absolute values
				mpfr_abs(absCoeff, absCoeff, GMP_RNDU);
				mpfr_add(sumAbsCoeff, sumAbsCoeff, absCoeff, GMP_RNDU);
				mpfr_clears(absCoeff, NULL);
			}
			// now sumAbsCoeff is the max value that the filter can take.
			double sumAbs = mpfr_get_d(sumAbsCoeff, GMP_RNDU);
			REPORT(INFO, "Scaling all the coefficients by 1/" << sumAbs);
			mpfr_clears(sumAbsCoeff, NULL);

			// now just replace each coeff with the scaled version
			for (int i=0; i< n; i++)	{
				ostringstream s;
				s << "(" << coeff[i] << ")/" << setprecision(20) << sumAbs;
				coeff[i] = s.str();
			}
		}


		fixSOPC = new FixSOPC(getTarget(), lsbInOut, lsbInOut, coeff);

		addSubComponent(fixSOPC);

		for(int i=0; i<n; i++) {
			inPortMap(fixSOPC, join("X",i), join("Y", i));
		}

		outPortMap(fixSOPC, "R", "Rtmp");

		vhdl << instance(fixSOPC, "fixSOPC");
		syncCycleFromSignal("Rtmp");

		addOutput("R", fixSOPC->msbOut - fixSOPC->lsbOut + 1,   true);
		vhdl << "R <= Rtmp;" << endl;

		// initialize stuff for emulate
		for(int i=0; i<=n; i++) {
			xHistory[i]=0;
		}
		currentIndex=0;
	};

	FixFIR::~FixFIR(){
	};



	void FixFIR::emulate(TestCase * tc){
		mpz_class sx = tc->getInputValue("X"); 		// get the input bit vector as an integer
		xHistory[currentIndex] = sx;

		// Not completely optimal in terms of object copies...
		vector<mpz_class> inputs;
		for (int i=0; i< n; i++)	{
			sx = xHistory[(currentIndex+n-i)%n];
			inputs.push_back(sx);
		}
		pair<mpz_class,mpz_class> results = fixSOPC-> computeSOPCForEmulate(inputs);

		tc->addExpectedOutput ("R", results.first);
		tc->addExpectedOutput ("R", results.second);
		currentIndex=(currentIndex+1)%n; //  circular buffer to store the inputs

#if 0
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

	OperatorPtr FixFIR::parseArguments(Target *target, vector<string> &args) {
		int lsbInOut;
		UserInterface::parseInt(args, "lsbInOut", &lsbInOut);
		bool rescale;
		UserInterface::parseBoolean(args, "rescale", &rescale);
		vector<string> input;
		string in;
		UserInterface::parseString(args, "coeff", &in);
		// tokenize a string, thanks Stack Overflow
		stringstream ss(in);
		while( ss.good() )	{
				string substr;
				getline( ss, substr, ':' );
				input.push_back( substr );
			}

		return new FixFIR(target, lsbInOut, input, rescale);
	}

	void FixFIR::registerFactory(){
		UserInterface::add("FixFIR", // name
											 "A fix-point Finite Impulse Filter generator.",
											 "FiltersEtc", // categories
											 "",
											 "lsbInOut(int): integer size in bits;\
                        rescale(bool)=false: If true, divides all coefficient by 1/sum(|coeff|);\
                        coeff(string): colon-separated list of real coefficients using Sollya syntax. Example: coeff=1.234567890123:sin(3*pi/8)",
											 "For more details, see <a href=\"bib/flopoco.html#DinIstoMas2014-SOPCJR\">this article</a>.",
											 FixFIR::parseArguments
											 ) ;
	}


}

