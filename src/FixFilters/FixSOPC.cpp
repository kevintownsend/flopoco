#include <iostream>
#include <sstream>
#include "gmp.h"
#include "mpfr.h"

#include "sollya.h"

#include "FixSOPC.hpp"

#include "ConstMult/FixRealKCM.hpp"

using namespace std;
namespace flopoco{

	const int veryLargePrec = 6400;  /*6400 bits should be enough for anybody */


	FixSOPC::FixSOPC(
			Target* target_, 
			int lsbIn_,
			int lsbOut_,
			vector<string> coeff_
		) : 
			Operator(target_),
			lsbOut(lsbOut_),
			coeff(coeff_),
			g(-1),
			computeMSBOut(true),
			computeGuardBits(true),
			addFinalRoundBit(true)
	{
		n = coeff.size();
		for (int i=0; i<n; i++) {
			msbIn.push_back(0);
			lsbIn.push_back(lsbIn_);
		}
		initialize();
	}


	FixSOPC::FixSOPC(
			Target* target_,
			vector<int> msbIn_,
			vector<int> lsbIn_,
			int msbOut_,
			int lsbOut_,
			vector<string> coeff_,
			int g_
		) :
			Operator(target_),
			msbIn(msbIn_),
			lsbIn(lsbIn_),
			msbOut(msbOut_),
			lsbOut(lsbOut_),
			coeff(coeff_),
			g(g_),
			computeMSBOut(false)
	{
		n = coeff.size();
		for(int i=0; i<n; i++)
			maxAbsX.push_back(intpow2(msbIn[i]));
		
		if (g==-1)
			computeGuardBits=true;

		addFinalRoundBit = (g==0 ? false : true);

		initialize();
	}


	FixSOPC::FixSOPC(
			Target* target_,
			vector<double> maxAbsX_,
			vector<int> lsbIn_,
			int msbOut_,
			int lsbOut_,
			vector<string> coeff_,
			int g_
		) :
			Operator(target_),
			maxAbsX(maxAbsX_),
			lsbIn(lsbIn_),
			msbOut(msbOut_),
			lsbOut(lsbOut_),
			coeff(coeff_),
			g(g_),
			computeMSBOut(false)
	{
		n = coeff.size();
		for(int i=0; i<n; i++)
			msbIn.push_back(ceil(log2(maxAbsX[i])));
		
		if (g==-1)
			computeGuardBits=true;

		addFinalRoundBit = (g==0 ? false : true);

		initialize();
	}




	FixSOPC::~FixSOPC()
	{
		for (int i=0; i<n; i++) {
			mpfr_clear(mpcoeff[i]);
		}
	}

	void FixSOPC::initialize()	{
		srcFileName="FixSOPC";

		ostringstream name;
		name << "FixSOPC_uid" << getNewUId(); 
		setName(name.str()); 
	
		setCopyrightString("Matei Istoan, Louis Besème, Florent de Dinechin (2013-2015)");
		
		for (int i=0; i< n; i++)
		{
			addInput(join("X",i), msbIn[i]-lsbIn[i]+1); 
		}

		//reporting on the filter
		ostringstream clist;
		for (int i=0; i< n; i++)
		{
			clist << "    " << coeff[i] << ", ";
		}
		REPORT(INFO, "Building a " << n << "-tap FIR; lsbOut=" << lsbOut << 
				" for coefficients " << clist.str());


		if(computeGuardBits) {
			// guard bits for a faithful result
			g = intlog2(n-1);
			REPORT(INFO, "g=" << g);
		}


		for (int i=0; i< n; i++)	{
				// parse the coeffs from the string, with Sollya parsing
				sollya_obj_t node;

				node = sollya_lib_parse_string(coeff[i].c_str());
				// If conversion did not succeed (i.e. parse error)
				if(node == 0)	{
						ostringstream error;
						error << srcFileName << ": Unable to parse string " << 
							coeff[i] << " as a numeric constant" << endl;
						throw error.str();
					}

				mpfr_init2(mpcoeff[i], veryLargePrec);
				sollya_lib_get_constant(mpcoeff[i], node);
				sollya_lib_clear_obj(node);
			}


		if(computeMSBOut) {
			mpfr_t sumAbsCoeff, absCoeff;
			mpfr_init2 (sumAbsCoeff, veryLargePrec);
			mpfr_init2 (absCoeff, veryLargePrec);
			mpfr_set_d (sumAbsCoeff, 0.0, GMP_RNDN);

			for (int i=0; i< n; i++)	{
			// Accumulate the absolute values
				mpfr_abs(absCoeff, mpcoeff[i], GMP_RNDU);
				mpfr_add(sumAbsCoeff, sumAbsCoeff, absCoeff, GMP_RNDU);
			}

			// now sumAbsCoeff is the max value that the SOPC can take.
			double sumAbs = mpfr_get_d(sumAbsCoeff, GMP_RNDU); // just to make the following loop easier
			REPORT(INFO, "sumAbs=" << sumAbs);
			msbOut=1;
			while(sumAbs>=2.0)		{
					sumAbs*=0.5;
					msbOut++;
				}
			while(sumAbs<1.0)	{
					sumAbs*=2.0;
					msbOut--;
				}
			REPORT(INFO, "Worst-case weight of MSB of the result is " << msbOut);
			mpfr_clears(sumAbsCoeff, absCoeff, NULL);
		}

		addOutput("R", msbOut-lsbOut+1);

		int sumSize = 1 + msbOut - lsbOut  + g ;
		REPORT(DETAILED, "Sum size is: "<< sumSize );

		//compute the guard bits from the KCM mulipliers, and take the max
		int guardBitsKCM = 0;
		int lsbOutKCM = lsbOut-g; // we want each KCM to be faithful to this ulp
		double targetUlpError = 1.0;

		for(int i=0; i<n; i++)		{
			int wInKCM = msbIn[i]-lsbIn[i]+1-g;	//p bits + 1 sign bit

			int temp = FixRealKCM::neededGuardBits(
					getTarget(), 
					wInKCM, 
					targetUlpError,
					coeff[i],
					lsbIn[i],
					lsbOut
				);			

			if(temp > guardBitsKCM)
				guardBitsKCM = temp;
		}

		sumSize += guardBitsKCM;
		REPORT(DETAILED, "Sum size with KCM guard bits is: "<< sumSize);

		if(!getTarget()->plainVHDL())
		{
			//create the bitheap that computes the sum
			bitHeap = new BitHeap(this, sumSize);

			for (int i=0; i<n; i++)	{
				// Multiplication: instantiating a KCM object. It will add bits also to the right of lsbOutKCM
				new FixRealKCM(
						this,				// the envelopping operator
						 getTarget(), 	// the target FPGA
						 getSignalByName(join("X",i)),
						 true, 		// signed
						 msbIn[i], 	
						 lsbIn[i], 		// input LSB weight
						 lsbOutKCM, 		// output LSB weight -- the output MSB is computed out of the constant
						 coeff[i], 	// pass the string unmodified
						 bitHeap,	// pass the reference to the bitheap that will accumulate the intermediary products
						 lsbOutKCM - guardBitsKCM 
					 );
			}

			//compress the bitheap
			bitHeap -> generateCompressorVHDL();
			
			vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << 
				range(sumSize-1, g+guardBitsKCM) << ";" << endl;

		}

		else
		{
			 THROWERROR("Sorry, plainVHDL doesn't work at the moment for FixSOPC. Somebody has to fix it and remove this message" );
			// Technically if you comment the line above it generates non-correct VHDL

			// All the KCMs in parallel
			for(int i=0; i< n; i++)	{
				FixRealKCM* mult = new FixRealKCM(getTarget(), 
						true, // signed
						msbIn[i]-1, // input MSB,TODO one sign bit will be added by KCM because it has a non-standard interface. To be fixed someday
						lsbIn[i], // input LSB weight
						lsbOut-g, // output LSB weight -- the output MSB is computed out of the constant
						coeff[i] // pass the string unmodified
					);
				addSubComponent(mult);
				inPortMap(mult,"X", join("X", i));
				outPortMap(mult, "R", join("P", i));
				vhdl << instance(mult, join("mult", i));
			}
			// Now advance to their output, and build a pipelined rake
			syncCycleFromSignal("P0");
			vhdl << tab << declare("S0", sumSize) << " <= " << zg(sumSize) << ";" << endl;
			for(int i=0; i< n; i++)		{
				//manage the critical path
				manageCriticalPath(getTarget()->adderDelay(sumSize));
				// Addition
				int pSize = getSignalByName(join("P", i))->width();
				vhdl << tab << declare(join("S", i+1), sumSize) << " <= " <<  join("S",i);
				vhdl << " + (" ;
				if(sumSize>pSize)
					vhdl << "("<< sumSize-1 << " downto " << pSize<< " => "<< join("P",i) << of(pSize-1) << ")" << " & " ;
				vhdl << join("P", i) << ");" << endl;
			}

			// Rounding costs one more adder to add the half-ulp round bit.
			// This could be avoided by pushing this bit in one of the KCM tables
			syncCycleFromSignal(join("S", n));
			manageCriticalPath(getTarget()->adderDelay(msbOut-lsbOut+1));
			vhdl << tab << declare("R_int", sumSize+1) << " <= " <<  join("S", n) << range(sumSize-1, sumSize-(msbOut-lsbOut+1)-1) << " + (" << zg(sumSize) << " & \'1\');" << endl;
			vhdl << tab << "R <= " <<  "R_int" << range(sumSize, 1) << ";" << endl;
		}


	};





	// Function that factors the work done by emulate() of FixFIR and the emulate() of FixSOPC
	pair<mpz_class,mpz_class> FixSOPC::computeSOPCForEmulate(vector<mpz_class> inputs) {
		// Not completely safe: we compute everything on veryLargePrec, and hope that rounding this result is equivalent to rounding the exact result
		mpfr_t t, s, rd, ru;
		mpfr_init2 (t, veryLargePrec);
		mpfr_init2 (s, veryLargePrec);
		mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0
		for (int i=0; i< n; i++)	{
			mpfr_t x;
			mpz_class sx = bitVectorToSigned(inputs[i], 1+msbIn[i]-lsbIn[i]); 						// convert it to a signed mpz_class
			mpfr_init2 (x, 1+msbIn[i]-lsbIn[i]);
			mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 				// convert this integer to an MPFR; this rounding is exact
			mpfr_mul_2si (x, x, lsbIn[i], GMP_RNDD); 						// multiply this integer by 2^-lsb to obtain a fixed-point value; this rounding is again exact

			mpfr_mul(t, x, mpcoeff[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter
			mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
			mpfr_clears (x, NULL);
		}

		// now we should have in s the (very accurate) sum
		// round it up and down

		// make s an integer -- no rounding here
		mpfr_mul_2si (s, s, -lsbOut, GMP_RNDN);

		mpfr_init2 (rd, 1+msbOut-lsbOut);
		mpfr_init2 (ru, 1+msbOut-lsbOut);

		mpz_class rdz, ruz;

		mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
		rdz=signedToBitVector(rdz, 1+msbOut-lsbOut);

		mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here
		ruz=signedToBitVector(ruz, 1+msbOut-lsbOut);

		mpfr_clears (t, s, rd, ru, NULL);

		return make_pair(rdz, ruz);
	}





	void FixSOPC::emulate(TestCase * tc) {
		vector<mpz_class> inputs;
		for (int i=0; i< n; i++)	{
			mpz_class sx = tc->getInputValue(join("X", i)); 		// get the input bit vector as an integer
			inputs.push_back(sx);
		}
		pair<mpz_class,mpz_class> results = computeSOPCForEmulate(inputs);

		tc->addExpectedOutput ("R", results.first);
		tc->addExpectedOutput ("R", results.second);
	}

	OperatorPtr FixSOPC::parseArguments(Target *target, vector<string> &args) {
		int lsbIn;
		UserInterface::parseInt(args, "lsbIn", &lsbIn);
		int lsbOut;
		UserInterface::parseInt(args, "lsbOut", &lsbOut);
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
		
		return new FixSOPC(target, lsbIn, lsbOut, input);
	}

	void FixSOPC::registerFactory(){
		UserInterface::add("FixSOPC", // name
											 "A fix-point Sum of Product by Constants.",
											 "FiltersEtc", // categories
											 "",
											 "lsbIn(int): input's last significant bit;\
                        lsbOut(int): output's last significant bit;\
                        coeff(string): colon-separated list of real coefficients using Sollya syntax. Example: coeff=1.234567890123:sin(3*pi/8);",
											 "",
											 FixSOPC::parseArguments
											 ) ;
	}


	// please fill me with regression tests or corner case tests
	void FixSOPC::buildStandardTestCases(TestCaseList * tcl) {

#if 0
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
#endif
	}

}
