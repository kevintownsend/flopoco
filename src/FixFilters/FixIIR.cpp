#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixIIR.hpp"

#include "sollya.h"

#include "ShiftReg.hpp"
#include "FixSOPC.hpp"

#include "ConstMult/FixRealKCM.hpp"


using namespace std;

namespace flopoco {

	FixIIR::FixIIR(Target* target, int msbOut_, int lsbOut_, double H_, vector<string> coeffb_, vector<string> coeffa_, map<string, double> inputDelays) : 
		Operator(target), msbOut(msbOut_), lsbOut(lsbOut_), H(H_), coeffb(coeffb_), coeffa(coeffa_)
	{
		srcFileName="FixIIR";
		setCopyrightString ( "Louis Beseme, Florent de Dinechin (2014)" );
		useNumericStd_Unsigned();


		ostringstream name;
		name << "FixIIR_"<< msbOut << "_"<< lsbOut << "_uid" << getNewUId();
		setNameWithFreq( name.str() );

		n = coeffb.size();
		m = coeffa.size();
		


		//manage the critical path
		setCriticalPath(getMaxInputDelays(inputDelays));

		addInput("X", 1+-lsbOut, true);


		// guard bits for a faithful result
		g= intlog2(2*H*(n+m)); 
		REPORT(INFO, "g=" << g);

		hugePrec = 10*(1+msbOut+-lsbOut+g);
		currentIndexA=0;
		currentIndexB=0;
		for (int i = 0; i<m; i++)
		{
			mpfr_init2 (yHistory[i], hugePrec);
			mpfr_set_d(yHistory[i], 0.0, GMP_RNDN);
		}
		for(int i=0; i<n; i++) {
			xHistory[i]=0;
		}


		for (int i=0; i< n; i++)
		{
			// parse the coeffs from the string, with Sollya parsing
			sollya_obj_t node;
			
			node = sollya_lib_parse_string(coeffb[i].c_str());
			// If conversion did not succeed (i.e. parse error)
			if(node == 0)
			{
				ostringstream error;
				error << srcFileName << ": Unable to parse string " << coeffb[i] << " as a numeric constant" << endl;
				throw error.str();
			}

			mpfr_init2(mpcoeffb[i], 10000);
			sollya_lib_get_constant(mpcoeffb[i], node);
			if(mpfr_get_d(mpcoeffb[i], GMP_RNDN) < 0)
				coeffsignb[i] = 1;
			else
				coeffsignb[i] = 0;
			
			mpfr_abs(mpcoeffb[i], mpcoeffb[i], GMP_RNDN);
				
		}

		for (int i=0; i< m; i++)
		{
			// parse the coeffs from the string, with Sollya parsing
			sollya_obj_t node;
			
			node = sollya_lib_parse_string(coeffa[i].c_str());
			// If conversion did not succeed (i.e. parse error)
			if(node == 0)
			{
				ostringstream error;
				error << srcFileName << ": Unable to parse string " << coeffa[i] << " as a numeric constant" << endl;
				throw error.str();
			}

			mpfr_init2(mpcoeffa[i], 10000);
			sollya_lib_get_constant(mpcoeffa[i], node);
			
			if(mpfr_get_d(mpcoeffa[i], GMP_RNDN) < 0)
				coeffsigna[i] = 1;
			else
				coeffsigna[i] = 0;
			
			mpfr_abs(mpcoeffa[i], mpcoeffa[i], GMP_RNDN);

		}
		

		wO = (msbOut - lsbOut) + 1; //1 + sign  ; 

		

		int size = wO + g ; 
		REPORT(INFO, "Sum size is: "<< size );

		
		//compute the guard bits from the KCM mulipliers
		int wInKCM_B = 1 -lsbOut;	//1 sign bit + p bit 
		int lsbOutKCM = lsbOut-g;
		double targetUlpError = 1.0;
		int guardBitsKCM_B = FixRealKCM::neededGuardBits(target, wInKCM_B, lsbOutKCM, targetUlpError);

		int wInKCM_A = 1+1+msbOut-lsbOut+g;
		int guardBitsKCM_A = FixRealKCM::neededGuardBits(target, wInKCM_A, lsbOutKCM, targetUlpError);

		int guardBitsKCM = max(guardBitsKCM_A, guardBitsKCM_B);

		// size += guardBitsKCM; // sign + overflow  bits on the left, guard bits + guard bits from KCMs on the right
		REPORT(INFO, "Sum size with KCM guard bits is: "<< size+guardBitsKCM);


		//Pour info
		REPORT(INFO, "guardBitsKCM part B = "<< guardBitsKCM_B);
		REPORT(INFO, "guardBitsKCM part A = "<< guardBitsKCM_A);

		//Shift register for the left part // TODO size
		ShiftReg *shiftRegB = new ShiftReg(target, 1-lsbOut , n, inputDelays);


		addSubComponent(shiftRegB);
		inPortMap(shiftRegB, "X", "X");

		for(int i = 0; i<n; i++) {
			outPortMap(shiftRegB, join("Xd", i), join("Yb", i));
		}

		vhdl << instance(shiftRegB, "shiftRegB");

		
		//Shift register for the right part
		setCycleFromSignal("X");

		ShiftReg *shiftRegA = new ShiftReg(target, wInKCM_A, m, inputDelays);

		addSubComponent(shiftRegA);
		declare("Rtmp", wO+g, false, Signal::registeredWithAsyncReset);
		nextCycle();
		inPortMap(shiftRegA, "X", "Rtmp");

		for(int i = 0; i<m; i++) {
			outPortMap(shiftRegA, join("Xd", i), join("Ya", i));
		}

		vhdl << instance(shiftRegA, "shiftRegA");
		setCycle(0);


		target->setNotPipelined(); //the following parts of the circuit will be combinatorial
		setCombinatorial();

		
		if (!target->plainVHDL())
		{
			//create the bitheap that computes the sum
			bitHeapB = new BitHeap(this, size+guardBitsKCM_B);
			bitHeapA = new BitHeap(this, size+guardBitsKCM_A);

			for (int i=0; i<n; i++) 
			{
				// TODO possible mem leak here? The pointer is lost, do we keep pointers to the subtables?
				// Multiplication: instantiating a KCM object. It will add bits also to the right of lsbOutKCM
				new FixRealKCM(this,				// the envelopping operator
														target, 	// the target FPGA
														getSignalByName(join("Yb",i)),
														true, 		// signed
														-1, 		// input MSB, but one sign bit will be added
														lsbOut, 		// input LSB weight
														lsbOutKCM, 		// output LSB weight -- the output MSB is computed out of the constant
														coeffb[i], 	// pass the string unmodified
														bitHeapB		// pass the reference to the bitheap that will accumulate the intermediary products
													);
			}


			for (int i=0; i<m; i++) 
			{
				// Multiplication: instantiating a KCM object. It will add bits also to the right of lsbOutKCM
				new FixRealKCM(this,				// the envelopping operator
														target, 	// the target FPGA
														getSignalByName(join("Ya",i)),
														true, 		// signed
														msbOut, 		// input MSB, but one sign bit will be added
														lsbOut-g, 		// input LSB weight
														lsbOutKCM, 		// output LSB weight -- the output MSB is computed out of the constant
														coeffa[i], 	// pass the string unmodified
														bitHeapA		// pass the reference to the bitheap that will accumulate the intermediary products
													);
			}


			//rounding to p+g - add 1/2 ulps
			if (guardBitsKCM_B>0)
			{
				bitHeapB->addConstantOneBit(guardBitsKCM-1);
			}
			if (guardBitsKCM_A>0)
			{
				bitHeapA->addConstantOneBit(guardBitsKCM-1);
			}

			//compress the bitheap
			bitHeapB -> generateCompressorVHDL();
			bitHeapA -> generateCompressorVHDL();

			vhdl << tab << "Rtmp" << " <= " << bitHeapB-> getSumName() << range(size+guardBitsKCM_B-1, guardBitsKCM_B) << " + " <<bitHeapA-> getSumName() << range(size+guardBitsKCM_A-1, guardBitsKCM_A) <<";" << endl;

			
		}
		else
		{
			setCycleFromSignal("Yb0");
			vhdl << tab << declare("S0", size, false, Signal::registeredWithAsyncReset) << " <= " << zg(size) << ";" << endl;
			REPORT(INFO, "cycle Yb0 is : "<<getCurrentCycle());
			for (int i=0; i<n; i++)
			{
				//manage the critical path
				setCycleFromSignal(join("Yb",i));

				// Multiplication: instantiating a KCM object. 
				FixRealKCM* mult = new FixRealKCM(target, 
												  true, // signed
												  -1, // input MSB, but one sign bit will be added
												  lsbOut, // input LSB weight
												  lsbOutKCM, // output LSB weight -- the output MSB is computed out of the constant
												  coeffb[i] // pass the string unmodified
												  );
				addSubComponent(mult);
				inPortMap(mult,"X", join("Yb", i));
				outPortMap(mult, "R", join("Pb", i));
				vhdl << instance(mult, join("mult", i));

				//manage the critical path
				syncCycleFromSignal(join("Pb", i));
				syncCycleFromSignal(join("S", i));
				manageCriticalPath(target->adderDelay(size));

				// Addition
				int pSize = getSignalByName(join("Pb", i))->width();
				vhdl << tab << declare(join("S", i+1), size, false, Signal::registeredWithAsyncReset) << " <= " <<  join("S",i);
				if(coeffsignb[i] == 1)
					vhdl << " - (" ;
				else
					vhdl << " + (" ;
				if(size>pSize) 
					vhdl << "("<< size-1 << " downto " << pSize<< " => "<< join("Pb",i) << of(pSize-1) << ")" << " & " ;
				vhdl << join("Pb", i) << ");" << endl;
			}

			for (int i=0; i<m; i++)
			{
				//manage the critical path. It's ugly but seems to be working
				setCycleFromSignal(join("Ya",i));
				previousCycle();
				getSignalByName(join("Ya",i))->setCycle(getCurrentCycle());

				// Multiplication: instantiating a KCM object. 
				FixRealKCM* mult = new FixRealKCM(target, 
												  true, // signed
												  msbOut, // input MSB, but one sign bit will be added
												  lsbOut-g, // input LSB weight
												  lsbOutKCM, // output LSB weight -- the output MSB is computed out of the constant
												  coeffa[i] // pass the string unmodified
												  );
				addSubComponent(mult);
				inPortMap(mult,"X", join("Ya", i));
				outPortMap(mult, "R", join("Pa", i));
				vhdl << instance(mult, join("mult", n+i));

				//manage the critical path
				manageCriticalPath(target->adderDelay(size));

				// Addition
				int pSize = getSignalByName(join("Pa", i))->width();
				vhdl << tab << declare(join("S", n+i+1), size, false, Signal::registeredWithAsyncReset) << " <= " <<  join("S",n+i);
				if(coeffsigna[i] == 1)
					vhdl << " - (" ;
				else
					vhdl << " + (" ;
				if(size>pSize) 
					vhdl << "("<< size-1 << " downto " << pSize<< " => "<< join("Pa",i) << of(pSize-1) << ")" << " & " ;
				vhdl << join("Pa", i) << ");" << endl;
			}

			//manage the critical path
			syncCycleFromSignal(join("S", n+m));
			manageCriticalPath(target->adderDelay(wO+1));


			// setSequential();
			
			// setCycleFromSignal(join("S", n+m));

			vhdl << tab << "Rtmp <= " << join("S", n+m) << ";" << endl;


			
		}
		setSequential();
		target->setPipelined();
		addOutput("R", wO, 2); 
		
		syncCycleFromSignal(join("Yb", n-1));
		nextCycle();

		//Adding one half ulp to obtain correct rounding
		vhdl << tab << declare("R_int", wO+1) << " <= " <<  "Rtmp_d1" << range(size-1, g - 1) << " + (" << zg(wO) << " & \'1\');" << endl;
		vhdl << tab << "R <= " <<  "R_int" << range(wO, 1) << ";" << endl;
		



	};

	FixIIR::~FixIIR(){

	};


	void FixIIR::emulate(TestCase * tc){

		mpz_class sx;

		sx = tc->getInputValue("X"); 		// get the input bit vector as an integer
		xHistory[currentIndexB] = sx;

		mpfr_t x, t, u, s;
		mpfr_init2 (x, 1-lsbOut);
		mpfr_init2 (t, hugePrec);
		mpfr_init2 (u, hugePrec);
		mpfr_init2 (s, hugePrec);	

		mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0

		for (int i=0; i< n; i++)
		{
			sx = xHistory[(currentIndexB+n-i)%n];		// get the input bit vector as an integer		
			sx = bitVectorToSigned(sx, 1-lsbOut); 						// convert it to a signed mpz_class		
			mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 				// convert this integer to an MPFR; this rounding is exact
			mpfr_div_2si (x, x, -lsbOut, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact

			mpfr_mul(t, x, mpcoeffb[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter

			if(coeffsignb[i]==1)
				mpfr_neg(t, t, GMP_RNDN); 

			mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
			
		}

		for (int i=0; i<m; i++)
		{

			mpfr_mul(u, yHistory[(currentIndexA+m-i-1)%m], mpcoeffa[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter

			if(coeffsigna[i]==1)
				mpfr_neg(u, u, GMP_RNDN); 

			mpfr_add(s, s, u, GMP_RNDN); 							// same comment as above

		}

		mpfr_set(yHistory[currentIndexA], s, GMP_RNDN);


		// now we should have in s the (exact in most cases) sum
		// round it up and down

		// make s an integer -- no rounding here 
		mpfr_mul_2si (s, s, -lsbOut, GMP_RNDN);

		
		// We are waiting until the first meaningful value comes out of the IIR

		mpz_class rdz, ruz;

		mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
		rdz=signedToBitVector(rdz, wO);
		tc->addExpectedOutput ("R", rdz);

		mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
		ruz=signedToBitVector(ruz, wO);
		tc->addExpectedOutput ("R", ruz);

		
		mpfr_clears (x, t, u, s, NULL);

		currentIndexB = (currentIndexB +1)%n; // We use a circular buffer to store the inputs
		currentIndexA = (currentIndexA +1)%m;


	};

	void FixIIR::buildStandardTestCases(TestCaseList* tcl){};

}
