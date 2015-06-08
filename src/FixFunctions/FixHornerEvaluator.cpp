/*

  This file is part of the FloPoCo project
  initiated by the Aric team at Ecole Normale Superieure de Lyon
  and developed by the Socrate team at Institut National des Sciences Appliquées de Lyon

  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2008-2014.
  All rights reserved.

*/
#include <iostream>

#include "FixHornerEvaluator.hpp"
#include "IntMult/FixMultAdd.hpp"

using namespace std;


	/*



			TODO: one LSB, or a vector of LSBs?
	 Analysis of the polynomial coefficients to determine the intermediate sizes.
	 This class knows the abs max of the coefficients to do a worst-case analysis.
		 This may be pessimistic, e.g. log has alternate coeffs and doesn't grow as far as worst case suggests.
		 Managing that properly is a TODO: it impacts only efficiency, not numerical quality.
		 Maybe it can be done with a bit more Sollya.

		 This is a simplified version of the computation in the ASAP 2010 paper, simplified because x is in [-1,1)

		 WRT the ASAP paper,
		 1/we want to use FixMultAdd operations:
		 \sigma_0  =  a_d
		 \sigma_j  =  a_{d-j} + x^T_i * \sigma_{j-1}     \forall j \in \{ 1...d\}
		 p(y)      =  \sigma_d

		 2/ we have only one value of g for all the steps.
		 Main advantage is that it is simpler to maintain.

		 We still should consider the DSP granularity to truncate x to the smallest DSP-friendly size larger than |sigma_{j-1}|
		 TODO

		 So all that remains is to compute the parameters of the FixMultAdd.
		 We need
		 * size of \sigma_d: MSB is that of  a_{d-j}, plus 1 for overflows (sign extended).
							 LSB is the common LSB: lsbOut-g
		 * size of x truncated x^T_i : since it is multiplied by \sigma_{j-1}, it should have the same size
		 * weight of the MSB of the product: since y \in [-1,1), this is the MSB of \sigma_{j-1}, i.e. the MSB of a_{d-j+1}
		 ** in case of signed coeffs, the max abs value is 2^(MSB-1).
		 **   however, since we have zero-extended x in this case there is some thinking TODO
		 * size of the truncated result:
		 * weight of the LSB of a_{d-j}
		 * weight of the MSB of a_{d-j}
		 */


#define LARGE_PREC 1000 // 1000 bits should be enough for everybody

namespace flopoco{

/** This builds an architecture such as eps_finalround < 2^(lsbOut-1) and eps_round<2^(lsbOut-2)
 */
  FixHornerEvaluator::FixHornerEvaluator(Target* target,
																				 int lsbIn_, int msbOut_, int lsbOut_,
																				 int degree_, vector<int> msbCoeff_, int lsbCoeff_, bool signedXandCoeffs_,
																				 bool finalRounding_, map<string, double> inputDelays)
	: Operator(target), degree(degree_), lsbIn(lsbIn_), msbOut(msbOut_), lsbOut(lsbOut_),
			msbCoeff(msbCoeff_), lsbCoeff(lsbCoeff_), signedXandCoeffs(signedXandCoeffs_),
			finalRounding(finalRounding_)
  {

	/* Generate unique name */
	{
	  std::ostringstream o;
	  o << "FixHornerEvaluator_" << getNewUId() << "_";
	  if(target->isPipelined())
				o << target->frequencyMHz() ;
	  else
				o << "comb";
	  uniqueName_ = o.str();
	}

	setCopyrightString("F. de Dinechin (2014)");
	srcFileName="FixHornerEvaluator";

		if(!signedXandCoeffs)
			REPORT(0,"signedXandCoeffs=false, this code has probably never been tested in this case. If it works, please remove this warning. If it doesn't, we deeply apologize and invite you to fix it.");

		// computing the coeff sizes
		for (int i=0; i<=degree; i++)
			coeffSize.push_back(msbCoeff[i]-lsbCoeff+1); // see FixConstant.hpp for the constant format

	// declaring inputs
		REPORT(0, "lsbIn=" << lsbIn);
		if(signedXandCoeffs){
			addInput("X"  , 0 - lsbIn+1);
			vhdl << tab << declareFixPoint("Xs", true, 0, lsbIn) << " <= signed(X);" << endl;
		}
		else{
			addInput("X"  , -lsbIn);
			vhdl << tab << declareFixPoint("Xs", false, -1, lsbIn) << " <= unsigned(X);" << endl;
		}
		for (int i=0; i<=degree; i++)
			addInput(join("A",i), coeffSize[i]);

	// declaring outputs
	addOutput("R", msbOut-lsbOut+1);
		setCriticalPath( getMaxInputDelays(inputDelays) + target->localWireDelay() );

		for(int i=0; i<=degree; i++) {
			vhdl << tab << declareFixPoint(join("As", i), signedXandCoeffs, msbCoeff[i], lsbCoeff)
					 << " <= " << (signedXandCoeffs?"signed":"unsigned") << "(" << join("A",i) << ");" <<endl;
		}


		// there are degree multiplications, with 2ulp errors: faithful mult, and truncation of x; the additions are exact. Hence,
		int g = intlog2(degree); // TODO: replace the +1 with a formula involving the  margin approx error.
		int lsbGlobal = lsbOut-2-g; // the shared internal lsb that will ensure epsilon_round < 2^(lsbOut-2)

		int msbSigma[1000], lsbSigma[1000], msbP[1000], lsbP[1000], lsbXTrunc[1000];

		// initialize the recurrence
		msbSigma[degree] = msbCoeff[degree];
		lsbSigma[degree] = lsbCoeff;

		for(int i=degree-1; i>=0; i--) {
			lsbXTrunc[i] = max(lsbIn, lsbGlobal-msbCoeff[i]);
			// P is the product of sigma_i+1 by xtrunc: sum of the MSBs+1, and sum the LSBs
			msbP[i] = msbSigma[i+1] + 0 + 1;
			lsbP[i] = lsbSigma[i+1] + lsbXTrunc[i];
			msbSigma[i] = max(msbP[i], msbCoeff[i]) + 1; // +1 for addition overflow, which is very rare
			lsbSigma[i] = lsbGlobal;
			REPORT(DEBUG, "i="<< i  << " lsbXtrunc=" << 	lsbXTrunc[i] << " msbP=" << msbP[i]<< " lsbP=" << lsbP[i]<< " msbSigma=" << 	msbSigma[i]<< " lsbSigma=" << 	lsbSigma[i]);
		}

		// Now generate the hardware
		vhdl << tab << declareFixPoint(join("Sigma", degree), true, msbSigma[degree], lsbSigma[degree])
				 << " <= " << join("As", degree)  << ";" << endl;

		for(int i=degree-1; i>=0; i--) {
			resizeFixPoint(join("XsTrunc", i), "Xs", 0, lsbXTrunc[i]);

			//  assemble faithful operators (either FixMultAdd, or truncated mult)
			if(target->plainVHDL()) {	// stupid pipelining here
				vhdl << tab << declareFixPoint(join("P", i), true, msbP[i],  lsbP[i])
						 <<  " <= "<< join("XsTrunc", i) <<" * Sigma" << i+1 << ";" << endl;

				// Align before addition
				resizeFixPoint(join("Ptrunc", i), join("P", i), msbSigma[i], lsbSigma[i]);
				resizeFixPoint(join("Aext", i), join("As", i), msbSigma[i], lsbSigma[i]);
				setCycle(getCurrentCycle() + target->plainMultDepth(1-lsbXTrunc[i], msbSigma[i]-lsbSigma[i]+1) );

				vhdl << tab << declareFixPoint(join("Sigma", i), true, msbSigma[i], lsbSigma[i])
						 << " <= " << join("Aext", i) << " + " << join("Ptrunc", i) << ";" << endl;
				nextCycle();
			}

			else { // using FixMultAdd
				//REPORT(DEBUG, "*** iteration " << i << );
				FixMultAdd::newComponentAndInstance(this,
																						join("Step",i),     // instance name
																						join("XsTrunc",i),  // x
																						join("Sigma", i+1), // y
																						join("As", i),       // a
																						join("Sigma", i),   // result
																						msbSigma[i], lsbSigma[i]
																						);
			}
			syncCycleFromSignal(join("Sigma", i));
			nextCycle();
		}
		if(finalRounding)
			resizeFixPoint("Ys", "Sigma0",  msbOut, lsbOut);

		vhdl << tab << "R <= " << "std_logic_vector(Ys);" << endl;
  }

	FixHornerEvaluator::~FixHornerEvaluator(){}

}
