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
#include "IntMult//FixMultAdd.hpp"

using namespace std;



#define LARGE_PREC 1000 // 1000 bits should be enough for everybody

namespace flopoco{

  FixHornerEvaluator::FixHornerEvaluator(Target* target, 
																				 int wX_, 
																				 int degree_, vector<int> coeffMSB_, vector<int> coeffLSB_, bool signedCoeffs_, 
																				 bool finalRounding_, bool plainStupidVHDL_, map<string, double> inputDelays)
    : Operator(target), wX(wX_), degree(degree_), coeffLSB(coeffLSB_), coeffMSB(coeffMSB_), signedCoeffs(signedCoeffs_), finalRounding(finalRounding_), plainStupidVHDL(plainStupidVHDL_)
  { 
    
    /* Generate unique name */
    {
      std::ostringstream o;
      o << "FixHornerEvaluator_" << wX ;
      for(int i=0; i<=degree; i++)
				o << "_" << coeffMSB[i]<< "_" << coeffLSB[i];
      o << "_" << (signedCoeffs?1:0);
      if(target->isPipelined()) 
				o << target->frequencyMHz() ;
      else
				o << "comb";
      uniqueName_ = o.str();
    }
    
    setCopyrightString("F. de Dinechin (2014)");
    srcFileName="FixHornerEvaluator";

		// computing the coeff sizes
		for (int i=0; i<coeffMSB.size(); i++)
			coeffSize.push_back(coeffMSB[i]-coeffLSB[i]+1); // see FixConstant.hpp for the constant format

    // declaring inputs
    addInput("X", wX);
		for (int i=0; i<coeffMSB.size(); i++)
			addInput(join("A",i), coeffSize[i]);
			
    // declaring outputs
    addOutput("R", coeffMSB[0]-coeffLSB[0]+1);


		/*
     Analysis of the polynomial coefficients to determine the intermediate sizes.
     This class knows the abs max of the coefficients to do a worst-case analysis.
		 This may be pessimistic, e.g. log has alternate coeffs and doesn't grow as far as worst case suggests.
		 Managing that properly is a TODO: it impacts only efficiency, not numerical quality.
		 Maybe it can be done with a bit more Sollya.
 
		 This is a simplified version of the computation in the ASAP 2010 paper, simplified because x is in [0,1) 

		 WRT the ASAP paper,	
		 1/we want to use FixMultAdd operations: 
		 \sigma_0  =  a_d
		 \sigma_j  =  a_{d-j} + x^T_i * \sigma_{j-1}     \forall j \in \{ 1...d\}
		 p(y)      =  \sigma_d
		 
		 2/ we integrate the multiplication guard bits (g^\pi_j in the paper) into the coeffs a_j themselves
		 so we truncate the product \pi_j to the LSB of a_{d-j} so that the addition is aligned
		 Philosophically, this should provide a better balance of approx versus rounding error.
		 But main advantage is that it is simpler to maintain.

		 We still should consider the DSP granularity to truncate y to the smallest DSP-friendly size larger than |sigma_{j-1}|
		 TODO

		 So all that remains is to compute the parameters of the FixMultAdd.
		 We need 
		 * size of \sigma_d: this is the size of a_{d-j}
		 * size of y truncated : see above
		 * weight of the MSB of the product: since y \in [0,1), this is the MSB of \sigma_{j-1}, i.e. the MSB of a_{d-j+1}
		 ** true in case of unsigned coeffs
		 ** in case of signed coeffs, the max abs value is 2^(MSB-1).
		 **   however, since we have zero-extended x in this case there is some thinking TODO 
		 * size of the truncated result: 
		 * weight of the LSB of a_{d-j}
		 * weight of the MSB of a_{d-j}
		 */

		vhdl << tab << tab << declare("h0", coeffSize[degree]) << " <= " << join("A", degree)  << endl;

		int wInX = wX+(signedCoeffs?1:0);
		vhdl << tab << tab << declare("Xse", wInX) << " <= " << (signedCoeffs?"\"0\"&":"") << "X;" << (signedCoeffs?"-- sign extension before signed mult":"") << endl;


		if(signedCoeffs) {
			// extend x in [0,1] with a sign bit
			int wInX=wX+1; //  
		}
		else{ 
			wInX=wX;
		}
		// anyway this extension is at the MSB. The LSB is still -wX.
		int lsbIn = -wX;

		// initialize the recurrence
		int sigmaMSB=coeffMSB[degree];
		int sigmaLSB=coeffLSB[degree];
    // Now assemble faithful FixMultAdd operators

		for(int i=degree-1; i>=0; i--) {

			int xTruncLSB = max(lsbIn, sigmaLSB-sigmaMSB);

			int pMSB=sigmaMSB+0 + 1;
			sigmaMSB = max(pMSB-1, coeffMSB[i]) +1; // +1 to absorb addition overflow
			sigmaLSB = coeffLSB[i];

			resizeFixPoint(join("XsTrunc", i), "Xs", 0, xTruncLSB);			

			if(plainStupidVHDL) {
				vhdl << tab << declareFixPoint(join("P", i), true, pMSB,  sigmaLSB  + xTruncLSB /*LSB*/) 
						 <<  " <= "<< join("XsTrunc", i) <<" * Sigma" << i+1 << ";" << endl;
				// However the bit of weight pMSB is a 0. We want to keep the bits from  pMSB-1
				resizeFixPoint(join("Ptrunc", i), join("P", i), sigmaMSB, sigmaLSB);
				resizeFixPoint(join("Aext", i), join("A", i), sigmaMSB, sigmaLSB);
				
				vhdl << tab << declareFixPoint(join("Sigma", i), true, sigmaMSB, sigmaLSB)   << " <= " << join("Aext", i) << " + " << join("Ptrunc", i) << ";" << endl;
			}

			else { // using FixMultAdd
				REPORT(0, " i=" << i);
				FixMultAdd::newComponentAndInstance(this,
																						join("Step",i),     // instance name
																						join("XsTrunc",i),  // x
																						join("Sigma", i+1), // y
																						join("A", i),       // a
																						join("Sigma", i),   // result 
																						true, sigmaMSB, sigmaLSB  // signed, outMSB, outLSB
																						);
			}
		}
    for (int i=1; i<=degree; i++) {
		
		}

  }


}
