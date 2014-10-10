/*

  A class that manages polynomial approximation for FloPoCo (and possibly later for metalibm). 

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL

  All rights reserved.

*/


/* 
	 The function is assumed to have inputs in [0,1]
*/
#include "PiecewisePolyApprox.hpp"
#include <sstream>
#include <limits.h>
#include <float.h>

namespace flopoco{

	PiecewisePolyApprox::PiecewisePolyApprox(FixFunction *f_, double targetAccuracy_, int degree_): 
		f(f_), targetAccuracy(targetAccuracy_), degree(degree_)
	{
		needToFreeF = false;
		srcFileName="PiecewisePolyApprox"; // should be somehow static but this is too much to ask me
		build();
	}


	PiecewisePolyApprox::PiecewisePolyApprox(string sollyaString_, double targetAccuracy_, int degree_):
	targetAccuracy(targetAccuracy_), degree(degree_)
	{
		//  parsing delegated to FixFunction
		f = new FixFunction(sollyaString_);
		needToFreeF = true;
		srcFileName="PiecewisePolyApprox"; // should be somehow static but this is too much to ask me
		build();
	}



	PiecewisePolyApprox::~PiecewisePolyApprox()
	{
		if(needToFreeF)	free(f);
	}


	/** a local function to build g_i(x) = f(2^-alpha*x + i*2^-alpha) */
	sollya_obj_t buildSubIntervalFunction(sollya_obj_t fS, int alpha, int i){
		stringstream s;
		s << "(1b-" << alpha << ")*x + ("<< i << "b-" << alpha << ")";
		string ss = s.str(); // do not use c_str directly on the stringstream, it is too transient (?) 
		sollya_obj_t newxS = sollya_lib_parse_string(ss.c_str());
		sollya_obj_t giS = sollya_lib_substitute(fS,newxS);
		sollya_lib_clear_obj(newxS);
		return giS;
	}

	// split into smaller and smaller intervals until the function can be approximated by a polynomial of degree degree.
	void PiecewisePolyApprox::build() {
		
		sollya_obj_t fS = f->getSollyaObj(); // no need to free this one
		int nbIntervals;

		// Limit alpha to 24, because alpha will be the number of bits input to a table
		// it will take too long before that anyway
		bool alphaOK;
		for (alpha=0; alpha<24; alpha++) {
			nbIntervals=1<<alpha;
			alphaOK=true;
			REPORT(DETAILED, " Testing alpha=" << alpha );
			for (int i=0; i<nbIntervals; i++) {
				// The worst case is typically on the left (i==0) or on the right (i==nbIntervals-1).
				// To test these two first, we do this small rotation of i
				int ii=(i+nbIntervals-1) & ((1<<alpha)-1);

				// First build g_i(x) = f(2^-alpha*x + i*2^-alpha)
				sollya_obj_t giS = buildSubIntervalFunction(fS, alpha, ii);

				if(DEBUG <= verbose)
					sollya_lib_printf("> PiecewisePolyApprox: alpha=%d, ii=%d, testing  %b \n", alpha, ii, giS);
				// Now what degree do we need to approximate gi?
				int degreeInf, degreeSup;
				BasicPolyApprox::guessDegree(giS, targetAccuracy, &degreeInf, &degreeSup);
				// REPORT(DEBUG, " guessDegree returned (" << degreeInf <<  ", " << degreeSup<<")" ); // no need to report, it is done by guessDegree()
				sollya_lib_clear_obj(giS);
				// For now we only consider degreeSup. Is this a TODO?
				if(degreeSup>degree) {
					REPORT(DEBUG, "   alpha=" << alpha << " failed." );
					alphaOK=false;
					break;
				}
			} // end for loop on i
			
			// Did we success?
			if (alphaOK)
				break;
		} // end for loop on alpha
		if (alphaOK)
			REPORT(DETAILED, " Found alpha=" << alpha << " OK");

		// Still have to do a while loop because we can't trust guessdegree, damn 
		bool success=false;
		while(!success) {
			// Now fill the vector of polynomials, computing the coefficient parameters along.
			//		LSB=INT_MAX; // very large
			// MSB=INT_MIN; // very small
			approxErrorBound = 0.0;
			BasicPolyApprox *p;

			REPORT(DETAILED, " Now computing the actual polynomials ");
			// initialize the vector of MSB weights
			for (int j=0; j<=degree; j++) {
				MSB.push_back(INT_MIN);
			}

			// Compute the LSB of each coefficient.
			// It should be at least floor(log2(targetAccuracy)); 
			// However we can get a bit extra accuracy/slack because the evaluation will use log2(degree) guard bits.
			// Adding these to the constants is almost for free: let's do it.  
			LSB = floor(log2(targetAccuracy/degree));
			REPORT(DEBUG, "To obtain target accuracy " << targetAccuracy << " with a degree-"<<degree <<" polynomial, we compute coefficients accurate to " << targetAccuracy/degree
						 << " (LSB="<<LSB<<")"); 
		
			for (int i=0; i<nbIntervals; i++) {
				REPORT(DETAILED, " ----------Interval " << i << "-------------");
				// Recompute the substitution. No big deal.
				sollya_obj_t giS = buildSubIntervalFunction(fS, alpha, i);

				p = new BasicPolyApprox(giS, degree, LSB);
				poly.push_back(p);
				if (approxErrorBound < p->approxErrorBound){ 
					REPORT(DEBUG, "   new approxErrorBound=" << p->approxErrorBound );
					approxErrorBound = p->approxErrorBound;
				}

				// Now compute the englobing MSB and LSB for each coefficient	
				for (int j=0; j<=degree; j++) {
					// if the coeff is zero, we can set its MSB to anything, so we exclude this case
					if (  (!p->coeff[j]->isZero())  &&  (p->coeff[j]->MSB > MSB[j])  )
						MSB[j] = p->coeff[j]->MSB;
				}
		
			} // end for loop on i


			if (approxErrorBound < targetAccuracy) {
				REPORT(INFO, " *** Success! Final approxErrorBound=" << approxErrorBound << "  is smaller than target accuracy: " << targetAccuracy  );
				success=true;
			}
			else {
				REPORT(INFO, " guessDegree mislead us, measured approx error:" << approxErrorBound << " is larger than target accuracy: " << targetAccuracy << ". Now increasing alpha and starting over, thank you for your patience");
				alpha++;
				nbIntervals=1<<alpha;
			}
		} // end while(!success)


		// Now we need to resize all the coefficients of degree i to the largest one
		for (int i=0; i<nbIntervals; i++) {
			for (int j=0; j<=degree; j++) {
				// REPORT(DEBUG "Resizing MSB of coeff " << j << " of poly " << i << " : from " << poly[i] -> coeff[j] -> MSB << " to " <<  MSB[j]);  
				poly[i] -> coeff[j] -> changeMSB(MSB[j]);
				// REPORT(DEBUG, "   Now  " << poly[i] -> coeff[j] -> MSB);  
			}
		}
		// TODO? In the previous loop we could also check if one of the coeffs is always positive or negative, and optimize generated code accordingly 

		// A bit of reporting
		REPORT(INFO,"Final report: ");
		REPORT(INFO,"  Degree=" << degree	<< "      maxApproxErrorBound=" << approxErrorBound);
 		int totalOutputSize=0;
		for (int j=0; j<=degree; j++) {
			int size = MSB[j]-LSB +1;
			totalOutputSize += size ;
			REPORT(INFO,"      MSB["<<j<<"] = " << MSB[j] << "  size=" << size);
		}
		REPORT(INFO, "  Total size of the table is " << nbIntervals << " x " << totalOutputSize << " bits");

	}
	
	mpz_class PiecewisePolyApprox::getCoeff(int i, int d){
		BasicPolyApprox* p = poly[i];
		FixConstant* c = p->coeff[d];
		return c->getBitVectorAsMPZ();
	} 

} //namespace
