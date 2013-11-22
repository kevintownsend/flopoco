/*

  A class that manages polynomial approximation for FloPoCo (and possibly later for metalibm). 

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL

  All rights reserved.

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


		// Now fill the vector of polynomials, computing the coefficient parameters along.
		//		LSB=INT_MAX; // very large
		// MSB=INT_MIN; // very small
		approxErrorBound = 0.0;
		BasicPolyApprox *p;


		for (int i=0; i<nbIntervals; i++) {
			// Recompute the substitution. No big deal.
			sollya_obj_t giS = buildSubIntervalFunction(fS, alpha, i);

			// TODO Bug here, BasicPolyApprox recomputes guessDegree, takes the min and tries that,whereas we want the sup and does all sorts of optimizations, 
			// define a simpler constructor
			p = new BasicPolyApprox(giS, targetAccuracy);
			poly.push_back(p);
			if (approxErrorBound < p->approxErrorBound){ 
				REPORT(DEBUG, "   new approxErrorBound=" << p->approxErrorBound );
				approxErrorBound = p->approxErrorBound;
			}
			
		} // end for loop on i

#if 0
		// First try alpha=0, because it is a bit special
		alpha=0; 
		p = new BasicPolyApprox(f, targetAccuracy);
		if (p->approxErrorBound <= targetAccuracy){
			// Success at first try! Store it and return
			return;
		}
#endif


	}
	

} //namespace
