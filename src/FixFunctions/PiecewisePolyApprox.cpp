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

	void PiecewisePolyApprox::build() {
		
		BasicPolyApprox *p;
		// First try alpha=0, because it is a bit special
		alpha=0; 
		p = new BasicPolyApprox(f, degree);
		if (p->approxErrorBound <= targetAccuracy){
			// Success at first try! Store it and return
			poly.push_back(p);
			approxErrorBound = p->approxErrorBound;
			return;
		}

		// Limit alpha to 24, because alpha will be the number of bits input to a table
		// it will take too long before that anyway
		for (alpha=1; alpha<24; alpha++) {
			int nbIntervals=1<<alpha;
			
		}

	}


} //namespace
