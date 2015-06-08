/*

  A class that manages fixed-point  polynomial approximation for FloPoCo (and possibly later for metalibm).

	At the moment, only on [0,1], no possibility to have it on [-1,1].
	Rationale: no point really as soon as we are not in the degenerate case alpha=0.

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
		degree(degree_), f(f_), targetAccuracy(targetAccuracy_)
	{
		needToFreeF = false;
		srcFileName="PiecewisePolyApprox"; // should be somehow static but this is too much to ask me
		build();
	}


	PiecewisePolyApprox::PiecewisePolyApprox(string sollyaString_, double targetAccuracy_, int degree_):
		degree(degree_), targetAccuracy(targetAccuracy_)
	{
		//  parsing delegated to FixFunction
		f = new FixFunction(sollyaString_, false /* on [0,1]*/);
		needToFreeF = true;
		srcFileName="PiecewisePolyApprox"; // should be somehow static but this is too much to ask me
		build();
	}



	PiecewisePolyApprox::~PiecewisePolyApprox()
	{
		if(needToFreeF)	free(f);
	}



#if 0
	/** a local function to build g_i(x) = f(2^(-alpha)*x + i*2^-alpha) defined on [0,1] */

	sollya_obj_t buildSubIntervalFunction(sollya_obj_t fS, int alpha, int i){
		stringstream s;
		s << "(1b-" << alpha << ")*x + ("<< i << "b-" << alpha << ")";
		string ss = s.str(); // do not use c_str directly on the stringstream, it is too transient (?)
		sollya_obj_t newxS = sollya_lib_parse_string(ss.c_str());
		sollya_obj_t giS = sollya_lib_substitute(fS,newxS);
		sollya_lib_clear_obj(newxS);
		return giS;
	}
#else
	/** a local function to build g_i(x) = f(2^(-alpha-1)*x + i*2^-alpha + 2^(-alpha-1)) defined on [-1,1] */
	sollya_obj_t buildSubIntervalFunction(sollya_obj_t fS, int alpha, int i){
		stringstream s;
		s << "(1b-" << alpha+1 << ")*x + ("<< i << "b-" << alpha << " + 1b-" << alpha+1 << ")";
		string ss = s.str(); // do not use c_str directly on the stringstream, it is too transient (?)
		sollya_obj_t newxS = sollya_lib_parse_string(ss.c_str());
		sollya_obj_t giS = sollya_lib_substitute(fS,newxS);
		sollya_lib_clear_obj(newxS);
		return giS;
	}
#endif




	// split into smaller and smaller intervals until the function can be approximated by a polynomial of degree degree.
	void PiecewisePolyApprox::build() {

		ostringstream cacheFileName;
		cacheFileName << "PiecewisePoly_"<<vhdlize(f->description) << "_" << degree << "_" << targetAccuracy << ".cache";

		// cerr << cacheFileName.str() <<endl<<endl;
		// Test existence of cache file
		fstream file;
		file.open(cacheFileName.str().c_str(), ios::in);

		// check for bogus .cache file
		if(file.is_open() && file.peek() == std::ifstream::traits_type::eof())
			{
				file.close();
				std::remove(cacheFileName.str().c_str());
				file.open(cacheFileName.str().c_str(), ios::in); // of course this one will fail
			}

		if(!file.is_open()){
			//********************** Do the work, then write the cache *********************
			sollya_obj_t fS = f->fS; // no need to free this one
			sollya_obj_t rangeS;

#if 0
			rangeS  = sollya_lib_parse_string("[0;1]");
#else
			rangeS  = sollya_lib_parse_string("[-1;1]");
#endif
			// TODO test with [-1,1] which is the whole point of current refactoring.
			// There needs to be a bit of logic here because rangeS should be [-1,1] by default to exploit signed arith,
			// except in the case when alpha=0 because then rangeS should be f->rangeS (and should not be freed)

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
					BasicPolyApprox::guessDegree(giS, rangeS, targetAccuracy, &degreeInf, &degreeSup);
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

			// Compute the LSB of each coefficient. Minimum value is:
			LSB = floor(log2(targetAccuracy));
			REPORT(DEBUG, "To obtain target accuracy " << targetAccuracy << " with a degree-"<<degree <<" polynomial, we compute coefficients accurate to " << targetAccuracy/degree
						 << " (LSB="<<LSB<<")");
			// It is pretty sure that adding intlog2(degree) bits is enough for FPMinimax.
			int lsbAttemptsMax = intlog2(degree);
			int lsbAttempts=0; // a counter of attempts to move the LSB down, caped by lsbAttemptsMax


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


				for (int i=0; i<nbIntervals; i++) {
					REPORT(DETAILED, " ... computing polynomial approx for interval " << i << " / "<< nbIntervals);
					// Recompute the substitution. No big deal.
					sollya_obj_t giS = buildSubIntervalFunction(fS, alpha, i);

					p = new BasicPolyApprox(giS, degree, LSB, true);
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
						REPORT(INFO, "So far measured approx error:" << approxErrorBound << " is larger than target accuracy: " << targetAccuracy << ". Thank you for your patience");
					//empty poly
					for (int i=0; i<nbIntervals; i++) {
						free(poly.back());
						poly.pop_back();
					}
					if(lsbAttempts<=lsbAttemptsMax) {
						lsbAttempts++;
						LSB--;
					}
					else {
						REPORT(INFO, "guessDegree mislead us, increasing alpha and starting over");
						LSB+=lsbAttempts;
						lsbAttempts=0;
						alpha++;
						nbIntervals=1<<alpha;
					}
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

			// Write the cache file
			REPORT(INFO, "Writing to cache file: " << cacheFileName.str());
			file.open(cacheFileName.str().c_str(), ios::out);
			file << "Polynomial data cache for " << cacheFileName.str() << endl;
			file << "Erasing this file is harmless, but do not try to edit it." <<endl;
			file << degree <<endl;
			file << alpha <<endl;
			file << LSB <<endl;
			for (int j=0; j<=degree; j++) {
				file << MSB[j] << endl;
			}
			file << approxErrorBound << endl;
			// now the coefficients themselves
			for (int i=0; i<(1<<alpha); i++) {
				for (int j=0; j<=degree; j++) {
					file <<  poly[i] -> coeff[j] -> getBitVectorAsMPZ() << endl;
				}
			}
			sollya_lib_clear_obj(rangeS);
		}

		else {
			REPORT(INFO, "Polynomial data cache found: " << cacheFileName.str());
			//********************** Just read the cache *********************
			string line;
			getline(file, line); // ignore the first line which is a comment
			getline(file, line); // ignore the second line which is a comment
			file >> degree;
			file >> alpha;
			file >> LSB;

			for (int j=0; j<=degree; j++) {
				int msb;
				file >> msb;
				MSB.push_back(msb);
			}
			file >> approxErrorBound;

			for (int i=0; i<(1<<alpha); i++) {
				vector<mpz_class> coeff;
				for (int j=0; j<=degree; j++) {
					mpz_class c;
					file >> c;
					coeff.push_back(c);
				}
				BasicPolyApprox* p = new BasicPolyApprox(degree,MSB,LSB,coeff);
				poly.push_back(p);
			}
		} // end if cache
	}




	mpz_class PiecewisePolyApprox::getCoeff(int i, int d){
		BasicPolyApprox* p = poly[i];
		FixConstant* c = p->coeff[d];
		return c->getBitVectorAsMPZ();
	}

} //namespace
