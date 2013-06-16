#ifndef _FIXFUNCTION_HPP_
#define _FIXFUNCTION_HPP_

#include <string>
#include <iostream>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp"

using namespace std;

/* Stylistic convention here: all the sollya_obj_t have names that end with a capital S */
namespace flopoco{
	
	/** A fixed-point constant is an arbitrary precision integer scaled by a constant value*/
	class FixCoefficient {
	public:
		FixCoefficient(sollya_obj_t x);
		virtual ~FixCoefficient();

	private:
		mpz_class intSignificand;
		int exponent;
	};






	/** The FixFunction objects holds a fixed-point function on an interval. 
			It provides an interface to Sollya services such as 
         parsing it,
				 evaluating it at arbitrary precision,
				 providing a polynomial approximation
	*/

	class FixFunction {
	public:

		/** A serious constructor with arbitrary precision xmin and xmax  */
		FixFunction(string sollyaString, mpfr_t xmin, mpfr_t xmax);

		/** A constructor with a simpler interface for everyday use, like [0,1] */
		FixFunction(string sollyaString, double xmin = 0, double xmax = 1);

		virtual ~FixFunction();

		string getDescription() const;
		double eval(double x) const;
		void eval(mpfr_t r, mpfr_t x) const;
		sollya_obj_t getSollyaObj() const;


		/** Attempts to build a polynomial approximation. 
				Returns degree if success, -1 if failure, but sometimes also never returns*/
		int buildPolyApprox(double targetAccuracy); 
	private:
		mpfr_t xmin;   /**< lower bound on the interval */
		mpfr_t xmax;   /**< upper bound on the interval */
 		mpfr_t xmaxabs; /**< maximum absolute value of x*/

		string description;
		sollya_obj_t fS;
		sollya_obj_t rangeS;

		string srcFileName; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		
		int verbose; /**< Flopoco-like verbosity level, from 0 to 4*/

		int polyApproxDegree;               /**< the degree of the polynomial approximation, if any */
		sollya_obj_t polynomialS;           /**< the Sollya version of the polynomial, if any */
		vector<FixCoefficient*> polyApprox;  /**< the civilized version of the polynomial approximation, if any */
		double polyApproxError;             /**< the guaranteed approximation error (computed using Sollya's supnorm) */

		void finishConstruction(string sollyaString) ;

	};

}
#endif // _FIXFUNCTION_HH_
