#ifndef _POLYAPPROX_HPP_
#define _POLYAPPROX_HPP_

#include <string>
#include <iostream>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FixFunction.hpp"

using namespace std;

/* Stylistic convention here: all the sollya_obj_t have names that end with a capital S */
namespace flopoco{
	
	/** The PolyApprox object builds and maintains a machine-efficient polynomial approximation to a fixed-point function over [0,1]
			Fixed point, hence only absolute errors/accuracy targets.
			Eventually there should be two constructors: 
			one that inputs target accuracy and computes degree and achieved accuracy,
			one that inputs degree and computes achieved accuracy.

			The two are needed in the typical case of a domain split: we will first call the first for each subinterval,
			then take the max of the degrees, and call the second with that.

			We then call guessDegree.
			target_accuracy implies the target LSB of the constant part of the polynomial.
			if  addGuardBitsToConstante, we add g=ceil(log2(degree+1)) bits to the LSB of the constant: 
			this provides a bit of freedom to fpminimax, for free in terms of evaluation.
			And we call fpminimax with that, and check what it provided.
			Since x is in [0,1], the MSB of coefficient i is then exactly the MSB of a_i.x^i
			For the same reason, all the coefficients should have the same target LSB
			
			To start with, this class is for "expert" users: 
			No domain splitting or other range reduction here: these should be managed in different classes
			No good handling of functions with zero coefficients for now. 
			- a function with zero constant should be transformed into a "proper" one outside this class. 
			   Example: sin, log(1+x)
		  - a function with odd or even Taylor should be transformed as per the Muller book.
			
			To implement a generic approximator we will need to lift these restrictions, but it is unclear that it needs to be in this class.
	*/

	class PolyApprox {
	public:

		/** A minimal constructor  
				@param addGuardBitsToConstant: 
				if >=0, add this number of bits to the LSB of the constant
				if -1, add the bits needed for a Horner evaluation based on faithful (truncated) multipliers

		 */
		PolyApprox(string sollyaString, double targetAccuracy, int addGuardBitsToConstant=0);

		virtual ~PolyApprox();

		string getDescription() const;
		int getDegree() const;

		// TODO do we need what follows? Is PolyApprox a FixFunction? Should it input a FixFunction, too? 
		double eval(double x) const;
		void eval(mpfr_t r, mpfr_t x) const;
		sollya_obj_t getSollyaPolynomial() const;

	private:

		string description;
		sollya_obj_t fS;
		sollya_obj_t polynomialS;

		string srcFileName; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		
		int degree;  /**< degree of the polynomial approximation */
		int lsb;     /**< lsb of the polynomial approximation, computed by adding enough bits to the lsb corresponding to targetAccuracy*/
		double approxError; /**< guaranteed upper bound on the approx error */

		vector<int> weightLSB;
		vector<int> weightMSB;
	};

}
#endif // _POLYAPPROX_HH_


// Garbage below


#if 0
	/** A fixed-point constant is an arbitrary precision integer scaled by a constant value*/
	class FixCoefficient {
	public:
		FixCoefficient(sollya_obj_t x);
		virtual ~FixCoefficient();

	private:
		mpz_class intSignificand;
		int exponent;
	};

		/** Attempts to build a polynomial approximation on an interval. 
				Returns degree if success, -1 if failure, but sometimes also never returns*/
		int buildPolyApprox(double targetAccuracy, mpfr_t xmin, mpfr_t xmax); 
		int polyApproxDegree;               /**< the degree of the polynomial approximation, if any */
		sollya_obj_t polynomialS;           /**< the Sollya version of the polynomial, if any */
		vector<FixCoefficient*> polyApprox;  /**< the civilized version of the polynomial approximation, if any */
		double polyApproxError;             /**< the guaranteed approximation error (computed using Sollya's supnorm) */
		sollya_obj_t rangeS;

		void finishConstruction(string sollyaString) ;
#endif
