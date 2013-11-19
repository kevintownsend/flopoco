#ifndef _PIECEWISEPOLYAPPROX_HPP_
#define _PIECEWISEPOLYAPPROX_HPP_

#include <string>
#include <iostream>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp" // mostly for reporting
#include "FixFunction.hpp"
#include "FixConstant.hpp"

using namespace std;

/* Stylistic convention here: all the sollya_obj_t have names that end with a capital S */
namespace flopoco{
	
	/** The PiecewisePolyApprox object builds and maintains a machine-efficient polynomial approximation to a fixed-point function over [0,1]
			Fixed point, hence only absolute errors/accuracy targets.
			Eventually there should be two constructors: 
			one that inputs target accuracy and computes degree and approximation error,
			one that inputs degree and computes approximation error.

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

	class PiecewisePolyApprox {
	public:

		/** A minimal constructor inputting target accuracy
				@param addGuardBitsToConstant: 
				if >=0, add this number of bits to the LSB of the constant
				if -1, add the bits needed for a Horner evaluation based on faithful (truncated) multipliers

		 */
		PiecewisePolyApprox(FixFunction* f, double targetAccuracy, int addGuardBitsToConstant=0);

		/** A minimal constructor that parses a sollya string, inputting target accuracy
				@param addGuardBitsToConstant: 
				if >=0, add this number of bits to the LSB of the constant
				if -1, add the bits needed for a Horner evaluation based on faithful (truncated) multipliers

		 */
		PiecewisePolyApprox(string sollyaString, double targetAccuracy, int addGuardBitsToConstant=0);


		/** A minimal constructor inputting degree
		 */
		PiecewisePolyApprox(FixFunction *f, int degree);


		virtual ~PiecewisePolyApprox();

		int getDegree() const;

		sollya_obj_t getSollyaPolynomial() const;

	private:
		FixFunction *f;                   /**< The function to be approximated */
		sollya_obj_t polynomialS;         /**< The polynomial approximating it */
		int degree;                       /**< degree of the polynomial approximation */
		vector<FixConstant*> coeff;  /**< polynomial coefficients in a hardware-ready form */
		int LSB;                          /**< weight of the LSB of the polynomial approximation */
		int constLSB;                     /**< weight of the LSB of the constant coeff, may be smaller than LSB for free */
		double approxError;               /**< guaranteed upper bound on the approx error */


		string srcFileName; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		bool needToFreeF;   /**< in an ideal world, this should not exist */
		void buildApproxFromTargetAccuracy(double targetAccuracy, int addGuardBitsToConstant); /** < constructor code factored out */
		void buildFixFormatVector(); /** < constructor code, factored out */


	};

}
#endif // _POLYAPPROX_HH_


// Garbage below
