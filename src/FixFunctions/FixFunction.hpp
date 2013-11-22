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
	
	/** The FixFunction objects holds a fixed-point function. 
			It provides an interface to Sollya services such as 
			parsing it,
			evaluating it at arbitrary precision,
			providing a polynomial approximation on an interval
	*/

	class FixFunction {
	public:

		FixFunction(string sollyaString);
		FixFunction(sollya_obj_t fS);

		virtual ~FixFunction();

		string getDescription() const;
		double eval(double x) const;
		void eval(mpfr_t r, mpfr_t x) const;
		sollya_obj_t getSollyaObj() const;

	private:

		string description;
		sollya_obj_t fS;

		string srcFileName; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		
		//		int verbose; /**< Flopoco-like verbosity level, from 0 to 4*/


	};

}
#endif // _FIXFUNCTION_HH_
