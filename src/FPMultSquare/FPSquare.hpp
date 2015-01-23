#ifndef FPSQUARER_HPP
#define FPSQUARER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "TestState.hpp"

#include "Operator.hpp"
#include "IntMult//IntMultiplier.hpp"

namespace flopoco{

	/** The FPSquare class */
	class FPSquare : public Operator
	{
	public:
	
		/**
		 * The FPMutliplier constructor
		 * @param[in]		target		the target device
		 * @param[in]		wEX			the the with of the exponent for the f-p number X
		 * @param[in]		wFX			the the with of the fraction for the f-p number X
		 * @param[in]		wFR			the the with of the fraction for the multiplication result
		 **/
		FPSquare(Target* target, int wEX, int wFX, int wFR);

		/**
		 * FPSquare destructor
		 */
		~FPSquare();

		/**
		 * Emulate the operator using MPFR.
		 * @param tc a TestCase partially filled with input values 
		 */
		void emulate(TestCase * tc);

		// Static method used for tests on the Operator
		static void nextTest ( TestState * ts );

	protected:
	
		int  wE_;                  /**< The width of the exponent for the input X */
		int  wFX_;                  /**< The width of the fraction for the input X */
		int  wFR_;                  /**< The width of the fraction for the output R */


	private:
 
	};

}
#endif