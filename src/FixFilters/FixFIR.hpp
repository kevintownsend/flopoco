#ifndef FIXFIR_HPP
#define FIXFIR_HPP

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco{ 

	
	class FixFIR : public Operator {
	  
		public:
			/* Constructor ; p must be greater than 6 and you must use bitheap in case of negative coefficient*/
			FixFIR(Target* target, int p_, vector<string> coeff_, bool useBitheap = true, map<string, double> inputDelays = emptyDelayMap); 

			/* Destructor */
			~FixFIR();

			// Below all the functions needed to test the operator
			/* the emulate function is used to simulate in software the operator
			  in order to compare this result with those outputed by the vhdl opertator */
			void emulate(TestCase * tc);

			/* function used to create Standard testCase defined by the developper */
			void buildStandardTestCases(TestCaseList* tcl);

	  	private:
			int p;							/**< The precision of inputs and outputs */ 
			int n;							/**< number of taps */
			vector<string> coeff;			/**< the coefficients as strings */

			bool useBitheap;

			int wO;

			mpfr_t mpcoeff[10000];			/**< the absolute values of the coefficients as MPFR numbers */
			bool coeffsign[10000];			/**< the signs of the coefficients */
	};

}

#endif
