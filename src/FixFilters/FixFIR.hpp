#ifndef FIXFIR_HPP
#define FIXFIR_HPP

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco{ 

	
	class FixFIR : public Operator {
	  
		public:
			/* Constructor  */
			FixFIR(Target* target, int p_, vector<string> coeff_, bool useBitheap = false, map<string, double> inputDelays = emptyDelayMap); 

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
	};

}

#endif
