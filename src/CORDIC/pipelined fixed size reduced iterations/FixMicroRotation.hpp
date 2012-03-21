#ifndef FIXMICROROTATION_HPP
#define FIXMICROROTATION_HPP

#include "../Operator.hpp"
#include "../IntAdder.hpp"
#include "../utils.hpp"

namespace flopoco{
 

	class FixMicroRotation : public Operator {
	  
	  public:
		int wIx, wFx, wIy, wFy, wIz, wFz, stage, xyIncrement;
		bool hasLargerZout;
	  
		// constructor, defined there with two parameters (default value 0 for each)
		FixMicroRotation(Target* target, int wIx_, int wFx_, int wIy_, int wFy_, int wIz_, int wFz_, int stage_, int xyIncrement_, map<string, double> inputDelays = emptyDelayMap);

		// destructor
		~FixMicroRotation() {};


		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
		  in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		/* function used to create Standard testCase defined by the developper */
		void buildStandardTestCases(TestCaseList* tcl);

		/* function used to generate n random test, where n is an argument of
		  the function , this function is DEPRECATED, use the new function below*/
		//void buildRandomTestCases(TestCaseList* tcl, int n);

		TestCase* buildRandomTestCase(int i);
		
		// definition of some function for the operator
		mpz_class fp2fix(mpfr_t x, int wI, int wF);
	};

}

#endif
