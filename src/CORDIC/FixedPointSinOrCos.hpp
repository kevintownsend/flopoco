#ifndef FixSinOrCos_HPP
#define FixSinOrCos_HPP

#include "../Operator.hpp"
#include "../utils.hpp"
#include "../IntMultiplier.hpp"
#include "../IntAdder.hpp"

#include "../FixedPointFunctions/FunctionEvaluator.hpp"

namespace flopoco{ 


	class FixedPointSinOrCos : public Operator {
	  
	  public:
		int w;
		gmp_randstate_t state;
	  
		// constructor, defined there with two parameters (default value 0 for each)
		FixedPointSinOrCos(Target* target, int w, map<string, double> inputDelays = emptyDelayMap);

		// destructor
		~FixedPointSinOrCos() {};


		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
		  in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		/* function used to create Standard testCase defined by the developper */
		void buildStandardTestCases(TestCaseList* tcl);


		/* function used to generate n random test, where n is an argument of
		  the function , this function is DEPRECATED, use the new function below*/
		TestCase* buildRandomTestCase(int n);

		// definition of some function for the operator
		std::string generateFixPointNumber(float x, int wI, int wF);
		std::string generateFixPointNumber(mpf_t x, int wI, int wF);
		
		std::string getParamName(std::string s, int number);
		
		mpz_class fp2fix(mpfr_t x, int wI, int wF);
	};

}

#endif
