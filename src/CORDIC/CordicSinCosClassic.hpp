#ifndef CordicSinCosClassicCLASSIC_HPP
#define CordicSinCosClassicCLASSIC_HPP

#include "../Operator.hpp"
#include "../utils.hpp"
#include "../IntMultiplier.hpp"
#include "../IntAdder.hpp"
#include "../ConstMult/FixRealKCM.hpp"

namespace flopoco{ 


	class CordicSinCosClassic : public Operator {
	  
	  public:
		int w;
		gmp_randstate_t state;
	  
		// constructor, defined there with two parameters (default value 0 for each)
		CordicSinCosClassic(Target* target, int w, map<string, double> inputDelays = emptyDelayMap);

		// destructor
		~CordicSinCosClassic() {};


		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
		  in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		/* function used to create Standard testCase defined by the developper */
		void buildStandardTestCases(TestCaseList* tcl);

		// definition of some function for the operator
		std::string generateFixPointNumber(float x, int wI, int wF);
		std::string generateFixPointNumber(mpf_t x, int wI, int wF);
		
		mpz_class fp2fix(mpfr_t x, int wI, int wF);
	};

}

#endif

