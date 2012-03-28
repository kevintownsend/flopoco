#ifndef CordicSinCosClassicCLASSIC_HPP
#define CordicSinCosClassicCLASSIC_HPP

#include "../Operator.hpp"
#include "../utils.hpp"
#include "../IntAdder.hpp"

namespace flopoco{ 


	class CordicSinCosClassic : public Operator {
	  
	  public:
		int w;
		gmp_randstate_t state;

		int guard;

		// constructor, defined there with two parameters (default value 0 for each)
		CordicSinCosClassic(Target* target, int w, map<string, double> inputDelays = emptyDelayMap);

		// destructor
		~CordicSinCosClassic() {};



		// definition of some function for the operator
		std::string generateFixPointNumber(float x, int wI, int wF);
		std::string generateFixPointNumber(mpf_t x, int wI, int wF);
		
		mpz_class fp2fix(mpfr_t x, int wI, int wF);
	};

}

#endif

