#ifndef FIXFIR_HPP
#define FIXFIR_HPP

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco{ 

	
	class FixFIR : public Operator {
	  
	public:
		/** normal constructor, building the FIR out of the coefficients */
		FixFIR(Target* target, int lsb_, vector<string> coeff_, map<string, double> inputDelays = emptyDelayMap); 
		
		/**empty constructor, to be called by subclasses */
		FixFIR(Target* target, int lsb_);

		/* Destructor */
		~FixFIR();

		/** The method that does the bulk of operator construction, isolated to enable sub-classes such as FixHalfSine etc */
		void buildVHDL();

		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
			 in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		int p;							/**< The precision (opposite of LSB weight) of inputs and outputs */ 
		int n;							/**< number of taps */
		vector<string> coeff;			/**< the coefficients as strings */

		int wO;

	private:
		mpfr_t mpcoeff[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsign[10000];			/**< the signs of the coefficients */
		mpz_class xHistory[10000]; // history of x used by emulate
		int currentIndex;
	};

}

#endif
