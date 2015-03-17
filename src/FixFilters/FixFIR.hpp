#ifndef FIXFIR_HPP
#define FIXFIR_HPP

#include "Operator.hpp"
#include "FixSOPC.hpp"
#include "utils.hpp"

namespace flopoco{ 

	
	class FixFIR : public Operator {
	  
	public:
		/** normal constructor, building the FIR out of the coefficients */
		FixFIR(Target* target, int lsbInOut, vector<string> coeff, map<string, double> inputDelays = emptyDelayMap); 
		
		/**empty constructor, to be called by subclasses */
		FixFIR(Target* target, int lsbInOut);

		/* Destructor */
		~FixFIR();

		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
			 in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

	private:
	protected:
		/** The method that does the bulk of operator construction, isolated to enable sub-classes such as FixHalfSine etc */
		void buildVHDL();

		int n;							/**< number of taps */
		vector<string> coeff;			  /**< the coefficients as strings */
		int lsbInOut;
		mpz_class xHistory[10000]; // history of x used by emulate
		int currentIndex;
		FixSOPC *fixSOPC; /**< most of the work done here */
	};

}

#endif
