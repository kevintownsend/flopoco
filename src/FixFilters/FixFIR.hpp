#ifndef FIXFIR_HPP
#define FIXFIR_HPP

#include "Operator.hpp"
#include "FixSOPC.hpp"
#include "utils.hpp"

namespace flopoco{


	class FixFIR : public Operator {

	public:
		/**
		 * @brief normal		constructor, building the FIR out of the coefficients.
		 * @param coeff			The coefficients are considered as real numbers, provided as string expresssions such as 0.1564565756768 or sin(3*pi/8)
		 * Input is assumed to be in [-1,1], with msb (sign bit) at position 0 and lsb at postion lsbInOut.
		 * @param lsbInOut		Output has lsb at position lsbInOut.
		 * @param rescale
		 *						If rescale=false, the msb of the output is computed so as to avoid overflow.
		 *						If rescale=true, all the coefficients are rescaled by 1/sum(|coeffs|).
		 * This way the output is also in [-1,1], output size is equal to input size, and the output signal makes full use of the output range.
		*/
		FixFIR(Target* target, int lsbInOut, vector<string> coeff, bool rescale=false, map<string, double> inputDelays = emptyDelayMap);

		/** @brief empty constructor, to be called by subclasses */
		FixFIR(Target* target, int lsbInOut, bool rescale=false);

		/** @brief Destructor */
		~FixFIR();

		// Below all the functions needed to test the operator
		/**
		 * @brief the emulate function is used to simulate in software the operator
		 * in order to compare this result with those outputed by the vhdl opertator
		 */
		void emulate(TestCase * tc);

		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(Target *target , vector<string> &args);
		static void registerFactory();

	private:
	protected:
		/** @brief The method that does the bulk of operator construction, isolated to enable sub-classes such as FixHalfSine etc */
		void buildVHDL();

		int n;							/**< number of taps */
		int lsbInOut;
		vector<string> coeff;			  /**< the coefficients as strings */
		bool rescale; /**< if true, the output is rescaled to [-1,1]  (to the same format as input) */
		mpz_class xHistory[10000]; // history of x used by emulate
		int currentIndex;
		FixSOPC *fixSOPC; /**< most of the work done here */
	};

}

#endif
