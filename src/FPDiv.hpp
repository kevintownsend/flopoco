#ifndef FPDIV_HPP
#define FPDIV_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FPNumber.hpp"
#include "division/SRT4Step.hpp"

namespace flopoco{

	/** The FPDiv class */
	class FPDiv : public Operator
	{
	public:
		/**
		 * The FPDiv constructor
		 * @param[in]		target		the target device
		 * @param[in]		wE			the the with of the exponent for the f-p number X
		 * @param[in]		wF			the the with of the fraction for the f-p number X
		 */
		FPDiv(Target* target, int wE, int wF);

		/**
		 * FPDiv destructor
		 */
		~FPDiv();

		/**
		 * Method belonging to the Operator class overloaded by the FPDiv class
		 * @param[in,out] o     the stream where the current architecture will be outputed to
		 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
		 **/
		void outputVHDL(std::ostream& o, std::string name);


		/**
		 * Emulate a correctly rounded division using MPFR.
		 * @param tc a TestCase partially filled with input values 
		 */
		void emulate(TestCase * tc);

		/* Overloading the Operator method */
		void buildStandardTestCases(TestCaseList* tcl);

	
	private:
		/** The width of the exponent for the input X */
		int wE; 
		/** The width of the fraction for the input X */
		int wF; 
		/** The number of iterations */
		int nDigit;

		/** A SRT4Step subcomponent */
		SRT4Step* srt4step;
		int stepDelay;
		int srtDelay;

	};
}
#endif //FPDIV_HPP
