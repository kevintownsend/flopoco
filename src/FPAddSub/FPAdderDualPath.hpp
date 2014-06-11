#ifndef FPADDERDUALPATH_HPP
#define FPADDERDUALPATH_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../Operator.hpp"
#include "../ShiftersEtc/LZOC.hpp"
#include "../ShiftersEtc/Shifters.hpp"
#include "../ShiftersEtc/LZOCShifterSticky.hpp"
#include "../TestBenches/FPNumber.hpp"
#include "../IntAddSubCmp/IntAdder.hpp"
#include "../IntAddSubCmp/IntDualSub.hpp"

namespace flopoco{

	/** The FPAdderDualPath class */
	class FPAdderDualPath : public Operator
	{
	public:
		/**
		 * The FPAdderDualPath constructor
		 * @param[in]		target		the target device
		 * @param[in]		wEX			the the with of the exponent for the f-p number X
		 * @param[in]		wFX			the the with of the fraction for the f-p number X
		 * @param[in]		wEY			the the with of the exponent for the f-p number Y
		 * @param[in]		wFY			the the with of the fraction for the f-p number Y
		 * @param[in]		wER			the the with of the exponent for the addition result
		 * @param[in]		wFR			the the with of the fraction for the addition result
		 */
		FPAdderDualPath(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR);

		/**
		 * FPAdderDualPath destructor
		 */
		~FPAdderDualPath();


		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);
		TestCase* buildRandomTestCase(int i);



	private:
		/** The width of the exponent for the input X */
		int wEX; 
		/** The width of the fraction for the input X */
		int wFX; 
		/** The width of the exponent for the input Y */
		int wEY; 
		/** The width of the fraction for the input Y */
		int wFY; 
		/** The width of the exponent for the output R */
		int wER; 
		/** The width of the fraction for the output R */
		int wFR;
		/** Signal if the output of the operator is to be or not normalized*/

		/** The combined leading zero counter and shifter for the close path */
		LZOCShifterSticky* lzocs; 
		/** The integer adder object for subtraction in the close path */
		IntAdder *fracSubClose; 
		/** The dual subtractor for the close path */
		IntDualSub *dualSubClose;
		/** The fraction adder for the far path */
		IntAdder *fracAddFar; 
		/** The adder that does the final rounding */
		IntAdder *finalRoundAdd; 
		/** The right shifter for the far path */
		Shifter* rightShifter;	


		int wF;
		int wE;
		int sizeRightShift;
	
	};

}

#endif 
