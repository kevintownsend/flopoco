#ifndef FIXATAN2_HPP
#define FIXATAN2_HPP
#include <iomanip>
#include <vector>
#include <sstream>
#include <math.h>
#include <gmp.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "Operator.hpp"

#include "../Tools/Point.hpp"
#include "../Tools/Plane.hpp"


namespace flopoco {

	class FixAtan2 : public Operator {

	public:
		/**
		 * The FixAtan2 generic constructor computes atan(y/x), faithful to outLSB.
		 * 		the inputs x and y are assumed to be x in [0, 1) an y in [0, 1)
		 * @param[in] target            target device
		 * @param[in] wIn               input width
		 * @param[in] wOut              output width
		 **/
		FixAtan2(Target* target, int wIn, int wOut, map<string, double> inputDelays = emptyDelayMap);


		/**
		 *  Destructor
		 */
		~FixAtan2();

		/**
		 * Generates a component, and produces VHDL code for the instance inside an operator.
		 * The parent operator uses the std_logic_vector type.
		 */
		static FixAtan2* newComponentAndInstance(Operator* op,
													int wIn,
													int wOut
												);

		/**
		 * Generates a component, and produces VHDL code for the instance inside an operator.
		 * The parent operator uses the signed/unsigned types.
		 */
		static FixAtan2* newComponentAndInstanceNumericStd(Operator* op,
																int wIn,
																int wOut
															);

		/**
		 * The emulate function.
		 * @param[in] tc               a test-case
		 */
		void emulate ( TestCase* tc );

		void buildStandardTestCases(TestCaseList* tcl);

		int wIn;                     /**< input width */
		int wOut;                    /**< output width */

	private:
		Target* target;

	};

}
#endif
