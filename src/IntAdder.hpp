#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco {
	
#define LOGIC      0
#define REGISTER   1
#define SLICE      2
#define LATENCY    3
	
	extern map<string, double> emptyDelayMap;
	/** The IntAdder class for experimenting with adders.
	 */
	class IntAdder : public Operator {
	public:
		/**
		 * The IntAdder constructor
		 * @param[in] target           the target device
		 * @param[in] wIn              the with of the inputs and output
		 * @param[in] inputDelays      the delays for each input
		 * @param[in] optimizeType     the type optimization we want for our adder.
		 *            0: optimize for logic (LUT/ALUT)
		 *            1: optimize register count
		 *            2: optimize slice/ALM count
		 * @param[in] srl              optimize for use of shift registers
		 **/
		IntAdder ( Target* target, int wIn, map<string, double> inputDelays = emptyDelayMap, int optimizeType = SLICE, bool srl = true, int implementation = -1 );
			

			
		/**
		 *  Destructor
		 */
		~IntAdder();
		
		/** Overloaded to do nothing */
		void outputVHDL(std::ostream& o, std::string name);

		void changeName(std::string operatorName);
				
		/**
		 * The emulate function.
		 * @param[in] tc               a list of test-cases
		 */
		void emulate ( TestCase* tc );
				
	protected:
		int wIn_;                         /**< the width for X, Y and R*/
			
	private:
		vector<Operator*> addImplementationList;
		int selectedVersion;
	};
	
}
#endif
