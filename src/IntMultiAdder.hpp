#ifndef IntMultiAdderS_HPP
#define IntMultiAdderS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco {
	
	extern map<string, double> emptyDelayMap;
	/** The IntMultiAdder class for experimenting with adders.
	 */
	class IntMultiAdder : public Operator {
	public:
	
		IntMultiAdder (Target* target, int wIn, int N, map<string, double> inputDelays, bool noAmbiguity): 
			Operator(target,inputDelays), wIn_(wIn){
		}
		
		/**
		 * The IntMultiAdder constructor
		 * @param[in] target           the target device
		 * @param[in] wIn              the with of the inputs and output
		 * @param[in] N                the number of operands
		 * @param[in] inputDelays      the delays for each input
		 **/
		IntMultiAdder ( Target* target, int wIn, int N, map<string, double> inputDelays = emptyDelayMap);

		/**
		 *  Destructor
		 */
		~IntMultiAdder();
		
		/**
		 * The emulate function.
		 * @param[in] tc               a list of test-cases
		 */
		void emulate ( TestCase* tc );
				
	protected:
		int wIn_; // the width
		int N_;   // the number of operands;
	};
	
}
#endif
