#ifndef IntNAdderS_HPP
#define IntNAdderS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "Operator.hpp"


namespace flopoco{
	extern map<string, double> emptyDelayMap;
	/** The IntNAdder class for experimenting with adders. 
	 */
	class IntNAdder : public Operator
	{
	public:
		/**
		 * The IntNAdder constructor
		 * @param[in] target the target device
		 * @param[in] wIn    the with of the inputs and output
		 * @param[in] inputDelays the delays for each input
		 **/
		IntNAdder(Target* target, int wIn, int N, map<string, double> inputDelays = emptyDelayMap);
		/*IntNAdder(Target* target, int wIn);
		  void cmn(Target* target, int wIn, map<string, double> inputDelays);*/
	
		/**
		 *  Destructor
		 */
		~IntNAdder();


		void emulate(TestCase* tc);

	protected:
		int wIn_;                         /**< the width for X, Y and R*/
		int N_;

	private:
		map<string, double> inputDelays_; /**< a map between input signal names and their maximum delays */
		int bufferedInputs;               /**< variable denoting an initial buffering of the inputs */
		double maxInputDelay;             /**< the maximum delay between the inputs present in the map*/
		int nbOfChunks;                   /**< the number of chunks that the addition will be split in */
		int chunkSize_;                   /**< the suggested chunk size so that the addition can take place at the objective frequency*/
		int *cSize;                       /**< array containing the chunk sizes for all nbOfChunks*/

	};
}
#endif
