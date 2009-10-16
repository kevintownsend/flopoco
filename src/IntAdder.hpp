#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"


namespace flopoco{

	extern map<string, double> emptyDelayMap;
	/** The IntAdder class for experimenting with adders. 
	 */
	class IntAdder : public Operator
	{
	public:
		/**
		 * The IntAdder constructor
		 * @param[in] target the target device
		 * @param[in] wIn    the with of the inputs and output
		 * @param[in] inputDelays the delays for each input
		 * @param[in] aType	the type of adder we want to instantiate. 
		 *			For now we have 2 versions available, 
		 *			0 = fastest CPA, 1 = previous version
		 **/
		IntAdder(Target* target, int wIn, map<string, double> inputDelays = emptyDelayMap, int aType = 0 );
	
		/**
		 *  Destructor
		 */
		~IntAdder();

		/**
			A method returning the size of last chunk. 
			Typical use: target->adderDelay(getLastChunkSize())  provides the output slack delay
		*/
		int getLastChunkSize() {
			return lastChunkSize;
		}


		void emulate(TestCase* tc);

	protected:
		int wIn_;                         /**< the width for X, Y and R*/

	private:
		map<string, double> inputDelays_; /**< a map between input signal names and their maximum delays */
		double maxInputDelay;             /**< the maximum delay between the inputs present in the map*/
		int nbOfChunks;                   /**< the number of chunks that the addition will be split in */
		int chunkSize_;                   /**< the suggested chunk size so that the addition can take place at the objective frequency*/
		int lastChunkSize;
		int *cSize;                       /**< array containing the chunk sizes for all nbOfChunks*/
		int *cIndex;                       /**< array containing the indexes for all Chunks*/

	};

}
#endif
