#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco{

#define LOGIC      0 
#define REGISTER   1
#define SLICE      2
#define LATENCY    3


#define CLASSICAL    0
#define ALTERNATIVE  1
#define SHORTLATENCY 2


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
		 * @param[in] the type optimization we want for our adder.
		 *            0: optimize for logic (LUT/ALUT)
		 *            1: optimize register count
		 *            2: optimize slice/ALM count
		 **/
		IntAdder(Target* target, int wIn, map<string, double> inputDelays = emptyDelayMap, int optimizeType = SLICE );
	
		int implementationSelector(Target* target, int wIn, map<string, double> inputDelays, int optimizeType);

		int getLutCostClassical(Target* target, int wIn, map<string, double> inputDelays);
		int getLutCostAlternative(Target* target, int wIn, map<string, double> inputDelays);
		int getLutCostShortLatency(Target* target, int wIn, map<string, double> inputDelays);
	
		int getRegCostClassical(Target* target, int wIn, map<string, double> inputDelays);
		int getRegCostAlternative(Target* target, int wIn, map<string, double> inputDelays);
		int getRegCostShortLatency(Target* target, int wIn, map<string, double> inputDelays);

		int getSliceCostClassical(Target* target, int wIn, map<string, double> inputDelays);
		int getSliceCostAlternative(Target* target, int wIn, map<string, double> inputDelays);
		int getSliceCostShortLatency(Target* target, int wIn, map<string, double> inputDelays);

		void updateParameters( Target* target, int &alpha, int &beta, int &k);
		void updateParameters( Target* target, map<string, double> inputDelays, int &alpha, int &beta, int &gamma, int &k);
		void updateParameters( Target* target, map<string, double> inputDelays, int &alpha, int &beta, int &k);
		
		bool tryOptimizedChunkSplittingShortLatency(Target* target, int wIn, int k);
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
		// new notations
		int alpha; /**< the chunk size */
		int beta;  /**< the last chunk size */
		int gamma; /**< the first chunk size when slack is considered */
		int k;     /**< the number of chunks */
		int w;     /**< the addition width */

		int selectedDesign;
		
		int classicalSlackVersion;
		int alternativeSlackVersion;
		
		int shortLatencyVersion;
		double objectivePeriod;
	};

}
#endif
