#ifndef IntCompressorTree_HPP
#define IntCompressorTree_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "Operator.hpp"

extern map<string, double> emptyDelayMap;
/** The IntCompressorTree class for experimenting with adders. 
*/
class IntCompressorTree : public Operator
{
public:
	/**
	 * The IntCompressorTree constructor
	 * @param[in] target the target device
	 * @param[in] wIn    the with of the inputs and output
	 * @param[in] inputDelays the delays for each input
	 **/
	IntCompressorTree(Target* target, int wIn, int N, map<string, double> inputDelays = emptyDelayMap);
	
	/**
	 *  Destructor
	 */
	~IntCompressorTree();


	void emulate(TestCase* tc);

protected:
	int wIn_;                         /**< the width for X_{0}..X_{n-1} and R*/
	int N_;                           /**< number of operands */

private:
	map<string, double> inputDelays_; /**< a map between input signal names and their maximum delays */
	int bufferedInputs;               /**< variable denoting an initial buffering of the inputs */
	double maxInputDelay;             /**< the maximum delay between the inputs present in the map*/
	int nbOfChunks;                   /**< the number of chunks that the addition will be split in */
	int chunkSize_;                   /**< the suggested chunk size so that the addition can take place at the objective frequency*/
	int *cSize;                       /**< array containing the chunk sizes for all nbOfChunks*/

};
#endif
