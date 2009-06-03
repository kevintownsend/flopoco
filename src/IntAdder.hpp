#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "Operator.hpp"

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
	 **/
	IntAdder(Target* target, int wIn, map<string, double> inputDelays = emptyDelayMap);
	/*IntAdder(Target* target, int wIn);
	void cmn(Target* target, int wIn, map<string, double> inputDelays);*/
	
	/**
	 *  Destructor
	 */
	~IntAdder();

	void emulate(TestCase* tc);

protected:
	int wIn_;                         /**< the width for X, Y and R*/

private:
	map<string, double> inputDelays_; /**< a map between input signal names and their maximum delays */
	double maxInputDelay;             /**< the maximum delay between the inputs present in the map*/
	int nbOfChunks;                   /**< the number of chunks that the addition will be split in */
	int chunkSize_;                   /**< the suggested chunk size so that the addition can take place at the objective frequency*/
	int *cSize;                       /**< array containing the chunk sizes for all nbOfChunks*/
	int *cIndex;                       /**< array containing the indexes for all Chunks*/

};
#endif
