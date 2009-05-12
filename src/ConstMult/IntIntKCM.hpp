#ifndef IntIntKCMS_HPP
#define IntIntKCMS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "../Operator.hpp"

extern map<string, double> emptyDelayMap;
/** The IntIntKCM class for experimenting with adders. 
*/
class IntIntKCM : public Operator
{
public:
	/**
	 * The IntIntKCM constructor
	 * @param[in] target the target device
	 * @param[in] wIn    the with of the inputs and output
	 * @param[in] inputDelays the delays for each input
	 **/
	IntIntKCM(Target* target, int wIn, mpz_class C, map<string, double> inputDelays = emptyDelayMap);

	/**
	 *  Destructor
	 */
	~IntIntKCM();

	void emulate(TestCase* tc);
	
	int getOutputWidth();             /**< returns the number of bits of the output */

protected:
	int wIn_;                         /**< the width for the input X*/
	mpz_class C_;                           /**< the constant to be used for the multiplication*/
	int wOut_;
	
private:
	map<string, double> inputDelays_; /**< a map between input signal names and their maximum delays */
//	int bufferedInputs;               /**< variable denoting an initial buffering of the inputs */
//	double maxInputDelay;             /**< the maximum delay between the inputs present in the map*/
//	int nbOfChunks;                   /**< the number of chunks that the addition will be split in */
//	int chunkSize_;                   /**< the suggested chunk size so that the addition can take place at the objective frequency*/
//	int *cSize;                       /**< array containing the chunk sizes for all nbOfChunks*/

};
#endif
