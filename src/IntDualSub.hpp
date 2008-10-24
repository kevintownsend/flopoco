#ifndef IntDualSub_HPP
#define IntDualSub_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "Operator.hpp"

extern map<string, double> emptyDelayMap;
/** The IntDualSub class for experimenting with adders. 
*/
class IntDualSub : public Operator
{
public:
	/**
	 * The IntDualSub constructor
	 * @param[in] target the target device
	 * @param[in] wIn    the with of the inputs and output
	 * @param[in] inputDelays the delays for each input
	 **/
	IntDualSub(Target* target, int wIn, map<string, double> inputDelays = emptyDelayMap);
	/*IntDualSub(Target* target, int wIn);
	void cmn(Target* target, int wIn, map<string, double> inputDelays);*/
	
	/**
	 *  Destructor
	 */
	~IntDualSub();

	/**
	 * Method belonging to the Operator class overloaded by the IntDualSub class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);

	/**
	 * Gets the signals which are interesting for TestCases.
	 * @see TestIOMap
	 */
	TestIOMap getTestIOMap();
	
	/**
	 * Gets the correct value associated to one or more inputs.
	 * @param a the array which contains both already filled inputs and
	 *          to be filled outputs in the order specified in getTestIOMap.
	 */
	void fillTestCase(mpz_class a[]);

	/** 
	 * Sets the default name of this operator
	 */
	void setOperatorName(); 
	 
protected:
	int wIn_;                         /**< the width for X, Y and R*/

private:
	map<string, double> inputDelays_; /**< a map between input signal names and their maximum delays */
	int bufferedInputs;               /**< variable denoting an initial buffering of the inputs */
	double maxInputDelay;             /**< the maximum delay between the inputs present in the map*/
	int nbOfChunks;                   /**< the number of chunks that the addition will be split in */
	int chunkSize_;                   /**< the suggested chunk size so that the addition can take place at the objective frequency*/
	int *cSize;                       /**< array containing the chunk sizes for all nbOfChunks*/

};
#endif
