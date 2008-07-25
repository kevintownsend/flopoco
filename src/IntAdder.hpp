#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "Operator.hpp"


/** The IntAdder class for experimenting with adders. 
*/
class IntAdder : public Operator
{
public:
	/**
	 * The IntAdder constructor
	 * @param[in] target the target device
	 * @param[in] wIn    the with of the inputs and output
	 * @param[in] p      if we want to force the pipeline //XXX Soon obsolete
	 **/
	IntAdder(Target* target, int wIn, const int p=0);
	
	/**
	 *  Destructor
	 */
	~IntAdder();

	/**
	 * Method belonging to the Operator class overloaded by the IntAdder class
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
	int wIn_;           /**< the width for X, Y and R*/

private:
	int chunkSize_;     /**< the maximum chunk size for an addition so that the requested frequency can theoretically be used*/
	int lastChunkSize_; /**< the last chunk size - the last one is slightly smaller*/
	int pipeLevels_;    /**< the number of pieline levels of the sequential operator*/ 
	int forcePipeline_; /**< force pipeline */
};
#endif
