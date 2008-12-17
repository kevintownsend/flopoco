#ifndef IntMult2_HPP
#define IntMult2_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
/** 
 * The Integer Multiplier class. Receives at input two numbers of 
 * wInX and wInY widths and outputs a result having the width wOut=wInX+wInY 
 **/
class IntMult2 : public Operator
{
public:
	/** 
	 * The constructor of the IntMult2 class
	 * @param target argument of type Target containing the data for which this operator will be optimized
	 * @param wInX integer argument representing the width in bits of the input X 
	 * @param wInY integer argument representing the width in bits of the input Y
	 **/
	IntMult2(Target* target, int wInX, int wInY);
	
	/** IntMult2 destructor */
	~IntMult2();
	
	/**
	 * Method belonging to the Operator class overloaded by the IntMult2 class
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

	int wInX_; /**< the width (in bits) of the input X  */
	int wInY_; /**< the width (in bits) of the input Y  */
	int wOut_; /**< the width (in bits) of the output R  */

private:

	int      partsX; 	          /**< The number of parts that the input X will be split in */
	int      partsY; 	          /**< The number of parts that the input Y will be split in */
	int      size;
};
#endif
