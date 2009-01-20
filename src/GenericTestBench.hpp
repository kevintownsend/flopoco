#ifndef GENERICTESTBENCH_HPP
#define GENERICTESTBENCH_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntAdder.hpp"

/** 
 * The Integer Multiplier class. Receives at input two numbers of 
 * wInX and wInY widths and outputs a result having the width wOut=wInX+wInY 
 **/
class GenericTestBench : public Operator
{
public:
	/** 
	 * The constructor of the GenericTestBench class
	 * @param target argument of type Target containing the data for which this operator will be optimized
	 * @param wInX integer argument representing the width in bits of the input X 
	 * @param wInY integer argument representing the width in bits of the input Y
	 **/
	GenericTestBench(Target* target, string name, int depth, int wIn);
	
	/** GenericTestBench destructor */
	~GenericTestBench();
	
	/**
	 * Method belonging to the Operator class overloaded by the GenericTestBench class
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

	int wIn_; /**< the width (in bits) of the input X  */
	string _name;

};
#endif
