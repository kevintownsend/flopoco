#ifndef GENERICTESTBENCH_HPP
#define GENERICTESTBENCH_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntAdder.hpp"

namespace flopoco{

/** 
 * The GenericTestBench class. Produces an empty black box 
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

}
#endif
