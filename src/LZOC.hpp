#ifndef LZOC_HPP
#define LZOC_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Operator.hpp"

extern map<string, double> emptyDelayMap;
/** The Leading zero/one counter class.  
 * Recursive structure with wOut stages. At most 2^wOut-1 leading zeros are counted.
 */
class LZOC : public Operator{

public:
	/** The LZOC constructor
	 * @param[in] target the target device for this operator
	 * @param[in] wIn the width of the input
	 * @param[in] wOut the width of the output 
	 */
	LZOC(Target* target, int wIn, map<string, double> inputDelays = emptyDelayMap);
	
	/** The LZOC destructor	*/
	~LZOC();

	/**
	 * Method belonging to the Operator class overloaded by the LZOC class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);

	/** 
	 * Sets the default name of this operator
	 */
	void setOperatorName();


	/**
	 * Emulate a correctly rounded division using MPFR.
	 * @param tc a TestCase partially filled with input values 
	 */
	void emulate(TestCase * tc);
	
protected:

	int wIn_;    /**< The width of the input */
	int wOut_;   /**< The width of the output */
	int p2wOut_; /**< The value of 2^wOut, which is computed as 1<<wOut */

};
#endif
