#ifndef LONGACC2FP_HPP
#define LONGACC2FP_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Operator.hpp"
#include "Shifters.hpp"
#include "LZOCShifterSticky.hpp"
#include "IntAdder.hpp"


/** Operator which converts the output of the long accumulator to the desired FP format
 */
class LongAcc2FP : public Operator
{
public:

	/** Constructor
	 * @param target the target device
	 * @param MaxMSBX the weight of the MSB of the expected exponent of X
	 * @param LSBA the weight of the least significand bit of the accumulator
	 * @param MSBA the weight of the most significand bit of the accumulator
 	 * @param wEOut the width of the output exponent 
	 * @param eFOut the width of the output fractional part
	 */ 
	LongAcc2FP(Target* target, int LSBA, int MSBA, int wEOut, int wFOut);

	/** Destructor */
	~LongAcc2FP();

	/**
	 * Method belonging to the Operator class overloaded by the LongAcc2FP class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(ostream& o, string name);


protected:
	int LSBA_;    /**< the weight of the least significand bit of the accumulator */
	int MSBA_;    /**< the weight of the most significand bit of the accumulator */
	int wEOut_;   /**< the width of the output exponent */
	int wFOut_;   /**< the width of the output fractional part */
	int extraPipeLevel;

private:
	IntAdder* adder_;
	LZOCShifterSticky* lzocShifterSticky_;   
	int      sizeAcc_;       /**< The size of the accumulator  = MSBA-LSBA+1; */
	int      expBias_;       /**< the exponent bias value */
	int      countWidth_;    /**< the number of bits that the leading zero/one conunter outputs the result on */

};
#endif
