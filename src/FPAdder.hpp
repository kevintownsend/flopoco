#ifndef FPADDER_HPP
#define FPADDER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "LZOC.hpp"
#include "Shifters.hpp"
#include "FPNumber.hpp"
#include "IntAdder.hpp"
#include "IntDualSub.hpp"
#include "LZOCShifterSticky.hpp"

/** The FPAdder class */
class FPAdder : public Operator
{
public:
	/**
	 * The FPAdder constructor
	 * @param[in]		target		the target device
	 * @param[in]		wEX			the the with of the exponent for the f-p number X
	 * @param[in]		wFX			the the with of the fraction for the f-p number X
	 * @param[in]		wEY			the the with of the exponent for the f-p number Y
	 * @param[in]		wFY			the the with of the fraction for the f-p number Y
	 * @param[in]		wER			the the with of the exponent for the addition result
	 * @param[in]		wFR			the the with of the fraction for the addition result
	 */
	FPAdder(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR);

	/**
	 * FPAdder destructor
	 */
	~FPAdder();

	/**
	 * Method belonging to the Operator class overloaded by the FPAdder class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);

	/** Method for setting the operator name
	*/
	void setOperatorName();	

	/**
	 * Emulate a correctly rounded addition using MPFR.
	 * @param tc a TestCase partially filled with input values 
	 */
	void emulate(TestCase * tc);
	

	void buildStandardTestCases(TestCaseList* tcl);



private:
	/** The width of the exponent for the input X */
	int wEX; 
	/** The width of the fraction for the input X */
	int wFX; 
	/** The width of the exponent for the input Y */
	int wEY; 
	/** The width of the fraction for the input Y */
	int wFY; 
	/** The width of the exponent for the output R */
	int wER; 
	/** The width of the fraction for the output R */
	int wFR;
	/** Signal if the output of the operator is to be or not normalized*/

	/** The integer adder object */
	/** The integer adder object */
	IntAdder *fracSubClose; 
	/** The integer adder object */
	IntAdder *complementAdderClose; 

	/** The dual subtractor for the close path */
	IntDualSub *dualSubClose;
	/** The integer adder object */
	IntAdder *fracAddFar; 
	/** The integer adder object */
	IntAdder *finalRoundAdd; 

	LZOC* leadingZeroCounter;
	Shifter* leftShifter;
	Shifter* rightShifter;	
	LZOCShifterSticky* lzocs; 

	int wF;
	int wE;
	int wOutLZC;
	int sizeRightShift;
	
	int swapDifferencePipelineDepth;
	int closePathDepth;
	int farPathDepth;
	int maxPathDepth;	
	int delaySDToRound;
};

#endif
