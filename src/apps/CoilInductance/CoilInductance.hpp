#ifndef CoilInductance_HPP
#define CoilInductance_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../../Operator.hpp"
#include "../../FPNumber.hpp"
#include "../../IntAdder.hpp"
#include "../../IntMultiplier.hpp"
#include "../../Fix2FP.hpp"
#include "../../IntSquarer.hpp"
#include "../../FPSqrt.hpp"



/** The FPAdder class */
class CoilInductance : public Operator
{
public:
	/**
	 * The FPAdder constructor
	 * @param[in]		target		the target device
	 * @param[in]		MSB			the MSB of the input number
	 * @param[in]		LSB			the LSB of the input number
	 * @param[in]		wER			the with of the exponent for the convertion result
	 * @param[in]		wFR			the with of the fraction for the convertion result
	 */
	CoilInductance(Target* target, int LSBI, int MSBI, int LSBO, int MSBO, char *filepath);

	/**
	 * FPAdder destructor
	 */
	~CoilInductance();


	void emulate(TestCase * tc);
	void buildStandardTestCases(TestCaseList* tcl);



private:
	/** The MSB for the input */
	int MSBI; 
	/** The LSB for the input */
	int LSBI; 
	/** The width of the exponent for the output R */
	int LSBO; 
	/** The width of the fraction for the output R */
	int MSBO;
	/** configuration file path from which to initialize the parameters of the coil*/
	char *filepath;

	
	/** The integer adder object for computing the first segment with coordinates X*/
	IntAdder* segment1X;
	/** The integer adder object for computing the second segment with coordinates X*/
	IntAdder* segment2X;
	/** The integer adder object for computing the first segment with coordinates Y*/
	IntAdder* segment1Y;
	/** The integer adder object for computing the second segment with coordinates Y*/
	IntAdder* segment2Y;
	/** The integer adder object for computing the first segment with coordinates Z*/
	IntAdder* segment1Z;
	/** The integer adder object for computing the second segment with coordinates Z*/
	IntAdder* segment2Z;


	/** The integer adder object for converting the difference x1-x0 into positive*/
	IntAdder* covertInitialXS1v1;
	/** The integer adder object for converting the difference x3-x2 into positive*/
	IntAdder* covertInitialXS2v1;
	/** The integer adder object for converting back the product of X segments*/
	IntAdder* covertFinalXv1;
	/** The integer multiplier object for computing the product of X segments*/
	IntMultiplier* multiplierXSegv1;
	
	/** The integer adder object for converting the difference y1-y0 into positive*/
	IntAdder* covertInitialYS1v1;
	/** The integer adder object for converting the difference y3-y2 into positive*/
	IntAdder* covertInitialYS2v1;
	/** The integer adder object for converting back the product of Y segments*/
	IntAdder* covertFinalYv1;
	/** The integer multiplier object for computing the product of Y segments*/
	IntMultiplier* multiplierYSegv1;
	
	
	/** The integer adder object for converting the difference z1-z0 into positive*/
	IntAdder* covertInitialZS1v1;
	/** The integer adder object for converting the difference z3-z2 into positive*/
	IntAdder* covertInitialZS2v1;
	/** The integer adder object for converting back the product of Z segments*/
	IntAdder* covertFinalZv1;
	/** The integer multiplier object for computing the product of Z segments*/
	IntMultiplier* multiplierZSegv1;
	
	/** The  2 integers for adding the products of the segments for computing variable1 */
	IntAdder* adderPartialResult4Var1;
	IntAdder* adderResult4Var1;
	/** Fix2FP object for variable 1 */
	Fix2FP* convert2FPv1;
	
	
	/** The squarer for the X part of the variable 2 */
	IntSquarer* sqrXv2;
	/** The squarer for the Y part of the variable 2 */
	IntSquarer* sqrYv2;
	/** The squarer for the Z part of the variable 2 */
	IntSquarer* sqrZv2;
	/** The  2 integers for adding the products of the segments for computing the the intermidiate sum for var2 */
	IntAdder* adderPartialResult4SQRTv2;
	IntAdder* adderResult4SQRTv2;
	/** Fix2FP object for variable 2 ; input to the sqrt*/
	Fix2FP* convert2FP4sqrtv2;
	/** The SQRT for the variable 2 */
	FPSqrt* sqrt4var2;

	
	
	
	
	int MSB;
	int LSB;
	int inputWidth;
	int inputWidthSegments;
	int outputWidth;
	int wE,wF;
	

};

#endif

