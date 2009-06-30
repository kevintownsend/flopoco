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

	
	
	/** The integer multiplier object for computing the product of X segments*/
	IntMultiplier* multiplierXSegv1;
	
	/** The integer multiplier object for computing the product of Y segments*/
	IntMultiplier* multiplierYSegv1;
	
	
	/** The integer multiplier object for computing the product of Z segments*/
	IntMultiplier* multiplierZSegv1;
	
	/** The  2 integers for adding the products of the segments for computing variable1 */
	IntNAdder* adder4var1;
	/** Fix2FP object for variable 1 */
	Fix2FP* convert2FPv1;
	
	
	/** The squarer for the X part of the variable 2 */
	IntSquarer* sqrXv2;
	/** The squarer for the Y part of the variable 2 */
	IntSquarer* sqrYv2;
	/** The squarer for the Z part of the variable 2 */
	IntSquarer* sqrZv2;
	/** The  2 integers for adding the products of the segments for computing the the intermidiate sum for var2 */
	IntNAdder* adder4SQRTv2;
	/** Fix2FP object for variable 2 ; input to the sqrt*/
	Fix2FP* convert2FP4sqrtv2;
	/** The SQRT for the variable 2 */
	FPSqrt* sqrt4var2;

	
	
	/** The multiplier for the (xo-x1)*t in variable 3*/
	IntMultiplier* multiplierXv3;
	/** The multiplier for the (yo-y1)*t in variable 3*/
	IntMultiplier* multiplierYv3;
	/** The multiplier for the (zo-z1)*t in variable 3*/
	IntMultiplier* multiplierZv3;
	
	/** The squarer for the X part of the variable 3 */
	IntSquarer* sqrXv3;
	/** The squarer for the Y part of the variable 3 */
	IntSquarer* sqrYv3;
	/** The squarer for the Z part of the variable 3 */
	IntSquarer* sqrZv3;
	/** The  2 integers for adding the products of the segments for computing the the intermidiate sum for var3 */
	IntNAdder* adder4SQRTv3;
	/** Fix2FP object for variable 3 ; input to the sqrt*/
	Fix2FP* convert2FP4sqrtv3;
	/** The SQRT for the variable 3 */
	FPSqrt* sqrt4var3;
	
	
	/** The multiplier for the (x1-x0)*t in variable 4*/
	IntMultiplier* multiplierXv4;
	/** The multiplier for the (y1-y0)*t in variable 4*/
	IntMultiplier* multiplierYv4;
	/** The multiplier for the (z1-z0)*t in variable 4*/
	IntMultiplier* multiplierZv4;
	
	/** The squarer for the X part of the variable 4 */
	IntSquarer* sqrXv4;
	/** The squarer for the Y part of the variable 4 */
	IntSquarer* sqrYv4;
	/** The squarer for the Z part of the variable 4 */
	IntSquarer* sqrZv4;
	/** The  2 integers for adding the products of the segments for computing the the intermidiate sum for var4 */
	IntNAdder* adder4SQRTv4;
	/** Fix2FP object for variable 4 ; input to the sqrt*/
	Fix2FP* convert2FP4sqrtv4;
	/** The SQRT for the variable 4 */
	FPSqrt* sqrt4var4;
	
	
	int MSB;
	int LSB;
	int inputWidth;
	int inputWidthSegments;
	int outputWidth;
	int wE,wF;
	int integratorWidth;
	int referenceCycle1;
	int referenceCycle2;

};

#endif

