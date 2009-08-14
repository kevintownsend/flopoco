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
#include "../../FPAdder.hpp"
#include "../../FPDiv.hpp"
#include "../../FPLog.hpp"
#include "../../FPMultiplier.hpp"
#include "../../LongAcc.hpp"
#include "../../ConstMult/IntConstMult.hpp"
#include "CoordinatesTableX.hpp"
#include "CoordinatesTableY.hpp"
#include "CoordinatesTableZ.hpp"



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
	CoilInductance(Target* target, int LSBI, int MSBI, int MaxMSBO,int LSBO, int MSBO, char *filepath);


	void outputVHDL(std::ostream& o, std::string name);
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
	/** The width of the MaxMSB for the output result*/
	int MaxMSBO;
	/** The width of the exponent for the output R */
	int LSBO; 
	/** The width of the fraction for the output R */
	int MSBO;
	/** configuration file path from which to initialize the parameters of the coil*/
	char *filepath;

	

	/** 3 input adder for variables */
	IntNAdder* adder4var;

	/** 3 input adder for each coordinate of var 5 */
	IntNAdder* adder4coordvar5;

	/** Fix2FP converter */
	Fix2FP* convert2FP;
	
	/** The  integer adder for the products of the segments for computing the the intermidiate sum for sqrt */
	IntNAdder* adder4SQRTv;
	
	/** Fix2FP object for variables  ; input to the sqrt*/
	Fix2FP* convert2FP4sqrtv;

	/** The SQRT for the variables  */
	FPSqrt* sqrt4var;

	/** The FPAdder object for creating the accumulator for the variables 3 & 4*/
	FPAdder* acc4var;
	
	/** The FP Division operator used in the final accumulator formula */
	FPDiv* div4Log;
	
	/** The FP Logarithm operator used in the final accumulator formula */
	FPLog* log4Acc;
	
	/** The FP Multiplier operator used in the final accumulator formula */
	FPMultiplier* mult4Acc;
	
	/** The FP Multiplier operator used in the final accumulator formula */
	LongAcc* finalAcc;
	
		
	/**The memories for the coordinates(X,Y,Z)*/
	CoordinatesTableX* memCoordX;
	CoordinatesTableY* memCoordY;
	CoordinatesTableZ* memCoordZ;
	
	
	

	int addressLength();
	
	int MSB;
	int LSB;
	int inputWidth;
	int inputWidthSegments;
	int outputWidth;
	int wE,wF;
	int integratorWidth;
	int addrWidth;

};

#endif

