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

	//~ /**The leading sign counter	*/
	//~ LZOCShifterSticky* lzocs; 
	//~ /**The leading zero counter	*/
	//~ LZOCShifterSticky* lzcs; 
	//~ /** The integer adder object for subtraction from the MSB the position of the leading 1, for shifting the number */
	//~ IntAdder* fractionConvert;
	//~ /** The integer adder object for adding 1 to the fraction part*/
	//~ IntAdder* roundingAdder;
	//~ /** The integer adder object for substracting 1 from the remainder of the fraction to establish if it is zero*/
	//~ IntAdder* oneSubstracter;
	//~ /** The integer adder object for transforming the Count of LZO to exponent*/
	//~ IntAdder* exponentConversion;
	//~ /** The integer adder object for adding the biass to the exponent*/
	//~ IntAdder* exponentFinal;
	//~ /** The integer adder object for zero detector of the input*/
	//~ IntAdder* zeroD;
	//~ /** The integer adder object for correcting the exponent*/
	//~ IntAdder* expCorrect;
	
	
	
	
	
	int MSB;
	int LSB;
	int inputWidth;
	

};

#endif

