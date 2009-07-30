#ifndef IntTillingMult_HPP
#define IntTillingMult_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

/** The IntTillingMult class */
class IntTillingMult : public Operator
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
	IntTillingMult(Target* target, int wInX, int wInY);

	/**
	 * IntTillingMult destructor
	 */
	~IntTillingMult();


	void emulate(TestCase * tc);
	void buildStandardTestCases(TestCaseList* tcl);



private:
	/** The width of first input*/
	int wInX; 
	/** The width of second input*/
	int wInY; 
	/** The width of output */
	int wOut;
	
	


};


#endif
