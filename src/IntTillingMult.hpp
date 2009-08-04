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
	IntTillingMult(Target* target, int wInX, int wInY,float ratio);

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
	/** The ratio between the number of DSPs and slices */
	float ratio;
	/* The working configuration of the tilling algorithm on DSPs */
	DSP** globalConfig;
	/* The best configuration of the after tilling DSPs */
	DSP** bestConfig;
	/* The target */
	Target * target;

	/* The number of estimated DSPs that will be used according to this parameter */
	int nrDSPs;

	int estimateDSPs();


	void tillingAlgorithm();
	


};


#endif
