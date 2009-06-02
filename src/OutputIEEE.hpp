#ifndef OUTPUTIEEE_HPP
#define OUTPUTIEEE_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FPNumber.hpp"

/** The OutputIEEE class */
class OutputIEEE : public Operator
{
public:
	/**
	 * The OutputIEEE constructor
	 * @param[in]		target		the target device
	 * @param[in]		wE			the the with of the exponent for the f-p number X
	 * @param[in]		wF			the the with of the fraction for the f-p number X
	 */
	OutputIEEE(Target* target, int wEI, int wFI, int wEO, int wFO);

	/**
	 * OutputIEEE destructor
	 */
	~OutputIEEE();



	/**
	 * Emulate a correctly rounded square root using MPFR.
	 * @param tc a TestCase partially filled with input values 
	 */
	void emulate(TestCase * tc);

	
private:
	/** The width of the exponent for the input X */
	int wEI; 
	/** The width of the fraction for the input X */
	int wFI; 
	/** The width of the exponent for the output X */
	int wEO; 
	/** The width of the fraction for the output X */
	int wFO; 
}
;

#endif //OUTPUTIEEE_HPP
