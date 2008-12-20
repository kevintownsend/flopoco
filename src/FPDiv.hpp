#ifndef FPDIV_HPP
#define FPDIV_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FPNumber.hpp"
#include "FPDiv/SRT4Step.hpp"

/** The FPDiv class */
class FPDiv : public Operator
{
public:
	/**
	 * The FPDiv constructor
	 * @param[in]		target		the target device
	 * @param[in]		wE			the the with of the exponent for the f-p number X
	 * @param[in]		wF			the the with of the fraction for the f-p number X
	 */
	FPDiv(Target* target, int wE, int wF);

	/**
	 * FPDiv destructor
	 */
	~FPDiv();

	/**
	 * Method belonging to the Operator class overloaded by the FPDiv class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);


	/**
	 * Gets the signals which are interesting for TestCases.
	 * @see TestIOMap
	 */
	TestIOMap getTestIOMap();

	/**
	 * Gets the correct value associated to one or more inputs.
	 * @param a the array which contains both already filled inputs and
	 *          to be filled outputs in the order specified in getTestIOMap.
	 */
	void fillTestCase(mpz_class a[]);
	
private:
	/** The width of the exponent for the input X */
	int wE; 
	/** The width of the fraction for the input X */
	int wF; 
	/** The number of iterations */
	int nDigit;
	
	/** A SRT4Step subcomponent */
	SRT4Step* srt4step;
};

#endif //FPDIV_HPP
