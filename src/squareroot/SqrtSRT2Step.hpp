#ifndef SQRTSRT2STEP_HPP
#define SQRTSRT2STEP_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../Operator.hpp"

/** The SqrtSRT2Step class */
class SqrtSRT2Step : public Operator
{
public:
	/**
	 * The SqrtSRT2Step constructor
	 * @param[in]		target		the target device
	 * @param[in]		wF			   the with of the fraction for the f-p number X
	 * @param[in]		step			the step
	 */
	SqrtSRT2Step(Target* target, int wF, int step);

	/**
	 * SqrtSRT2Step destructor
	 */
	~SqrtSRT2Step();

	
private:
	/** The width of the fraction for the input X */
	int wF; 

	/** The step of the digit recurrence */
	int step; 

};

#endif //SQRTSRT2STEP_HPP
