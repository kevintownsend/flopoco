#ifndef __HOTBM_HPP
#define __HOTBM_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

class HOTBMInstance;
class Function;

/**
 * Implements an Operator around HOTBM. Acts like a wrapper around
 * HOTBM classes.
 */
class HOTBM : public Operator
{
public:
	HOTBM(Target* target, string func, int wI, int wO, int n);
	~HOTBM();

	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);

	TestCaseList generateStandardTestCases(int n);
	TestCaseList generateRandomTestCases(int n);

private:
	HOTBMInstance *inst;
	Function &f;

	int wI, wO;
};

#endif
