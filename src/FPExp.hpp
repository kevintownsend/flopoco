#ifndef __FPEXP_HPP
#define __FPEXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

class Fragment;
class FloFP;

class FPExp : public Operator
{
public:
	FPExp(Target* target, int wE, int wF);
	~FPExp();

	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);

	TestCaseList generateStandardTestCases(int n);
	TestCaseList generateRandomTestCases(int n);

private:
	/**
	 * Adds a test case for the FPExp operator.
	 * @param tcl The TestCaseList to which the TestCase should be added.
	 * @param x The input value for which the TestCase should be generated.
	 */
	void addTestCase(TestCaseList& tcl, FloFP x);
	
	int wE, wF;
	Fragment *f;
	int result_length, g;
	double area, max_error;
};

#endif
