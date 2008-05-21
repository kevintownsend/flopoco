#ifndef __FPEXP_HPP
#define __FPEXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

class Fragment;

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
	int wE, wF;
	Fragment *f;
};

#endif
