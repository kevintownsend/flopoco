#ifndef __EXP_HPP
#define __EXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

class Fragment;

class Exp : public Operator
{
public:
	Exp(Target* target, int wE, int wF);
	~Exp();

	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);

	TestCaseList generateStandardTestCases(int n);
	TestCaseList generateRandomTestCases(int n);

private:
	int wE, wF;
	Fragment *f;
};

#endif
