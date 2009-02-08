#ifndef __FPEXP_HPP
#define __FPEXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

class Fragment;
class FPNumber;

class FPExp : public Operator
{
public:
	FPExp(Target* target, int wE, int wF);
	~FPExp();

	// Overloading the virtual functions of Operator
	void outputVHDL(std::ostream& o, std::string name);

	void fillTestCase(mpz_class a[]);

private:
	int wE, wF;
	Fragment *f;
	int result_length, g;
	double area, max_error;
};

#endif
