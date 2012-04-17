#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "IntPower.hpp"
#include "FormalBinaryProduct.hpp"

//std::string IntPower::operatorInfo = "UserDefinedInfo param0 param1 <options>";

using namespace flopoco;

IntPower::IntPower(Target* target,
                   size_t wIn, size_t n,
		   std::map<std::string,double> inputDelays)
	:GenericBinaryPolynomial(target,
			ProductIR::identity(wIn).toPow(n),
			inputDelays) {
};

	
void IntPower::emulate(TestCase * tc) {
}

void IntPower::buildStandardTestCases(TestCaseList * tcl) {
}

void IntPower::buildRandomTestCases(TestCaseList *  tcl, int n) {
}

TestCase* IntPower::buildRandomTestCases(int i) {
  TestCase* tc = new TestCase(this);
  return tc;
}
