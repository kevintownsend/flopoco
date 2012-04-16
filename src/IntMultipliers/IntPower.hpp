#include "Operator.hpp"
#include "utils.hpp"
#include "GenericBinaryPolynomial.hpp"

namespace flopoco {

class IntPower : public GenericBinaryPolynomial {
  public:
    IntPower(Target* target, size_t wIn, size_t n, std::map<string,double> inputDelays = emptyMap);

    ~IntPower() {};

    void emulate(TestCase * tc);
    void buildStandardTestCases(TestCaseList* tcl);
    void buildRandomTestCases(TestCaseList* tcl, int n);
    TestCase* buildRandomTestCases(int i);
  protected:
};

}

