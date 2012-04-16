#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco {

// TODO: it's a generic class, it should be put in a separate file
// or eliminated completely in favor of using directly std::auto_ptr
template<typename T> class Option {
	public:
		Option (T x)
			:empty(false),value(std::auto_ptr<T> (new T(x)))
		{
		}
		Option ()
			:empty(true),value(std::auto_ptr<T> (0))
		{
		}
		bool is_empty () {
			return empty;
		}
		T get_value () {
			if (empty)
				throw "Option object is empty";
			return *value;
		}
	protected:
		bool empty;
		std::auto_ptr<typename T> data;
};

class GenericBinaryPolynomial : public Operator {
  public:
    //static std::string operatorInfo;
    Product p;

  public:
    GenericBinaryPolynomial(Target* target, Product p, std::map<string,double> inputDelays = emptyMap);

    ~GenericBinaryPolynomial() {};

    void emulate(TestCase * tc);
    void buildStandardTestCases(TestCaseList* tcl);
    void buildRandomTestCases(TestCaseList* tcl, int n);
    TestCase* buildRandomTestCases(int i);
  protected:
    IntMultiAdder* ima;
};

}

