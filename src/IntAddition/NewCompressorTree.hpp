#include "Operator.hpp"
#include "utils.hpp"

#include <vector>

using namespace flopoco;


class NewCompressorTree : public Operator {
	public:
		static string operatorInfo;
		unsigned w; // length of the result
		std::vector<unsigned> vert_operands;
		// at index i [little endian], contains
		// the number of bits of weight i to be added


	public:

		NewCompressorTree(Target* target,
		                  std::vector<unsigned> vert_ops);

		~NewCompressorTree() {};


		void emulate(TestCase * tc);

		void buildStandardTestCases(TestCaseList* tcl);

		void buildRandomTestCases(TestCaseList* tcl, int n);

		TestCase* buildRandomTestCases(int i);
};
