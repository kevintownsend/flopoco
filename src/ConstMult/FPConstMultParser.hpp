#ifndef FPCONSTMULTPARSER_HPP
#define FPCONSTMULTPARSER_HPP
#include "../Operator.hpp"

namespace flopoco{

	class FPConstMultParser : public FPConstMult
	{
	public:
		FPConstMultParser(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, int wF_C, string constant);
		~FPConstMultParser();

		int cst_width;

		string constant;

		// No need for the other usual methods of Operators, use those of FPConstMult

		// except for test case generation: TODO

	};

}
#endif
