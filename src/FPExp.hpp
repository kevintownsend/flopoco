#ifndef __FPEXP_HPP
#define __FPEXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

class Fragment;
class FPNumber;


namespace flopoco{

	class FPExp : public Operator
	{
	public:
		FPExp(Target* target, int wE, int wF);
		~FPExp();
		
		// Overloading the virtual functions of Operator
		// void outputVHDL(std::ostream& o, std::string name);
		
		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);
		
	private:
		int wE, wF;
		int p; // Size of the address bits for the first table
		int result_length, g;
		// Fragment *f;
		double area, max_error;
	};

}
#endif
