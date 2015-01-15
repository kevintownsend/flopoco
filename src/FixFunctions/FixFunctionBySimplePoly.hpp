#ifndef FixFunctionBySimplePoly_HPP
#define FixFunctionBySimplePoly_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../Operator.hpp"
#include "FixFunction.hpp"
#include "BasicPolyApprox.hpp"

namespace flopoco{


	/** The FixFunctionBySimplePoly class */
	class FixFunctionBySimplePoly : public Operator
	{
	public:
		/**
		 * The FixFunctionBySimplePoly constructor
			 @param[string] func    a string representing the function, input range should be [0,1]
			 @param[int]    lsbIn   input LSB weight
			 @param[int]    msbOut  output MSB weight, used to determine wOut
			 @param[int]    lsbOut  output LSB weight
			 @param[int]    lsbOut  output LSB weight
			 @param[bool]   finalRounding: if false, the operator outputs its guard bits as well, saving the half-ulp rounding error. 
			                 This makes sense in situations that further process the result with further guard bits.
			 @param[bool]   plainStupidVHDL: if true, generate * and +; if false, use BitHeap-based FixMultAdd
			 @param[float]   DSPThreshold between 0 and 1, 0 means logic only, 1 means use DSPs even for very small multiplications
			 
			 One could argue that MSB weight is redundant, as it can be deduced from an analysis of the function. 
			 This would require quite a lot of work for non-trivial functions (isolating roots of the derivative etc).
			 So this is currently left to the user.
		 */
		FixFunctionBySimplePoly(Target* target, string func, int lsbIn, int msbOut, int lsbOut, bool finalRounding = true, bool plainStupidVHDL=false, float DSPThreshold=0.7,  map<string, double> inputDelays = emptyDelayMap);

		/**
		 * FixFunctionBySimplePoly destructor
		 */
		~FixFunctionBySimplePoly();
		
		void emulate(TestCase * tc);

		void buildStandardTestCases(TestCaseList* tcl);

	private:
		FixFunction *f; 
		BasicPolyApprox *poly;
		bool finalRounding;

		vector<int> coeffMSB; 
		vector<int> coeffLSB; 
		vector<int> coeffSize; 
	};

}

#endif
