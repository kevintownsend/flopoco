#ifndef FixSOPC_HPP
#define FixSOPC_HPP

#include "Operator.hpp"
#include "utils.hpp"

#include "BitHeap/BitHeap.hpp"

/*  All flopoco operators and utility functions are declared within
  the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
  functions.
*/

namespace flopoco{

	class FixSOPC : public Operator {
	public:

		/** simplest constructor for inputs in the fixed-point format (0, lsbIn), computing msbOut out of the coeffs */
		FixSOPC(Target* target, int lsbIn, int lsbOut, vector<string> coeff);

		/** constructor for inputs in various formats, computing msbOut  */
		FixSOPC(Target* target, vector<int> msbIn, vector<int> lsbIn, int lsbOut, vector<string> coeff);

		/** constructor for inputs in various formats, msbOut provided  */
		FixSOPC(Target* target, vector<int> msbIn, vector<int> lsbIn, int msbOut, int lsbOut, vector<string> coeff_);

		/** destructor */
		~FixSOPC() {};

		/** The method that does most of operator construction for the two constructors */
		void initialize();

		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
		   in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		/* function used to create Standard testCase defined by the developper */
		void buildStandardTestCases(TestCaseList* tcl);


	private:
	protected:
		int n;							        /**< number of products, also size of the vectors coeff, msbIn and lsbIn */
		vector<string> coeff;			  /**< the coefficients as strings */
		vector<mpfr_t> mpcoeff;			/**< the coefficients as MPFR numbers */
		vector<int> msbIn;			    /**< MSB weights of the inputs */
		vector<int> lsbIn;			    /**< LSB weights of the inputs */
		bool computeMsbOut;         /**< if true, initialize should compute msbOut. Otherwise the caller has provided it */
		int msbOut;							    /**< MSB weight of the output, may be computed out of the constants (depending on the constructor used) */
		int lsbOut;							    /**< LSB weight of the output */
		BitHeap* bitHeap;    			/**< The heap of weighted bits that will be used to do the additions */

	};


}

#endif
