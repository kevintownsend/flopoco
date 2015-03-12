#ifndef FIXIIR_HPP
#define FIXIIR_HPP

#include "Operator.hpp"
#include "utils.hpp"
#include "BitHeap/BitHeap.hpp"

namespace flopoco{ 

	class FixIIR : public Operator {

	public:
		/* Constructor ; you must use bitheap in case of negative coefficient*/
		FixIIR(Target* target, int msbOut_, int lsbOut_, double H_, vector<string> coeffb_, vector<string> coeffa_, map<string, double> inputDelays = emptyDelayMap); 

		/* Destructor */
		~FixIIR();

		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
		  in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		/* function used to create Standard testCase defined by the developper */
		void buildStandardTestCases(TestCaseList* tcl);

  	private:
		int p;							/**< The precision of inputs and outputs */ 
		int n;							/**< number of taps */
		int m;							/**< number of taps */
		vector<string> coeffa;			/**< the coefficients as strings */
		vector<string> coeffb;			/**< the coefficients as strings */

		BitHeap* bitHeapA;    			/**< The heap of weighted bits that will be used to do the additions */
		BitHeap* bitHeapB;    			/**< The heap of weighted bits that will be used to do the additions */

		int g;							/**< number of guard bits */
		int wO;							/**< width of the result */
		int msbOut;					/**< weight of the MSB in the result */
		double H;						/**< Worst case peak gain */

	private:
		int hugePrec;

		mpfr_t mpcoeffb[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsignb[10000];			/**< the signs of the coefficients */

		mpfr_t mpcoeffa[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsigna[10000];			/**< the signs of the coefficients */

		mpz_class xHistory[10000]; // history of x used by emulate
		int currentIndexA;
		int currentIndexB;
		mpfr_t yHistory[10000]; // history of y (result) used by emulate


	};

}

#endif
