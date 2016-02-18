#ifndef FIXIIR_HPP
#define FIXIIR_HPP

#include "Operator.hpp"
#include "utils.hpp"
#include "BitHeap/BitHeap.hpp"

namespace flopoco{

	class FixIIR : public Operator {

	public:
		/** @brief Constructor ; you must use bitheap in case of negative coefficient*/
		FixIIR(Target* target, int lsbIn, int msbOut, int lsbOut, vector<string> coeffb, vector<string> coeffa, double H=0.0);

		/** @brief Destructor */
		~FixIIR();

		// Below all the functions needed to test the operator
		/**
		 * @brief the emulate function is used to simulate in software the operator
		 * in order to compare this result with those outputed by the vhdl opertator
		 */
		void emulate(TestCase * tc);

		/** @brief function used to create Standard testCase defined by the developper */
		void buildStandardTestCases(TestCaseList* tcl);

		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(Target *target , vector<string> &args);
		static void registerFactory();

	private:
		int lsbIn;					/**< weight of the LSB in the input, considered as a signed number in (-1,1) */
		int msbOut;					/**< weight of the MSB in the result */
		int lsbOut;					/**< weight of the LSB in the result */
		vector<string> coeffb;			/**< the b_i coefficients as strings */
		vector<string> coeffa;			/**< the a_i coefficients as strings */
		double H;						/**< Worst case peak gain of this filter */
		double Heps;						/**< Worst case peak gain of the error filter */
		int n;							/**< number of taps on the numerator */
		int m;							/**< number of taps on the denominator */

		BitHeap* bitHeapB;    			/**< The bit heap for the FIR part */
		BitHeap* bitHeapA;    			/**< The bit heap for the recursive part */

		int g;							/**< number of guard bits */
		int wO;							/**< width of the result */

	private:
		int hugePrec;
		// TODO All arrays or all mallocs below,

		mpfr_t mpcoeffb[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsignb[10000];			/**< the signs of the coefficients */
		mpfr_t mpcoeffa[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsigna[10000];			/**< the signs of the coefficients */
		double *coeffb_d;           /**< version of coeffb as C-style arrays of double, because WCPG needs it this way */
		double *coeffa_d;           /**< version of coeffa as C-style arrays of double, because WCPG needs it this way */

		mpz_class xHistory[10000]; // history of x used by emulate
		int currentIndexA;
		int currentIndexB;
		mpfr_t yHistory[10000]; // history of y (result) used by emulate


	};

}

#endif
