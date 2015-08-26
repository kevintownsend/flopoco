#ifndef FIXCOMPLEXKCM_HPP
#define FIXCOMPLEXKCM_HPP
#include <string>
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"

#include "gmp.h"
#include "mpfr.h"

namespace flopoco {

	/**
	 * @brief FixComplexKCM : operator that compute the product between a 
	 * complex number and a complex constant.
	 */
	class FixComplexKCM : public Operator {
		public:

			/**
			 * @brief FixComplexKCM: constructor
			 * @param target : the target FPGA
			 * @param msb_in : the power of two of the weight associated to the
			 * 				   input msb
			 * @param lsb_in : the power of two of the weight associated to the
			 * 				   input lsb
			 * @param signedInput : whether the input is to be considered as
			 * 					    signed.
			 * 						The input width is always msb_in - lsb_in +1
			 * 						even when signedInput is set;
			 * @param lsb_out : the required precision
			 * @param constant : constant by which the input is multiplied
			 */
			FixComplexKCM(
					Target* target, 
					bool signedInput,
					int msb_in, 
					int lsb_in, 
					int lsb_out,
					string constant_re,
					string constant_im
				);

			virtual ~FixComplexKCM();



			// Below all the functions needed to test the operator
			/* the emulate function is used to simulate in software the operator
			   in order to compare this result with those outputed by the vhdl
			   opertator */
			void emulate(TestCase * tc);

			/* function used to create Standard testCase defined by the developper */
			void buildStandardTestCases(TestCaseList* tcl);
		
		private:
			void init();

			bool signedInput;
			int msb_in;
			int lsb_in;
			int lsb_out;
			int msbre_out;
			int msbim_out;
			int input_width;
			int outputre_width;
			int outputim_width;
			int constantReMsb;
			int constantImMsb;
			bool constantReNeg;
			bool constantImNeg;

			string constant_re;
			string constant_im;

			//Computed values for testBenchs
			mpfr_t mpfr_constant_re;
			mpfr_t mpfr_constant_im;

	};


}//namespace
#endif
