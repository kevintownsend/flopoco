#ifndef FIXCOMPLEXKCM_HPP
#define FIXCOMPLEXKCM_HPP
#include <string>
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"


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
			 * @param lsb_out : the required precision
			 * @param constant : constant by which the input is multiplied
			 */
			FixComplexKCM(
					Target* target, 
					int msb_in, 
					int lsb_in, 
					int lsb_out,
					string constant
				);


			// Below all the functions needed to test the operator
			/* the emulate function is used to simulate in software the operator
			   in order to compare this result with those outputed by the vhdl
			   opertator */
			void emulate(TestCase * tc);

			/* function used to create Standard testCase defined by the developper */
			void buildStandardTestCases(TestCaseList* tcl);
		
		private:
			int msb_in;
			int lsb_in;
			int lsb_out;
			string constant;

	};


}//namespace
#endif
