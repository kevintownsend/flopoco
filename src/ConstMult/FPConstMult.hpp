#ifndef FPCONSTMULT_HPP
#define FPCONSTMULT_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../Operator.hpp"
#include "IntConstMult.hpp"


namespace flopoco{

	class FPConstMult : public Operator
	{
	public:
		/** The generic constructor */
		FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, int cst_sgn, int cst_exp, mpz_class cst_sig);
		
		/** An empty constructor,  used by CRFPConstMult */
		FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out);
		
#ifdef HAVE_SOLLYA
		/** A constructor that parses an expression for the constant */
		FPConstMult(Target* target, int wE_in, int wF_in, int wE_out, int wF_out, int wF_C, string constant);
#endif //HAVE_SOLLYA
		~FPConstMult();

		int wE_in; 
		int wF_in; 
		int wE_out; 
		int wF_out; 
		int cst_sgn;
		int cst_exp_when_mantissa_1_2;
		int cst_exp_when_mantissa_int;
		int cst_width;
		mpz_class cst_sig;
		string cst_name;

		mpfr_t mpfr_cst; 
		mpfr_t mpfr_cst_sig; // between 1 and 2
		mpfr_t mpfr_xcut_sig; // between 1 and 2

		mpz_class xcut_sig_rd; // an int on wF_in+1 bits, which is mpfr_xcut_sig rounded down 

		bool mantissa_is_one; /**< is the mantissa equal to 1? */
		IntConstMult *icm;
		int icm_depth;

		/** The method that sets up all the attributes, called by the various constructors to avoid code duplication */
		void setup();

		/** The method that declares all the signals and sets up the pipeline.
			 It is called by the constructors of FPConstMult and CRFPConstMult to avoid code duplication */
		void buildVHDL();

		void emulate(TestCase *tc);

		/* The value of the constant multiplicand */
		mpfr_t mpY;

		void fillTestCase(mpz_class a[]);
	};

}
#endif
