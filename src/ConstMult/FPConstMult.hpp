#ifndef FPCONSTMULT_HPP
#define FPCONSTMULT_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../Operator.hpp"
#include "IntConstMult.hpp"

class FPConstMult : public Operator
{
public:
	FPConstMult(Target* target, int wE_in, int wF_in,int wE_out, int wF_out, int cst_sgn, int cst_exp, mpz_class cst_sig);
	~FPConstMult();

	int wE_in; 
	int wF_in; 
	int wE_out; 
	int wF_out; 
	int cst_sgn;
	int cst_exp_when_mantissa_1_2;
	int cst_exp_when_mantissa_int;
	mpz_class cst_sig;
	string cst_name;

	mpfr_t mpfr_cst_sig; // between 1 and 2
	mpfr_t mpfr_xcut_sig; // between 1 and 2

	mpz_class xcut_sig_rd; // an int on wF_in+1 bits, which is mpfr_xcut_sig rounded down 

	IntConstMult *icm;

	// Overloading the virtual functions of Operator
	//  void output_vhdl_component(ostream& o, string name);
	void output_vhdl(ostream& o, string name);

	virtual TestCaseList generateStandardTestCases(int n);
	virtual TestCaseList generateRandomTestCases(int n);
};


#endif
