#ifndef INTCONSTMULT_HPP
#define INTCONSTMULT_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../Operator.hpp"
#include "ShiftAddOp.hpp"
#include "ShiftAddDag.hpp"

/**
	Integer constant multiplication.

	See also ShiftAddOp, ShiftAddDag
	ShiftAddOp defines a shift-and-add operation for IntConstMult.

	ShiftAddDag deals with the implementation of an IntConstMult as a
	vector of ShiftAddOp. It defines the intermediate variables with
	their bit sizes and provide methods for evaluating the cost of an
	implementation.

*/


class IntConstMult : public Operator
{
public:
	IntConstMult(Target* target, int xsize, mpz_class n);
	~IntConstMult();

	mpz_class n;
	int nsize;
	int xsize;
	int rsize;
	int* bits;
	int* BoothCode;
	int nonZeroInBoothCode;

	ShiftAddDag* implementation;


	// Overloading the virtual functions of Operator
	void outputVHDL(std::ostream& o, std::string name);

	TestIOMap getTestIOMap();
	void fillTestCase(mpz_class a[]);

private:
	int build_pipeline(ShiftAddOp* sao, double& delay);
	void recodeBooth();
	void printBoothCode();
	void buildMultBooth();
	void buildMultBoothTree();
	void showShiftAddDag();
	void optimizeLefevre(const vector<mpz_class>& constants);

};

#endif
