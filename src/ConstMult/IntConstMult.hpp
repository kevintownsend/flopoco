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

  void recodeBooth();
  void printBoothCode();

  //  int computeLUTCost(ShiftAddDag);
  void buildMultBooth();
  void buildMultBoothTree();
  void showShiftAddDag();

  void optimizeLefevre(const vector<mpz_class>& constants);

  // Overloading the virtual functions of Operator
  void output_vhdl(std::ostream& o, std::string name);

  virtual TestCaseList generateStandardTestCases(int n);
  virtual TestCaseList generateRandomTestCases(int n);
};

#endif
