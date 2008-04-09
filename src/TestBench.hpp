#ifndef TESTBENCH_HPP
#define TESTBENCH_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "TestBench.hpp"

class TestBench : public Operator
{
public:
  // The operator to wrap
  Operator* op;

  // The test case list 
  vector<TestCase> test_case_list;

  TestBench(Target* target, Operator* op, int n);
  ~TestBench();
  void output_vhdl(ostream& o, string name);
};


#endif
