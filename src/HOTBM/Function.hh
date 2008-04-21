#ifndef _FUNCTION_HH_
#define _FUNCTION_HH_

#include <string>
#include <iostream>

#include "Util.hh"
#include "sollya.h"

using namespace std;



class Function {
public:
  Function(string name_);
  virtual ~Function();

  string getName() const;
  double eval(double x) const;
  void eval(mpfr_t r, mpfr_t x) const;
  sollya_node_t getSollyaNode() const;

private:
  string name;
  sollya_node_t node;
};

#endif // _FUNCTION_HH_
