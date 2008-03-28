#ifndef _SIN_HH_
#define _SIN_HH_

#include "Function.hh"

using namespace std;



class Sin : public Function {
public:
  Sin();

  double eval(double x, int n = 0) const;
  void mpEval(mpfr_t mpR, mpfr_t mpX, int n = 0) const;
};

#endif // _SIN_HH_
