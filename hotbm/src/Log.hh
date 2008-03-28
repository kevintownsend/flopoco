#ifndef _LOG_HH_
#define _LOG_HH_

#include "Function.hh"

using namespace std;



class Log : public Function {
public:
  Log();

  double eval(double x, int n = 0) const;
  void mpEval(mpfr_t mpR, mpfr_t mpX, int n = 0) const;
};

#endif // _LOG_HH_
