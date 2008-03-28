#ifndef _FUNCTION_HH_
#define _FUNCTION_HH_

#include <string>
#include <iostream>

#include "Util.hh"

using namespace std;



class Function {
public:
  Function(string name_, string shortName_);
  virtual ~Function();

  string getName() const;
  string getShortName() const;

  virtual double eval(double x, int n = 0) const = 0;
  virtual void mpEval(mpfr_t mpR, mpfr_t mpX, int n = 0) const;

private:
  string name, shortName;
};

#endif // _FUNCTION_HH_
