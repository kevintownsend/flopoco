#include "Function.hh"



Function::Function(string name_, string shortName_)
  : name(name_), shortName(shortName_)
{
}

Function::~Function()
{
}

string Function::getName() const
{
  return name;
}

string Function::getShortName() const
{
  return shortName;
}

void Function::mpEval(mpfr_t mpR, mpfr_t mpX, int n) const
{
  double x = mpfr_get_d(mpX, GMP_RNDN);
  double r = eval(x, n);
  mpfr_set_d(mpR, r, GMP_RNDN);
}
