#include "Sin.hh"



Sin::Sin()
  : Function(string("sin(x)"), string("sin"))
{
}

#define PIo4 0.785398163397448309615660845820

double Sin::eval(double x, int n) const
{
  return (n%4 < 2 ? 1 : -1) * (n%2 ? cos(x*PIo4) : sin(x*PIo4)) * pow(PIo4, n);
}

void Sin::mpEval(mpfr_t mpR, mpfr_t mpX, int n) const
{
  mpfr_t mpTmp;
  mpfr_init(mpTmp);
  mpfr_set_d(mpTmp, PIo4, GMP_RNDN);
  mpfr_mul(mpR, mpX, mpTmp, GMP_RNDN);
  if (n%2)
    mpfr_cos(mpR, mpR, GMP_RNDN);
  else
    mpfr_sin(mpR, mpR, GMP_RNDN);
  mpfr_pow_ui(mpTmp, mpTmp, n, GMP_RNDN);
  mpfr_mul(mpR, mpR, mpTmp, GMP_RNDN);
  mpfr_mul_si(mpR, mpR, n%4 < 2 ? 1 : -1, GMP_RNDN);
  mpfr_clear(mpTmp);
}
