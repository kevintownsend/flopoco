#include "Log.hh"



Log::Log()
  : Function(string("log2(1+x)"), string("log"))
{
}

double Log::eval(double x, int n) const
{
  return !n ? log2(1+x) : ((n%2?1:-1) * fact(n-1) / (pow(1+x, n) * log(2)));
}

void Log::mpEval(mpfr_t mpR, mpfr_t mpX, int n) const
{
  mpfr_t mpXp1;

  mpfr_init(mpXp1);
  mpfr_add_ui(mpXp1, mpX, 1, GMP_RNDN);

  if (!n)
    mpfr_log2(mpR, mpXp1, GMP_RNDN);
  else {
    mpfr_t mpTmp;
    mpfr_init(mpTmp);
    mpfr_pow_ui(mpR, mpXp1, n, GMP_RNDN);
    mpfr_set_ui(mpTmp, 2, GMP_RNDN);
    mpfr_log(mpTmp, mpTmp, GMP_RNDN);
    mpfr_mul(mpR, mpR, mpTmp, GMP_RNDN);
    mpfr_fac_ui(mpTmp, n-1, GMP_RNDN);
    mpfr_div(mpR, mpTmp, mpR, GMP_RNDN);
    mpfr_mul_si(mpR, mpR, n%2 ? 1 : -1, GMP_RNDN);
    mpfr_clear(mpTmp);
  }

  mpfr_clear(mpXp1);
}
