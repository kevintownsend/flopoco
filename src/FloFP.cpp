#include "FloFP.hpp"
#include "utils.hpp"

FloFP::FloFP(int wE, int wF) : wE(wE), wF(wF)
{
	if (wE > 30)
		throw "FloFP::FloFP: Using exponents larger than 30 bits is not supported.";
}

FloFP::FloFP(int wE, int wF, mpfr_t m) : wE(wE), wF(wF)
{
	if (wE > 30)
		throw "FloFP::FloFP: Using exponents larger than 30 bits is not supported.";
	operator=(m);
}

mpz_class FloFP::getMantissaSignalValue() { return mantissa; }
mpz_class FloFP::getExceptionSignalValue() { return exception; }
mpz_class FloFP::getSignSignalValue() { return sign; }
mpz_class FloFP::getExponentSignalValue() { return exponent; }
mpz_class FloFP::getFractionSignalValue() { return mantissa + (mpz_class(1)<<wF); }

FloFP FloFP::operator*(FloFP fp)
{
	mpfr_t x, y, r;
	mpfr_init2(x, wF+1);
	mpfr_init2(y, fp.wF+1);
	mpfr_init2(r, wF+fp.wF+2);
	getMPFR(x);
	fp.getMPFR(y);
	mpfr_mul(r, x, y, GMP_RNDN);
	FloFP flofp(max(wE, fp.wE) + 1, wF + fp.wF + 2, r);
	mpfr_clears(r, x, y, 0);
	return flofp;
}

FloFP FloFP::operator*(mpfr_t mpX)
{
	FloFP fpX(24, mpfr_get_prec(mpX)+32, mpX);

	return *this * fpX;
}

void FloFP::getMPFR(mpfr_t mp)
{
	/* NaN */
	if (exception == 3)
	{
		mpfr_set_nan(mp);
		return;
	}

	/* Infinity */
	if (exception == 2)
	{
		mpfr_set_inf(mp, (sign == 1) ? -1 : 1);
		return;
	}

	/* Zero */
	if (exception == 0)
	{
		mpfr_set_d(mp, (sign == 1) ? -0.0 : +0.0, GMP_RNDN);
		return;
	}
	
	/* „Normal” numbers
	 * mp = (-1) * (1 + (mantissa / 2^wF)) * 2^unbiased_exp
	 * unbiased_exp = exp - (1<<(wE-1)) + 1
	 */
	mpfr_set_z(mp, mantissa.get_mpz_t(), GMP_RNDN);
	mpfr_div_2si(mp, mp, wF, GMP_RNDN);
	mpfr_add_ui(mp, mp, 1, GMP_RNDN);
	
	mp_exp_t exp = exponent.get_si();
	exp -= ((1<<(wE-1))-1);
	mpfr_mul_2si(mp, mp, exp, GMP_RNDN);

	/* Sign */
	if (sign == 1)
		mpfr_neg(mp, mp, GMP_RNDN);
}

FloFP& FloFP::operator=(mpfr_t mp_)
{
	mpfr_t mp;
	mpfr_init2(mp, mpfr_get_prec(mp_));
	mpfr_set(mp, mp_, GMP_RNDN);

	/* NaN */
	if (mpfr_nan_p(mp))
	{
		exception = 3;
		sign = 0;
		exponent = 0;
		mantissa = 0;
		return *this;
	}

	/* Inf */
	if (mpfr_inf_p(mp))
	{
		exception = 2;
		sign = mpfr_sgn(mp) > 0 ? 0 : 1;
		exponent = 0;
		mantissa = 0;
		return *this;
	}

	/* Zero */
	if (mpfr_zero_p(mp))
	{
		exception = 0;
		sign = mpfr_sgn(mp) > 0 ? 0 : 1;
		exponent = 0;
		mantissa = 0;
		return *this;
	}

	/* Normal numbers */
	exception = 1;
	sign = mpfr_sgn(mp) > 0 ? 0 : 1;
	mpfr_abs(mp, mp, GMP_RNDN);

	/* Get exponent
	 * mpfr_get_exp() return exponent for significant in [1/2,1)
	 * but we require [1,2). Hence the -1.
	 */
	mp_exp_t exp = mpfr_get_exp(mp)-1;

	/* Extract mantissa */
	mpfr_div_2si(mp, mp, exp, GMP_RNDN);
	mpfr_sub_ui(mp, mp, 1, GMP_RNDN);
	mpfr_mul_2si(mp, mp, wF, GMP_RNDN);
	mpfr_get_z(mantissa.get_mpz_t(), mp, GMP_RNDN);

	// Due to rounding, the mantissa might overflow (i.e. become bigger
	// then we expect).
	if (mantissa == mpz_class(1) << wF)
	{
		mantissa = 0;
		exp++;
	}

	if (mantissa >= mpz_class(1) << wF)
		throw std::string("Mantissa is to big after conversion to VHDL signal.");
	if (mantissa < 0)
		throw std::string("Mantissa is negative after conversion to VHDL signal.");

	/* Bias and store exponent */
	exp += ((1<<(wE-1))-1);
	exponent = exp;

	/* Handle underflow */
	if (exponent < 0)
	{
		exception = 0;
		exponent = 0;
		mantissa = 0;
	}

	/* Handle overflow */
	if (exponent >= (1<<(wE)))
	{
		exception = 2;
		exponent = 0;
		mantissa = 0;
	}

	mpfr_clear(mp);

	return *this;
}

FloFP& FloFP::operator=(mpz_class s)
{
	mantissa = s & ((mpz_class(1) << wF) - 1); s = s >> wF;
	exponent = s & ((mpz_class(1) << wE) - 1); s = s >> wE;
	sign = s & mpz_class(1); s = s >> 1;
	exception = s & mpz_class(3); s = s >> 2;

	if (s != 0)
		throw std::string("FloFP::operator= s is bigger than expected.");

	return *this;
}

mpz_class FloFP::getSignalValue()
{
	/* Sanity checks */
	if ((sign != 0) && (sign != 1))
		throw std::string("FloFP::getSignal: sign is invalid.");
	if ((exception < 0) || (exception > 3))
		throw std::string("FloFP::getSignal: exception is invalid.");
	if ((exponent < 0) || (exponent >= (1<<wE)))
		throw std::string("FloFP::getSignal: exponent is invalid.");
	if ((mantissa < 0) || (mantissa >= (mpz_class(1)<<wF)))
		throw std::string("FloFP::getSignal: mantissa is invalid.");
	return (((((exception << 1) + sign) << wE) + exponent) << wF) + mantissa;
}

FloFP& FloFP::operator=(FloFP fp)
{
	/* Pass this through MPFR to lose precision */
	mpfr_t mp;
	mpfr_init2(mp, wF+1);
	fp.getMPFR(mp);
	operator=(mp);
	mpfr_clear(mp);

	return *this;
}

