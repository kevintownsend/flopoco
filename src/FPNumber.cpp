#include "FPNumber.hpp"
#include "utils.hpp"

/* Exponent of 2, used to represent signed zero */
#define ZERO_EXPONENT -1000000000

// All this is broken if wf_out is different from wf_in
// Replace all the operators with an emulate() method of Operator.

/* SwW: Slipper when Wet
 * When the FPMultiplier is set not to normalise results,
 * it outputs significants and exponents which have little 
 * logic when viewed from a software side. Unfortunately,
 * in order to have a good test bench, we have to emulate
 * those behaviours in software. Whenever you see this tag
 * expect hard-to-understand or ilogic code.
 *
 * What happens is that the binary fraction comma is placed
 * after the second bit. Exponents are
 * simply added together. Sometimes the first bit
 * is 1 (the result is „overnormalised”), other times
 * it is 0 (the result is normalized). We will detect this
 * „condition” by seeing if the exponent of the normalized result
 * is higher than the sum of the exponents of the operands.
 *
 * We will do the following. We will always store numbers as
 * normalized numbers (1.mantissa), exponents like FPMultiplier
 * and when necessary (i.e. mustunnormalize==true)), we will do
 * rounding with one bit „faster” and shift the significant in
 * getFractionSignalValue().
 */

FPNumber::FPNumber(int wE, int wF, bool normalise)
	: wE(wE), wF(wF), normalise(normalise), mustAddLeadingZero(0), mustRoundDown(0)
{
	if (wE > 30)
		throw "FPNumber::FPNumber: Using exponents larger than 30 bits is not supported.";
}

FPNumber::FPNumber(int wE, int wF, mpfr_t m, bool normalise)
	: wE(wE), wF(wF), normalise(normalise), mustAddLeadingZero(0), mustRoundDown(0)
{
	if (wE > 30)
		throw "FPNumber::FPNumber: Using exponents larger than 30 bits is not supported.";
	operator=(m);
}

mpz_class FPNumber::getMantissaSignalValue()
{
	if (mustRoundDown)
		throw std::string("This FPNumber does not have the correct mantissa value.");
	return mantissa;
}

mpz_class FPNumber::getExceptionSignalValue() { return exception; }

mpz_class FPNumber::getSignSignalValue() { return sign; }
mpz_class FPNumber::getExponentSignalValue()
{
	if (mustRoundDown)
		throw std::string("This FPNumber does not have the correct exponent value.");
	return exponent;
}

mpz_class FPNumber::getFractionSignalValue()
{
	/* SwW: Add a leading zero */
	if (!normalise && mustAddLeadingZero)
		return (mantissa + (mpz_class(1)<<wF)) >> 1;

	return mantissa + (mpz_class(1)<<wF);
}




FPNumber FPNumber::operator*(FPNumber fp)
{
	mpfr_t x, y, r;
	mpfr_init2(x, wF+1);
	mpfr_init2(y, fp.wF+1);
	mpfr_init2(r,   wF + fp.wF + 2); // r will hold an exact product
	getMPFR(x);
	fp.getMPFR(y);
	mpfr_mul(r, x, y, GMP_RNDN);
	FPNumber flofp(max(wE, fp.wE) + 1, wF + fp.wF + 2, r);

	/* SwW: Detect the „condition” */
	if (mpfr_get_exp(x) + mpfr_get_exp(y) != mpfr_get_exp(r))
		flofp.mustAddLeadingZero = true;

	mpfr_clears(r, x, y, 0, NULL);
	return flofp;
}

FPNumber FPNumber::operator+(FPNumber fp)
{
	mpfr_t x, y, r;
	mpfr_init2(x, 1+wF);
	mpfr_init2(y, 1+fp.wF);
	mpfr_init2(r, wF+fp.wF+3); // FIXME double rounding here
	getMPFR(x);
	fp.getMPFR(y);
	mpfr_add(r, x, y, GMP_RNDN);
	FPNumber flofp(wE, wF, r);

	return flofp;
}

FPNumber FPNumber::operator/(FPNumber fp)
{
	mpfr_t x, y, r;
	mpfr_init2(x, 1+wF);
	mpfr_init2(y, 1+fp.wF);
	mpfr_init2(r, wF+fp.wF+3);
	getMPFR(x);
	fp.getMPFR(y);
	mpfr_div(r, x, y, GMP_RNDN);
	FPNumber flofp(wE, wF, r);

	return flofp;
}

void FPNumber::getMPFR(mpfr_t mp, bool withFakeZero)
{
	if (!normalise)
		throw "FPNumber::getMPFR: Non-normalised case not implemented.";

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
		if (withFakeZero)
		{
			/* MPFR does NOT have the concept of +/- zero due to its wide range of
			 * exponent values. As FloPoCo does have +/- zero (as a result of different
			 * underflows), we must somehow simulate it in MPFR. We will do this
			 * by storing zero as a really small number, what we won't ever encounder in
			 * FloPoCo */
			mpfr_set_d(mp, (sign == 1) ? -1 : +1, GMP_RNDN);
			mpfr_mul_2si(mp, mp, ZERO_EXPONENT, GMP_RNDN);
			return;
		}
		else
		{
			mpfr_set_d(mp, 0, GMP_RNDN);
			return;
		}
	}
	
	/* „Normal” numbers
	 * mp = (-1) * (1 + (mantissa / 2^wF)) * 2^unbiased_exp
	 * unbiased_exp = exp - (1<<(wE-1)) + 1
	 */
	mpfr_set_prec(mp, wF+2);
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

FPNumber& FPNumber::operator=(mpfr_t mp_)
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
		sign = 0;	/* MPFR does not have a sign for zero */
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
	if (!normalise && mustAddLeadingZero)
	{
		/* SwW: we need to round with one bit earlier */ 
		mpfr_mul_2si(mp, mp, wF-1, GMP_RNDN);
		mpfr_get_z(mantissa.get_mpz_t(), mp, GMP_RNDN);
		mantissa = mantissa << 1;
	}
	else
	{
		mpfr_mul_2si(mp, mp, wF, GMP_RNDN);
		mpfr_get_z(mantissa.get_mpz_t(), mp, mustRoundDown ? GMP_RNDD : GMP_RNDN);
	}

	/* SwW: exponent is smaller in the „normalised” case */
	if (!normalise && !mustAddLeadingZero)
		exp--;

	// Due to rounding, the mantissa might overflow (i.e. become bigger
	// then we expect).
	if (mantissa == mpz_class(1) << wF)
	{
		/* SwW: don't add leading zero when mantissa has overflown */
		if (!normalise && mustAddLeadingZero)
		{
			mustAddLeadingZero = false;
			exp--;
		}
		exp++;
		mantissa = 0;
	}

	if (mantissa >= mpz_class(1) << wF)
		throw std::string("Mantissa is too big after conversion to VHDL signal.");
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

FPNumber& FPNumber::operator=(mpz_class s)
{
	mantissa = s & ((mpz_class(1) << wF) - 1); s = s >> wF;
	exponent = s & ((mpz_class(1) << wE) - 1); s = s >> wE;
	sign = s & mpz_class(1); s = s >> 1;
	exception = s & mpz_class(3); s = s >> 2;

	if (s != 0)
		throw std::string("FPNumber::operator= s is bigger than expected.");

	return *this;
}

mpz_class FPNumber::getSignalValue()
{
	if (mustRoundDown)
		throw std::string("This FPNumber does not have the correct value.");

	/* Sanity checks */
	if ((sign != 0) && (sign != 1))
		throw std::string("FPNumber::getSignal: sign is invalid.");
	if ((exception < 0) || (exception > 3))
		throw std::string("FPNumber::getSignal: exception is invalid.");
	if ((exponent < 0) || (exponent >= (1<<wE)))
		throw std::string("FPNumber::getSignal: exponent is invalid.");
	if ((mantissa < 0) || (mantissa >= (mpz_class(1)<<wF)))
		throw std::string("FPNumber::getSignal: mantissa is invalid.");
	return (((((exception << 1) + sign) << wE) + exponent) << wF) + mantissa;
}

FPNumber& FPNumber::operator=(FPNumber fp)
{
	mustAddLeadingZero = fp.mustAddLeadingZero;
	mustRoundDown = fp.mustRoundDown;

	/* Pass this through MPFR to lose precision */
	mpfr_t mp;
	mpfr_init(mp);	// XXX: Precision set in operator=
	fp.getMPFR(mp);
	operator=(mp);
	mpfr_clear(mp);

	return *this;
}

FPNumber FPNumber::exp()
{
	/* Compute exponential using MPFR */
	mpfr_t mpR, mpX;
	mpfr_init2(mpR, wF+3);	// XXX: is this enough?
	mpfr_init(mpX);	// XXX: precision set in getMPFR()
	getMPFR(mpX);
	mpfr_exp(mpR, mpX, GMP_RNDD);
	
	/* Create FPNumber */
	FPNumber ret(wE, wF);
	ret.mustRoundDown = true;
	ret = mpR;

	/* Cleanup */
	mpfr_clears(mpX, mpR, 0, NULL);

	return ret;
}

FPNumber FPNumber::log()
{
	/* Compute logarithm using MPFR */
	mpfr_t mpR, mpX;
	mpfr_init2(mpR, wF+3);	// XXX: is this enough?
	mpfr_init(mpX);	// XXX: precision set in getMPFR()
	getMPFR(mpX, false);
	mpfr_log(mpR, mpX, GMP_RNDD);
	
	/* Create FPNumber */
	FPNumber ret(wE, wF);
	ret.mustRoundDown = true;
	ret = mpR;

	/* Cleanup */
	mpfr_clears(mpX, mpR, 0, NULL);

	return ret;
}

mpz_class FPNumber::getRoundedDownSignalValue()
{
	if (!mustRoundDown)
		throw std::string("Only correct value is stored.");
	return (((((exception << 1) + sign) << wE) + exponent) << wF) + mantissa;
}

mpz_class FPNumber::getRoundedUpSignalValue()
{
	if (!mustRoundDown)
		throw std::string("Only correct value is stored.");
	return (((((exception << 1) + sign) << wE) + exponent) << wF) + mantissa + 1;
}

void FPNumber::getPrecision(int &wE, int &wF)
{
	wE = this->wE;
	wF = this->wF;
}

FPNumber& FPNumber::operator=(double x)
{
	mpfr_t mp;
	mpfr_init2(mp, 53);
	mpfr_set_d(mp, x, GMP_RNDN);

	operator=(mp);

	mpfr_clear(mp);

	return *this;
}

FPNumber& FPNumber::operator--(int)
{
	/* Decrementation from special values is hard-coded */
	if (exception == 3)
	{
		if (sign == 0) /* before +NaN follows +inf */
			exception = 2;
		else /* before -NaN follows +NaN */
			sign = 0;
	}
	else if (exception == 2) /* inf */
	{
		if (sign == 0)	/* before +inf comes max */
		{
			exception = 1;
			exponent = (mpz_class(1) << wE) - 1;
			mantissa = (mpz_class(1) << wF) - 1;
		}
		else /* before -inf comes -NaN */
			exception = 3;
	}
	else if (exception == 0) /* zero */
	{
		if (sign == 0)	/* before +0 comes -0 */
			sign = 1;
		else /* before -0 comes -min */
		{
			exception = 1;
			exponent = 0;
			mantissa = 0;
		}
	}
	else /* other numbers */
	{
		/* Combine exception, exponent & mantissa, for nice trick */
		mpz_class tz = (((exception << wE) + exponent) << wF) + mantissa;
		if (sign == 0)
			tz--;
		else
			tz++;
		mantissa = tz & ((mpz_class(1) << wF) - 1); tz = tz >> wF;
		exponent = tz & ((mpz_class(1) << wE) - 1); tz = tz >> wE;
		exception = tz & mpz_class(3);
	}

	return *this;
}

FPNumber& FPNumber::operator++(int)
{
	/* Negate number */
	sign = 1 - sign;
	/* Decrement it */
	operator--(0);
	/* Negate it again */
	sign = 1 - sign;

	return *this;
}

FPNumber& FPNumber::operator-=(int x)
{
	for ( ; x < 0; x++)
		operator++(0);
	for ( ; x > 0; x--)
		operator--(0);
	return *this;
}

FPNumber& FPNumber::operator+=(int x)
{
	return operator-=(-x);
}

