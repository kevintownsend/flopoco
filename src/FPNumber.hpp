#ifndef __FPNumber_HPP
#define __FPNumber_HPP

#include <gmpxx.h>
#include <mpfr.h>

/**
 * Flopoco internal Floating Point. Defines an
 * abstraction on which arithmetic operations can easily be applied
 * but at the same times can easily be converted to VHDL signals.
 * Used for TestBench generation.
 */
class FPNumber
{
public:
	/**
	 * Constructs a new FPNumber.
	 * @param wE the width of the exponent
	 * @param wF the width of the significant
	 */
	FPNumber(int wE, int wF, bool normalise = true);

	/**
	 * Constructs a new initialised FPNumber.
	 * @param wE the width of the exponent
	 * @param wF the width of the significant
	 * @param m the initial value.
	 */
	FPNumber(int wE, int wF, mpfr_t m, bool normalise = true);

	/**
	 * Retrieves the significant.
	 * @return Returns an mpz_class, representing the
	 * VHDL signal of the mantissa, without leading 1.
	 */
	mpz_class getMantissaSignalValue();

	/**
	 * Retrieves the fraction.
	 * @return An mpz_class, representing the VHDL
	 * signal of the mantissa, plus leading 1.
	 */
	mpz_class getFractionSignalValue();

	/**
	 * Retrieves the two exception bits.
	 * @return the two exception bits as VHDL signals.
	 */
	mpz_class getExceptionSignalValue();

	/**
	 * Retrives the sign.
	 * @return the sign as a VHDL signal.
	 */
	mpz_class getSignSignalValue();

	/**
	 * Retrieves the exponent.
	 * @return the exponent as a VHDL signal.
	 */
	mpz_class getExponentSignalValue();

	/**
	 * Multiplies two FPNumbers using MPFR.
	 * @return a FPNumber representing the result of the multiplication.
	 */
	FPNumber operator*(FPNumber);

	/**
	 * Divides two FPNumbers using MPFR.
	 * @return a FPNumber which is the correctly rounded quotient.
	 */
	FPNumber operator/(FPNumber);

	/**
	 * Adds two FPNumbers using MPFR.
	 * @return a FPNumber representing the result of the addition.
	 */
	FPNumber operator+(FPNumber);
	


	/**
	 * Converts the currently stored FPNumber to an mpfr_t
	 * @param[out] m a preinitialized mpfr_t where to store the floating point
	 * @param[in] withFakeZero use a very small number to represent zero, so that the sign is preserved
	 */
	void getMPFR(mpfr_t m, bool withFakeZero = true);

	/**
	 * Stores an mpfr_t as an internal representation of Flopoco.
	 * @param m the mpfr_t to store.
	 */
	FPNumber &operator=(mpfr_t m);

	/**
	 * Assignes a signal value. Converts the signal value to the
	 * relevant FPNumber fields.
	 * @param s the signal value to assign.
	 */
	FPNumber &operator=(mpz_class s);

	/**
	 * Retrieved the VHDL signal representation of this floating point.
	 * @return a VHDL signal stored as mpz_class.
	 */
	mpz_class getSignalValue();

	/**
	 * Equality operator. Everything does through MPFR to make sure
	 * correct rounding occurs.
	 */
	FPNumber &operator=(FPNumber fp);

	/**
	 * Returns the exponential of the current FPNumber.
	 * @return a FPNumber storing the exponential.
	 */
	FPNumber exp();

	/**
	 * Returns the natural logarithm of the current FPNumber.
	 * @return a FPNumber storing the natural logarithm.
	 */
	FPNumber log();

	/**
	 * Returns the whole signal rounded down.
	 */
	mpz_class getRoundedDownSignalValue();

	/**
	 * Returns the whole signal rounded up.
	 */
	mpz_class getRoundedUpSignalValue();

	/**
	 * Returns wE and wF.
	 * @param[out] wE exponent precision
	 * @param[out] wF fraction precision
	 */
	void getPrecision(int &wE, int &wF);

	/**
	 * Assigns a double.
	 */
	FPNumber &operator=(double x);

	/**
	 * Changes this FPNumber, so that it is decremented with
	 * one ULP. Postfix form.
	 */
	FPNumber &operator--(int);

	/**
	 * Changes this FPNumber, so that it is incremented with
	 * one ULP. Postfix form.
	 */
	FPNumber &operator++(int);

	/**
	 * Changes this FPNumber, so that it is decremented with
	 * the given number of ULPs.
	 * @param x how many ULPs to decrement it with.
	 */
	FPNumber &operator-=(int x);

	/**
	 * Changes this FPNumber, so that it is incremented with
	 * the given number of ULPs.
	 * @param x how many ULPs to increment it with.
	 */
	FPNumber &operator+=(int x);

private:
	/** The width of the exponent */
	int wE;

	/** The width of the significant (without leading zero) */
	int wF;

	/** The value of the sign */
	mpz_class sign;
	
	/** The value of the exception */
	mpz_class exception;

	/** The value of the exponent */
	mpz_class exponent;

	/** The value of the mantissa (without leading one) */
	mpz_class mantissa;

	/** Should we operate like a normalising or unnormalising
	 * FPMultiplier */
	bool normalise;

	/** Should we shift the fraction with one bit?
	 * See the comment at the beginig of FPNumber.cpp */
	bool mustAddLeadingZero;

	/** Should we store all values rounded down? */
	bool mustRoundDown;
};

#endif

