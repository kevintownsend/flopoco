#ifndef __FLOFP_HPP
#define __FLOFP_HPP

#include <gmpxx.h>
#include <mpfr.h>

/**
 * Flopoco internal Floating Point. Defines an
 * abstraction on which arithmetic operations can easily be applied
 * but at the same times can easily be converted to VHDL signals.
 * Used for TestBench generation.
 */
class FloFP
{
public:
	/**
	 * Constructs a new FloFP.
	 * @param wE the width of the exponent
	 * @param wF the width of the significant
	 */
	FloFP(int wE, int wF);

	/**
	 * Constructs a new initialised FloFP.
	 * @param wE the width of the exponent
	 * @param wF the width of the significant
	 * @param m the initial value.
	 */
	FloFP(int wE, int wF, mpfr_t m);

	/**
	 * Retrieves the significant.
	 * @return Returns an mpz_class, representing the
	 * VHDL signal of the mantissa, without leading 1.
	 */
	mpz_class getMantissa();

	/**
	 * Retrieves the two exception bits.
	 * @return the two exception bits as VHDL signals.
	 */
	mpz_class getException();

	/**
	 * Retrives the sign.
	 * @return the sign as a VHDL signal.
	 */
	mpz_class getSign();

	/**
	 * Retrieves the exponent.
	 * @return the exponent as a VHDL signal.
	 */
	mpz_class getExponent();

	/**
	 * Multiplies to FloFP using MPFR.
	 * @return a FloFP representing the result of the multiplication.
	 */
	FloFP operator*(FloFP);

	/**
	 * Converts the currently stored FloFP to an mpfr_t
	 * @param[out] m a preinitialized mpfr_t where to store the floating point
	 */
	void getMPFR(mpfr_t m);

	/**
	 * Stores an mpfr_t as an internal representation of Flopoco.
	 * @param m the mpfr_t to store.
	 */
	FloFP &operator=(mpfr_t m);

	/**
	 * Assignes a signal value. Converts the signal value to the
	 * relevant FloFP fields.
	 * @param s the signal value to assign.
	 */
	FloFP &operator=(mpz_class s);

	/**
	 * Retrieved the VHDL signal representation of this floating point.
	 * @return a VHDL signal stored as mpz_class.
	 */
	mpz_class getSignal();

	/**
	 * Equality operator. Everything does through MPFR to make sure
	 * correct rounding occurs.
	 */
	FloFP &operator=(FloFP fp);

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
};

#endif

