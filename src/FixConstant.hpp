#ifndef __FIXCONSTANT_HPP
#define __FIXCONSTANT_HPP

#include <iostream>
#include <sstream>
#include <gmpxx.h>
#include <mpfr.h>

namespace flopoco{

	/**
	 * A class representing a fixed-point constant. 
	 The fixed-point format is represented as two weights: (wMSB, wLSB). 
	 These two weights are included: the size of the signal is wMSB-wLSB+1.
	 If the constant is signed, the bit of weight wMSB has a negative weight. 

	 By default the constructor has to provide MSB and LSB. There will be a method to resize/minimize a constant at some point.
	 */
	class FixConstant
	{
	public:

		/** A constructor with explicit wLSB and wMSB. May perform rounding, or have MSB or LSB zeroes */
		FixConstant(const int wMSB, const int wLSB, const bool isSigned, const mpfr_t val);

#if 0
		/** A constructor that deduces wLSB and wMSB from the value of the mpfr passed */
		FixConstant(const bool isSigned, const mpfr_t val);
#endif

		~FixConstant();

		/** Reports if the constant should be interpreted as a signed constant		 */		
		bool isSigned();

		int getMSB(); 
		int getLSB();

		std::string getBitVector(int margins=0);     /**< return the textual version of the constant. See utils.hpp for margins  */ 

	private:
		bool isSignedFormat;        /**< if the constant is a two's complement one */
		int MSB;                   /**< weight of the MSB*/
		int LSB;                   /**< weight of the LSB*/
		int width;                   /**< wMSB - wLSB +1*/
		mpz_class intValue;         /**< the bit vector in an mpz_class */
		mpfr_t fpValue;             /**< the mpfr_t value */
	};

}

#endif

