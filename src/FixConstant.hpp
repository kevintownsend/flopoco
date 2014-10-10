#ifndef __FIXCONSTANT_HPP
#define __FIXCONSTANT_HPP

#include <iostream>
#include <sstream>
#include <gmpxx.h>
#include <mpfr.h>

namespace flopoco{

	/**
	 * A class representing a fixed-point constant. 
	 The fixed-point format is represented as two bit positions: (MSB, LSB). 
	 These two positions are included: the size of the signal is wMSB-wLSB+1.
	 If signed, the sign bit, in position wMSB, has weight -2^wMSB. 

	 We can have a positive constant in a signed Fixconstant object, just as we can have leading zeroes.

	 By default the constructor has to provide MSB and LSB. 
	 TODO convention for zero in terms of LSB and MSB.


	 */
	class FixConstant
	{
	public:

		/** A constructor with explicit wLSB and wMSB. May perform rounding, or have MSB or LSB zeroes */
		FixConstant(const int wMSB, const int wLSB, const bool isSigned, const mpfr_t val);

#if 0
		/** A constructor that deduces wLSB and wMSB from the value of the mpfr passed 
				This is where we need a convention for zero */
		FixConstant(const bool isSigned, const mpfr_t val);
#endif

		~FixConstant();
		bool isSigned;        /**< if the constant is a two's complement one */
		int MSB;              /**< weight of the MSB*/
		int LSB;              /**< weight of the LSB*/
		int width;            /**< wMSB - wLSB +1*/
		mpfr_t fpValue;       /**< the value of this constant as an mpfr_t */

		mpz_class getBitVectorAsMPZ(); /**< return the constant as a bit vector in a positive MPZ.   */ 
		std::string getBitVector(int margins=0);     /**< return the textual version of the constant. See utils.hpp for margins  */ 
		bool isZero();        /**< as name suggests */
		void addRoundBit(int weight);   /**< updates the value to add a bit of a certain weight */
		void changeMSB(int newMSB);
		void changeLSB(int newMSB);
		std::string report();
		
	};

}

#endif

