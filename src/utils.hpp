#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#ifdef _WIN32
  #include "pstdint.h"
#else
  #include <inttypes.h>
#endif


using namespace std;

/** Returns under the form of a string of given size, the unsigned binary represenation of an integer.
 * @param x the number to be represented in unsigned binary
 * @param size the size of the output string
 * @return the string binary representation of x
 */
string unsignedBinary(mpz_class x, int size);

/** Return the binary representation of a floating point number in the
 * FPLibrary/FloPoCo format
 * @param x the number to be represented
 * @param wE the width (in bits) of the result exponent
 * @param wF the width (in bits) of the result fraction
 */
string fp2bin(mpfr_t x, int wE, int wF);

/** Prints the binary representation of a integer on size bits
 * @param o the output stream
 * @param number [uint64_t] the number to be represented
 * @param size the number of bits of the output
 */
void printBinNum(ostream& o, uint64_t x, int size);

/** Prints the binary representation of a integer on size bits
 * @param o the output stream
 * @param number the number to be represented
 * @param size the number of bits of the output
 */
void printBinNumGMP(ostream& o, mpz_class number, int size);

/** Prints the binary representation of a positive integer on size bits
 * @param o the output stream
 * @param number the number to be represented
 * @param size the number of bits of the output
 */
void printBinPosNumGMP(ostream& o, mpz_class number, int size);

/** Function which rounds a FP on a given number of bits 
 * @param number the numer to be rounded
 * @param bits the number of bits of the result
 * @return the rounded number on "bits" bits 
 */
double iround(double number, int bits);

/** Function which truncates a FP on a given number of bits 
 * @param number the numer to be truncated
 * @param bits the number of bits of the result
 * @return the truncated number on "bits" bits 
 */
double itrunc(double number, int bits);

/** Function which returns floor a FP on a given number of bits 
 * @param number the numer to be floored
 * @param bits the number of bits of the result
 * @return the floored number on "bits" bits 
 */
double ifloor(double number, int bits);

/** Function which returns the maximum exponent on wE bits
 * @param wE the number of bits
 * @return the maximum exponent on wE bits
 */
mpz_class maxExp(int wE);

/** Function which returns the minimum exponent on wE bits
 * @param wE the number of bits
 * @return the minimum exponent on wE bits
 */
mpz_class minExp(int wE);

/** Function which returns the bias for an exponent on wE bits
 * @param wE the number of bits of the exponent
 * @return the bias corresponding to this exponent
 */
mpz_class bias(int wE);

/** 2 to the power function.
 * @param power the power at which 2 is raised
 * @return 2^power
 */
double intpow2(int power);

/** 2 ^ (- minusPower). Exact, no round
 * @param minusPower - the power at which 2 will be raised
 * @return 2 ^ (- minusPower)
 */
double invintpow2(int minusPower);

/** How many bits does it take to write number. 
 * @param number the number to be represented (floating point)
 * @return the number of bits
 */
int intlog2(double number);

/** computes a logarithm in a given base
 * @param base base of the logarithm
 * @param number the number to be represented (floating point)
 * @return the result is bits (ceil)
 */
int intlog(mpz_class base, mpz_class number);


/** How many bits does it take to write number. 
 * @param number the number to be represented (mpz_class)
 * @return the number of bits
 */
int intlog2(mpz_class number);



/** Maximum.
 * @param[double] x first number 
 * @param[double] y second number
 * @return maximum between x and y
 */
inline double max(double x, double y) {return (x > y ? x : y);}

/** Minimum.
 * @param[double] x first number 
 * @param[double] y second number
 * @return minimum between x and y
 */
inline double min(double x, double y) {return (x < y ? x : y);}

/** Maximum.
 * @param[int] x first number 
 * @param[int] y second number
 * @return maximum between x and y
 */
inline int max(int x, int y) {return (x > y ? x : y);}

/** Minimum.
 * @param[int] x first number 
 * @param[int] y second number
 * @return minimum between x and y
 */
inline int min(int x, int y) {return (x < y ? x : y);}

/** Minimum.
 * @param[int] x first number 
 * @param[int] y second number
 * @return maximum between x and y
 */
inline mpz_class max(mpz_class x, mpz_class y) {return (x > y ? x : y);}

/** Minimum.
 * @param[mpz_class] x first number 
 * @param[mpz_class] y second number
 * @return minimum between x and y
 */
inline mpz_class min(mpz_class x, mpz_class y) {return (x < y ? x : y);}

/**
 * Generate a very big random number.
 * Due to rereusage of a PRNG, this function might be suboptimal.
 * @param n bit-width of the target random number.
 * @return an mpz_class representing the random number.
 */
mpz_class getLargeRandom(int n);

/**
 * A zero generator method which takes as input two arguments and returns a string of zeros with quotes as stated by the second argurment
 * @param[in] n		    integer argument representing the number of zeros on the output string
 * @param[in] margins	integer argument determining the position of the quotes in the output string. The options are: -2= no quotes; -1=left quote; 0=both quotes 1=right quote
 * @return returns a string of zeros with the corresonding quotes given by margins
 **/
string zeroGenerator(int n, int margins);

/**
 * Turns an arbitrary string (e.g. Sollya expression or FP number) to
 * part of a valid VHDL identifier. May (and usually will) loose information.
 * Looks ugly.
 * May begin with a digit.
 * @param[in] expr		expression to convert
 **/
string vhdlize(string const & expr);


string vhdlize(double num);

string mpz2string(mpz_class x);


#endif
