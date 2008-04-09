#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <inttypes.h>

using namespace std;


void mpfr_shift_left(mpfr_t& x, int s);
void mpfr_shift_right(mpfr_t& x, int s);
void printBinNumGMP(ostream& o, mpz_class number, int size);
string unsigned_binary(mpz_class x, int size);
void printBinPosNumGMP(ostream& o, mpz_class number, int size);
void printBinNum(ostream& o, uint64_t x, int size);
double iround(double number, int bits);
double itrunc(double number, int bits);
double ifloor(double number, int bits);
double ifloor_strict(double number, int bits);

//  2 ^ power
double intpow2(int power);

//  2 ^ -minusPower. Exact, no round
double invintpow2(int minusPower);

// How many bits does it take to write number ?
int intlog2(double number);
int intlog2(mpz_class number);

inline double max(double x, double y) {return (x > y ? x : y);}
inline double min(double x, double y) {return (x < y ? x : y);}
inline int max(int x, int y) {return (x > y ? x : y);}
inline int min(int x, int y) {return (x < y ? x : y);}

#endif
