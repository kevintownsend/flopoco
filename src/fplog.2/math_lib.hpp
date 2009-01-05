#ifndef MATH_LIB_H
#define MATH_LIB_H

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <inttypes.h>

using namespace std;


void mpfr_shift_left(mpfr_t& x, int s);
void mpfr_shift_right(mpfr_t& x, int s);
void printBinNumGMP(ostream& o, mpz_class number, int size);
void printBinPosNumGMP(ostream& o, mpz_class number, int size);
void printBinNum(ostream& o, uint64_t x, int size);
double iround(double number, int bits);
double itrunc(double number, int bits);
double ifloor(double number, int bits);
double ifloor_strict(double number, int bits);
double powOf2(int power);
double negPowOf2(int minusPower);
int intlog2(double number);
int intpow2(int number);

inline double max(double x, double y) {return (x > y ? x : y);}
inline double min(double x, double y) {return (x < y ? x : y);}
inline int max(int x, int y) {return (x > y ? x : y);}
inline int min(int x, int y) {return (x < y ? x : y);}

#endif
