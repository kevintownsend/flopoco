/*
 * utility functions for FloPoCo
 *
 * Author : Florent de Dinechin
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <iostream>
#include "utils.hpp"
#include "math.h"

using namespace std;

// Print out binary writing of an integer on "size" bits

void printBinNum(ostream& o, uint64_t x, int size)
{
  uint64_t po2 = ((uint64_t) 1)<<size; 
  char bit;

  if(size>=64){
    cerr << "\n printBinNum: size larger than 64" << endl;
    exit(1);    
  }

  if ((x<0) || (x >= po2) ) {
    cerr << "\n printBinNum: input out of range" << endl;
    exit(1);
  }
  for (int i = 0; i < size ; i++) {
    po2 = po2>>1;
    if (x >= po2) {
      bit = '1';
      x -= po2;
    }
    else {
      bit = '0';
    }
    o << bit;
  }
}


void mpfr_shift_left(mpfr_t& x, int s) {
  mpfr_t two;
  int i;

  mpfr_init2(two, 2);
  mpfr_set_ui(two, 2, GMP_RNDD);
  for(i=0; i<s; i++) {
    mpfr_mul(x, x, two, GMP_RNDD);
  }
  mpfr_clear(two);
}

void mpfr_shift_right(mpfr_t& x, int s) {
  mpfr_t two;
  int i;

  mpfr_init2(two, 2);
  mpfr_set_ui(two, 2, GMP_RNDD);
  for(i=0; i<s; i++) {
    mpfr_div(x, x, two, GMP_RNDD);
  }
  mpfr_clear(two);
}



void printBinNumGMP(ostream& o, mpz_class x, int size)
{
  mpz_class px;  
  if(x<0) {
    o<<"-";
    px=-x;
  }
  else {
    o<<" ";
    px=x;
  }
  printBinPosNumGMP(o, px, size);
}




string unsigned_binary(mpz_class x, int size){
  string s;
  mpz_class po2, number;
  char bit;

  if(x<0) {
    cerr<<"printBinPosNumGMP: Positive number expected, got x="<<x.get_d()<<endl;
    exit(EXIT_FAILURE);
  }
  po2 = ((mpz_class) 1)<<size;
  number=x;
    
  for (int i = 0; i < size ; i++) {
    po2 = po2>>1;
    if (number >= po2) {
      bit = '1';
      number -= po2;
    }
    else {
      bit = '0';
    }
    s +=  bit;
  }
  return s;
}

void printBinPosNumGMP(ostream& o, mpz_class x, int size)
{
#if 0 // to be removed some day
  mpz_class po2, number;
  char bit;

  if(x<0) {
    cerr<<"printBinPosNumGMP: Positive number expected, got x="<<x.get_d()<<endl;
    exit(1);
  }
  po2 = ((mpz_class) 1)<<size;
  number=x;
    
  for (int i = 0; i < size ; i++) {
    po2 = po2>>1;
    if (number >= po2) {
      bit = '1';
      number -= po2;
    }
    else {
      bit = '0';
    }
    o << bit;
  }
#else
  o << unsigned_binary(x,size);
#endif
}



double iround(double number, int bits)
{
  double shift = intpow2(bits);
  double x = number * shift;
  return floor(x + 0.5) / shift;
}

double ifloor(double number, int bits)
{
  double shift = intpow2(bits);
  return floor(number * shift) / shift;
}

//  2 ^ power
double intpow2(int power)
{
  double x = 1;
  for (int i = 0; i < power; i++)
    x *= 2;
  return x;
}

//  2 ^ -minusPower. Exact, no round
double invintpow2(int minusPower)
{
  double x = 1;
  for (int i = 0; i < minusPower; i++)
    x /= 2;
  return x;
}

// How many bits does it take to write number ?
int intlog2(double number)
{
  double po2 = 1.0; int result = 0;
  while (po2 <= number) {
    po2 *= 2;
    result++;
  }
  return result;
}

int intlog2(mpz_class number)
{
  mpz_class po2 = 1; 
  int result = 0;
  while (po2 <= number) {
    po2 *= 2;
    result++;
  }
  return result;
}

