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
#include <sstream>
#include "utils.hpp"
#include "math.h"

using namespace std;

/** return a string representation of an mpz_class on a given number of bits */
string unsignedBinary(mpz_class x, int size){
	string s;
	mpz_class po2, number;
	char bit;

	if(x<0) {
		cerr<<"unsigned_binary: Positive number expected, got x="<<x.get_d()<<endl;
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

/** return the binary representation of a floating point number in the
		FPLibrary/FloPoCo format */
string fp2bin(mpfr_t x, int wE, int wF){
	mpfr_t mpx, one, two;
	ostringstream s;

	// copy the input
	mpfr_init2 (mpx, wF+1);
	mpfr_set (mpx, x, GMP_RNDN);


	// exception bits
	if(mpfr_nan_p (mpx)) {
		s << "11";
		for(int i=0; i<wE+wF+1; i++)
			s<< "0";
		return s.str();
	}

	if(mpfr_zero_p (mpx)) {
		s << "00";
		for(int i=0; i<wE+wF+1; i++)
			s<< "0";
		return s.str();
	}

	if(mpfr_inf_p (mpx)) {
		s << "10";
		for(int i=0; i<wE+wF+1; i++)
			s<< "0";
		return s.str();
	}

	// otherwise normal number
	s << "01";

	// sign bit
	int sign = mpfr_sgn(mpx);

	if(sign<0) {
		mpfr_neg(mpx, mpx, GMP_RNDN);
		s << "1";
	}
	else
		s << "0";

	// compute exponent and mantissa
	mpz_class exponent = 0;
	mpz_class biased_exponent;

	mpfr_init2(one, 2);
	mpfr_set_d(one, 1.0, GMP_RNDN);
	mpfr_init2(two, 2);
	mpfr_set_d(two, 2.0, GMP_RNDN);

	while(mpfr_less_p(mpx,one)) {
		mpfr_mul(mpx, mpx, two, GMP_RNDN);
		exponent --;
	}
	while(mpfr_greaterequal_p(mpx, two)) {
		mpfr_div(mpx, mpx, two, GMP_RNDN);
		exponent ++;
	}

	// add exponent bias
	biased_exponent = exponent + (mpz_class(1)<<(wE-1)) - 1;

	if ( biased_exponent<0 || biased_exponent>=(mpz_class(1)<<wE) )  {
		cerr << "Exponent out of range in fp2bin, exiting"<<endl;
		exit(EXIT_FAILURE);
	}

	// exponent
	s << unsignedBinary(biased_exponent, wE);
		
	// significand
	
	mpfr_sub(mpx, mpx, one, GMP_RNDN);
	for (int i=0; i<wF; i++) {
		mpfr_mul(mpx, mpx, two, GMP_RNDN);
		if(mpfr_greaterequal_p(mpx, one)) {
			s << "1";
			mpfr_sub(mpx, mpx, one, GMP_RNDN);
		}
		else
			s << "0";
	}

	mpfr_clear(mpx);
	mpfr_clear(one);
	mpfr_clear(two);
	return s.str();
}

/** Print out binary writing of an integer on "size" bits */
// TODO remove this function
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

void printBinNumGMP(ostream& o, mpz_class x, int size){
	mpz_class px;  
	if(x<0) {
		o<<"-";
		px=-x;
	}
	else {
		//		o<<" "; removed because VHDL " 110" (with a space) is bad
		px=x;
	}
	printBinPosNumGMP(o, px, size);
}

void printBinPosNumGMP(ostream& o, mpz_class x, int size)
{
	o << unsignedBinary(x,size);
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

int intlog(int base, double number)
{
double n1,n2,r;
	n1 = intlog2(number);
	n2 = intlog2(base);
	r = ceil(n1/n2);
	return int(r);
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

mpz_class maxExp(int wE){
	return mpz_class(1)<<(wE-1);
}

mpz_class minExp(int wE){
	return 1 - (mpz_class(1)<<(wE-1));
}

mpz_class bias(int wE){
	return (mpz_class(1)<<(wE-1)) - 1;
}

mpz_class getLargeRandom(int n)
{
	static gmp_randstate_t state;
	static bool init = false;
	if (init == false)
	{
		init = true;
		gmp_randinit_default(state);
	}
	mpz_class o;
	mpz_urandomb(o.get_mpz_t(), state, n);
	return o;
}

string zeroGenerator(int n, int margins){
	ostringstream left,full, right, zeros;
	int i;

	for (i=1; i<=n;i++)
		zeros<<"0";

	left<<"\""<<zeros.str();
	full<<left.str()<<"\"";
	right<<zeros.str()<<"\"";

	switch(margins){
		case -2: return zeros.str(); break;
		case -1: return left.str(); break;
		case  0: return full.str(); break;
		case  1: return right.str(); break;
		default: return full.str();
	}
}
