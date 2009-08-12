#include <iostream>
#include <math.h>
#include <cstdlib>
#include "../utils.hpp"
#include "ExpLogFragmentExpTable.hpp"
using namespace std;


ExpLogFragmentExpTable::ExpLogFragmentExpTable(Target* target, Fragment* fragment) : 
	Table(target, 
			fragment->reallength, // wIn
			fragment->accuracy - fragment->start), //wOut
	fragment(fragment)
 {
	ostringstream name; 
	name <<"ExpLogFragmentExpTable_"<<fragment->accuracy<<"_"<<fragment->start<<"_"<<fragment->end;
	setName(name.str());
	setCopyrightString("X. Pujol (2007), C. Klein  (2008), F. de Dinechin (2009)");

}

ExpLogFragmentExpTable::~ExpLogFragmentExpTable() {}



mpz_class ExpLogFragmentExpTable::function(int x) {
	mpfr_t value, result;
	mpz_t r;
	mpz_class rr;

	mpfr_init2(value, fragment->accuracy + 1);
	mpfr_init2(result, fragment->accuracy + (x >= 0)); // TODO check that

	//	if(negative) x--; // was  delta = negPowOf2(exp_accuracy); if (negative) x -= delta; but exp_accuracy is set to fragment->end

	// convert input to real:

	mpfr_set_si(value, x, GMP_RNDD);            // exact
	mpfr_mul_2si(value, value, -fragment->end, GMP_RNDN);


	//evaluate function
	mpfr_exp(result, value, GMP_RNDD);         // round down.  This seems the only difference with SimpleFragmentExpTable
	mpfr_sub(result, result, value, GMP_RNDD); // -x, should be exact
	mpfr_set_d(value, 1.0, GMP_RNDD);          // 
	mpfr_sub(result, result, value, GMP_RNDD); // -1

	// convert function value to mpz
	
	mpfr_mul_2si(result, result, fragment->accuracy, GMP_RNDN);
	mpz_init2(r,400); // 400 bits should be enough for anybody
	mpfr_get_z(r, result, GMP_RNDN);
	rr=mpz_class(r);

	// free memory
	mpfr_clear(value);
	mpfr_clear(result);
	mpz_clear(r);
	return  rr;
}

