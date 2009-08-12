#include <iostream>
#include <math.h>
#include <cstdlib>
#include "../utils.hpp"
#include "SimpleFragmentTable.hpp"
using namespace std;


/** This table holds e^x-x-1*/

SimpleFragmentTable::SimpleFragmentTable(Target* target, Fragment* fragment) : 
	Table(target, 
			fragment->reallength, // wIn
			fragment->accuracy - 2*fragment->start -1),  //wOut: e^x-x-1 is smaller than x^2/2
	fragment(fragment)
 {
	ostringstream name; 
	name <<"SimpleFragmentTable_"<<fragment->accuracy<<"_"<<fragment->start<<"_"<<fragment->end;
	setName(name.str());
	setCopyrightString("X. Pujol (2007), C. Klein  (2008), F. de Dinechin (2009)");
}

SimpleFragmentTable::~SimpleFragmentTable() {}




mpz_class SimpleFragmentTable::function(int x) {
	mpfr_t value, result;
	mpz_t r;
	mpz_class rr;

	mpfr_init2(value, fragment->accuracy + 1);
	mpfr_init2(result, fragment->accuracy + (x >= 0)); // TODO check that

	// convert input to real:
	mpfr_set_si(value, x, GMP_RNDD);            // exact
	mpfr_mul_2si(value, value, -fragment->end, GMP_RNDN);


	//evaluate function
	mpfr_exp(result, value, GMP_RNDN);         // arrondi au plus proche
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





#if 0

int    SimpleFragmentTable::double2input(double x){
  int result;
  throw string("??? SimpleFragmentTable::double2input not yet implemented ");
  //  exit(1);
  //return result;
}


double SimpleFragmentTable::input2double(int x) {
throw}





mpz_class    SimpleFragmentTable::double2output(double y){
  mpz_class z; 
  return z; 
}



double SimpleFragmentTable::output2double(mpz_class x) {
  cerr<<" SimpleFragmentTable::output2double TODO"; exit(1);
  //double y=((long double)x) /  ((long double)(1<<outputPrecision));
  //return(y);
}


#endif
