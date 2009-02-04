#include <iostream>
#include "math_lib.hpp"
#include "math.h"

using namespace std;

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

