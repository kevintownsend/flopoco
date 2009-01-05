#include <iostream>
#include <math.h>
#include "math_lib.hpp"
#include <cstdlib>
#include "FirstInvTable.hpp"
#include "FirstLogTable.hpp"
using namespace std;


FirstLogTable::FirstLogTable(int wIn, int wOut, FirstInvTable* fit) : 
  Table(wIn, wOut), fit(fit)
 {
   minIn = 0;
   maxIn = (1<<wIn) -1;
   if (wIn!=fit->wIn) {
     cerr<< "FirstLogTable::FirstLogTable: Please use same wIn as FirstInvTable"<<endl;
     exit(1);
   }
 }

FirstLogTable::~FirstLogTable() {}

  

int    FirstLogTable::double2input(double x){
  int result;
  cerr << "??? FirstLogTable::double2input not yet implemented ";
  exit(1);
//   result = (int) floor(x*((double)(1<<(wIn-1))));
//   if( result < minIn || result > maxIn) {
//     cerr << "??? FirstLogTable::double2input:  Input "<< result <<" out of range ["<<minIn<<","<<maxIn<<"]";
//     exit(1);
//  }
  return result;
}

double FirstLogTable::input2double(int x) {
  return(fit->input2double(x));
}

mpz_class    FirstLogTable::double2output(double y){
  //Here y is between -0.5 and 0.5 strictly, whatever wIn.  Therefore
  //we multiply y by 2 before storing it, so the table actually holds
  //2*log(1/m)
  double z = floor(2*y*((double)(1<<(wOut-1)))); 

  // otherwise, signed arithmetic on wOut bits
  if(z>=0)
    return (mpz_class) z;
  else
    return (z + (double)(1<<wOut));
    
}

double FirstLogTable::output2double(mpz_class x) {
  cerr<<" FirstLogTable::output2double TODO"; exit(1);
  //  double y=((double)x) /  ((double)(1<<(wOut-1)));
  //return(y);
}


mpz_class FirstLogTable::function(int x)
{ 
  mpz_class result;
  double apprinv;
  mpfr_t i,l;
  mpz_t r;

  mpfr_init(i);
  mpfr_init(l);
  mpz_init2(r,400);
  apprinv = fit->output2double(fit->function(x));;
  // result = double2output(log(apprinv));
  mpfr_set_d(i, apprinv, GMP_RNDN);
  mpfr_log(l, i, GMP_RNDN);



#if 0
  if (x>>(wIn-1)) { // the log will be positive
    mpfr_shift_left(l, wOut); 
    mpfr_get_z(r, l, GMP_RNDD);
    result=mpz_class(r);
  }
  else { // the log will be negative -- code it in two's complement
    mpfr_neg(l, l, GMP_RNDN);
    mpfr_shift_left(l, wOut); 
    mpfr_get_z(r, l, GMP_RNDD);
    result=mpz_class(r);
  };

#else

  if (x>>(wIn-1)) { // the log will be negative  -- code it in two's complement
    mpfr_shift_left(l, wOut); 
    mpfr_get_z(r, l, GMP_RNDD);
    result=mpz_class(r);
    // code in two's complement
    mpz_class t = mpz_class(1) << wOut;
    result = t-result;
  }
  else { // the log will be positive
    mpfr_neg(l, l, GMP_RNDN);
    mpfr_shift_left(l, wOut); 
    mpfr_get_z(r, l, GMP_RNDD);
    result=mpz_class(r);
  };
#endif  


  //  cout << "x="<<x<<" apprinv="<<apprinv<<" logapprinv="<<log(apprinv)<<" result="<<result<<endl;
  mpfr_clear(i);
  mpfr_clear(l);
  mpz_clear(r);
  return  result;
}



//void FirstLogTable::check_accuracy() {}
