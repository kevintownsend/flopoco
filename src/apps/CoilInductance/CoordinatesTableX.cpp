#include <iostream>
#include <math.h>
#include <cstdlib>
#include "../../utils.hpp"
#include "CoordinatesTableX.hpp"
#include <stdio.h>
#include <mpfr.h>


using namespace std;

CoordinatesTableX::CoordinatesTableX(Target* target, int wIn, int LSBI,int MSBI) : 
   Table(target, wIn, (MSBI-LSBI)), MSBI(MSBI), LSBI(LSBI), adrWidth(wIn) , wOutm(MSBI-LSBI)
   {
	ostringstream name; 
	name <<"CoordinatesX_"<<wIn<<"_"<<MSBI<<"_"<<LSBI;
	setName(name.str());
   }

CoordinatesTableX::~CoordinatesTableX() {}
	

int    CoordinatesTableX::double2input(double x){
  int result;
  cerr << "??? PolynomialTable::double2input not yet implemented ";
  exit(1);
  return result;
}


double CoordinatesTableX::input2double(int x) {
  double y;
  cerr << "??? PolynomialTable::double2input not yet implemented ";
  exit(1);
  return(y);
}

mpz_class CoordinatesTableX::double2output(double x){
  cerr << "??? PolynomialTable::double2input not yet implemented ";
  exit(1);
  return 0;
}

double CoordinatesTableX::output2double(mpz_class x) {
  double y;
  cerr << "??? PolynomialTable::double2input not yet implemented ";
  exit(1);
  
  return(y);
}


mpz_class CoordinatesTableX::function(int x)
{
	return mpz_class(0);
}