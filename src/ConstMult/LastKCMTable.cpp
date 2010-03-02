/*
  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License, 2008-2010.
*/


#include <iostream>
#include <math.h>
#include <cstdlib>
#include "../utils.hpp"
#include "LastKCMTable.hpp"
using namespace std;

namespace flopoco{

	// A table for the KCM Multiplication
	// the input of the table will be the number of inputs/LUT
	LastKCMTable::LastKCMTable(Target* target, int wIn, int wOut, mpz_class C) : 
		Table(target, wIn, wOut),  C_(C), wIn_(wIn), wOut_(wOut)
	{
		ostringstream name; 
		name <<"KCMLastTable_"<<wIn<<"_"<<wOut;
		setName(name.str());
	}

	LastKCMTable::~LastKCMTable() {}

	mpz_class LastKCMTable::function(int x) {
		mpz_class result;
		//  int n = intlog2(x); // the number of bits of x;
  
		mpz_class maxPositive = ( mpz_class(1)  << (wIn_-1)) -1;
  
		mpz_class mx = x;
  
		if ( mx > maxPositive) //this means that x is negative;
			mx = mx - ( mpz_class(1) << wIn_ );
  
		//determine if x is positive or negative
    
		result = mpz_class(mx) * C_;
		if (result < 0 )
			result = (mpz_class(1)<<(wOut_)) + result; //this is how we pass negatives to 2's complement
		return  result;
	}

}
