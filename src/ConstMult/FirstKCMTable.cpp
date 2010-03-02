/*
  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License,
  2008-2010.
*/

#include <iostream>
#include <math.h>
#include <cstdlib>
#include "../utils.hpp"
#include "FirstKCMTable.hpp"
using namespace std;


namespace flopoco{

	// A table for the KCM Multiplication
	// the input of the table will be the number of inputs/LUT
	FirstKCMTable::FirstKCMTable(Target* target, int wIn, int wOut, mpz_class C) : 
		Table(target, wIn, wOut),  C_(C), wOut_(wOut)
	{
		ostringstream name; 
		name <<"KCMFirstTable_"<<wIn<<"_"<<wOut;
		setName(name.str());
	}

	FirstKCMTable::~FirstKCMTable() {}

	mpz_class FirstKCMTable::function(int x) {
		mpz_class result;
  
		result = mpz_class(x) * C_;

		if (result < 0 )
			result = (mpz_class(1)<<(wOut_)) + result; //this is how we pass negatives to 2's complement
		return  result;
	}

}
