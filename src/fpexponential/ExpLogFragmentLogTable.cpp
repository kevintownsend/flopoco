#include <iostream>
#include <math.h>
#include <cstdlib>
#include "../utils.hpp"
#include "ExpLogFragmentLogTable.hpp"
using namespace std;


ExpLogFragmentLogTable::ExpLogFragmentLogTable(Target* target, Fragment* fragment) : 
	Table(target, 
			fragment->reallength, // wIn
			fragment->accuracy - fragment->start), //wOut
	fragment(fragment)
 {
	ostringstream name; 
	name <<"ExpLogFragmentLogTable_"<<fragment->accuracy<<"_"<<fragment->start<<"_"<<fragment->end;
	setName(name.str());
	setCopyrightString("X. Pujol (2007), C. Klein  (2008), F. de Dinechin (2009)");
}

ExpLogFragmentLogTable::~ExpLogFragmentLogTable() {}



mpz_class ExpLogFragmentLogTable::function(int x) {
	mpz_class r;
	return r;
}

