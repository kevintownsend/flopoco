/*
 * A generic class for tables of values
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
#include <cstdlib>
#include "utils.hpp"
#include "Table.hpp"

using namespace std;




int Table::double2input(double x){
	throw string("Error, double2input is being used and has not been overriden");
}

double Table::input2double(int x) {
	throw string("Error, input2double is being used and has not been overriden");
}

mpz_class Table::double2output(double x){
	throw string("Error, double2output is being used and has not been overriden");
}

double Table::output2double(mpz_class x) {
	throw string("Error, output2double is being used and has not been overriden");
}

#if 0 // TODO some day
mpz_class Table::mpfr2output(mpfr_t x){
	throw string("Error, mpfr2output is being used and has not been overriden");
}

void Table::output2mpfr(mpz_class x, mpfr_t y) {
	throw string("Error, output2mpfr is being used and has not been overriden");
}
#endif



Table::Table(Target* target, int _wIn, int _wOut, int _minIn, int _maxIn) : 
	Operator(target),
	wIn(_wIn), wOut(_wOut), minIn(_minIn), maxIn(_maxIn)
	{
	setCopyrightString("Florent de Dinechin (2007)");

	// Set up the IO signals
	addInput ("X"  , wIn);
	addOutput ("Y"  , wOut);
	if(maxIn==-1) maxIn=(1<<wIn)-1;
	if(minIn<0) {
		cerr<<"ERROR in Table::Table, minIn<0\n";
		exit(EXIT_FAILURE);
	}
	if(maxIn>=(1<<wIn)) {
		cerr<<"ERROR in Table::Table, maxIn too large\n";
		exit(EXIT_FAILURE);
	}
	if((minIn==0) && (maxIn==(1<<wIn)-1)) 
		full=true;
	else
		full=false;
	if (wIn > 10)
		cerr << "WARNING : FloPoCo is building a table with " << wIn << " input bits, it will be large." << endl;
	}


// We have to define this method because the constructor of Table cannot use the (pure virtual) function()
void Table::outputVHDL(std::ostream& o, std::string name) {
	int i,x;
	mpz_class y;

	vhdl	<< "  with X select  Y <= " << endl;
	for (x = minIn; x <= maxIn; x++) {
		y=function(x);
		//    cout << x <<"  "<< y << endl;
		vhdl 	<< tab << "\"" << unsignedBinary(y, wOut) << "\" when \"" << unsignedBinary(x, wIn) << "\"," << endl;
	}
	vhdl << tab << "\"";
	for (i = 0; i < wOut; i++) 
		vhdl << "-";
	vhdl <<  "\" when others;" << endl;

	Operator::outputVHDL(o,  name);
}


int Table::size_in_LUTs() {
	return wOut*int(intpow2(wIn-target_->lutInputs()));
}
