/*
 * Floating Point Divider for FloPoCo
 *
 * Author : Jeremie Detrey, Florent de Dinechin
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


// TODO Test even and odd significands
// use IntAdder for normalization and SRT4Step

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "FPSqrt.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPSqrt::FPSqrt(Target* target, int wE, int wF) :
	Operator(target), wE(wE), wF(wF) {

	ostringstream name;

	name<<"FPSqrt_"<<wE<<"_"<<wF;

	uniqueName_ = name.str(); 

	if(wF!=23)
		throw "Only wF=23 at the moment";

	// -------- Parameter set up -----------------

	addFPInput ("X", wE, wF);
	addFPOutput("R", wE, wF);

	int coeff_msb[3];
	coeff_msb[0]=1;
	coeff_msb[1]=-2;
	coeff_msb[2]=-9;

	int res_lsb=-26;
	int y_lsb=-wF;
	int y_msb=-8;

	
	vhdl << tab << declare("fX", wF) << " <= X" << range(wF-1, 0) << ";"  << endl; 
	vhdl << tab << declare("expX", wE) << " <= X" << range(wE+wF-1, wF) << ";"  << endl; 
	vhdl << tab << declare("OddExp") << " <= expX(0);"  << endl;  
	

	vhdl << tab << declare("address", 0-y_msb) << " <= X" << range(wF, wF+y_msb+1) << ";"  << endl; 
	vhdl << tab << declare("Y", 16) << " <= X" << range(wF+y_msb, wF+y_lsb) << ";"  << endl;
	
	// --------- Sub-components ------------------


}

FPSqrt::~FPSqrt() {
}






void FPSqrt::emulate(TestCase * tc)
{
}



