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

		
	// -------- Parameter set up -----------------

	addFPInput ("X", wE, wF);
	addFPOutput("R", wE, wF);
	
	// --------- Sub-components ------------------


}

FPSqrt::~FPSqrt() {
}






void FPSqrt::emulate(TestCase * tc)
{
}



