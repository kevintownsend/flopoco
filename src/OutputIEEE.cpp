/*
 * Conversion from  FloPoCo format to IEEE-like compact floating-point format
 *
 * Author : 
 * Florent de Dinechin
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
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "OutputIEEE.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0

OutputIEEE::OutputIEEE(Target* target, int wEI, int wFI, int wEO, int wFO) :
	Operator(target), wEI(wEI), wFI(wFI), wEO(wEO), wFO(wFO) {

	ostringstream name;

	name<<"OutputIEEE_"<<wEI<<"_"<<wFI<<"_to_"<<wEO<<"_"<<wFO;

	uniqueName_ = name.str(); 

	// -------- Parameter set up -----------------

	addFPInput ("X", wEI, wFI);
	addFPOutput("R", wEO, wFO);

	
	vhdl << tab << "R <= X; " << endl; 

}

OutputIEEE::~OutputIEEE() {
}






void OutputIEEE::emulate(TestCase * tc)
{
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("X");

	/* Compute correct value */
	FPNumber fpx(wEI, wFI);
	fpx = svX;
	mpfr_t x, r;
	mpfr_init2(x, 1+wFI);
	mpfr_init2(r, 1+wFO); 
	fpx.getMPFR(x);

	mpfr_set(r, x, GMP_RNDN); ///TODO probably not enough
	FPNumber  fpr(wEO, wFO, r);

	/* Set outputs */
	mpz_class svr= fpr.getSignalValue();
	tc->addExpectedOutput("R", svr);

	mpfr_clears(x, r, NULL);
}



