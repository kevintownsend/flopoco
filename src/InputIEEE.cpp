/*
 * Conversion from IEEE-like compact floating-point format to FloPoCo format
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
#include "InputIEEE.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0

InputIEEE::InputIEEE(Target* target, int wEI, int wFI, int wEO, int wFO) :
	Operator(target), wEI(wEI), wFI(wFI), wEO(wEO), wFO(wFO) {

	ostringstream name;

	name<<"InputIEEE_"<<wEI<<"_"<<wFI<<"_to_"<<wEO<<"_"<<wFO;

	uniqueName_ = name.str(); 

	// -------- Parameter set up -----------------

	addFPInput ("X", wEI, wFI);
	addFPOutput("R", wEO, wFO);

	vhdl << tab << declare("expX", wEI) << "  <= X" << range(wEI+wFI-1, wFI) << ";" << endl;
	vhdl << tab << declare("fracX", wFI) << "  <= X" << range(wFI-1, 0) << ";" << endl;
	vhdl << tab << declare("sX") << "  <= X(" << wEI+wFI << ");" << endl;
	vhdl << tab << declare("expZero") << "  <= '1' when expX = " << rangeAssign(wEI-1,0, "'0'") << " else '0';" << endl;
	vhdl << tab << declare("expInfty") << "  <= '1' when expX = " << rangeAssign(wEI-1,0, "'1'") << " else '0';" << endl;
	vhdl << tab << declare("fracNotZero") << " <= '0' when fracX = " << rangeAssign(wFI-1,0, "'0'") << " else '1';" << endl;
	
	
	if(wEI==wEO){ // copy. Works for zero and infty since in such cases the exn bits prevail
		vhdl << tab << declare("expR", wEO) << " <= expX;" << endl;
		// FLUSH TO ZERO, TODO some day: implement subnormals
		vhdl << tab << declare("exnR",2) << " <= (expInfty) &  (fracNotZero or not expZero); " << endl;
		cout << "Warning: subnormal inputs must be flushed to zero" << endl;
		// TODO we may still represent some subnormals with our minimal exponent
		vhdl << tab << declare("fracR",wFO) << " <= fracX;" << endl;
	}
	else if (wEI<wEO) { // No overflow possible. Subnormal inputs need to be normalized
		cout << "Warning: subnormal inputs would be representable in the destination format,\n   but will be flushed to zero anyway (TODO)" << endl;
	}
	else { // wEI>wEO: some exponents will lead to over/underflow. Subnormals are flushed to zero
		cout << "Warning: subnormal inputs must be flushed to zero" << endl;
	}

	vhdl << tab << "R <= exnR & sX & expR & fracR; " << endl; 

}

InputIEEE::~InputIEEE() {
}






void InputIEEE::emulate(TestCase * tc)
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



