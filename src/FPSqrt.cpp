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

	correctRounding=false; //TODO set by the constructor
	useDSP=true; //TODO set by the constructor
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

	if(useDSP) {
		vhdl << tab << declare("fX", wF) << " <= X" << range(wF-1, 0) << ";"  << endl; 
		vhdl << tab << declare("expX", wE) << " <= X" << range(wE+wF-1, wF) << ";"  << endl; 
		vhdl << tab << declare("OddExp") << " <= expX(0);"  << endl;  
	

		vhdl << tab << declare("address", 0-y_msb) << " <= X" << range(wF, wF+y_msb+1) << ";"  << endl; 
		vhdl << tab << declare("Y", 16) << " <= X" << range(wF+y_msb, wF+y_lsb) << ";"  << endl;
	
		vhdl << tab << declare("exnsX", 3) << " <= X" << range(wE+wF+2, wE+wF) << ";"  << endl;



		vhdl << tab << "-- sign/exception handling" << endl;
		vhdl << tab << "with " << use("exnsX") << " select" <<endl
			  << tab << tab <<  declare("exnR", 3) << " <= " << "\"01\" when \"010\", -- positive, normal number (TODO underflow) " << endl
			  << tab << tab << use("exnsX") << range(2, 1) << " when \"001\" | \"000\" | \"100\", " << endl
			  << tab << tab << "\"11\" when others;"  << endl;

		vhdl << tab << "R <= exnR & '0' & computed; -- TODO" << endl; 
	}

	else {		// Digit-recurrence implementation recycled from FPLibrary

		vhdl << tab << declare("fX", wF) << " <= X" << range(wF-1, 0) << ";"  << endl; 
		vhdl << tab << declare("fA0", wF+2) << " <= \"1\" & fX when X(wF) = '0' else"
			  << tab << "\"01\" & fX;" << endl;

		vhdl << tab << declare(join("w",wF+3), wF+4) << " <= \"111\" & fX & \"0\" when X(" << wF << ") = '0' else" << endl
			  << tab << tab << "\"1101\" & fX;" << endl;
		vhdl << tab << declare(join("d",wF+3)) << " <= '0';" << endl;

		for(int i=wF+2; i>0; i--) {

		}
#if 0
  sqrt : FPSqrt_Sqrt
    generic map ( wF => wF )
    port map ( fA => fA0,
               fR => fRn0 );
  with fRn0(wF+3) select
    fRn1(wF+2 downto 0) <= fRn0(wF+3 downto 2) & (fRn0(1) or fRn0(0)) when '1',
                              fRn0(wF+2 downto 0)                        when others;

  eRn0 <= "00" & nA(wE+wF-1 downto wF+1);
  eRn1 <= eRn0 + ("000" & (wE-3 downto 0 => '1')) + nA(wF);

  round : FP_Round
    generic map ( wE => wE,
                  wF => wF )
    port map ( eA => eRn1,
               fA => fRn1,
               eR => eRn,
               fR => fRn );

  sRn <= nA(wE+wF);

  format : FP_Format
    generic map ( wE => wE,
                  wF => wF )
    port map ( sA => sRn,
               eA => eRn,
               fA => fRn,
               nR => nRn );

  xsA <= nA(wE+wF+2 downto wE+wF);

  with xsA select
    nR(wE+wF+2 downto wE+wF+1) <= nRn(wE+wF+2 downto wE+wF+1) when "010",
                                            xsA(2 downto 1)                       when "001" | "000" | "100",
                                            "11"                                  when others;
  nR(wE+wF downto 0) <= nRn(wE+wF downto 0);

#endif
	}


}

FPSqrt::~FPSqrt() {
}






void FPSqrt::emulate(TestCase * tc)
{
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("X");

	/* Compute correct value */
	FPNumber fpx(wE, wF);
	fpx = svX;
	mpfr_t x, r;
	mpfr_init2(x, 1+wF);
	mpfr_init2(r, 1+wF); 
	fpx.getMPFR(x);

	if(correctRounding) {
		mpfr_sqrt(r, x, GMP_RNDN);
		FPNumber  fpr(wE, wF, r);
		/* Set outputs */
		mpz_class svr= fpr.getSignalValue();
		tc->addExpectedOutput("R", svr);
	}
	else { // faithful rounding 
		mpfr_sqrt(r, x, GMP_RNDU);
		FPNumber  fpru(wE, wF, r);
		mpz_class svru = fpru.getSignalValue();
		tc->addExpectedOutput("R", svru);

		mpfr_sqrt(r, x, GMP_RNDD);
		FPNumber  fprd(wE, wF, r);
		mpz_class svrd = fprd.getSignalValue();
		/* Set outputs */
		tc->addExpectedOutput("R", svrd);
	}

	mpfr_clears(x, r, NULL);
}



