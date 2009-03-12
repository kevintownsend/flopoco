/*
 * A toy example so support the FPL paper 


This is a collision detector: it inputs 3 FP coordinates X,Y and Z and
the square of a radius, R2, and computes a boolean predicate which is
true iff X^2 + Y^2 + Z^2 < R2

There are two versions, selectable by changing the value of
USE_FP_OPERATORS below. One combines existing FloPoCo floating-point
operators, and the other one is a specific datapath designed in the
FloPoCo philosophy. 

Both versions have its "grey area", situations where the predicate may
be wrong with respect to the true result. This happens when X^2+Y^2+Z^2 is very close to R2. 

The grey area of the combination of FP operators is about 2.5 units in the last place of R2. 
The pure FloPoCo version is slightly more accurate, with a grey area smaller than 1 ulp of R2.

It is also a lot smaller and faster, of course.
 
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
#include <fstream>
#include <sstream>
#include <math.h>	// for NaN
#include "Collision.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"


using namespace std;

extern vector<Operator *> oplist;


#define USE_FP_OPERATORS 1


Collision::Collision(Target* target, int wE, int wF)
	: Operator(target), wE(wE), wF(wF)
{
	setCopyrightString("F. de Dinechin (2009)");

	ostringstream o;
	o << "Collision_" << wE << "_" << wF;
	setName(o.str());

	addFPInput("X", wE, wF);
	addFPInput("Y", wE, wF);
	addFPInput("Z", wE, wF);
	addFPInput("R2", wE, wF);
	addOutput("P"); 

#if USE_FP_OPERATORS
	vhdl << tab << declare("R2pipe", wE+wF) << " <= R2(" << wE+wF-1 << " downto 0);"  << endl; 
	//////////////////////////////////////////////////////////////////:
	//            A version that assembles FP operators             //

	FPMultiplier* mult = new FPMultiplier(target, wE, wF, wE, wF, wE, wF, 1);
	oplist.push_back(mult);
	FPAdder* add =  new FPAdder(target, wE, wF, wE, wF, wE, wF);
	oplist.push_back(add);

	inPortMap (mult, "X", "X");
	inPortMap (mult, "Y", "X");
	outPortMap(mult, "R", "X2");
	vhdl << instance(mult, "multx");

	inPortMap (mult, "X", "Y");
	inPortMap (mult, "Y", "Y");
	outPortMap(mult, "R", "Y2");
	vhdl << instance(mult, "multy");

	inPortMap (mult, "X", "Z");
	inPortMap (mult, "Y", "Z");
	outPortMap(mult, "R", "Z2");
	vhdl << instance(mult, "multz");

	syncCycleFromSignal("Z2", false);
	nextCycle(); 

	inPortMap (add, "X", "X2");
	inPortMap (add, "Y", "Y2");
	outPortMap(add, "R", "X2PY2");
	vhdl << instance(add, "add1");

	syncCycleFromSignal("X2PY2", false);
	nextCycle(); 

	inPortMap (add, "X", "X2PY2");
	inPortMap (add, "Y", "Z2");
	outPortMap(add, "R", "X2PY2PZ2");
	vhdl << instance(add, "add2");

	syncCycleFromSignal("X2PY2PZ2", false);
	nextCycle(); 

	// TODO An IntAdder here
	vhdl << tab << declare("diff", wE+wF+1) << " <= ('0'& " << use("R2pipe") << ") - ('0'& " << use("X2PY2PZ2") << ");" << endl;
	vhdl << tab <<  "P <= " << use("diff") << "(" << wE+wF << ");"  << endl ;
	
#else

	// Error analysis
	// 3 ulps in the multiplier truncation
	// 2 half-ulps in the addition: we add one carry in that turns the two truncations into two roundings 
	// Total 5 ulps, or 3 bits

	// guard bits for a faithful result
	int g=4; // so we're safe 


	// extract the three biased exponents. 
	vhdl << declare("EX", wE) << " <=  X" << range(wE+wF-1, wF)  << ";" << endl;
	vhdl << declare("EY", wE) << " <=  Y" << range(wE+wF-1, wF) << ";" << endl;
	vhdl << declare("EZ", wE) << " <=  Z" << range(wE+wF-1, wF) << ";" << endl;

	// determine the max of the exponents
	vhdl << declare("DEXY", wE+1) << " <=   ('0' & EX) - ('0' & EY);" << endl;
	vhdl << declare("DEYZ", wE+1) << " <=   ('0' & EY) - ('0' & EZ);" << endl;
	vhdl << declare("DEXZ", wE+1) << " <=   ('0' & EX) - ('0' & EZ);" << endl;
	vhdl << declare("XltY") << " <=   DEXY("<< wE<<");" << endl;
	vhdl << declare("YltZ") << " <=   DEYZ("<< wE<<");" << endl;
	vhdl << declare("XltZ") << " <=   DEXZ("<< wE<<");" << endl;

	// rename the exponents  to A,B,C with A>=(B,C)
  	vhdl << tab << declare("EA", wE)  << " <= " << endl
		  << tab << tab << use("EZ") << " when (" << use("XltZ") << "='1') and (" << use("YltZ") << "='1')  else " << endl
		  << tab << tab << use("EY") << " when (" << use("XltY") << "='1') and (" << use("YltZ") << "='0')  else " << endl
		  << tab << tab << use("EX") << " when others; " << endl;
  	vhdl << tab << declare("EB", wE)  << " <= " << endl
		  << tab << tab << use("EX") << " when (" << use("XltZ") << "='1') and (" << use("YltZ") << "='1')  else " << endl
		  << tab << tab << use("EZ") << " when (" << use("XltY") << "='1') and (" << use("YltZ") << "='0')  else " << endl
		  << tab << tab << use("EY") << " when others; " << endl;
  	vhdl << tab << declare("EC", wE)  << " <= " << endl
		  << tab << tab << use("EY") << " when (" << use("XltZ") << "='1') and (" << use("YltZ") << "='1')  else " << endl
		  << tab << tab << use("EX") << " when (" << use("XltY") << "='1') and (" << use("YltZ") << "='0')  else " << endl
		  << tab << tab << use("EZ") << " when others; " << endl;

	nextCycle();

	// Now recompute our two shift values -- they were already computed at cycle 0 but it is cheaper this way.
	vhdl << declare("ShiftValB", wE) << " <=    " << use("EA") << " - " << use("EB") << "; -- positive result, no overflow " << endl;
	vhdl << declare("ShiftValC", wE) << " <=    " << use("EA") << " - " << use("EC") << "; -- positive result, no overflow" << endl;
	


	// Back to cycle 0 for the significand datapath
	setCycle(0);
	// Square the significands TODO an IntSquarer here some day
	IntMultiplier* mult = new IntMultiplier(target, 1+ wF, 1+ wF);
	oplist.push_back(mult);

	vhdl << tab << declare("mX", wF+1)  << " << '1' & X" << range(wF-1, 0) << "; " << endl;

	inPortMap (mult, "X", "mX");
	inPortMap (mult, "Y", "mX");
	outPortMap(mult, "R", "mX2");
	vhdl << instance(mult, "multx");
	
	vhdl << tab << declare("mY", wF+1)  << " <= '1' & Y" << range(wF-1, 0) << "; " << endl;

	inPortMap (mult, "X", "mY");
	inPortMap (mult, "Y", "mY");
	outPortMap(mult, "R", "mY2");
	vhdl << instance(mult, "multy");

	vhdl << tab << declare("mZ", wF+1)  << " <= '1' & Z" << range(wF-1, 0) << "; " << endl;
	
	inPortMap (mult, "X", "mZ");
	inPortMap (mult, "Y", "mZ");
	outPortMap(mult, "R", "mZ2");
	vhdl << instance(mult, "multz");

	syncCycleFromSignal("mZ2", false);

	// truncate the three results to wF+g+1
  	vhdl << tab << declare("X2t", wF+g+1)  << " <= " << use("mX2") << range(wF+g, 0) << "; " << endl;
  	vhdl << tab << declare("Y2t", wF+g+1)  << " <= " << use("mY2") << range(wF+g, 0) << "; " << endl;
  	vhdl << tab << declare("Z2t", wF+g+1)  << " <= " << use("mZ2") << range(wF+g, 0) << "; " << endl;

	nextCycle(); 

	
	// Now we have our three FP squares, we rename them to A,B,C with A>=(B,C) 
	// only 3 3-muxes

  	vhdl << tab << declare("MA", wF+g+1)  << " <= " << endl
		  << tab << tab << use("mZ2") << " when (" << use("XltZ") << "='1') and (" << use("YltZ") << "='1')  else " << endl
		  << tab << tab << use("mY2") << " when (" << use("XltY") << "='1') and (" << use("YltZ") << "='0')  else " << endl
		  << tab << tab << use("mX2") << " when others; " << endl;
  	vhdl << tab << declare("MB", wF+g+1)  << " <= " << endl
		  << tab << tab << use("mX2") << " when (" << use("XltZ") << "='1') and (" << use("YltZ") << "='1')  else " << endl
		  << tab << tab << use("mZ2") << " when (" << use("XltY") << "='1') and (" << use("YltZ") << "='0')  else " << endl
		  << tab << tab << use("mY2") << " when others; " << endl;
  	vhdl << tab << declare("MC", wF+g+1)  << " <= " << endl
		  << tab << tab << use("mY2") << " when (" << use("XltZ") << "='1') and (" << use("YltZ") << "='1')  else " << endl
		  << tab << tab << use("mX2") << " when (" << use("XltY") << "='1') and (" << use("YltZ") << "='0')  else " << endl
		  << tab << tab << use("mZ2") << " when others; " << endl;


	//Synchronize exponent and significand datapath
	setCycleFromSignal("MA", false); // useless here but harmless, too
	syncCycleFromSignal("ShiftValB", false);
	nextCycle();

	Shifter* rightShifter = new Shifter(target,wF+g+1, wF+g+1,Right); 
	oplist.push_back(rightShifter);

	inPortMap  (rightShifter, "X", "MB");
	inPortMap  (rightShifter, "S", "ShiftValB");
	outPortMap (rightShifter, "R","shiftedB");
	vhdl << instance(rightShifter, "ShifterForB");

	inPortMap  (rightShifter, "X", "MC");
	inPortMap  (rightShifter, "S", "ShiftValC");
	outPortMap (rightShifter, "R","shiftedC");
	vhdl << instance(rightShifter, "ShifterForC");

	// superbly ignore the bits that are shifted out
	syncCycleFromSignal("shiftedB", false);
	int shiftedB_size = getSignalByName("shiftedB")->width();
  	vhdl << tab << declare("alignedB", wF+g+1)  << " <= " << use("shiftedB") << range(shiftedB_size-1, shiftedB_size -(wF+g)) << "; " << endl;
  	vhdl << tab << declare("alignedC", wF+g+1)  << " <= " << use("shiftedC") << range(shiftedB_size-1, shiftedB_size -(wF+g)) << "; " << endl;

	nextCycle();

	IntAdder* fracAddFar = new IntAdder(target,wF+4);
	fracAddFar->changeName(getName()+"_fracAddFar");
	oplist.push_back(fracAddFar);
	inPortMap  (fracAddFar, "X", "fracXfar");
	inPortMap  (fracAddFar, "Y", "fracYfarXorOp");
	inPortMap  (fracAddFar, "Cin", "cInAddFar");
	outPortMap (fracAddFar, "R","fracResultfar0");
	vhdl << instance(fracAddFar, "fracAdderFar");


#endif

}	

Collision::~Collision()
{
}







void Collision::emulate(TestCase * tc)
{
#if 0
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("X");

	/* Compute correct value */
	FPNumber fpx(wE, wF);
	fpx = svX;

	mpfr_t x, ru,rd;
	mpfr_init2(x,  1+wF);
	mpfr_init2(ru, 1+wF);
	mpfr_init2(rd, 1+wF); 
	fpx.getMPFR(x); 
	mpfr_log(rd, x, GMP_RNDD);
	mpfr_log(ru, x, GMP_RNDU);
#if 0
	mpfr_out_str (stderr, 10, 30, x, GMP_RNDN); cerr << " ";
	mpfr_out_str (stderr, 10, 30, rd, GMP_RNDN); cerr << " ";
	mpfr_out_str (stderr, 10, 30, ru, GMP_RNDN); cerr << " ";
	cerr << endl;
#endif
	FPNumber  fprd(wE, wF, rd);
	FPNumber  fpru(wE, wF, ru);
	mpz_class svRD = fprd.getSignalValue();
	mpz_class svRU = fpru.getSignalValue();
	tc->addExpectedOutput("R", svRD);
	tc->addExpectedOutput("R", svRU);
	mpfr_clears(x, ru, rd, NULL);
#endif
}
 

