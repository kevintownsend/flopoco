/*
 * A toy example so support the FPL paper
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

#if 1
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

	// Square the mantissas 
	IntMultiplier* mult = new IntMultiplier(target, 1+ wF, 1+ wF);
	oplist.push_back(mult);

	vhdl << tab << declare("mX", wF+1)  << " << '1' & X" << range(wF-1, 0) << "; " << endl;

	inPortMap (mult, "X", "mX");
	inPortMap (mult, "Y", "mX");
	outPortMap(mult, "R", "mX2");
	vhdl << instance(mult, "multx");

	vhdl << tab << declare("mY", wF+1)  << " << '1' & Y" << range(wF-1, 0) << "; " << endl;

	inPortMap (mult, "X", "mY");
	inPortMap (mult, "Y", "mY");
	outPortMap(mult, "R", "mY2");
	vhdl << instance(mult, "multy");

	vhdl << tab << declare("mZ", wF+1)  << " << '1' & Z" << range(wF-1, 0) << "; " << endl;
	
	inPortMap (mult, "X", "mZ");
	inPortMap (mult, "Y", "mZ");
	outPortMap(mult, "R", "mZ2");
	vhdl << instance(mult, "multz");

	syncCycleFromSignal("X2PY2PZ2", false);
	nextCycle(); 

	// Shift the three squares to align them to R2
	
	vhdl << declare("EX2", wE+1) << " <=  X" << range(wE+wF-1, wF) << " & '0';" << endl;
	vhdl << declare("EY2", wE+1) << " <=  Y" << range(wE+wF-1, wF) << " & '0';" << endl;
	vhdl << declare("EZ2", wE+1) << " <=  Z" << range(wE+wF-1, wF) << " & '0';" << endl;
	int bias_val= (1<<(wEX_-1))-1;
	vhdl << tab << declare("bias", wE+1) << " <= CONV_STD_LOGIC_VECTOR("<<bias_val<<","<<wE_+1<<");"<<endl; 
	vhdl << declare("ER2", wE+1) << " <= '0' & R2" << range(wE+wF-1, wF) << ";" << endl;
	vhdl << tab << 
			//substract the bias value from the exponents' sum
			o<<tab<<"Exponents_Sum_Post_Bias_Substraction <= Exponents_Sum_Pre_Bias_Substraction_d - bias;"<<endl;//wEX_ + 2   

	Shifter* 
	
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
 

