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

#include "FPDiv.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPDiv::FPDiv(Target* target, int wE, int wF) :
	Operator(target), wE(wE), wF(wF) {

	int i, j;
	ostringstream name;

	name<<"FPDiv_"<<wE<<"_"<<wF; 
	uniqueName_ = name.str(); 

		
	// -------- Parameter set up -----------------
	nDigit = int(ceil( double(wF+6)/2));

	addFPInput ("X", wE, wF);
	addFPInput ("Y", wE, wF);
	addFPOutput("R", wE, wF);
	
	// --------- Sub-components ------------------

	srt4step = new SRT4Step(target, wF);
	oplist.push_back(srt4step);



	setOperatorType();
	
	vhdl << tab << lhs("fX",wF+1) << " <= \"1\" & X(" << wF-1 << " downto 0);" << endl;
	vhdl << tab << lhs("fY",wF+1) << " <= \"1\" & Y(" << wF-1 << " downto 0);" << endl;

	vhdl << tab << "-- exponent difference, sign and exception combination computed early, to have less bits to pipeline" << endl;
	 
	vhdl << tab << lhs("expR0", wE+2) << " <= (\"00\" & X(" << wE+wF-1 << " downto " << wF << ")) - (\"00\" & Y(" << wE+wF-1 << " downto " << wF<< "));" << endl;
	vhdl << tab << lhs("sR") << " <= X(" << wE+wF << ") xor Y(" << wE+wF<< ");" << endl;
	vhdl << tab << "-- early exception handling " <<endl;
	vhdl << tab << lhs("exnXY",4) << " <= X(" << wE+wF+2 << " downto " << wE+wF+1  << ") & Y(" << wE+wF+2 << " downto " << wE+wF+1 << ");" <<endl;
	vhdl << tab << "with exnXY select" <<endl;
	vhdl << tab << tab << lhs("exnR0", 2) << " <= " << endl;
	vhdl << tab << tab << tab << "\"01\"  when \"0101\",                   -- normal" <<endl;
	vhdl << tab << tab << tab << "\"00\"  when \"0001\" | \"0010\" | \"0110\", -- zero" <<endl;
	vhdl << tab << tab << tab << "\"10\"  when \"0100\" | \"1000\" | \"1001\", -- overflow" <<endl;
	vhdl << tab << tab << tab << "\"11\"  when others;                   -- NaN" <<endl;
	vhdl << tab << " -- compute 3Y" << endl;
	vhdl << tab << lhs("fYTimes3",wF+3) << " <= (\"00\" & fY) + (\"0\" & fY & \"0\");" << endl; // TODO an IntAdder here
	
	ostringstream wInit;
	wInit << "w"<<nDigit-1;
	vhdl << tab << lhs(wInit.str(), wF+3) <<" <=  \"00\" & fX;" << endl;

	nextCycle();/////////////////////////////////////////////////////////////

	for(i=nDigit-1; i>=1; i--) {
		ostringstream wi, qi, wim1;
		wi << "w" << i;
		qi << "q" << i;
		wim1 << "w" << i-1;
		//		vhdl << tab << "-- SRT4 step "; 
		vhdl << tab << "step" << i << ": " << srt4step->getOperatorName();
		if (isSequential()) 
		vhdl << "  -- pipelineDepth="<< srt4step->getPipelineDepth();
		vhdl << endl;
		vhdl << tab << tab << "port map ( x       => " << rhs(wi.str()) << "," << endl;
		vhdl << tab << tab << "           d       => " << rhs("fY") << "," << endl;
		vhdl << tab << tab << "           dtimes3 => " << rhs("fYTimes3") << "," << endl;
		if(isSequential()) {
			vhdl << tab << tab << "           clk  => clk, " << endl;
			vhdl << tab << tab << "           rst  => rst, " << endl;
		}
		vhdl << tab << tab << "              q  => " << lhs(qi.str(),3) << "," << endl;
		vhdl << tab << tab << "              w  => " << lhs(wim1.str(), wF+3) << "     );" <<endl;
		nextCycle();///////////////////////////////////////////////////////////////////////

	}
 
 
	
  	vhdl << tab << lhs("q0",3) << "(2 downto 0) <= \"000\" when  " << rhs("w0") << " = (" << wF+2 << " downto 0 => '0')" << endl;
	vhdl << tab << "             else " << rhs("w0") << "(" << wF+2 << ") & \"10\";" << endl;

	for(i=nDigit-1; i>=1; i--) {
		ostringstream qi, qPi, qMi;
		qi << "q" << i;
		qPi << "qP" << i;
		qMi << "qM" << i;
		vhdl << tab << lhs(qPi.str(), 2) <<" <=      " << rhs(qi.str()) << "(1 downto 0);" << endl;
		vhdl << tab << lhs(qMi.str(), 2)<<" <=      " << rhs(qi.str()) << "(2) & \"0\";" << endl;
	}

	vhdl << tab << lhs("qP0", 2) << " <= " << rhs("q0") << "(1 downto 0);" << endl;
	vhdl << tab << lhs("qM0", 2) << " <= " << rhs("q0") << "(2)  & \"0\";" << endl;

	vhdl << tab << lhs("qP", 2*nDigit) << " <= qP" << nDigit-1;
	for (i=nDigit-2; i>=0; i--)
		vhdl << " & qP" << i;
	vhdl << ";" << endl;

	vhdl << tab << lhs("qM", 2*nDigit) << " <= qM" << nDigit-1 << "(0)";
	for (i=nDigit-2; i>=0; i--)
		vhdl << " & qM" << i;
	vhdl << " & \"0\";" << endl;


	// TODO an IntAdder here
  	vhdl << tab << lhs("fR0", 2*nDigit) << " <= " << rhs("qP") << " - " << rhs("qM") << ";" << endl;

	nextCycle();///////////////////////////////////////////////////////////////////////
	
	vhdl << tab << lhs("fR", wF+4) << " <= "; 
	if (1 == (wF & 1) ) // odd wF
    	vhdl << rhs("fR0") << "(" << 2*nDigit-1 << " downto 1);  -- odd wF" << endl;
	else 
    	vhdl << rhs("fR0") << "(" << 2*nDigit-1 << " downto 3)  & (fR0(2) or fR0(1));  -- even wF, fixing the round bit" << endl;


	vhdl << tab << "-- normalisation" << endl;
	vhdl << tab << "with " << rhs("fR") << "(" << wF+3 << ") select" << endl;

	vhdl << tab << tab << lhs("fRn1", wF+2) << " <= " << rhs("fR") << "(" << wF+2 << " downto 2) & (" << rhs("fR") << "(1) or " << rhs("fR") << "(0)) when '1'," << endl;
	vhdl << tab << tab << "        " << rhs("fR") << "(" << wF+1 << " downto 0)                    when others;" << endl;

	vhdl << tab << lhs("expR1", wE+2) << " <= "<< rhs("expR0") 
		  << " + (\"000\" & (" << wE-2 << " downto 1 => '1') & " << rhs("fR") <<"(" << wF+3 << ")); -- add back bias" << endl;



	vhdl << tab << lhs("round") << " <= " << rhs("fRn1") << "(1) and (" << rhs("fRn1") << "(2) or " << rhs("fRn1") << "(0)); -- fRn1(0) is the sticky bit" << endl;

	nextCycle();///////////////////////////////////////////////////////////////////////
	vhdl << tab << "-- final rounding" <<endl;
	vhdl << tab <<  lhs("expfrac", wE+wF+2) << " <= " 
		  << rhs("expR1") << " & " << rhs("fRn1") << "(" << wF+1 << " downto 2) ;" << endl;
	vhdl << tab << lhs("expfracR", wE+wF+2) << " <= " 
		  << rhs("expfrac") <<" + ((" << wE+wF+1 << " downto 1 => '0') & " << rhs("round") << ");" << endl;
	vhdl << tab <<  lhs("exnR", 2) << " <=      \"00\"  when " << rhs("expfracR") << "(" << wE+wF+1 << ") = '1'   -- underflow" <<endl;
	vhdl << tab << "        else \"10\"  when  " << rhs("expfracR") << "(" << wE+wF+1 << " downto " << wE+wF << ") =  \"01\" -- overflow" <<endl;
	vhdl << tab << "        else \"01\";      -- 00, normal case" <<endl;


	vhdl << tab << "with " << rhs("exnR0") << " select" <<endl;
	vhdl << tab << tab << lhs("exnRfinal", 2) << " <= " <<endl;
	vhdl << tab << tab << tab << rhs("exnR") << "   when \"01\", -- normal" <<endl;
	vhdl << tab << tab << tab << rhs("exnR0") << "  when others;" <<endl;
	vhdl << tab << "R <= " << rhs("exnRfinal") << " & " << rhs("sR") << " & " 
		  << rhs("expfracR") << "(" << wE+wF-1 << " downto 0);" <<endl;

}

FPDiv::~FPDiv() {
}


void FPDiv::outputVHDL(std::ostream& o, std::string name) {
	licence(o,"Jeremie Detrey, Florent de Dinechin (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);	
	srt4step->outputVHDLComponent(o);
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);
	o << buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);
}




TestIOMap FPDiv::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("X"));
	tim.add(*getSignalByName("Y"));
	tim.add(*getSignalByName("R"));
	return tim;
}

void FPDiv::fillTestCase(mpz_class a[])
{
	/* Get I/O values */
	mpz_class& svX = a[0];
	mpz_class& svY = a[1];
	mpz_class& svR = a[2];

	/* Compute correct value */
	FPNumber x(wE, wF), y(wE, wF), r(wE, wF);
	x = svX;
	y = svY;

	r = x/y; // TODO
	/* Set outputs */
	svR = r.getSignalValue();
	
}

