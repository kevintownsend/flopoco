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
	nDigit = (wF+6) >> 1; 

	addFPInput ("X", wE, wF);
	addFPInput ("Y", wE, wF);
	addFPOutput("R", wE, wF);
	
	// --------- Sub-components ------------------

	srt4step = new SRT4Step(target, wF);
	oplist.push_back(srt4step);



	setOperatorType();
	
	vhdl << tab << declare("fX",wF+1) << " <= \"1\" & X(" << wF-1 << " downto 0);" << endl;
	vhdl << tab << declare("fY",wF+1) << " <= \"1\" & Y(" << wF-1 << " downto 0);" << endl;

	vhdl << tab << "-- exponent difference, sign and exception combination computed early, to have less bits to pipeline" << endl;
	 
	vhdl << tab << declare("expR0", wE+2) << " <= (\"00\" & X(" << wE+wF-1 << " downto " << wF << ")) - (\"00\" & Y(" << wE+wF-1 << " downto " << wF<< "));" << endl;
	vhdl << tab << declare("sR") << " <= X(" << wE+wF << ") xor Y(" << wE+wF<< ");" << endl;
	vhdl << tab << "-- early exception handling " <<endl;
	vhdl << tab << declare("exnXY",4) << " <= X(" << wE+wF+2 << " downto " << wE+wF+1  << ") & Y(" << wE+wF+2 << " downto " << wE+wF+1 << ");" <<endl;
	vhdl << tab << "with exnXY select" <<endl;
	vhdl << tab << tab << declare("exnR0", 2) << " <= " << endl;
	vhdl << tab << tab << tab << "\"01\"  when \"0101\",                   -- normal" <<endl;
	vhdl << tab << tab << tab << "\"00\"  when \"0001\" | \"0010\" | \"0110\", -- zero" <<endl;
	vhdl << tab << tab << tab << "\"10\"  when \"0100\" | \"1000\" | \"1001\", -- overflow" <<endl;
	vhdl << tab << tab << tab << "\"11\"  when others;                   -- NaN" <<endl;
	vhdl << tab << " -- compute 3Y" << endl;
	vhdl << tab << declare("fYTimes3",wF+3) << " <= (\"00\" & fY) + (\"0\" & fY & \"0\");" << endl; // TODO an IntAdder here
	
	ostringstream wInit;
	wInit << "w"<<nDigit-1;
	vhdl << tab << declare(wInit.str(), wF+3) <<" <=  \"00\" & fX;" << endl;

	nextCycle();/////////////////////////////////////////////////////////////

	for(i=nDigit-1; i>=1; i--) {
		ostringstream wi, qi, wim1, stepi;
		wi << "w" << i;
		qi << "q" << i;
		wim1 << "w" << i-1;
		stepi << "step" << i; // the instance name
	
#if 1
	
		inPortMap(srt4step, "x", wi.str());
		inPortMap(srt4step, "d", "fY");
		inPortMap(srt4step, "dtimes3", "fYTimes3");
		outPortMap(srt4step, "q", qi.str());
		outPortMap(srt4step, "w", wim1.str());
		vhdl << instance(srt4step, stepi.str());

#else
		//		vhdl << tab << "-- SRT4 step "; 
		vhdl << tab << "step" << i << ": " << srt4step->getOperatorName();
		if (isSequential()) 
		vhdl << "  -- pipelineDepth="<< srt4step->getPipelineDepth();
		vhdl << endl;
		vhdl << tab << tab << "port map ( x       => " << use(wi.str()) << "," << endl;
		vhdl << tab << tab << "           d       => " << use("fY") << "," << endl;
		vhdl << tab << tab << "           dtimes3 => " << use("fYTimes3") << "," << endl;
		if(isSequential()) {
			vhdl << tab << tab << "           clk  => clk, " << endl;
			vhdl << tab << tab << "           rst  => rst, " << endl;
		}
		vhdl << tab << tab << "              q  => " << declare(qi.str(),3) << "," << endl;
		vhdl << tab << tab << "              w  => " << declare(wim1.str(), wF+3) << "     );" <<endl;

#endif
		nextCycle();///////////////////////////////////////////////////////////////////////

	}
 
 
	
  	vhdl << tab << declare("q0",3) << "(2 downto 0) <= \"000\" when  " << use("w0") << " = (" << wF+2 << " downto 0 => '0')" << endl;
	vhdl << tab << "             else " << use("w0") << "(" << wF+2 << ") & \"10\";" << endl;

	for(i=nDigit-1; i>=1; i--) {
		ostringstream qi, qPi, qMi;
		qi << "q" << i;
		qPi << "qP" << i;
		qMi << "qM" << i;
		vhdl << tab << declare(qPi.str(), 2) <<" <=      " << use(qi.str()) << "(1 downto 0);" << endl;
		vhdl << tab << declare(qMi.str(), 2)<<" <=      " << use(qi.str()) << "(2) & \"0\";" << endl;
	}

	vhdl << tab << declare("qP0", 2) << " <= " << use("q0") << "(1 downto 0);" << endl;
	vhdl << tab << declare("qM0", 2) << " <= " << use("q0") << "(2)  & \"0\";" << endl;

	vhdl << tab << declare("qP", 2*nDigit) << " <= qP" << nDigit-1;
	for (i=nDigit-2; i>=0; i--)
		vhdl << " & qP" << i;
	vhdl << ";" << endl;

	vhdl << tab << declare("qM", 2*nDigit) << " <= qM" << nDigit-1 << "(0)";
	for (i=nDigit-2; i>=0; i--)
		vhdl << " & qM" << i;
	vhdl << " & \"0\";" << endl;


	// TODO an IntAdder here
  	vhdl << tab << declare("fR0", 2*nDigit) << " <= " << use("qP") << " - " << use("qM") << ";" << endl;

	nextCycle();///////////////////////////////////////////////////////////////////////
	

	vhdl << tab << declare("fR", wF+4) << " <= "; 
	if (1 == (wF & 1) ) // odd wF
    	vhdl << use("fR0") << "(" << 2*nDigit-1 << " downto 1);  -- odd wF" << endl;
	else 
    	vhdl << use("fR0") << "(" << 2*nDigit-1 << " downto 3)  & (" << use("fR0") << "(2) or " << use("fR0") << "(1));  -- even wF, fixing the round bit" << endl;


	vhdl << tab << "-- normalisation" << endl;
	vhdl << tab << "with " << use("fR") << "(" << wF+3 << ") select" << endl;

	vhdl << tab << tab << declare("fRn1", wF+2) << " <= " << use("fR") << "(" << wF+2 << " downto 2) & (" << use("fR") << "(1) or " << use("fR") << "(0)) when '1'," << endl;
	vhdl << tab << tab << "        " << use("fR") << "(" << wF+1 << " downto 0)                    when others;" << endl;

	vhdl << tab << declare("expR1", wE+2) << " <= "<< use("expR0") 
		  << " + (\"000\" & (" << wE-2 << " downto 1 => '1') & " << use("fR") <<"(" << wF+3 << ")); -- add back bias" << endl;



	vhdl << tab << declare("round") << " <= " << use("fRn1") << "(1) and (" << use("fRn1") << "(2) or " << use("fRn1") << "(0)); -- fRn1(0) is the sticky bit" << endl;

	nextCycle();///////////////////////////////////////////////////////////////////////
	vhdl << tab << "-- final rounding" <<endl;
	vhdl << tab <<  declare("expfrac", wE+wF+2) << " <= " 
		  << use("expR1") << " & " << use("fRn1") << "(" << wF+1 << " downto 2) ;" << endl;
	vhdl << tab << declare("expfracR", wE+wF+2) << " <= " 
		  << use("expfrac") <<" + ((" << wE+wF+1 << " downto 1 => '0') & " << use("round") << ");" << endl;
	vhdl << tab <<  declare("exnR", 2) << " <=      \"00\"  when " << use("expfracR") << "(" << wE+wF+1 << ") = '1'   -- underflow" <<endl;
	vhdl << tab << "        else \"10\"  when  " << use("expfracR") << "(" << wE+wF+1 << " downto " << wE+wF << ") =  \"01\" -- overflow" <<endl;
	vhdl << tab << "        else \"01\";      -- 00, normal case" <<endl;


	vhdl << tab << "with " << use("exnR0") << " select" <<endl;
	vhdl << tab << tab << declare("exnRfinal", 2) << " <= " <<endl;
	vhdl << tab << tab << tab << use("exnR") << "   when \"01\", -- normal" <<endl;
	vhdl << tab << tab << tab << use("exnR0") << "  when others;" <<endl;
	vhdl << tab << "R <= " << use("exnRfinal") << " & " << use("sR") << " & " 
		  << use("expfracR") << "(" << wE+wF-1 << " downto 0);" <<endl;

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




void FPDiv::emulate(TestCase * tc)
{
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");

	/* Compute correct value */
	FPNumber fpx(wE, wF), fpy(wE, wF);
	fpx = svX;
	fpy = svY;
	mpfr_t x, y, r;
	mpfr_init2(x, 1+wF);
	mpfr_init2(y, 1+wF);
	mpfr_init2(r, 1+wF); 
	fpx.getMPFR(x);
	fpy.getMPFR(y);
	mpfr_div(r, x, y, GMP_RNDN);
	FPNumber  fpr(wE, wF, r);

	/* Set outputs */
	mpz_class svR = fpr.getSignalValue();
	tc->addExpectedOutput("R", svR);
	mpfr_clears(x, y, r, NULL);
}


