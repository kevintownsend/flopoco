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
	nDigit = int(ceil((wF+6)/2));

	addFPInput ("X", wE, wF);
	addFPInput ("Y", wE, wF);
	addFPOutput("R", wE, wF);
	
	// --------- Sub-components ------------------

	srt4step = new SRT4Step(target, wF);
	oplist.push_back(srt4step);
	stepDelay= srt4step->getPipelineDepth()+1;
	srtDelay=nDigit*stepDelay;

	mult3AdderDelay=1; // TODO replace with an IntAdder
	conversionAdderDelay=1;
	finalRoundAdderDelay=1;

	// -------- Pipeline setup--------------------
	
	setOperatorType();
	
	if(isSequential())		
		setPipelineDepth(nDigit+mult3AdderDelay+conversionAdderDelay+finalRoundAdderDelay);
	else
		setPipelineDepth(0);

	


	// Signals-------------------------------
	addDelaySignal("fX", wF+1, mult3AdderDelay);
	addDelaySignal("fY", wF+1, srtDelay + mult3AdderDelay);
	addDelaySignal("fYTimes3", wF+3, srtDelay);
	for(i=0; i<=nDigit-1; i++){
		ostringstream w;
		w << "w" << i;
		addDelaySignal(w.str(), wF+3, 1);
	}
	for(i=0; i<nDigit; i++){
		ostringstream  q, qP, qM;
		q << "q" << i;
		addDelaySignal(q.str(), 3, i);
		qP << "qP" << i;
		addSignal(qP.str(), 2);
		qM << "qM" << i;
		addSignal(qM.str(), 2);
	}
	addDelaySignal("qP", 2*nDigit,1);
	addDelaySignal("qM", 2*nDigit,1); 
	addDelaySignal("fR0", 2*nDigit, 1);
	addSignal("fR", wF+4); // significand plus ov bit on the left and round and sticky bit on the right  
	addSignal("fRn1", wF+2); // normd significand, without implicit 1, plus round and sticky bit on the right  
	// The computation is done early, to have less bits to pipeline. Same for sR
	addDelaySignal("eRn0", wE+2,     // two bits more on the left to detect overflow and underflow conditions
				 srtDelay + mult3AdderDelay + 1 + conversionAdderDelay); 
	addDelaySignal("eRn1", wE+2, 1); // idem
	addSignal("round");
	addSignal("expfrac", wE+wF+2); // two bits more to detect overflow and underflow
	addSignal("expfracR", wE+wF+2); // idem
	addDelaySignal("sR", 1, srtDelay + mult3AdderDelay + conversionAdderDelay + finalRoundAdderDelay);
	addSignal("exnXY ", 4);
	addDelaySignal("exnR0", 2, srtDelay + mult3AdderDelay + conversionAdderDelay + finalRoundAdderDelay);
	addSignal("exnR ", 2);
	addSignal("exnRfinal ", 2);
}

FPDiv::~FPDiv() {
}


void FPDiv::outputVHDL(std::ostream& o, std::string name) {
	
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;
	int i;

	licence(o,"Jeremie Detrey, Florent de Dinechin (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);	
	srt4step->outputVHDLComponent(o);
	outputVHDLSignalDeclarations(o);	  
	beginArchitecture(o);
	outputVHDLRegisters(o); o<<endl;

	o << tab << "fX <= \"1\" & X(" << wF-1 << " downto 0);" << endl;
	o << tab << "fY <= \"1\" & Y(" << wF-1 << " downto 0);" << endl;

	o << tab << "fYTimes3 <= (\"00\" & fY) + (\"0\" & fY & \"0\");" << endl; // TODO an IntAdder here
	o << tab << "w"<<nDigit-1<<" <=  \"00\" & " << delaySignal("fX", mult3AdderDelay) << ";" << endl;

	for(i=nDigit-1; i>=1; i--) {
		ostringstream wi;
		wi << "w" << i;
		//		o << tab << "-- SRT4 step "; 
		o << tab << "step" << i << ": " << srt4step->getOperatorName();
		if (isSequential()) 
		o << "  -- pipelineDepth="<< srt4step->getPipelineDepth();
		o << endl;
		o << tab << tab << "port map ( x       => " << delaySignal(wi.str(),1) << "," << endl;
		o << tab << tab << "           d       => " 
		  << delaySignal("fY", nDigit -i + mult3AdderDelay) << "," << endl;
		o << tab << tab << "           dtimes3 => " 
		  << delaySignal("fYTimes3",nDigit-i) << "," << endl;
		if(isSequential()) {
			o << tab << tab << "           clk  => clk, " << endl;
			o << tab << tab << "           rst  => rst, " << endl;
		}
		o << tab << tab << "              q  => q" << i << "," << endl;
		o << tab << tab << "              w  => w" << i-1 << "     );" <<endl;
	}
 
 
	
  	o << tab << "q0(2 downto 0) <= \"000\" when  " << delaySignal("w0", 1) << " = (" << wF+2 << " downto 0 => '0')" << endl;
	o << tab << "             else " << delaySignal("w0", 1) << "(" << wF+2 << ") & \"10\";" << endl;

	for(i=nDigit-1; i>=1; i--) {
		ostringstream qi;
		qi << "q" << i;
		o << tab << "qP" << i <<" <=      " << delaySignal(qi.str(), i) << "(1 downto 0);" << endl;
		o << tab << "qM" << i <<" <=      " << delaySignal(qi.str(), i) << "(2) & \"0\";" << endl;
	}

	o << tab << "qP0 <= q0(1 downto 0);" << endl;
	o << tab << "qM0 <= q0(2)  & \"0\";" << endl;

	o << tab << "qP <= qP" << nDigit-1;
	for (i=nDigit-2; i>=0; i--)
		o << " & qP" << i;
	o << ";" << endl;

	o << tab << "qM <= qM" << nDigit-1 << "(0)";
	for (i=nDigit-2; i>=0; i--)
		o << " & qM" << i;
	o << " & \"0\";" << endl;

	// TODO an IntAdder here
  	o << tab << "fR0 <= " << delaySignal("qP",1) << " - " << delaySignal("qM", 1) << ";" << endl;

	if (1 == (wF & 1) ) // odd wF
    	o << tab << "fR <= "<< delaySignal("fR0",1) << "(" << 2*nDigit-1 << " downto 1);  -- odd wF" << endl;
	else 
    	o << tab << "fR <= "<< delaySignal("fR0",1) << "(" << 2*nDigit-1 << " downto 3)  & (fR0(2) or fR0(1));  -- even wF, fixing the round bit" << endl;


	o << tab << "-- normalisation" << endl;
	o << tab << "with fR(" << wF+3 << ") select" << endl;

	o << tab << tab << "fRn1 <= fR(" << wF+2 << " downto 2) & (fR(1) or fR(0)) when '1'," << endl;
	o << tab << tab << "        fR(" << wF+1 << " downto 0)                    when others;" << endl;

	o << tab << "-- exponent difference" << endl;
	// The computation is done early, to have less bits to pipeline. Same for sR
	o << tab << "eRn0 <= (\"00\" & X(" << wE+wF-1 << " downto " << wF << ")) - (\"00\" & Y(" << wE+wF-1 << " downto " << wF<< "));" << endl;
	o << tab << "eRn1 <= "<< delaySignal("eRn0",  mult3AdderDelay + srtDelay + 1 + conversionAdderDelay) << " + (\"000\" & (" << wE-2 << " downto 1 => '1') & fR(" << wF+3 << ")); -- add back bias" << endl;

	o << tab << "-- final rounding" <<endl;
	o << tab << "round <= fRn1(1) and (fRn1(2) or fRn1(0)); -- fRn1(0) is the sticky bit" << endl;
	o << tab << "expfrac <= eRn1 & fRn1(" << wF+1 << " downto 2) ;" << endl;
	o << tab << "expfracR <= expfrac + ((" << wE+wF+1 << " downto 1 => '0') & round);" << endl;
	o << tab << "exnR <=      \"00\"  when expfracR(" << wE+wF+1 << ") = '1'   -- underflow" <<endl;
	o << tab << "        else \"10\"  when expfracR(" << wE+wF+1 << " downto " << wE+wF << ") =  \"01\" -- overflow" <<endl;
	o << tab << "        else \"01\";      -- 00, normal case" <<endl;

	o << tab << "-- exception handling" <<endl;
	// Computed at early stage so that less signals are pipelined through
	o << tab << "exnXY <= X(" << wE+wF+2 << " downto " << wE+wF+1  << ") & Y(" << wE+wF+2 << " downto " << wE+wF+1 << ");" <<endl;
	o << tab << "with exnXY select" <<endl;
	o << tab << tab << "exnR0 <= " << endl;
	o << tab << tab << tab << "\"01\"  when \"0101\",                   -- normal" <<endl;
	o << tab << tab << tab << "\"00\"  when \"0001\" | \"0010\" | \"0110\", -- zero" <<endl;
	o << tab << tab << tab << "\"10\"  when \"0100\" | \"1000\" | \"1001\", -- overflow" <<endl;
	o << tab << tab << tab << "\"11\"  when others;                   -- NaN" <<endl;

	o << tab << "with " << delaySignal("exnR0", getPipelineDepth()) << " select" <<endl;
	o << tab << tab << "exnRfinal <= " <<endl;
	o << tab << tab << tab << "exnR   when \"01\", -- normal" <<endl;
	o << tab << tab << tab << "exnR0  when others;" <<endl;
	o << tab << "sR <= X(" << wE+wF << ") xor Y(" << wE+wF<< ");" << endl;
	o << tab << "R <= exnRfinal & " 
	  << delaySignal("sR", getPipelineDepth()) 
	  << " & expfracR(" << wE+wF-1 << " downto 0);" <<endl;
	endArchitecture(o);
	checkDelays();
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

