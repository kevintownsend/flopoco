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

#include "FPDiv.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPDiv::FPDiv(Target* target, int wE, int wF) :
	Operator(target), wE(wE), wF(wF) {

	int i, j;
	ostringstream name, synch, synch2;

	setOperatorName();
	setOperatorType();
		
	//parameter set up
	nDigit = int(ceil((wF+6)/2));

	addFPInput ("X", wE, wF);
	addFPInput ("Y", wE, wF);
	addFPOutput("R", wE, wF);
	


	// -------- Pipeline setup--------------------




	if(isSequential())		
		setPipelineDepth(1);
	else
		setPipelineDepth(0);

	


	// Signals-------------------------------
  addSignal("fA0", wF+1);
  addSignal("fB0", wF+1);
  addSignal("eRn0", wE+1);
  addSignal("fRn0", wF+4);
  addSignal("eRn1", wE+1);
  addSignal("fRn1", wF+3);
  addSignal("sRn");
  addSignal("eRn ", wE+1);
  addSignal("fRn ", wF+1);
  addSignal("nRn ", wE+wF+3);
  addSignal("xAB", 4);

  addSignal("fB3", wF+3);
  for(i=0; i<=nDigit; i++){
	  // q0->q, qd -> qtimesd, xo->xx, w0-> ww
	  ostringstream w, sel, q, qtimesd, xx, ww;
	  w << "w" << i;
	  q << "q" << i;
	  sel << "sel" << i;
	  qtimesd << "q"<< i << "TimesD" ;

	  xx << "xx" << i;
	  ww << "ww" << i;
	  addSignal(w.str(), wF+4);
	  addSignal(sel.str(), 5);
	  addSignal(q.str(), 3);
	  addSignal(qtimesd.str(), wF+4);
	  addSignal(xx.str(), wF+4);
	  addSignal(ww.str(), wF+4); // Beware was wF+4 downto 1
  }
  addSignal("q", nDigit*3);
  addSignal("qP", 2*nDigit);
  addSignal("qM", 2*nDigit); // Beware was 2*nDigit downto 0

  addSignal("fR0", 2*nDigit);

	//	addDelaySignal("",wF+3, 1);
	//addSignal("",wE+wF+3);
					
	//		cout<<"signals pass";	
}

FPDiv::~FPDiv() {
}

void FPDiv::setOperatorName(){
	ostringstream name;
	/* The name has the format: FPDiv_wE_wF : 
 */
	name<<"FPDiv_"<<wE<<"_"<<wF; 
	uniqueName_ = name.str(); 
}

void FPDiv::outputVHDL(std::ostream& o, std::string name) {
  
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;
	int i;

	licence(o,"Jeremie Detrey, Florent de Dinechin (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);	
			
	outputVHDLSignalDeclarations(o);	  
	beginArchitecture(o);
	outputVHDLRegisters(o); o<<endl;

	cout << tab << "fA0 <= \"1\" & X(" << wF-1 << " downto 0);" << endl;
	cout << tab << "fB0 <= \"1\" & Y(" << wF-1 << " downto 0);" << endl;

	cout << tab << "-- significand division " << endl; 

	cout << tab << "fB3 <= (\"00\" & fB) + (\"0\" & fB & \"0\");" << endl; // TODO pipeline
	cout << tab << "w"<<nDigit<<" <=  \"00\" & fA0" << endl;

	for(i=nDigit-1; i>=1; i--) {
	cout << tab << "-- SRT4 step "; 
	ostringstream wi, wip1, sel, qi, qiTimesD, xx, ww;
	wi << "w" << i;
	wip1 << "w" << i+1;
	qi << "q" << i;
	sel << "sel" << i;
	qiTimesD << "q"<< i << "TimesD" ;
	xx << "xx" << i;
	ww << "ww" << i;


	cout << tab << sel << " <= " << wip1 << " & Y(" << wF-1 << ");" << endl; 
  	cout << tab << "with "<< sel << " select" << endl;
   cout << tab << qi << "<= " << endl;
	cout << tab << tab << "\"001\" when \"00010\" | \"00011\"," << endl;
	cout << tab << tab << "\"010\" when \"00100\" | \"00101\" | \"00111\"," << endl;
	cout << tab << tab << "\"011\" when \"00110\" | \"01000\" | \"01001\" | \"01010\" | \"01011\" | \"01101\" | \"01111\"," << endl;
	cout << tab << tab << "\"101\" when \"11000\" | \"10110\" | \"10111\" | \"10100\" | \"10101\" | \"10011\" | \"10001\"," << endl;
	cout << tab << tab << "\"110\" when \"11010\" | \"11011\" | \"11001\"," << endl;
	cout << tab << tab << "\"111\" when \"11100\" | \"11101\"," << endl;
	cout << tab << tab << "\"000\" when \"00000\" | \"00001\" | \"11110\" | \"11111\"," << endl;
	cout << tab << tab << "\"---\" when others;" << endl;
	cout << tab << "	with "<< qi << " select" << endl;
   cout << tab << tab << qiTimesD << " <= "<< endl ;
	cout << tab << tab << tab << "\"000\" & fBO            when \"001\" | \"111\"," << endl;
	cout << tab << tab << tab << "\"00\" & fB0 & \"0\"     when \"010\" | \"110\"," << endl;
	cout << tab << tab << tab << "\"0\" & d3             when \"011\" | \"101\"," << endl;
	cout << tab << tab << tab << "(wF+3 downto 0 => '0') when \"000\"," << endl;
	cout << tab << tab << tab << "(wF+3 downto 0 => '-') when others;" << endl;

	cout << tab << xx << " <= fA0 & \"0\"" << endl;
	cout << tab << "with " << qi << "(2) select" << endl;
   cout << tab << ww << " <= " << xx << " - " << qiTimesD << " when '0'," << endl;
	cout << tab << "       " << xx << " + "  << qiTimesD << " when others;" << endl;

	//   q <= q0; // TODO
  	cout << tab << wi << " <= " << ww << "(" << wF+2 << " downto 1) & \"0\"" << endl;
	
	} 
 
	cout << tab << " <= ;" << endl; 
	


	/*
  srt : for i in nDigit-1 downto 1 generate
    step : FPDiv_SRT4_Step
      generic map ( wF => wF )
      port map ( x  => w((i+1)*(wF+3)-1 downto i*(wF+3)),
                 d  => fB,
                 d3 => fB3,
                 q  => q((i+1)*3-1 downto i*3),
                 w  => w(i*(wF+3)-1 downto (i-1)*(wF+3)) );
    qP(2*i+1 downto 2*i)   <= q(i*3+1 downto i*3);
    qM(2*i+2 downto 2*i+1) <= q(i*3+2) & "0";
  end generate;

-----------------------------------------------------------------
  sel <= x(wF+2 downto wF-1) & d(wF-1);
  with sel select
    q0 <= "001" when "00010" | "00011",
          "010" when "00100" | "00101" | "00111",
          "011" when "00110" | "01000" | "01001" | "01010" | "01011" | "01101" | "01111",
          "101" when "11000" | "10110" | "10111" | "10100" | "10101" | "10011" | "10001",
          "110" when "11010" | "11011" | "11001",
          "111" when "11100" | "11101",
          "000" when "00000" | "00001" | "11110" | "11111",
          "---" when others;

  with q0 select
    qd <= "000" & d                 when "001" | "111",
          "00" & d & "0"            when "010" | "110",
          "0" & d3                  when "011" | "101",
          (wF+3 downto 0 => '0') when "000",
          (wF+3 downto 0 => '-') when others;

  x0 <= x & "0";
  with q0(2) select
    w0 <= x0 - qd when '0',
          x0 + qd when others;

  q <= q0;
  w <= w0(wF+2 downto 1) & "0";
--------------------------------------

  q(2 downto 0) <= "000" when w(wF+2 downto 0) = (wF+1 downto 0 => '0') else
                   w(wF+2) & "10";
  qP(1 downto 0) <= q(1 downto 0);
  qM(2 downto 1) <= q(2) & "0";

  fR0 <= qP - (qM(2*nDigit-1 downto 1) & "0");

  round_odd : if wF mod 2 = 1 generate
    fR <= fR0(2*nDigit-1 downto 1);
  end generate;
  round_even : if wF mod 2 = 0 generate
    fR <= fR0(2*nDigit-1 downto 3) & (fR0(2) or fR0(1));
  end generate;


  with fRn0(wF+3) select
    fRn1(wF+2 downto 0) <= fRn0(wF+3 downto 2) & (fRn0(1) or fRn0(0)) when '1',
                              fRn0(wF+2 downto 0)                        when others;

  eRn0 <= ("0" & nA(wE+wF-1 downto wF)) - ("0" & nB(wE+wF-1 downto wF));
  eRn1 <= eRn0 + ("00" & (wE-2 downto 1 => '1') & fRn0(wF+3));

  round : FP_Round
    generic map ( wE => wE,
                  wF => wF )
    port map ( eA => eRn1,
               fA => fRn1,
               eR => eRn,
               fR => fRn );

  sRn <= nA(wE+wF) xor nB(wE+wF);

  format : FP_Format
    generic map ( wE => wE,
                  wF => wF )
    port map ( sA => sRn,
               eA => eRn,
               fA => fRn,
               nR => nRn );

  xAB <= nA(wE+wF+2 downto wE+wF+1) & nB(wE+wF+2 downto wE+wF+1);

  with xAB select
    nR(wE+wF+2 downto wE+wF+1) <= nRn(wE+wF+2 downto wE+wF+1) when "0101",
                                            "00"                                  when "0001" | "0010" | "0110",
                                            "10"                                  when "0100" | "1000" | "1001",
                                            "11"                                  when others;

  nR(wE+wF downto 0) <= nRn(wE+wF downto 0);
	*/
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
	r = x+y; // TODO
	/* Set outputs */
	svR = r.getSignalValue();
	
}

