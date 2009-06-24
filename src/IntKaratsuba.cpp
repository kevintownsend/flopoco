/*
 * An integer multiplier for FloPoCo using the Karatsuba method
 *
 * Authors : Bogdan Pasca
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
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"
#include "IntKaratsuba.hpp"

using namespace std;
extern vector<Operator*> oplist;

IntKaratsuba:: IntKaratsuba(Target* target, int wIn) :
	Operator(target), wIn_(wIn), wOut_(2*wIn){

	ostringstream name;
	name << "IntKaratsuba_" << wIn_<<"_f"<<target->frequencyMHz();
	setName(name.str());
	
	/* Set up the IO signals
	  * X and Y have wInX_ and wInY_ bits respectively 
	  * R has wOut_ bits where wOut_ = (wInX_ + wInY_) bits
	  */
	addInput ("X", wIn_);
	addInput ("Y", wIn_);
	addOutput("R", wOut_);
	
	//the 17 should be replaced by a generic multiplierWidth()
	
	int chunks = ( wIn % 17 ==0 ? wIn/17 : ceil( double(wIn)/double(17)) );

	if (chunks == 1)
		cerr << " The TwoWayKaratsuba and the ThreeWayKaratsuba are implemented. Your 17 < wIn <= 51" << endl;
	else if (chunks == 2){
		//pad inputs to 34 bits
		vhdl << tab << declare ("sX", 34) << " <= " << use("X") << " & " << zeroGenerator(34-wIn, 0) << ";" << endl;
		vhdl << tab << declare ("sY", 34) << " <= " << use("Y") << " & " << zeroGenerator(34-wIn, 0) << ";" << endl;
		//chunk splitting
		vhdl << tab << declare ("x0", 18) << " <= " << " \"0\" & " << use("sX") << range(16,0)  << ";" << endl;
		vhdl << tab << declare ("x1", 18) << " <= " << " \"0\" & " << use("sX") << range(33,17) << ";" << endl;
		vhdl << tab << declare ("y0", 18) << " <= " << " \"0\" & " << use("sY") << range(16,0)  << ";" << endl;
		vhdl << tab << declare ("y1", 18) << " <= " << " \"0\" & " << use("sY") << range(33,17) << ";" << endl;
		nextCycle();////////////////////////////////////////////////////////////
		//precomputing
		vhdl << tab << declare ("dx", 18) << " <= " << use("x1") << " - " << use("x0") << ";" << endl;
		vhdl << tab << declare ("dy", 18) << " <= " << use("y1") << " - " << use("y0") << ";" << endl;
		//computing
		vhdl << tab << declare ("r0", 36) << " <= " << use("x0") << " * " << use("y0") << ";" << endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare ("r1", 36) << " <= " << use("r0") << " - " << "(" << use("dx") << " * " << use("dy") << ");"<<endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare ("r2", 36) << " <= " << use("r1") << " + " << "(" << use("x1") << " * " << use("y1") << ");"<<endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare ("t0", 36) << " <= " << use("r0") << ";" << endl;
		vhdl << tab << declare ("t1", 36) << " <= " << use("r2") << ";" << endl;		
		vhdl << tab << declare ("t2", 36) << " <= " << use("r2") << " - " << use("r1") << ";" << endl;		
		vhdl << tab << declare ("a1", 51) << " <= " << use("t2")<<range(33,0) << " & " << use("t0")<<range(33,17) << ";" << endl;
		vhdl << tab << declare ("b1", 51) << " <= " << zeroGenerator(15,0) << " & " << use("r2") << ";" << endl;
		if (68 - 2*wIn < 17 )
			vhdl << tab << "R <= " << "(" << use("a1") << " + " << use("b1") << ") & " << use("t0")<<range(16, (68 - 2*wIn)) << ";" << endl;
		else{
			vhdl << tab << declare("a1b1sum",51) <<  " <= " << use("a1") << " + " << use("b1") << ";" << endl;
			vhdl << tab << "R <= " << use("a1b1sum") << range(50, 68-2*wIn-17) << ";" << endl;
		}
	}else if (chunks == 3){
		//pad inputs to 51 bits
		vhdl << tab << declare ("sX", 51) << " <= " << use("X") << " & " << zeroGenerator(51-wIn, 0) << ";" << endl;
		vhdl << tab << declare ("sY", 51) << " <= " << use("Y") << " & " << zeroGenerator(51-wIn, 0) << ";" << endl;
		//chunk splitting
		vhdl << tab << declare ("x0", 18) << " <= " << " \"0\" & " << use("sX") << range(16,0)  << ";" << endl;
		vhdl << tab << declare ("x1", 18) << " <= " << " \"0\" & " << use("sX") << range(33,17) << ";" << endl;
		vhdl << tab << declare ("x2", 18) << " <= " << " \"0\" & " << use("sX") << range(50,34) << ";" << endl;
		vhdl << tab << declare ("y0", 18) << " <= " << " \"0\" & " << use("sY") << range(16,0)  << ";" << endl;
		vhdl << tab << declare ("y1", 18) << " <= " << " \"0\" & " << use("sY") << range(33,17) << ";" << endl;
		vhdl << tab << declare ("y2", 18) << " <= " << " \"0\" & " << use("sY") << range(50,34) << ";" << endl;
		nextCycle();////////////////////////////////////////////////////////////
		//precomputing
		vhdl << tab << declare ("dx10", 18) << " <= " << use("x1") << " - " << use("x0") << ";" << endl;
		vhdl << tab << declare ("dx20", 18) << " <= " << use("x2") << " - " << use("x0") << ";" << endl;
		vhdl << tab << declare ("dx21", 18) << " <= " << use("x2") << " - " << use("x1") << ";" << endl;
		vhdl << tab << declare ("dy10", 18) << " <= " << use("y1") << " - " << use("y0") << ";" << endl;
		vhdl << tab << declare ("dy20", 18) << " <= " << use("y2") << " - " << use("y0") << ";" << endl;
		vhdl << tab << declare ("dy21", 18) << " <= " << use("y2") << " - " << use("y1") << ";" << endl;
		//computing
		vhdl << tab << declare("p00",36) << " <= " << use("x0") << " * " << use("y0") << ";" << endl;
		vhdl << tab << declare("p11",36) << " <= " << use("x1") << " * " << use("y1")<<";"<<endl;
		vhdl << tab << declare("p22",36) << " <= " << use("x2") << " * " << use("y2")<<";"<<endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare("tk_0",36) << " <= " << use("p00") << " - " << "(" << use("dx10") << " * " << use("dy10") << ");"<<endl;
		vhdl << tab << declare("t2k_0",36) << " <= " << use("p11") << " - " << "(" << use("dx20") << " * " << use("dy20") << ");"<<endl;
		vhdl << tab << declare("t3k_0",36) << " <= " << use("p22") << " - " << "(" << use("dx21") << " * " << use("dy21") << ");"<<endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare("tk",36) << " <= " << use("tk_0") << " + " << use("p11") << ";" << endl;
		vhdl << tab << declare("t2k_1",36) << " <= " << use("t2k_0") << " + " << use("p00") << ";" << endl;
		vhdl << tab << declare("t3k",36) << " <= " << use("t3k_0") << " + " << use("p11") << ";" << endl;
		vhdl << tab << declare("s1_1",34) << " <= " << use("p22")<<range(33,0)<< " + (\"0\" & " << use("t3k")<<range(34,17)<<");"<<endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare("t2k",36) << " <= " << use("t2k_1") << " + " << use("p22") << ";" << endl;
		vhdl << tab << declare("s0",35) << " <= ( \"0\" & " << use("p00")<<range(33,17)<<") + " << use("tk")<<range(34,0) << ";" << endl;
		vhdl << tab << declare("s1",51) << " <= " << use("s1_1") << " & " << use("t3k")<<range(16,0) << ";" << endl;
		vhdl << tab << declare("finalSumLow_0",17) << " <= " << use("p00") << range(16,0) << ";" << endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare("finalSumLow_1",17) << " <= " << use("s0") << range(16,0) << ";" << endl;
		vhdl << tab << declare("l5_0",18) << " <= " << "( \"0\" & " << use("s0")<<range(33,17)<<") + (\"0\" & " << use("t2k")<<range(16,0)<<");"<<endl;
		vhdl << tab << declare("l5_1",27) << " <= " << "( \"0\" & " << use("s1")<<range(25,0)<<") + ("<< zeroGenerator(7,0) << " & " << use("t2k")<<range(35,17)<< ") + " << use("s0") << of(34) <<";"<<endl;
		vhdl << tab << declare("l5_2",25) << " <= " << use("s1") << range(50,26) << ";" << endl;
		nextCycle();////////////////////////////////////////////////////////////
		vhdl << tab << declare("finalSumLow_2",17) << " <= " << use("l5_0") << range(16,0) << ";" << endl;
		vhdl << tab << declare("l6_0",27) << " <= " << use("l5_1")<< " + (\"0\" & " << use("l5_0")<<range(17,17) << ");" << endl;		
		vhdl << tab << declare("l6_1",25) << " <= " << use("l5_2") << ";" << endl;
		
				
		vhdl << tab << declare("finalSumLow_3",26) << " <= " << use("l6_0")<<range(25,0) << ";" << endl;
		vhdl << tab << declare("finalSumLow_4",25) << " <= " << use("l6_1") << " + " << use("l6_0")<<of(26) << ";" << endl;

		
		if (102-2*wIn < 17){                        
			vhdl << tab << " R <= " << use("finalSumLow_4") << " & " 
				                    << use("finalSumLow_3") << " & " 
				                    << use("finalSumLow_2") << " & " 
				                    << use("finalSumLow_1") << " & " 
				                    << use("finalSumLow_0")<<range(16, 102-2*wIn) << ";" << endl;
		}else{
			vhdl << tab << " R <= " << use("finalSumLow_4") << " & " 
				                    << use("finalSumLow_3") << " & " 
				                    << use("finalSumLow_2") << " & " 
				                    << use("finalSumLow_1") << range(16, 102-2*wIn-17) << ";" << endl;
		}
	}
}

void IntKaratsuba::outputVHDL(std::ostream& o, std::string name) {
  
	licence(o);
	o << "library ieee; " << endl;
	o << "use ieee.std_logic_1164.all;" << endl;
	o << "use ieee.std_logic_arith.all;" << endl;
	o << "use ieee.std_logic_signed.all;" << endl;
	o << "library work;" << endl;
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLComponentDeclarations();	
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);		
	o<<buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);

}


IntKaratsuba::~IntKaratsuba() {
}

void IntKaratsuba::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");

	mpz_class svR = svX * svY;

	tc->addExpectedOutput("R", svR);
}

