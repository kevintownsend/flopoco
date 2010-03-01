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
#include <cstdlib>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"
#include "IntKaratsuba.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	IntKaratsuba:: IntKaratsuba(Target* target, int wIn) :
		Operator(target), wIn_(wIn), wOut_(2*wIn){

		ostringstream name;
		name << "IntKaratsuba_" << wIn_<<"_f"<<target->frequencyMHz();
		setCopyrightString("Bogdan Pasca (2008-2009)");
		setName(name.str());
	
		/* Set up the IO signals
		 * X and Y have wInX_ and wInY_ bits respectively 
		 * R has wOut_ bits where wOut_ = (wInX_ + wInY_) bits
		 */
		addInput ("X", wIn_);
		addInput ("Y", wIn_);
		addOutput("R", wOut_);
	
		//TODO replace 17 with a generic multiplierWidth()
	
		int chunks = ( wIn % 17 ==0 ? wIn/17 : ceil( double(wIn)/double(17)) );

		if (chunks > 4){
			REPORT(INFO, " The TwoWayKaratsuba and the ThreeWayKaratsuba are implemented. wIn <= 68 "); 
			exit ( EXIT_FAILURE );	
		}else if (chunks == 1){
			/* no need for karatsuba here */
			vhdl << tab << "R <= X * Y;" << endl;
		}else if (chunks == 2){
			//pad inputs to 34 bits
			vhdl << tab << declare ("sX", 34) << " <= X & " << zg(34-wIn, 0) << ";" << endl;
			vhdl << tab << declare ("sY", 34) << " <= Y & " << zg(34-wIn, 0) << ";" << endl;
			//chunk splitting
			vhdl << tab << declare ("x0", 18) << " <= \"0\" & sX" << range(16,0)  << ";" << endl;
			vhdl << tab << declare ("x1", 18) << " <= \"0\" & sX" << range(33,17) << ";" << endl;
			vhdl << tab << declare ("y0", 18) << " <= \"0\" & sY" << range(16,0)  << ";" << endl;
			vhdl << tab << declare ("y1", 18) << " <= \"0\" & sY" << range(33,17) << ";" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			//precomputing
			vhdl << tab << declare ("dx", 18) << " <= x1 - x0;" << endl;
			vhdl << tab << declare ("dy", 18) << " <= y1 - y0;" << endl;
			//computing
			vhdl << tab << declare ("r0", 36) << " <= x0 * y0;" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare ("r1", 36) << " <= r0 - (dx * dy);"<<endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare ("r2", 36) << " <= r1 + (x1 * y1);"<<endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare ("t0", 36) << " <= r0;" << endl;
			vhdl << tab << declare ("t1", 36) << " <= r2;" << endl;		
			vhdl << tab << declare ("t2", 36) << " <= r2 - r1;" << endl;		
			vhdl << tab << declare ("a1", 51) << " <= t2"<<range(33,0) << " & t0" <<range(33,17) << ";" << endl;
			vhdl << tab << declare ("b1", 51) << " <= " << zg(15,0) << " & r2;" << endl;
			if (68 - 2*wIn < 17 )
				vhdl << tab << "R <= (a1 + b1) & t0"<<range(16, (68 - 2*wIn)) << ";" << endl;
			else{
				vhdl << tab << declare("a1b1sum",51) <<  " <= a1 + b1;" << endl;
				vhdl << tab << "R <= a1b1sum" << range(50, 68-2*wIn-17) << ";" << endl;
			}
		}else if (chunks == 3){
			//pad inputs to 51 bits
			vhdl << tab << declare ("sX", 51) << " <= X & " << zg(51-wIn, 0) << ";" << endl;
			vhdl << tab << declare ("sY", 51) << " <= Y & " << zg(51-wIn, 0) << ";" << endl;
			//chunk splitting
			vhdl << tab << declare ("x0", 18) << " <= \"0\" & sX" << range(16,0)  << ";" << endl;
			vhdl << tab << declare ("x1", 18) << " <= \"0\" & sX" << range(33,17) << ";" << endl;
			vhdl << tab << declare ("x2", 18) << " <= \"0\" & sX" << range(50,34) << ";" << endl;
			vhdl << tab << declare ("y0", 18) << " <= \"0\" & sY" << range(16,0)  << ";" << endl;
			vhdl << tab << declare ("y1", 18) << " <= \"0\" & sY" << range(33,17) << ";" << endl;
			vhdl << tab << declare ("y2", 18) << " <= \"0\" & sY" << range(50,34) << ";" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			//precomputing
			vhdl << tab << declare ("dx10", 18) << " <= x1 - x0;" << endl;
			vhdl << tab << declare ("dx20", 18) << " <= x2 - x0;" << endl;
			vhdl << tab << declare ("dx21", 18) << " <= x2 - x1;" << endl;
			vhdl << tab << declare ("dy10", 18) << " <= y1 - y0;" << endl;
			vhdl << tab << declare ("dy20", 18) << " <= y2 - y0;" << endl;
			vhdl << tab << declare ("dy21", 18) << " <= y2 - y1;" << endl;
			//computing
			vhdl << tab << declare("p00",36) << " <= x0 * y0;" << endl;
			vhdl << tab << declare("p11",36) << " <= x1 * y1;" << endl;
			vhdl << tab << declare("p22",36) << " <= x2 * y2;" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("tk_0", 36) << " <= p00 - (dx10 * dy10);"<<endl;
			vhdl << tab << declare("t2k_0",36) << " <= p11 - (dx20 * dy20);"<<endl;
			vhdl << tab << declare("t3k_0",36) << " <= p22 - (dx21 * dy21);"<<endl;
			nextCycle();////////////////////////////////////////////////////////////
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("tk"   ,36) << " <= tk_0  + p11;" << endl;
			vhdl << tab << declare("t2k_1",36) << " <= t2k_0 + p00;" << endl;
			vhdl << tab << declare("t3k"  ,36) << " <= t3k_0 + p11;" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("s1_1",34) << " <= p22"<<range(33,0)<< " + (\"0\" & t3k"<<range(34,17)<<");"<<endl;
			vhdl << tab << declare("t2k" ,36) << " <= t2k_1 + p22;" << endl;
			vhdl << tab << declare("s0"  ,35) << " <= ( \"0\" & p00"<<range(33,17)<<") + tk"<<range(34,0) << ";" << endl;
			vhdl << tab << declare("s1"  ,51) << " <= s1_1 & t3k"<<range(16,0) << ";" << endl;
			vhdl << tab << declare("finalSumLow_0",17) << " <= p00" << range(16,0) << ";" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("finalSumLow_1",17) << " <= s0" << range(16,0) << ";" << endl;
			vhdl << tab << declare("l5_0",18) << " <= ( \"0\" & s0"<<range(33,17)<<") + (\"0\" & t2k"<<range(16,0)<<");"<<endl;
			vhdl << tab << declare("l5_1",27) << " <= ( \"0\" & s1"<<range(25,0)<<") + ("<< zg(7,0) << " & t2k"<<range(35,17)<< ") + s0"<< of(34) <<";"<<endl;
			vhdl << tab << declare("l5_2",25) << " <= s1" << range(50,26) << ";" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("finalSumLow_2",17) << " <= l5_0" << range(16,0) << ";" << endl;
			vhdl << tab << declare("l6_0",27) << " <= l5_1 + (\"0\" & l5_0"<<range(17,17) << ");" << endl;		
			vhdl << tab << declare("l6_1",25) << " <= l5_2;" << endl;
			nextCycle();////////////////////////////////////////////////////////////
				
			vhdl << tab << declare("finalSumLow_3",26) << " <= l6_0"<<range(25,0) << ";" << endl;
			vhdl << tab << declare("finalSumLow_4",25) << " <= l6_1 + l6_0"<<of(26) << ";" << endl;
		
			if (102-2*wIn < 17){                        
				vhdl << tab << " R <= finalSumLow_4 & finalSumLow_3 & finalSumLow_2 & finalSumLow_1 & finalSumLow_0"<<range(16, 102-2*wIn) << ";" << endl;
			}else{
				vhdl << tab << " R <= finalSumLow_4 & finalSumLow_3 & finalSumLow_2 & finalSumLow_1" << range(16, 102-2*wIn-17) << ";" << endl;
			}
		}else if (chunks == 4){
			//pad inputs to 68 bits
			vhdl << tab << declare ("sX", 68) << " <= X & " << zg(68-wIn, 0) << ";" << endl;
			vhdl << tab << declare ("sY", 68) << " <= Y & " << zg(68-wIn, 0) << ";" << endl;
			//chunk splitting
			vhdl << tab << declare ("x0", 18) << " <= \"0\" & sX" << range(16,0)  << ";" << endl;
			vhdl << tab << declare ("x1", 18) << " <= \"0\" & sX" << range(33,17) << ";" << endl;
			vhdl << tab << declare ("x2", 18) << " <= \"0\" & sX" << range(50,34) << ";" << endl;
			vhdl << tab << declare ("x3", 18) << " <= \"0\" & sX" << range(67,51) << ";" << endl;
			vhdl << tab << declare ("y0", 18) << " <= \"0\" & sY" << range(16,0)  << ";" << endl;
			vhdl << tab << declare ("y1", 18) << " <= \"0\" & sY" << range(33,17) << ";" << endl;
			vhdl << tab << declare ("y2", 18) << " <= \"0\" & sY" << range(50,34) << ";" << endl;
			vhdl << tab << declare ("y3", 18) << " <= \"0\" & sY" << range(67,51) << ";" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			//precomputing
			vhdl << tab << declare ("dx10", 18) << " <= x1 - x0;" << endl;
			vhdl << tab << declare ("dx20", 18) << " <= x2 - x0;" << endl;
			vhdl << tab << declare ("dx21", 18) << " <= x2 - x1;" << endl;
			vhdl << tab << declare ("dx30", 18) << " <= x3 - x0;" << endl;
			vhdl << tab << declare ("dx31", 18) << " <= x3 - x1;" << endl;
			vhdl << tab << declare ("dx32", 18) << " <= x3 - x2;" << endl;

			vhdl << tab << declare ("dy10", 18) << " <= y1 - y0;" << endl;
			vhdl << tab << declare ("dy20", 18) << " <= y2 - y0;" << endl;
			vhdl << tab << declare ("dy21", 18) << " <= y2 - y1;" << endl;
			vhdl << tab << declare ("dy30", 18) << " <= y3 - y0;" << endl;
			vhdl << tab << declare ("dy31", 18) << " <= y3 - y1;" << endl;
			vhdl << tab << declare ("dy32", 18) << " <= y3 - y2;" << endl;

			nextCycle();////////////////////////////////////////////////////////////

			//computing
			vhdl << tab << declare("p00",36) << " <= x0 * y0;" << endl;
			vhdl << tab << declare("p10",36) << " <= dx10 * dy10;" << endl;
			vhdl << tab << declare("p20",36) << " <= dx20 * dy20;" << endl;
			vhdl << tab << declare("p21",36) << " <= dx21 * dy21;" << endl;
			vhdl << tab << declare("p31",36) << " <= dx31 * dy31;" << endl;
			vhdl << tab << declare("p32",36) << " <= dx32 * dy32;" << endl;

			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("t1k_0",36) << " <= p00 + (x1 * y1);"<<endl;
			vhdl << tab << declare("p30",36)   << " <= p21 + (dx30 * dy30);" << endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("t2k_0",36) << " <= t1k_0 + (x2 * y2);"<<endl;
			vhdl << tab << declare("t1k_1",36) << " <= t1k_0 - p10;"<<endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("t3k_0",36) << " <= t2k_0 + (x3 * y3);"<<endl;
			vhdl << tab << declare("t2k_1",36) << " <= t2k_0 - p20;"<<endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("t3k_1",36) << " <= t3k_0 - p30;"<<endl;
			vhdl << tab << declare("t4k_0",36) << " <= t3k_0 - p00;"<<endl;
			vhdl << tab << declare("t5k_0",36) << " <= t3k_0 - t1k_0;"<<endl;
			vhdl << tab << declare("t6k_0",36) << " <= t3k_0 - t2k_0;"<<endl;
			nextCycle();////////////////////////////////////////////////////////////
			vhdl << tab << declare("t4k_1",36) << " <= t4k_0 - p31;"<<endl;
			vhdl << tab << declare("t5k_1",36) << " <= t5k_0 - p32;"<<endl;

			nextCycle();////////////////////////////////////////////////////////////

			IntCompressorTree* c = new IntCompressorTree(target, 119, 7); //68*2 - 17
			oplist.push_back(c);

			vhdl << tab << declare("op0", 119) << " <= " << zg(100,0) << " & p00"<<range(35, 17) << ";" << endl;
			vhdl << tab << declare("op1", 119) << " <= " << zg(83,0) << " & t1k_1"<<range(35, 0) << ";" << endl;
			vhdl << tab << declare("op2", 119) << " <= " << zg(66,0) << " & t2k_1"<<range(35, 0) << " & " << zg(17,0) << ";" << endl;
			vhdl << tab << declare("op3", 119) << " <= " << zg(49,0) << " & t3k_1"<<range(35, 0) << " & " << zg(34,0) << ";" << endl;
			vhdl << tab << declare("op4", 119) << " <= " << zg(32,0) << " & t4k_1"<<range(35, 0) << " & " << zg(51,0) << ";" << endl;
			vhdl << tab << declare("op5", 119) << " <= " << zg(15,0) << " & t5k_1"<<range(35, 0) << " & " << zg(68,0) << ";" << endl;
			vhdl << tab << declare("op6", 119) << " <= " << "t6k_0"<<range(33, 0) << " & " << zg(85,0) << ";" << endl;

			inPortMap( c, "X0", "op0");
			inPortMap( c, "X1", "op1");
			inPortMap( c, "X2", "op2");
			inPortMap( c, "X3", "op3");
			inPortMap( c, "X4", "op4");
			inPortMap( c, "X5", "op5");
			inPortMap( c, "X6", "op6");
			outPortMap(c, "R", "addRes");
			
			vhdl << instance(c, "finalAdder");
			syncCycleFromSignal("addRes");
//			nextCycle();///////////////////////////////////////////////////////

			if (136-2*wIn < 17){                        
				vhdl << tab << " R <= addRes & p00"<<range(16, 136-2*wIn) << ";" << endl;
			}else{
				vhdl << tab << " R <= addRes" << range(118, 136-2*wIn-17) << ";" << endl;
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


	void IntKaratsuba::emulate(TestCase* tc){
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svR = svX * svY;
		tc->addExpectedOutput("R", svR);
	}

}
