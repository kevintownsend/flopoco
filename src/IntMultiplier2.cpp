/*
 * An integer multiplier for FloPoCo
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
#include "IntMultiplier2.hpp"

using namespace std;
extern vector<Operator*> oplist;

IntMultiplier2:: IntMultiplier2(Target* target, int wInX, int wInY) :
	Operator(target), wInX_(wInX), wInY_(wInY), wOut_(wInX + wInY){
 
	ostringstream name;

	/* Name Setup procedure
	 *  The name has the format: IntMultiplier2_wInX__wInY_
	 *  wInX_ = width of the X input
	 *  wInY_ = width of the Y input
	 */  
	name <<"IntMultiplier2_"<<wInX_<<"_"<<wInY_;
	setName(name.str());
	
	addInput ("X", wInX_);
	addInput ("Y", wInY_);
	addOutput("R", wOut_); /* wOut_ = wInX_ + wInY_ */

	if ((wInX==wInY) && (wInX>17) && (wInX<=34)){
		//sigle precision optimizations
		vhdl << tab << declare("sX",34) << " <= " << zeroGenerator(34-wInX,0) << " & " << "X" << ";" << endl;
		vhdl << tab << declare("sY",34) << " <= " << zeroGenerator(34-wInY,0) << " & " << "Y" << ";" << endl;
		
		vhdl << tab << declare("x0",17) << " <= " << use("sX") << range(16,0)  << ";" << endl;
		vhdl << tab << declare("x1",17) << " <= " << use("sX") << range(33,17) << ";" << endl;
		vhdl << tab << declare("y0",17) << " <= " << use("sY") << range(16,0)  << ";" << endl;
		vhdl << tab << declare("y1",17) << " <= " << use("sY") << range(33,17) << ";" << endl;


		vhdl << tab << declare("p00",34) << " <= " << use("x0") << " * " << use("y0") << ";" << endl; 
		nextCycle();
		vhdl << tab << declare("p01",34) << " <= " << use("x0") << " * " << use("y1") << ";" << endl;
		vhdl << tab << declare("ps01",35) << " <= " << "(\"0\" & " << use("p01") << ")" <<  " + " << use("p00")<<range(33,17) << ";" << endl; 
		nextCycle();
		vhdl << tab << declare("p11",34) << " <= " << use("x1") << " * " << use("y1") << ";" << endl;
		vhdl << tab << declare("ps11",34) << " <= "<< use("p11") << " + " << use("ps01")<<range(34,17) << ";" << endl; 
		vhdl << tab << declare("p10",34) << " <= " << use("x1") << " * " << use("y0") << ";" << endl;
		nextCycle();
		vhdl << tab << declare("operand1", 51) << " <= " << use("ps11") << " & " <<  use("ps01") << range(16,0) << ";" <<endl;
		vhdl << tab << declare("operand2", 51) << " <= " << zeroGenerator(51-34,0) << " & " << use("p10") << ";" <<endl;

		IntAdder* add =  new IntAdder(target, 51);
		oplist.push_back(add);
		
		inPortMap (add, "X", "operand1");
		inPortMap (add, "Y", "operand2");
		inPortMapCst(add, "Cin", "'0'");
		outPortMap(add, "R", "addResult");
		vhdl << instance(add, "addition");
		
		syncCycleFromSignal("addResult");
		vhdl << tab << "R <=" << use("addResult")<<range(2*wInX-1-17,0) << " & " << use("p00") << range(16,0) << ";" << endl;  
	}
	if ((wInX==wInY) && (wInX>51) && (wInX<=68)){
		//a 16 DSP version
		//sigle precision optimizations
		vhdl << tab << declare("sX",68) << " <= " << zeroGenerator(68-wInX,0) << " & " << "X" << ";" << endl;
		vhdl << tab << declare("sY",68) << " <= " << zeroGenerator(68-wInY,0) << " & " << "Y" << ";" << endl;

		vhdl << tab << declare("x0",17) << " <= " << "sX" << range(16,0 ) << ";" << endl;
		vhdl << tab << declare("x1",17) << " <= " << "sX" << range(33,17) << ";" << endl;
		vhdl << tab << declare("x2",17) << " <= " << "sX" << range(50,34) << ";" << endl;
		vhdl << tab << declare("x3",17) << " <= " << "sX" << range(67,51) << ";" << endl;
	
		vhdl << tab << declare("y0",17) << " <= " << "sY" << range(16,0 ) << ";" << endl;
		vhdl << tab << declare("y1",17) << " <= " << "sY" << range(33,17) << ";" << endl;
		vhdl << tab << declare("y2",17) << " <= " << "sY" << range(50,34) << ";" << endl;
		vhdl << tab << declare("y3",17) << " <= " << "sY" << range(67,51) << ";" << endl;

		for (int i=0; i<=3; i++){
			ostringstream x, p0, p1, p2 ,p3, tp1, tp2, tp3;
			x << "x" << i;
			p0 << "p" << x.str() << "y0";
			p1 << "p" << x.str() << "y1";
			p2 << "p" << x.str() << "y2";
			p3 << "p" << x.str() << "y3";
			tp1 << "t" << p1.str();
			tp2 << "t" << p2.str();
			tp3 << "t" << p3.str();

			vhdl << tab << declare(p0.str(),34) << " <= " << use(x.str()) << " * " << use("y0") << ";" << endl;
			nextCycle();
			vhdl << tab << declare(tp1.str(),34) << " <= " << use(x.str()) << " * " << use("y1")  << ";" << endl; 
			vhdl << tab << declare(p1.str(),35) << " <= " << "( \"0\" & " << use(tp1.str()) << ")" << " + " << use(p0.str())<<range(33,17) << ";" << endl; 
			nextCycle();
			vhdl << tab << declare(tp2.str(),34) << " <= " << use(x.str()) << " * " << use("y2")  << ";" << endl; 
			vhdl << tab << declare(p2.str(),35) << " <= " << "( \"0\" & " << use(tp2.str()) << ")" << " + " << use(p1.str())<<range(33,17) << ";" << endl; 
			nextCycle();
			vhdl << tab << declare(tp3.str(),34) << " <= " << use(x.str()) << " * " << use("y3")  << ";" << endl; 
			vhdl << tab << declare(p3.str(),35) << " <= " << "( \"0\" & " << use(tp3.str()) << ")" << " + " << use(p2.str())<<range(33,17) << ";" << endl; 
						
			if (i<3) setCycle(0); 
		}
		
		nextCycle();
		vhdl << tab << declare("sum0",85) << " <= " << use("px0y3") << range(33, 0) << " & " << 
		                                               use("px0y2") << range(16, 0) << " & " << 
		                                               use("px0y1") << range(16, 0) << " & " << 
		                                               use("px0y0") << range(16, 0) << ";" << endl;
		vhdl << tab << declare("sum0Low",17) << "<=" << use("sum0")<<range(16,0) << ";" << endl;
		vhdl << tab << declare("sum1",85) << " <= " << use("px1y3") << range(33, 0) << " & " << 
		                                               use("px1y2") << range(16, 0) << " & " << 
		                                               use("px1y1") << range(16, 0) << " & " << 
		                                               use("px1y0") << range(16, 0) << ";" << endl;
		vhdl << tab << declare("sum2",85) << " <= " << use("px2y3") << range(33, 0) << " & " << 
		                                               use("px2y2") << range(16, 0) << " & " << 
		                                               use("px2y1") << range(16, 0) << " & " << 
		                                               use("px2y0") << range(16, 0) << ";" << endl;
		vhdl << tab << declare("sum3",85) << " <= " << use("px3y3") << range(33, 0) << " & " << 
		                                               use("px3y2") << range(16, 0) << " & " << 
		                                               use("px3y1") << range(16, 0) << " & " << 
		                                               use("px3y0") << range(16, 0) << ";" << endl;
		IntNAdder* add =  new IntNAdder(target, 119, 4);
		oplist.push_back(add);
		
		vhdl << tab << declare ("addOp0",119) << " <= " << zeroGenerator(51,0) << " & " << use("sum0") << range(84,17) << ";" <<endl;
		vhdl << tab << declare ("addOp1",119) << " <= " << zeroGenerator(34,0) << " & " << use("sum1") << range(84,0) << ";" <<endl;
		vhdl << tab << declare ("addOp2",119) << " <= " << zeroGenerator(17,0) << " & " << use("sum2") << range(84,0) << " &  " << zeroGenerator(17,0) << ";" <<endl;
		vhdl << tab << declare ("addOp3",119) << " <= " << use("sum3")<< range(84,0)  << " & " <<  zeroGenerator(34,0) << ";" <<endl;
		
		inPortMap (add, "X0", "addOp0");
		inPortMap (add, "X1", "addOp1");
		inPortMap (add, "X2", "addOp2");
		inPortMap (add, "X3", "addOp3");
		inPortMapCst(add, "Cin", "'0'");
		outPortMap(add, "R", "addRes");
		vhdl << instance(add, "adder4");

		syncCycleFromSignal("addRes");
		
		vhdl << tab << "R <= " << use("addRes")<<range(wInX+wInY-1-17,0) << " & " <<  use("sum0Low") << ";" << endl;
	
	}
	if (0==1){
	// to be general version (parametrized etc) 
			//a 16 DSP version
		//sigle precision optimizations
		vhdl << tab << declare("sX",68) << " <= " << zeroGenerator(68-wInX,0) << " & " << "X" << ";" << endl;
		vhdl << tab << declare("sY",68) << " <= " << zeroGenerator(68-wInY,0) << " & " << "Y" << ";" << endl;

		vhdl << tab << declare("x0",17) << " <= " << "sX" << range(16,0 ) << ";" << endl;
		vhdl << tab << declare("x1",17) << " <= " << "sX" << range(33,17) << ";" << endl;
		vhdl << tab << declare("x2",17) << " <= " << "sX" << range(50,34) << ";" << endl;
		vhdl << tab << declare("x3",17) << " <= " << "sX" << range(67,51) << ";" << endl;
	
		vhdl << tab << declare("y0",17) << " <= " << "sY" << range(16,0 ) << ";" << endl;
		vhdl << tab << declare("y1",17) << " <= " << "sY" << range(33,17) << ";" << endl;
		vhdl << tab << declare("y2",17) << " <= " << "sY" << range(50,34) << ";" << endl;
		vhdl << tab << declare("y3",17) << " <= " << "sY" << range(67,51) << ";" << endl;

		for (int i=0; i<=3; i++){
			ostringstream x, p0, p1, p2 ,p3, tp1, tp2, tp3;
			x << "x" << i;
			p0 << "p" << x.str() << "y0";
			p1 << "p" << x.str() << "y1";
			p2 << "p" << x.str() << "y2";
			p3 << "p" << x.str() << "y3";
			tp1 << "t" << p1.str();
			tp2 << "t" << p2.str();
			tp3 << "t" << p3.str();

			vhdl << tab << declare(p0.str(),34) << " <= " << use(x.str()) << " * " << use("y0") << ";" << endl;
			nextCycle();
			vhdl << tab << declare(tp1.str(),34) << " <= " << use(x.str()) << " * " << use("y1")  << ";" << endl; 
			vhdl << tab << declare(p1.str(),35) << " <= " << "( \"0\" & " << use(tp1.str()) << ")" << " + " << use(p0.str())<<range(33,17) << ";" << endl; 
			nextCycle();
			vhdl << tab << declare(tp2.str(),34) << " <= " << use(x.str()) << " * " << use("y2")  << ";" << endl; 
			vhdl << tab << declare(p2.str(),35) << " <= " << "( \"0\" & " << use(tp2.str()) << ")" << " + " << use(p1.str())<<range(33,17) << ";" << endl; 
			nextCycle();
			vhdl << tab << declare(tp3.str(),34) << " <= " << use(x.str()) << " * " << use("y3")  << ";" << endl; 
			vhdl << tab << declare(p3.str(),35) << " <= " << "( \"0\" & " << use(tp3.str()) << ")" << " + " << use(p2.str())<<range(33,17) << ";" << endl; 
						
			if (i<3) setCycle(0); 
		}
		
		nextCycle();
		vhdl << tab << declare("sum0",85) << " <= " << use("px0y3") << range(33, 0) << " & " << 
		                                               use("px0y2") << range(16, 0) << " & " << 
		                                               use("px0y1") << range(16, 0) << " & " << 
		                                               use("px0y0") << range(16, 0) << ";" << endl;
		vhdl << tab << declare("sum0Low",17) << "<=" << use("sum0")<<range(16,0) << ";" << endl;
		vhdl << tab << declare("sum1",85) << " <= " << use("px1y3") << range(33, 0) << " & " << 
		                                               use("px1y2") << range(16, 0) << " & " << 
		                                               use("px1y1") << range(16, 0) << " & " << 
		                                               use("px1y0") << range(16, 0) << ";" << endl;
		vhdl << tab << declare("sum2",85) << " <= " << use("px2y3") << range(33, 0) << " & " << 
		                                               use("px2y2") << range(16, 0) << " & " << 
		                                               use("px2y1") << range(16, 0) << " & " << 
		                                               use("px2y0") << range(16, 0) << ";" << endl;
		vhdl << tab << declare("sum3",85) << " <= " << use("px3y3") << range(33, 0) << " & " << 
		                                               use("px3y2") << range(16, 0) << " & " << 
		                                               use("px3y1") << range(16, 0) << " & " << 
		                                               use("px3y0") << range(16, 0) << ";" << endl;
		IntNAdder* add =  new IntNAdder(target, 119, 4);
		oplist.push_back(add);
		
		vhdl << tab << declare ("addOp0",119) << " <= " << zeroGenerator(51,0) << " & " << use("sum0") << range(84,17) << ";" <<endl;
		vhdl << tab << declare ("addOp1",119) << " <= " << zeroGenerator(34,0) << " & " << use("sum1") << range(84,0) << ";" <<endl;
		vhdl << tab << declare ("addOp2",119) << " <= " << zeroGenerator(17,0) << " & " << use("sum2") << range(84,0) << " &  " << zeroGenerator(17,0) << ";" <<endl;
		vhdl << tab << declare ("addOp3",119) << " <= " << use("sum3")<< range(84,0)  << " & " <<  zeroGenerator(34,0) << ";" <<endl;
		
		inPortMap (add, "X0", "addOp0");
		inPortMap (add, "X1", "addOp1");
		inPortMap (add, "X2", "addOp2");
		inPortMap (add, "X3", "addOp3");
		inPortMapCst(add, "Cin", "'0'");
		outPortMap(add, "R", "addRes");
		vhdl << instance(add, "adder4");

		syncCycleFromSignal("addRes");
		
		vhdl << tab << "R <= " << use("addRes")<<range(wInX+wInY-1-17,0) << " & " <<  use("sum0Low") << ";" << endl;
	

	
	
	
	}

	  
}

IntMultiplier2::~IntMultiplier2() {
}

void IntMultiplier2::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");

	mpz_class svR = svX * svY;

	tc->addExpectedOutput("R", svR);
}

