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

	int chunksX, chunksY;
	chunksX =  int(ceil( ( double(wInX) / 17.0 ) ));
	chunksY =  int(ceil( ( double(wInY) / 17.0 ) ));
	
	if (verbose)
		cout << "X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; " << endl;


	if ((wInX>17) && (wInX<=34) && (wInY>17) && (wInY<=34) ){ //this case is for a 2 x 2 chunk multiplier scheme
		//sigle precision optimizations
		vhdl << tab << declare("sX",34) << " <= " << zeroGenerator(34-wInX,0) << " & " << "X" << ";" << endl;
		vhdl << tab << declare("sY",34) << " <= " << zeroGenerator(34-wInY,0) << " & " << "Y" << ";" << endl;
		
		vhdl << tab << declare("x0",17) << " <= " << use("sX") << range(16,0)  << ";" << endl;
		vhdl << tab << declare("x1",17) << " <= " << use("sX") << range(33,17) << ";" << endl;
		vhdl << tab << declare("y0",17) << " <= " << use("sY") << range(16,0)  << ";" << endl;
		vhdl << tab << declare("y1",17) << " <= " << use("sY") << range(33,17) << ";" << endl;
		
		nextCycle();///////////

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
	}else
		if (chunksX + chunksY > 2) { // up to 17 x 17 bit on Virtex4 can be written as an "*" @ 400++ MHz 
		// to be general version (parametrized etc) 
			//TODO swap X and Y under certain conditions

			//NOT TOO SMART
			bool swap=false;
			if (chunksX > chunksY) 
				swap=true;
			
			if (swap){
				int tmp = chunksX;
				chunksX = chunksY;
				chunksY = tmp; 
			}
			
			if (verbose)
				cout << "Perform swapping = " << swap << endl;
			
			if (swap){
				vhdl << tab << declare("sX",17*chunksX) << " <= " << zeroGenerator(17*chunksX-wInY,0) << " & " << "Y" << ";" << endl;
				vhdl << tab << declare("sY",17*chunksY) << " <= " << zeroGenerator(17*chunksY-wInX,0) << " & " << "X" << ";" << endl;
			}else{
				vhdl << tab << declare("sX",17*chunksX) << " <= " << zeroGenerator(17*chunksX-wInX,0) << " & " << "X" << ";" << endl;
				vhdl << tab << declare("sY",17*chunksY) << " <= " << zeroGenerator(17*chunksY-wInY,0) << " & " << "Y" << ";" << endl;
			}
			////////////////////////////////////////////////////
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++){
				ostringstream dname;
				dname << "x"<<k;
				vhdl << tab << declare(dname.str(),17) << " <= " << "sX" << range((k+1)*17-1,k*17) << ";" << endl;
			}
			for (int k=0; k<chunksY ; k++){
				ostringstream dname;
				dname << "y"<<k;
				vhdl << tab << declare(dname.str(),17) << " <= " << "sY" << range((k+1)*17-1,k*17) << ";" << endl;
			}

			//MULTIPLICATIONS WITH SOME ACCUMULATIONS
			for (int i=0; i<chunksX; i++){ 
				for (int j=0; j<chunksY; j++){ 
					ostringstream currP; //currentProduct
					ostringstream currPS;//currentProductSum
					ostringstream prevPS;//previous Product/ProductSum
					currP  << "tp" << "x" << i << "y" << j;
					currPS <<  "p" << "x" << i << "y" << j;
					prevPS << "p" << "x" << i << "y" << j-1;
					if (j==0){ // @ the first the operation is only multiplication, not MAC
						vhdl << tab << declare(currPS.str(),34) << " <= " << use(join("x",i)) << " * " << use(join("y",j)) << ";" << endl;
					}else{
						vhdl << tab << declare(currP.str(),34) << " <= " << use(join("x",i)) << " * " << use(join("y",j))  << ";" << endl; 
						vhdl << tab << declare(currPS.str(),35) << " <= " << "( \"0\" & " << use(currP.str()) << ")" << " + " << use(prevPS.str())<<range(33,17) << ";" << endl; 
					}
					nextCycle();
				}
				if (i<chunksX-1) setCycle(0);
			}
			//FORM THE INTERMEDIARY PRODUCTS
			vhdl << tab << declare ("sum0Low", 17) << " <= " << use("px0y0")<<range(16,0) << ";" << endl;
			for (int i=0; i<chunksX ; i++){
				vhdl << tab << declare(join("sum",i),17*(chunksY+1)) << " <= "; 
				for (int j=chunksY-1;j>=0; j--){
					ostringstream uname;
					uname << "px" << i << "y"<<j;
					vhdl << use(uname.str()) << range ( (j==chunksY-1?33:16) ,0) << (j==0 ? ";" : " & ");
				}
				vhdl << endl;
			}

			if (chunksX>1){
				IntNAdder* add =  new IntNAdder(target, 17*(chunksX+chunksY-1), chunksX);
				oplist.push_back(add);
		
				for (int i=0; i< chunksX; i++){
					if (i==0) vhdl << tab << declare (join("addOp",i),17*(chunksX+chunksY-1)) << " <= " << zeroGenerator(17*(chunksX-1),0) << " & " << use(join("sum",i)) << range(17*(chunksY+1)-1,17) << ";" <<endl;
					if (i==1) vhdl << tab << declare (join("addOp",i),17*(chunksX+chunksY-1)) << " <= " << zeroGenerator(17*(chunksX-2),0) << " & " << use(join("sum",i)) << range(17*(chunksY+1)-1,0) << ";" <<endl;
					if (i> 1) 
						vhdl << tab << declare (join("addOp",i),17*(chunksX+chunksY-1)) << " <= " << zeroGenerator(17*(chunksX-(i+1)),0) << " & " << use(join("sum",i)) << range(17*(chunksY+1)-1,0) << " & "
								                                                                  << zeroGenerator(17*(i-1),0) << ";" <<endl;
				}
				for (int i=0; i< chunksX; i++)
					inPortMap (add, join("X",i) , join("addOp",i));
		
				inPortMapCst(add, "Cin", "'0'");
				outPortMap(add, "R", "addRes");
				vhdl << instance(add, "adder");

				syncCycleFromSignal("addRes");
		
				vhdl << tab << "R <= " << use("addRes")<<range(wInX+wInY-1-17,0) << " & " <<  use("sum0Low") << ";" << endl;
			}else{
				vhdl << tab << "R <= " << use(join("sum",0))<<range(wInX+wInY-1,0) << ";" << endl;
			}
		}
		else 
			vhdl << tab << "R <= X * Y ;" <<endl;
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

