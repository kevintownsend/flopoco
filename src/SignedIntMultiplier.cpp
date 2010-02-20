/*
 * A signed integer multiplier for FloPoCo
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

#include <typeinfo>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "SignedIntMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	SignedIntMultiplier:: SignedIntMultiplier(Target* target, int wInX, int wInY) :
		Operator(target), wInX_(wInX), wInY_(wInY), wOut_(wInX + wInY){
 		srcFileName = "SignedIntMultiplier";
 
		ostringstream name;
		name <<"SignedIntMultiplier_"<<wInX_<<"_"<<wInY_;
		setName(name.str());
	
		setCopyrightString("Bogdan Pasca (2010)");
	
		addInput ("X", wInX_);
		addInput ("Y", wInY_);
		addOutput("R", wOut_); /* wOut_ = wInX_ + wInY_ */
	
	
		if ((target->getUseHardMultipliers()) && (target->getNumberOfDSPs()>0))
		{
			REPORT(DETAILED, "The target is: " << target->getID());
			if ((target->getID()=="Virtex4")||(target->getID()=="Spartan3")) 
			{
				int x, y, xS, yS, cOp1, cOp2;

				string op1("Y"), op2("X");
				int wOp1(wInY), wOp2(wInX);				

				target->getDSPWidths(xS, yS, true);
				target->getDSPWidths(x,  y,  false);
				 
				cOp1 = getChunkNumber( wOp1, x, xS);
				cOp2 = getChunkNumber( wOp2, y, yS);
				
				REPORT(DEBUG, "cOp1="<<cOp1<<" cOp2="<<cOp2<< " x="<<x<<" xS="<<xS<< " y="<<y<<" yS="<<yS);

				if (cOp1 < cOp2){
					/* swap */
					op1 = string("X"); op2 = string("Y");
					wOp1 = wInX; wOp2 = wInY;
					int tmp = cOp1; cOp1 = cOp2; cOp2 = tmp; 
				}
	
				if (cOp1 + cOp2 > 2) { 

					vhdl << tab << declare("sY", y*(cOp1-1)+yS) << " <= " << op1 << " & " << zg(x*(cOp1-1)+xS -wOp1,0) << ";" << endl;
					vhdl << tab << declare("sX", x*(cOp2-1)+xS) << " <= " << op2 << " & " << zg(y*(cOp2-1)+yS -wOp2,0) << ";" << endl; //pad to the right with 0 (xst)

					//SPLITTINGS
					for (int k=0; k < cOp2 ; k++)
						vhdl << tab << declare(join("x",k), xS) << " <= " << (k<cOp2-1?"\"0\" & ":"") << "sX" << range((k<cOp2-1?(k+1)*x:x*(cOp2-1)+xS)-1,k*x) << ";" << endl;
					for (int k=0; k < cOp1 ; k++)
						vhdl << tab << declare(join("y",k), yS) << " <= " << (k<cOp1-1?"\"0\" & ":"") << "sY" << range((k<cOp1-1?(k+1)*y:y*(cOp1-1)+yS)-1,k*y) << ";" << endl;
				
					//MULTIPLICATIONS WITH SOME ACCUMULATIONS
					for (int i=0; i < cOp2; i++){ 
						for (int j=0; j < cOp1; j++){ 
							if (j==0){ // @ the first the operation is only multiplication, not MAC
								vhdl << tab << declare(join("px",i,"y",j), xS+yS) << " <= " << join("x",i) << " * " << join("y",j) << ";" << endl;
							}else{
								vhdl << tab << declare(join("tpx",i,"y",j),xS+yS) << " <= " << join("x",i) << " * " << join("y",j) << ";" << endl; 
								vhdl << tab << declare(join("px",i,"y",j), xS+yS) << " <= " << join("tpx",i,"y",j) << " + " 
								                                                            << join("px",i,"y",j-1) << range(xS+yS-1,y) << ";" << endl; //sign extend TODO
							}
							if (!((j==cOp1-1) && (i<cOp2-1))) nextCycle();
						}
						if (i<cOp2-1) setCycle(0); //reset cycle
					}
			
					//FORM THE INTERMEDIARY PRODUCTS
					for (int i=0; i < cOp2 ; i++){
						vhdl << tab << declare(join("sum",i), y*(cOp1-1) + xS + yS) << " <= "; 
						for (int j=cOp1-1;j>=0; j--)
							vhdl << join("px",i,"y",j) << (j<cOp1-1?range(x-1,0):"") << (j==0 ? ";" : " & ");
						vhdl << endl;
					}

					if (cOp2 > 1){
						vhdl << tab << declare ("sum0Low", x) << " <= " << use("sum0")<<range(x-1,0) << ";" << endl;

						IntNAdder* add =  new IntNAdder(target, y*(cOp1-1)+yS +  x*(cOp2-1)+xS - x, cOp2);
						oplist.push_back(add);
		
						/* prepare operands */
						for (int i=0; i < cOp2; i++){
							if (i==0){ 
								vhdl << tab << declare (join("addOp",i), (cOp2-1)*x + y*(cOp1-1) + xS + yS - x ) << " <= "
								            << rangeAssign( (cOp2-1-i)*x-1, 0, join("sum",i)+of(y*(cOp1-1) + xS + yS - 1 ))  << " & " 
								            << join("sum",i) << range(y*(cOp1-1) + xS + yS - 1,x) << ";" <<endl;
							}else if (i==1){ 
								vhdl << tab << declare (join("addOp",i),(cOp2-1)*x + y*(cOp1-1) + xS + yS - x) << " <= " 
								            << rangeAssign( (cOp2-1-i)*x-1, 0, join("sum",i)+of(y*(cOp1-1) + xS + yS - 1 ))  << ((cOp2-1-i)*x-1>=0?" & ":"") 
								            << join("sum",i) << ";" <<endl;
							}else if ((i > 1) && ( i!= cOp2-1)){ 
								vhdl << tab << declare (join("addOp",i),(cOp2-1)*x + y*(cOp1-1) + xS + yS - x) << " <= " 
								            << rangeAssign( (cOp2-1-i)*x-1, 0, join("sum",i)+of(y*(cOp1-1) + xS + yS - 1 ))  << ((cOp2-1-i)*x-1>=0?" & ":"") 
								            << join("sum",i)  << " & "
								            << zg( (i-1)*x, 0) << ";" <<endl;
							}else if (i == cOp2-1){
								vhdl << tab << declare (join("addOp",i),(cOp2-1)*x + y*(cOp1-1) + xS + yS - x) << " <= " 
								            << join("sum",i)  << " & "
								            << zg( (i-1)*x, 0) << ";" <<endl;
							}
						}//TODO compact ifs

						for (int i=0; i< cOp2; i++)
							inPortMap (add, join("X",i) , join("addOp",i));
		
						inPortMapCst(add, "Cin", "'0'");
						outPortMap(add, "R", "addRes");
						vhdl << instance(add, "adder");

						syncCycleFromSignal("addRes");

						if ( (y*(cOp1-1)+yS + x*(cOp2-1)+xS ) - (wInX + wInY) < x )
							vhdl << tab << "R <= addRes & sum0Low "<<range(x-1, (y*(cOp1-1)+yS + x*(cOp2-1)+xS ) - (wInX + wInY) )<< ";" << endl;
						else						
							vhdl << tab << "R <= addRes"<<range(y*(cOp1-1)+yS +  x*(cOp2-1)+xS - x -1, (y*(cOp1-1)+yS + x*(cOp2-1)+xS ) - (wInX + wInY) + x )<< ";" << endl;
					}else{
						vhdl << tab << "R <= sum0"<<range( y*(cOp1-1) + xS + yS - 1, y*(cOp1-1) + xS + yS - (wInX + wInY ) ) << ";" << endl;
					}
				}
//				else 
//					vhdl << tab << "R <= X * Y ;" <<endl;
			}else{
				cerr << " Sorry: only Virtex4 and Spartan3 are supported at this time. We'll be back soon" << endl;
				throw "Error";
			}
		}else
			vhdl << tab << "R <= X * Y ;" <<endl;
	}

	SignedIntMultiplier::~SignedIntMultiplier() {
	}
	
	void SignedIntMultiplier::outputVHDL(std::ostream& o, std::string name) {
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
	

	void SignedIntMultiplier::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
		mpz_class big1 = (mpz_class(1) << (wInX_));
		mpz_class big1P = (mpz_class(1) << (wInX_-1));
		mpz_class big2 = (mpz_class(1) << (wInY_));
		mpz_class big2P = (mpz_class(1) << (wInY_-1));

		if ( svX >= big1P)
			svX = svX-big1;

		if ( svY >= big2P)
			svY = svY -big2;
			
		mpz_class svR = svX * svY;
		if ( svR < 0){
			mpz_class tmpSUB = (mpz_class(1) << (wOut_));
			svR = tmpSUB + svR; 
		}

		tc->addExpectedOutput("R", svR);
	}
}
