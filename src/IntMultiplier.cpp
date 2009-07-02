/*
 * An integer multiplier for FloPoCo
 *
 * Authors : Bogdan Pasca, Sebastian Banescu
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
#include "IntMultiplier.hpp"

using namespace std;
extern vector<Operator*> oplist;

IntMultiplier:: IntMultiplier(Target* target, int wInX, int wInY) :
	Operator(target), wInX_(wInX), wInY_(wInY), wOut_(wInX + wInY){
 
	ostringstream name;

	/* Name Setup procedure
	 *  The name has the format: IntMultiplier_wInX__wInY_
	 *  wInX_ = width of the X input
	 *  wInY_ = width of the Y input
	 */  
	name <<"IntMultiplier_"<<wInX_<<"_"<<wInY_;
	setName(name.str());
	
	setCopyrightString("Bogdan Pasca, Sebastian Banescu (2008-2009)");
	
	addInput ("X", wInX_);
	addInput ("Y", wInY_);
	addOutput("R", wOut_); /* wOut_ = wInX_ + wInY_ */
	
	if (target->getUseHardMultipliers())
	{
	int chunksX, chunksY;
	int x, y;
	target->suggestSubmultSize(x, y, wInX, wInY);
	int chunkSize_ = x;
	chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
	chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
	
	if (verbose)
		cout << "X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; " << endl;

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
				vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= " << "Y" << " & " << zeroGenerator(chunkSize_*chunksX-wInY,0) << ";" << endl;
				vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= " << "X" << " & " << zeroGenerator(chunkSize_*chunksY-wInX,0) << ";" << endl;
			}else{
				vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= " << "X" << " & " << zeroGenerator(chunkSize_*chunksX-wInX,0) << ";" << endl;
				vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= " << "Y" << " & " << zeroGenerator(chunkSize_*chunksY-wInY,0) << ";" << endl;
			}
			////////////////////////////////////////////////////
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++){
				ostringstream dname;
				dname << "x"<<k;
				vhdl << tab << declare(dname.str(),chunkSize_) << " <= " << "sX" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
			}
			for (int k=0; k<chunksY ; k++){
				ostringstream dname;
				dname << "y"<<k;
				vhdl << tab << declare(dname.str(),chunkSize_) << " <= " << "sY" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
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
						vhdl << tab << declare(currPS.str(),2*chunkSize_) << " <= " << use(join("x",i)) << " * " << use(join("y",j)) << ";" << endl;
					}else{
						vhdl << tab << declare(currP.str(),2*chunkSize_) << " <= " << use(join("x",i)) << " * " << use(join("y",j))  << ";" << endl; 
						vhdl << tab << declare(currPS.str(),2*chunkSize_+1) << " <= " << "( \"0\" & " << use(currP.str()) << ")" << " + " << use(prevPS.str())<<range(2*chunkSize_-1,chunkSize_) << ";" << endl; 
					}
					nextCycle();
				}
				if (i<chunksX-1) setCycle(0);
			}
			//FORM THE INTERMEDIARY PRODUCTS
			vhdl << tab << declare ("sum0Low", chunkSize_) << " <= " << use("px0y0")<<range(chunkSize_-1,0) << ";" << endl;
			for (int i=0; i<chunksX ; i++){
				vhdl << tab << declare(join("sum",i),chunkSize_*(chunksY+1)) << " <= "; 
				for (int j=chunksY-1;j>=0; j--){
					ostringstream uname;
					uname << "px" << i << "y"<<j;
					vhdl << use(uname.str()) << range ( (j==chunksY-1?(2*chunkSize_-1):chunkSize_-1) ,0) << (j==0 ? ";" : " & ");
				}
				vhdl << endl;
			}

			if (chunksX>1){
				IntNAdder* add =  new IntNAdder(target, chunkSize_*(chunksX+chunksY-1), chunksX);
				oplist.push_back(add);
		
				for (int i=0; i< chunksX; i++){
					if (i==0) vhdl << tab << declare (join("addOp",i),chunkSize_*(chunksX+chunksY-1)) << " <= " << zeroGenerator(chunkSize_*(chunksX-1),0) << " & " << use(join("sum",i)) << range(chunkSize_*(chunksY+1)-1,chunkSize_) << ";" <<endl;
					if (i==1) vhdl << tab << declare (join("addOp",i),chunkSize_*(chunksX+chunksY-1)) << " <= " << zeroGenerator(chunkSize_*(chunksX-2),0) << " & " << use(join("sum",i)) << range(chunkSize_*(chunksY+1)-1,0) << ";" <<endl;
					if (i> 1) 
						vhdl << tab << declare (join("addOp",i),chunkSize_*(chunksX+chunksY-1)) << " <= " << zeroGenerator(chunkSize_*(chunksX-(i+1)),0) << " & " << use(join("sum",i)) << range(chunkSize_*(chunksY+1)-1,0) << " & "
								                                                                  << zeroGenerator(chunkSize_*(i-1),0) << ";" <<endl;
				}
				for (int i=0; i< chunksX; i++)
					inPortMap (add, join("X",i) , join("addOp",i));
		
				inPortMapCst(add, "Cin", "'0'");
				outPortMap(add, "R", "addRes");
				vhdl << instance(add, "adder");

				syncCycleFromSignal("addRes");
		
				if ((chunkSize_*(chunksX+chunksY)) - wInX - wInY < chunkSize_){
					vhdl << tab << "R <= " << use("addRes")<<range((chunkSize_*(chunksX+chunksY-1))-1,0) << " & " <<  use("sum0Low")<<range(chunkSize_-1, chunkSize_-1 + ((chunkSize_*(chunksX+chunksY-1)) - wInX - wInY)+1) << ";" << endl;
				}else{
					vhdl << tab << "R <= " << use("addRes")<<range((chunkSize_*(chunksX+chunksY-1))-1, ((chunkSize_*(chunksX+chunksY-1)) - wInX - wInY)) << ";" << endl;
				
				}
			}else{
				vhdl << tab << "R <= " << use(join("sum",0))<<range(chunkSize_*(chunksX+chunksY)-1, chunkSize_*(chunksX+chunksY) -( wInX+wInY)) << ";" << endl;
			}
		}
		else 
			vhdl << tab << "R <= X * Y ;" <<endl;
	}
	else
	{
		int chunkSize_ = target->lutInputs()/2;
	
	int chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
	int chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
		
	if (verbose)
        cout << "X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; " << endl;
	if (chunksX + chunksY > 2) {
	int widthX = wInX_;
	int widthY = wInY_;	
	
	if (chunksX > chunksY)
	{
		int tmp = chunksX;
		chunksX = chunksY;
		chunksY = tmp;
		
		tmp = widthX;
		widthX = widthY;
		widthY = tmp;
	}
	
	vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= " << "X" << " & " << zeroGenerator(chunkSize_*chunksX-widthX,0) << ";" << endl;
    vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= " << "Y" << " & " << zeroGenerator(chunkSize_*chunksY-widthY,0) << ";" << endl;
    
	////////////////////////////////////////////////////
    //SPLITTINGS
    for (int k=0; k<chunksX ; k++)
	{
		ostringstream dname;
        dname << "x"<<k;
        vhdl << tab << declare(dname.str(),chunkSize_) << " <= " << "sX" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
    }
    for (int k=0; k<chunksY ; k++)
	{
		ostringstream dname;
        dname << "y"<<k;
        vhdl << tab << declare(dname.str(),chunkSize_) << " <= " << "sY" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
    }

	// COMPUTE PARTIAL PRODUCTS
	for (int i=0; i<chunksY; i++)
	{
		for (int j=0; j<chunksX; j++)
		{
			ostringstream partialProd;
            partialProd  << "p" << "x" << j << "y" << i;
			vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("x",j)) << " * " << use(join("y",i)) << ";" << endl;
		}
	}		
	
	int adderWidth = (chunksX+chunksY)*chunkSize_;
	
	// CONCATENATE PARTIAL PRODUCTS
	for (int i=0; i<chunkSize_; i++)
	{
		for (int j=0; j<chunksX; j++)
		{
			ostringstream concatPartialProd;
            concatPartialProd  << "cp" << i << j;
			vhdl << tab << declare(concatPartialProd.str(),adderWidth) << " <= ";
			
			int startIdx = chunksY-1-i;
			int paddWidth = adderWidth - 2*chunkSize_*(startIdx/chunkSize_+1);
			int endPaddWidth = chunkSize_*(j+startIdx%chunkSize_);
			vhdl << zeroGenerator(paddWidth-endPaddWidth,0); 
			
			for (int k=startIdx; k>=0; k-=chunkSize_)
			{
				ostringstream partialProd;
				partialProd  << "p" << "x" << j << "y" << k;
				vhdl << " & " << use(partialProd.str());
			}
			vhdl << " & " << zeroGenerator(endPaddWidth,0) << ";" << endl;
		}
	}
	
	IntNAdder* add =  new IntNAdder(target, adderWidth, chunksX*chunkSize_);
    oplist.push_back(add);
	
	for (int i=0; i<chunkSize_; i++)
		for (int j=0; j<chunksX; j++)
		{
			ostringstream concatPartialProd;
            concatPartialProd  << "cp" << i << j;
			
			inPortMap (add, join("X",i*chunksX+j) , concatPartialProd.str());
		}	
			
	inPortMapCst(add, "Cin", "'0'");
    outPortMap(add, "R", "addRes");
    vhdl << instance(add, "adder");

    syncCycleFromSignal("addRes");
	
	vhdl << tab << "R <= " << use("addRes")<<range(adderWidth-1,adderWidth-wInX_-wInY_) << ";" << endl;		
	}
		else 
			vhdl << tab << "R <= X * Y ;" <<endl;
	}
}

IntMultiplier::~IntMultiplier() {
}

void IntMultiplier::emulate(TestCase* tc)
{
	mpz_class svX = tc->getInputValue("X");
	mpz_class svY = tc->getInputValue("Y");

	mpz_class svR = svX * svY;

	tc->addExpectedOutput("R", svR);
}

