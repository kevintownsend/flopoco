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
		if (target->lutInputs() == 4) // then the target is Virtex 4 (TODO: change this to something like: if (target instanceof Virtex4))
		{
			int chunksX, chunksY;
			int x, y;
			target->suggestSubmultSize(x, y, wInX, wInY);
			int chunkSize_ = x;
			chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
			chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
			
			if (verbose)
				cerr << "> IntMultiplier:   X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; " << endl;

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
						cerr << "> IntMultiplier:  Perform swapping = " << swap << endl;
					
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
		else // the target is a StratixII
		{
			double f = target->frequencyMHz();
			bool quadMultiply = false; // TRUE if we cannot use double the chunckSize width multipliers

				int zeroPad18 = ceil((double)wInX/18)*18 - wInX + ceil((double)wInY/18)*18 - wInY;
				int zeroPad9 = ceil((double)wInX/9)*9 - wInX + ceil((double)wInY/9)*9 - wInY;
				int chunkSize_ = 18;
			
				if ((zeroPad9 < zeroPad18) || 
					(f > 367))
					chunkSize_ = 9;
			
				if ((f > 250 && chunkSize_ == 18) || 	// don't use 36x36
					(f > 367 && chunkSize_ == 9)) 	// don't use 18x18
						quadMultiply = true;
				
				int chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
				int chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
			
				int widthX = wInX_;
				int widthY = wInY_;	
			
				if (chunksX > chunksY) // the tiling is done on a vertical rectangle
				{
					int tmp = chunksX;
					chunksX = chunksY;
					chunksY = tmp;
				
					tmp = widthX;
					widthX = widthY;
					widthY = tmp;
				
					vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= " << "Y" << " & " << zeroGenerator(chunkSize_*chunksX-widthX,0) << ";" << endl;
					vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= " << "X" << " & " << zeroGenerator(chunkSize_*chunksY-widthY,0) << ";" << endl;
				}
				else
				{
					vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= " << "X" << " & " << zeroGenerator(chunkSize_*chunksX-widthX,0) << ";" << endl;
					vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= " << "Y" << " & " << zeroGenerator(chunkSize_*chunksY-widthY,0) << ";" << endl;
				}
				
				if (verbose)
					cerr << "> IntMultiplier:  X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; The chunk size is " << chunkSize_ << ";" << endl;

				int chX = chunksX; 	// number of chunks of X and decreases at each iteration of the while loop below
				int level = 0;		// index of the current iteration of the while loop
				int height;			// height of the tiling after the while loop
				
				int adderWidth = (chunksX+chunksY)*chunkSize_;	// width of the final adder
				int opCount = 0;								// current index of the operands of the final adder
				
				ostringstream partialProd;	// temporary holder of a partial product
				ostringstream partialProd2; // temporary holder of a partial product
				ostringstream sum;			// temporary holder of a sum of partial products
				ostringstream operands[2];	// operands of the final summation that hold the concatenation of partial products
				ostringstream carrys;		// operand of the final summation that holds the concatenation of carry bits of the partial products

			// COMPUTE PARTIAL PRODUCTS
			if (quadMultiply)
			{
				
				////////////////////////////////////////////////////
				//SPLITTINGS
				for (int k=0; k<chunksX ; k++)
				{
					ostringstream dname;
					dname << "x"<<k;
					vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
				}
				for (int k=0; k<chunksY ; k++)
				{
					ostringstream dname;
					dname << "y"<<k;
					vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
				}
				
				if (verbose)	
					cerr << "> IntMultiplier: Cannot use double chunckSize_ multipliers" << endl;
					
				while (chX/4 > 0)
				{
					// top-right tile
					setCycle(1);
					partialProd.str("");
					partialProd << "px0y" << 4*level;
					vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y", 4*level)) << ";" << endl;
					
					setCycle(3);
					operands[0] << use(partialProd.str()) << " & " << zeroGenerator(chunkSize_*4*level,0) << ";" << endl;
					
					
					// top-right diagonal pair
					setCycle(1);
					partialProd2.str("");
					partialProd2 << "px0y" << 4*level+1;
					vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",4*level+1)) << ";" << endl;
					partialProd.str("");
					partialProd << "px1y" << 4*level;
					vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y",4*level)) << ";" << endl;
					
					setCycle(2);
					sum.str("");
					sum << "addOp" << level << "_shift_" << 4*level+1;
					vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
					
					setCycle(3);
					operands[1] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << zeroGenerator((4*level+1)*chunkSize_,0) << ";" << endl;
					carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zeroGenerator((4*level+3)*chunkSize_,0) << ";" << endl;
					
					
					// top-right diagonal triplet
					setCycle(1);
					// first compute the 3 partial products
					for (int j=0; j<3; j++)
					{
						partialProd.str("");
						partialProd << "px" << j << "y" << 4*level+2-j;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",4*level+2-j)) << ";" << endl; 
					}
					setCycle(2);
					// add the 3 partial products
					sum.str("");
					sum << "addOp" << level << "_shift_" << 4*level+2;
					vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
					
					for (int j=0; j<3; j++)
					{
						partialProd.str("");
						partialProd << "px" << j << "y" << 4*level+2-j;
						if (j == 0)
							vhdl << use(partialProd.str());
						else
							vhdl << ") + (\"00\" & " << use(partialProd.str());
					}
					vhdl << ");" << endl;
					
					setCycle(3);
					operands[0].seekp(ios_base::beg);
					operands[0] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[0].str();
					carrys.seekp(ios_base::beg);
					carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-1,0) << " & " << carrys.str();
					
						
					//setCycle(0);
						
					// compute the sum of four tiles (diagonal) starting from top-right going left
					for (int i=3; i<chX-1; i++)
					{
						setCycle(1);
						// first compute the 4 partial products
						for (int j=0; j<4; j++)
						{
							partialProd.str("");
							partialProd << "px" << i-j << "y" << j+4*level;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",i-j)) << " * " << use(join("y",j+4*level)) << ";" << endl;
						}
						setCycle(2);
						// add the 4 partial products
						sum.str("");
						sum << "addOp" << level << "_shift_" << i+4*level;
						vhdl << tab << declare(sum.str(),2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
						
						for (int j=0; j<4; j++)
						{
							partialProd.str("");
							partialProd << "px" << i-j << "y" << j+4*level;
							if (j == 0)
								vhdl << use(partialProd.str());
							else
								vhdl << ") + (\"00\" & " << use(partialProd.str());
						}
						vhdl << ");" << endl;
						
						setCycle(3);
						operands[i%2].seekp(ios_base::beg);
						operands[i%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[i%2].str();  
						carrys.seekp(ios_base::beg);
						carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
					}
					
					// compute the sum of four tiles (diagonal) going down to the bottom-left corner
					for (int j=4*level; j<chunksY-3; j++)
					{
						setCycle(1);
						// first compute the 4 partial products
						for (int i=0; i<4; i++)
						{
							partialProd.str("");
							partialProd << "px" << chX-i-1 << "y" << i+j;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-i-1)) << " * " << use(join("y",i+j)) << ";" << endl; 
						}
						setCycle(2);
						// add the 4 partial products
						sum.str("");
						sum << "addOp" << level << "_shift_" << chX+j-1;
						vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
						
						for (int i=0; i<4; i++)
						{
							partialProd.str("");
							partialProd << "px" << chX-i-1 << "y" << i+j;
							if (i == 0)
								vhdl << use(partialProd.str());
							else
								vhdl << ") + (\"00\" & " << use(partialProd.str());
						}
						vhdl << ");" << endl;
						
						setCycle(3);
						operands[(j+chX-1)%2].seekp(ios_base::beg);
						operands[(j+chX-1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(j+chX-1)%2].str(); 
						carrys.seekp(ios_base::beg);
						carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
					} 
					
					// bottom-left diagonal triplet
					setCycle(1);
					// first compute the 3 partial products
					for (int j=0; j<3; j++)
					{
						partialProd.str("");
						partialProd << "px" << chX-3+j << "y" << chunksY-1-j;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-3+j)) << " * " << use(join("y",chunksY-1-j)) << ";" << endl; 
					}
					setCycle(2);
					// add the 3 partial products
					sum.str("");
					sum << "addOp" << level << "_shift_" << chunksY+chX-4;
					vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
					
					for (int j=0; j<3; j++)
					{
						partialProd.str("");
						partialProd << "px" << chX-3+j << "y" << chunksY-1-j;
						if (j == 0)
							vhdl << use(partialProd.str());
						else
							vhdl << ") + (\"00\" & " << use(partialProd.str());
					}
					vhdl << ");" << endl;
					
					setCycle(3);
					operands[(chX+chunksY-4)%2].seekp(ios_base::beg);
					operands[(chX+chunksY-4)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(chX+chunksY-4)%2].str();
					carrys.seekp(ios_base::beg);
					carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
					
					// bottom-left diagonal pair
					setCycle(1);
					partialProd2.str("");
					partialProd2 << "px" << chX-1 << "y" << chunksY-2;
					vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-1)) << " * " << use(join("y", chunksY-2)) << ";" << endl;
					partialProd.str("");
					partialProd << "px" << chX-2 << "y" << chunksY-1;
					vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x", chX-2)) << " * " << use(join("y", chunksY-1)) << ";" << endl;
					
					setCycle(2);
					sum.str("");
					sum << "addOp" << level << "_shift_" << chX+chunksY-3;
					vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
					
					setCycle(3);
					operands[(chX+chunksY-3)%2].seekp(ios_base::beg);
					operands[(chX+chunksY-3)%2] << zeroGenerator((4*level+1)*chunkSize_,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(chX+chunksY-3)%2].str();
					carrys.seekp(ios_base::beg);
					carrys << zeroGenerator((4*level+1)*chunkSize_-1,0) << " & " << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
					
					// bottom-left tile
					setCycle(1);
					partialProd.str("");
					partialProd << "px"<< chX-1 << "y" << chunksY-1;
					vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-1)) << " * " << use(join("y", chunksY-1)) << ";" << endl;
					
					setCycle(3);
					operands[(chX+chunksY-2)%2].seekp(ios_base::beg);
					operands[(chX+chunksY-2)%2] << zeroGenerator(4*level*chunkSize_,0) << " & " << use(partialProd.str()) << " & " << operands[(chX+chunksY-2)%2].str();
					
					// reinitialize the concatenation buffers for the next iteration of the while loop
					for (int i=0; i<2; i++)
					{
						setCycle(3);
						vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[i].str() << endl;
						opCount++;
						operands[i].str("");
					}
					vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
					opCount++;	
					carrys.str("");
					
					chX -= 4;
					level++;
				}
				
				int i = level*4;	
				bool oneCol = true;		
				switch (chX)
				{
					case 1: // remaining tiles are in a single column
						if (verbose)
							cerr << ">IntMultiplier: Case 1: remaining tiles are in a single column" << endl;
						
						operands[i%2] << zeroGenerator(level*4*chunkSize_,0) << ";";
						if (chunksY-i > 1) // then there are more than one tiles in the column
						{
							if (verbose)
								cerr << ">IntMultiplier: More than one tile in the column" << endl;
							operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << ";";
							oneCol = false;
						}
						
						for (i=level*4; i<chunksY; i++)
						{
							setCycle(1);
							partialProd.str("");
							partialProd << "px0y" << i;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i)) << ";" << endl;
							
							setCycle(3);
							operands[i%2].seekp(ios_base::beg);
							operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
							
						}
						
						if (!oneCol) // then there are more than one tiles in the column
						{
							operands[i%2].seekp(ios_base::beg);
							operands[i%2] << zeroGenerator((level*4+1)*chunkSize_,0) << " & " << operands[i%2].str();
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[i%2].str() << endl;
							opCount++;
						}
						
						operands[(i+1)%2].seekp(ios_base::beg);
						operands[(i+1)%2] << zeroGenerator(level*4*chunkSize_,0) << " & " << operands[(i+1)%2].str(); 
						vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(i+1)%2].str() << endl;
						opCount++;
							
						break;
					case 2: // remaining tiles are in 2 columns
						if (verbose)
							cerr << ">IntMultiplier: Case 2: remaining tiles are in 2 columns" << endl;
						setCycle(1);
						// top right tile
						partialProd.str("");
						partialProd << "px0y" << i;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y", i)) << ";" << endl;
						setCycle(3);
						operands[i%2] << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << zeroGenerator((4*level)*chunkSize_,0) << ";";
						operands[(i+1)%2] << zeroGenerator((4*level+1)*chunkSize_,0) << ";";
						carrys << zeroGenerator((4*level+2)*chunkSize_+1,0) << ";";
						 
						while (i<chunksY-1) // diagonal pairs
						{
							setCycle(1);
							partialProd2.str("");
							partialProd2 << "px0y" << (i+1);
							vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y", i+1)) << ";" << endl;
							partialProd.str("");
							partialProd << "px1y" << i;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y", i)) << ";" << endl;
							setCycle(2);
							sum.str("");
							sum << "addOp" << level << "_shift_" << (i+1);
							vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
							setCycle(3);
							operands[(i+1)%2].seekp(ios_base::beg);
							operands[(i+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
							carrys.seekp(ios_base::beg);
							carrys << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-1,0) << " & " << carrys.str();
					
							i++;			
						}
						
						// bottom-left tile
						setCycle(1);
						partialProd.str("");
						partialProd << "px1y" << i;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y", i)) << ";" << endl;
						
						setCycle(3);
						operands[(i+1)%2].seekp(ios_base::beg);
						operands[(i+1)%2] << zeroGenerator((4*level)*chunkSize_,0) << " & " << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
						
						operands[i%2].seekp(ios_base::beg);
						operands[i%2] << zeroGenerator((4*level+1)*chunkSize_,0) << " & " << operands[i%2].str();
						carrys.seekp(ios_base::beg);
						carrys << zeroGenerator((4*level+1)*chunkSize_-1,0) << " & " << carrys.str();
						
						// reinitialize the concatenation buffers for the next iteration of the while loop
						for (int j=0; j<2; j++)
						{
							setCycle(3);
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[j].str() << endl;
							opCount++;
							operands[j].str("");
						}
						vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
						opCount++;	
						carrys.str("");
						break;
					case 3: // Remaining tiles are on 3 columns
						if (verbose)
							cerr << ">IntMultiplier: Case 3: Remaining tiles are on 3 columns" << endl;
						
						height = chunksY - i;
						
						switch (height)
						{
							case 0: break;
							case 1: // only one row
								if (verbose)
									cerr << ">IntMultiplier:" << tab << " Subcase 1: there is one row" << endl;
							
								operands[i%2] << zeroGenerator(level*4*chunkSize_,0) << ";";
								operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << ";";
						
								for (int j=0; j<3; j++) // each tile is computed and padded with zeros becoming an operand of IntNAdder
								{
									setCycle(1);
									partialProd.str("");
									partialProd << "px" << j << "y" << i;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",i)) << ";" << endl;
									setCycle(3);
									operands[(j+i)%2].seekp(ios_base::beg);
									operands[(j+i)%2] << use(partialProd.str()) << " & " << operands[(j+i)%2].str();
								}
								
								operands[(3+i)%2].seekp(ios_base::beg);
								operands[(3+i)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << " & " << operands[(3+i)%2].str();
								vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(3+i)%2].str() << endl;
								opCount++;
								
								operands[(i+4)%2].seekp(ios_base::beg);
								operands[(i+4)%2] << zeroGenerator(level*4*chunkSize_,0) << " & " << operands[(i+4)%2].str(); 
								vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(i+4)%2].str() << endl;
								opCount++;
								break;
							case 2: // two rows
								if (verbose)
									cerr << ">IntMultiplier:" << tab << " Subcase 2: there are two rows" << endl;
								
								operands[i%2] << zeroGenerator(level*4*chunkSize_,0) << ";";
								operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << ";";
								carrys << zeroGenerator(level*4*chunkSize_,0) << ";";
								// top-right tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px0y" << i;
								vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i)) << ";" << endl;
								operands[i%2].seekp(ios_base::beg);
								operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
								
								// diagonal tiles
								for (int j=0; j<2; j++)
								{
									setCycle(1);
									partialProd.str("");
									partialProd << "px"<< j+1 << "y" << i;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j+1)) << " * " << use(join("y",i)) << ";" << endl;
									partialProd2.str("");
									partialProd2 << "px"<< j << "y" << i+1;
									vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",i+1)) << ";" << endl;
									
									setCycle(2);
									sum.str("");
									sum << "addOpd" << level << "_shift_" << i+j+1;
									vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
									
									setCycle(3);
									operands[(i+j+1)%2].seekp(ios_base::beg);
									operands[(i+j+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) <<" & " << operands[(i+j+1)%2].str();
									carrys.seekp(ios_base::beg);
									carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zeroGenerator(chunkSize_,0) << "&" << carrys.str();
								}
								
								// bottom-left tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px2y" << i+2;
								vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x2") << " * " << use(join("y",i+2)) << ";" << endl;
								setCycle(3);
								operands[(i+4)%2].seekp(ios_base::beg);
								operands[(i+4)%2] << zeroGenerator(level*4*chunkSize_,0) << " & " << use(partialProd.str()) << " & " << operands[(i+4)%2].str();	
								operands[(i+5)%2].seekp(ios_base::beg);
								operands[(i+5)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << " & " << operands[(i+5)%2].str(); 
								carrys.seekp(ios_base::beg);
								carrys << zeroGenerator((level*4+1)*chunkSize_-1,0) << " & " << carrys.str();
								// reinitialize the concatenation buffers for the next iteration of the while loop
								for (int j=0; j<2; j++)
								{
									vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[j].str() << endl;
									opCount++;
									operands[j].str("");
								}
								vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
								opCount++;	
								carrys.str(""); 
								break;
							default: // more than 2 tiles
								if (verbose)
									cerr << ">IntMultiplier:" << tab << " Subcase 3: there are more than 2 rows" << endl;
								
								operands[i%2] << zeroGenerator(level*4*chunkSize_,0) << ";";
								operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << ";";
								carrys << zeroGenerator((level*4+3)*chunkSize_,0) << ";";
								
								// top-right tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px0y" << i;
								vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i)) << ";" << endl;
								setCycle(3);
								operands[i%2].seekp(ios_base::beg);
								operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
								
								// top-right diagonal pair
								setCycle(1);
								partialProd.str("");
								partialProd << "px1y" << i;
								vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y",i)) << ";" << endl;
								partialProd2.str("");
								partialProd2 << "px0y" << i+1;
								vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i+1)) << ";" << endl;
								
								setCycle(2);
								sum.str("");
								sum << "addOpd" << level << "_shift_" << i+1;
								vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
								
								setCycle(3);
								operands[(i+1)%2].seekp(ios_base::beg);
								operands[(i+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
								carrys.seekp(ios_base::beg);
								carrys << "\"0\" & " << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << carrys.str(); 
								// sum of diagonal triplets from top to bottom
								while(i < chunksY-2)
								{
									setCycle(1);
									// first compute the 3 partial products
									for (int j=0; j<3; j++)
									{
										partialProd.str("");
										partialProd << "px" << j << "y" << i+2-j;
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",i+2-j)) << ";" << endl; 
									}
									setCycle(2);
									// add the 3 partial products
									sum.str("");
									sum << "addOp" << level << "_shift_" << i+2;
									vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
									
									for (int j=0; j<3; j++)
									{
										partialProd.str("");
										partialProd << "px" << j << "y" << i+2-j;
										if (j == 0)
											vhdl << use(partialProd.str());
										else
											vhdl << ") + (\"00\" & " << use(partialProd.str());
									}
									vhdl << ");" << endl;
									setCycle(3);
									operands[(i+2)%2].seekp(ios_base::beg);
									operands[(i+2)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+2)%2].str();
									carrys.seekp(ios_base::beg);
									carrys << use(sum.str()) << range(2*chunkSize_+1,2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
							
									i++;
								}
								
								// bottom-left diagonal pair
								setCycle(1);
								partialProd.str("");
								partialProd << "px2y" << chunksY-2;
								vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x2") << " * " << use(join("y",chunksY-2)) << ";" << endl;
								partialProd2.str("");
								partialProd2 << "px1y" << chunksY-1;
								vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y",chunksY-1)) << ";" << endl;
								
								setCycle(2);
								sum.str("");
								sum << "addOpd" << level << "_shift_" << chunksY;
								vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
								
								setCycle(3);
								operands[(i+2)%2].seekp(ios_base::beg);
								operands[(i+2)%2] << zeroGenerator((level*4+1)*chunkSize_,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+2)%2].str();
								carrys.seekp(ios_base::beg);
								carrys << zeroGenerator((level*4+1)*chunkSize_-1,0) << " & " << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str(); 
								// bottom-left tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px2y" << chunksY-1;
								vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x2") << " * " << use(join("y",chunksY-1)) << ";" << endl;
								setCycle(3);
								operands[(i+3)%2].seekp(ios_base::beg);
								operands[(i+3)%2] << zeroGenerator((level*4)*chunkSize_,0) << " & " << use(partialProd.str()) << " & " << operands[(i+3)%2].str();
								
								// create the operands for the final summation
								for (int j=0; j<2; j++)
								{
									vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[j].str() << endl;
									opCount++;
									operands[j].str("");
								}
								vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
								opCount++;	
								carrys.str(""); 
								break;
						}
						break;
				}
				setCycle(3);
			}
		/**********************************************************************************************************************/
			else // not(quadMultiply) => can use double chunckSize_ multipliers
			{	
				if (chunkSize_ == 18)
				{
					if (verbose)
						cerr << ">IntMultiplier: Block by block adition using double chunkSize multipliers" << endl;
					int i, j, k;	
					
					////////////////////////////////////////////////////
					//SPLITTINGS
					for (k=0; k<chunksX-1; k+=2)
					{
						ostringstream dname;
						dname << "x"<< k;
						vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range((k+2)*chunkSize_-1,k*chunkSize_) << ";" << endl;
					}
					
					/* The following splitting for y uses a small trick, i.e. if we have an even number of chunks for y
					 * nothing changes and the splitting procedes normally. Otherwise we skip over the first chunk y(17:0)
					 * and start counting the chunks from then on. At the end the last indexed chunk will be y(17:0).
					 * This was done to avoid major changes in the partial product computation and concatenation algorithm.*/
					int start = 0;
					
					if (chunksY%2 == 1)
						start = 1;
						
					for (k=0; k<chunksY-1 ; k+=2)
					{
						ostringstream dname;
						dname << "y"<< k;
						vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range((k+2+start)*chunkSize_-1,(k+start)*chunkSize_) << ";" << endl;
					}
					
					if (start == 1)
					{
						ostringstream dname;
						dname << "y"<< chunksY-1;
						vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range(chunkSize_-1,0) << ";" << endl;
					}
					
					nextCycle();
					
					// COMPUTE PARTIAL PRODUCTS
					for (i=0; i<chunksX-1; i+=2)
					{
						operands[0] << zeroGenerator((i+start)*chunkSize_,0) << ";";
						operands[1] << zeroGenerator((i+2+start)*chunkSize_,0) << ";";
						
						// compute and concatenate partial products from right to left
						for (k=0; k<chunksX-i-1; k+=2)
						{
							setCycleFromSignal(join("x",k));
							nextCycle();
							ostringstream partialProd;
							partialProd  << "p" << "x" << k << "y" << i;
							vhdl << tab << declare(partialProd.str(),4*chunkSize_) << " <= " << use(join("x",k)) << " * " << use(join("y",i)) << ";" << endl;
							
							nextCycle();
							operands[(k/2)%2].seekp(ios_base::beg);
							operands[(k/2)%2] << use(partialProd.str()) << " & " << operands[(k/2)%2].str();
						}
						
						// compute and concatenate partial products from top to bottom
						for (j=i+2; j<chunksY-1; j+=2)
						{
							setCycleFromSignal(join("x",k-2));
							nextCycle();
							ostringstream partialProd;
							partialProd  << "p" << "x" << k-2 << "y" << j;
							vhdl << tab << declare(partialProd.str(),4*chunkSize_) << " <= " << use(join("x",k-2)) << " * " << use(join("y",j)) << ";" << endl;
						
							nextCycle();
							operands[((j-i+k-2)/2)%2].seekp(ios_base::beg);
							operands[((j-i+k-2)/2)%2] << use(partialProd.str()) << " & " << operands[((j-i+k-2)/2)%2].str();
						}
						
						operands[((i+j+k-2)/2)%2].seekp(ios_base::beg);
						operands[((i+j+k-2)/2)%2] << zeroGenerator((i+2+chunksX%2)*chunkSize_,0) << " & " << operands[((i+j+k-2)/2)%2].str();	
						operands[((i+j+k-2)/2+1)%2].seekp(ios_base::beg);
						operands[((i+j+k-2)/2+1)%2] << zeroGenerator((i+chunksX%2)*chunkSize_,0) << " & " << operands[((i+j+k-2)/2+1)%2].str();	
						
						// reinitialize the concatenation buffers for the next iteration of the loop
						for (k=0; k<2; k++)
						{
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[k].str() << endl;
							opCount++;
							operands[k].str("");
						}
					}
					
					operands[0] << "\"\";";
					operands[1] << zeroGenerator(chunkSize_,0) << ";";
					
					bool halfPadded = false; /* TRUE when there is either 
												a chunkSize width column on the left side of the tiling or 
												a chunkSize width row on the bottom side of the tiling */	
					
					if (start == 1) // then there is a chunkSize width row on the bottom side of the tiling
					{
						halfPadded = true;
						
						// compute and concatenate partial products from right to left
						for (k=0; k<chunksX; k++)
						{
							setCycleFromSignal("sX");
							
							ostringstream dname;
							dname << "xx" << k;
							vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
							
							nextCycle();
							
							ostringstream partialProd;
							partialProd  << "p" << "xx" << k << "y" << chunksY-1;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("xx",k)) << " * " << use(join("y",chunksY-1)) << ";" << endl;
							
							nextCycle();
							
							operands[k%2].seekp(ios_base::beg);
							operands[k%2] << use(partialProd.str()) << " & " << operands[k%2].str();
						}
					}
					else // padd operands with zeros
					{
							operands[0].seekp(ios_base::beg);
							operands[0] << zeroGenerator((chunksX-1)*2*chunkSize_,0) << " & " << operands[0].str();
							operands[1].seekp(ios_base::beg);
							operands[1] << zeroGenerator((chunksX-1)*2*chunkSize_,0) << " & " << operands[1].str();
					}
					
					if (i == chunksX-1) // then there is a chunkSize width column on the left side of the tiling
					{
						halfPadded = true;
						setCycleFromSignal("sX");
						
						ostringstream dname;
						dname << "x" << chunksX-1;
						vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range(chunksX*chunkSize_-1,(chunksX-1)*chunkSize_) << ";" << endl;
						
						for (k=start; k<chunksY; k++)
						{
							setCycleFromSignal("sY");
							dname.str("");
							dname << "yy" << k;
							vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
						
							nextCycle();
							ostringstream partialProd;
							partialProd  << "p" << "x" << chunksX-1 << "yy" << k;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("x",chunksX-1)) << " * " << use(join("yy",k)) << ";" << endl;
						
							nextCycle();
							operands[(k+chunksX+1)%2].seekp(ios_base::beg);
							operands[(k+chunksX+1)%2] << use(partialProd.str()) << " & " << operands[(k+chunksX+1)%2].str();
						}
							
					}
					else // padd operands with zeros
					{
						operands[0].seekp(ios_base::beg);
						operands[0] << zeroGenerator((chunksY-2)*2*chunkSize_,0) << " & " << operands[0].str();
						operands[1].seekp(ios_base::beg);
						operands[1] << zeroGenerator((chunksY-2)*2*chunkSize_,0) << " & " << operands[1].str();	
					}
					
					if (halfPadded) // then we add the two operands to the sum
					{
						if ((chunksX+chunksY-1)%2 == 0)
						{
							operands[0].seekp(ios_base::beg); 
							operands[0] << zeroGenerator(chunkSize_,0) << " & " << operands[0].str(); 	
						}
						else
						{
							operands[1].seekp(ios_base::beg); 
							operands[1] << zeroGenerator(chunkSize_,0) << " & " << operands[1].str(); 	
						}
						
						for (k=0; k<2; k++)
						{
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[k].str() << endl;
							opCount++;
						}
					}
				}
				else // then chunkSize is 9 and we can use double chunkSize multipliers
				{
					if (verbose)
						cerr << ">IntMultiplier: Using double chunkSize_ multipliers" << endl;
					
					int i, j, k;	
					
					////////////////////////////////////////////////////
					//SPLITTINGS
					for (k=0; k<chunksX-1; k+=2)
					{
						ostringstream dname;
						dname << "x"<< k/2;
						vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range((k+2)*chunkSize_-1,k*chunkSize_) << ";" << endl;
					}
					
					/* The following splitting for y uses a small trick, i.e. if we have an even number of chunks for y
					 * nothing changes and the splitting procedes normally. Otherwise we skip over the first chunk y(17:0)
					 * and start counting the chunks from then on. At the end the last indexed chunk will be y(17:0).
					 * This was done to avoid major changes in the partial product computation and concatenation algorithm.*/
					int start = 0;
					
					if (chunksY%2 == 1)
						start = 1;
						
					for (k=0; k<chunksY-1 ; k+=2)
					{
						ostringstream dname;
						dname << "y"<< k/2;
						vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range((k+2+start)*chunkSize_-1,(k+start)*chunkSize_) << ";" << endl;
					}
					
					if (start == 1)
					{
						ostringstream dname;
						dname << "y"<< k/2;
						vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range(chunkSize_-1,0) << ";" << endl;
					}
					int startX = 0; // added this variable to apply the same trick for zeroGenerators as with the start variable that adds zeros behind if chunksY is odd 
					bool leftColumn = false; // TRUE if there are an odd number of slices in X and there is a chunkSize width column to the left of the tiling that has no pair
					if (chunksX % 2 == 1)
					{
						leftColumn = true;
						startX = 1;
					}
					
					/* modifing the quantities so that we can use the same algorithm as for the case where we have no
					 * double chunkSize multipliers */
					chX /= 2;
					chunksY /= 2;
					chunksX /= 2;
					chunkSize_ *=2;
					
					while (chX/4 > 0)
					{
						
						// top-right tile
						setCycle(1);
						partialProd.str("");
						partialProd << "px0y" << 4*level;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y", 4*level)) << ";" << endl;
						setCycle(3);
						operands[0] << use(partialProd.str()) << " & " << zeroGenerator(chunkSize_*4*level + start*chunkSize_/2,0) << ";" << endl;
						
						
						// top-right diagonal pair
						setCycle(1);
						partialProd2.str("");
						partialProd2 << "px0y" << 4*level+1;
						vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",4*level+1)) << ";" << endl;
						partialProd.str("");
						partialProd << "px1y" << 4*level;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y",4*level)) << ";" << endl;
						setCycle(2);
						sum.str("");
						sum << "addOp" << level << "_shift_" << 4*level+1;
						vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
						setCycle(3);
						operands[1] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << zeroGenerator((4*level+1)*chunkSize_ + start*chunkSize_/2,0) << ";" << endl;
						carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zeroGenerator((4*level+3)*chunkSize_ + start*chunkSize_/2,0) << ";" << endl;
						
						
						// top-right diagonal triplet
						setCycle(1);
						// first compute the 3 partial products
						for (int j=0; j<3; j++)
						{
							partialProd.str("");
							partialProd << "px" << j << "y" << 4*level+2-j;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",4*level+2-j)) << ";" << endl; 
						}
						setCycle(2);
						// add the 3 partial products
						sum.str("");
						sum << "addOp" << level << "_shift_" << 4*level+2;
						vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
						
						for (int j=0; j<3; j++)
						{
							partialProd.str("");
							partialProd << "px" << j << "y" << 4*level+2-j;
							if (j == 0)
								vhdl << use(partialProd.str());
							else
								vhdl << ") + (\"00\" & " << use(partialProd.str());
						}
						vhdl << ");" << endl;
						setCycle(3);
						operands[0].seekp(ios_base::beg);
						operands[0] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[0].str();
						carrys.seekp(ios_base::beg);
						carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-1,0) << " & " << carrys.str();
						
							
						//setCycle(0);
							
						// compute the sum of four tiles (diagonal) starting from top-right going left
						for (int i=3; i<chX-1; i++)
						{
							setCycle(1);
							// first compute the 4 partial products
							for (int j=0; j<4; j++)
							{
								partialProd.str("");
								partialProd << "px" << i-j << "y" << j+4*level;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",i-j)) << " * " << use(join("y",j+4*level)) << ";" << endl;
							}
							setCycle(2);
							// add the 4 partial products
							sum.str("");
							sum << "addOp" << level << "_shift_" << i+4*level;
							vhdl << tab << declare(sum.str(),2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
							
							for (int j=0; j<4; j++)
							{
								partialProd.str("");
								partialProd << "px" << i-j << "y" << j+4*level;
								if (j == 0)
									vhdl << use(partialProd.str());
								else
									vhdl << ") + (\"00\" & " << use(partialProd.str());
							}
							vhdl << ");" << endl;
							setCycle(3);
							operands[i%2].seekp(ios_base::beg);
							operands[i%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[i%2].str();  
							carrys.seekp(ios_base::beg);
							carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
						}
						
						// compute the sum of four tiles (diagonal) going down to the bottom-left corner
						for (int j=4*level; j<chunksY-3; j++)
						{
							setCycle(1);
							// first compute the 4 partial products
							for (int i=0; i<4; i++)
							{
								partialProd.str("");
								partialProd << "px" << chX-i-1 << "y" << i+j;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-i-1)) << " * " << use(join("y",i+j)) << ";" << endl; 
							}
							setCycle(2);
							// add the 4 partial products
							sum.str("");
							sum << "addOp" << level << "_shift_" << chX+j-1;
							vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
							
							for (int i=0; i<4; i++)
							{
								partialProd.str("");
								partialProd << "px" << chX-i-1 << "y" << i+j;
								if (i == 0)
									vhdl << use(partialProd.str());
								else
									vhdl << ") + (\"00\" & " << use(partialProd.str());
							}
							vhdl << ");" << endl;
							setCycle(3);
							operands[(j+chX-1)%2].seekp(ios_base::beg);
							operands[(j+chX-1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(j+chX-1)%2].str(); 
							carrys.seekp(ios_base::beg);
							carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
						} 
						
						// bottom-left diagonal triplet
						setCycle(1);
						// first compute the 3 partial products
						for (int j=0; j<3; j++)
						{
							partialProd.str("");
							partialProd << "px" << chX-3+j << "y" << chunksY-1-j;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-3+j)) << " * " << use(join("y",chunksY-1-j)) << ";" << endl; 
						}
						setCycle(2);
						// add the 3 partial products
						sum.str("");
						sum << "addOp" << level << "_shift_" << chunksY+chX-4;
						vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
						
						for (int j=0; j<3; j++)
						{
							partialProd.str("");
							partialProd << "px" << chX-3+j << "y" << chunksY-1-j;
							if (j == 0)
								vhdl << use(partialProd.str());
							else
								vhdl << ") + (\"00\" & " << use(partialProd.str());
						}
						vhdl << ");" << endl;
						setCycle(3);
						operands[(chX+chunksY-4)%2].seekp(ios_base::beg);
						operands[(chX+chunksY-4)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(chX+chunksY-4)%2].str();
						carrys.seekp(ios_base::beg);
						carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
						
						// bottom-left diagonal pair
						setCycle(1);
						partialProd2.str("");
						partialProd2 << "px" << chX-1 << "y" << chunksY-2;
						vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-1)) << " * " << use(join("y", chunksY-2)) << ";" << endl;
						partialProd.str("");
						partialProd << "px" << chX-2 << "y" << chunksY-1;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x", chX-2)) << " * " << use(join("y", chunksY-1)) << ";" << endl;
						setCycle(2);
						sum.str("");
						sum << "addOp" << level << "_shift_" << chX+chunksY-3;
						vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
						setCycle(3);
						operands[(chX+chunksY-3)%2].seekp(ios_base::beg);
						operands[(chX+chunksY-3)%2] << zeroGenerator((4*level+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(chX+chunksY-3)%2].str();
						carrys.seekp(ios_base::beg);
						carrys << zeroGenerator((4*level+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
						
						// bottom-left tile
						setCycle(1);
						partialProd.str("");
						partialProd << "px"<< chX-1 << "y" << chunksY-1;
						vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-1)) << " * " << use(join("y", chunksY-1)) << ";" << endl;
						setCycle(3);
						operands[(chX+chunksY-2)%2].seekp(ios_base::beg);
						operands[(chX+chunksY-2)%2] << zeroGenerator(4*level*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << " & " << operands[(chX+chunksY-2)%2].str();
						
						// reinitialize the concatenation buffers for the next iteration of the while loop
						for (int i=0; i<2; i++)
						{
							setCycle(3);
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[i].str() << endl;
							opCount++;
							operands[i].str("");
						}
						vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
						opCount++;	
						carrys.str("");
						
						chX -= 4;
						level++;
					}
					
					i = level*4;	
					bool oneCol = true;		
					switch (chX)
					{
						case 1: // remaining tiles are in a single column
							if (verbose)
								cerr << ">IntMultiplier: Case 1: remaining tiles are in a single column" << endl;
							
							operands[i%2] << zeroGenerator(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
							if (chunksY-i > 1) // then there are more than one tiles in the column
							{
								if (verbose)
									cerr << ">IntMultiplier: More than one tile in the column" << endl;
								operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
								oneCol = false;
							}
							
							for (i=level*4; i<chunksY; i++)
							{
								setCycle(1);
								partialProd.str("");
								partialProd << "px0y" << i;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i)) << ";" << endl;
								
								setCycle(3);
								operands[i%2].seekp(ios_base::beg);
								operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
								
							}
							
							if (!oneCol) // then there are more than one tiles in the column
							{
								operands[i%2].seekp(ios_base::beg);
								operands[i%2] << zeroGenerator((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[i%2].str();
								vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[i%2].str() << endl;
								opCount++;
							}
							
							operands[(i+1)%2].seekp(ios_base::beg);
							operands[(i+1)%2] << zeroGenerator(level*4*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(i+1)%2].str(); 
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(i+1)%2].str() << endl;
							opCount++;
								
							break;
						case 2: // remaining tiles are in 2 columns
							if (verbose)
								cerr << ">IntMultiplier: Case 2: remaining tiles are in 2 columns" << endl;
							setCycle(1);
							// top right tile
							partialProd.str("");
							partialProd << "px0y" << i;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y", i)) << ";" << endl;
							setCycle(3);
							operands[i%2] << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << zeroGenerator((4*level)*chunkSize_ + start*chunkSize_/2,0) << ";";
							operands[(i+1)%2] << zeroGenerator((4*level+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
							carrys << zeroGenerator((4*level+2)*chunkSize_+1 + start*chunkSize_/2,0) << ";";
							 
							while (i<chunksY-1) // diagonal pairs
							{
								setCycle(1);
								partialProd2.str("");
								partialProd2 << "px0y" << (i+1);
								vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y", i+1)) << ";" << endl;
								partialProd.str("");
								partialProd << "px1y" << i;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y", i)) << ";" << endl;
								setCycle(2);
								sum.str("");
								sum << "addOp" << level << "_shift_" << (i+1);
								vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
								setCycle(3);
								operands[(i+1)%2].seekp(ios_base::beg);
								operands[(i+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
								carrys.seekp(ios_base::beg);
								carrys << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zeroGenerator(chunkSize_-1,0) << " & " << carrys.str();
						
								i++;			
							}
							
							// bottom-left tile
							setCycle(1);
							partialProd.str("");
							partialProd << "px1y" << i;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y", i)) << ";" << endl;
							
							setCycle(3);
							operands[(i+1)%2].seekp(ios_base::beg);
							operands[(i+1)%2] << zeroGenerator((4*level)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
							
							operands[i%2].seekp(ios_base::beg);
							operands[i%2] << zeroGenerator((4*level+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[i%2].str();
							carrys.seekp(ios_base::beg);
							carrys << zeroGenerator((4*level+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << carrys.str();
						
							// create operands for final sumation
							for (int j=0; j<2; j++)
							{
								setCycle(3);
								vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[j].str() << endl;
								opCount++;
							}
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
							opCount++;	
							break;
						case 3: // Remaining tiles are on 3 columns
							if (verbose)
								cerr << ">IntMultiplier: Case 3: Remaining tiles are on 3 columns" << endl;
							
							height = chunksY - i;
							
							switch (height)
							{
								case 0: break;
								case 1: // only one row
									if (verbose)
										cerr << ">IntMultiplier: " << tab << "Subcase 1: there is one row" << endl;
								
									operands[i%2] << zeroGenerator(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
									operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
							
									for (int j=0; j<3; j++) // each tile is computed and concatenated with an coresponding operand
									{
										setCycle(1);
										partialProd.str("");
										partialProd << "px" << j << "y" << i;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",i)) << ";" << endl;
										setCycle(3);
										operands[(j+i)%2].seekp(ios_base::beg);
										operands[(j+i)%2] << use(partialProd.str()) << " & " << operands[(j+i)%2].str();
									}
									
									operands[(3+i)%2].seekp(ios_base::beg);
									operands[(3+i)%2] << zeroGenerator((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(3+i)%2].str();
									vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(3+i)%2].str() << endl;
									opCount++;
									
									operands[(i+4)%2].seekp(ios_base::beg);
									operands[(i+4)%2] << zeroGenerator(level*4*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(i+4)%2].str(); 
									vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(i+4)%2].str() << endl;
									opCount++;
									break;
								case 2: // two rows
									if (verbose)
										cerr << ">IntMultiplier: " << tab << "Subcase 2: there are two rows" << endl;
									
									operands[i%2] << zeroGenerator(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
									operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
									carrys << zeroGenerator(level*4*chunkSize_,0) << ";";
									// top-right tile
									setCycle(1);
									partialProd.str("");
									partialProd << "px0y" << i;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i)) << ";" << endl;
									operands[i%2].seekp(ios_base::beg);
									operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
									
									// diagonal tiles
									for (int j=0; j<2; j++)
									{
										setCycle(1);
										partialProd.str("");
										partialProd << "px"<< j+1 << "y" << i;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j+1)) << " * " << use(join("y",i)) << ";" << endl;
										partialProd2.str("");
										partialProd2 << "px"<< j << "y" << i+1;
										vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",i+1)) << ";" << endl;
										
										setCycle(2);
										sum.str("");
										sum << "addOpd" << level << "_shift_" << i+j+1;
										vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
										
										setCycle(3);
										operands[(i+j+1)%2].seekp(ios_base::beg);
										operands[(i+j+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) <<" & " << operands[(i+j+1)%2].str();
										carrys.seekp(ios_base::beg);
										carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zeroGenerator(chunkSize_,0) << "&" << carrys.str();
									}
									
									// bottom-left tile
									setCycle(1);
									partialProd.str("");
									partialProd << "px2y" << i+2;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x2") << " * " << use(join("y",i+2)) << ";" << endl;
									setCycle(3);
									operands[(i+4)%2].seekp(ios_base::beg);
									operands[(i+4)%2] << zeroGenerator(level*4*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << " & " << operands[(i+4)%2].str();	
									operands[(i+5)%2].seekp(ios_base::beg);
									operands[(i+5)%2] << zeroGenerator((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(i+5)%2].str(); 
									carrys.seekp(ios_base::beg);
									carrys << zeroGenerator((level*4+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << carrys.str();
									// create the operands for the final summation
									for (int j=0; j<2; j++)
									{
										vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[j].str() << endl;
										opCount++;
									}
									vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
									opCount++;	
									break;
								default: // more than 2 tiles
									if (verbose)
										cerr << ">IntMultiplier: " << tab << "Subcase 3: there are more than 2 rows" << endl;
									
									operands[i%2] << zeroGenerator(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
									operands[(i+1)%2] << zeroGenerator((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
									carrys << zeroGenerator((level*4+3)*chunkSize_ + start*chunkSize_/2,0) << ";";
									
									// top-right tile
									setCycle(1);
									partialProd.str("");
									partialProd << "px0y" << i;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i)) << ";" << endl;
									setCycle(3);
									operands[i%2].seekp(ios_base::beg);
									operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
									
									// top-right diagonal pair
									setCycle(1);
									partialProd.str("");
									partialProd << "px1y" << i;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y",i)) << ";" << endl;
									partialProd2.str("");
									partialProd2 << "px0y" << i+1;
									vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x0") << " * " << use(join("y",i+1)) << ";" << endl;
									
									setCycle(2);
									sum.str("");
									sum << "addOpd" << level << "_shift_" << i+1;
									vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
									
									setCycle(3);
									operands[(i+1)%2].seekp(ios_base::beg);
									operands[(i+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
									carrys.seekp(ios_base::beg);
									carrys << "\"0\" & " << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << carrys.str(); 
									// sum of diagonal triplets from top to bottom
									while(i < chunksY-2)
									{
										setCycle(1);
										// first compute the 3 partial products
										for (int j=0; j<3; j++)
										{
											partialProd.str("");
											partialProd << "px" << j << "y" << i+2-j;
											vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",j)) << " * " << use(join("y",i+2-j)) << ";" << endl; 
										}
										setCycle(2);
										// add the 3 partial products
										sum.str("");
										sum << "addOp" << level << "_shift_" << i+2;
										vhdl << tab << declare(sum.str(), 2*chunkSize_+2,true,Signal::registeredWithAsyncReset) << " <= (\"00\" & ";
										
										for (int j=0; j<3; j++)
										{
											partialProd.str("");
											partialProd << "px" << j << "y" << i+2-j;
											if (j == 0)
												vhdl << use(partialProd.str());
											else
												vhdl << ") + (\"00\" & " << use(partialProd.str());
										}
										vhdl << ");" << endl;
										setCycle(3);
										operands[(i+2)%2].seekp(ios_base::beg);
										operands[(i+2)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+2)%2].str();
										carrys.seekp(ios_base::beg);
										carrys << use(sum.str()) << range(2*chunkSize_+1,2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str();
								
										i++;
									}
									
									// bottom-left diagonal pair
									setCycle(1);
									partialProd.str("");
									partialProd << "px2y" << chunksY-2;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x2") << " * " << use(join("y",chunksY-2)) << ";" << endl;
									partialProd2.str("");
									partialProd2 << "px1y" << chunksY-1;
									vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x1") << " * " << use(join("y",chunksY-1)) << ";" << endl;
									
									setCycle(2);
									sum.str("");
									sum << "addOpd" << level << "_shift_" << chunksY;
									vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
									
									setCycle(3);
									operands[(i+2)%2].seekp(ios_base::beg);
									operands[(i+2)%2] << zeroGenerator((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+2)%2].str();
									carrys.seekp(ios_base::beg);
									carrys << zeroGenerator((level*4+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zeroGenerator(chunkSize_-2,0) << " & " << carrys.str(); 
									// bottom-left tile
									setCycle(1);
									partialProd.str("");
									partialProd << "px2y" << chunksY-1;
									vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use("x2") << " * " << use(join("y",chunksY-1)) << ";" << endl;
									setCycle(3);
									operands[(i+3)%2].seekp(ios_base::beg);
									operands[(i+3)%2] << zeroGenerator((level*4)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << " & " << operands[(i+3)%2].str();
									
									// create the operands for the final summation
									for (int j=0; j<2; j++)
									{
										vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[j].str() << endl;
										opCount++;
										operands[j].str("");
									}
									vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << carrys.str() << endl;
									opCount++;	
									carrys.str("");
									break;
							}
							break;
					}
					
					// REESTABLISH THE INITIAL VALUES FOR THESE VARIABLES
					chunksY = chunksY*2 + start;
					chunksX = chunksX*2 + leftColumn;
					chunkSize_ /=2;
					// HANDLE THE CHUNK SIZE PADDING IF NECESSARY
					operands[0].str("");
					operands[0] << "\"\";";
					operands[1].str("");
					operands[1] << zeroGenerator(chunkSize_,0) << ";";
					
					bool halfPadded = false; /* TRUE when there is either 
												a chunkSize width column on the left side of the tiling or 
												a chunkSize width row on the bottom side of the tiling */	
					
					if (start == 1) // then there is a chunkSize width row on the top side of the tiling
					{
						halfPadded = true;
						
						// compute and concatenate partial products from right to left
						for (k=0; k<chunksX; k++)
						{
							//setCycleFromSignal("sX");
							setCycle(0);
							ostringstream dname;
							dname << "xx" << k;
							vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
							
							//nextCycle();
							setCycle(1);
							ostringstream partialProd;
							partialProd  << "p" << "xx" << k << "y" << chunksY/2;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("xx",k)) << " * " << use(join("y",chunksY/2)) << ";" << endl;
							
							//nextCycle();
							setCycle(3);
							operands[k%2].seekp(ios_base::beg);
							operands[k%2] << use(partialProd.str()) << " & " << operands[k%2].str();
						}
						
					}
					else // padd operands with zeros
					{
							operands[0].seekp(ios_base::beg);
							operands[0] << zeroGenerator((chunksX-1)*chunkSize_,0) << " & " << operands[0].str();
							operands[1].seekp(ios_base::beg);
							operands[1] << zeroGenerator((chunksX-1)*chunkSize_,0) << " & " << operands[1].str();
					}
					
					if (leftColumn) // then there is a chunkSize width column on the left side of the tiling
					{
						halfPadded = true;
						//setCycleFromSignal("sX");
						setCycle(0);
						ostringstream dname;
						dname << "x" << chunksX-1;
						vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sX" << range(chunksX*chunkSize_-1,(chunksX-1)*chunkSize_) << ";" << endl;
						
						for (k=start; k<chunksY; k++)
						{
							//setCycleFromSignal("sY");
							setCycle(0);
							dname.str("");
							dname << "yy" << k;
							vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << "sY" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
						
							//nextCycle();
							setCycle(1);
							ostringstream partialProd;
							partialProd  << "p" << "x" << chunksX-1 << "yy" << k;
							vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("x",chunksX-1)) << " * " << use(join("yy",k)) << ";" << endl;
						
							//nextCycle();
							setCycle(3);
							operands[(k+chunksX+1)%2].seekp(ios_base::beg);
							operands[(k+chunksX+1)%2] << use(partialProd.str()) << " & " << operands[(k+chunksX+1)%2].str();
						}
							
					}
					else // padd operands with zeros
					{
						operands[0].seekp(ios_base::beg);
						operands[0] << zeroGenerator((chunksY-1)*chunkSize_,0) << " & " << operands[0].str();
						operands[1].seekp(ios_base::beg);
						operands[1] << zeroGenerator((chunksY-1)*chunkSize_,0) << " & " << operands[1].str();	
					}
					
					//setCycle(3);
					
					if (halfPadded) // then we add the two operands to the sum
					{
						if ((chunksX+chunksY-1)%2 == 0)
						{
							operands[0].seekp(ios_base::beg); 
							operands[0] << zeroGenerator(chunkSize_,0) << " & " << operands[0].str(); 	
						}
						else
						{
							operands[1].seekp(ios_base::beg); 
							operands[1] << zeroGenerator(chunkSize_,0) << " & " << operands[1].str(); 	
						}
						
						for (k=0; k<2; k++)
						{
							vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[k].str() << endl;
							opCount++;
						}
					}
				}
			}	
			
			IntNAdder* add =  new IntNAdder(target, adderWidth, opCount);
			//IntCompressorTree* add =  new IntCompressorTree(target, adderWidth, opCount);
			oplist.push_back(add);


			for (int j=0; j<opCount; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOp" << j;
				
				inPortMap (add, join("X",j) , concatPartialProd.str());
			}	
					
			inPortMapCst(add, "Cin", "'0'");
			outPortMap(add, "R", "addRes");
			vhdl << instance(add, "adder");

			syncCycleFromSignal("addRes");
			
			vhdl << tab << "R <= " << use("addRes")<<range(adderWidth-1,adderWidth-wInX_-wInY_) << ";" << endl;	
			
		}
	}
	else
	{
		int chunkSize_ = target->lutInputs()/2;
	
	int chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
	int chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
		
	if (verbose)
        cerr << "> IntMultiplier:  X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; " << endl;
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

