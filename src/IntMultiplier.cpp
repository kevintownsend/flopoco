/*
  An integer multiplier for FloPoCo
  
  Authors: Bogdan Pasca, Sebastian Banescu

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

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
#include "IntMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	IntMultiplier:: IntMultiplier(Target* target, int wInX, int wInY, map<string, double> inputDelays, bool sign) :
		Operator(target, inputDelays), wInX_(wInX), wInY_(wInY), wOut_(wInX + wInY), sign_(sign){
 
		ostringstream name;
		double target_freq = target->frequency();
		double maxDSPFrequency = int(floor(1.0/ (target->DSPMultiplierDelay() + 1.0e-10)));
		if (target_freq > maxDSPFrequency)
			target->setFrequency(maxDSPFrequency);
			

		/* Name Setup procedure
		 *  The name has the format: IntMultiplier_wInX_wInY
		 *  wInX = width of the X input
		 *  wInY = width of the Y input
		 */  
		name <<"IntMultiplier_"<<wInX_<<"_"<<wInY_;
		if (sign_)
			name << "_signed";
		srcFileName="IntMultiplier";
		setName(name.str());
		
		REPORT(INFO, "################################################################################");
		REPORT(INFO, " ------------------- wInX="<<wInX<<" wInY="<<wInY<<" sign="<<sign<<" -------------------");
	
		setCopyrightString("Bogdan Pasca, Sebastian Banescu (2008-2009)");
	
		addInput ("X", wInX_,true);
		addInput ("Y", wInY_,true);
		addOutput("R", wOut_,1, true); /* wOut_ = wInX_ + wInY_ */
	
		double delay = 0.0;
	
		if ((target->getUseHardMultipliers()) && (target->getNumberOfDSPs()>0)){
			if (verbose)
				cerr << "> IntMultiplier: The target is " << target->getID() << endl;
	
			if (target->getVendor()=="Xilinx"){
				int chunksX, chunksY;
				int x, y;
				int score1, score2;

				int explicitSwap;
				target->suggestSubmultSize(x, y, wInX, wInY);

				score1 = int(ceil((double(wInX)/double(x)))*ceil((double(wInY)/double(y))));
				score2 = int(ceil((double(wInY)/double(x)))*ceil((double(wInX)/double(y))));

				explicitSwap = false;
				if (score1 > score2)
					explicitSwap = true;

				int extension = (y<x)?(-x):(-y); // for asymmetric sized multipliers (e.g. 24x17 on Virtex5)  

				if (score1 <= score2){
					chunksX =  int(ceil( ( double(wInX) / (double) x) ));
					chunksY =  int(ceil( ( double(wInY) / (double) y) ));
				}
				else{
					chunksX =  int(ceil( ( double(wInX) / (double) y) ));
					chunksY =  int(ceil( ( double(wInY) / (double) x) ));
				}

				if (chunksX + chunksY > 2) { // up to 17 x 17 bit on Virtex4 can be written as an "*" @ 400++ MHz 
					// to be general version (parametrized etc) 
					//TODO swap X and Y under certain conditions

					//NOT TOO SMART
					bool swap=false;
					if (((chunksX > chunksY) && (x == y)) || (explicitSwap == true))
						swap=true;
		
					if (swap){
						int tmp = chunksX;
						chunksX = chunksY;
						chunksY = tmp; 
					}
		
					if (verbose)
						cerr << "> IntMultiplier: Perform swapping = " << (swap==1?"true":"false") << endl;
					
					REPORT(INFO, "x=" << x << " y= " << y << " chunksX="<<chunksX << " chunksY="<<chunksY);
		
					if (swap){
						vhdl << tab << declare("sX",x*chunksX) << " <= Y & " << zg(x*chunksX-wInY,0) << ";" << endl;
						vhdl << tab << declare("sY",y*chunksY) << " <= X & " << zg(y*chunksY-wInX,0) << ";" << endl;
					}else{
						vhdl << tab << declare("sX",x*chunksX) << " <= X & " << zg(x*chunksX-wInX,0) << ";" << endl;
						vhdl << tab << declare("sY",y*chunksY) << " <= Y & " << zg(y*chunksY-wInY,0) << ";" << endl;
					}
					////////////////////////////////////////////////////
					//SPLITTINGS
				
					setCriticalPath( getMaxInputDelays(inputDelays));
				
					for (int k=0; k<chunksX ; k++)
						vhdl << tab << declare(join("x",k),x) << " <= sX" << range((k+1)*x-1,k*x) << ";" << endl;
					for (int k=0; k<chunksY ; k++)
						vhdl << tab << declare(join("y",k),y) << " <= sY" << range((k+1)*y-1,k*y) << ";" << endl;

					//MULTIPLICATIONS WITH SOME ACCUMULATIONS
					for (int i=0; i<chunksX; i++){ 
						setCriticalPath ( getMaxInputDelays(inputDelays) );
						for (int j=0; j<chunksY; j++){ 
							if (j==0){ // @ the first the operation is only multiplication, not MAC
								manageCriticalPath(target->DSPMultiplierDelay());
								vhdl << tab << declare(join("px",i,"y",j),x+y) << " <= " << use(join("x",i)) << " * " << use(join("y",j)) << ";" << endl;
							}else{
								manageCriticalPath(target->DSPCascadingWireDelay() + target->DSPAdderDelay());
								vhdl << tab << declare(join("tpx",i,"y",j),x+y) << " <= " << use(join("x",i)) << " * " << use(join("y",j))  << ";" << endl; 
								vhdl << tab << declare(join("px",i,"y",j),x+y+1) << " <= ( \"0\" & " << use(join("tpx",i,"y",j)) << ") + " << use(join("px",i,"y",j-1))<<range(x+y-1,y) << ";" << endl; 
							}
	//									if (!((j==chunksY-1) && (i<chunksX-1))) nextCycle();
						}
						if (i<chunksX-1) setCycle(0); //reset cycle
					}

					manageCriticalPath( 2*target->DSPToLogicWireDelay() + target->DSPCascadingWireDelay() );
		
					//FORM THE INTERMEDIARY PRODUCTS
					for (int i=0; i<chunksX ; i++){
						vhdl << tab << declare(join("sum",i),y*chunksY+x) << " <= "; 
						for (int j=chunksY-1;j>=0; j--)
							vhdl << use(join("px",i,"y",j)) << range ( (j==chunksY-1?(x+y-1):y-1) ,0) << (j==0 ? ";" : " & ");
						vhdl << endl;
					}
					vhdl << tab << declare ("sum0Low", x) << " <= sum0" << range(x-1,0) << ";" << endl;
		
					if (chunksX>1){
						manageCriticalPath(target->DSPToLogicWireDelay());
				
						REPORT(DEBUG, "delay at adder input " << getCriticalPath() ); 
	//								IntNAdder* add =  new IntNAdder(target, x*chunksX+y*chunksY+extension, chunksX, inDelayMap("X0",getCriticalPath()));
						IntCompressorTree* add =  new IntCompressorTree(target, x*chunksX+y*chunksY+extension, chunksX, inDelayMap("X0",getCriticalPath()));

						oplist.push_back(add);
	
						for (int i=0; i< chunksX; i++){
							if (i==0) vhdl << tab << declare (join("addOp",i),x*chunksX+y*chunksY+extension) << " <= " << zg(x*chunksX+extension,0) << " & " << use(join("sum",i)) << range(y*chunksY+x-1,x) << ";" <<endl;
							if (i==1) vhdl << tab << declare (join("addOp",i),x*chunksX+y*chunksY+extension) << " <= " << zg(x*chunksX+extension-x,0) << " & " << use(join("sum",i)) << range(y*chunksY+x-1,0) << ";" <<endl;
							if (i> 1) 
								vhdl << tab << declare (join("addOp",i),x*chunksX+y*chunksY+extension) << " <= " << zg(x*chunksX-x*i+extension,0) << " & " << use(join("sum",i)) << range(y*chunksY+x-1,0) << " & " << zg(x*i-x,0) << ";" <<endl;
						}
						for (int i=0; i< chunksX; i++)
							inPortMap (add, join("X",i) , join("addOp",i));
	
	//								inPortMapCst(add, "Cin", "'0'");
						outPortMap(add, "R", "addRes");
						vhdl << instance(add, "adder");

						syncCycleFromSignal("addRes");

						outDelayMap["R"] = (add->getOutDelayMap())["R"];
									
						if ((x*chunksX+y*chunksY) - wInX - wInY < -extension){
							vhdl << tab << "R <= addRes" << range((x*chunksX+y*chunksY+extension)-1,0) << " & sum0Low" << range(x-1, x-1 + ((x*chunksX+y*chunksY+extension) - wInX - wInY)+1) << ";" << endl;
						}else{
							vhdl << tab << "R <= addRes" << range((x*chunksX+y*chunksY+extension)-1, ((x*chunksX+y*chunksY+extension) - wInX - wInY)) << ";" << endl;
	//									outDelayMap["R"] = 0.0;
						}
					}else{
						vhdl << tab << "R <= " << use(join("sum",0))<<range(y*chunksY+x-1, x*chunksX+y*chunksY -( wInX+wInY)) << ";" << endl;
						outDelayMap["R"] = getCriticalPath();
					}
				}
				else{
					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPMultiplierDelay());
					vhdl << tab << "R <= X * Y ;" <<endl;
				
					outDelayMap["R"] = 1.0/ target->frequency();
				}
		}else if (target->getVendor()=="Altera"){
			//Bogdan version of the altera multiplier.
			string opX,opY;
			int tmp;
			if (wInX < wInY){
				opX = "Y";
				opY = "X";
				tmp = wInX;
				wInX = wInY;
				wInY = tmp;
			}else{
				opX = "X";
				opY = "Y";
			}
			
			if ((wInX <= 18) && (wInY <= 18)){
				setCriticalPath( getMaxInputDelays(inputDelays));
				manageCriticalPath( target->DSPToLogicWireDelay());
				manageCriticalPath( target->DSPMultiplierDelay());
				vhdl << tab << "R <= X * Y;" << endl;
				manageCriticalPath( target->DSPToLogicWireDelay() );
				outDelayMap["R"]=getCriticalPath();
			}else if ((wInX <= 36) && (wInY <= 36)){				
				setCriticalPath( getMaxInputDelays(inputDelays));
				manageCriticalPath( target->DSPToLogicWireDelay());
				manageCriticalPath( target->DSPMultiplierDelay());
//						manageCriticalPath( target->DSPAdderDelay());
//						manageCriticalPath( target->DSPAdderDelay());
				vhdl << tab << "R <= X * Y;" << endl;
				manageCriticalPath( target->DSPToLogicWireDelay() );
				outDelayMap["R"]=getCriticalPath();
			}else if ((wInX <= 54) && (wInY <= 36)){
				//first multiplication
//							setCriticalPath( getMaxInputDelays(inputDelays));
//							manageCriticalPath( target->DSPToLogicWireDelay());
//							manageCriticalPath( target->DSPMultiplierDelay());
				
					vhdl << tab << declare("x0",36) << " <= "<< opX << range(35,0) << ";" << endl;
					vhdl << tab << declare("y0",36) << " <= "<< opY << range(35,0) << ";" << endl;

					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPMultiplierDelay() + target->DSPToLogicWireDelay());
//							manageCriticalPath( target->DSPAdderDelay());
//							manageCriticalPath( target->DSPAdderDelay());
					vhdl << tab << declare("r0",72) << "<= x0 * y0;" << endl;
				
				setCycle(0); //reset for second multiplication
					vhdl << tab << declare("x1",wInX-36) << " <= "<< opX << range(wInX-1,36) << ";" << endl;
					vhdl << tab << declare("y1",18) << " <= "<< opY << range(17,0) << ";" << endl;
					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPMultiplierDelay() + target->DSPToLogicWireDelay());
					vhdl << tab << declare("r1",18 + (wInX-36)) << "<= x1 * y1;" << endl;
					
				setCycle(0); //reset for second multiplication
					vhdl << tab << declare("y2",wInY-18) << " <= "<< opY << range(wInY-1,18) << ";" << endl;
					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPMultiplierDelay() + target->DSPToLogicWireDelay());
					vhdl << tab << declare("r2",(wInY-18) + (wInX-36)) << "<= x1 * y2;" << endl;
					
				//syncronization
					syncCycleFromSignal("r0");	
					syncCycleFromSignal("r1");
					syncCycleFromSignal("r2");
					
//							manageCriticalPath( target->DSPToLogicWireDelay() );
					nextCycle();///////////////
					vhdl << tab << declare("addOp1",18 + (wInY-18) + (wInX-36)) << " <= "
					            << zg(18 + (wInY-18) + (wInX-36) - 36,0) << " & " << "r0" << range(71, 36) << ";" << endl;
					vhdl << tab << declare("addOp2",18 + (wInY-18) + (wInX-36)) << " <= "
					            << zg(18 + (wInY-18) + (wInX-36) - (18 + (wInX-36)),0) << " & " << "r1;" << endl;
					vhdl << tab << declare("addOp3",18 + (wInY-18) + (wInX-36)) << " <= "
					            << "r2 & " << zg(18,0) <<";" << endl;
					
					manageCriticalPath( target->adderDelay(18 + (wInY-18) + (wInX-36)) );							
					vhdl << tab << declare("addRes", 18 + (wInY-18) + (wInX-36)) << " <= addOp1 + addOp2 + addOp3;"<<endl;
					
					vhdl << tab << "R <= addRes & r0"<<range(35,0) << ";" << endl;
					outDelayMap["R"]=getCriticalPath();
											
			}	else if ((wInX <= 54) && (wInY <= 54)){
				//first multiplication
					vhdl << tab << declare("x0",36) << " <= "<< opX << range(35,0) << ";" << endl;
					vhdl << tab << declare("y0",36) << " <= "<< opY << range(35,0) << ";" << endl;

					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPToLogicWireDelay() + target->DSPMultiplierDelay());
					vhdl << tab << declare("r00",72) << "<= x0 * y0;" << endl;

				///////////////////////////////////////////////////////						
				setCycle(0); //reset for second multiplication
					vhdl << tab << declare("x1",wInX-36) << " <= "<< opX << range(wInX-1,36) << ";" << endl;
					vhdl << tab << declare("y1",18) << " <= "<< opY << range(17,0) << ";" << endl;
					vhdl << tab << declare("x3",18) << " <= "<< opX << range(17,0) << ";" << endl;
					vhdl << tab << declare("y3",wInY-36) << " <= "<< opY << range(wInY-1,36) << ";" << endl;
					vhdl << tab << declare("y2",18) << " <= "<< opY << range(35,18) << ";" << endl;
					vhdl << tab << declare("x2",18) << " <= "<< opX << range(35,18) << ";" << endl;

					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPMultiplierDelay() + target->DSPToLogicWireDelay());

					vhdl << tab << declare("r11",18 + (wInX-36)) << "<= x1 * y1;" << endl;
					vhdl << tab << declare("r33",18 + (wInY-36)) << "<= x3 * y3;" << endl;
					vhdl << tab << declare("r12",18 + (wInX-36)) << "<= x1 * y2;" << endl;
					vhdl << tab << declare("r23",18 + (wInY-36)) << "<= x2 * y3;" << endl;
					
					manageCriticalPath( target->DSPAdderDelay() );
					int maxOperandWidth = max (18 + (wInX-36), 18 + (wInY-36) );
					vhdl << tab << declare ("sum_r11_r33", maxOperandWidth+1 )  << " <= (" << zg(maxOperandWidth+1 - (18 + (wInX-36)),0) << " & r11) + "
					                                                            << "  (" << zg(maxOperandWidth+1 - (18 + (wInY-36)),0) << " & r33);" << endl;
					vhdl << tab << declare ("sum_r12_r23", maxOperandWidth+1 )  << " <= (" << zg(maxOperandWidth+1 - (18 + (wInX-36)),0) << " & r12) + "
					                                                            << "  (" << zg(maxOperandWidth+1 - (18 + (wInY-36)),0) << " & r23);" << endl;

				///////////////////////////////////////////////////////

				setCycle(0); //reset for the third multiplication
					setCriticalPath( getMaxInputDelays(inputDelays));
					manageCriticalPath( target->DSPMultiplierDelay() + target->DSPToLogicWireDelay());
					vhdl << tab << declare("r13",(wInY-36) + (wInX-36)) << "<= x1 * y3;" << endl;
					
				//syncronization
					syncCycleFromSignal("sum_r11_r33");
					syncCycleFromSignal("r00");	
					
//							manageCriticalPath( target->DSPToLogicWireDelay() );
//							nextCycle();/////////////// DSP out

					vhdl << tab << declare("addOp1",36 + (wInY-36) + (wInX-36)) << " <= r13  & r00" << range(71, 36) << ";" << endl;
					vhdl << tab << declare("addOp2",36 + (wInY-36) + (wInX-36)) << " <= "
					            << zg(36 + (wInY-36) + (wInX-36)  - (maxOperandWidth+1)   ,0) << " & " << "sum_r11_r33;" << endl;
					vhdl << tab << declare("addOp3",36 + (wInY-36) + (wInX-36)) << " <= "
					            << zg(36 + (wInY-36) + (wInX-36)  - (maxOperandWidth+1+18) ,0) << " & " << "sum_r12_r23 & " << zg(18,0) <<";" << endl;
					
					manageCriticalPath( target->adderDelay(36 + (wInY-36) + (wInX-36)) );							
					vhdl << tab << declare("addRes", 36 + (wInY-36) + (wInX-36)) << " <= addOp1 + addOp2 + addOp3;"<<endl;
					
					vhdl << tab << "R <= addRes & r00"<<range(35,0) << ";" << endl;
					outDelayMap["R"]=getCriticalPath();
			} else{
				cout << "ooooups!! 2" << endl;
				exit(-1);
			
			}
			
		
		
		}	
		else if (0) // the target is a Stratix
			{
				bool quadMultiply = false; // TRUE if we cannot use double the chunckSize width multipliers

				int chunkSize_ = 1, x, y;
	
				quadMultiply = !target->suggestSubmultSize(x, y, wInX, wInY);
				chunkSize_ = max(x, y);
		
//						setCriticalPath( getMaxInputDelays(inputDelays)); //FIXME remove hardcoded register insertions
//						manageCriticalPath( target->DSPMultiplierDelay() );
				delay += target->distantWireDelay((int)target->frequencyMHz());
				
				
				
		
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
		
						vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= Y & " << zg(chunkSize_*chunksX-widthX,0) << ";" << endl;
						vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= X & " << zg(chunkSize_*chunksY-widthY,0) << ";" << endl;
					}
				else
					{
						vhdl << tab << declare("sX",chunkSize_*chunksX) << " <= X & " << zg(chunkSize_*chunksX-widthX,0) << ";" << endl;
						vhdl << tab << declare("sY",chunkSize_*chunksY) << " <= Y & " << zg(chunkSize_*chunksY-widthY,0) << ";" << endl;
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
								vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
							}
						for (int k=0; k<chunksY ; k++)
							{
								ostringstream dname;
								dname << "y"<<k;
								vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range((k+1)*chunkSize_-1,k*chunkSize_) << ";" << endl;
							}
		
						if (verbose)	
							cerr << "> IntMultiplier: Cannot use double chunckSize_ multipliers" << endl;
			
						while (chX/4 > 0)
							{
								// top-right tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px0y" << 4*level;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y", 4*level)) << ";" << endl;
			
								setCycle(3);
								operands[0] << use(partialProd.str()) << " & " << zg(chunkSize_*4*level,0) << ";" << endl;
			
			
								// top-right diagonal pair
								setCycle(1);
								partialProd2.str("");
								partialProd2 << "px0y" << 4*level+1;
								vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",4*level+1)) << ";" << endl;
								partialProd.str("");
								partialProd << "px1y" << 4*level;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y",4*level)) << ";" << endl;
			
								setCycle(2);
								sum.str("");
								sum << "addOp" << level << "_shift_" << 4*level+1;
								vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
			
								setCycle(3);
								operands[1] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << zg((4*level+1)*chunkSize_,0) << ";" << endl;
								carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zg((4*level+3)*chunkSize_,0) << ";" << endl;
			
			
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
								carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-1,0) << " & " << carrys.str();
			
				
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
										carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
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
										carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
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
								carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
			
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
								operands[(chX+chunksY-3)%2] << zg((4*level+1)*chunkSize_,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(chX+chunksY-3)%2].str();
								carrys.seekp(ios_base::beg);
								carrys << zg((4*level+1)*chunkSize_-1,0) << " & " << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
			
								// bottom-left tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px"<< chX-1 << "y" << chunksY-1;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-1)) << " * " << use(join("y", chunksY-1)) << ";" << endl;
			
								setCycle(3);
								operands[(chX+chunksY-2)%2].seekp(ios_base::beg);
								operands[(chX+chunksY-2)%2] << zg(4*level*chunkSize_,0) << " & " << use(partialProd.str()) << " & " << operands[(chX+chunksY-2)%2].str();
			
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
				
								operands[i%2] << zg(level*4*chunkSize_,0) << ";";
								if (chunksY-i > 1) // then there are more than one tiles in the column
									{
										if (verbose)
											cerr << ">IntMultiplier: More than one tile in the column" << endl;
										operands[(i+1)%2] << zg((level*4+1)*chunkSize_,0) << ";";
										oneCol = false;
									}
				
								for (i=level*4; i<chunksY; i++)
									{
										setCycle(1);
										partialProd.str("");
										partialProd << "px0y" << i;
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i)) << ";" << endl;
					
										setCycle(3);
										operands[i%2].seekp(ios_base::beg);
										operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
					
									}
				
								if (!oneCol) // then there are more than one tiles in the column
									{
										operands[i%2].seekp(ios_base::beg);
										operands[i%2] << zg((level*4+1)*chunkSize_,0) << " & " << operands[i%2].str();
										vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[i%2].str() << endl;
										opCount++;
									}
				
								operands[(i+1)%2].seekp(ios_base::beg);
								operands[(i+1)%2] << zg(level*4*chunkSize_,0) << " & " << operands[(i+1)%2].str(); 
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
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y", i)) << ";" << endl;
								setCycle(3);
								operands[i%2] << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << zg((4*level)*chunkSize_,0) << ";";
								operands[(i+1)%2] << zg((4*level+1)*chunkSize_,0) << ";";
								carrys << zg((4*level+2)*chunkSize_+1,0) << ";";
				 
								while (i<chunksY-1) // diagonal pairs
									{
										setCycle(1);
										partialProd2.str("");
										partialProd2 << "px0y" << (i+1);
										vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y", i+1)) << ";" << endl;
										partialProd.str("");
										partialProd << "px1y" << i;
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y", i)) << ";" << endl;
										setCycle(2);
										sum.str("");
										sum << "addOp" << level << "_shift_" << (i+1);
										vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
										setCycle(3);
										operands[(i+1)%2].seekp(ios_base::beg);
										operands[(i+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
										carrys.seekp(ios_base::beg);
										carrys << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zg(chunkSize_-1,0) << " & " << carrys.str();
			
										i++;			
									}
				
								// bottom-left tile
								setCycle(1);
								partialProd.str("");
								partialProd << "px1y" << i;
								vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y", i)) << ";" << endl;
				
								setCycle(3);
								operands[(i+1)%2].seekp(ios_base::beg);
								operands[(i+1)%2] << zg((4*level)*chunkSize_,0) << " & " << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
				
								operands[i%2].seekp(ios_base::beg);
								operands[i%2] << zg((4*level+1)*chunkSize_,0) << " & " << operands[i%2].str();
								carrys.seekp(ios_base::beg);
								carrys << zg((4*level+1)*chunkSize_-1,0) << " & " << carrys.str();
				
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
					
										operands[i%2] << zg(level*4*chunkSize_,0) << ";";
										operands[(i+1)%2] << zg((level*4+1)*chunkSize_,0) << ";";
				
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
										operands[(3+i)%2] << zg((level*4+1)*chunkSize_,0) << " & " << operands[(3+i)%2].str();
										vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(3+i)%2].str() << endl;
										opCount++;
						
										operands[(i+4)%2].seekp(ios_base::beg);
										operands[(i+4)%2] << zg(level*4*chunkSize_,0) << " & " << operands[(i+4)%2].str(); 
										vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(i+4)%2].str() << endl;
										opCount++;
										break;
									case 2: // two rows
										if (verbose)
											cerr << ">IntMultiplier:" << tab << " Subcase 2: there are two rows" << endl;
						
										operands[i%2] << zg(level*4*chunkSize_,0) << ";";
										operands[(i+1)%2] << zg((level*4+1)*chunkSize_,0) << ";";
										carrys << zg(level*4*chunkSize_,0) << ";";
										// top-right tile
										setCycle(1);
										partialProd.str("");
										partialProd << "px0y" << i;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i)) << ";" << endl;
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
												carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zg(chunkSize_,0) << "&" << carrys.str();
											}
						
										// bottom-left tile
										setCycle(1);
										partialProd.str("");
										partialProd << "px2y" << i+2;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x2 * " << use(join("y",i+2)) << ";" << endl;
										setCycle(3);
										operands[(i+4)%2].seekp(ios_base::beg);
										operands[(i+4)%2] << zg(level*4*chunkSize_,0) << " & " << use(partialProd.str()) << " & " << operands[(i+4)%2].str();	
										operands[(i+5)%2].seekp(ios_base::beg);
										operands[(i+5)%2] << zg((level*4+1)*chunkSize_,0) << " & " << operands[(i+5)%2].str(); 
										carrys.seekp(ios_base::beg);
										carrys << zg((level*4+1)*chunkSize_-1,0) << " & " << carrys.str();
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
										//delay += target->distantWireDelay((int)target->frequencyMHz());
						
										operands[i%2] << zg(level*4*chunkSize_,0) << ";";
										operands[(i+1)%2] << zg((level*4+1)*chunkSize_,0) << ";";
										carrys << zg((level*4+3)*chunkSize_,0) << ";";
						
										// top-right tile
										setCycle(1);
										partialProd.str("");
										partialProd << "px0y" << i;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i)) << ";" << endl;
										setCycle(3);
										operands[i%2].seekp(ios_base::beg);
										operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
						
										// top-right diagonal pair
										setCycle(1);
										partialProd.str("");
										partialProd << "px1y" << i;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y",i)) << ";" << endl;
										partialProd2.str("");
										partialProd2 << "px0y" << i+1;
										vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i+1)) << ";" << endl;
						
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
												carrys << use(sum.str()) << range(2*chunkSize_+1,2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
					
												i++;
											}
						
										// bottom-left diagonal pair
										setCycle(1);
										partialProd.str("");
										partialProd << "px2y" << chunksY-2;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x2 * " << use(join("y",chunksY-2)) << ";" << endl;
										partialProd2.str("");
										partialProd2 << "px1y" << chunksY-1;
										vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y",chunksY-1)) << ";" << endl;
						
										setCycle(2);
										sum.str("");
										sum << "addOpd" << level << "_shift_" << chunksY;
										vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
						
										setCycle(3);
										operands[(i+2)%2].seekp(ios_base::beg);
										operands[(i+2)%2] << zg((level*4+1)*chunkSize_,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+2)%2].str();
										carrys.seekp(ios_base::beg);
										carrys << zg((level*4+1)*chunkSize_-1,0) << " & " << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str(); 
										// bottom-left tile
										setCycle(1);
										partialProd.str("");
										partialProd << "px2y" << chunksY-1;
										vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x2 * " << use(join("y",chunksY-1)) << ";" << endl;
										setCycle(3);
										operands[(i+3)%2].seekp(ios_base::beg);
										operands[(i+3)%2] << zg((level*4)*chunkSize_,0) << " & " << use(partialProd.str()) << " & " << operands[(i+3)%2].str();
						
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
				//**********************************************************************************************************************
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
										vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range((k+2)*chunkSize_-1,k*chunkSize_) << ";" << endl;
									}
			
								// * The following splitting for y uses a small trick, i.e. if we have an even number of chunks for y
								// * nothing changes and the splitting procedes normally. Otherwise we skip over the first chunk y(17:0)
								// * and start counting the chunks from then on. At the end the last indexed chunk will be y(17:0).
								// * This was done to avoid major changes in the partial product computation and concatenation algorithm.
								int start = 0;
			
								if (chunksY%2 == 1)
									start = 1;
				
								for (k=0; k<chunksY-1 ; k+=2)
									{
										ostringstream dname;
										dname << "y"<< k;
										vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range((k+2+start)*chunkSize_-1,(k+start)*chunkSize_) << ";" << endl;
									}
			
								if (start == 1)
									{
										ostringstream dname;
										dname << "y"<< chunksY-1;
										vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range(chunkSize_-1,0) << ";" << endl;
									}
			
								nextCycle();
			
								// COMPUTE PARTIAL PRODUCTS
								for (i=0; i<chunksX-1; i+=2)
									{
										operands[0] << zg((i+start)*chunkSize_,0) << ";";
										operands[1] << zg((i+2+start)*chunkSize_,0) << ";";
				
										// compute and concatenate partial products from right to left
										for (k=0; k<chunksX-i-1; k+=2)
											{
												setCycleFromSignal(join("x",k));
												nextCycle();
												ostringstream partialProd;
												partialProd  << "px" << k << "y" << i;
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
												partialProd  << "px" << k-2 << "y" << j;
												vhdl << tab << declare(partialProd.str(),4*chunkSize_) << " <= " << use(join("x",k-2)) << " * " << use(join("y",j)) << ";" << endl;
				
												nextCycle();
												operands[((j-i+k-2)/2)%2].seekp(ios_base::beg);
												operands[((j-i+k-2)/2)%2] << use(partialProd.str()) << " & " << operands[((j-i+k-2)/2)%2].str();
											}
				
										operands[((i+j+k-2)/2)%2].seekp(ios_base::beg);
										operands[((i+j+k-2)/2)%2] << zg((i+2+chunksX%2)*chunkSize_,0) << " & " << operands[((i+j+k-2)/2)%2].str();	
										operands[((i+j+k-2)/2+1)%2].seekp(ios_base::beg);
										operands[((i+j+k-2)/2+1)%2] << zg((i+chunksX%2)*chunkSize_,0) << " & " << operands[((i+j+k-2)/2+1)%2].str();	
				
										// reinitialize the concatenation buffers for the next iteration of the loop
										for (k=0; k<2; k++)
											{
												vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[k].str() << endl;
												opCount++;
												operands[k].str("");
											}
									}
			
								operands[0] << "\"\";";
								operands[1] << zg(chunkSize_,0) << ";";
			
								bool halfPadded = false; // TRUE when there is either 
								//	a chunkSize width column on the left side of the tiling or 
								//	a chunkSize width row on the bottom side of the tiling 	
			
								if (start == 1) // then there is a chunkSize width row on the bottom side of the tiling
									{
										halfPadded = true;
				
										// compute and concatenate partial products from right to left
										for (k=0; k<chunksX; k++)
											{
												setCycleFromSignal("sX");
					
												ostringstream dname;
												dname << "xx" << k;
												vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
					
												nextCycle();
					
												ostringstream partialProd;
												partialProd  << "pxx" << k << "y" << chunksY-1;
												vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("xx",k)) << " * " << use(join("y",chunksY-1)) << ";" << endl;
					
												nextCycle();
					
												operands[k%2].seekp(ios_base::beg);
												operands[k%2] << use(partialProd.str()) << " & " << operands[k%2].str();
											}
									}
								else // padd operands with zeros
									{
										operands[0].seekp(ios_base::beg);
										operands[0] << zg((chunksX-1)*2*chunkSize_,0) << " & " << operands[0].str();
										operands[1].seekp(ios_base::beg);
										operands[1] << zg((chunksX-1)*2*chunkSize_,0) << " & " << operands[1].str();
									}
			
								if (i == chunksX-1) // then there is a chunkSize width column on the left side of the tiling
									{
										halfPadded = true;
										setCycleFromSignal("sX");
				
										ostringstream dname;
										dname << "x" << chunksX-1;
										vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range(chunksX*chunkSize_-1,(chunksX-1)*chunkSize_) << ";" << endl;
				
										for (k=start; k<chunksY; k++)
											{
												setCycleFromSignal("sY");
												dname.str("");
												dname << "yy" << k;
												vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
				
												nextCycle();
												ostringstream partialProd;
												partialProd  << "px" << chunksX-1 << "yy" << k;
												vhdl << tab << declare(partialProd.str(),2*chunkSize_) << " <= " << use(join("x",chunksX-1)) << " * " << use(join("yy",k)) << ";" << endl;
				
												nextCycle();
												operands[(k+chunksX+1)%2].seekp(ios_base::beg);
												operands[(k+chunksX+1)%2] << use(partialProd.str()) << " & " << operands[(k+chunksX+1)%2].str();
											}
					
									}
								else // padd operands with zeros
									{
										operands[0].seekp(ios_base::beg);
										operands[0] << zg((chunksY-2)*2*chunkSize_,0) << " & " << operands[0].str();
										operands[1].seekp(ios_base::beg);
										operands[1] << zg((chunksY-2)*2*chunkSize_,0) << " & " << operands[1].str();	
									}
			
								if (halfPadded) // then we add the two operands to the sum
									{
										if ((chunksX+chunksY-1)%2 == 0)
											{
												operands[0].seekp(ios_base::beg); 
												operands[0] << zg(chunkSize_,0) << " & " << operands[0].str(); 	
											}
										else
											{
												operands[1].seekp(ios_base::beg); 
												operands[1] << zg(chunkSize_,0) << " & " << operands[1].str(); 	
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
			
								int i, k;	
			
								////////////////////////////////////////////////////
								//SPLITTINGS
								for (k=0; k<chunksX-1; k+=2)
									{
										ostringstream dname;
										dname << "x"<< k/2;
										vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range((k+2)*chunkSize_-1,k*chunkSize_) << ";" << endl;
									}
			
								//* The following splitting for y uses a small trick, i.e. if we have an even number of chunks for y
								//* nothing changes and the splitting procedes normally. Otherwise we skip over the first chunk y(17:0)
								//* and start counting the chunks from then on. At the end the last indexed chunk will be y(17:0).
								//* This was done to avoid major changes in the partial product computation and concatenation algorithm.
								int start = 0;
			
								if (chunksY%2 == 1)
									start = 1;
				
								for (k=0; k<chunksY-1 ; k+=2)
									{
										ostringstream dname;
										dname << "y"<< k/2;
										vhdl << tab << declare(dname.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range((k+2+start)*chunkSize_-1,(k+start)*chunkSize_) << ";" << endl;
									}
			
								if (start == 1)
									{
										ostringstream dname;
										dname << "y"<< k/2;
										vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range(chunkSize_-1,0) << ";" << endl;
									}
								int startX = 0; // added this variable to apply the same trick for zgs as with the start variable that adds zeros behind if chunksY is odd 
								bool leftColumn = false; // TRUE if there are an odd number of slices in X and there is a chunkSize width column to the left of the tiling that has no pair
								if (chunksX % 2 == 1)
									{
										leftColumn = true;
										startX = 1;
									}
			
								//* modifing the quantities so that we can use the same algorithm as for the case where we have no
								//* double chunkSize multipliers 
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
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y", 4*level)) << ";" << endl;
										setCycle(3);
										operands[0] << use(partialProd.str()) << " & " << zg(chunkSize_*4*level + start*chunkSize_/2,0) << ";" << endl;
				
				
										// top-right diagonal pair
										setCycle(1);
										partialProd2.str("");
										partialProd2 << "px0y" << 4*level+1;
										vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",4*level+1)) << ";" << endl;
										partialProd.str("");
										partialProd << "px1y" << 4*level;
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y",4*level)) << ";" << endl;
										setCycle(2);
										sum.str("");
										sum << "addOp" << level << "_shift_" << 4*level+1;
										vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
										setCycle(3);
										operands[1] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << zg((4*level+1)*chunkSize_ + start*chunkSize_/2,0) << ";" << endl;
										carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zg((4*level+3)*chunkSize_ + start*chunkSize_/2,0) << ";" << endl;
				
				
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
										carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-1,0) << " & " << carrys.str();
				
					
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
												carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
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
												carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
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
										carrys << use(sum.str()) << range(2*chunkSize_+1, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
				
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
										operands[(chX+chunksY-3)%2] << zg((4*level+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(chX+chunksY-3)%2].str();
										carrys.seekp(ios_base::beg);
										carrys << zg((4*level+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
				
										// bottom-left tile
										setCycle(1);
										partialProd.str("");
										partialProd << "px"<< chX-1 << "y" << chunksY-1;
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= " << use(join("x",chX-1)) << " * " << use(join("y", chunksY-1)) << ";" << endl;
										setCycle(3);
										operands[(chX+chunksY-2)%2].seekp(ios_base::beg);
										operands[(chX+chunksY-2)%2] << zg(4*level*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << " & " << operands[(chX+chunksY-2)%2].str();
				
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
					
										operands[i%2] << zg(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
										if (chunksY-i > 1) // then there are more than one tiles in the column
											{
												if (verbose)
													cerr << ">IntMultiplier: More than one tile in the column" << endl;
												operands[(i+1)%2] << zg((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
												oneCol = false;
											}
					
										for (i=level*4; i<chunksY; i++)
											{
												setCycle(1);
												partialProd.str("");
												partialProd << "px0y" << i;
												vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i)) << ";" << endl;
						
												setCycle(3);
												operands[i%2].seekp(ios_base::beg);
												operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
						
											}
					
										if (!oneCol) // then there are more than one tiles in the column
											{
												operands[i%2].seekp(ios_base::beg);
												operands[i%2] << zg((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[i%2].str();
												vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[i%2].str() << endl;
												opCount++;
											}
					
										operands[(i+1)%2].seekp(ios_base::beg);
										operands[(i+1)%2] << zg(level*4*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(i+1)%2].str(); 
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
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y", i)) << ";" << endl;
										setCycle(3);
										operands[i%2] << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << zg((4*level)*chunkSize_ + start*chunkSize_/2,0) << ";";
										operands[(i+1)%2] << zg((4*level+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
										carrys << zg((4*level+2)*chunkSize_+1 + start*chunkSize_/2,0) << ";";
					 
										while (i<chunksY-1) // diagonal pairs
											{
												setCycle(1);
												partialProd2.str("");
												partialProd2 << "px0y" << (i+1);
												vhdl << tab << declare(partialProd2.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y", i+1)) << ";" << endl;
												partialProd.str("");
												partialProd << "px1y" << i;
												vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y", i)) << ";" << endl;
												setCycle(2);
												sum.str("");
												sum << "addOp" << level << "_shift_" << (i+1);
												vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
												setCycle(3);
												operands[(i+1)%2].seekp(ios_base::beg);
												operands[(i+1)%2] << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
												carrys.seekp(ios_base::beg);
												carrys << use(sum.str()) << range(2*chunkSize_, 2*chunkSize_) << " & " << zg(chunkSize_-1,0) << " & " << carrys.str();
				
												i++;			
											}
					
										// bottom-left tile
										setCycle(1);
										partialProd.str("");
										partialProd << "px1y" << i;
										vhdl << tab << declare(partialProd.str(),2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y", i)) << ";" << endl;
					
										setCycle(3);
										operands[(i+1)%2].seekp(ios_base::beg);
										operands[(i+1)%2] << zg((4*level)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+1)%2].str();
					
										operands[i%2].seekp(ios_base::beg);
										operands[i%2] << zg((4*level+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[i%2].str();
										carrys.seekp(ios_base::beg);
										carrys << zg((4*level+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << carrys.str();
				
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
						
												operands[i%2] << zg(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
												operands[(i+1)%2] << zg((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
					
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
												operands[(3+i)%2] << zg((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(3+i)%2].str();
												vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(3+i)%2].str() << endl;
												opCount++;
							
												operands[(i+4)%2].seekp(ios_base::beg);
												operands[(i+4)%2] << zg(level*4*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(i+4)%2].str(); 
												vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[(i+4)%2].str() << endl;
												opCount++;
												break;
											case 2: // two rows
												if (verbose)
													cerr << ">IntMultiplier: " << tab << "Subcase 2: there are two rows" << endl;
							
												operands[i%2] << zg(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
												operands[(i+1)%2] << zg((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
												carrys << zg(level*4*chunkSize_,0) << ";";
												// top-right tile
												setCycle(1);
												partialProd.str("");
												partialProd << "px0y" << i;
												vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i)) << ";" << endl;
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
														carrys << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zg(chunkSize_,0) << "&" << carrys.str();
													}
							
												// bottom-left tile
												setCycle(1);
												partialProd.str("");
												partialProd << "px2y" << i+2;
												vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x2 * " << use(join("y",i+2)) << ";" << endl;
												setCycle(3);
												operands[(i+4)%2].seekp(ios_base::beg);
												operands[(i+4)%2] << zg(level*4*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << " & " << operands[(i+4)%2].str();	
												operands[(i+5)%2].seekp(ios_base::beg);
												operands[(i+5)%2] << zg((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << operands[(i+5)%2].str(); 
												carrys.seekp(ios_base::beg);
												carrys << zg((level*4+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << carrys.str();
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
							
												operands[i%2] << zg(level*4*chunkSize_ + start*chunkSize_/2,0) << ";";
												operands[(i+1)%2] << zg((level*4+1)*chunkSize_ + start*chunkSize_/2,0) << ";";
												carrys << zg((level*4+3)*chunkSize_ + start*chunkSize_/2,0) << ";";
							
												// top-right tile
												setCycle(1);
												partialProd.str("");
												partialProd << "px0y" << i;
												vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i)) << ";" << endl;
												setCycle(3);
												operands[i%2].seekp(ios_base::beg);
												operands[i%2] << use(partialProd.str()) << " & " << operands[i%2].str();
							
												// top-right diagonal pair
												setCycle(1);
												partialProd.str("");
												partialProd << "px1y" << i;
												vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y",i)) << ";" << endl;
												partialProd2.str("");
												partialProd2 << "px0y" << i+1;
												vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x0 * " << use(join("y",i+1)) << ";" << endl;
							
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
														carrys << use(sum.str()) << range(2*chunkSize_+1,2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str();
						
														i++;
													}
							
												// bottom-left diagonal pair
												setCycle(1);
												partialProd.str("");
												partialProd << "px2y" << chunksY-2;
												vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x2 * " << use(join("y",chunksY-2)) << ";" << endl;
												partialProd2.str("");
												partialProd2 << "px1y" << chunksY-1;
												vhdl << tab << declare(partialProd2.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x1 * " << use(join("y",chunksY-1)) << ";" << endl;
							
												setCycle(2);
												sum.str("");
												sum << "addOpd" << level << "_shift_" << chunksY;
												vhdl << tab << declare(sum.str(), 2*chunkSize_+1,true,Signal::registeredWithAsyncReset) << " <= (\"0\" & " << use(partialProd.str()) << ") + (\"0\" & " << use(partialProd2.str()) << ");" << endl; 
							
												setCycle(3);
												operands[(i+2)%2].seekp(ios_base::beg);
												operands[(i+2)%2] << zg((level*4+1)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_-1,0) << " & " << operands[(i+2)%2].str();
												carrys.seekp(ios_base::beg);
												carrys << zg((level*4+1)*chunkSize_-1 + startX*chunkSize_/2,0) << " & " << use(sum.str()) << range(2*chunkSize_,2*chunkSize_) << " & " << zg(chunkSize_-2,0) << " & " << carrys.str(); 
												// bottom-left tile
												setCycle(1);
												partialProd.str("");
												partialProd << "px2y" << chunksY-1;
												vhdl << tab << declare(partialProd.str(), 2*chunkSize_,true,Signal::registeredWithAsyncReset) << " <= x2 * " << use(join("y",chunksY-1)) << ";" << endl;
												setCycle(3);
												operands[(i+3)%2].seekp(ios_base::beg);
												operands[(i+3)%2] << zg((level*4)*chunkSize_ + startX*chunkSize_/2,0) << " & " << use(partialProd.str()) << " & " << operands[(i+3)%2].str();
							
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
								operands[1] << zg(chunkSize_,0) << ";";
			
								bool halfPadded = false; // TRUE when there is either 
								//	a chunkSize width column on the left side of the tiling or 
								//	a chunkSize width row on the bottom side of the tiling 	
			
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
												vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
					
												//nextCycle();
												setCycle(1);
												ostringstream partialProd;
												partialProd  << "pxx" << k << "y" << chunksY/2;
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
										operands[0] << zg((chunksX-1)*chunkSize_,0) << " & " << operands[0].str();
										operands[1].seekp(ios_base::beg);
										operands[1] << zg((chunksX-1)*chunkSize_,0) << " & " << operands[1].str();
									}
			
								if (leftColumn) // then there is a chunkSize width column on the left side of the tiling
									{
										halfPadded = true;
										//setCycleFromSignal("sX");
										setCycle(0);
										ostringstream dname;
										dname << "x" << chunksX-1;
										vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sX" << range(chunksX*chunkSize_-1,(chunksX-1)*chunkSize_) << ";" << endl;
				
										for (k=start; k<chunksY; k++)
											{
												//setCycleFromSignal("sY");
												setCycle(0);
												dname.str("");
												dname << "yy" << k;
												vhdl << tab << declare(dname.str(),chunkSize_,true,Signal::registeredWithAsyncReset) << " <= sY" << range((k+1)*chunkSize_-1, k*chunkSize_) << ";" << endl;
				
												//nextCycle();
												setCycle(1);
												ostringstream partialProd;
												partialProd  << "px" << chunksX-1 << "yy" << k;
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
										operands[0] << zg((chunksY-1)*chunkSize_,0) << " & " << operands[0].str();
										operands[1].seekp(ios_base::beg);
										operands[1] << zg((chunksY-1)*chunkSize_,0) << " & " << operands[1].str();	
									}
			
								//setCycle(3);
			
								if (halfPadded) // then we add the two operands to the sum
									{
										if ((chunksX+chunksY-1)%2 == 0)
											{
												operands[0].seekp(ios_base::beg); 
												operands[0] << zg(chunkSize_,0) << " & " << operands[0].str(); 	
											}
										else
											{
												operands[1].seekp(ios_base::beg); 
												operands[1] << zg(chunkSize_,0) << " & " << operands[1].str(); 	
											}
				
										for (k=0; k<2; k++)
											{
												vhdl << tab << declare(join("addOp", opCount), adderWidth) << " <= " << operands[k].str() << endl;
												opCount++;
											}
									}
							}
					}	
	
				map<string, double> inMap;
				inMap["X0"] = delay;

				IntNAdder* add =  new IntNAdder(target, adderWidth, opCount, inMap);
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
	
				vhdl << tab << "R <= addRes" << range(adderWidth-1,adderWidth-wInX_-wInY_) << ";" << endl;	
	
			}
		}
		/***************************************************************************
		 ***************************************************************************/
		else{ //IntMultiplier Version -- multiplication performed in LUTs. Works for both Xilinx and Altera
			int chunkSize_ = target->lutInputs()/2;
			int chunksX =  int(ceil( ( double(wInX) / (double) chunkSize_) ));
			int chunksY =  int(ceil( ( double(wInY) / (double) chunkSize_) ));
		
			setCriticalPath( getMaxInputDelays(inputDelays) );

			REPORT(DEBUG, "X splitted in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
			if (chunksX + chunksY > 2) { //we do more than 1 subproduct
				int widthX = wInX_;
				int widthY = wInY_;	
				bool swap = false;

				if (chunksX > chunksY){ //interchange X with Y
					int tmp = chunksX;
					chunksX = chunksY;
					chunksY = tmp;
		
					tmp = widthX;
					widthX = widthY;
					widthY = tmp;
					swap = true;
				}		

				vhdl<<tab<<declare("sX",chunkSize_*chunksX)<<" <= "<< (swap? "Y":"X")<<" & "<<zg( chunkSize_ * chunksX -widthX,0)<< ";"<<endl;
				vhdl<<tab<<declare("sY",chunkSize_*chunksY)<<" <= "<< (swap? "X":"Y")<<" & "<<zg(chunkSize_*chunksY-widthY,0)<< ";"<<endl;	
				//SPLITTINGS
				for (int k=0; k<chunksX ; k++)
					vhdl<<tab<<declare(join("x",k),chunkSize_)<<" <= sX"<<range((k+1)*chunkSize_-1,k*chunkSize_)<<";"<<endl;
				for (int k=0; k<chunksY ; k++)
					vhdl<<tab<<declare(join("y",k),chunkSize_)<<" <= sY"<<range((k+1)*chunkSize_-1,k*chunkSize_)<<";"<<endl;

				manageCriticalPath( target->lutDelay() + target->lutInputs()*target->localWireDelay() );
			
				if (sign_){
					/* ----------------------------------------------------*/
					//COMPUTE PARTIAL PRODUCTS
					for (int i=0; i<chunksY; i++)
						for (int j=0; j<chunksX; j++){
							if ((i<chunksY-1) && (j<chunksX-1)){
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_+2)<<" <= ( \"0\" & " << join("x",j) <<") * " 
		                                            << "( \"0\" & " << join("y",i)<< ");" << endl;
							} else if ((i==chunksY-1) && (j<chunksX-1)){
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_+1)<<" <= ( \"0\" & " << join("x",j) <<") * " 
		                                            << "(" << join("y",i)<< ");" << endl;
							} else if ((i<chunksY-1) && (j==chunksX-1)){	
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_+1)<<" <= ( " << join("x",j) <<") * " 
		                                            << "( \"0\" & " << join("y",i)<< ");" << endl;
							} else{
								vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_)<<" <= "<< join("x",j)<<" * " << join("y",i)<< ";" << endl;
							}								
						}		

					map<string, int> availableZeros;
					map<string, int> neededSignExtensions;

					int adderWidth = (chunksX+chunksY)*chunkSize_;
					// CONCATENATE PARTIAL PRODUCTS
					for (int i=0; i<2; i++){
						for (int j=0; j<chunksX; j++){
							int startIdx = chunksY-1-i;
							int paddWidth = adderWidth - chunkSize_*2*(startIdx/2+1);
							int endPaddWidth = chunkSize_*(j+startIdx%2);

							if (j!=chunksX-1){ //the last chunk makes us problems 
								if (i!=0)
									availableZeros[join("cp",i,j)] = paddWidth-endPaddWidth; 
							}else{
								//each of these subproducts requires sign extension
								for (int k=startIdx; k>=0; k-=2)
									if (k!=startIdx)
										neededSignExtensions[join("px",j,"y",k)]= adderWidth - (j+k+2)*chunkSize_;
							}
						}
					}

					map<string,int>::iterator it, it2;

					REPORT(DEBUG, "available ");						
					for (it = availableZeros.begin(); it!=availableZeros.end(); it++)
						REPORT(DEBUG, it->first << " -> " << it->second);

					REPORT(DEBUG, "needed");
					for (it = neededSignExtensions.begin(); it!=neededSignExtensions.end(); it++)
						REPORT(DEBUG, it->first << " -> " << it->second);
				
					map< string, pair<string,int> > reallocations;
					map< string, int> unallocated; 
				
					for (it=neededSignExtensions.begin(); it!= neededSignExtensions.end(); it++){
						bool reallocationsStatus = false;
						REPORT(DEBUG, "trying to allocate " << it->first << "," << it->second);
						for (it2=availableZeros.begin(); it2!= availableZeros.end() && (!reallocationsStatus); it2++){
							if (it2->second >= it->second){ //we have enough zeros to allocate
								reallocations[ it2->first ] = pair<string,int>(it->first, it->second); 
								reallocationsStatus = true;
								availableZeros.erase(it2->first); //delete the entry.
							}
						}
						if (!reallocationsStatus) // we could not realocate
							unallocated[it->first] = it->second;
					}
				
					map< string, pair<string,int> >::iterator it3;
				
					REPORT(DEBUG, "reallocated ");
					for (it3 = reallocations.begin(); it3!=reallocations.end(); it3++)
						REPORT(DEBUG, it3->first << " --> " << it3->second.first << ", " << it3->second.second);

					REPORT(DEBUG, "unallocated");
					for (it = unallocated.begin(); it!=unallocated.end(); it++)
						REPORT(DEBUG, it->first << " -> " << it->second);
				
					vector<string> operandNames;
				
					for (int i=0; i<2; i++){
						for (int j=0; j<chunksX; j++){
							int startIdx = chunksY-1-i;
							int paddWidth = adderWidth - chunkSize_*2*(startIdx/2+1);
							int endPaddWidth = chunkSize_*(j+startIdx%2);

							operandNames.push_back(join("cp",i,j));
							vhdl << tab << declare(join("cp",i,j),adderWidth) << " <= ";

							if (j!=chunksX-1){ //the last chunk makes us problems 
								if (i==0){ //multiply with the last chunk of y, so sign extend
									vhdl << rangeAssign(paddWidth-endPaddWidth-1, 0, join("px",j,"y",startIdx)+of(2*chunkSize_-1));
								}else{
									//we check the reallocations. Some zeros might be reallocated
									it3 = reallocations.find(join("cp",i,j));
	
									if (it3!=reallocations.end()){ // some reallocation happended here
										REPORT(DEBUG, "reallocated " << it3->second.second << " sign bits over " << paddWidth-endPaddWidth << " zero bits");
										vhdl << rangeAssign( it3->second.second-1, 0, it3->second.first+of(2*chunkSize_-1)) << " & " << zg(paddWidth-endPaddWidth-it3->second.second,0) ; 	
									}else{
										vhdl << zg(paddWidth-endPaddWidth,0) ; 	
									}
								}
								for (int k=startIdx; k>=0; k-=2)
									vhdl << " & " << use(join("px",j,"y",k))<<range(2*chunkSize_-1,0);
								vhdl << " & " << zg((startIdx<0?2:endPaddWidth),0) << ";" << endl;
							}else{
								//the ms one does't need any reallocations
								vhdl << rangeAssign(paddWidth-endPaddWidth-1, 0, join("px",j,"y",startIdx)+of(2*chunkSize_-1));
								for (int k=startIdx; k>=0; k-=2){
									if (k==startIdx){
										vhdl << " & " << join("px",j,"y",k)<<range(2*chunkSize_-1,0);
									}else{
										//check if the sign extension was reallocated
										bool reallocationsStatus = false;
										for(it3 = reallocations.begin(); it3!=reallocations.end() && (!reallocationsStatus); it3++){
											if (it3->second.first == join("px",j,"y",k)){
												reallocationsStatus = true;
											}
										}
									
										if (reallocationsStatus){ //if the sign extension was reallocated we simply continue
											vhdl << " & " << join("px",j,"y",k)<<range(2*chunkSize_-1,0);
										}else{
											//not reallocated.
											//we need to finish the current cp and start a new one
											vhdl << " & " << zg( (j+1 + k+1)*chunkSize_, 0) << ";" << endl;
											operandNames.push_back(join("cp",i,j,k));
											vhdl << tab << declare(join("cp",i,j,k),adderWidth) << " <= ";
											vhdl << rangeAssign(adderWidth - (j+k+2)*chunkSize_ -1 , 0, join("px",j,"y",k)+of(2*chunkSize_-1))
												 << " & " << join("px",j,"y",k) << range(2*chunkSize_-1,0);
										}											
									}
									if (k-2 <0) //we need to finish
										vhdl << " & " << zg( (j + k)*chunkSize_, 0) << ";--" << endl;
								}
							}
						}
					}

					IntCompressorTree* add =  new IntCompressorTree(target, adderWidth, operandNames.size(), inDelayMap("X0", getCriticalPath()));
					oplist.push_back(add);

					for (int i=0; i< signed(operandNames.size()); i++)
						inPortMap (add, join("X",i) , operandNames[i]);
	
					outPortMap(add, "R", "addRes");
					vhdl << instance(add, "adder");
					syncCycleFromSignal("addRes");
					setCriticalPath(add->getOutputDelay("R"));

					vhdl<<tab<<"R<=addRes" << range(adderWidth-1,adderWidth-wInX_-wInY_) << ";" << endl;		
				}else{ //unsigned version
				/* ----------------------------------------------------*/
				//COMPUTE PARTIAL PRODUCTS
					for (int i=0; i<chunksY; i++)
						for (int j=0; j<chunksX; j++)
							vhdl<<tab<<declare(join("px",j,"y",i),2*chunkSize_)<<" <= "<< join("x",j)<< " * " << join("y",i) << ";" << endl;
					
					int adderWidth = (chunksX+chunksY)*chunkSize_;
					// CONCATENATE PARTIAL PRODUCTS
					for (int i=0; i<2; i++){
						for (int j=0; j<chunksX; j++){
							vhdl << tab << declare(join("cp",i,j),adderWidth) << " <= ";
	
							int startIdx = chunksY-1-i;
							int paddWidth = adderWidth - chunkSize_*2*(startIdx/2+1);
							int endPaddWidth = chunkSize_*(j+startIdx%2);
							vhdl << zg(paddWidth-endPaddWidth,0); 
	
							for (int k=startIdx; k>=0; k-=2)
								vhdl << " & " << use(join("px",j,"y",k));
				
							vhdl << " & " << zg((startIdx<0?2:endPaddWidth),0) << ";" << endl;
						}
					}

					IntCompressorTree* add =  new IntCompressorTree(target, adderWidth, 2*chunksX, inDelayMap("X0", getCriticalPath()));
					oplist.push_back(add);

					for (int i=0; i<2; i++)
						for (int j=0; j<chunksX; j++)
							inPortMap (add, join("X",i*chunksX+j) , join("cp",i,j));
	
					outPortMap(add, "R", "addRes");
					vhdl << instance(add, "adder");
					syncCycleFromSignal("addRes");
					setCriticalPath(add->getOutputDelay("R"));
					vhdl<<tab<<"R<=addRes" << range(adderWidth-1,adderWidth-wInX_-wInY_) << ";" << endl;		
				}
		}else 
			vhdl << tab << "R <= X * Y ;" <<endl;
	}

		if (target_freq > maxDSPFrequency)
			target->setFrequency( target_freq );
	}

	IntMultiplier::~IntMultiplier() {
	}

	void IntMultiplier::outputVHDL(std::ostream& o, std::string name) {
		licence(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.std_logic_arith.all;" << endl;
		if (sign_)
			o << "use ieee.std_logic_signed.all;" << endl;
		else
			o << "use ieee.std_logic_unsigned.all;" << endl;
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


	void IntMultiplier::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
		if (! sign_){

			mpz_class svR = svX * svY;

			tc->addExpectedOutput("R", svR);
		}else{
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
}
