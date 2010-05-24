/*
 * An constant multiplier for FloPoCo using the KCM method
  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author :  Bogdan.Pasca, Florent.de.Dinechin, both @ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License, 2008-2010.

 */

#include <iostream>
#include <sstream>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "../IntNAdder.hpp"
#include "IntIntKCM.hpp"
#include "KCMTable.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator *> oplist;

	IntIntKCM::IntIntKCM(Target* target, int wIn, mpz_class C, bool inputTwosComplement, map<string, double> inputDelays):
		Operator(target), wIn_(wIn), signedInput_(inputTwosComplement), C_(C), inputDelays_(inputDelays) 
	{
		// Set up the IO signals
		addInput ("X" , wIn_);

		if(C<0)
			throw string("IntIntKCM: only positive constants are supported");

		wOut_ = intlog2(C) + wIn_; 
				
		addOutput("R" , wOut_);

		ostringstream name;
		name << "IntIntKCM_" << wIn_ << "_" << C << (signedInput_?"_signed":"_unsigned");
		setName(name.str());



		// if(0*inputTwosComplement){

 		// 	// cout << "C="<<unsignedBinary( (C<0? (0-C) :C) , intlog2( (C<0? (0-C) :C) ) +1  ) << endl;

		// 	// if (verbose){
		// 	// 	cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
		// 	// }
	
		// 	int constantWidth = intlog2( (C < 0?0-C:C) ) + 1; //the constant is represented in 2's complement
		// 	int lutWidth = target->lutInputs();
		// 	chunkSize_ = lutWidth;
		// 	int nbOfTables = int ( ceil( double(wIn)/double(lutWidth)) );
		// 	int lastLutWidth = (wIn%lutWidth==0?lutWidth: wIn%lutWidth);
		// 	double delay = 0.0;
	
		// 	KCMTable* t1; 
		// 	LastKCMTable* t2; //will containt the signed multiplication result 
	
		// 	// if (verbose){
		// 	// 	cout << "The width of the constant in 2's complement is = " << constantWidth<<endl;
		// 	// 	cout << "The number of inputs / LUT is = " << lutWidth << endl;
		// 	// 	cout << "The number of tables needed is = " << nbOfTables << endl;	
		// 	// }
	
		// 	if (nbOfTables>1) {
		// 		//let us build one table 
		// 		t1 = new KCMTable(target, lutWidth, constantWidth + lutWidth, C);
		// 		oplist.push_back(t1);
		// 		if (target->getVendor()=="Xilinx") 
		// 			{
		// 				addAttribute("rom_extract", "string", t1->getName()+": component", "yes");
		// 				addAttribute("rom_style", "string", t1->getName()+": component", "distributed");
		// 			}
		// 		else if (strncmp(typeid(*target).name(), "7Virtex", 7) == 0) // then the target is Virtex
		// 			addAttribute("altera_attribute", "string", t1->getName()+": component", "-name ALLOW_ANY_ROM_SIZE_FOR_RECOGNITION ON");
		// 	}
		
		// 	t2 = new LastKCMTable(target, lastLutWidth , constantWidth + lastLutWidth, C);
		// 	oplist.push_back(t2);

		// 	if (target->getVendor()=="Xilinx") 
		// 		{
		// 			addAttribute("rom_extract", "string", t2->getName()+": component", "yes");
		// 			addAttribute("rom_style", "string", t2->getName()+": component", "distributed");
		// 		}
		// 	else if (target->getVendor()=="Altera") 
		// 		addAttribute("altera_attribute", "string", t2->getName()+": component", "-name ALLOW_ANY_ROM_SIZE_FOR_RECOGNITION ON");
	
		// 	//first split the input X into digits having lutWidth bits -> this is as generic as it gets :)
		// 	for (int i=0; i<nbOfTables; i++)
		// 		if (i < nbOfTables-1)
		// 			vhdl << tab << declare( join("d",i), lutWidth ) << " <= X" << range(lutWidth*(i+1)-1, lutWidth*i ) << ";" <<endl;
		// 		else
		// 			vhdl << tab << declare( join("d",i), wIn -  lutWidth*i ) << " <= X" << range( wIn-1 , lutWidth*i ) << ";" <<endl;
	
		// 	// cout << "Generating the maps on the tables ... "<<endl;
		// 	//perform nbOfTables multiplications
		// 	for ( int i=nbOfTables-1; i >= 0; i--){
		// 		if (i==(nbOfTables-1)){
		// 			//instantiate the KCMLastTable
		// 			inPortMap (t2, "X", join("d",i));
		// 			outPortMap(t2, "Y", join("pp",i));
		// 			vhdl << instance(t2, join("LastKCMTable_SignedMul",i));
		// 			// cout << "Generated map for last table " << endl;
		// 		}else{
		// 			//instantiate the KCMLastTable
		// 			inPortMap (t1, "X", join("d",i));
		// 			outPortMap(t1, "Y", join("pp",i));
		// 			vhdl << instance(t1, join("KCMTable_UnsignedMul",i));
		// 			// cout << "Generated map for first table " << endl;
		// 		}
		// 	}
		
		// 	delay += target->lutDelay();
		
		// 	//determine the addition operand size
		// 	int addOpSize = (nbOfTables - 2) * lutWidth + (constantWidth + 	lastLutWidth);
		// 	//if (verbose)
		// 	// cout << "The addition operand size is: " << addOpSize << endl;
		
		// 	for (int i=0; i<nbOfTables; i++){
		// 		vhdl << tab << declare( join("addOp",i), addOpSize ) << " <= ";
		// 		if (i!=nbOfTables-1){ //if not the last table
		// 			for (int j=addOpSize-1; j>= (constantWidth + lutWidth) + (i-1)*lutWidth ; j--) //sign extension
		// 				vhdl << join("pp",i)<<of(constantWidth + lutWidth -1) << " & ";
					
		// 			vhdl << join("pp",i) << range(constantWidth + lutWidth -1, (i==0?lutWidth:0)) << " & " << zg((i-1)*lutWidth,0) << ";" << endl;
		// 		}
		// 		else{
		// 			for (int j=addOpSize-1; j>= (constantWidth + lastLutWidth)+ (i-1)*lutWidth ; j--)
		// 				vhdl << join("pp",i)<<range(constantWidth + lastLutWidth -1,constantWidth + lastLutWidth -1) << " & ";
					
		// 			vhdl << join("pp",i)<<range(constantWidth + lastLutWidth -1, (i==0?lutWidth:0))<< " & " << zg((i-1)*lutWidth,0) << ";" << endl;
		// 		}
		// 	}
		
		// 	map<string, double> inMap;
		// 	inMap["X0"] = delay;
		
		// 	IntNAdder* adder = new IntNAdder(target, addOpSize, nbOfTables, inMap);
		// 	oplist.push_back(adder);
		// 	for (int i=0; i<nbOfTables; i++)
		// 		inPortMap (adder, join("X",i) , join("addOp",i));
		
		// 	inPortMapCst(adder, "Cin", "'0'");
		// 	outPortMap(adder, "R", "OutRes");
		// 	vhdl << instance(adder, "Result_Adder");
		
		// 	syncCycleFromSignal("OutRes");
		
		// 	vhdl << tab << "R <= OutRes & pp0"<<range(lutWidth-1,0) << ";" <<endl;
		// }





		// else {

			// cout << "C="<<unsignedBinary( (C<0? (0-C) :C) , intlog2( (C<0? (0-C) :C) ) +1  ) << endl;

			// if (verbose==2){
			// 	cerr <<"> IntIntKCM:\t delay for X is   "<< inputDelays["X"]<<endl;	
			// }
	
			int constantWidth = intlog2( C ); 
			int lutWidth = target->lutInputs();
			chunkSize_ = lutWidth;
			int nbOfTables = int ( ceil( double(wIn)/double(lutWidth)) );
			int lastLutWidth = (wIn%lutWidth==0?lutWidth: wIn%lutWidth);
			double delay = 0.0;
	
	
			// if (verbose){
			// 	cerr << "> IntIntKCM:\t The width of the constant is = " << constantWidth<<endl;
			// 	cerr << "> IntIntKCM:\t The number of inputs / LUT is = " << lutWidth << endl;
			// 	cerr << "> IntIntKCM:\t The number of tables needed is = " << nbOfTables << endl;
			// }


			//first split the input X into digits having lutWidth bits -> this is as generic as it gets :)
			for (int i=0; i<nbOfTables; i++)
				if (i < nbOfTables-1)
					vhdl << tab << declare( join("d",i), lutWidth ) << " <= X" << range(lutWidth*(i+1)-1, lutWidth*i ) << ";" <<endl;
				else {
					vhdl << tab << declare( join("d",i), lutWidth ) << " <= " ;
					if(lutWidth*(i+1)>wIn){	
						if(signedInput_) {
							vhdl << "(" << lutWidth*(i+1)-1 << " downto " <<  wIn << "=> X(" << wIn-1 << ")) & ";
						}
						else
							vhdl << rangeAssign(lutWidth*(i+1)-1, wIn, "'0'") << " & ";
					}
				vhdl << "X" << range( wIn-1 , lutWidth*i ) << ";" <<endl;
				}
	
			KCMTable *t1, *t2; 
 			t1 = new KCMTable(target, lutWidth, constantWidth + lutWidth, C, false);
			oplist.push_back(t1);
 			if (target->getVendor() == "Xilinx") 
				{
					addAttribute("rom_extract", "string", t1->getName()+": component", "yes");
					addAttribute("rom_style", "string", t1->getName()+": component", "distributed");
				}
			if (target->getVendor() == "Altera") 
				addAttribute("altera_attribute", "string", t1->getName()+": component", "-name ALLOW_ANY_ROM_SIZE_FOR_RECOGNITION ON");


			if(signedInput_) {
				t2 = new KCMTable(target, lutWidth, constantWidth + lutWidth, C, true);
				oplist.push_back(t2);

				if (target->getVendor() == "Xilinx") 
					{
						addAttribute("rom_extract", "string", t2->getName()+": component", "yes");
						addAttribute("rom_style", "string", t2->getName()+": component", "distributed");
					}
				if (target->getVendor() == "Altera") 
					addAttribute("altera_attribute", "string", t2->getName()+": component", "-name ALLOW_ANY_ROM_SIZE_FOR_RECOGNITION ON");
			}


			//perform nbOfTables multiplications
 			for ( int i=0; i<nbOfTables; i++){

				if(signedInput_ && (i==nbOfTables-1)) {
					inPortMap (t2, "X", join("d",i));
					outPortMap(t2, "Y", join("pp",i));
					vhdl << instance(t2, join("KCMTable_",i));
				}
				else {
					inPortMap (t1, "X", join("d",i));
					outPortMap(t1, "Y", join("pp",i));
					vhdl << tab << "-- last table multiplies by signed input" << endl;
					vhdl << instance(t1, join("KCMTable_",i));
				}
			}
		
			delay += target->lutDelay();
		
			//determine the addition operand size
			int addOpSize = (nbOfTables - 2) * lutWidth + (constantWidth +  lastLutWidth);
			// if (verbose)
			// 	cerr << "> IntIntKCM:\t The addition operand size is: " << addOpSize << endl;
		
			for (int i=0; i<nbOfTables; i++){
				vhdl << tab << declare( join("addOp",i), addOpSize ) << " <= ";
				if (i!=nbOfTables-1){ //if not the last table
					vhdl << rangeAssign(addOpSize-1, constantWidth + i*lutWidth, "'0'") << " & " 
						  <<  join("pp",i) << range(constantWidth + lutWidth -1, (i==0?lutWidth:0)) << " & " 
						  << zg((i-1)*lutWidth,0) << ";" << endl;
				}
				else{
					vhdl << join("pp",i)<<range(constantWidth + lastLutWidth -1, (i==0?lutWidth:0))<< " & " << zg((i-1)*lutWidth,0) << ";" << endl;
				}
			}
		
			map<string, double> inMap;
			inMap["X0"] = delay;
		
			IntNAdder* adder = new IntNAdder(target, addOpSize, nbOfTables, inMap);
			oplist.push_back(adder);
			for (int i=0; i<nbOfTables; i++)
				inPortMap (adder, join("X",i) , join("addOp",i));
		
			inPortMapCst(adder, "Cin", "'0'");
			outPortMap(adder, "R", "OutRes");
			vhdl << instance(adder, "Result_Adder");
		
			syncCycleFromSignal("OutRes");
		
			vhdl << tab << "R <= OutRes & pp0" << range(lutWidth-1,0) << ";" <<endl;
	}

	IntIntKCM::~IntIntKCM() {
	}

	int IntIntKCM::getOutputWidth(){
		return wOut_;
	}

	void IntIntKCM::emulate(TestCase* tc)
	{
		mpz_class svX =  tc->getInputValue("X");
		// cout << "-------------------------"<<endl;
		// cout << "X="<<unsignedBinary(svX, wIn_);
		bool xneg = false;
		//x is in 2's complement, so it's value is
		if(signedInput_) {
			if ( svX > ( (mpz_class(1)<<(wIn_-1))-1) ){
				// cout << "X is negative" << endl;
				xneg = true;
				svX = svX - (mpz_class(1)<<wIn_);
			}
		}
		mpz_class svR;
	
		svR = svX * C_;

		if ( svR < 0)
			svR = (mpz_class(1)<<wOut_) + svR;
	
		// cout << "R="<<unsignedBinary(svR, wOut_);

		tc->addExpectedOutput("R", svR);
	}


}
