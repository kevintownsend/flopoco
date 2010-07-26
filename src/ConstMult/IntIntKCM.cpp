/*
 * An constant multiplier for FloPoCo using the KCM method
  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author :  Bogdan.Pasca, Florent.de.Dinechin, both @ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License, 2008-2010.

 */

// TODO if LUTsize=5 and wIn=11, we should have 2 tables, not 3

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
		Operator(target, inputDelays), wIn_(wIn), signedInput_(inputTwosComplement), C_(C), inputDelays_(inputDelays) 
	{
	
		setCopyrightString("Bogdan Pasca, Florent de Dinechin (2009,2010)");		
		// Set up the IO signals
		addInput ("X" , wIn_);
		

		if(C<0)
			throw string("IntIntKCM: only positive constants are supported");

		wOut_ = intlog2(C) + wIn_; 
				
		addOutput("R" , wOut_);

		ostringstream name;
		name << "IntIntKCM_" << wIn_ << "_" << C << (signedInput_?"_signed":"_unsigned");
		setName(name.str());

		setCriticalPath( getMaxInputDelays(inputDelays) );

		int constantWidth = intlog2( C ); 
		int lutWidth = target->lutInputs();
		chunkSize_ = lutWidth;
		int nbOfTables = int ( ceil( double(wIn)/double(lutWidth)) );
		int lastLutWidth = (wIn%lutWidth==0?lutWidth: wIn%lutWidth);

#if 0 // TODO need to hack table generation, too
		// Better to double an existing table than adding one more table and one more addition.
		if (lastLutWidth==1){ 
		  nbOfTables--;
		  lastLutWidth=lutWidth + 1;
		}
#endif
//			double delay = 0.0;
	
	
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
			useSoftRAM(t1);


			if(signedInput_) {
				t2 = new KCMTable(target, lutWidth, constantWidth + lutWidth, C, true);
				oplist.push_back(t2);

				useSoftRAM(t2);
			}

			manageCriticalPath(target->lutDelay() + 2*target->localWireDelay());

			//perform nbOfTables multiplications
 			for ( int i=0; i<nbOfTables; i++){

				if(signedInput_ && (i==nbOfTables-1)) {
					inPortMap (t2, "X", join("d",i));
					outPortMap(t2, "Y", join("pp",i));
					vhdl << tab << "-- last table multiplies by signed input" << endl;
					vhdl << instance(t2, join("KCMTable_",i));
				}
				else {
					inPortMap (t1, "X", join("d",i));
					outPortMap(t1, "Y", join("pp",i));
					vhdl << instance(t1, join("KCMTable_",i));
				}
			}
			
			for ( int i=0; i<nbOfTables; i++)
				syncCycleFromSignal(join("pp",i));
		
		
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
		
//			if(wIn_>32 && target->normalizedFrequency()>=0.5) { // TODO a real test, or fix the inMap?
//				map<string, double> inMap;
//				inMap["X0"] = delay;
		
				IntCompressorTree* adder = new IntCompressorTree(target, addOpSize, nbOfTables, inDelayMap("X0",getCriticalPath()));
				oplist.push_back(adder);
				for (int i=0; i<nbOfTables; i++)
					inPortMap (adder, join("X",i) , join("addOp",i));
		
//				inPortMapCst(adder, "Cin", "'0'");
				outPortMap(adder, "R", "OutRes");
				vhdl << instance(adder, "Result_Adder");
				syncCycleFromSignal("OutRes");
//			}

//			else{
//				if(target->normalizedFrequency()>0.5)
//					nextCycle();
//				// stupid addition
//				vhdl << tab << declare("OutRes", addOpSize) << " <= ";
//				for (int i=0; i<nbOfTables; i++) {
//					vhdl <<  join("addOp",i);
//					if(i<nbOfTables-1)
//						vhdl << " + ";
//				}
//				vhdl << ";" << endl;
//			}

			outDelayMap["R"] = adder->getOutputDelay("R");
		
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
