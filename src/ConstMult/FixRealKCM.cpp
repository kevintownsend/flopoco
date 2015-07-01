/*
 * A faithful multiplier by a real constant, using a variation of the KCM method
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 * 
 * Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr
 *
 * Initial software.
 * Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
 * 2008-2011.
 * All rights reserved.
 */

// TODO the first table should have lutSize input bits!

#include "../Operator.hpp"

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <sollya.h>
#include "../utils.hpp"
#include "FixRealKCM.hpp"
#include "../IntAddSubCmp/IntAdder.hpp"

using namespace std;

namespace flopoco{

	/**
	* @brief init : all operator initialisation stuff goes here
	*/
	void FixRealKCM::init()
	{
	srcFileName="FixRealKCM";

		if(lsbIn>msbIn) 
			throw string("FixRealKCM: Error, lsbIn>msbIn");
    
		if(targetUlpError > 1.0)
			THROWERROR("FixRealKCM: Error, targetUlpError="<<
					targetUlpError<<">1.0. Should be in ]0.5 ; 1].");
		//Error on final rounding is er <= 2^{lsbout - 1} = 1/2 ulp so 
		//it's impossible to reach 0.5 ulp of precision
		if(targetUlpError <= 0.5) 
			THROWERROR("FixRealKCM: Error, targetUlpError="<<
					targetUlpError<<"<0.5. Should be in ]0.5 ; 1]");
		
		int signBit=0;
		if(signedInput)
			signBit=1;
		wIn+=signBit;

		// Convert the input string into a sollya evaluation tree
		sollya_obj_t node;
		node = sollya_lib_parse_string(constant.c_str());	
		/* If  parse error throw an exception */
		if (sollya_lib_obj_is_error(node))
		{
				ostringstream error;
				error << srcFileName << ": Unable to parse string "<< 
					constant << " as a numeric constant" <<endl;
				throw error.str();
		}

		mpfr_init2(mpC, 10000);
		mpfr_init2(absC, 10000);
		sollya_lib_get_constant(mpC, node);

		//if negative constant, then set negativeConstant and remake the
		//constant positive
		negativeConstant = false;
		if(mpfr_cmp_si(mpC, 0) < 0)
		{
			//throw string("FixRealKCMBH: only positive constants are supported");
			negativeConstant = true;
		}

		mpfr_abs(absC, mpC, GMP_RNDN);

		REPORT(DEBUG, "Constant evaluates to " << mpfr_get_d(mpC, GMP_RNDN));

		// build the name
		ostringstream name; 
		name <<"FixRealKCM_" << vhdlize(lsbIn)  << "_" << vhdlize(msbIn) << 
			"_" << vhdlize(lsbOut) << "_" << vhdlize(constant)  << 
			(signedInput  ?"_signed" : "_unsigned");
		setName(name.str());

		mpfr_t log2C;
		mpfr_init2(log2C, 100); // should be enough for anybody
		mpfr_log2(log2C, absC, GMP_RNDN);
		msbC = mpfr_get_si(log2C, GMP_RNDU);
		mpfr_clears(log2C, NULL);

		msbOut = msbC + msbIn;
		wOut = msbOut + signBit - lsbOut+1;

		REPORT(DEBUG, "msbConstant=" << msbC  << "   (msbIn,lsbIn)=("<< 
				msbIn << "," << lsbIn << ")   wIn=" << wIn << 
				"   (msbOut,lsbOut)=(" << msbOut << "," << lsbOut <<
				")   wOut=" << wOut
			);

		g = neededGuardBits(getTarget(), wIn, targetUlpError);
	}
	
	int FixRealKCM::computeTableNumbers(
			Target* target,
			int wIn,
			int** disize_target
		)
	{

#ifndef WIP_FORGET
#define WIP_FORGET
#endif
#ifdef WIP_FORGET
		int lutWidth = target->lutInputs(); 
		int* diSize = new int[17*42];		
		int currentSize;
		currentSize = diSize[0] = (wIn < lutWidth) ? wIn : lutWidth;
		diSize[1] = 0;
		int nbOfTables = 1;
		size_t idx = 1;

		int remainingBits = (wIn - diSize[0]) % (lutWidth - 1);

		size_t i;
		for(i = idx ; currentSize < wIn - remainingBits ; ++i)
		{
			currentSize += diSize[i] = ((wIn - currentSize) < (lutWidth - 1))?
				wIn - currentSize : lutWidth - 1;
			nbOfTables++;
		}

		if(remainingBits <= lutWidth/2)
		{
			diSize[1] += remainingBits;
			if(nbOfTables < 2)
			{
				nbOfTables++;
			}
		}
		else
		{
			diSize[idx] = remainingBits;
			nbOfTables++;
		}

#else
		int lutWidth = target->lutInputs(); 
	
		// First set up all the sizes. Table 0 is the most one
		int nbOfTables = 0;
		int* diSize = new int[17*42];		

		diSize[0] = lutWidth;
		int currentSize = diSize[0];
		int counter=1;

		while(currentSize < wIn) {
			// -1 because the tools are able to pack LUT + addition in one LUT 
			diSize[counter] = lutWidth-1; 

			currentSize += diSize[counter];
			counter++;
		}

		nbOfTables = counter;
		counter--;

		//Last lut
		diSize[counter] = wIn - (currentSize - diSize[counter]);

		//Better to move the remaining bits to the first tables, than to have
		//them in a new table
		if (diSize[counter] <= lutWidth/2)
		{
			diSize[1] += diSize[counter];
			
			counter=2;
			currentSize = diSize[0] + diSize[1];
			while(currentSize < wIn) 
			{
				diSize[counter] = lutWidth-1;
				currentSize += diSize[counter];
				counter++;	
			}
			nbOfTables = counter;
			counter--;
			diSize[counter] = wIn - (currentSize - diSize[counter]);
		}
#endif
		if(disize_target != nullptr)
		{
			*disize_target = diSize;
		}
		else
		{
			delete[] diSize;
		}
		return nbOfTables;
	}

	//standalone operator
	FixRealKCM::FixRealKCM(
				Target* target, 
				bool signedInput_, 
				int msbIn_, 
				int lsbIn_, 
				int lsbOut_, 
				string constant_, 
				double targetUlpError_,
				map<string, double> inputDelays 
			):Operator(target, inputDelays), 
			signedInput(signedInput_),
			msbIn(msbIn_), 
			lsbIn(lsbIn_), 
			wIn(msbIn_-lsbIn_+1), 
			lsbOut(lsbOut_), 
			constant(constant_), 
			targetUlpError(targetUlpError_)
	{
		init();		
		
		int* diSize;
		int nbOfTables = computeTableNumbers(target, wIn, &diSize);
		int lutWidth = target->lutInputs();
		
		REPORT(INFO, "Constant multiplication in "<< nbOfTables << " tables");
		//manage the critical path
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		addInput("X", wIn);
		addOutput("R", wOut);

#define WIP_FORGET 
#ifdef WIP_FORGET
		int* doSize;

		FixRealKCMTable** t = createTables(
				target,
				diSize,
				nbOfTables,
				&doSize
			);


		setCriticalPath(getMaxInputDelays(inputDelays));

		if(!target->plainVHDL())
		{
			//create the bitheap
			bitHeap = new BitHeap(this, wOut+g);
		}

		manageCriticalPath(target->localWireDelay() + target->lutDelay());

		for(int i = 0; i < nbOfTables; i++)
		{
			REPORT(DEBUG, "Adding bits for table " << i)
				//manage the critical path
				setCycleFromSignal(join("d",i));

			inPortMap (t[i] , "X", join("d",i));
			outPortMap(t[i] , "Y", join("pp",i));
			vhdl << instance(t[i] , join("KCMTable_",i));

			//manage the critical path
			syncCycleFromSignal(join("pp",i));

			//add the bits to the bit heap
			for(int w=0; w<=doSize[i]-1; w++)
			{
				stringstream s;

				//manage the critical path
				manageCriticalPath(target->lutDelay());

				if(negativeConstant != (w==(doSize[i]-1) && i==(nbOfTables-1)))
					s << "not(pp" << i << of(w) << ")";
				else
					s << "pp" << i << of(w);

				bitHeap->addBit(w, s.str());
			}
		}

		//compress the bitheap and produce the result
		bitHeap->generateCompressorVHDL();

		//manage the critical path
		syncCycleFromSignal(bitHeap->getSumName());

		//because of final add in bit heap, add one more bit to the result
		vhdl << declare("OutRes", wOut+g) << " <= " << 
			bitHeap->getSumName() << range(wOut+g-1, 0) << ";" << endl;	

#else
		
		if(wIn <= lutWidth+1) //Marche bien 
		{
			///////  multiplication using 1 table only ////////////////////////
			REPORT(INFO, 
				"Constant multiplication in a single table,"
				"will be correctly rounded");

			FixRealKCMTable *t;
			t = new FixRealKCMTable(
					target, 
					this, 
					0, 
					0, 
					wIn, 
					wOut, 
					signedInput, 
					false
				);

			addSubComponent(t);
			useSoftRAM(t);

			//manage the critical path
			manageCriticalPath(target->localWireDelay() + target->lutDelay());
      
			inPortMap (t , "X", "X");
			outPortMap(t , "Y", "Y");
			vhdl << instance(t , "KCMTable");
			
			//manage the critical path
			syncCycleFromSignal("Y");

			//negate the result if necessary
			if(negativeConstant)
			{
				//manage the critical path
				manageCriticalPath(target->localWireDelay() + target->lutDelay());
				
				vhdl << tab << declare("Y_xored", wOut) << " <= Y xor " << 
					og(wOut) << ";" << endl;
				vhdl << tab << declare("Y_negated", wOut) << 
					"<= Y_xored + (" << zg(wOut-1) << " & \'1\');" << endl;
				
				vhdl << tab << "R <= Y_negated;" << endl;
			}
			else
			{
				vhdl << tab << "R <= Y;" << endl;
			}
			
		  	outDelayMap["R"] = getCriticalPath();
		}
		else
		{
			//////////////////   Generic Case  ///////////////////////
			// All the tables are read in parallel
			setCriticalPath(getMaxInputDelays(inputDelays));

			int ppiSize[42*17]; // should be more than enough for everybody
			FixRealKCMTable *t[17*42]; 
		
			//first split the input X into digits having lutWidth bits -> this
			//is as generic as it gets :)
			bool tableSigned=false;
			bool last;
			int highBit = wIn;
			int ppI = wOut+g;
			
			for (int i=nbOfTables-1; i>=0; i--)
			{
				// The bit width of the output of this table
				// The last table has to have wOut+g  bits.
				// The previous one wOut+g-lastLutWidth
				// the previous one wOut+g-lastLutWidth-lutWidth etc

				vhdl << tab << declare( join("d",i), diSize[i] ) << " <= X" << 
					range(highBit-1,   highBit - diSize[i]) << ";" << endl;
				highBit -= diSize[i];
				ppiSize[i] = ppI;
				ppI -= diSize[i];
				if(i < nbOfTables-1)
				{
					tableSigned=false;
					last=false;
				}
				else 
				{
					// last table is a bit special
					if(signedInput)
						tableSigned=true;
					last=true;
				}

				REPORT(DEBUG, "Table i=" << i << ", input size=" << 
						diSize[i] << ", output size=" << ppiSize[i]);

				// Now produce the VHDL
				
				// already updated 
				t[i] = new FixRealKCMTable(
							target, 
							this, 
							i, 
							highBit, 
							diSize[i], 
							ppiSize[i], 
							tableSigned, 
							last, 
							1
						);
				useSoftRAM(t[i]);
				addSubComponent(t[i]);
            	

#define USE_MADDER 0
#if USE_MADDER // size= 190
				vhdl << tab << declare( join("addOp",i), wOut+g ) << " <= ";
				if (i!=nbOfTables-1) //if not the last table
					vhdl << rangeAssign(wOut+g-1, ppiSize[i], "'0'") << " & " ;	
				vhdl << join("pp",i) << ";" << endl;
#endif
			}
			
			Operator* adder;
			
#if USE_MADDER
			if(nbOfTables>2)
			{
				adder = new IntMultiAdder(target, wOut+g, nbOfTables, inDelayMap("X0",target->localWireDelay() + getCriticalPath()));
				addSubComponent(adder);
				for (int i=0; i<nbOfTables; i++)
					inPortMap (adder, join("X",i) , join("addOp",i));
			}
			else
			{
				// 2 tables only
				adder = new IntAdder(target, wOut+g, inDelayMap("X",target->localWireDelay() + getCriticalPath()));
				addSubComponent(adder);
				inPortMap (adder, "X" , join("addOp",0));
				inPortMap (adder, "Y" , join("addOp",1));
				inPortMapCst(adder, "Cin" , "'0'");
			}
			outPortMap(adder, "R", "OutRes");
			vhdl << instance(adder, "Result_Adder");
			syncCycleFromSignal("OutRes");
			setCriticalPath( adder->getOutputDelay("R") );

#else // Back to the rake!
			if(!target->plainVHDL())
			{
				//create the bitheap
				bitHeap = new BitHeap(this, wOut+g);
			}

			if(nbOfTables>2)
			{
				if(!target->plainVHDL())
				{
					//manage the pipeline
					manageCriticalPath(target->localWireDelay() + target->lutDelay());
					
					for(int i = 0; i < nbOfTables; i++)
					{
						REPORT(DEBUG, "Adding bits for table " << i)
						//manage the critical path
						setCycleFromSignal(join("d",i));
						
						inPortMap (t[i] , "X", join("d",i));
						outPortMap(t[i] , "Y", join("pp",i));
						vhdl << instance(t[i] , join("KCMTable_",i));
						
						//manage the critical path
						syncCycleFromSignal(join("pp",i));

						//add the bits to the bit heap
						for(int w=0; w<=ppiSize[i]-1; w++)
						{
							stringstream s;
							
							//manage the critical path
							manageCriticalPath(target->lutDelay());
							
							if(negativeConstant != (w==(ppiSize[i]-1) && i==(nbOfTables-1)))
								s << "not(pp" << i << of(w) << ")";
							else
								s << "pp" << i << of(w);
							
							bitHeap->addBit(w, s.str());
						}
						
						// Sign extension
						if(negativeConstant   ||   (i==(nbOfTables-1)  && !negativeConstant) )
						{
							for(int w=( 
 										(i==(nbOfTables-1)) ? 
											ppiSize[i]-1 : 
											ppiSize[i]
										); 
									w < wOut+g; 
									w++)
							{
								bitHeap->addConstantOneBit(w);
							}
							
							if(negativeConstant){}
								bitHeap->addConstantOneBit(0);
						}
					}
					//add one bit at weight g-1, for faithfull rounding
					bitHeap->addConstantOneBit(g-1);				
					
					//compress the bitheap and produce the result
					bitHeap->generateCompressorVHDL();
					
					//manage the critical path
					syncCycleFromSignal(bitHeap->getSumName());
						
					//because of final add in bit heap, add one more bit to the result
					vhdl << declare("OutRes", wOut+g) << " <= " << 
						bitHeap->getSumName() << range(wOut+g-1, 0) << ";" << endl;
				}
				else // Target->plainVHDL() is true 
				{
					bool pipeinit;
					bool pipe[17*42];
					// First evaluate the pipeline by hand, because
					// PipelineMadeEasy(TM) overestimates; 
					double slack = 1/target->frequency() - getCriticalPath();

					if (slack<0)
					{
						slack = 1/target->frequency() - getTarget()->ffDelay();
						pipeinit = true;
					}
					else
						pipeinit = false;
					 
					slack -=  target->localWireDelay(ppiSize[0]) + target->lutDelay();
					if (slack<0)
					{
						slack = 1/target->frequency() - getTarget()->ffDelay();
						pipe[0] = true;
					}
					else
						pipe[0] = false;

					// Second table + addition costs 1 adder delay
					slack -= target->localWireDelay(ppiSize[1]) + 
						target->adderDelay(ppiSize[1]);
					if (slack<0)
					{
						slack = 1/target->frequency()- getTarget()->ffDelay();
						pipe[1] = true;
					}
					else
						pipe[1] = false;
					
					// Following tables add just local routing and the delay of
					// one small addition (this includes the additional LUT
					// delay)
					for (int i=2; i<nbOfTables; i++)
					{
						slack -= target->adderDelay(target->lutInputs()-1) + 
							target->localWireDelay(ppiSize[i]);
						if (slack<0)
						{
							slack = 1/target->frequency()- getTarget()->ffDelay();
							pipe[i] = true;
						}
						else
							pipe[i] = false;
					}

					// Now we may build the pipeline
					if(pipeinit)
						nextCycle();
					
					int i=0;
					inPortMap (t[i] , "X", join("d",i));
					outPortMap(t[i] , "Y", join("pp",i));
					vhdl << instance(t[i] , join("KCMTable_",i));
					
					if(pipe[0])
						nextCycle();
					
					vhdl << tab << declare("sum0", ppiSize[0]) << " <= pp0;" << endl;
					
					for (i=1; i<nbOfTables; i++)
					{
						if(pipe[i])
							nextCycle();
						
						inPortMap (t[i] , "X", join("d",i));
						outPortMap(t[i] , "Y", join("pp",i));
						vhdl << instance(t[i] , join("KCMTable_",i));
						
						// TODO prove that this addition never overflows, or add
						// +1 
						vhdl << tab << declare(join("sum",i), ppiSize[i]) 
							 << " <= " << join("pp",i) << " + " << 
							 join("sum", i-1) << ";" << endl;
					}
					vhdl << declare("OutRes", ppiSize[nbOfTables-1]) << 
						" <= " << join("sum", nbOfTables-1) << ";" << endl;
					setCriticalPath( 
							1.0/getTarget()->frequency() - 
							getTarget()->ffDelay() - 
							slack 
						);
				}
			}
			else 
			{
				// 2 tables only
				
				//manage the pipeline
				manageCriticalPath(target->localWireDelay() + target->lutDelay());
				
				for(int i=0; i<nbOfTables; i++) {
					inPortMap (t[i] , "X", join("d",i));
					outPortMap(t[i] , "Y", join("pp",i));
					vhdl << instance(t[i] , join("KCMTable_",i));
				}
				
				//manage the critical path
				syncCycleFromSignal(join("pp", 1));
				manageCriticalPath(target->adderDelay(wOut+g));
				
				vhdl << tab << declare("addOp0", wOut+g ) << " <= " << 
					rangeAssign(wOut+g-1, ppiSize[0], "'0'") << " & pp0;" << endl;
				
				adder = new IntAdder(
						target, 
						wOut+g, 
						inDelayMap(
							"X",  
							target->localWireDelay() + 
								getCriticalPath() + 
								target->localWireDelay(ppiSize[1]) + 
								target->lutDelay()  
							)
					);
				addSubComponent(adder);
				
				inPortMap (adder, "X" , "addOp0");
				inPortMap (adder, "Y" , "pp1");
				inPortMapCst(adder, "Cin" , "'0'");
				if(negativeConstant)
				{
					outPortMap(adder, "R", "OutRes_int");
				}
				else
					outPortMap(adder, "R", "OutRes");
				vhdl << instance(adder, "Result_Adder");
				
				if(negativeConstant)
				{
					syncCycleFromSignal("OutRes_int");
					
					vhdl << tab << declare("OutRes_int_xored", wOut+g) << 
						" <= OutRes_int xor " << og(wOut+g) << ";" << endl;
					
					inPortMap	 (adder, "X" , "OutRes_int_xored");
					inPortMapCst (adder, "Y" , "\'0\'");
					inPortMapCst (adder, "Cin" , "\'1\'");
					outPortMap(adder, "R", "OutRes");
				}
				
				syncCycleFromSignal("OutRes");
			}

			
#endif
			vhdl << tab << "R <= OutRes" << range(wOut+g-1, g) << ";" << endl;
			outDelayMap["R"] = getCriticalPath();
		}
#endif
		vhdl << tab << "R <= OutRes" << range(wOut+g-1, g) << ";" << endl;
		outDelayMap["R"] = getCriticalPath();
		delete[] diSize;
	}
	
	
	//operator incorporated into a global compression
	//	for use as part of a bigger operator
	FixRealKCM::FixRealKCM(
			Operator* parentOp, 
			Target* target, 
			Signal* multiplicandX, 
			bool signedInput_, 
			int msbIn_, 
			int lsbIn_, 
			int lsbOut_, 
			string constant_, 
			BitHeap* bitHeap_, 
			double targetUlpError_, 
			map<string, double> inputDelays
		):	
			Operator(target, inputDelays),
			signedInput(signedInput_), 
			msbIn(msbIn_),
			lsbIn(lsbIn_), 
			wIn(msbIn_-lsbIn_+1), 
			lsbOut(lsbOut_), 
			constant(constant_), 
			targetUlpError(targetUlpError_), 
			bitHeap(bitHeap_)
	{
		
		init();
		
		// -1 because the tools are able to pack LUT + addition in one LUT
		int lutWidth = target->lutInputs(); 

		// First set up all the sizes
		int *diSize;
		int nbOfTables = computeTableNumbers(target, wIn, &diSize);
		

		REPORT(INFO, "Constant multiplication in "<< nbOfTables << " tables");		
		
		//manage pipeline
		parentOp->syncCycleFromSignal(multiplicandX->getName());

		if(wIn <= lutWidth+1)
		{
			//////////////// multiplication using 1 table only ////////////////
			REPORT(INFO, "Constant multiplication in a single table, will be correctly rounded");
			int bitheapSize = bitHeap->getMaxWeight() - bitHeap->getMinWeight()  + 1;

			FixRealKCMTable *t;
			t = new FixRealKCMTable(target, this, 0, 0, wIn, wOut, signedInput, false);
			parentOp->addSubComponent(t);
			parentOp->useSoftRAM(t);

			//manage critical path
			parentOp->manageCriticalPath(
					target->localWireDelay() + target->lutDelay()
				);
			
			parentOp->inPortMap (t , "X", multiplicandX->getName());
			parentOp->outPortMap(t , "Y", join("Y_kcmMult_", getuid()));
			parentOp->vhdl << 
				parentOp->instance(t , join("KCMTable_kcmMult_", getuid()));
			
			//manage pipeline
			parentOp->syncCycleFromSignal(join("Y_kcmMult_", getuid()));

			//add the resulting bits to the bit heap
			int ySize;
			
			ySize = parentOp->getSignalByName(join("Y_kcmMult_", getuid()))->width();
			//all but the msb, which is handled separately
			for(int w=0; w<=ySize-2; w++)
			{
				stringstream s;
				
				if(negativeConstant)
					s << "not(" << join("Y_kcmMult_", getuid()) << of(w) << ")";
				else
					s << join("Y_kcmMult_", getuid()) << of(w);
				
				bitHeap->addBit(w, s.str());
				
				if(negativeConstant)
				{
					for(int w2=w; w2<bitheapSize; w2++)
						bitHeap->addConstantOneBit(w2);
				}
			}
			
			//add the msb, complemented
			{
				stringstream s;
				
				s << "not(" << join("Y_kcmMult_", getuid()) << of(ySize-1) << ")";
				bitHeap->addBit(ySize-1, s.str());
			}
			
			//add the rest of the constant bits
			for(int w=ySize-1; w<bitheapSize; w++)
				bitHeap->addConstantOneBit(w);
		}
		else
		{
			///////////////////////   Generic Case  ///////////////////////////

			//manage pipeline
			parentOp->syncCycleFromSignal(multiplicandX->getName());

			int ppiSize[42*17]; // should be more than enough for everybody
			FixRealKCMTable *t[17*42];
			//first split the input X into digits having lutWidth bits -> this
			//is as generic as it gets :)
			bool tableSigned=false;
			bool last;
			int highBit = wIn;
			int ppI = wOut+g;
			for (int i=nbOfTables-1; i>=0; i--) 
			{
				// The bit width of the output of this table
				// The last table has to have wOut+g  bits.
				// The previous one wOut+g-lastLutWidth
				// the previous one wOut+g-lastLutWidth-lutWidth etc

				parentOp->vhdl << tab << 
					parentOp->declare(join("d", i, "_kcmMult_", getuid()), diSize[i]) 
					<< " <= " << multiplicandX->getName() << 
					range(highBit-1, highBit-diSize[i]) << ";" << endl;
				highBit -= diSize[i];
				ppiSize[i] = ppI;
				ppI -= diSize[i];
				if (i < nbOfTables-1){
					tableSigned=false;
					last=false;
				}
				else {// last table is a bit special
					if(signedInput)
						tableSigned=true;
					last=true;
				}

				REPORT(DEBUG, "Table i=" << i << ", input size=" << diSize[i] << 
						", output size=" << ppiSize[i]);

				// Now produce the VHDL
				t[i] = new FixRealKCMTable(
								target, 
								this, 
								i, 
								highBit, // already updated
								diSize[i], 
								ppiSize[i], 
								tableSigned, 
								last, 
								1
							);

				parentOp->useSoftRAM(t[i]);
				parentOp->addSubComponent(t[i]);
			}
			
			int bitheapSize = bitHeap->getMaxWeight()-bitHeap->getMinWeight();
			
			REPORT(DEBUG, "working with bitheap of size=" << bitheapSize);
			
			for(int i=0; i < nbOfTables; i++)
			{
				//manage pipeline
				parentOp->setCycleFromSignal(join("d", i, "_kcmMult_", getuid()));
				parentOp->manageCriticalPath(target->lutDelay());
				
				parentOp->inPortMap (t[i] , "X", join("d", i, "_kcmMult_", getuid()));
				parentOp->outPortMap(t[i] , "Y", join("pp", i, "_kcmMult_", getuid()));
				parentOp->vhdl << parentOp->instance(t[i] , 
						join("KCMTable_", i, "_kcmMult_", getuid()));
				
				//manage pipeline
				parentOp->syncCycleFromSignal(join("pp", i, "_kcmMult_", getuid()));
				
				//add the bits to the bit heap
				for(int w=0; w<=ppiSize[i]-1; w++)
				{
					stringstream s;
					
					if(negativeConstant != (w==(ppiSize[i]-1) && i==(nbOfTables-1)))
					{
						s << "not(pp" << i << "_kcmMult_" << getuid() << of(w) << ")";
					}
					else
						s << "pp" << i << "_kcmMult_" << getuid() << of(w);
					
					bitHeap->addBit(w, s.str());
				}
				
				if(negativeConstant || (i==(nbOfTables-1)  && !negativeConstant))
				{
					for(
							int w=(((negativeConstant && i==(nbOfTables-1)) ||
									(i==(nbOfTables-1))) ? 
										ppiSize[i]-1 :
										ppiSize[i]); 
							w<bitheapSize; 
							w++
					)
					{
						bitHeap->addConstantOneBit(w);
					}//End for
					
					if(negativeConstant)	
						bitHeap->addConstantOneBit(0);
				}
				
			}
		}
		delete[] diSize;
	}

	FixRealKCMTable** FixRealKCM::createTables(
			Target* target,
			int* diSize,
			int nbOfTables,
			int** doSize_target	
		)
	{
		FixRealKCMTable** t = new FixRealKCMTable*[nbOfTables]; 
		int* doSize = new int[17*42];

		//first split the input X into digits having lutWidth bits -> this
		//is as generic as it gets :)
		bool signedOutput = negativeConstant || signedInput;
		bool last = true;
		int highBit = wIn;
		int tableDo = wOut+g;
		
		for (int i=nbOfTables-1; i>=0; i--)
		{
			// The last table has to have wOut+g  bits.
			// The previous one wOut+g-lastLutWidth
			// the previous one wOut+g-anteLastLutWidth-lastlutWidth etc
			
			vhdl << tab << declare( join("d",i), diSize[i] ) << " <= X" << 
				range(highBit-1,   highBit - diSize[i]) << ";" << endl;
			highBit -= diSize[i];
			doSize[i] = tableDo;
			if(!(last && signedInput))
			{
				doSize[i]++;
			}
			tableDo -= diSize[i];

			REPORT(DEBUG, "Table i=" << i << ", input size=" << 
					diSize[i] << ", output size=" << doSize[i] << 
					", weight= " << highBit);

			// Now produce the VHDL

			// already updated 
			t[i] = new FixRealKCMTable(
					target, 
					this, 
					i, 
					highBit, 
					diSize[i], 
					doSize[i], 
					signedOutput, 
					last 
				);
			useSoftRAM(t[i]);
			addSubComponent(t[i]);

			last = false;
		}


		if(doSize_target != nullptr)
		{
			*doSize_target = doSize;
		}
		else
		{
			delete[] doSize;
		}

		return t;

	}

	// To have MPFR work in fix point, we perform the multiplication in very
	// large precision using RN, and the RU and RD are done only when converting
	// to an int at the end.
	void FixRealKCM::emulate(TestCase* tc)
	{
		// Get I/O values
		mpz_class svX = tc->getInputValue("X");
		bool negativeInput = false;
		
		// get rid of two's complement
		if(signedInput)
		{
			if ( svX > ( (mpz_class(1)<<(wIn-1))-1) )
			{
				svX = svX - (mpz_class(1)<<wIn);
				negativeInput = true;
			}
		}
		
		// Cast it to mpfr 
		mpfr_t mpX; 
		mpfr_init2(mpX, msbIn-lsbIn+2);	
		mpfr_set_z(mpX, svX.get_mpz_t(), GMP_RNDN); // should be exact
		
		// scale appropriately: multiply by 2^lsbIn
		mpfr_mul_2si(mpX, mpX, lsbIn, GMP_RNDN); //Exact
		
		// prepare the result
		mpfr_t mpR;
		mpfr_init2(mpR, 10*wOut);
		
		// do the multiplication
		mpfr_mul(mpR, mpX, absC, GMP_RNDN);
		
		// scale back to an integer
		mpfr_mul_2si(mpR, mpR, -lsbOut, GMP_RNDN); //Exact
		mpz_class svRu, svRd;
		
		mpfr_get_z(svRd.get_mpz_t(), mpR, GMP_RNDD);
		mpfr_get_z(svRu.get_mpz_t(), mpR, GMP_RNDU);

		if(negativeInput != negativeConstant)
		{
			svRd = (mpz_class(1) << wOut) - svRd;
			svRu = (mpz_class(1) << wOut) - svRu;
		}

		//Border cases
		if(svRd > (mpz_class(1) << wOut) - 1 )
		{
			svRd = 0;
		}

		if(svRu > (mpz_class(1) << wOut) - 1 )
		{
			svRu = 0;
		}

		tc->addExpectedOutput("R", svRd);
		tc->addExpectedOutput("R", svRu);

		// clean up
		mpfr_clears(mpX, mpR, NULL);
	}

	int FixRealKCM::neededGuardBits(
			Target* target, 
			int wIn, 
			double targetUlpError
		)
	{
		int guardBits;
		int nbOfTables = computeTableNumbers(target, wIn, nullptr);

		if(targetUlpError > 1.0 || targetUlpError <= 0.5)
		{
			cout << "ERREUR !!!!" << endl;
			return -1; 
		}
				
		if((nbOfTables <= 2 && targetUlpError==1.0) || nbOfTables == 1)
		{
			// specific case: two CR table make up a faithful sum
			guardBits = 0; 
		}
		else
		{
			/* Since we want the output error to be strictly less than 
			 * 2^(lsbout) * targetUlpError and we have 
			 */ 
			guardBits = ceil(log2((double)nbOfTables/(2*targetUlpError - 1)));
		}
			
		return guardBits;
	}



	/************************** The FixRealKCMTable class ********************/


	FixRealKCMTable::FixRealKCMTable(
			Target* target, 
			FixRealKCM* mother, 
			int i, 
			int weight, 
			int wIn, 
			int wOut, 
			bool signedOutput, 
			bool last, 
			int pipeline
		):
			Table(target, wIn, wOut, 0, -1, pipeline), 
			mother(mother), 
			index(i), 
			weight(weight), 
			signedOutput(signedOutput), 
			last(last)
	{
		ostringstream name; 
		srcFileName="FixRealKCM";
		name << mother->getName() << "_Table_" << index;
		setName(name.str());
	}
  
	mpz_class FixRealKCMTable::function(int x0)
	{
#ifdef WIP_FORGET
#pragma message("Using WIP_FORGET FixRealKCMTable function")
		int x;
		bool negativeInput = false;
		
		// get rid of two's complement
		x = x0;
		//Only the last "digit" has a negative weight
		if(mother->signedInput && last)
		{
			if ( x0 > ((1<<(wIn-1))-1) )
			{
				x = (1<<wIn) - x;
				negativeInput = true;
			}
		}

		mpz_class result;
		mpfr_t mpR, mpX;

		mpfr_init2(mpR, 10*wOut);
		mpfr_init2(mpX, wIn);

		if(x0 == 0)
		{
			mpfr_set_si(mpX, 0, GMP_RNDN); // should be exact
		}
		else
		{
			mpfr_set_si(mpX, x, GMP_RNDN); // should be exact
			//Getting real weight
			mpfr_mul_2si(mpX, mpX, weight, GMP_RNDN); //Exact

			// do the mult in large precision
			mpfr_mul(mpR, mpX, mother->absC, GMP_RNDN);
			
			int rescaleShift = (last && mother->signedInput) ? 0 : 1;

			// Result is integer*C, which is more or less what we need: just
			// scale to add g bits.
			mpfr_mul_2si(
					mpR, 
					mpR,
					wOut - mother->msbC - weight - wIn - rescaleShift + mother->g,
					GMP_RNDN
				); //Exact

			// Here is when we do the rounding
			mpfr_get_z(result.get_mpz_t(), mpR, GMP_RNDN); // Should be exact

				//Get the real sign
			if(negativeInput != mother->negativeConstant)
			{
				result = (mpz_class(1) << wOut) - result;
			}
		}

#else
		mpz_class result;
		//If the table contains just two values
		//To test
		if(wIn < 2)
		{
			// Now we want to compute the product correctly rounded to 
			// lsbOut-g but we have to coerce MPFR into rounding to this
			// fixed-point format.

			mpfr_t mpR;
			
			mpfr_init2(mpR, 10*wOut);
			// do the mult in large precision
			if(x0 == 0)
				mpfr_set_d(mpR, 0, GMP_RNDN);
			else
			{
				mpfr_set(mpR, mother->absC, GMP_RNDN);
			}

			// Result is integer*C, which is more or less what we need: just
			// scale to add g bits.
			mpfr_mul_2si(mpR, mpR, mother->wOut - mother->wIn - mother->msbC +
					mother->g, GMP_RNDN); //Exact

			// Here is when we do the rounding
			mpfr_get_z(result.get_mpz_t(), mpR, GMP_RNDN); // Should be exact
		}
		else
		{
			int x;
			// get rid of two's complement
			x = x0;
			if(signedInput)
			{
				if ( x0 > ((1<<(wIn-1))-1) )
					x = (1<<wIn) - x;
			}
			// Cast x to mpfr 
			mpfr_t mpX;
			
			mpfr_init2(mpX, wIn);	
			
			mpfr_set_si(mpX, x, GMP_RNDN); // should be exact
			mpfr_mul_2si(mpX, mpX, weight, GMP_RNDN); //Exact
			// now mpX is the integer radix-LUTinput digit, with its proper
			// weight 

			// Now we want to compute the product correctly rounded to LSB
			// lsbOut-g but we have to coerce MPFR into rounding to this
			// fixed-point format.
			mpfr_t mpR;
			
			mpfr_init2(mpR, 10*wOut);	
			// do the mult in large precision
			mpfr_mul(mpR, mpX, mother->absC, GMP_RNDN);

			// Result is integer*C, which is more or less what we need: just
			// scale to add g bits.
			mpfr_mul_2si(
					mpR, 
					mpR, 
					mother->wOut - mother->wIn - mother->msbC + mother->g,
					GMP_RNDN
				); //Exact

			// and add the half-ulp of the result that turns truncation into
			// rounding if g=0 (meaning either one table, or two tables+faithful
			// rounding), we don't need to add this bit
			if(last && mother->g!=0) 
				mpfr_add_si(mpR, mpR, 1<<(mother->g-1), GMP_RNDN);

			// Here is when we do the rounding
			mpfr_get_z(result.get_mpz_t(), mpR, GMP_RNDN); // Should be exact

						
		}
			//Change number sign only if input XOR cste is negative
			if(
						(
					 			(signedInput && (x0 > (1<<(wIn-1))-1)) 
							!=
								mother->negativeConstant
						)
					&& 
					 	result != 0
				) 
			{
					result = (mpz_class(1)<<wOut) - result;
			}
			return  result;
#endif
			return result;
	}
}



