/*
 * A faithful multiplier by a real constant, using a variation of the KCM method
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 * 
 * Authors : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr
 * 			 3IF-Dev-Team-2015
 *
 * Initial software.
 * Copyright © ENS-Lyon, INRIA, CNRS, UCBL, 
 * 2008-2011.
 * All rights reserved.
 */

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

#define WIP_FORGET
#ifdef WIP_LFORGET
#pragma message("Version du FixRealKCM en cours de développement")
#pragma message("Statut : Fonctionne globalement bien mais n'est pas encore optimisé")
#endif

namespace flopoco{

	/**
	* @brief init : all operator initialisation stuff goes here
	*/
	void FixRealKCM::init()
	{
		srcFileName="FixRealKCM";

		setCopyrightString("Florent de Dinechin (2007-2011-?), 3IF Dev Team 2015");

		if(lsbIn>msbIn) 
			throw string("FixRealKCM: Error, lsbIn>msbIn");
    
		if(targetUlpError > 1.0)
			THROWERROR("FixRealKCM: Error, targetUlpError="<<
					targetUlpError<<">1.0. Should be in ]0.5 ; 1].");
		//Error on final rounding is er <= 2^{lsbout - 1} = 1/2 ulp so 
		//it's impossible to reach 0.5 ulp of precision if cste is real
		if(targetUlpError <= 0.5) 
			THROWERROR("FixRealKCM: Error, targetUlpError="<<
					targetUlpError<<"<0.5. Should be in ]0.5 ; 1]");
		
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
		if(!signedInput && negativeConstant)
		{
			msbOut++; //Result would be signed
		}

		wOut = msbOut - lsbOut +1;

		REPORT(DEBUG, "msbConstant=" << msbC  << "   (msbIn,lsbIn)=("<< 
				msbIn << "," << lsbIn << ")   wIn=" << wIn << 
				"   (msbOut,lsbOut)=(" << msbOut << "," << lsbOut <<
				")   wOut=" << wOut
			);

		g = neededGuardBits(
				getTarget(), 
				wIn, 
				targetUlpError, 
				constant, 
				lsbIn, 
				lsbOut
			);
	}
	
	int FixRealKCM::guardBitsFromTableNumber(
			int nbOfTables,
			double targetUlpError
		)
	{
		int guardBits;
		if(targetUlpError > 1.0 || targetUlpError <= 0.5)
		{
			cerr << "WARNING : Target ulp error should be in ]0.5 ; 1]. " <<
					"Value provided : " << targetUlpError << 
					"Will be considered as 1.0. Please, considere this"
					" warning as a bug.";
			targetUlpError = 1.0; 
		}
				
		if((nbOfTables <= 2 && targetUlpError==1.0) || nbOfTables == 1)
		{
			// specific case: two CR table make up a faithful sum
			guardBits = 0; 
		}
		else
		{
			guardBits = ceil(log2((double)nbOfTables/(2*targetUlpError - 1)));
		}
			
		return guardBits;
	}

	int FixRealKCM::computeTableNumbers(
			Target* target,
			int wIn,
			int msbC,
			int lsbIn,
			int lsbOut,
			double targetUlpError,
			int** disize_target
		)
	{
		/** will be target->lutInputs() or target->lutInputs()-1  */
		int optimalTableInputWidth = target->lutInputs()-1;
		int* diSize = new int[17*42];		
		int nbOfTables, guardBits;
		int newWIn = wIn;
		int newLsbIn = lsbIn;

		//The loop is here to prevent neglictible input bits from being
		//tabulated.
		do
		{ 
			wIn = newWIn;
			lsbIn = newLsbIn;

			int offset = 0;
			int nbTablesEntieres = wIn / optimalTableInputWidth;
			int remainingBits = wIn % optimalTableInputWidth;
			nbOfTables = nbTablesEntieres;

			if(remainingBits != 0)
			{ 
				//On each case we will need to handle first table in width size
				//separately
				offset++;

				int guardBits_extendedTable = 
					guardBitsFromTableNumber(nbTablesEntieres, targetUlpError);
				int guardBits_extraTable =
					guardBitsFromTableNumber(nbTablesEntieres + 1, targetUlpError);
				//Cost for extended table is nb  of extra tables * output width
				int lutCost_extendedTable = ((1 << remainingBits) - 1) * (
						optimalTableInputWidth + remainingBits -lsbOut);
				//TODO measure bitheap impact 
				//(here cost is only extraguardbits * nbOfTables) 
				int lutCost_extraTable = 
					(guardBits_extendedTable - guardBits_extraTable) * 
					(nbTablesEntieres + 1);
				if(lutCost_extraTable < lutCost_extendedTable)
				{
					nbOfTables++;
					diSize[0] = remainingBits;
				}
				else
				{
					diSize[0] = remainingBits + optimalTableInputWidth;
				}
			}
			cout << "Nb tables : " << nbOfTables << endl;
			for(int i = offset ; i < nbOfTables ; diSize[i++] = optimalTableInputWidth);
			

			guardBits = guardBitsFromTableNumber(nbOfTables, targetUlpError);
			newLsbIn = lsbOut - guardBits - msbC;
			newWIn = wIn - (newLsbIn - lsbIn);
		}while(newLsbIn > lsbIn);

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
	

	void FixRealKCM::connectBitHeap(
				FixRealKCMTable** t,
				int* doSize,
				int nbOfTables,
				Operator* op
			)
	{
		Target* target = getTarget();

		for(int i = 0; i < nbOfTables; i++)
		{
			REPORT(DEBUG, "Adding bits for table " << i)
				//manage the critical path
			op->setCycleFromSignal(join("d", i, "_kcmMult_", getuid()));
			op->manageCriticalPath(target->lutDelay());

			op->inPortMap (t[i] , "X", join("d",i, "_kcmMult_", getuid()));
			op->outPortMap(t[i] , "Y", join("pp",i, "_kcmMult_", getuid()));
			op->vhdl << op->instance(t[i] , join("KCMTable_",i, "_kcmMult_", getuid()));

			//manage the critical path
			op->syncCycleFromSignal(join("pp",i,"_kcmMult_", getuid()));

			int bitheapSize = bitHeap->getMaxWeight()-bitHeap->getMinWeight();
			//add the bits to the bit heap
			int w;
			for(w=0; w < doSize[i]; w++)
			{
				stringstream s;

				//manage the critical path
				manageCriticalPath(target->lutDelay());

				s << join("pp",i, "_kcmMult_", getuid()) << of(w);

				bitHeap->addBit(w, s.str());
			} // w = table.msb + 1
	
			//Negative subproduct sign extension :
			//As fast sign extension was enabled, we only need to add
			//1 to each weight from table.msb to wOut - 1
			if(t[i]->negativeSubproduct)
			{	
				for(; w < bitheapSize ; w++)
				{
					bitHeap->addConstantOneBit(w);
				}
				bitHeap->addConstantOneBit(0);
			} 
			else if (i == nbOfTables - 1 && signedInput)
			{
				for(w-- ; w < bitheapSize ; w++)
				{
					bitHeap->addConstantOneBit(w);
				}
			}
		}

		if(g > 0)
		{
			REPORT(INFO, "Adding one half ulp to transform truncation into"
					" faithful rounding");
			bitHeap->addConstantOneBit(g-1);
		}

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
		int nbOfTables = computeTableNumbers(
				target, 
				wIn, 
				msbC, 
				lsbIn, 
				lsbOut, 
				targetUlpError, 
				&diSize
			);
		
		REPORT(INFO, "Constant multiplication in "<< nbOfTables << " tables." <<
				g << "guards bits are used.");
		//manage the critical path
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		addInput("X", wIn);
		addOutput("R", wOut);

		int* doSize;

		FixRealKCMTable** t = createTables(
				target,
				diSize,
				nbOfTables,
				&doSize,
				this
			);


		setCriticalPath(getMaxInputDelays(inputDelays));

		//create the bitheap
		bitHeap = new BitHeap(this, wOut+g);
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		connectBitHeap(t, doSize, nbOfTables, this);
		vhdl << tab << "R <= OutRes" << range(wOut+g-1, g) << ";" << endl;
		outDelayMap["R"] = getCriticalPath();

		//compress the bitheap and produce the result
		bitHeap->generateCompressorVHDL();

		//manage the critical path
		syncCycleFromSignal(bitHeap->getSumName());

		//because of final add in bit heap, add one more bit to the result
		vhdl << declare("OutRes", wOut+g) << " <= " << 
			bitHeap->getSumName() << range(wOut+g-1, 0) << ";" << endl;	

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

		// First set up all the sizes
		int *diSize;
		int nbOfTables = computeTableNumbers(
				target, 
				wIn, 
				msbC, 
				lsbIn, 
				lsbOut, 
				targetUlpError, 
				&diSize
			);

		REPORT(INFO, "Constant multiplication in "<< nbOfTables << " tables");		
		

		int* doSize;

		FixRealKCMTable** t = createTables(
				target,
				diSize,
				nbOfTables,
				&doSize,
				parentOp,
				multiplicandX->getName()
			);

		connectBitHeap(t, doSize, nbOfTables, parentOp);
		vhdl << tab << "R <= OutRes" << range(wOut+g-1, g) << ";" << endl;
		outDelayMap["R"] = getCriticalPath();

		delete[] diSize;
	}

	FixRealKCMTable** FixRealKCM::createTables(
			Target* target,
			int* diSize,
			int nbOfTables,
			int** doSize_target,
			Operator* op,
			string inputSignalName
		)
	{
		FixRealKCMTable** t = new FixRealKCMTable*[nbOfTables]; 
		int* doSize = new int[17*42];

		//Useful bit entry width (we doesn't care of too precise bits)
		int effectiveWin;
		int i;
		for(effectiveWin = i = 0; i < nbOfTables ; effectiveWin += diSize[i++]);

		//first split the input X into digits having lutWidth bits -> this
		//is as generic as it gets :)
		bool last = true;
		//Will be current table lsb weight
		int highBit = msbIn + 1;
		int offset = wIn;
		//out width of current table
		int tableDo = wOut+g;

		for (int i=nbOfTables-1; i>=0; i--)
		{
			// The last table has to have wOut+g  bits.
			// The previous one wOut+g-lastLutWidth
			// the previous one wOut+g-anteLastLutWidth-lastlutWidth etc
			
			op-> vhdl << tab << 
				op->declare(join("d",i,"_kcmMult_", getuid()), diSize[i]) << 
				" <= " << 
				inputSignalName << range(offset - 1, offset - diSize[i]) <<  
				";" << endl;


			highBit -= diSize[i];
			offset -= diSize[i];
			
			// already updated 
			t[i] = new FixRealKCMTable(
					target, 
					this, 
					i, 
					highBit, 
					diSize[i], 
					tableDo, 
					negativeConstant && !last, 
					last 
				);

			doSize[i] = tableDo;
			REPORT(DEBUG, "Table i=" << i << ", input size=" << 
					diSize[i] << ", output size=" << doSize[i] << 
					", weight= " << highBit);

			op->useSoftRAM(t[i]);
			op->addSubComponent(t[i]);

			tableDo -= diSize[i];
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
				svX = (mpz_class(1)<<wIn) - svX;
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
			double targetUlpError,
			string constant,
			int lsbIn, 
			int lsbOut
		)
	{
		// Convert the input string into a sollya evaluation tree
		sollya_obj_t node;
		node = sollya_lib_parse_string(constant.c_str());	
		/* If  parse error throw an exception */
		if (sollya_lib_obj_is_error(node))
		{
				ostringstream error;
				error << " neededGuardBits : Unable to parse string "<< 
					constant << " as a numeric constant" <<endl;
				throw error.str();
		}

		mpfr_t mpC, absC;

		mpfr_init2(mpC, 10000);
		mpfr_init2(absC, 10000);
		sollya_lib_get_constant(mpC, node);

		mpfr_abs(absC, mpC, GMP_RNDN);

		mpfr_t log2C;
		mpfr_init2(log2C, 100); // should be enough for anybody
		mpfr_log2(log2C, absC, GMP_RNDN);
		int msbC = mpfr_get_si(log2C, GMP_RNDU);
		mpfr_clears(log2C, absC, mpC, NULL);

		int nbOfTables = computeTableNumbers(
				target, 
				wIn, 
				msbC, 
				lsbIn, 
				lsbOut,
				targetUlpError, 
				nullptr
			);
		return guardBitsFromTableNumber(nbOfTables, targetUlpError);
	}



	/************************** The FixRealKCMTable class ********************/


	FixRealKCMTable::FixRealKCMTable(
			Target* target, 
			FixRealKCM* mother, 
			int i, 
			int weight, 
			int wIn, 
			int wOut, 
			bool negativeSubproduct, 
			bool last, 
			int pipeline
		):
			Table(target, wIn, wOut, 0, -1, pipeline), 
			mother(mother), 
			index(i), 
			weight(weight), 
			negativeSubproduct(negativeSubproduct), 
			last(last)
	{
		ostringstream name; 
		srcFileName="FixRealKCM";
		name << mother->getName() << "_Table_" << index;
		setName(name.str());
		setCopyrightString("Florent de Dinechin (2007-2011-?), 3IF Dev Team"); 
	}
  
	mpz_class FixRealKCMTable::function(int x0)
	{
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
		mpfr_init2(mpX, 2*wIn); //To avoid mpfr bug if wIn = 1

		if(x == 0)
		{
			result = mpz_class(0);
		}
		else
		{
			mpfr_set_si(mpX, x, GMP_RNDN); // should be exact
			//Getting real weight
			mpfr_mul_2si(mpX, mpX, weight, GMP_RNDN); //Exact


			// do the mult in large precision
			mpfr_mul(mpR, mpX, mother->absC, GMP_RNDN);
			
			// Result is integer*C, which is more or less what we need: just
			// scale to add g bits.
			mpfr_mul_2si(
					mpR, 
					mpR,
					-mother->lsbOut + mother->g,	
					GMP_RNDN
				); //Exact

			// Here is when we do the rounding
			mpfr_get_z(result.get_mpz_t(), mpR, GMP_RNDN); // Should be exact

			//Get the real sign
			if(negativeInput != mother->negativeConstant && result != 0)
			{
				result = (mpz_class(1) << wOut) - result;
			}
		}

		//Result is result -1 in order to avoid to waste a sign bit
		if(negativeSubproduct)
		{
			if(result == 0)
				result = (mpz_class(1) << wOut) - 1;
			else
				result -= 1;
		}

		//In case of global bitHeap, we need to invert msb for signedInput
		//last bit for fast sign extension.
		if(last && mother->signedInput)
		{
			mpz_class shiffter = mpz_class(1) << (wOut - 1);
			if(result < shiffter )
			{
				result += shiffter;
			}
			else
			{
				result -= shiffter;
			}
		}


		return result;
	}
}



