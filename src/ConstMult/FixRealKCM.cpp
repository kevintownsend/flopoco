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

namespace flopoco{

	/**
	* @brief init : all operator initialisation stuff goes here
	*/
	void FixRealKCM::init()
	{
		useNumericStd();

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


	bool FixRealKCM::handleSpecialConstant(Operator* op, string inputName)
	{
		if(mpfr_cmp_ui_2exp(absC, 1, msbC) == 0)
		{
			REPORT(DEBUG, "Constant is a power of two. " 
					"Simple shift will be used instead of tables");

			int shiffterWidth = msbC + msbIn - lsbOut + 1;
			int copyLength = min(wIn, shiffterWidth) - 1;
			int shifStart = shiffterWidth - 1;
			int shifStop = shifStart - copyLength;

			op->setCycleFromSignal(inputName);
			op->vhdl << tab << 
				op->declare(join("kcm_", getuid(), "_shiffter"), shiffterWidth)<<	
				range(shifStart, shifStop) << " <= " << 
				inputName << range(wIn - 1, wIn - 1 - copyLength) << ";" << endl;
			if(wIn < shiffterWidth)
			{
				op->vhdl << tab << join("kcm_", getuid(), "_shiffter") << 
					range (shiffterWidth - wIn - 1 , 0) << " <= " << 
					zg(shiffterWidth - wIn, 0) << ";" << endl;
			}

			if(signedInput)
			{
				if(negativeConstant)
				{
					bitHeap->subtractSignedBitVector(
							lsbOut - bitheaplsb,
							join("kcm_", getuid(), "_shiffter"),
							shiffterWidth
						);
				}
				else
				{
					bitHeap->addSignedBitVector(
							lsbOut - bitheaplsb,
							join("kcm_", getuid(), "_shiffter"),
							shiffterWidth	
						);
				}
			}
			else
			{
				if(negativeConstant)
				{
					bitHeap->subtractUnsignedBitVector(
							lsbOut - bitheaplsb,
							join("kcm_", getuid(), "_shiffter"),
							shiffterWidth	
						);
				}
				else
				{
					bitHeap->addUnsignedBitVector(
							lsbOut - bitheaplsb,
							join("kcm_", getuid(), "_shiffter"),
							shiffterWidth	
						);
				}
			}
			return true;
		}
		return false;
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

	int FixRealKCM::computeTablesNumber(
			Target* target,
			int wIn,
			int msbC,
			int lsbIn,
			int lsbOut,
			double targetUlpError,
			int** disize_target
		)
	{
		int oldWIn = wIn;
		/** will be target->lutInputs() or target->lutInputs()-1  */
		int optimalTableInputWidth = target->lutInputs()-1;
		int* diSize = nullptr;		
		int nbOfTables, guardBits;
		int wOut = wIn + lsbIn + msbC - lsbOut;
		int newWIn = wOut + msbC + 1;
		int newGuardBits = 0;

		if(wIn <= newWIn)
		{
			int nbTablesEntieres = wIn / optimalTableInputWidth;
			int remainingBits = wIn % optimalTableInputWidth;
			nbOfTables = nbTablesEntieres;

			if(remainingBits != 0)
			{ 
				int guardBits_extendedTable = 
					guardBitsFromTableNumber(nbTablesEntieres, targetUlpError);
				int guardBits_extraTable =
					guardBitsFromTableNumber(nbTablesEntieres + 1, targetUlpError);

				//TODO : compute more accurately costs and compare them
				if	(	guardBits_extraTable == guardBits_extendedTable )
				{
					nbOfTables++;
					diSize = new int[nbOfTables];
					for(int i = 0 ; i + 1 < nbOfTables ; diSize[i++] = optimalTableInputWidth );
					diSize[nbOfTables - 1] = remainingBits;
				}
				else
				{
					diSize = new int[nbOfTables];
					diSize[0] = remainingBits + optimalTableInputWidth;
					for(int i = 1 ; i < nbOfTables ; diSize[i++] = optimalTableInputWidth);
				}
			}
			else
			{
				diSize = new int[nbOfTables];	
				for(int i = 0 ; i < nbOfTables ; diSize[i++] = optimalTableInputWidth);
			}

		}
		else {
			cerr << "Input precision higher than required. Trying to optimize" << endl;
			//The loop is here to prevent neglictible input bits from being
			//tabulated.
			do { 
				if(diSize != nullptr) {
					delete diSize;
					diSize = nullptr;
				}
				wIn = newWIn;
				guardBits = newGuardBits;

				int nbTablesEntieres = wIn / optimalTableInputWidth;
				int remainingBits = wIn % optimalTableInputWidth;
				nbOfTables = nbTablesEntieres;

				if (remainingBits != 0) { 
					int guardBits_extendedTable = 
						guardBitsFromTableNumber(nbTablesEntieres, targetUlpError);
					int guardBits_extraTable =
						guardBitsFromTableNumber(nbTablesEntieres + 1, targetUlpError);
					
					//TODO : compute more accurately costs and compare them
					if	(	guardBits_extraTable == guardBits_extendedTable ) {
						nbOfTables++;
						diSize = new int[nbOfTables];
						for(int i = 0 ; i + 1 < nbOfTables ; diSize[i++] = optimalTableInputWidth);
						diSize[nbOfTables - 1] = remainingBits;
					} else {
						diSize = new int[nbOfTables];
						diSize[0] = remainingBits + optimalTableInputWidth;
						for(int i = 1 ; i < nbOfTables ; diSize[i++] = optimalTableInputWidth );
					}
				} else {
					diSize = new int[nbOfTables];
					for(int i = 0 ; i < nbOfTables ; diSize[i++] = optimalTableInputWidth);
				}
				
				newGuardBits = guardBitsFromTableNumber(nbOfTables, targetUlpError);
				newWIn = wIn + newGuardBits - guardBits;
			}while(newWIn > wIn && wIn <= oldWIn && newWIn <= oldWIn);
		}

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
	

	void FixRealKCM::connectTablesToBitHeap(
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
			int offset = lsbOut-g-bitheaplsb;
			int w;
			for(w=0; w < doSize[i]; w++)
			{
				stringstream s;

				//manage the critical path
				manageCriticalPath(target->lutDelay());

				s << join("pp",i, "_kcmMult_", getuid()) << of(w);
				bitHeap->addBit(w+offset, s.str());

			} // w = table.msb + 1
	
			//Negative subproduct sign extension :
			//As fast sign extension was enabled, we only need to add
			//1 to each weight from table.msb to wOut - 1
			if(t[i]->negativeSubproduct)
			{	
				REPORT(DEBUG, "Negative subproduct sign extension for table "<<i);
				for(; w < bitheapSize ; w++)
				{
					bitHeap->addConstantOneBit(w);
				}
				bitHeap->addConstantOneBit(0);
			} 
			else if (i == nbOfTables - 1 && signedInput)
			{
				REPORT(DEBUG, "Sign extension for signed input msb table");
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
			bitHeap->addConstantOneBit(lsbOut-bitheaplsb-1);
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
		bitheaplsb = lsbOut - g;
		addInput("X", wIn);
		addOutput("R", wOut);
	
		//Zero constant
		if(mpfr_zero_p(mpC) != 0)
		{
			vhdl << tab << "R" << range(wOut - 1, 0) << " <= " << zg(wOut, 0) <<
				";" << endl;
		} else { //NonZero constant

			//create the bitheap
			bitHeap = new BitHeap(this, wOut+g);

			if(!handleSpecialConstant(this))
			{
				int* diSize;
				int nbOfTables = computeTablesNumber(
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

				int* doSize;

				setCriticalPath(getMaxInputDelays(inputDelays));
				manageCriticalPath(target->localWireDelay() + target->lutDelay());

				FixRealKCMTable** t = createTables(
						diSize,
						nbOfTables,
						&doSize,
						this
						);

				connectTablesToBitHeap(t, doSize, nbOfTables, this);

				delete[] diSize;
			}

			//compress the bitheap and produce the result
			bitHeap->generateCompressorVHDL();

			//manage the critical path
			syncCycleFromSignal(bitHeap->getSumName());

			//because of final add in bit heap, add one more bit to the result
			vhdl << declare("OutRes", wOut+g) << " <= " << 
				bitHeap->getSumName() << range(wOut+g-1, 0) << ";" << endl;	

			vhdl << tab << "R <= OutRes" << range(wOut+g-1, g) << ";" << endl;

		}

				outDelayMap["R"] = getCriticalPath();
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
			int bitheapLsb,
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
			bitHeap(bitHeap_),
			bitheaplsb(bitheapLsb)
	{
		
		init();

		//If constant is not zero or a power of two, then create tables and
		//branc them to bitHeap
		if	(	mpfr_zero_p(mpC) == 0 && 
				!handleSpecialConstant(parentOp, multiplicandX->getName())
			)
		{
			// First set up all the sizes
			int *diSize;
			int nbOfTables = computeTablesNumber(
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
					diSize,
					nbOfTables,
					&doSize,
					parentOp,
					multiplicandX->getName()
					);

			connectTablesToBitHeap(t, doSize, nbOfTables, parentOp);

			delete[] diSize;
		}

	}

	FixRealKCMTable** FixRealKCM::createTables(
			int* diSize,
			int nbOfTables,
			int** doSize_target,
			Operator* op,
			string inputSignalName
		)
	{
		Target* target = getTarget();
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
					negativeConstant && (!last || !signedInput), 
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

		//Special cases : constant is 0 or a power of two
		if(mpfr_cmp_ui_2exp(absC, 1, msbC) == 0 || mpfr_zero_p(mpC) != 0)
		{
			mpfr_clears(log2C, absC, mpC, NULL);
			return 0;
		}
		mpfr_clears(log2C, absC, mpC, NULL);

		int nbOfTables = computeTablesNumber(
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

	OperatorPtr FixRealKCM::parseArguments(
			Target* target, 
			std::vector<std::string> &args
		)
	{
		int lsbIn, lsbOut, msbIn;
		bool signedInput;
		double targetUlpError;
		string constant;
		UserInterface::parseInt(args, "lsbIn", &lsbIn);
		UserInterface::parseString(args, "constant", &constant);
		UserInterface::parseInt(args, "lsbOut", &lsbOut);
		UserInterface::parseInt(args, "msbIn", &msbIn);
		UserInterface::parseBoolean(args, "signedInput", &signedInput);
		UserInterface::parseFloat(args, "targetUlpError", &targetUlpError);	
		return new FixRealKCM(
				target, 
				signedInput,
				msbIn,
				lsbIn,
				lsbOut,
				constant, 
				targetUlpError
			);
	}

	void FixRealKCM::registerFactory()
	{
		UserInterface::add(
				"FixRealKCM",
				"Table based real multiplier. Output size is computed",
				"ConstMultDiv",
				"",
				"signedInput(bool): 0=unsigned, 1=signed; \
         msbIn(int): weight associated to most significant bit (including sign bit);\
         lsbIn(int): weight associated to least significant bit;\
         lsbOut(int): weight associated to output least significant bit; \
         constant(string): constant given in arbitrary-precision decimal, or as a Sollya expression, e.g \"log(2)\"; \
         targetUlpError(real)=1.0: required precision on last bit. Should be strictly greater than 0.5 and lesser than 1;",
				"This variant of Ken Chapman's Multiplier is briefly described in <a href=\"bib/flopoco.html#DinIstoMas2014-SOPCJR\">this article</a>.<br> Special constants, such as 0 or powers of two, are handled efficiently.",
				FixRealKCM::parseArguments		
			);
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





