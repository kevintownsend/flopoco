
/*
   A multiply-and-add in a single bit heap

Author:  Florent de Dinechin, Matei Istoan

This file is part of the FloPoCo project
developed by the Arenaire team at Ecole Normale Superieure de Lyon

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
2012-2013.
All rights reserved.
*/



#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "FixMultAdd.hpp"
#include "IntMultiplier.hpp"

using namespace std;

namespace flopoco {


	// The constructor for a stand-alone operator
	FixMultAddBitheap::FixMultAddBitheap(Target* target, Signal* x_, Signal* y_, Signal* a_, int outMSB_, int outLSB_,
												 float ratio_, bool enableSuperTiles_, map<string, double> inputDelays_):
		Operator ( target, inputDelays_ ),
		x(x_), y(y_), a(a_),
		wOut(outMSB_- outLSB_ + 1),
		outMSB(outMSB_),
		outLSB(outLSB_),
		ratio(ratio_),
		enableSuperTiles(enableSuperTiles_)
{

		srcFileName="FixMultAdd";
		setCopyrightString ( "Florent de Dinechin, 2012-2014" );

		// Set up the VHDL library style
		useNumericStd();

		wX = x->MSB() - x->LSB() +1;
		wY = y->MSB() - y->LSB() +1;
		wA = a->MSB() - a->LSB() +1;

		signedIO = (x->isSigned() && y->isSigned());
		// TODO: manage the case when one is signed and not the other.
		if((x->isSigned() && !y->isSigned()) || (!x->isSigned() && y->isSigned()))
		{
			THROWERROR("One operator signed and the other unsigned is currently not supported.")
		}

		// Set the operator name
		{
			ostringstream name;
			name <<"FixMultAdd_";
			name << wX << "x" << wY << "p" << wA << "r" << wOut << "" << (signedIO?"signed":"unsigned");
			name << Operator::getNewUId();
			setName(name.str());
			REPORT(DEBUG, "Building " << name.str());
		}

		// Set up the IO signals
		xname="X";
		yname="Y";
		aname="A";
		rname="R";

		//create the inputs and the outputs of the operator
		addInput (xname,  wX);
		addInput (yname,  wY);
		addInput (aname,  wA);
		addOutput(rname,  wOut, possibleOutputs);

		// Determine the msb and lsb of the full product X*Y
		pMSB = x->MSB() + y->MSB() + 1;
		pLSB = x->LSB() + y->LSB();

		// Determine the actual msb and lsb of the product,
		// from the output's msb and lsb, and the (possible) number of guard bits

		//lsb
		if(pLSB < outLSB)
		{
			//the result of the multiplication will be truncated
			workPLSB = outLSB;
			g = IntMultiplier::neededGuardBits(wX, wY, pMSB-pLSB+1-(outLSB-pLSB));
			workPLSB += g;

			possibleOutputs = 2;
			REPORT(DETAILED, "Faithfully rounded architecture");

			ostringstream dbgMsg;
			dbgMsg << "Multiplication of " << wX << " by " << wY << " bits, with the result truncated to weight " << workPLSB;
			REPORT(DETAILED, dbgMsg.str());
		}
		else
		{
			//the result of the multiplication will be truncated
			workPLSB = outLSB;
			g = 0;

			possibleOutputs = 1; // No faithful rounding
			REPORT(DETAILED, "Exact architecture");

			ostringstream dbgMsg;
			dbgMsg << "Full multiplication of " << wX << " by " << wY;
			REPORT(DETAILED, dbgMsg.str());
		}
		//msb
		if(pMSB < outMSB)
		{
			//all msbs of the product are used
			workPMSB = pMSB;
		}
		else
		{
			//not all msbs of the product are used
			workPMSB = outMSB;
		}

		// Determine the actual msb and lsb of the addend,
		// from the output's msb and lsb

		//lsb
		if(a->LSB() < outLSB)
		{
			//truncate the addend
			workALSB = (outLSB+g < a->LSB()) ? a->LSB() : outLSB+g;
		}
		else
		{
			//keep the full addend
			workALSB = a->LSB();
		}
		//msb
		if(a->MSB() < outMSB)
		{
			//all msbs of the addend are used
			workAMSB = a->MSB();
		}
		else
		{
			//not all msbs of the product are used
			workAMSB = outMSB;
		}

		//create the bit heap
		{
			ostringstream dbgMsg;
			dbgMsg << "Using " << g << " guard bits" << endl;
			dbgMsg << "Creating bit heap of size " << wOut+g << ", out of which " << g << " guard bits";
			REPORT(DETAILED, dbgMsg.str());
		}
		bitHeap = new BitHeap(this,								//parent operator
							  wOut+g,							//size of the bit heap
							  enableSuperTiles);				//whether super-tiles are used
		bitHeap->setSignedIO(signedIO);

		//create the multiplier
		//	this is a virtual operator, which uses the bit heap to do all its computations
		mult = new IntMultiplier(this,							//parent operator
								 bitHeap,						//the bit heap that performs the compression
								 getSignalByName(xname),		//first input to the multiplier (a signal)
								 getSignalByName(yname),		//second input to the multiplier (a signal)
								 wX,							//width of the first operator
								 wY,							//width of the second operator
								 pMSB-pLSB+1-(outLSB-pLSB),		//width of the result
								 g,								//offset in the bit heap (info provided by neededGuardBits() method)
								 false /*negate*/,				//whether to subtract the result of the multiplication from the bit heap
								 signedIO,						//signed/unsigned operator
								 ratio);						//DSP ratio

		//add the addend to the bit heap
		int addendWeight;

		//if the addend is truncated, no shift is needed
		//	else, compute the shift from the output lsb
		addendWeight = (a->LSB() < outLSB) ? 0 : outLSB - a->LSB();
		if(signedIO)
		{
			bitHeap->addSignedBitVector(addendWeight,			//weight of signal in the bit heap
										aname,					//name of the signal
										workAMSB-workALSB+1,	//size of the signal added
										a->LSB()-workALSB);		//index of the lsb in the bit vector from which to add the bits of the addend
		}
		else
		{
			bitHeap->addUnsignedBitVector(addendWeight,			//weight of signal in the bit heap
										  aname,				//name of the signal
										  workAMSB-workALSB+1,	//size of the signal added
										  a->LSB()-workAMSB,	//index of the msb in the actual bit vector from which to add the bits of the addend
										  a->LSB()-workALSB);	//index of the lsb in the actual bit vector from which to add the bits of the addend
		}

		//compress the bit heap
		bitHeap -> generateCompressorVHDL();

		vhdl << tab << rname << " <= " << bitHeap->getSumName() << range(wOut+g-1, g) << ";" << endl;
	}




	FixMultAddBitheap::~FixMultAddBitheap() {
		if(mult)
			free(mult);
	}


	FixMultAddBitheap* FixMultAddBitheap::newComponentAndInstance(Operator* op,
																	string instanceName,
																	string xSignalName,
																	string ySignalName,
																	string aSignalName,
																	string rSignalName,
																	int rMSB,
																	int rLSB
																)
	{
		FixMultAdd* f = new FixMultAdd(op->getTarget(),
										op->getSignalByName(xSignalName),
										op->getSignalByName(ySignalName),
										op->getSignalByName(aSignalName),
										rMSB, rLSB);
		op->addSubComponent(f);
		op->inPortMap(f, "X", xSignalName);
		op->inPortMap(f, "Y", ySignalName);
		op->inPortMap(f, "A", aSignalName);

		op->outPortMap(f, "R", join(rSignalName, "_slv"));

		op->vhdl << op->instance(f, instanceName);
		op->vhdl << tab << op->declareFixPoint(rSignalName,f->signedIO, rMSB, rLSB) << " <= " <<  "signed(" << (join(rSignalName, "_slv")) << ");" << endl;
	}


	//FIXME: is this right? emulate function needs to be checked
	//		 it makes no assumptions on the signals (only on their relative alignments),
	//		 so it should work for all combinations of the inputs
	void FixMultAddBitheap::emulate ( TestCase* tc )
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svA = tc->getInputValue("A");

		mpz_class svP, svR, svRAux;
		mpz_class twoToWR = (mpz_class(1) << (wOut));
		mpz_class twoToWRm1 = (mpz_class(1) << (wOut-1));

		if(!signedIO)
		{
			int outShift;

			svP = svX * svY;

			//align the product and the addend
			if(a->LSB() < pLSB)
			{
				svP = svP << pLSB-a->LSB();
				outShift = outLSB - a->LSB();
			}
			else
			{
				svA = svA >> a->LSB()-pLSB;
				outShift = outLSB - pLSB;
			}

			svR = svP + svA;

			//align the multiply-and-add with the output format
			if(outShift > 0)
				svR = svR >> outShift;
			else
				svR = svR << (-outShift);

			//add only the bits corresponding to the output format
			svRAux = svR & (twoToWR -1);

			tc->addExpectedOutput("R", svRAux);
			if(possibleOutputs==2)
			{
				svR++;
				svR &= (twoToWR -1);
				tc->addExpectedOutput("R", svR);
			}
		}
		else
		{
			//TODO

			// Manage signed digits
			mpz_class twoToWX = (mpz_class(1) << (wX));
			mpz_class twoToWXm1 = (mpz_class(1) << (wX-1));
			mpz_class twoToWY = (mpz_class(1) << (wY));
			mpz_class twoToWYm1 = (mpz_class(1) << (wY-1));
			mpz_class twoToWA = (mpz_class(1) << (wA));
			mpz_class twoToWAm1 = (mpz_class(1) << (wA-1));

			if (svX >= twoToWXm1)
				svX -= twoToWX;

			if (svY >= twoToWYm1)
				svY -= twoToWY;

			if (svA >= twoToWAm1)
				svA -= twoToWA;

			svP = svX * svY; //signed
			if(lsbPfull>=0)
			{
				// no truncation
				svR = (svA << lsbA) + (svP << lsbPfull); // signed
				// manage two's complement at output
				if ( svR < 0)
					svR += twoToWR;
				tc->addExpectedOutput("R", svR);
				return;
			}
			else
			{
				//lsbPfull<0
				int shift=-lsbPfull;
				// fully accurate result, product-anchored
				svR = (svA<<(lsbA + shift)) + svP;
				// manage its two's complement
				if ( svR < 0)
					svR += (mpz_class(1) << (wOut+shift));
				// shift back to place: truncation
				svR = svR >> shift;
				tc->addExpectedOutput("R", svR); // this was rounded down
				svR++;
				svR &= (twoToWR -1);
				tc->addExpectedOutput("R", svR);
			}
		}
	}



	void FixMultAddBitheap::buildStandardTestCases(TestCaseList* tcl)
	{
		//TODO
	}




}
