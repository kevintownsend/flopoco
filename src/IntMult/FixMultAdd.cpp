
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
	FixMultAdd::FixMultAdd(Target* target, Signal* x_, Signal* y_, Signal* a_, int outMSB_, int outLSB_,
												 float ratio_, bool enableSuperTiles_, map<string, double> inputDelays_):
		Operator ( target, inputDelays_ ),
		x(x_), y(y_), a(a_),
		outMSB(outMSB_),
		outLSB(outLSB_),
		wOut(outMSB_- outLSB_ + 1),
		ratio(ratio_),
		enableSuperTiles(enableSuperTiles_)
{

		srcFileName="FixMultAdd";
		setCopyrightString ( "Matei Istoan, Florent de Dinechin, 2012-2014" );

		// Set up the VHDL library style
		//useNumericStd();

		wX = x->MSB() - x->LSB() +1;
		wY = y->MSB() - y->LSB() +1;
		wA = a->MSB() - a->LSB() +1;

		signedIO = (x->isFixSigned() && y->isFixSigned());
		// TODO: manage the case when one is signed and not the other.
		if((x->isFixSigned() && !y->isFixSigned()) || (!x->isFixSigned() && y->isFixSigned()))
			THROWERROR("One operator signed and the other unsigned is currently not supported.")

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

		// Determine the msb and lsb of the full product X*Y
		pMSB = x->MSB() + y->MSB() + 1;
		pLSB = x->LSB() + y->LSB();

		// Determine the actual msb and lsb of the product,
		// from the output's msb and lsb, and the (possible) number of guard bits


		// workPLSB is the actual weight of the bit that has weight 0 in the bit heap
		int PlsbWeightInBitHeap;
		//lsb
		if(pLSB < outLSB)	{
			possibleOutputs = 2;
			REPORT(DEBUG, "pLSB="<<pLSB<<" < outLSB="<< outLSB<< ", the result of the multiplication will be truncated and the architecture will be faithful");
			// workPLSB = outLSB;
			g = IntMultiplier::neededGuardBits(target, wX, wY, workPMSB-workPLSB+1);
			// workPLSB -= g;
			//			REPORT(DETAILED,  "Multiplication of " << wX << " by " << wY << " bits, with the result truncated to weight " << workPLSB);
		}
		else { // pLSB>=outLSB: the result of the multiplication will not be truncated
			//workPLSB = pLSB;
      g=0;
			//			PlsbWeightInBinHeap=
			possibleOutputs = 1; // No faithful rounding
			REPORT(DETAILED, "Exact architecture: Full multiplication of " << wX << " by " << wY);
		}

		//msb
		if(pMSB <= outMSB)	{	    //all msbs of the product are used
			workPMSB = pMSB;
		}
		else	{                   //not all msbs of the product are used
			workPMSB = outMSB;
		}

		// Determine the actual msb and lsb of the addend, from the output's msb and lsb

		//lsb
		if(a->LSB() < outLSB) {	  //truncate the addend
			workALSB = (outLSB-g < a->LSB()) ? a->LSB() : outLSB-g;
		}
		else	{                   //keep the full addend
			workALSB = a->LSB();
		}
		//msb
		if(a->MSB() <= outMSB)	{	//all msbs of the addend are used
			workAMSB = a->MSB();
		}
		else	{	                  //not all msbs of the addend are used
			workAMSB = outMSB;
		}

		//create the inputs and the outputs of the operator
		addInput (xname,  wX);
		addInput (yname,  wY);
		addInput (aname,  wA);
		addOutput(rname,  wOut, possibleOutputs);

		//create the bit heap
		REPORT(DETAILED, "Creating bit heap of size " << wOut+g << ", out of which " << g << " guard bits");
		
		bitHeap = new BitHeap(this,								//parent operator
													wOut+g,							//size of the bit heap
													enableSuperTiles);	//whether super-tiles are used
		bitHeap->setSignedIO(signedIO);

		//create the multiplier
		//	this is a virtual operator, which uses the bit heap to do all its computations
		mult = new IntMultiplier(this,							//parent operator
														 bitHeap,						//the bit heap that performs the compression
														 getSignalByName(xname),		//first input to the multiplier (a signal)
														 getSignalByName(yname),		//second input to the multiplier (a signal)
														 workPLSB-g,			                  //offset of the LSB of the multiplier in the bit heap
														 false /*negate*/,				//whether to subtract the result of the multiplication from the bit heap
														 signedIO,						//signed/unsigned operator
														 ratio);						//DSP ratio
		
		//add the addend to the bit heap
		int addendWeight;

		//if the addend is truncated, no shift is needed
		//	else, compute the shift from the output lsb
		//the offset for the signal in the bitheap can be either positive (signal's lsb < than bitheap's lsb), or negative
		//	the case of a negative offset is treated with a correction term, as it makes more sense in the context of a bitheap
		//addendWeight = (a->LSB() < outLSB) ? 0 : outLSB - a->LSB();
		addendWeight = a->LSB()-(outLSB-g);

		if(signedIO)		{
			bitHeap->addSignedBitVector(addendWeight,			//weight of signal in the bit heap
																	aname,					//name of the signal
																	workAMSB-workALSB+1,	//size of the signal added
																	workALSB-a->LSB(),		//index of the lsb in the bit vector from which to add the bits of the addend
																	(addendWeight<0));		//if we are correcting the index in the bit vector with a negative weight
		}
		else		{
			bitHeap->addUnsignedBitVector(addendWeight,			//weight of signal in the bit heap
																		aname,				//name of the signal
																		workAMSB-workALSB+1,	//size of the signal added
																		(a->MSB()-a->LSB())-(workAMSB-a->MSB()),
																		//index of the msb in the actual bit vector from which to add the bits of the addend
																		workALSB-a->LSB(),	//index of the lsb in the actual bit vector from which to add the bits of the addend
																		(addendWeight<0));	//if we are correcting the index in the bit vector with a negative weight
		}

		//final rounding, if needed
		if(g >=1 )	{
			bitHeap -> addConstantOneBit(g-1);
		}

		//compress the bit heap
		bitHeap -> generateCompressorVHDL();

		vhdl << tab << rname << " <= " << bitHeap-> getSumName() << range(wOut+g-1, g) << ";" << endl;
}




	FixMultAdd::~FixMultAdd()
	{
		if(mult)
			free(mult);
		if(plotter)
			free(plotter);
	}


	FixMultAdd* FixMultAdd::newComponentAndInstance(Operator* op,
																	string instanceName,
																	string xSignalName,
																	string ySignalName,
																	string aSignalName,
																	string rSignalName,
																	int rMSB,
																	int rLSB
																)
	{
		Signal* x = op->getSignalByName(xSignalName);
		Signal* y = op->getSignalByName(ySignalName);
		Signal* a = op->getSignalByName(aSignalName);
		FixMultAdd* f = new FixMultAdd(op->getTarget(),
																	 x, y, a,
																	 rMSB, rLSB);

		op->vhdl << tab << "-- FixMultAdd: X("<<x->MSB()<<", "<<x->LSB()<<")  Y("<<y->MSB()<<", "<<y->LSB()<<") A("<<a->MSB()<<", "<<a->LSB()<<") R("<<rMSB<<", "<<rLSB<<")" << endl;
		op->addSubComponent(f);
		op->inPortMap(f, "X", xSignalName);
		op->inPortMap(f, "Y", ySignalName);
		op->inPortMap(f, "A", aSignalName);

		op->outPortMap(f, "R", join(rSignalName, "_slv"));

		op->vhdl << op->instance(f, instanceName);
		op->vhdl << tab << op->declareFixPoint(rSignalName,f->signedIO, rMSB, rLSB) << " <= " <<  "signed(" << (join(rSignalName, "_slv")) << ");" << endl;

		return f;
	}


	//FIXME: is this right? emulate function needs to be checked
	//		 it makes no assumptions on the signals, except their relative alignments,
	//		 so it should work for all combinations of the inputs
	void FixMultAdd::emulate ( TestCase* tc )
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
				svP = svP << (pLSB-a->LSB());
				outShift = outLSB - a->LSB();
			}
			else
			{
				svA = svA << (a->LSB()-pLSB);
				outShift = outLSB - pLSB;
			}

			svR = svP + svA;

			//align the multiply-and-add with the output format
			if(outShift > 0)
			{
				svR = svR >> outShift;
				possibleOutputs = 2;
			}
			else
			{
				svR = svR << (-outShift);
				possibleOutputs = 1;
			}

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
			int outShift;

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

			//align the product and the addend
			if(a->LSB() < pLSB)
			{
				svP = svP << (pLSB-a->LSB());
				outShift = outLSB - a->LSB();
			}
			else
			{
				svA = svA >> (a->LSB()-pLSB);
				outShift = outLSB - pLSB;
			}

			svR = svP + svA;

			//align the multiply-and-add with the output format
			if(outShift > 0)
			{
				svR = svR >> outShift;
				possibleOutputs = 2;
			}
			else
			{
				svR = svR << (-outShift);
				possibleOutputs = 1;
			}

			// manage two's complement at output
			if(svR < 0)
				svR += twoToWR;

			//add only the bits corresponding to the output format
			svRAux = svR & (twoToWR -1);

			tc->addExpectedOutput("R", svRAux);
			if(possibleOutputs == 2)
			{
				svR++;
				svR &= (twoToWR -1);
				tc->addExpectedOutput("R", svR);
			}
		}
	}



	void FixMultAdd::buildStandardTestCases(TestCaseList* tcl)
	{
		//TODO
	}




}
