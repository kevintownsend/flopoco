/*
  An integer squarer for FloPoCo, using bit heaps

  Author: Bogdan Pasca, Matei Istoan

  This file is part of the FloPoCo project

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2014.
  All rights reserved.
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"
#include "IntSquarerBitheap.hpp"

#include "BitHeap.hpp"

//only one of these constants can be 1 at any given time
#define COMPRESSED_SQR_ARCH 1
#define FULL_SQR_ARCH 0

using namespace std;

namespace flopoco{

		IntSquarerBitheap::IntSquarerBitheap(Target* target, int wIn, map<string, double> inputDelays):
		Operator(target, inputDelays), wIn_(wIn), inputDelays_(inputDelays)
	{
		ostringstream name;
		name << "IntSquarerBitheap_" << wIn_ << "_uid" << getNewUId();
		setName(name.str());
		setCopyrightString("Bogdan Pasca, Matei Istoan (2014)");

		srcFileName = "IntSquarerBitheap";

		// Set up the IO signals
		addInput ("X"  , wIn_);
		addOutput("R"  , 2*wIn_);

		setCriticalPath( getMaxInputDelays(inputDelays) );

#if FULL_SQR_ARCH
		if (wIn <= 17 )
		{
			//no point in using a bit heap here
			//	the only thing needed is an addition

			vhdl << tab << declare( "sX", wIn) << " <= X;" << endl;
			vhdl << tab << declare( "sY", wIn) << " <= X;" << endl;

			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );
			vhdl << tab << "R <= sX * sY;" << endl;
			outDelayMap["R"] = getCriticalPath();
		}
		else if ((wIn>17) && (wIn<=34))
		{
			//create the bit heap
			bitHeap = new BitHeap(this, 2*wIn);

			//prepare the groups of bits
			vhdl << tab << declare("x0_16",18) << " <= \"0\" & X" << range(16,0) << ";" << endl;
			vhdl << tab << declare("x17_32",18) << " <= " << (wIn==34 ? zg(1, 0) : zg(34-wIn, 0)) << " & X" << range(wIn-1,17) << ";" << endl;

			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );
			vhdl << tab << declare("p0",36) << " <= x0_16 * x0_16;" <<endl;
			//add the bits to the bit heap
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "p0" << of(i);

				bitHeap->addBit(i, s.str());
			}

			vhdl << tab << declare("p1_x2",36) << " <= x17_32 * x0_16;" <<endl;
			//add the bits to the bit heap
			int bitHeapOffset = 17+1;				//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "p1_x2" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}

			manageCriticalPath( target->DSPCascadingWireDelay() + target->DSPAdderDelay() );
			vhdl << tab << declare("p2",36) << " <= x17_32 * x17_32;" <<endl;
			//add the bits to the bit heap
			bitHeapOffset = 2*17;				//each of the multiplicands is offset
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "p2" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}

			//compress the bit heap and produce the result
			bitHeap->generateCompressorVHDL();

			//manage the pipeline
			syncCycleFromSignal(bitHeap->getSumName());

			//because of final add in bit heap, add one more bit to the result
			vhdl << declare("Rint", 2*wIn) << " <= " << bitHeap->getSumName() << range(2*wIn-1, 0) << ";" << endl;

			vhdl << tab << "R <= Rint;" << endl;

			outDelayMap["R"] = getCriticalPath();

		}
		else if ((wIn>34) && (wIn<=51))
		{
			int bitHeapOffset;
			int sliceOffset = 17;		//the size of a slice of the input

			//create the bit heap
			bitHeap = new BitHeap(this, 2*wIn);

			if (wIn<51)
				vhdl << tab << declare("sigX",51) << "<= " << zg(51-wIn,0) << " & X;" << endl;
			else
				vhdl << tab << declare("sigX",51) << "<= X;" << endl;

			vhdl << tab << declare("x0_16",  18) << " <= \"0\" & sigX" << range(16,  0) << ";" << endl;
			vhdl << tab << declare("x17_33", 18) << " <= \"0\" & sigX" << range(33, 17) << ";" << endl;
			vhdl << tab << declare("x34_50", 18) << " <= \"0\" & sigX" << range(50, 34) << ";" << endl;

			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );
			vhdl << tab << declare("x0_16_sqr", 36)     << "<=  x0_16  * x0_16;" << endl;
			vhdl << tab << declare("x17_33_sqr", 36)    << "<= x17_33 * x17_33;" << endl;
			vhdl << tab << declare("x34_50_sqr", 36)    << "<= x34_50 * x34_50;" << endl;
			vhdl << tab << declare("x0_16_x17_33", 36)  << "<= x0_16 * x17_33;"  << endl;
			vhdl << tab << declare("x0_16_x34_50", 36)  << "<= x0_16 * x34_50;"  << endl;
			vhdl << tab << declare("x17_33_x34_50", 36) << "<= x17_33 * x34_50;" << endl;


			//add the bits of x0_16_sqr to the bit heap
			bitHeapOffset = 0;
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of x17_33_sqr to the bit heap
			bitHeapOffset = 2*sliceOffset;
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x17_33_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of x34_50_sqr to the bit heap
			bitHeapOffset = 2*(2*sliceOffset);
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x34_50_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}


			//add the bits of 2*x0_16_x17_33 to the bit heap
			bitHeapOffset = 0+sliceOffset+1;							//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_x17_33" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x0_16_x34_50 to the bit heap
			bitHeapOffset = 0+2*sliceOffset+1;							//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_x34_50" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x17_33_x34_50 to the bit heap
			bitHeapOffset = sliceOffset+2*sliceOffset+1;				//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x17_33_x34_50" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}

			//compress the bit heap and produce the result
			bitHeap->generateCompressorVHDL();

			//manage the pipeline
			syncCycleFromSignal(bitHeap->getSumName());

			//because of final add in bit heap, add one more bit to the result
			vhdl << declare("Rint", 2*wIn) << " <= " << bitHeap->getSumName() << range(2*wIn-1, 0) << ";" << endl;

			vhdl << tab << "R <= Rint;" << endl;

			outDelayMap["R"] = getCriticalPath();
		}
		else if ((wIn>51) && (wIn<=68))
		{
			int bitHeapOffset;
			int sliceOffset = 17;		//the size of a slice of the input

			//create the bit heap
			bitHeap = new BitHeap(this, 2*wIn);

			if (wIn<68)
				vhdl << tab << declare("sigX",68) << "<= " << zg(68-wIn,0) << " & X;" << endl;
			else
				vhdl << tab << declare("sigX",68) << "<= X;" << endl;

			vhdl << tab << declare("x0_16",  18) << " <= \"0\" & sigX" << range(16,  0) << ";" << endl;
			vhdl << tab << declare("x17_33", 18) << " <= \"0\" & sigX" << range(33, 17) << ";" << endl;
			vhdl << tab << declare("x34_50", 18) << " <= \"0\" & sigX" << range(50, 34) << ";" << endl;
			vhdl << tab << declare("x51_67", 18) << " <= \"0\" & sigX" << range(67, 51) << ";" << endl;

			manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );
			vhdl << tab << declare("x0_16_sqr", 36)     << "<=  x0_16  * x0_16;" << endl;
			vhdl << tab << declare("x17_33_sqr", 36)    << "<= x17_33 * x17_33;" << endl;
			vhdl << tab << declare("x34_50_sqr", 36)    << "<= x34_50 * x34_50;" << endl;
			vhdl << tab << declare("x51_67_sqr", 36)    << "<= x51_67 * x51_67;" << endl;

			vhdl << tab << declare("x0_16_x17_33", 36)  << "<= x0_16 * x17_33;"  << endl;
			vhdl << tab << declare("x0_16_x34_50", 36)  << "<= x0_16 * x34_50;"  << endl;
			vhdl << tab << declare("x0_16_x51_67", 36)  << "<= x0_16 * x51_67;"  << endl;
			vhdl << tab << declare("x17_33_x34_50", 36) << "<= x17_33 * x34_50;" << endl;
			vhdl << tab << declare("x17_33_x51_67", 36) << "<= x17_33 * x51_67;" << endl;
			vhdl << tab << declare("x34_50_x51_67", 36) << "<= x34_50 * x51_67;" << endl;


			//add the bits of x0_16_sqr to the bit heap
			bitHeapOffset = 0;
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of x17_33_sqr to the bit heap
			bitHeapOffset = 2*sliceOffset;
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x17_33_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of x34_50_sqr to the bit heap
			bitHeapOffset = 2*(2*sliceOffset);
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x34_50_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of x51_67_sqr to the bit heap
			bitHeapOffset = 2*(3*sliceOffset);
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x51_67_sqr" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}


			//add the bits of 2*x0_16_x17_33 to the bit heap
			bitHeapOffset = 0+sliceOffset+1;							//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_x17_33" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x0_16_x34_50 to the bit heap
			bitHeapOffset = 0+2*sliceOffset+1;							//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_x34_50" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x0_16_x51_67 to the bit heap
			bitHeapOffset = 0+3*sliceOffset+1;							//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x0_16_x51_67" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x17_33_x34_50 to the bit heap
			bitHeapOffset = sliceOffset+2*sliceOffset+1;				//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x17_33_x34_50" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x17_33_x51_67 to the bit heap
			bitHeapOffset = sliceOffset+3*sliceOffset+1;				//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x17_33_x51_67" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}
			//add the bits of 2*x34_50_x51_67 to the bit heap
			bitHeapOffset = 2*sliceOffset+3*sliceOffset+1;				//+1 because the sub-product appears twice, so is multiplied by 2
			for(int i=0; i<36; i++)
			{
				stringstream s;

				s << "x34_50_x51_67" << of(i);

				bitHeap->addBit(i+bitHeapOffset, s.str());
			}

			//compress the bit heap and produce the result
			bitHeap->generateCompressorVHDL();

			//manage the pipeline
			syncCycleFromSignal(bitHeap->getSumName());

			//because of final add in bit heap, add one more bit to the result
			vhdl << declare("Rint", 2*wIn) << " <= " << bitHeap->getSumName() << range(2*wIn-1, 0) << ";" << endl;

			vhdl << tab << "R <= Rint;" << endl;

			outDelayMap["R"] = getCriticalPath();
		}
		else
		{
			cerr << " For the moment IntSquarer does not support inputs larger than 68 bits. " << endl;
			exit (EXIT_FAILURE);
		}
#endif

#if COMPRESSED_SQR_ARCH

		int nbSlices, sliceSize, bitHeapOffset;

		//create the bit heap
		bitHeap = new BitHeap(this, 2*wIn);

		//determine the number of parts into which to split the input
		bitHeapOffset = 0;
		sliceSize = 17;				//the size of a DSP on Xilinx architectures
		nbSlices = wIn/17 + 1;
		if(wIn%17 == 0)
			nbSlices--;

		//pad the input, for easier processing
		if (wIn<sliceSize*nbSlices)
			vhdl << tab << declare("sigX", sliceSize*nbSlices) << " <= " << zg(sliceSize*nbSlices-wIn, 0) << " & X;" << endl;
		else
			vhdl << tab << declare("sigX", sliceSize*nbSlices) << "<= X;" << endl;

		//split the input into slices
		for(int i=0; i<nbSlices; i++)
		{
			vhdl << tab << declare(join("x", i*sliceSize, "_", (i+1)*sliceSize-1),  sliceSize+1)
					<< " <= \"0\" & sigX" << range((i+1)*sliceSize-1,  i*sliceSize) << ";" << endl;
		}

		manageCriticalPath( target->LogicToDSPWireDelay() + target->DSPMultiplierDelay() );

		//create the squares
		for(int i=0; i<nbSlices; i++)
		{
			vhdl << tab << declare(join("x", i*sliceSize, "_", (i+1)*sliceSize-1, "_sqr"), 2*(sliceSize+1))
					<< "<= " << join("x", i*sliceSize, "_", (i+1)*sliceSize-1)
					<< " * " << join("x", i*sliceSize, "_", (i+1)*sliceSize-1) << ";" << endl;
		}

		//create the rest of the sub-products
		for(int i=0; i<nbSlices; i++)
			for(int j=i+1; j<nbSlices; j++)
			{
				vhdl << tab << declare(join("x", i*sliceSize, "_", (i+1)*sliceSize-1, "_x", j*sliceSize, "_", (j+1)*sliceSize-1), 2*(sliceSize+1))
						<< "<= " << join("x", i*sliceSize, "_", (i+1)*sliceSize-1)
						<< " * " << join("x", j*sliceSize, "_", (j+1)*sliceSize-1) << ";" << endl;
			}

		//add the bits of the squares to the bit heap
		for(int i=0; i<nbSlices; i++)
		{
			bitHeapOffset = 2*i*sliceSize;

			for(int j=0; j<2*sliceSize; j++)
			{
				stringstream s;

				s << "x" << i*sliceSize << "_" << (i+1)*sliceSize-1 << "_sqr" << of(j);

				bitHeap->addBit(j+bitHeapOffset, s.str());
			}
		}

		//add the bits of the sub-products to the bit heap
		for(int i=0; i<nbSlices; i++)
			for(int j=i+1; j<nbSlices; j++)
			{
				bitHeapOffset = (i+j)*sliceSize+1;					//+1 because the sub-product is added twice, so the multiplication result is multiplied by 2

				for(int k=0; k<2*sliceSize; k++)
				{
					stringstream s;

					s << "x" << i*sliceSize << "_" << (i+1)*sliceSize-1 << "_x" << j*sliceSize << "_" << (j+1)*sliceSize-1 << of(k);

					bitHeap->addBit(k+bitHeapOffset, s.str());
				}
			}

		//compress the bit heap and produce the result
		bitHeap->generateCompressorVHDL();

		//manage the pipeline
		syncCycleFromSignal(bitHeap->getSumName());

		//because of final add in bit heap, add one more bit to the result
		vhdl << declare("Rint", 2*wIn) << " <= " << bitHeap->getSumName() << range(2*wIn-1, 0) << ";" << endl;

		vhdl << tab << "R <= Rint;" << endl;

		outDelayMap["R"] = getCriticalPath();

#endif


	}

	IntSquarerBitheap::~IntSquarerBitheap() {
	}


	void IntSquarerBitheap::outputVHDL(std::ostream& o, std::string name) {
		ostringstream signame;
		licence(o);
		pipelineInfo(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.std_logic_arith.all;" << endl;
		if ((wIn_>17) && (wIn_<34)) {
			o << "use ieee.std_logic_signed.all;" << endl;
		}else
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




	void IntSquarerBitheap::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svR = svX * svX ;
		tc->addExpectedOutput("R", svR);
	}


}
