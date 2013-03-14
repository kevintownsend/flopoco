/*
 * Operator computing (X^3)/6, for X being a signed two's complement number
 
 Assumptions: 0<=X<1

 This file is part of the FloPoCo project developed by the Arenaire/ARIC
 team at Ecole Normale Superieure de Lyon
  
 Author : Florent de Dinechin, Matei Istoan

 Initial software.
 Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
 2008-2013.
  All rights reserved.
*/

#include "FixSinPoly.hpp"

#ifdef HAVE_SOLLYA
#include "../sollya.h"

using namespace std;

namespace flopoco{


	//FIXME: for now, the output width (wOut) is computed inside the operator, 
	//		and not from the parameters given to the constructor
	
	//standalone operator 
	FixSinPoly::FixSinPoly(Target* target, int msbIn_, int lsbIn_, int msbOut_, int lsbOut_, bool signedInput_, map<string, double> inputDelays) :
		Operator(target, inputDelays), msbIn(msbIn_), lsbIn(lsbIn_), msbOut(msbOut_), lsbOut(lsbOut_), signedInput(signedInput_)
	{
		srcFileName="FixSinPoly";

		if(lsbIn > msbIn)
			throw string("FixSinPoly: Error, lsbIn should not be greater than msbIn");
    
		if(lsbIn > msbIn)
			throw string("FixSinPoly: Error, lsbOut should not be greater than msbOut");
		
		wIn = msbIn - lsbIn + 1;
		wOut = (msbIn<0 ? msbIn : 3*msbIn) - (lsbIn<0 ? 3*lsbIn : lsbIn) + 1 + (lsbIn<0 ? 1 : 0);				//one more bit to cover all the range of produced bits
		
		// build the name
		ostringstream name; 
		name << "FixSinPoly_" << vhdlize(wIn) << "_" << vhdlize(wOut) << (signedInput ? "_signed" : "_unsigned");
		setName(name.str());
		
		//create the input and the output
		addInput("X", wIn);
		addOutput("R", wOut);
		
		//create the bitheap that computes the sum
		bitHeap = new BitHeap(this, wOut);
		
		REPORT(DEBUG, "Adding the bits for X");
		
		//add the bits corresponding to sum_{i=imin}^imax(2^i*x_i)
		for(int i=lsbIn; i<=msbIn; i++)
		{
			stringstream s;
			
			vhdl << tab << declare(join("X_orig_", (lsbIn<0 ? i-lsbIn : i))) 
				<< " <= X" << of(lsbIn<0 ? i-lsbIn : i) << ";" << endl;
						
			s << join("X_orig_", (lsbIn<0 ? i-lsbIn : i));
			
			if(msbIn<0 && lsbIn<0)
			{
				bitHeap->addBit(wOut-1-(msbIn-i), s.str());
			}
			else if(lsbIn<0)
			{
				bitHeap->addBit((i-lsbIn)-2*lsbIn+1, s.str());
			}
			else
			{
				bitHeap->addBit(i, s.str());
			}
		}
		
		REPORT(DEBUG, "Adding the bits for the first sum");
		
		//add the terms corresponding to sum_{i=imin}^imax(2^(3i-1)*x_i)
		//	negated
		IntConstDiv3 *divider;
		divider = new IntConstDiv3(target, wIn, 3, -1, 2);
		addSubComponent(divider);

		//manageCriticalPath(divider->getOutputDelay("Q"));
  
		inPortMap (divider , "X", "X");
		outPortMap(divider , "Q", "XZeroIntDiv3");
		vhdl << instance(divider , "Divider");
		
		for(int i=0; i<=((msbIn-lsbIn+1)*3-2)-1; i++)
		{
			stringstream s;
			
			vhdl << tab << declare(join("XZeroIntDiv3_inverted_", i)) << " <= not XZeroIntDiv3" << of(i) << ";" << endl;
			
			s << "XZeroIntDiv3_inverted_" << i;
			
			bitHeap->addBit(i, s.str());
			
			for(int j=i; j<wOut; j++)
				bitHeap->addConstantOneBit(j);
		}
		
		REPORT(DEBUG, "Adding the bits for the second sum");
		
		//add the terms corresponding to sum_i_j_imin^imax(2^(i+2j-1)*x_i*x_j)
		//	negated
		for(int i=lsbIn; i<=msbIn; i++)
			for(int j=lsbIn; j<=msbIn; j++)
			{
				//add the bit only if i != j
				if(i == j)
					continue;
				
				stringstream s;
				
				vhdl << tab << declare(join("X_temp_", (lsbIn<0 ? i-lsbIn : i), "_", (lsbIn<0 ? j-lsbIn : j))) 
					<< " <= not (X" << of(lsbIn<0 ? i-lsbIn : i) << " and X" << of(lsbIn<0 ? j-lsbIn : j) << ");" << endl;
							
				s << "X_temp_" << (lsbIn<0 ? i-lsbIn : i) << "_" << (lsbIn<0 ? j-lsbIn : j);
				
				bitHeap->addBit((lsbIn<0 ? (i-lsbIn)+2*(j-lsbIn)-1 : i+2*j-1) + 1, s.str());
				
				for(int k=(lsbIn<0 ? (i-lsbIn)+2*(j-lsbIn)-1 : i+2*j-1) + 1; k<wOut; k++)
					bitHeap->addConstantOneBit(k);
			}
		
		REPORT(DEBUG, "Adding the bits for the third sum");
		
		//add the terms corresponding to sum_i_j_k_imin^imax(2^(i+j+k)*x_i*x_j*x_k)
		//	negated
		for(int i=lsbIn; i<=msbIn; i++)
			for(int j=i+1; j<=msbIn; j++)
				for(int k=j+1; k<=msbIn; k++)
				{
					//i < j < k
					stringstream s;					
					
					vhdl << tab << declare(join("X_temp2_", (lsbIn<0 ? i-lsbIn : i), "_", (lsbIn<0 ? j-lsbIn : j), "_", (lsbIn<0 ? k-lsbIn : k))) 
						<< " <= not (X" << of(lsbIn<0 ? i-lsbIn : i) << " and X" << of(lsbIn<0 ? j-lsbIn : j) << " and X" << of(lsbIn<0 ? k-lsbIn : k) << ");" << endl;
					
					s << "X_temp2_" << (lsbIn<0 ? i-lsbIn : i) << "_" << (lsbIn<0 ? j-lsbIn : j) << "_" << (lsbIn<0 ? k-lsbIn : k);
					
					bitHeap->addBit((lsbIn<0 ? (i-lsbIn)+(j-lsbIn)+(k-lsbIn) : i+j+k) + 1, s.str());
					
					for(int l=(lsbIn<0 ? (i-lsbIn)+(j-lsbIn)+(k-lsbIn) : i+j+k) + 1; l<wOut; l++)
						bitHeap->addConstantOneBit(l);
				}
				
		//compress the bitheap
		bitHeap -> generateCompressorVHDL();
		
		//generate the final result
		vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << range(wOut-1, 0) << ";" << endl;
						
	}
	
	
	//operator incorporated into a global compression
	//	for use as part of a bigger operator
	FixSinPoly::FixSinPoly(Operator* parentOp_, Target* target, Signal* multiplicandX, int msbIn_, int lsbIn_, int msbOut_, int lsbOut_, BitHeap* bitHeap_, bool signedInput_, map<string, double> inputDelays) :
		Operator(target, inputDelays), msbIn(msbIn_), lsbIn(lsbIn_), msbOut(msbOut_), lsbOut(lsbOut_), signedInput(signedInput_),
		wIn(msbIn_-lsbIn+1), wOut(msbOut_-lsbOut+1),
		parentOp(parentOp_), bitHeap(bitHeap_) 
	{
		srcFileName="FixSinPoly";

		// build the name
		ostringstream name; 
		name <<"FixSinPoly_" << vhdlize(wIn) << "_" << vhdlize(wOut) << "_" << (signedInput?"_signed":"_unsigned");
		setName(name.str()); 

		
	}




	FixSinPoly::~FixSinPoly()
	{
		// TODO 
	}


	void FixSinPoly::emulate(TestCase* tc)
	{
		// get I/O values
		mpz_class svX = tc->getInputValue("X");
		mpz_class svR;
		
		//compute the value of the result
		svR = svX*svX*svX;
		svR = svR/3;
		svR = (svX << (lsbIn-3*lsbIn + 1)) - svR;
		
		// add the result
		tc->addExpectedOutput("R", svR);
	}

}




#endif //HAVE_SOLLYA
