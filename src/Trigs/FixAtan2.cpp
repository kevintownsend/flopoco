/*
   An arctangent(y/x) implementation

	Author:  Florent de Dinechin, Matei Istoan

	This file is part of the FloPoCo project

	Initial software.
	Copyright © INSA Lyon, INRIA, CNRS, UCBL,
	2014.
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

#include "../utils.hpp"
#include "Operator.hpp"
#include "FixAtan2.hpp"

using namespace std;

namespace flopoco {


	// The constructor for a stand-alone operator
	FixAtan2::FixAtan2(Target* target_, int wIn_, int wOut_, int architectureType_, map<string, double> inputDelays_):
		Operator(target_, inputDelays_),
		target(target_), wIn(wIn_), wOut(wOut_), architectureType(architectureType_)
	{

		srcFileName = "FixAtan2";
		setCopyrightString("Florent de Dinechin, Matei Istoan, 2014");

		//set the VHDL generation style
		plainVHDL =  true;

		//set the ratio for the multiplications
		ratio = 0.95;

		// build the name
		ostringstream name;
		name <<"FixAtan2_" << vhdlize(wIn) << "_" << vhdlize(wOut) << "_archType_" << vhdlize(architectureType_);
		setName(name.str());

		//set the libraries to use
		//useNumericStd();
		useNumericStd_Unsigned();

		//initialize global variables
		maxValA = -1;		//this is the absolute value, so should be always non-negative, once used
		maxValB = -1;
		maxValC = -1;
		maxValD = -1;
		maxValE = -1;
		maxValF = -1;

		//declare the inputs and the outputs
		addInput ("X",  wIn);
		addInput ("Y",  wIn);
		addOutput("R",  wOut, 2 /*number of possible output values*/);

		//arguments for the range reduction
		bool negateByComplement = true;

		/////////////////////////////////////////////////////////////////////////////
		//
		//    First range reduction
		//
		/////////////////////////////////////////////////////////////////////////////
		vhdl << tab << declare("sgnX") << " <= X" << of(wIn-1) << ";" << endl;
		vhdl << tab << declare("sgnY") << " <= Y" << of(wIn-1) << ";" << endl;

		// TODO: replace the following with LUT-based comparators
		// and first maybe experiment with synthesis tools
		vhdl << tab << declare("XmY", wIn+1) << " <= std_logic_vector(sgnX & X) - std_logic_vector(sgnY & Y);" << endl;
		vhdl << tab << declare("XpY", wIn+1) << " <= (sgnX & X)+(sgnY & Y);" << endl;
		vhdl << tab << declare("XltY") << " <= XmY" << of(wIn) << ";" << endl;
		vhdl << tab << declare("mYltX") << " <= not XpY" << of(wIn) <<";" << endl;
		// Range reduction: we define 4 quadrants, each centered on one axis (these are not just the sign quadrants)
		// Then each quadrant is decomposed in its positive and its negative octant.
		vhdl << tab << "-- quadrant will also be the angle to add at the end" << endl;
		vhdl << tab << declare("quadrant", 2) << " <= " << endl;
		vhdl << tab << tab << "\"00\"  when (not sgnX and not XltY and     mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"01\"  when (not sgnY and     XltY and     mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"10\"  when (    sgnX and     XltY and not mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"11\";"    << endl;

		if(negateByComplement)
		{
			vhdl << tab << declare("pX", wIn) << " <=      X;" << endl;
			vhdl << tab << declare("pY", wIn) << " <=      Y;" << endl;
			vhdl << tab << declare("mX", wIn) << " <= (not X);  -- negation by not, implies one ulp error." << endl;
			vhdl << tab << declare("mY", wIn) << " <= (not Y);  -- negation by not, implies one ulp error. " << endl;
		}else
		{
			vhdl << tab << declare("pX", wIn) << " <= X;" << endl;
			vhdl << tab << declare("pY", wIn) << " <= Y;" << endl;
			vhdl << tab << declare("mX", wIn) << " <= (" << zg(wIn) << " - X);" << endl;
			vhdl << tab << declare("mY", wIn) << " <= (" << zg(wIn) << " - Y);" << endl;
		}

		//no need for sign bit any longer
		vhdl << tab << declare("XR", wIn-1) << " <= " << endl;
		vhdl << tab << tab << "pX" << range(wIn-2, 0) << " when quadrant=\"00\"   else " << endl;
		vhdl << tab << tab << "pY" << range(wIn-2, 0) << " when quadrant=\"01\"   else " << endl;
		vhdl << tab << tab << "mX" << range(wIn-2, 0) << " when quadrant=\"10\"   else " << endl;
		vhdl << tab << tab << "mY" << range(wIn-2, 0) << ";"    << endl;

		//no need for sign bit any longer
		vhdl << tab << declare("YR", wIn-1) << " <= " << endl;
		vhdl << tab << tab << "pY" << range(wIn-2, 0) << " when quadrant=\"00\" and sgnY='0'  else " << endl;
		vhdl << tab << tab << "mY" << range(wIn-2, 0) << " when quadrant=\"00\" and sgnY='1'  else " << endl;
		vhdl << tab << tab << "pX" << range(wIn-2, 0) << " when quadrant=\"01\" and sgnX='0'  else " << endl;
		vhdl << tab << tab << "mX" << range(wIn-2, 0) << " when quadrant=\"01\" and sgnX='1'  else " << endl;
		vhdl << tab << tab << "pY" << range(wIn-2, 0) << " when quadrant=\"10\" and sgnY='0'  else " << endl;
		vhdl << tab << tab << "mY" << range(wIn-2, 0) << " when quadrant=\"10\" and sgnY='1'  else " << endl;
		vhdl << tab << tab << "pX" << range(wIn-2, 0) << " when quadrant=\"11\" and sgnX='0'  else "    << endl;
		vhdl << tab << tab << "mX" << range(wIn-2, 0) << " ;"    << endl;

		vhdl << tab << declare("finalAdd") << " <= " << endl;
		vhdl << tab << tab << "'1' when (quadrant=\"00\" and sgnY='0') or(quadrant=\"01\" and sgnX='1') or (quadrant=\"10\" and sgnY='1') or (quadrant=\"11\" and sgnX='0')" << endl;
		vhdl << tab << tab << " else '0';  -- this information is sent to the end of the pipeline, better compute it here as one bit" << endl;
		/////////////////////////////////////////////////////////////////////////////
		//
		//    End of first range reduction
		//
		/////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////
		//
		//    Second range reduction, scaling
		//
		////////////////////////////////////////////////////////////////////////////
		vhdl << tab << declare("XorY", wIn-2) << " <= XR" << range(wIn-2,1) << " or YR" << range(wIn-2,1) << ";" << endl;
		// The LZC
		LZOC* lzc = new	LZOC(target, wIn-2);
		addSubComponent(lzc);

		inPortMap(lzc, "I", "XorY");
		inPortMapCst(lzc, "OZB", "'0'");
		outPortMap(lzc, "O", "S");
		vhdl << instance(lzc, "lzc");
		//setCycleFromSignal("lzo");
		//		setCriticalPath( lzc->getOutputDelay("O") );

		// The two shifters are two instance of the same component
		Shifter* lshift = new Shifter(target, wIn-1, wIn-2, Shifter::Left);
		addSubComponent(lshift);

		inPortMap(lshift, "S", "S");

		inPortMap(lshift, "X", "XR");
		outPortMap(lshift, "R", "XRSfull");
		vhdl << instance(lshift, "Xshift");
		vhdl << tab << declare("XRS", wIn-1) << " <=  XRsfull " << range(wIn-2,0) << ";" << endl;
		//eliminate the MSB, which is a constant '1'
		vhdl << tab << declare("XRS_short", wIn-2) << " <=  XRS " << range(wIn-3,0) << ";" << endl;

		inPortMap(lshift, "X", "YR");
		outPortMap(lshift, "R", "YRSfull");
		vhdl << instance(lshift, "Yshift");
		vhdl << tab << declare("YRS", wIn-1) << " <=  YRsfull " << range(wIn-2,0) << ";" << endl;
		////////////////////////////////////////////////////////////////////////////
		//
		//    End of second range reduction, scaling
		//
		////////////////////////////////////////////////////////////////////////////


		//build the architecture
		if((architectureType == 0) || (architectureType == 1))
		{
			//based on the plane's equation or using a first order polynomial
			//	(same architecture, but the parameters are generated differently)

			//check if the selected type of architecture can achieve the required accuracy
			//	also, determine the size of the parts of the input
			//	signals to use as an address for the table
			int k;
			int guardBitsSum, guardBitsApprox;

			guardBitsSum = 2;				// the error when computing Ax+By+C is at most 2.5 ulps, so we'll need 2 extra guard bits
			guardBitsApprox = 2;			// determined when computing the parameters for the tables
			g = guardBitsApprox + guardBitsSum;

			k = checkArchitecture(architectureType);
			kSize = k;

			//determine the size of the constants to store in the table
			// maxValA, maxValB and maxValC should have been set by the call to checkArchitecture()
			int msbA = intlog2(maxValA)+1;
			int msbB = intlog2(maxValB)+1;
			int msbC = intlog2(maxValC)+1;
			int maxMSB = maxInt(3, msbA, msbB, msbC);

			if(plainVHDL)
			{
				//split the input signals, and create the address signal for the table

				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(XRS_short" << range(wIn-3, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(YRS" << range(wIn-2, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;

				/*
				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(X" << range(wIn-3, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(Y" << range(wIn-2, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;
				*/

				//create the table for atan(y/x)
				Atan2Table *table = new Atan2Table(target, 2*k-1, msbA+msbB+msbC+3*(wOut-1+g), architectureType, msbA, msbB, msbC);

				//add the table to the operator
				addSubComponent(table);
				useSoftRAM(table);
				//useHardRAM(table);
				inPortMap (table , "X", "atan2TableInput");
				outPortMap(table , "Y", "atan2TableOutput");
				vhdl << instance(table , "KCMTable");

				//split the output of the table to the corresponding parts A, B and C
				vhdl << tab << declareFixPoint("C", true, msbC-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbC+(wOut-1+g)-1, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("B", true, msbB-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbB+msbC+2*(wOut-1+g)-1, msbC+(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("A", true, msbA-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbA+msbB+msbC+3*(wOut-1+g)-1, msbB+msbC+2*(wOut-1+g)) << ");" << endl;

				vhdl << tab << declareFixPoint("XLow", true, -k, -wIn+1) << " <= signed('0' & X" << range(wIn-1-k-1, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("YLow", true, -k, -wIn+1) << " <= signed('0' & Y" << range(wIn-1-k-1, 0) << ");" << endl;

				//create A*X_low and B*Y_low
				vhdl << tab << declareFixPoint("AXLow", true, msbA-k, -wOut+1-g-wIn+1) << " <= A * XLow;" << endl;
				vhdl << tab << declareFixPoint("BYLow", true, msbB-k, -wOut+1-g-wIn+1) << " <= B * YLow;" << endl;
				//align A*X_low and B*Y_low to the output format
				resizeFixPoint("AXLow_sgnExtended", "AXLow", maxMSB-1, -wOut+1-g);
				resizeFixPoint("BYLow_sgnExtended", "BYLow", maxMSB-1, -wOut+1-g);

				//add everything up
				vhdl << tab << declareFixPoint("AXLowAddBYLow", true, maxMSB, -wOut+1-g)
						<< " <= (AXLow_sgnExtended(AXLow_sgnExtended'HIGH) & AXLow_sgnExtended) + (BYLow_sgnExtended(BYLow_sgnExtended'HIGH) & BYLow_sgnExtended);" << endl;
				resizeFixPoint("C_sgnExtended", "C", maxMSB+1, -wOut+1-g);
				vhdl << tab << declareFixPoint("AXLowAddBYLowAddC", true, maxMSB+1, -wOut+1-g) << " <= (AXLowAddBYLow(AXLowAddBYLow'HIGH) & AXLowAddBYLow) + C_sgnExtended;" << endl;

				//extract the final result
				resizeFixPoint("Rtmp", "AXLowAddBYLowAddC", 1, -wOut);
				vhdl << tab << declareFixPoint("Rtmp_rndCst", true, 1, -wOut) << " <= signed(std_logic_vector\'(\"" << zg(wOut+1, -2) << "1\"));" << endl;
				vhdl << tab << declareFixPoint("Rtmp_rnd", true, 1, -wOut) << " <= Rtmp + Rtmp_rndCst;" << endl;

				//return the result
				/*
				resizeFixPoint("Rtmp_stdlv", "Rtmp_rnd", 0, -wOut+1);
				vhdl << tab << declare("R_int", wOut) << " <= std_logic_vector(Rtmp_stdlv);" << endl;
				*/
				resizeFixPoint("Rtmp_stdlv", "Rtmp_rnd", -2, -wOut+1);
				vhdl << tab << declare("R_int", wOut-2) << " <= std_logic_vector(Rtmp_stdlv);" << endl;

				//vhdl << tab << "R <= std_logic_vector(Rtmp_stdlv);" << endl;

			}else
			{
				//create the bitheap
				bitHeap = new BitHeap(this, (maxMSB+2)+wOut+g);

				//create the input signals for the table
				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(X" << range(wIn-3, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(Y" << range(wIn-2, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;

				//create the table for atan(y/x)
				Atan2Table *table = new Atan2Table(target, 2*k-1, msbA+msbB+msbC+3*(wOut-1+g), architectureType, msbA, msbB, msbC);

				//add the table to the operator
				addSubComponent(table);
				useSoftRAM(table);
				//useHardRAM(table);
				inPortMap (table , "X", "atan2TableInput");
				outPortMap(table , "Y", "atan2TableOutput");
				vhdl << instance(table , "KCMTable");

				//extract signals A, B and C
				vhdl << tab << declare("C", msbC+wOut-1+g) << " <= atan2TableOutput"
						<< range(msbC+(wOut-1+g)-1, 0) << ";" << endl;
				vhdl << tab << declare("B", msbB+wOut-1+g) << " <= atan2TableOutput"
						<< range(msbB+msbC+2*(wOut-1+g)-1, msbC+(wOut-1+g)) << ";" << endl;
				vhdl << tab << declare("A", msbA+wOut-1+g) << " <= atan2TableOutput"
						<< range(msbA+msbB+msbC+3*(wOut-1+g)-1, msbB+msbC+2*(wOut-1+g)) << ";" << endl;

				//create Ax and By
				vhdl << tab << declare("XLow", wIn-k) << " <= '0' & X" << range(wIn-1-k-1, 0) << ";" << endl;
				vhdl << tab << declare("YLow", wIn-k) << " <= '0' & Y" << range(wIn-1-k-1, 0) << ";" << endl;

				IntMultiplier* multAx;
				multAx = new IntMultiplier(this,								//parent operator
											 bitHeap,							//the bit heap that performs the compression
											 getSignalByName("XLow"),			//first input to the multiplier (a signal)
											 getSignalByName("A"),				//second input to the multiplier (a signal)
											 -wIn+1,							//offset of the LSB of the multiplier in the bit heap
											 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
											 true,								//signed/unsigned operator
											 ratio);							//DSP ratio
				IntMultiplier* multBy;
				multBy = new IntMultiplier(this,								//parent operator
											 bitHeap,							//the bit heap that performs the compression
											 getSignalByName("YLow"),			//first input to the multiplier (a signal)
											 getSignalByName("B"),				//second input to the multiplier (a signal)
											 -wIn+1,							//offset of the LSB of the multiplier in the bit heap
											 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
											 true,								//signed/unsigned operator
											 ratio);							//DSP ratio

				bitHeap->addSignedBitVector(0,									//weight of signal in the bit heap
											"C",								//name of the signal
											msbC+wOut-1+g,						//size of the signal added
											0,									//index of the lsb in the bit vector from which to add the bits of the addend
											false);								//if we are correcting the index in the bit vector with a negative weight

				//add the rounding bit - take into consideration the final alignment
				bitHeap->addConstantOneBit(g-1);

				//compress the bit heap
				bitHeap -> generateCompressorVHDL();

				//extract the result - take into consideration the final alignment
				vhdl << tab << "R <= " << bitHeap->getSumName() << range(wOut+g-1, g) << ";" << endl;
			}

		}else if(architectureType == 2)
		{
			//based on an order 2 Taylor approximating polynomial

			//check if the selected type of architecture can achieve the required accuracy
			//	also, determine the size of the parts of the input
			//	signals to use as an address for the table
			int k;
			int guardBitsSum, guardBitsApprox;

			guardBitsSum = 2;				// the error when computing Ax+By+C+Dx^2+Ey^2+Fxy is at most 2.5 ulps, so we'll need 2 extra guard bits
			guardBitsApprox = 2;			// determined when computing the parameters for the tables
			g = guardBitsApprox + guardBitsSum;

			k = checkArchitecture(architectureType);
			kSize = k;

			//determine the size of the constants to store in the table
			// maxValA, maxValB and maxValC should have been set by the call to checkArchitecture()
			int msbA = intlog2(maxValA)+1;
			int msbB = intlog2(maxValB)+1;
			int msbC = intlog2(maxValC)+1;
			int msbD = intlog2(maxValD)+1;
			int msbE = intlog2(maxValE)+1;
			int msbF = intlog2(maxValF)+1;
			int maxMSB = maxInt(6, msbA, msbB, msbC, msbD, msbE, msbF);

			if(plainVHDL)
			{
				//split the input signals, and create the address signal for the table
				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(X" << range(wIn-3, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(Y" << range(wIn-2, wIn-1-k) << ");" << endl;
				vhdl << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;
				vhdl << endl;

				//create the table for atan(y/x)
				Atan2Table *table = new Atan2Table(target, 2*k-1, msbA+msbB+msbC+msbD+msbE+msbF+6*(wOut-1+g),
													architectureType, msbA, msbB, msbC, msbD, msbE, msbF);

				//add the table to the operator
				addSubComponent(table);
				useSoftRAM(table);
				//useHardRAM(table);
				inPortMap (table , "X", "atan2TableInput");
				outPortMap(table , "Y", "atan2TableOutput");
				vhdl << instance(table , "KCMTable");
				vhdl << endl;

				//split the output of the table to the corresponding parts A, B, C, D, E and F
				vhdl << tab << declareFixPoint("F", true, msbF-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbF+1*(wOut-1+g)-1, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("E", true, msbE-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbE+msbF+2*(wOut-1+g)-1, msbF+1*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("D", true, msbD-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbD+msbE+msbF+3*(wOut-1+g)-1, msbE+msbF+2*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("C", true, msbC-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbC+msbD+msbE+msbF+4*(wOut-1+g)-1, msbD+msbE+msbF+3*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("B", true, msbB-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbB+msbC+msbD+msbE+msbF+5*(wOut-1+g)-1, msbC+msbD+msbE+msbF+4*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("A", true, msbA-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbA+msbB+msbC+msbD+msbE+msbF+6*(wOut-1+g)-1, msbB+msbC+msbD+msbE+msbF+5*(wOut-1+g)) << ");" << endl;
				vhdl << endl;

				vhdl << tab << declareFixPoint("DeltaX", true, -k-1,  -wIn+1) << " <= signed(not(X(" << of(wIn-k-2) << ")) & X" << range(wIn-k-3, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("DeltaY", true, -k-1,  -wIn+1) << " <= signed(not(Y(" << of(wIn-k-2) << ")) & Y" << range(wIn-k-3, 0) << ");" << endl;
				vhdl << endl;

				//create A*DeltaX and B*DeltaY
				vhdl << tab << declareFixPoint("A_DeltaX", true, msbA-k-1,  -wOut+1-g-wIn+1) << " <= A * DeltaX;" << endl;
				vhdl << tab << declareFixPoint("B_DeltaY", true, msbB-k-1,  -wOut+1-g-wIn+1) << " <= B * DeltaY;" << endl;
				vhdl << endl;

				//create DeltaX^2, DeltaY^2 and DeltaX*DeltaY
				vhdl << tab << declareFixPoint("DeltaX2", true, -2*k-1,  -2*wIn+2) << " <= DeltaX * DeltaX;" << endl;
				vhdl << tab << declareFixPoint("DeltaY2", true, -2*k-1,  -2*wIn+2) << " <= DeltaY * DeltaY;" << endl;
				vhdl << tab << declareFixPoint("DeltaX_DeltaY", true, -2*k-1,  -2*wIn+2) << " <= DeltaX * DeltaY;" << endl;
				vhdl << endl;

				//align the products, discard the extra lsb-s
				resizeFixPoint("DeltaX2_short", "DeltaX2", -2*k-1, -wOut+1-g);
				resizeFixPoint("DeltaY2_short", "DeltaY2", -2*k-1, -wOut+1-g);
				resizeFixPoint("DeltaX_DeltaY_short", "DeltaX_DeltaY", -2*k-1, -wOut+1-g);
				vhdl << endl;

				//create D*DeltaX^2, E*DeltaY^2 and F*DeltaX*DeltaY
				/*
				vhdl << tab << declareFixPoint("D_DeltaX2", true, msbD-2*k-1,  -2*wIn-wOut-g) << " <= D * DeltaX2;" << endl;
				vhdl << tab << declareFixPoint("E_DeltaY2", true, msbE-2*k-1,  -2*wIn-wOut-g) << " <= E * DeltaY2;" << endl;
				vhdl << tab << declareFixPoint("F_DeltaX_DeltaY", true, msbF-2*k-1,  -2*wIn-wOut-g) << " <= F * DeltaX_DeltaY;" << endl;
				*/
				vhdl << tab << declareFixPoint("D_DeltaX2", true, msbD-2*k-1,  -2*(wOut-1+g)) << " <= D * DeltaX2_short;" << endl;
				vhdl << tab << declareFixPoint("E_DeltaY2", true, msbE-2*k-1,  -2*(wOut-1+g)) << " <= E * DeltaY2_short;" << endl;
				vhdl << tab << declareFixPoint("F_DeltaX_DeltaY", true, msbF-2*k-1,  -2*(wOut-1+g)) << " <= F * DeltaX_DeltaY_short;" << endl;
				vhdl << endl;

				//align the signals to the output format to the output format
				resizeFixPoint("A_DeltaX_sgnExt", "A_DeltaX", maxMSB-1, -wOut+1-g);
				resizeFixPoint("B_DeltaY_sgnExt", "B_DeltaY", maxMSB-1, -wOut+1-g);
				resizeFixPoint("C_sgnExt", "C", maxMSB-1, -wOut+1-g);
				resizeFixPoint("D_DeltaX2_sgnExt", "D_DeltaX2", maxMSB-1, -wOut+1-g);
				resizeFixPoint("E_DeltaY2_sgnExt", "E_DeltaY2", maxMSB-1, -wOut+1-g);
				resizeFixPoint("F_DeltaX_DeltaY_sgnExt", "F_DeltaX_DeltaY", maxMSB-1, -wOut+1-g);
				vhdl << endl;

				//add everything up
				vhdl << tab << declareFixPoint("Sum1", true, maxMSB, -wOut+1-g)
						<< " <= (A_DeltaX_sgnExt(A_DeltaX_sgnExt'HIGH) & A_DeltaX_sgnExt) + (B_DeltaY_sgnExt(B_DeltaY_sgnExt'HIGH) & B_DeltaY_sgnExt);" << endl;
				vhdl << tab << declareFixPoint("Sum2", true, maxMSB, -wOut+1-g)
						<< " <= (C_sgnExt(C_sgnExt'HIGH) & C_sgnExt) + (D_DeltaX2_sgnExt(D_DeltaX2_sgnExt'HIGH) & D_DeltaX2_sgnExt);" << endl;
				vhdl << tab << declareFixPoint("Sum3", true, maxMSB, -wOut+1-g)
						<< " <= (E_DeltaY2_sgnExt(E_DeltaY2_sgnExt'HIGH) & E_DeltaY2_sgnExt) + (F_DeltaX_DeltaY_sgnExt(F_DeltaX_DeltaY_sgnExt'HIGH) & F_DeltaX_DeltaY_sgnExt);" << endl;
				vhdl << tab << declareFixPoint("Sum4", true, maxMSB+1, -wOut+1-g)
						<< " <= (Sum1(Sum1'HIGH) & Sum1) + (Sum2(Sum2'HIGH) & Sum2);" << endl;
				vhdl << tab << declareFixPoint("Sum5", true, maxMSB+2, -wOut+1-g)
						<< " <= (Sum4(Sum4'HIGH) & Sum4) + (Sum3(Sum3'HIGH) & Sum3(Sum3'HIGH) & Sum3);" << endl;
				vhdl << endl;

				//extract the final result
				resizeFixPoint("Rtmp", "Sum5", 1, -wOut);
				vhdl << tab << declareFixPoint("Rtmp_rndCst", true, 1, -wOut) << " <= signed(std_logic_vector\'(\"" << zg(wOut+1, -2) << "1\"));" << endl;
				vhdl << tab << declareFixPoint("Rtmp_rnd", true, 1, -wOut) << " <= Rtmp + Rtmp_rndCst;" << endl;

				//return the result
				resizeFixPoint("Rtmp_stdlv", "Rtmp_rnd", 0, -wOut+1);
				vhdl << tab << "R <= std_logic_vector(Rtmp_stdlv);" << endl;
			}else
			{
				//create the bitheap
				bitHeap = new BitHeap(this, (maxMSB+2)+wOut+g);

				//create the input signals for the table
				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(X" << range(wIn-3, wIn-1-k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(Y" << range(wIn-2, wIn-1-k) << ");" << endl;
				vhdl << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;
				vhdl << endl;

				//create the table for atan(y/x)
				Atan2Table *table = new Atan2Table(target, 2*k-1, msbA+msbB+msbC+msbD+msbE+msbF+6*(wOut-1+g),
													architectureType, msbA, msbB, msbC, msbD, msbE, msbF);

				//add the table to the operator
				addSubComponent(table);
				useSoftRAM(table);
				//useHardRAM(table);
				inPortMap (table , "X", "atan2TableInput");
				outPortMap(table , "Y", "atan2TableOutput");
				vhdl << instance(table , "KCMTable");
				vhdl << endl;

				//split the output of the table to the corresponding parts A, B, C, D, E and F
				vhdl << tab << declareFixPoint("F", true, msbF-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbF+1*(wOut-1+g)-1, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("E", true, msbE-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbE+msbF+2*(wOut-1+g)-1, msbF+1*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("D", true, msbD-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbD+msbE+msbF+3*(wOut-1+g)-1, msbE+msbF+2*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("C", true, msbC-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbC+msbD+msbE+msbF+4*(wOut-1+g)-1, msbD+msbE+msbF+3*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("B", true, msbB-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbB+msbC+msbD+msbE+msbF+5*(wOut-1+g)-1, msbC+msbD+msbE+msbF+4*(wOut-1+g)) << ");" << endl;
				vhdl << tab << declareFixPoint("A", true, msbA-1,  -wOut+1-g) << " <= signed(atan2TableOutput"
						<< range(msbA+msbB+msbC+msbD+msbE+msbF+6*(wOut-1+g)-1, msbB+msbC+msbD+msbE+msbF+5*(wOut-1+g)) << ");" << endl;
				vhdl << endl;

				//create DeltaX and DeltaY
				vhdl << tab << declareFixPoint("DeltaX", true, -k-1,  -wIn+1) << " <= signed(not(X(" << of(wIn-k-2) << ")) & X" << range(wIn-k-3, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("DeltaY", true, -k-1,  -wIn+1) << " <= signed(not(Y(" << of(wIn-k-2) << ")) & Y" << range(wIn-k-3, 0) << ");" << endl;
				vhdl << endl;

				//create A*DeltaX
				//	convert the terms to std_logic_vector
				vhdl << tab << declare("DeltaX_stdlv", wIn-1-k) << " <= std_logic_vector(DeltaX);" << endl;
				vhdl << tab << declare("A_stdlv", msbA+wOut-1+g) << " <= std_logic_vector(A);" << endl;
				//	multiply
				IntMultiplier* multADeltaX;
				multADeltaX = new IntMultiplier(this,							//parent operator
											 bitHeap,							//the bit heap that performs the compression
											 getSignalByName("DeltaX_stdlv"),	//first input to the multiplier (a signal)
											 getSignalByName("A_stdlv"),		//second input to the multiplier (a signal)
											 -wIn+1,							//offset of the LSB of the multiplier in the bit heap
											 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
											 true,								//signed/unsigned operator
											 ratio);							//DSP ratio
				vhdl << endl;

				//create B*DeltaY
				//	convert the terms to std_logic_vector
				vhdl << tab << declare("DeltaY_stdlv", wIn-1-k) << " <= std_logic_vector(DeltaY);" << endl;
				vhdl << tab << declare("B_stdlv", msbB+wOut-1+g) << " <= std_logic_vector(B);" << endl;
				//	multiply
				IntMultiplier* multBDeltaY;
				multBDeltaY = new IntMultiplier(this,								//parent operator
												 bitHeap,							//the bit heap that performs the compression
												 getSignalByName("DeltaY_stdlv"),	//first input to the multiplier (a signal)
												 getSignalByName("B_stdlv"),		//second input to the multiplier (a signal)
												 -wIn+1,							//offset of the LSB of the multiplier in the bit heap
												 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
												 true,								//signed/unsigned operator
												 ratio);							//DSP ratio
				vhdl << endl;

				//Add C
				bitHeap->addSignedBitVector(0,									//weight of signal in the bit heap
											"C",								//name of the signal
											msbC+wOut-1+g,						//size of the signal added
											0,									//index of the lsb in the bit vector from which to add the bits of the addend
											false);								//if we are correcting the index in the bit vector with a negative weight
				vhdl << endl;

				//create DeltaX^2, DeltaY^2 and DeltaX*DeltaY
				vhdl << tab << declareFixPoint("DeltaX2", true, -2*k-1,  -2*wIn+2) << " <= DeltaX * DeltaX;" << endl;
				vhdl << tab << declareFixPoint("DeltaY2", true, -2*k-1,  -2*wIn+2) << " <= DeltaY * DeltaY;" << endl;
				vhdl << tab << declareFixPoint("DeltaX_DeltaY", true, -2*k-1,  -2*wIn+2) << " <= DeltaX * DeltaY;" << endl;
				vhdl << endl;

				//align the products, discard the extra lsb-s
				resizeFixPoint("DeltaX2_short", "DeltaX2", -2*k-1, -wOut+1-g);
				resizeFixPoint("DeltaY2_short", "DeltaY2", -2*k-1, -wOut+1-g);
				resizeFixPoint("DeltaX_DeltaY_short", "DeltaX_DeltaY", -2*k-1, -wOut+1-g);
				vhdl << endl;

				//create D*DeltaX^2, E*DeltaY^2 and F*DeltaX*DeltaY

				//create D*DeltaX^2
				//	convert the terms to std_logic_vector
				vhdl << tab << declare("DeltaX2_short_stdlv", wOut-1+g-2*k) << " <= std_logic_vector(DeltaX2_short);" << endl;
				vhdl << tab << declare("D_stdlv", msbD+wOut-1+g) << " <= std_logic_vector(D);" << endl;
				//	multiply
				IntMultiplier* multDDeltaX2;
				multDDeltaX2 = new IntMultiplier(this,								//parent operator
												 bitHeap,							//the bit heap that performs the compression
												 getSignalByName("DeltaX2_short_stdlv"),	//first input to the multiplier (a signal)
												 getSignalByName("D_stdlv"),		//second input to the multiplier (a signal)
												 -wOut+1-g,							//offset of the LSB of the multiplier in the bit heap
												 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
												 true,								//signed/unsigned operator
												 ratio);							//DSP ratio
				vhdl << endl;

				//	create E*DeltaY^2
				//	convert the terms to std_logic_vector
				vhdl << tab << declare("DeltaY2_short_stdlv", wOut-1+g-2*k) << " <= std_logic_vector(DeltaY2_short);" << endl;
				vhdl << tab << declare("E_stdlv", msbE+wOut-1+g) << " <= std_logic_vector(E);" << endl;
				//	multiply
				IntMultiplier* multEDeltaY2;
				multEDeltaY2 = new IntMultiplier(this,								//parent operator
												 bitHeap,							//the bit heap that performs the compression
												 getSignalByName("DeltaY2_short_stdlv"),	//first input to the multiplier (a signal)
												 getSignalByName("E_stdlv"),		//second input to the multiplier (a signal)
												 -wOut+1-g,							//offset of the LSB of the multiplier in the bit heap
												 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
												 true,								//signed/unsigned operator
												 ratio);							//DSP ratio
				vhdl << endl;

				//	create F*DeltaX_DeltaY
				//	convert the terms to std_logic_vector
				vhdl << tab << declare("DeltaX_DeltaY_short_stdlv", wOut-1+g-2*k) << " <= std_logic_vector(DeltaX_DeltaY_short);" << endl;
				vhdl << tab << declare("F_stdlv", msbF+wOut-1+g) << " <= std_logic_vector(F);" << endl;
				//	multiply
				IntMultiplier* multFDeltaXDeltaY;
				multFDeltaXDeltaY = new IntMultiplier(this,								//parent operator
													 bitHeap,							//the bit heap that performs the compression
													 getSignalByName("DeltaX_DeltaY_short_stdlv"),	//first input to the multiplier (a signal)
													 getSignalByName("F_stdlv"),		//second input to the multiplier (a signal)
													 -wOut+1-g,							//offset of the LSB of the multiplier in the bit heap
													 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
													 true,								//signed/unsigned operator
													 ratio);							//DSP ratio
				vhdl << endl;

				//add the rounding bit - take into consideration the final alignment
				bitHeap->addConstantOneBit(g-1);

				//compress the bit heap
				bitHeap -> generateCompressorVHDL();

				//extract the result - take into consideration the final alignment
				vhdl << tab << "R <= " << bitHeap->getSumName() << range(wOut+g-1, g) << ";" << endl;
			}

		}else if(architectureType == 3)
		{
			//based on rotating using a table (for sine and cosine)
			//	and then using the Taylor series for (1/x)
			THROWERROR("in FixAtan2 constructor: architecture not yet implemented");
		}else
		{
			THROWERROR("in FixAtan2 constructor: invalid value for the requested architecture type");
		}

		////////////////////////////////////////////////////////////////////////////
		//
		//                            reconstruction
		//
		////////////////////////////////////////////////////////////////////////////

		vhdl << tab << declare("qangle", wOut) << " <= (quadrant & " << zg(wOut-2) << ");" << endl;
		vhdl << tab << declare("Rfinal", wOut) << " <= \"00\" & R_int; -- sign-extended and rounded" << endl;
		vhdl << tab << "R <= "
				<< tab << tab << "     qangle + Rfinal  when finalAdd='1'" << endl
				<< tab << tab << "else qangle - Rfinal;" << endl;

	}


	FixAtan2::~FixAtan2()
	{

	}

	FixAtan2* FixAtan2::newComponentAndInstanceNumericStd(Operator* op,
															int wIn,
															int wOut)
	{
		//TODO
		return 0;
	}

	int FixAtan2::checkArchitecture(int archType)
	{
		int k = (wIn/4);

		double error, errorMax, errorMin;
		bool errorSatisfied;
		double errorLimit = 1.0/(1 << (wOut-1));

		cout << "Beginning computations for error checking method based on the "
				<< (archType==0 ? "plane's equation" : (archType==1 ? "Taylor polynomial of order 1" : "Taylor polynomial of order 2"))
				<< " Atan2, with wIn=" << wIn << " and wOut=" << wOut << endl;

		cout << tab << "Computing the parameters for the a, b and d tables" << endl;
		errorSatisfied = false;
		while(!errorSatisfied)
		{
			bool runLoop = true;
			cout << tab << tab << "trying k=" << k << endl;

			errorMax = 0.0;
			errorMin = 1.0;

			for(int j=0; (j<(1<<k)-1 && runLoop); j++)
			{
				for(int i=(1<<(k-1)); (i<((1<<k)-1) && runLoop); i++)
				{
					Point *P1, *P2, *P3, *P4, *P5;
					Plane *plane;
					double a, b, c, d, e, f;
					double distance, distanceMin;
					double valI, valJ, centerValI, centerValJ, increment, halfIncrement;
					double distanceFromCenter;

					//create the actual values for the points in which we compute the atan2
					valI = 1.0*i/(1<<k);
					valJ = 1.0*j/(1<<k);
					increment = 1.0/(1<<k);
					halfIncrement = increment/2.0;
					centerValI = valI + halfIncrement;
					centerValJ = valJ + halfIncrement;
					distanceFromCenter = 1.0/(1<<(wIn+1));

					/**
					 * version 1: compute the values stored in the tables using the equation of the plane
					 */
					if(archType == 0)
					{
						//create the points
						P1 = new Point(centerValI,						centerValJ,						atan2(centerValJ, 						centerValI						));
						P2 = new Point(centerValI,						centerValJ+distanceFromCenter,	atan2(centerValJ+distanceFromCenter, 	centerValI						));
						P3 = new Point(centerValI+distanceFromCenter,	centerValJ, 					atan2(centerValJ, 						centerValI+distanceFromCenter	));
						P4 = new Point(centerValI, 						centerValJ-distanceFromCenter, 	atan2(centerValJ-distanceFromCenter, 	centerValI			  			));
						P5 = new Point(centerValI-distanceFromCenter,	centerValJ, 					atan2(centerValJ, 						centerValI-distanceFromCenter	));

						//create all the planes given by all the possible combinations of the 4 points
						// and determine the minimum distance from the plane given by three of the points
						// to the remaining point
						//the plane given by P2, P3 and P4
						plane = new Plane(P2, P3, P4);
						distance = P1->getZ() - ((-P1->getX()*plane->getA() - P1->getY()*plane->getB() - plane->getD())/plane->getC());
						if(distance<0)
							distance *= -1;
						distanceMin = distance;
						a = plane->getA();
						b = plane->getB();
						c = plane->getC();
						d = plane->getD();
						delete plane;

						//the plane given by P2, P3 and P5
						plane = new Plane(P2, P3, P5);
						distance = P1->getZ() - ((-P1->getX()*plane->getA() - P1->getY()*plane->getB() - plane->getD())/plane->getC());
						if(distance<0)
							distance *= -1;
						if(distance<distanceMin)
						{
							distanceMin = distance;
							a = plane->getA();
							b = plane->getB();
							c = plane->getC();
							d = plane->getD();
						}
						delete plane;

						//the plane given by P2, P4 and P5
						plane = new Plane(P2, P4, P5);
						distance = P1->getZ() - ((-P1->getX()*plane->getA() - P1->getY()*plane->getB() - plane->getD())/plane->getC());
						if(distance<0)
							distance *= -1;
						if(distance<distanceMin)
						{
							distanceMin = distance;
							a = plane->getA();
							b = plane->getB();
							c = plane->getC();
							d = plane->getD();
						}
						delete plane;

						//the plane given by P3, P4 and P5
						plane = new Plane(P3, P4, P5);
						distance = P1->getZ() - ((-P1->getX()*plane->getA() - P1->getY()*plane->getB() - plane->getD())/plane->getC());
						if(distance<0)
							distance *= -1;
						if(distance<distanceMin)
						{
							distanceMin = distance;
							a = plane->getA();
							b = plane->getB();
							c = plane->getC();
							d = plane->getD();
						}
						delete plane;

						//update the variables for use in the error measurement
						a = -(a/c);
						b = -(b/c);
						d = -(d/c);

						delete P1;
						delete P2;
						delete P3;
						delete P4;
						delete P5;

						//re-adjust the plane coordinates so as to minimize
						//	the error at the corners and in the center
						double errorP1, errorP2, errorP3, errorP4;
						double errorCenter;
						double errorMaxx, errorMinn, errorMean, errorAdjust;

						errorCenter = atan2(centerValJ, centerValI) - (a*centerValI + b*centerValJ + d);
						errorP1 = atan2(valJ, valI) - (a*valI + b*valJ + d);
						errorP2 = atan2(valJ, valI+increment) - (a*(valI+increment) + b*valJ + d);
						errorP3 = atan2(valJ+increment, valI) - (a*valI + b*(valJ+increment) + d);
						errorP4 = atan2(valJ+increment, valI+increment) - (a*(valI+increment) + b*(valJ+increment) + d);

						errorMaxx = max(5, errorCenter, errorP1, errorP2, errorP3, errorP4);
						errorMinn = min(5, errorCenter, errorP1, errorP2, errorP3, errorP4);

						errorMean = (fabs(errorMaxx)+fabs(errorMinn))/2.0;
						if((errorMaxx>=0) && (errorMinn>=0))
						{
							errorAdjust = errorMinn+(errorMaxx-errorMinn)/2.0;
						}else if((errorMaxx>=0) && (errorMinn<0))
						{
							errorAdjust = errorMaxx-errorMean;
						}else if((errorMaxx<0) && (errorMinn<0))
						{
							errorAdjust = errorMaxx+(errorMinn-errorMaxx)/2.0;
						}

						//equalize the maximum and the minimum error
						errorCenter -= errorAdjust;
						errorP1 -= errorAdjust;
						errorP2 -= errorAdjust;
						errorP3 -= errorAdjust;
						errorP4 -= errorAdjust;

						//test if any of the errors are above the error limit
						if((errorCenter>errorLimit) || (errorP1>errorLimit) || (errorP2>errorLimit) || (errorP2>errorLimit) || (errorP2>errorLimit))
						{
							cout << tab << tab << tab << "not worth continuing: distance from corners after adjustments is greater than the error limit" << endl;
							cout << tab << tab << tab << tab << "distance from center=" << errorCenter << " P1=" << errorP1 << " P2=" << errorP2
									<< " P3=" << errorP3 << " P4=" << errorP4 << " error limit=" << errorLimit << endl;

							runLoop = false;
							errorMax = errorMaxx-errorAdjust;
							errorMin = errorMinn-errorAdjust;
							break;
						}

						//compute the new equation of the plane
						P1 = new Point(valI, valJ, (a*valI + b*valJ + d)-errorAdjust);
						P2 = new Point(valI+increment, valJ, (a*(valI+increment) + b*valJ + d)-errorAdjust);
						P3 = new Point(valI, valJ+increment, (a*valI + b*(valJ+increment) + d)-errorAdjust);
						plane = new Plane(P1, P2, P3);
						a = plane->getA();
						b = plane->getB();
						c = plane->getC();
						d = plane->getD();
						a = -(a/c);
						b = -(b/c);
						d = -(d/c);
						delete plane;
						delete P1;
						delete P2;
						delete P3;
					}

					/**
					 * version 2: compute the values stored in the tables using a Taylor polynomial (of order 1)
					 * 	then atan(y/x) = a*x + b*y + c
					 */
					else if(archType == 1)
					{
						a = -(1.0*centerValJ)/(centerValI*centerValI+centerValJ*centerValJ);
						b = (1.0*centerValI)/(centerValI*centerValI+centerValJ*centerValJ);
						c = atan2(centerValJ, centerValI);

						//test if any of the errors are above the error limit
						double errorCenter 	= (a*centerValI + b*centerValJ + c) - atan2(centerValJ, centerValI);
						if(errorCenter<0)
							errorCenter = -errorCenter;
						double errorP1 		= (a*valI + b*valJ + c) - atan2(valJ, valI);
						if(errorP1<0)
							errorP1 = -errorP1;
						double errorP2 		= (a*(valI+increment) + b*valJ + c) - atan2(valJ, (valI+increment));
						if(errorP2<0)
							errorP2 = -errorP2;
						double errorP3 		= (a*valI + b*(valJ+increment) + c) - atan2((valJ+increment), valI);
						if(errorP3<0)
							errorP3 = -errorP3;
						double errorP4 		= (a*(valI+increment) + b*(valJ+increment) + c) - atan2((valJ+increment), (valI+increment));
						if(errorP4<0)
							errorP4 = -errorP4;
						if((errorCenter>errorLimit) || (errorP1>errorLimit) || (errorP2>errorLimit)
								|| (errorP3>errorLimit) || (errorP4>errorLimit))
						{
							cout << tab << tab << tab << "not worth continuing: distance from corners is greater than the error limit" << endl;
							cout << tab << tab << tab << tab << "distance from center=" << errorCenter << " P1=" << errorP1 << " P2=" << errorP2
									<< " P3=" << errorP3 << " P4=" << errorP4 << " error limit=" << errorLimit << endl;

							runLoop = false;
							errorMax = max(4, errorP1, errorP2, errorP3, errorP4);
							errorMin = min(5, errorP1, errorP2, errorP3, errorP4);
							break;
						}
					}


					/**
					 * version 3: compute the values stored in the tables using a Taylor polynomial (of order 2)
					 * 	then atan(x/y) = a*x + b*y + c + d*x^2 + e*y^2 + f*x*y
					 */
					else if(archType == 2)
					{
						double denominator = centerValI*centerValI + centerValJ*centerValJ;
						double denominatorSqr = denominator*denominator;

						a = -(2.0*centerValJ)/(denominator);
						b = (2.0*centerValI)/(denominator);
						c = atan2(centerValJ, centerValI);
						d = (1.0*centerValI*centerValJ)/(denominatorSqr);
						e = -(1.0*centerValI*centerValJ)/(denominatorSqr);
						f = (1.0*centerValJ*centerValJ-1.0*centerValI*centerValI)/(denominatorSqr);

						//test if any of the errors are above the error limit
						double errorCenter 	= (a*centerValI + b*centerValJ + c + d*centerValI*centerValI + e*centerValJ*centerValJ + f*centerValI*centerValJ)
																- atan2(centerValJ, centerValI);
						if(errorCenter<0)
							errorCenter = -errorCenter;
						double errorP1 		= (a*valI + b*valJ + c + d*valI*valI + e*valJ*valJ + f*valI*valJ)
																- atan2(valJ, valI);
						if(errorP1<0)
							errorP1 = -errorP1;
						double errorP2 		= (a*(valI+increment) + b*valJ + c + d*(valI+increment)*(valI+increment) + e*valJ*valJ + f*(valI+increment)*valJ)
																- atan2(valJ, (valI+increment));
						if(errorP2<0)
							errorP2 = -errorP2;
						double errorP3 		= (a*valI + b*(valJ+increment) + c + d*valI*valI + e*(valJ+increment)*(valJ+increment) + f*valI*(valJ+increment))
																- atan2((valJ+increment), valI);
						if(errorP3<0)
							errorP3 = -errorP3;
						double errorP4 		= (a*(valI+increment) + b*(valJ+increment) + c
								+ d*(valI+increment)*(valI+increment) + e*(valJ+increment)*(valJ+increment) + f*(valI+increment)*(valJ+increment))
																- atan2((valJ+increment), (valI+increment));
						if(errorP4<0)
							errorP4 = -errorP4;
						if((errorCenter>errorLimit) || (errorP1>errorLimit) || (errorP2>errorLimit)
								|| (errorP3>errorLimit) || (errorP4>errorLimit))
						{
							cout << tab << tab << tab << "not worth continuing: distance from corners is greater than the error limit" << endl;
							cout << tab << tab << tab << tab << "distance from center=" << errorCenter << " P1=" << errorP1 << " P2=" << errorP2
									<< " P3=" << errorP3 << " P4=" << errorP4 << " error limit=" << errorLimit << endl;

							runLoop = false;
							errorMax = max(4, errorP1, errorP2, errorP3, errorP4);
							errorMin = min(4, errorP1, errorP2, errorP3, errorP4);
							break;
						}
					}
					else
					{
						THROWERROR("Error: unknown type of architecture to check");
					}

					//determine the maximum values of A, B and D (or C, depending on the architecture)
					if((archType == 0) || (archType == 1) || (archType == 2))
					{
						double auxA, auxB, auxC, auxD, auxE, auxF;

						// maximum value of A/(Pi/2)
						auxA = ((a < 0) ? -a : a);
						auxA = auxA / (M_PI/2.0);
						if(auxA > maxValA)
							maxValA = auxA;
						// maximum value of B
						auxB = ((b < 0) ? -b : b);
						auxB = auxB / (M_PI/2.0);
						if(auxB > maxValB)
							maxValB = auxB;
						if(archType == 0)
						{
							// maximum value of D
							auxD = d + a*valI + b*valJ;
							auxD = ((auxD < 0) ? -auxD : auxD);
							auxD = auxD / (M_PI/2.0);
							if(auxD > maxValC)
								maxValC = auxD;
						}else if(archType == 1)
						{
							// maximum value of C
							auxC = c + a*valI + b*valJ;
							auxC = ((auxC < 0) ? -auxC : auxC);
							auxC = auxC / (M_PI/2.0);
							if(auxC > maxValC)
								maxValC = auxC;
						}else if(archType == 2)
						{
							// maximum value of C
							auxC = ((c < 0) ? -c : c);
							auxC = auxC / (M_PI/2.0);
							if(auxC > maxValC)
								maxValC = auxC;
							// maximum value of D
							auxD = ((d < 0) ? -d : d);
							auxD = auxD / (M_PI/2.0);
							if(auxD > maxValD)
								maxValD = auxD;
							// maximum value of E
							auxE = ((e < 0) ? -e : e);
							auxE = auxE / (M_PI/2.0);
							if(auxE > maxValE)
								maxValE = auxE;
							// maximum value of F
							auxF = ((f < 0) ? -f : f);
							auxF = auxF / (M_PI/2.0);
							if(auxF > maxValF)
								maxValF = auxF;
						}
					}


					//now check the error against all the points in the plane, at the given resolution
					for(int n=0; (n<(1<<(wIn-k)) && runLoop); n++)
						for(int m=0; (m<(1<<(wIn-k)) && runLoop); m++)
						{
							double valIPrime, valJPrime;

							//create the actual values for the atan2
							valIPrime = 1.0*((i<<(wIn-k))+n)/(1<<wIn);
							valJPrime = 1.0*((j<<(wIn-k))+m)/(1<<wIn);

							double referenceValue = atan2(valJPrime, valIPrime);
							double computedValue;

							if(archType == 0){
								// v1
								computedValue  = a*valIPrime + b*valJPrime + d;
							}else if(archType == 1){
								// v2
								computedValue  = a*valIPrime + b*valJPrime + c;
							}else{
								// v3
								computedValue  = a*valIPrime + b*valJPrime + c + d*valIPrime*valIPrime + e*valJPrime*valJPrime + f*valIPrime*valJPrime;
							}

							error = referenceValue-computedValue;
							if(error < 0)
								error *= -1;

							if(error > errorMax)
								errorMax = error;
							if((error < errorMin) && !((n==0)&&(m==0)))
								errorMin = error;

							if(errorMax > errorLimit)
							{
								cout << tab << tab << tab << "not worth continuing: error=" << errorMax << " error limit=" << errorLimit << endl;
								cout << tab << tab << tab << tab << "i=" << i << " j=" << j << " n=" << n << " m=" << m << endl;

								runLoop = false;
								break;
							}
						}
				}

				cout << "\r                                                                                                                   \r";
				cout <<  tab << tab << tab << std::setprecision(4) << ((1.0*j/((1<<k)-1))*100) << "% done"
						<< tab << "current maximum error=" << errorMax << " current minimum error=" << errorMin << " error limit=" << errorLimit;
				cout.flush();
			}

			cout << endl << tab << tab << tab << "computed the maximum error" << endl;

			//check to see if the current precision k is enough
			//	if not, increase k and compute the maximum error again
			if(errorMax < errorLimit)
			{
				errorSatisfied = true;
			}else
			{
				cout << tab << tab << tab << "error limit of " << errorLimit << " not satisfied at k=" << k
						<< ", with wIn=" << wIn << " and wOut=" << wOut << endl;
				cout << tab << tab << tab << tab << "actual error errorMax=" << errorMax << endl;

				k++;
			}
		}

		cout << "Error limit of " << (1.0/(1 << (wOut-1))) << " satisfied at k=" << k
				<< ", with wIn=" << wIn << " and wOut=" << wOut << endl;
		if((archType == 0) || (archType == 1) || (archType == 2))
		{
			cout << tab << "Computed the maximum values of the table parameters (divided by Pi/2): maxValA="
					<< maxValA << " maxValB=" << maxValB;
			if(archType == 0)
				cout << " maxValD=" << maxValC << endl;
			else if(archType == 1)
				cout << " maxValC=" << maxValC << endl;
			else
				cout << " maxValC=" << maxValC << " maxValD=" << maxValD
					<< " maxValE=" << maxValE << " maxValF=" << maxValF << endl;
		}

		return k;
	}

//	void FixAtan2::emulate(TestCase * tc)
//	{
//		mpfr_t x,y,a, constPi;
//		mpfr_init2(x, 10*wIn);
//		mpfr_init2(y, 10*wIn);
//		mpfr_init2(a, 10*wIn);
//		mpfr_init2(constPi, 10*wIn);
//		mpfr_const_pi( constPi, GMP_RNDN);
//
//		mpz_class az;
//
//		/* Get I/O values */
//		mpz_class svX = tc->getInputValue("X");
//		mpz_class svY = tc->getInputValue("Y");
//
//		// interpret as signed two'ss complement
//		if (1==(svX >> (wIn-1))) // sign bit
//			svX -= (1<<wIn);
//		if (1==(svY >> (wIn-1))) // sign bit
//			svY -= (1<<wIn);
//		/* Compute correct value */
//
//		mpfr_set_z (x, svX.get_mpz_t(), GMP_RNDN); //  exact
//		mpfr_set_z (y, svY.get_mpz_t(), GMP_RNDN); //  exact
//
//		mpfr_atan2(a, y, x, GMP_RNDN); // a between -pi and pi
//		mpfr_div(a, a, constPi, GMP_RNDN); // a between -1 and 1
//
//		// Now convert a to fix point
//		// Align to fix point by adding 6 -- if we just add 4 there is a 1-bit shift in case a<0
//		mpfr_add_d(a, a, 6.0, GMP_RNDN);
//		mpfr_mul_2si (a, a, wOut-1, GMP_RNDN); // exact scaling
//
//		mpz_class mask = (mpz_class(1)<<wOut) -1;
//
//		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDD); // there can be a real rounding here
//		az -= mpz_class(6)<<(wOut-1);
//		az &= mask;
//		tc->addExpectedOutput ("A", az);
//
//		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDU); // there can be a real rounding here
//		az -= mpz_class(6)<<(wOut-1);
//		az &= mask;
//		tc->addExpectedOutput ("A", az);
//
//		// clean up
//		mpfr_clears (x,y,a, constPi, NULL);
//	}


	void FixAtan2::emulate(TestCase* tc)
	{
		mpfr_t mpX, mpY, mpR, pi_mpfr;
		mpz_class svRu, svRd, aux, temp;

		/// Get I/O values
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		//get the true value of X
		/*
		temp = mpz_class(1);
		temp = temp << (wIn-2);
		svX = svX + temp;
		*/

		// interpret as signed two's complement
		if(1 == (svX >> (wIn-1))) // sign bit
			svX -= (1<<wIn);
		if(1 == (svY >> (wIn-1))) // sign bit
			svY -= (1<<wIn);

		mpfr_inits2(10000, mpX, mpY, mpR, pi_mpfr, (mpfr_ptr)0);

		// Cast X and Y to mpfr
		mpfr_set_z(mpX, svX.get_mpz_t(), GMP_RNDN);
		mpfr_set_z(mpY, svY.get_mpz_t(), GMP_RNDN);

		// scale appropriately, on the input format: multiply by 2^(-wIn+1)
		mpfr_div_2si(mpX, mpX, wIn-1, GMP_RNDN);
		mpfr_div_2si(mpY, mpY, wIn-1, GMP_RNDN);

		// do the computations
		mpfr_atan2(mpR, mpY, mpX, GMP_RNDN);

		//divide by Pi/2
		mpfr_const_pi(pi_mpfr, GMP_RNDN);
		mpfr_div_si(pi_mpfr, pi_mpfr, 2, GMP_RNDN);
		mpfr_div(mpR, mpR, pi_mpfr, GMP_RNDN);

		// scale back to an integer, on the output format: multiply by 2^(wOut-1)
		mpfr_mul_2si(mpR, mpR, wOut-1, GMP_RNDN);

		// extract the result
		//	rounded down
		mpfr_get_z(svRd.get_mpz_t(), mpR, GMP_RNDD);
		if(mpz_class(svRd) < 0)
		{
			aux = mpz_class(1);
			aux = aux << (wOut);
			svRd = aux + svRd;
			cout << svRd << " ";
		}
		tc->addExpectedOutput("R", svRd);
		//	rounded up
		mpfr_get_z(svRu.get_mpz_t(), mpR, GMP_RNDU);
		if(mpz_class(svRu) < 0)
		{
			aux = mpz_class(1);
			aux = aux << (wOut);
			svRu = aux + svRu;
			cout << " " << svRu << endl;
		}
		tc->addExpectedOutput("R", svRu);



		//----------------------------------------------------------------------
		//-------- Generate the signals of the architecture using software
		/*
		mpfr_t x, y, deltax, deltay, deltax2, deltay2, deltaxDeltay,
				aDeltax, bDeltay, dDeltax2, eDeltay2, fDeltaxDeltay,
				sum1, sum2, sum3, sum4, sum5, result,
				tmp;
		mpfr_t a, b, c, d, e, f;
		mpz_class tmpMpz;

		mpfr_inits2(10000, x, y, deltax, deltay, deltax2, deltay2, deltaxDeltay, aDeltax, bDeltay, dDeltax2, eDeltay2, fDeltaxDeltay, sum1, sum2, sum3, sum4, sum5, result, tmp, (mpfr_ptr)0);
		mpfr_inits2(10000, a, b, c, d, e, f, (mpfr_ptr)0);

		//generate the parameters
		generateTaylorOrder2Parameters(svX.get_si(), svY.get_si(), a, b, c, d, e, f);

		//create x and y
		mpfr_set(x, mpX, GMP_RNDN);
		mpfr_set(y, mpY, GMP_RNDN);

		//create deltax and deltay
		//	create deltax
		tmpMpz = mpz_class(svX);
		tmpMpz = tmpMpz >> (wIn-kSize);
		tmpMpz = tmpMpz << (wIn-kSize);
		tmpMpz = svX - tmpMpz;
		mpfr_set_z(deltax, tmpMpz.get_mpz_t(), GMP_RNDN);
		mpfr_div_2si(deltax, deltax, wIn, GMP_RNDN);
		mpfr_set_si(tmp, 1, GMP_RNDN);
		mpfr_div_2si(tmp, tmp, kSize+1, GMP_RNDN);
		mpfr_sub(deltax, deltax, tmp, GMP_RNDN);
		//	create deltay
		tmpMpz = mpz_class(svY);
		tmpMpz = tmpMpz >> (wIn-kSize);
		tmpMpz = tmpMpz << (wIn-kSize);
		tmpMpz = svY - tmpMpz;
		mpfr_set_z(deltay, tmpMpz.get_mpz_t(), GMP_RNDN);
		mpfr_div_2si(deltay, deltay, wIn, GMP_RNDN);
		mpfr_set_si(tmp, 1, GMP_RNDN);
		mpfr_div_2si(tmp, tmp, kSize+1, GMP_RNDN);
		mpfr_sub(deltay, deltay, tmp, GMP_RNDN);

		//create deltax2 and deltay2 and deltaxDeltay
		mpfr_sqr(deltax2, deltax, GMP_RNDN);
		mpfr_sqr(deltay2, deltay, GMP_RNDN);
		mpfr_mul(deltaxDeltay, deltax, deltay, GMP_RNDN);

		//create aDeltax, bDeltay
		mpfr_mul(aDeltax, a, deltax, GMP_RNDN);
		mpfr_mul(bDeltay, b, deltay, GMP_RNDN);

		//create dDeltax2, eDeltay2, fDeltaxDeltay
		mpfr_mul(dDeltax2, d, deltax2, GMP_RNDN);
		mpfr_mul(eDeltay2, e, deltay2, GMP_RNDN);
		mpfr_mul(fDeltaxDeltay, f, deltaxDeltay, GMP_RNDN);

		//create the sums
		mpfr_add(sum1, aDeltax, bDeltay, GMP_RNDN);
		mpfr_add(sum2, c, dDeltax2, GMP_RNDN);
		mpfr_add(sum3, eDeltay2, fDeltaxDeltay, GMP_RNDN);
		mpfr_add(sum4, sum1, sum2, GMP_RNDN);
		mpfr_add(sum5, sum3, sum4, GMP_RNDN);

		//extract the values out of the mpfr variables
		double x_d, y_d, deltax_d, deltay_d, deltax2_d, deltay2_d, deltaxDeltay_d,
				aDeltax_d, bDeltay_d, dDeltax2_d, eDeltay2_d, fDeltaxDeltay_d,
				sum1_d, sum2_d, sum3_d, sum4_d, sum5_d;

		x_d = mpfr_get_d(x, GMP_RNDN);
		y_d = mpfr_get_d(y, GMP_RNDN);
		deltax_d = mpfr_get_d(deltax, GMP_RNDN);
		deltay_d = mpfr_get_d(deltay, GMP_RNDN);
		deltax2_d = mpfr_get_d(deltax2, GMP_RNDN);
		deltay2_d = mpfr_get_d(deltay2, GMP_RNDN);
		deltaxDeltay_d = mpfr_get_d(deltaxDeltay, GMP_RNDN);
		aDeltax_d = mpfr_get_d(aDeltax, GMP_RNDN);
		bDeltay_d = mpfr_get_d(bDeltay, GMP_RNDN);
		dDeltax2_d = mpfr_get_d(dDeltax2, GMP_RNDN);
		eDeltay2_d = mpfr_get_d(eDeltay2, GMP_RNDN);
		fDeltaxDeltay_d = mpfr_get_d(fDeltaxDeltay, GMP_RNDN);
		sum1_d = mpfr_get_d(sum1, GMP_RNDN);
		sum2_d = mpfr_get_d(sum2, GMP_RNDN);
		sum3_d = mpfr_get_d(sum3, GMP_RNDN);
		sum4_d = mpfr_get_d(sum4, GMP_RNDN);
		sum5_d = mpfr_get_d(sum5, GMP_RNDN);


		mpfr_clears(x, y, deltax, deltay, deltax2, deltay2, deltaxDeltay, aDeltax, bDeltay, dDeltax2, eDeltay2, fDeltaxDeltay, sum1, sum2, sum3, sum4, sum5, result, tmp, (mpfr_ptr)0);
		mpfr_clears(a, b, c, d, e, f, (mpfr_ptr)0);
		*/
		//----------------------------------------------------------------------


		// clean up
		mpfr_clears(mpX, mpY, mpR, pi_mpfr, (mpfr_ptr)0);
	}

	void FixAtan2::generateTaylorOrder2Parameters(int x, int y, mpfr_t &fa, mpfr_t &fb, mpfr_t &fc, mpfr_t &fd, mpfr_t &fe, mpfr_t &ff)
	{
		mpfr_t centerValI, centerValJ, increment, temp, tempSqr, temp2;
		int k = kSize;

		mpfr_inits2(10000, centerValI, centerValJ, increment, temp, tempSqr, temp2, (mpfr_ptr)0);

		mpfr_set_si(increment, 1, GMP_RNDN);
		mpfr_div_2si(increment, increment, k+1, GMP_RNDN);
		mpfr_set_si(centerValI, x, GMP_RNDN);
		mpfr_div_2si(centerValI, centerValI, k, GMP_RNDN);
		mpfr_add(centerValI, centerValI, increment, GMP_RNDN);
		mpfr_set_si(centerValJ, y, GMP_RNDN);
		mpfr_div_2si(centerValJ, centerValJ, k, GMP_RNDN);
		mpfr_add(centerValJ, centerValJ, increment, GMP_RNDN);

		mpfr_sqr(temp, centerValI, GMP_RNDN);
		mpfr_sqr(temp2, centerValJ, GMP_RNDN);
		mpfr_add(temp, temp, temp2, GMP_RNDN);
		mpfr_sqr(tempSqr, temp, GMP_RNDN);
		//create A
		mpfr_set(fa, centerValJ, GMP_RNDN);
		mpfr_div(fa, fa, temp, GMP_RNDN);
		mpfr_neg(fa, fa, GMP_RNDN);
		//create B
		mpfr_set(fb, centerValI, GMP_RNDN);
		mpfr_div(fb, fb, temp, GMP_RNDN);
		//create C
		mpfr_atan2(fc, centerValJ, centerValI, GMP_RNDN);
		//create D
		mpfr_set(fd, centerValI, GMP_RNDN);
		mpfr_mul(fd, fd, centerValJ, GMP_RNDN);
		mpfr_div(fd, fd, tempSqr, GMP_RNDN);
		//create E
		mpfr_set(fe, fd, GMP_RNDN);
		mpfr_neg(fe, fe, GMP_RNDN);
		//create F
		mpfr_sqr(ff, centerValJ, GMP_RNDN);
		mpfr_sqr(temp2, centerValI, GMP_RNDN);
		mpfr_sub(ff, ff, temp2, GMP_RNDN);
		mpfr_div(ff, ff, tempSqr, GMP_RNDN);

		//clean-up
		mpfr_clears(centerValI, centerValJ, increment, temp, tempSqr, temp2, (mpfr_ptr)0);
	}



	void FixAtan2::buildStandardTestCases(TestCaseList* tcl)
	{
		//TODO
	}




}
