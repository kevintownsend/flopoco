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
		plainVHDL = false;

		// build the name
		ostringstream name;
		name <<"FixAtan2_" << vhdlize(wIn) << "_" << vhdlize(wOut) << "_archType_" << vhdlize(architectureType_);
		setName(name.str());

		//set the libraries to use
		useNumericStd();

		//initialize global variables
		maxValA = -1;		//this is the absolute value, so should be always non-negative, once used
		maxValB = -1;
		maxValC = -1;

		//declare the inputs and the outputs
		addInput ("X",  wIn-1);
		addInput ("Y",  wIn);
		addOutput("R",  wOut, 2 /*number of possible output values*/);

		//build the architecture
		if(architectureType == 0)
		{
			//based on the plane's equation

			//check if the selected type of architecture can achieve the required accuracy
			//	also, determine the size of the parts of the input
			//	signals to use as an address for the table
			int k;
			int guardBitsSum, guardBitsApprox;

			guardBitsSum = 2;				// the error when computing Ax+By+C is at most 2.5 ulps, so we'll need 2 extra guard bits
			guardBitsApprox = 2;			// determined when computing the parameters for the tables
			g = guardBitsApprox + guardBitsSum;

			//set the ratio for the multiplications
			ratio = 0.95;

			k = checkArchitecture(0);

			//determine the size of the constants to store in the table
			// maxValA, maxValB and maxValC should have been set by the call to checkArchitecture()
			int msbA = intlog2(maxValA)+1;
			int msbB = intlog2(maxValB)+1;
			int msbC = intlog2(maxValC)+1;
			int maxMSB = max(3, msbA, msbB, msbC);

			//create the bitheap
			bitHeap = new BitHeap(this, (maxMSB+2)+wOut+g);

			if(plainVHDL)
			{
				//split the input signals, and create the address signal for tha table
				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(X" << range(2*k-2, k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(Y" << range(2*k-1, k) << ");" << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;

				//create the table for atan(y/x)
				Atan2Table *table = new Atan2Table(target, 2*k-1, msbA+msbB+msbC+3*(wOut+g), 0, msbA, msbB, msbC);

				//add the table to the operator
				addSubComponent(table);
				useSoftRAM(table);
				//useHardRAM(table);
				inPortMap (table , "X", "atan2TableInput");
				outPortMap(table , "Y", "atan2TableOutput");
				vhdl << instance(table , "KCMTable");

				//split the output of the table to the corresponding parts A, B and C
				vhdl << tab << declareFixPoint("C", true, msbC-1,  -wOut-g) << " <= signed(atan2TableOutput"
						<< range(msbC+wOut+g-1, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("B", true, msbB-1,  -wOut-g) << " <= signed(atan2TableOutput"
						<< range(msbB+wOut+g+msbC+wOut+g-1, msbC+wOut+g) << ");" << endl;
				vhdl << tab << declareFixPoint("A", true, msbA-1,  -wOut-g) << " <= signed(atan2TableOutput"
						<< range(msbA+wOut+g+msbB+wOut+g+msbC+wOut+g-1, msbB+wOut+g+msbC+wOut+g) << ");" << endl;

				vhdl << tab << declareFixPoint("XLow", true, -k,  -2*k) << " <= signed('0' & X" << range(k-1, 0) << ");" << endl;
				vhdl << tab << declareFixPoint("YLow", true, -k,  -2*k) << " <= signed('0' & Y" << range(k-1, 0) << ");" << endl;

				//create A*X_low and B*Y_low
				vhdl << tab << declareFixPoint("AXLow", true, msbA-k,  -wOut-g-2*k) << " <= A * XLow;" << endl;
				vhdl << tab << declareFixPoint("BYLow", true, msbB-k,  -wOut-g-2*k) << " <= B * YLow;" << endl;
				//align A*X_low and B*Y_low to the output format
				resizeFixPoint("AXLow_sgnExtended", "AXLow", maxMSB-1, -wOut-g);
				resizeFixPoint("BYLow_sgnExtended", "BYLow", maxMSB-1, -wOut-g);

				//add everything up
				vhdl << tab << declareFixPoint("AXLowAddBYLow", true, maxMSB,  -wOut-g)
						<< " <= (AXLow_sgnExtended(AXLow_sgnExtended'HIGH) & AXLow_sgnExtended) + (BYLow_sgnExtended(BYLow_sgnExtended'HIGH) & BYLow_sgnExtended);" << endl;
				resizeFixPoint("C_sgnExtended", "C", maxMSB+1, -wOut-g);
				vhdl << tab << declareFixPoint("AXLowAddBYLowAddC", true, maxMSB+1,  -wOut-g) << " <= (AXLowAddBYLow(AXLowAddBYLow'HIGH) & AXLowAddBYLow) + C_sgnExtended;" << endl;

				//extract the final result
				resizeFixPoint("Rint", "AXLowAddBYLowAddC", 1, -wOut+1);
				vhdl << tab << declareFixPoint("Rint_rndCst", true, 1, -wOut+1) << " <= signed(std_logic_vector\'(\"" << zg(wOut, -2) << "1\"));" << endl;
				vhdl << tab << declareFixPoint("Rint_rnd", true, 1, -wOut+1) << " <= Rint + Rint_rndCst;" << endl;

				//return the result
				resizeFixPoint("Rint_stdlv", "Rint_rnd", 1, -wOut+2);
				vhdl << tab << "R <= std_logic_vector(Rint_stdlv);" << endl;

			}else
			{
				//create the input signals for the table
				vhdl << tab << declare("XHigh", k-1) << " <= std_logic_vector(X" << range(2*k-2, k) << ");" << endl;
				vhdl << tab << declare("YHigh", k)   << " <= std_logic_vector(Y" << range(2*k-1, k) << ");" << endl;
				vhdl << tab << declare("atan2TableInput", 2*k-1) << " <= std_logic_vector(XHigh) & std_logic_vector(YHigh);" << endl;

				//create the table for atan(y/x)
				Atan2Table *table = new Atan2Table(target, 2*k-1, msbA+msbB+msbC+3*(wOut+g), 0, msbA, msbB, msbC);

				//add the table to the operator
				addSubComponent(table);
				useSoftRAM(table);
				//useHardRAM(table);
				inPortMap (table , "X", "atan2TableInput");
				outPortMap(table , "Y", "atan2TableOutput");
				vhdl << instance(table , "KCMTable");

				//extract signals A, B and C
				vhdl << tab << declare("C", msbC+wOut+g) << " <= atan2TableOutput" << range(msbC+wOut+g-1, 0) << ";" << endl;
				vhdl << tab << declare("B", msbB+wOut+g) << " <= atan2TableOutput" << range(msbB+wOut+g+msbC+wOut+g-1, msbC+wOut+g) << ";" << endl;
				vhdl << tab << declare("A", msbA+wOut+g) << " <= atan2TableOutput" << range(msbA+wOut+g+msbB+wOut+g+msbC+wOut+g-1, msbB+wOut+g+msbC+wOut+g) << ";" << endl;

				//create Ax and By
				vhdl << tab << declare("XLow", k+1) << " <= '0' & X" << range(k-1, 0) << ";" << endl;
				vhdl << tab << declare("YLow", k+1) << " <= '0' & Y" << range(k-1, 0) << ";" << endl;

				IntMultiplier* multAx;
				multAx = new IntMultiplier(this,								//parent operator
											 bitHeap,							//the bit heap that performs the compression
											 getSignalByName("XLow"),			//first input to the multiplier (a signal)
											 getSignalByName("A"),				//second input to the multiplier (a signal)
											 -2*k,								//offset of the LSB of the multiplier in the bit heap
											 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
											 true,								//signed/unsigned operator
											 ratio);							//DSP ratio
				IntMultiplier* multBy;
				multBy = new IntMultiplier(this,								//parent operator
											 bitHeap,							//the bit heap that performs the compression
											 getSignalByName("YLow"),			//first input to the multiplier (a signal)
											 getSignalByName("B"),				//second input to the multiplier (a signal)
											 -2*k,								//offset of the LSB of the multiplier in the bit heap
											 false /*negate*/,					//whether to subtract the result of the multiplication from the bit heap
											 true,								//signed/unsigned operator
											 ratio);							//DSP ratio

				bitHeap->addSignedBitVector(0,									//weight of signal in the bit heap
											"C",								//name of the signal
											msbC+wOut+g,						//size of the signal added
											0,									//index of the lsb in the bit vector from which to add the bits of the addend
											false);								//if we are correcting the index in the bit vector with a negative weight

				//add the rounding bit - take into consideration the final alignment
				bitHeap->addConstantOneBit(g-1+2);

				//compress the bit heap
				bitHeap -> generateCompressorVHDL();

				//extract the result - take into consideration the final alignment
				vhdl << tab << "R <= " << bitHeap->getSumName() << range(wOut+g-1+2, g+2) << ";" << endl;
			}

		}else if(architectureType == 1)
		{
			//based on an order 1 Taylor approximating polynomial
			THROWERROR("in FixAtan2 constructor: architecture not yet implemented");
		}else if(architectureType == 2)
		{
			//based on an order 2 Taylor approximating polynomial
			THROWERROR("in FixAtan2 constructor: architecture not yet implemented");
		}else if(architectureType == 3)
		{
			//based on rotating using a table (for sine and cosine)
			//	and then using the Taylor series for (1/x)
			THROWERROR("in FixAtan2 constructor: architecture not yet implemented");
		}else
		{
			THROWERROR("in FixAtan2 constructor: invalid value for the requested architecture type");
		}
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
		double errorLimit = 1.0/(1 << wOut);

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
						a = -(2.0*centerValJ)/(centerValI*centerValI+centerValJ*centerValJ);
						b = (2.0*centerValI)/(centerValI*centerValI+centerValJ*centerValJ);
						c = atan2(centerValJ, centerValI);
						d = (1.0*centerValI*centerValJ)
												/((centerValI*centerValI+centerValJ*centerValJ)*(centerValI*centerValI+centerValJ*centerValJ));
						e = -(1.0*centerValI*centerValJ)
												/((centerValI*centerValI+centerValJ*centerValJ)*(centerValI*centerValI+centerValJ*centerValJ));
						f = (1.0*centerValJ*centerValJ-1.0*centerValI*centerValI)
												/((centerValI*centerValI+centerValJ*centerValJ)*(centerValI*centerValI+centerValJ*centerValJ));

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
					if((archType == 0) || (archType == 1))
					{
						double auxA, auxB, auxD;

						auxA = ((a < 0) ? -a : a);
						if(auxA > maxValA)
							maxValA = auxA;
						auxB = ((b < 0) ? -b : b);
						if(auxB > maxValB)
							maxValB = auxB;
						auxD = d + a*valI + b*valJ;
						auxD = ((auxD < 0) ? -auxD : auxD);
						if(auxD > maxValC)
							maxValC = auxD;
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

		cout << "Error limit of " << (1 >> wOut) << " satisfied at k=" << k
				<< ", with wIn=" << wIn << " and wOut=" << wOut << endl;
		if((archType == 0) || (archType == 1))
		{
			cout << tab << "Computed the maximum values of the table parameters: maxValA="
					<< maxValA << " maxValB=" << maxValB << " maxValD=" << maxValC << endl;
		}

		return k;
	}


	void FixAtan2::emulate(TestCase* tc)
	{
		mpfr_t mpX, mpY, mpR;
		mpz_class svRu, svRd, aux, temp;

		/// Get I/O values
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		//get the true value of X
		temp = mpz_class(1);
		temp = temp << (wIn-1);
		svX = svX + temp;

		mpfr_inits2(10000, mpX, mpY, mpR, (mpfr_ptr)0);

		// Cast X and Y to mpfr
		mpfr_set_z(mpX, svX.get_mpz_t(), GMP_RNDN);
		mpfr_set_z(mpY, svY.get_mpz_t(), GMP_RNDN);

		// scale appropriately, on the input format: multiply by 2^(-wIn)
		mpfr_div_2si(mpX, mpX, wIn, GMP_RNDN);
		mpfr_div_2si(mpY, mpY, wIn, GMP_RNDN);

		// do the computations
		mpfr_atan2(mpR, mpY, mpX, GMP_RNDN);

		// scale back to an integer, on the output format: multiply by 2^(wOut-2)
		mpfr_mul_2si(mpR, mpR, wOut-2, GMP_RNDN);

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

		// clean up
		mpfr_clears(mpX, mpY, mpR, (mpfr_ptr)0);
	}



	void FixAtan2::buildStandardTestCases(TestCaseList* tcl)
	{
		//TODO
	}




}
