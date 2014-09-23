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
	FixAtan2::FixAtan2(Target* target_, int wIn_, int wOut_, map<string, double> inputDelays_):
		Operator(target_, inputDelays_),
		wIn(wIn_), wOut(wOut_)
	{

		srcFileName = "FixAtan2";
		setCopyrightString("Florent de Dinechin, Matei Istoan, 2014");

		int k = (wIn/4);//(wIn+1)/2;

		double error, errorMax, errorMin;
		bool errorSatisfied;
		double errorLimit = 1.0/(1 << wOut);

		cout << "Beginning computations for Atan2, with wIn=" << wIn << " and wOut=" << wOut << endl;

		cout << tab << "Computing the parameters for the a, b and d tables" << endl;
		errorSatisfied = false;
		while(!errorSatisfied)
		{
			bool runLoop = true;
			cout << tab << tab << "trying k=" << k << endl;

			errorMax = 0.0;
			errorMin = 1.0;

			for(int i=0; (i<(1<<k)-1 && runLoop); i++)
			{
				for(int j=(1<<(k-1)); (j<((1<<k)-1) && runLoop); j++)
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
					distanceFromCenter = 1.0/(1<<(wIn-1));

					/**
					 * version 1: compute the values stored in the tables using the equation of the plane
					 */
					/*
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
					*/


					/**
					 * version 2: compute the values stored in the tables using a Taylor polynomial (of order 1)
					 * 	then atan(y/x) = a*x + b*y + c
					 */
					/*
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
					*/


					/**
					 * version 3: compute the values stored in the tables using a Taylor polynomial (of order 2)
					 * 	then atan(x/y) = a*x + b*y + c + d*x^2 + e*y^2 + f*x*y
					 */
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


					//now check the error against all the points in the plane, at the given resolution
					for(int n=0; (n<(1<<(wIn-k)) && runLoop); n++)
						for(int m=0; (m<(1<<(wIn-k)) && runLoop); m++)
						{
							double valIPrime, valJPrime;

							//create the actual values for the atan2
							valIPrime = 1.0*((i<<(wIn-k))+n)/(1<<wIn);
							valJPrime = 1.0*((j<<(wIn-k))+m)/(1<<wIn);

							double referenceValue = atan2(valJPrime, valIPrime);
							// v1
							//double computedValue  = a*valIPrime + b*valJPrime + d;
							// v2
							//double computedValue  = a*valIPrime + b*valJPrime + c;
							// v3
							double computedValue  = a*valIPrime + b*valJPrime + c + d*valIPrime*valIPrime + e*valJPrime*valJPrime + f*valIPrime*valJPrime;

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
				cout <<  tab << tab << tab << std::setprecision(4) << ((1.0*i/((1<<k)-1))*100) << "% done"
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
	}


	FixAtan2::~FixAtan2()
	{

	}

	FixAtan2* FixAtan2::newComponentAndInstanceNumericStd(Operator* op,
															int wIn,
															int wOut
														)
	{
		return 0;
	}


	void FixAtan2::emulate ( TestCase* tc )
	{
		//TODO
	}



	void FixAtan2::buildStandardTestCases(TestCaseList* tcl)
	{
		//TODO
	}




}
