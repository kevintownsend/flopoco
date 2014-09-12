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

		/*
		//test the plane equation
		cout << "Testing the plane equations" << endl;

		Point *P1, *P2, *P3;
		Plane* plane;

		P1 = new Point(1, 2, 3);
		P2 = new Point(3, 9, 8);
		P3 = new Point(2, 4, 7);

		plane = new Plane(P1, P2, P3);

		cout << "Created the plane. Plane equation is: "
				"(" << plane->getA() << ")*x + (" << plane->getB() << ")*y + (" << plane->getC() << ")*z + (" << plane->getD() << ") = 0" << endl << endl;

		//test the distance from a point to plane
		cout << "Testing the distance from a point to a plane" << endl;

		Point *P12;
		Plane* plane2;

		P12 = new Point(9, 2, 3);

		plane2 = new Plane(2, 5, 7, 9);

		cout << "Created a point P(" << P12->getX() << ", " << P12->getY() << ", " << P12->getZ() << ")"
				<< " and a plane of equation (" << plane2->getA() << ")*x + (" << plane2->getB() << ")*y + (" << plane2->getC() << ")*z + (" << plane2->getD() << ") = 0" << endl
				<< tab << tab << "the distance from the point to the plane is: " << plane2->distanceToPlane(P12) << endl;
		*/

		int k = (wIn+1)/2;
		//double aTable[1<<k][1<<(k-1)];
		double **aTable;
		double **bTable;
		double **dTable;
		double error, errorMax;
		bool errorSatisfied;

		cout << "Beginning computations for Atan2, with wIn=" << wIn << " and wOut=" << wOut << endl;

		cout << tab << "Computing the parameters for the a, b and d tables" << endl;
		errorSatisfied = false;
		while(!errorSatisfied)
		{
			cout << tab << tab << "trying k=" << k << endl;

			cout << tab << tab << "computing the tables" << endl;

			//allocate the tables
			aTable = (double**)malloc((1<<k)*(1<<(k-1))*sizeof(double));
			for(int i=0; i<(1<<k); i++)
				aTable[i] = (double*)malloc((1<<(k-1))*sizeof(double));
			bTable = (double**)malloc((1<<k)*(1<<(k-1))*sizeof(double));
			for(int i=0; i<(1<<k); i++)
				bTable[i] = (double*)malloc((1<<(k-1))*sizeof(double));
			dTable = (double**)malloc((1<<k)*(1<<(k-1))*sizeof(double));
			for(int i=0; i<(1<<k); i++)
				dTable[i] = (double*)malloc((1<<(k-1))*sizeof(double));

			for(int i=0; i<(1<<k)-1; i++)
				for(int j=(1<<(k-1)); j<(1<<k)-1; j++)
				{
					Point *P1, *P2, *P3, *P4;
					Plane *plane;
					double a, b, c, d;
					double distance, distanceMin;
					double valI, valJ, increment;

					//create the actual values for the atan2
					valI = 1.0*i/(1<<k);
					valJ = 1.0*j/(1<<k);
					increment = 1.0/(1<<k);

					//create the points
					P1 = new Point(valI, 			valJ, 			atan2(valJ, 			valI		  ));
					P2 = new Point(valI+increment, 	valJ, 			atan2(valJ, 			valI+increment));
					P3 = new Point(valI, 			valJ+increment, atan2(valJ+increment, 	valI		  ));
					P4 = new Point(valI+increment, 	valJ+increment, atan2(valJ+increment, 	valI+increment));

					//create all the planes given by all the possible combinations of the 4 points
					// and determine the minimum distance from the plane given by three of the points
					// to the remaining point
					//the plane given by P1, P2 and P3
					plane = new Plane(P1, P2, P3);
					distance = plane->distanceToPlane(P4);
					if(distance<0)
						distance *= -1;
					distanceMin = distance;
					a = plane->getA();
					b = plane->getB();
					c = plane->getC();
					d = plane->getD();

					//the plane given by P1, P2 and P4
					plane = new Plane(P1, P2, P4);
					distance = plane->distanceToPlane(P3);
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

					//the plane given by P2, P3 and P4
					plane = new Plane(P2, P3, P4);
					distance = plane->distanceToPlane(P1);
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

					aTable[i][j-(1<<(k-1))] = -(a/c);
					bTable[i][j-(1<<(k-1))] = -(b/c);
					dTable[i][j-(1<<(k-1))] = -(d/c);

					cout << "\r                                  \r";
					cout << tab << tab << tab << std::setprecision(4)
						<< (((i<<(k-1))+j-(1<<(k-1)))/(1.0*((1<<(2*k-1)) + (1<<(k-1)) - 2))*100) << "% done";
				}

			cout << endl << tab << tab << "computed the tables" << endl;

			cout << tab << tab << "measuring the error" << endl;

			//measure the error
			errorMax = 0;
			for(int i=0; i<(1<<k)-1; i++)
			{
				for(int j=(1<<(k-1)); j<(1<<k)-1; j++)
				{
					for(int n=0; n<(1<<(wIn-k)); n++)
						for(int m=0; m<(1<<(wIn-k)); m++)
						{
							double valI, valJ, valIPrime, valJPrime;

							//create the actual values for the atan2
							valI = 1.0*((i<<(wIn-k))+n)/(1<<wIn);
							valJ = 1.0*((j<<(wIn-k))+m)/(1<<wIn);
							valIPrime = 1.0*i/(1<<k);
							valJPrime = 1.0*j/(1<<k);

							double referenceValue = atan2(valJ, valI);
							double computedValue  = aTable[i][j-(1<<(k-1))]*valIPrime + bTable[i][j-(1<<(k-1))]*valJPrime + dTable[i][j-(1<<(k-1))];

							error = fabs(referenceValue-computedValue);

							if(error > errorMax)
								errorMax = error;
						}

					/*
					cout << "\r                                  \r";
					cout << tab << tab << tab << std::setprecision(4)
					<< (((i<<(k-1))+j-(1<<(k-1)))/(1.0*((1<<(2*k-1)) + (1<<(k-1)) - 2))*100) << "% done";
					*/
				}

				cout << "\r                                  \r";
				cout << tab << tab << tab << std::setprecision(4) << ((1.0*i/((1<<k)-1))*100) << "% done";
			}

			cout << endl << tab << tab << "computed the maximum error" << endl;

			//check to see if the current precision k is enough
			//	if not, increase k and compute the maximum error again
			double errorLimit = 1.0/(1 << wOut);
			if(errorMax < errorLimit)
			{
				errorSatisfied = true;
			}else
			{
				cout << tab << tab << "error limit of " << errorLimit << " not satisfied at k=" << k
						<< ", with wIn=" << wIn << " and wOut=" << wOut << endl;
				cout << tab << tab << tab << "actual error errorMax=" << errorMax << endl;
				k++;

				//deallocate the space used for the tables, as a new try is needed
				for(int i=0; i<(1<<k); i++)
					free(aTable[i]);
				free(aTable);
				for(int i=0; i<(1<<k); i++)
					free(bTable[i]);
				free(bTable);
				for(int i=0; i<(1<<k); i++)
					free(dTable[i]);
				free(dTable);
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
