/*
	A plane in the 3D space

	Author: Matei Istoan

	This file is part of the FloPoCo project

	Initial software.
	Copyright © INSA Lyon, INRIA, CNRS, UCBL,
	2012-2014.
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
#include "../Operator.hpp"
#include "Plane.hpp"

using namespace std;

namespace flopoco {

	Plane::Plane(Point* P1_, Point* P2_, Point* P3_):
			P1(P1_), P2(P2_), P3(P3_),
			x1(P1_->getX()), y1(P1_->getY()), z1(P1_->getZ()),
			x2(P2_->getX()), y2(P2_->getY()), z2(P2_->getZ()),
			x3(P3_->getX()), y3(P3_->getY()), z3(P3_->getZ())
	{
		P1 = new Point(P1_);
		P2 = new Point(P2_);
		P3 = new Point(P3_);

		a = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1);
		b = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1);
		c = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
		d = x1*((z2-z1)*(y3-y1) - (y2-y1)*(z3-z1))
				+ y1*((x2-x1)*(z3-z1) - (z2-z1)*(x3-x1))
				+ z1*((y2-y1)*(x3-x1) - (x2-x1)*(y3-y1));
	}

	Plane::Plane(double a_, double b_, double c_, double d_)
	{
		a = a_;
		b = b_;
		c = c_;
		d = d_;

		P1 = NULL;
		P2 = NULL;
		P3 = NULL;

		x1 = 0; y1 = 0; z1 = 0;
		x2 = 0; y2 = 0; z2 = 0;
		x3 = 0; y3 = 0; z3 = 0;
	}

	bool Plane::pointsCollinear(Point* P1, Point* P2, Point* P3)
	{
		double x1=P1->getX(), y1=P1->getY(), z1=P1->getZ();
		double x2=P2->getX(), y2=P2->getY(), z2=P2->getZ();
		double x3=P3->getX(), y3=P3->getY(), z3=P3->getZ();
		double xRatio, yRatio, zRatio;

		//sanity checks
		if(x1==x3)
			return (x1==x2);
		if(y1==y3)
			return (y1==y2);
		if(z1==z3)
			return (z1==z2);

		xRatio = (x2-x1)/(x3-x1);
		yRatio = (y2-y1)/(y3-y1);
		zRatio = (z2-z1)/(z3-z1);

		if((xRatio==yRatio) && (yRatio==zRatio) && (xRatio==zRatio))
			return true;
		else
			return false;
	}

	double Plane::distanceToPlane(Point* P)
	{
		double x0=P->getX(), y0=P->getY(), z0=P->getZ();

		return ((a*x0+b*y0+c*z0+d) / sqrt(a*a+b*b+c*c));
	}

	double Plane::getA()
	{
		return a;
	}

	double Plane::getB()
	{
		return b;
	}

	double Plane::getC()
	{
		return c;
	}

	double Plane::getD()
	{
		return d;
	}

	Plane::~Plane()
	{
		if(P1)
			delete P1;
		if(P2)
			delete P2;
		if(P3)
			delete P3;
	}

}
