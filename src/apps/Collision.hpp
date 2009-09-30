/*
 * An FP collision operator for FloPoCo

This is mostly an example of a coarse-grain, nonstandard operator.

A 3D collision detector inputs 3 FP coordinates X,Y and Z and
the square of a radius, R2, and computes a boolean predicate which is
true iff X^2 + Y^2 + Z^2 < R2

There are two versions, selectable by the useFPOperators parameter.
One combines existing FloPoCo floating-point operators, and the other
one is a specific datapath designed in the FloPoCo philosophy.

As this is a floating-point operator, each versions has its "grey
area", when X^2+Y^2+Z^2 is very close to R2.  In this case the
predicate may be wrong with respect to the infinitely accurate result.

The grey area of the combination of FP operators is about 2.5 units in
the last place of R2.  The pure FloPoCo version (which is a lot
smaller and faster) is more accurate, with a grey area smaller than 1
ulp of R2.

 *
 * Author : Florent de Dinechin
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/
#ifndef __COLLISION_HPP
#define __COLLISION_HPP
#include <vector>
#include <sstream>

#include "../Operator.hpp"
#include "../Shifters.hpp"
#include "../LZOC.hpp"
#include "../LZOCShifterSticky.hpp"
#include "../IntAdder.hpp"
#include "../IntMultiplier.hpp"
#include "../IntSquarer.hpp"
#include "../FPMultiplier.hpp"
#include "../FPAdder.hpp"

namespace flopoco{

	class Collision : public Operator
	{
	public:
		Collision(Target* target, int wE, int wF, int optimize);
		~Collision();

		void emulate(TestCase * tc);

		void buildRandomTestCases(TestCaseList* tcl, int n);

		int wE, wF;

	};
}
#endif
