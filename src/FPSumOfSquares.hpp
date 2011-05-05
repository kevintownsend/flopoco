/*
 * An FP FPSumOfSquares operator for FloPoCo
 *
 * the operator computes R <= X^2 + Y^2 + Z^2 
 * 
 * There are two versions, selectable by the useFPOperators parameter.
 * One combines existing FloPoCo floating-point operators, and the other
 * one is a specific datapath designed in the FloPoCo philosophy.
 * 
 * Author : Florent de Dinechin, Bogdan Pasca
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
#ifndef __FPSumOfSquares_HPP
#define __FPSumOfSquares_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "Shifters.hpp"
#include "LZOC.hpp"
#include "LZOCShifterSticky.hpp"
#include "IntAdder.hpp"
#include "IntNAdder.hpp"
#include "IntMultiplier.hpp"
#include "IntSquarer.hpp"
#include "FPMultiplier.hpp"
#include "FPAdderSinglePath.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"

namespace flopoco{

	class FPSumOfSquares : public Operator
	{
	public:
		FPSumOfSquares(Target* target, int wE, int wF, int optimize);
		~FPSumOfSquares();

		void emulate(TestCase * tc);

		void buildRandomTestCases(TestCaseList* tcl, int n);

		int wE, wF;

	};
}
#endif
