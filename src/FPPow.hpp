/*
 * An FP power for FloPoCo
 *
 * Author : Pedro Echevarria, Florent de Dinechin
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
#ifndef __FPPOW_HPP
#define __FPPOW_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "LZOC.hpp"
#include "LZOCShifterSticky.hpp"
#include "Shifters.hpp"


namespace flopoco{


	class FPPow : public Operator
	{
	public:
		FPPow(Target* target, int wE, int wF, int logTableSize, int expTableSize, int expDegree, int expG, int logG );
		~FPPow();

		void compute_error(mpfr_t & r, mpfr_t &epsE, mpfr_t& epsM, mpfr_t& epsL );

		//		Overloading the virtual functions of Operator
		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);
		/**Overloading the function of Operator with a function that tests only positive FP numbers (full range)*/
		void buildRandomTestCases(TestCaseList* tcl, int n); 

		int wE, wF;
	};
}
#endif
