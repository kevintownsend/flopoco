/*
 * An FP logarithm for FloPoCo
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
#ifndef __FPLOG_HPP
#define __FPLOG_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "LZOC.hpp"
#include "LZOCShifterSticky.hpp"
#include "Shifters.hpp"
#include "fplogarithm/LogRangeRed.hpp"

class FPNumber;

#define MAXRRSTAGES 2000 // 4000 bits of accuracy should be enough for anybody

class FPLog : public Operator
{
public:
	FPLog(Target* target, int wE, int wF);
	~FPLog();

	//		Overloading the virtual functions of Operator
	void outputVHDL(std::ostream& o, std::string name);
	void emulate(TestCase * tc);
	void buildStandardTestCases(TestCaseList* tcl);
	/**Overloading the function of Operator with a function that tests only positive FP numbers (full range)*/
	void buildRandomTestCases(TestCaseList* tcl, int n); 

	int wE, wF;
	// The input sizes to the successive tables
	int a[MAXRRSTAGES]; 
	// The intermediate precision: at step i, the exp part is bounded by
	//    1 <= m < 1+2^-p[i]
	int p[MAXRRSTAGES]; 

	// The size of the product, and hence of the subword of Z[i] input to the mult
	int psize[MAXRRSTAGES]; 

	// The total size of non-zero bits, will be limited to wF+g
	int s[MAXRRSTAGES]; 

	// The numbers of bits to truncate to limit the size to wF+g
	int t[MAXRRSTAGES];

	// size before truncation, should be equal to s[i]+t[i]
	int sbt[MAXRRSTAGES];

	// The max number of stages
	int stages;

	int sfinal; // equal to s[stages+1]

	int pfinal;  // equal to s[stages+1]

	// The guard bits 
	int gLog;

	// The range reduction box
	LogRangeRed* rr;

	// This is the size of the virtual mantissa:
	// 1.000(p[i] zeroes)00000xxxxxxx
	// the xxx which will be actually computed have p[i] bits less 
	int sfullZ[MAXRRSTAGES]; 

	// A global boolean flag disabling the simulation of the fullsize mantissa
	// as soon as it would be more than 64 bits
	int fullSizeSim;

	// The target precision: numbers may be truncated so that their LSB has weight -target_prec 
	int target_prec;  

	// Various subcomponents
	Shifter *ao_lshift;   // ao stands for "almost one"
	Shifter *ao_rshift;
	LZOC *lzoc;
	LogRangeRed *rrbox;
	LZOCShifterSticky *final_norm;
};

#endif
