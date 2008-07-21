#ifndef LZOCShifterSticky_HPP
#define LZOCShifterSticky_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"

#include "Operator.hpp"

/*
 * A leading zero/one counter + shifter + sticky bit evaluator for FloPoCo
 *
 * Authors : Florent de Dinechin, Bogdan Pasca
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

/** 
 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
 */ 
class LZOCShifterSticky : public Operator
{
public:
	typedef enum {generic, specific} entityType_t;

	/* Constructor(s) and Destructors */
	LZOCShifterSticky(Target* target, int wIn, int wOut, bool compute_sticky, const int countType=-1);
	~LZOCShifterSticky();

	/*  accessor methods */
	void setEntityType(entityType_t eType);
	int getCountWidth() const;
	
	
	/* Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);
	void set_operator_name(std::string prefix, std::string postfix);
	
	/* Overloading the virtual functions of TestBench */
	TestIOMap getTestIOMap();
	void      fillTestCase(mpz_class a[]);

private:
	/** The number of bits of the input */
	int wIn;
	/** The number of bits of the count */
	int wCount;
	/** The number of bits of the shifted output */
	int wOut;
	/** If true, compute the sticky bit. If false, save this hardware */
	bool computeSticky; 
	/** -1|0|1. If -1 is present then generic LZOC is instatiated */
	int countType;
	/** The names of the signals, just to make code more readable */
	string level[42];  
	/** Same but possibly delayed  */
	string leveld[42]; 
	/** Their size. Do we need to count more than 2^42 bits in FloPoCo? */ 
	int size[42];      
	/** Entity type. Can be either generic or specific */
	entityType_t entityType;
	/** if boolean true, the corresponding level signal is registered*/ 
	bool level_registered [128];
	/** the depths for the levels of the architecture */	
	int countDepth[42];
	/** utilitary var */		
	mpz_class maxValue;
};

#endif
