#ifndef LZOCShifterSticky_HPP
#define LZOCShifterSticky_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

/*
 * A leading zero/one counter + shifter +sticky bit evaluator for FloPoCo
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



class LZOCShifterSticky : public Operator
{
public:
	LZOCShifterSticky(Target* target, int wIn, int wOut, bool compute_sticky);
	~LZOCShifterSticky();

	int wIn;
	int wOut;

	/** The number of bits of the count */
	int wCount;


	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);

private:
	/** if true, compute the sticky bit. If false, save this hardware */
	bool compute_sticky; 
	string level[42]; // the names of the signals, just to make code more readable 
	string leveld[42]; // same but possibly delayed 
	int size[42]; // Their size. Do we need to count more than 2^42 bits in FloPoCo? 
};



#endif
