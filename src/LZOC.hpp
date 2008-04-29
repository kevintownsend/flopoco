#ifndef LZOC_HPP
#define LZOC_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

/*
 * A leading zero/one counter for FloPoCo
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


/* 
 * Recursive structure with wOut stages. At most 2^wOut-1 leading zeros are counted.
 */


class LZOC : public Operator
{
public:
	LZOC(Target* target, int wIn, int wOut);
	~LZOC();

	int wIn;
	int wOut;
	int p2wOut;


	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);
};



#endif
