/*
 * A step of SRT2 square root for FloPoCo
 *
 * Author : Jeremie Detrey, Florent de Dinechin
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


#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "../utils.hpp"

#include "SqrtSRT2Step.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


SqrtSRT2Step::SqrtSRT2Step(Target* target, int wF, int step) :
	Operator(target), wF(wF), step(step) {

	ostringstream name;

	setCopyrightString("Jeremie Detrey, Florent de Dinechin (2002-2009)");
	name<<"SqrtSRT2Step_"<<wF; 
	setName(name.str()); 

	addInput ("x", wF+4);
	addInput ("s", step);
	addOutput ("d");
	addOutput ("w", wF+4);
	

	vhdl << tab << declare("d0") << " <= x("<< wF+3<<");" << endl;
	vhdl << tab << declare("x0",wF+5) << " <= x & \"0\";" << endl;
	vhdl << tab << declare("s0",step+1) << " <= \"0\" & s;" << endl;
	vhdl << tab << declare("ds",step+3) << " <= s0" << range(step,1) << "& (not d0) & d0 & \"1\";" << endl;
	vhdl << tab << declare("x1",step+3) << " <= x0" << range(wF+4, wF+2-step) << ";" << endl;
  	vhdl << tab << "with d0 select" << endl
		  << tab << tab <<  declare("w1", step+3) << " <= x1 - ds when '0'," << endl
        << tab << tab << "      x1 + ds when others;" << endl;
	vhdl << tab << declare("w0",wF+5) << " <= w1";
	if(step <= wF+1) 
		vhdl << " & x0" << range(wF+1-step, 0) << ";" << endl;  
	else
		vhdl << ";" << endl; 

	vhdl << tab << "d <= d0;" << endl;
	vhdl << tab << "w <= w0" << range(wF+3, 0) << ";" << endl;
}

SqrtSRT2Step::~SqrtSRT2Step() {
}


