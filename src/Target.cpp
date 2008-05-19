/*
 * The abstract class that models different chips. Classes for real
 * chips inherit from this one.
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

#include <iostream>
#include <sstream>
#include "Target.hpp"
using namespace std;

extern int verbose;

void Target::set_pipelined() {
	_pipeline=true;
}

void Target::set_not_pipelined() {
	_pipeline=false;
}

bool Target::is_pipelined() {
	return _pipeline;
}

int Target::lut_inputs() {
	return _lut_inputs;
}

double Target::frequency(){
	return _frequency;
};

void Target::set_frequency(double f){
	_frequency=f;
};

void Target::set_use_hard_multipliers(bool v){
				_use_hard_multipliers = v;  
};

bool Target::get_use_hard_multipliers(){
				return _use_hard_multipliers;
}; 
