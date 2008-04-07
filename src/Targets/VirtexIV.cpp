/*
 * A model of FPGA that works well enough for Altera and Virtex chips. 
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
#include "VirtexIV.hpp"
#include <iostream>
#include <sstream>


double VirtexIV::adder_delay(int size) {
  return _lut_delay  +  size * _fastcarry_delay; 
};

double VirtexIV::carry_propagate_delay() {
  return  _fastcarry_delay; 
};

double VirtexIV::local_wire_delay(){
  return _lut2lut_delay;
};

double VirtexIV::distant_wire_delay(int n){
  return n*_elem_wire_delay;
};

double VirtexIV::lut_delay(){
  return _lut_delay;
};


bool VirtexIV::suggest_submult_size(int &x, int &y, int wInX, int wInY){
	if (get_use_hard_multipliers()){
		x=17;
		y=17;
		return true;	//TODO
	}else{
		return true;    //TODO
	}
};	 
	 
	 
bool VirtexIV::suggest_subadd_size(int &x, int wIn){

	int chunk_size = (int)floor( (1./frequency() - lut_delay()) / carry_propagate_delay()); // 1 if no need for pipeline
	x = chunk_size;		
	
	if (x>0) return true;
	else {
		x=1;		
		return false;
	} 
};
  
