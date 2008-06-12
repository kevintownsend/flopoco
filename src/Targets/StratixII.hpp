/*
 * A model of Stratix II FPGA 
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

#ifndef STRATIXII_HPP
#define  STRATIXII_HPP
#include "../Target.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

class StratixII : public Target
{
 public:


  // generic constructor
  StratixII() : Target()
  {
    // all these values are set more or less randomly, to match  Stratix II more or less
    _fastcarry_delay = 3.4e-11; // s    
    _elem_wire_delay = 0.3e-11;
    _lut2lut_delay = 1.5e-10;
    _lut_delay =  1.5e-9; 
    _mult_x_inputs=18;
    _mult_y_inputs=18;
  }

  // constructor for real chips: TODO

  virtual ~StratixII() {}

  // overloading the virtual functions of Target

  double carry_propagate_delay();

  double adder_delay(int size);

  double local_wire_delay();

  double lut_delay();

  double distant_wire_delay(int n);
  
  bool suggest_submult_size(int &x, int &y, int wInX, int wInY);
  
  bool suggest_subadd_size(int &x, int wIn);


private:

  double _fastcarry_delay; // in seconds

  double _lut2lut_delay; // in seconds

  double _elem_wire_delay; // in seconds

  double _lut_delay; // in seconds
};


#endif
