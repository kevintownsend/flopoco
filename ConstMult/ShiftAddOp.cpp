/*
 * Shift-and-add operators for integer constant multiplication
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
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "IntConstMult.hpp"
#include "ShiftAddOp.hpp"
#include "ShiftAddDag.hpp"

using namespace std;




ShiftAddOp::ShiftAddOp(ShiftAddDag* impl, OpType op, ShiftAddOp* i, int s, ShiftAddOp* j) :
  impl(impl), op(op), i(i), s(s), j(j) {
  n=impl->computeConstant(op, i, s, j);
  // now compute the size according to this constant
  if (n >= 0)
    size = intlog2(n-1) + impl->icm->xsize; // -1 so that powers of two are handled properly
  else
    size = intlog2(-n-1) + impl->icm->xsize + 1; // we need a sign bit

  // build the variable name
  ostringstream o; 
  if(n==1)         o <<"X"; 
  else  if(n>=0)   o<<"P"<<n<<"X";  
  else             o<<"M"<<-n<<"X";
  name =  o.str();

  // and add this object to the dictionary of the implementation
  if (op!=X)
    impl->saolist.push_back(this);
}


 
