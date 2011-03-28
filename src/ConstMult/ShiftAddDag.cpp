/*
 * Shift-and-add DAG for integer constant multiplication
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
#include "../Operator.hpp"
#include "ShiftAddOp.hpp"
#include "ShiftAddDag.hpp"
#include "IntConstMult.hpp"

using namespace std;

namespace flopoco{

	mpz_class ShiftAddDag::computeConstant(ShiftAddOpType op, ShiftAddOp* i, int s, ShiftAddOp* j) {
		switch(op) {
		case X:
			if(i!=NULL || j!=NULL) {
				throw string("ERROR in ShiftAddDag::computeConstant: Unexpected non-null pointer in X\n");
			}
			return mpz_class(1);
			break;

		case Add: 
			if(i==NULL || j==NULL) {
				throw string("ERROR in ShiftAddDag::computeConstant: Unexpected non-null pointer in Add\n");
			}
			return (i->n << s) +  j->n;
			break;

		case Sub: 
			if(i==NULL || j==NULL) {
				throw string("ERROR in ShiftAddDag::computeConstant: Unexpected non-null pointer in Sub\n");
			}
			return (i->n << s) -  j->n;
			break;

		case RSub: 
			if(i==NULL || j==NULL) {
				throw string("ERROR in ShiftAddDag::computeConstant: Unexpected non-null pointer in RSub\n");
			}
			return  j->n - (i->n << s) ;
			break;

		case Shift:
			if(i==NULL || j!=NULL) {
				throw string("ERROR in ShiftAddDag::computeConstant: Unexpected pointer in Shift\n");
			}
			return i->n << s ;
			break;

		case Neg:
			if(i==NULL) {
				throw string("ERROR in ShiftAddDag::computeConstant: Unexpected null pointer in Neg\n");
			}
			return - i->n;
			break;
		}

		return 0; 	// should never be reached, added to suppress a warning
	}



	// This method looks up in the current Dag if the requireed op
	// exists, and either returns a pointer to it, or creates the
	// corresponding node.
	ShiftAddOp* ShiftAddDag::provideShiftAddOp(ShiftAddOpType op, ShiftAddOp* i, int s, ShiftAddOp* j){
		mpz_class n=this->computeConstant(op, i, s, j);
		if(n==1)
			return PX;
		else {
			for(unsigned int ii=0; ii<this->saolist.size(); ii++) {
				if (n==saolist[ii]->n) 
					return saolist[ii];
			}
			return new ShiftAddOp(this, op, i, s, j);
		}
	}



}
