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

mpz_class ShiftAddDag::computeConstant(OpType op, ShiftAddOp* i, int s, ShiftAddOp* j) {
	switch(op) {
	case X:
		if(i!=NULL || j!=NULL) {
			cerr << "ERROR Unexpected non-null pointer in computeConstant(X,..). \n";
			exit(EXIT_FAILURE);
		}
		return 1;

	case Add: 
		if(i==NULL || j==NULL) {
			cerr << "ERROR Unexpected null pointer in computeConstant(Add, ...), exiting. \n";
			exit(EXIT_FAILURE);
		}
		// the constant by which this variable multiplies 
		return (i->n << s) +  j->n;

	case Shift:
		if(i==NULL || j!=NULL) {
			cerr << "ERROR In inputs to computeConstant(Shift,...), exiting. \n";
			exit(EXIT_FAILURE);
		}
		return i->n << s ;

	case Neg:
		if(i==NULL) {
			cerr << "ERROR In ShiftAddOp, cannot construct such a Neg, exiting. \n";
			exit(EXIT_FAILURE);
		}
		return - i->n;
	}

}



// This method looks up in the current Dag if the requireed op
// exists, and either returns a pointer to it, or creates the
// corresponding node.
ShiftAddOp* ShiftAddDag::provideShiftAddOp(OpType op, ShiftAddOp* i, int s, ShiftAddOp* j){
	mpz_class n=this->computeConstant(op, i, s, j);
	for(int ii=0; ii<this->saolist.size(); ii++) {
		if (n==saolist[ii]->n) 
			return saolist[ii];
	}
	return new ShiftAddOp(this, op, i, s, j);
}







