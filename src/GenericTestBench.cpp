/*
 * Generic Test Bench generator for FloPoCo
 *
 * Authors : Bogdan Pasca
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
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "GenericTestBench.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	GenericTestBench:: GenericTestBench(Target* target, string name, int depth, int wIn) :
		Operator(target), _name(name){
 
		setOperatorName();
	
		if(isSequential()) 
			{
				setPipelineDepth(depth);
			}
	
	
	
		addInput ("X", wIn);
		addInput ("Y", wIn);
		addOutput("R", 2*wIn); /* wOut_ = wInX_ + wInY_ */
  
	}

	GenericTestBench::~GenericTestBench() {
	}

	void GenericTestBench::setOperatorName(){	
		/* Name Setup procedure
		 *  The name has the format: GenericTestBench_wInX__wInY_
		 *  wInX_ = width of the X input
		 *  wInY_ = width of the Y input
		 */  
		ostringstream name;
		name <<_name;
		uniqueName_ = name.str(); 
	}

	void GenericTestBench::fillTestCase(mpz_class a[])
	{
		mpz_class &x = a[0];
		mpz_class &y = a[1];
		mpz_class &r = a[2];

		r = x * y;
	}

}
