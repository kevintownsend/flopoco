/*
 * The range reduction box for the FP logarithm
 * 
 * Author : Florent de Dinechin
 *
 * For a description of the algorithm, see the Arith17 paper. 
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
#ifndef LOGRANGERED_H
#define LOGRANGERED_H

#include "../Operator.hpp"
#include "FirstInvTable.hpp"
#include "FirstLogTable.hpp"
#include "SecondInvTable.hpp"
#include "OtherLogTable.hpp"

class FPLog;

class LogRangeRed : public Operator
{
 public:

	LogRangeRed(Target *target,
					FPLog *fplog)  ;

		
  virtual ~LogRangeRed() ;

  // The table objects
  FirstInvTable* it0;
  SecondInvTable* it1;
  FirstLogTable* lt0;
  OtherLogTable* lt[2000];


  // This is the size of the virtual mantissa:
  // 1.000(p[i] zeroes)00000xxxxxxx
  // the xxx which will be actually computed have p[i] bits less 
	//  int sfullZ[42]; 
private:

	/** The FPLog instantiating this LogRangeRed, if any */
	FPLog *fplog;

	/** Equal to fplog->stages, increases readibility */
	int stages;
	/** Equal to fplog->a, increases readibility */
	int *a;
	/** Equal to fplog->p, increases readibility */
	int *p;
	/** Equal to fplog->psize, increases readibility */
	int *psize;
	/** Equal to fplog->s, increases readibility */
	int *s;

};


#endif
