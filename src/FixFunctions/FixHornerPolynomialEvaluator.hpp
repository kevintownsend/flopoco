/*
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2008-2014.
  All rights reserved.

*/
#ifndef __FIXHORNERPOLYNOMIALEVALUATOR_HPP
#define __FIXHORNERPOLYNOMIALEVALUATOR_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "Table.hpp"
#include "DualTable.hpp"


namespace flopoco{


  class FixHornerPolynomialEvaluator : public Operator
  {
  public:
    /** The constructor with manual control of all options
     * @param wE exponent size
     */

    FixHornerPolynomialEvaluator(Target* target, 
				 int wX, 
				 int degree, 
				 int LSB, vector<int> MSB, 
				 bool signedCoeffs=true, 
				 map<string, double> inputDelays = emptyDelayMap);
    ~FixHornerPolynomialEvaluator();
		

  private:
    int wX;                           /**< size of X in bits; X is unsigned, in [0,1) */
    int degree;                       /**< degree of the polynomial */
    int LSB;                          /**< common weight of the LSBs of the polynomial approximations */
    bool signedCoeffs;                /**< if false, all the coeffs are unsigned and the operator may use unsigned arithmetc. 
					 Usually true unless known Taylor etc */
    vector<int> MSB;                  /**< vector of MSB weights for each coefficient */

  };

}
#endif
