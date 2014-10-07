/*
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2008-2014.
  All rights reserved.

*/
#ifndef __FIXHORNEREVALUATOR_HPP
#define __FIXHORNEREVALUATOR_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "Table.hpp"
#include "DualTable.hpp"


namespace flopoco{

	/** An Horner polynomial evaluator computing just right.
	 It assumes the input X is an unsigned number in [0, 1[ so msbX=-wX.
	*/

  class FixHornerEvaluator : public Operator
  {
  public:
    /** The constructor with manual control of all options.
     * @param wE exponent size
     */

    FixHornerEvaluator(Target* target, 
											 int wX, 
											 int degree, 
											 vector<int> coeffMSB, 
											 vector<int> coeffLSB, 
											 bool signedCoeffs=true, 
											 bool finalRounding=true,
											 bool oplainStupidVHDL=false,
											 map<string, double> inputDelays = emptyDelayMap);
    ~FixHornerEvaluator();
		

  private:
    int wX;                           /**< size of X in bits; X is unsigned, in [0,1) */
    int degree;                       /**< degree of the polynomial */
    vector<int> coeffMSB;                  /**< vector of MSB weights for each coefficient */
    vector<int> coeffLSB;                  /**< vector of LSB weights for each coefficient */
    bool signedCoeffs;                /**< if false, all the coeffs are unsigned and the operator may use unsigned arithmetc. 
																				 Usually true unless known Taylor etc */
		bool finalRounding;               /** If true, the operator returns a rounded result (i.e. add the half-ulp then truncate)
																					If false, the operator returns the full, unrounded results including guard bits */
		bool plainStupidVHDL;             /** If true, generate portable VHDL with multiplications as "*" and additions as "+"
																					If false, use FloPoCo FixMultAdd operator*/
    vector<int> coeffSize;            /**< vector of the sizes of the coefficients, computed out of MSB and LSB. See FixConstant.hpp for the constant format */

  };

}
#endif
