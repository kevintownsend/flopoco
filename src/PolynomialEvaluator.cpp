/*
 * A polynomial evaluator for FloPoCo
 *
 * Author : Bogdan Pasca
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
#include <stdio.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "PolynomialEvaluator.hpp"

using namespace std;

namespace flopoco{

	PolynomialEvaluator::PolynomialEvaluator(Target* target, vector<FixedPointCoefficient*> coef, YVar* y, int targetPrec):
		Operator(target), coef_(coef), y_(y), targetPrec_(targetPrec){
		srcFileName = "PolynomialEvaluator";

		int degree = coef.size()-1;
		setName( join("PolynomialEvaluator_d",degree) );
		
		REPORT(DETAILED, "Polynomial to evaluate: " << printPolynomial(coef, y));
		
		/* FIXME most of the parameters fixed here should belong to the constructor */
//		int rWidth = 10;
		
		int g = 0; /* the number of guard bits to the right of the lsb of the coeff. 
		              this number might increase during iterations 
		            */
			
		mpfr_t roundingError;                      /* keep track of the rounding error */
		mpfr_init2(roundingError, 100);            /* initialize the precision */

		/* ================== I/O declarations ===============================*/
		addInput("X", y->getSize());
		for (uint32_t i=0; i <= unsigned(degree); i++)
			addInput(join("C",i), coef[i]->getSize() );

		/* initialize first output level */
		LevelSignal* outputSignal     =  new LevelSignal(join("level",0), coef[degree]->getSize(), coef[degree]->getWeight());
		LevelSignal* firstOutputSignal = new LevelSignal(outputSignal);

		LevelSignal* inputSignal;
		
		int realPrecOld = 0;
		int realPrecCurrent = 0;
		
		vector<LevelSignal*> levelSignalList;
		
		while ( realPrecCurrent > targetPrec){
			levelSignalList.push_back(firstOutputSignal);
			mpfr_set_si(roundingError, 0 , GMP_RNDN);  /* rounding is initialized to 0 */
			
			REPORT( DETAILED,  " ======================================================= " );
			REPORT( DETAILED,  " ====================== g=" << g << " ================================= " );
			REPORT( DETAILED,  "OutputSignal name=" <<firstOutputSignal->getName() << " size=" << firstOutputSignal->getSize() << " weight = " << firstOutputSignal->getWeight());
			for(uint32_t i=0; i < unsigned(2*degree); i++) {
				REPORT( DETAILED,  " ------" );
				if (i!=0)
					inputSignal  = new LevelSignal(outputSignal); /* the input level at this iteration is the output level at the previous one */
				else
					inputSignal  = new LevelSignal(firstOutputSignal);
				
				if ( i % 2 == 0 ){
					/* the evaluation path is composed of multiplications an additions. 
					   This branch handles multiplication by X */
					mpfr_t mError;
					mpfr_init2(  mError, 300);
					mpfr_set_ui( mError, 2, GMP_RNDN);
					mpfr_pow_si( mError, mError, y->getWeight(), GMP_RNDN);
					mpfr_mul( roundingError, roundingError, mError, GMP_RNDN);
					int truncationULP = inputSignal->getWeight() + y->getWeight() - (inputSignal->getSize()   + y->getSize());
					int roundingULP = (mpfr_get_d(roundingError,GMP_RNDN) > 0 ? mpfr_get_exp(roundingError): MINF);
					int realULP =  max ( truncationULP, roundingULP);
					int realSize = inputSignal->getWeight() + y->getWeight() - realULP;

					outputSignal = new LevelSignal( join("level",i+1), 
						                            realSize,  
						                            inputSignal->getWeight() + y->getWeight());
					
					levelSignalList.push_back(outputSignal); 			

//					REPORT( DETAILED,  "InputSignal  name="  << inputSignal->getName() << " size=" << inputSignal->getSize() << " weight = " << inputSignal->getWeight() << endl;
//					REPORT( DETAILED,  "y size=" << y->getSize() << " weight="<<y->getWeight() << endl;
					REPORT( DETAILED,  "OutputSignal name=" << outputSignal->getName() << " size=" << outputSignal->getSize() << " weight = " << outputSignal->getWeight() );
					REPORT( DETAILED,  " -- " << "Error produced at step " << i << " is " << (mpfr_get_d(roundingError,GMP_RNDN)>0?mpfr_get_exp(roundingError):0) );

				}else{
					if (i != unsigned(2*degree-1)){ /* except last addition a0 + prevLevel */
						int cIndex = coef.size()-(i+3)/2;
						int outputLevelShift = max(inputSignal->getWeight(),coef[cIndex]->getWeight())+1;
						int lsbCoef = coef[cIndex]->getWeight() - coef[cIndex]->getSize() - g;
						int lsbLev  = inputSignal->getWeight()  - inputSignal->getSize();
						int lowerBound = max (lsbCoef, lsbLev);
						int outputLevelSize = outputLevelShift - lowerBound;
					
//						REPORT( DETAILED,  "Input level weight " << inputSignal->getWeight() << " of size " << inputSignal->getSize() << endl; 
//						REPORT( DETAILED,  "Input coeff weight " << coef[cIndex]->getWeight() << " of size " << coef[cIndex]->getSize() << endl;
//						REPORT( DETAILED,  "LSB coef = " << lsbCoef << " LSB Lev = " << lsbLev << endl;
//						REPORT( DETAILED,  "outputLevelShift " << outputLevelShift << " with outputLevelSize = " << outputLevelSize << endl;
					
						outputSignal = new LevelSignal( join("level", i+1),
							                            outputLevelSize,
							                            outputLevelShift); 
					
						levelSignalList.push_back(outputSignal);
						
						mpfr_t aError;
						mpfr_init2( aError, 300);
						mpfr_set_ui( aError, 2, GMP_RNDN);
						mpfr_pow_si( aError, aError, lowerBound, GMP_RNDN);
						mpfr_max( roundingError, roundingError, aError , GMP_RNDN);

						REPORT( DETAILED,  "OutputSignal name=" <<outputSignal->getName() << " size="  << outputSignal->getSize() << " weight = " << outputSignal->getWeight() );
						REPORT( DETAILED,  " -- " << "Error produced at step " << i << " is " << (mpfr_get_d(roundingError,GMP_RNDN)>0?mpfr_get_exp(roundingError):0));
					}
				}
	//			vhdl << " --" << " Result size = " << outputLevelSize << " Result Shift = " << outputLevelShift << endl;
	////			double d = mpfr_get_d( roundingError, GMP_RNDN);
	//			vhdl << " --" << " The rounding error at this step is = " << mpfr_get_exp(roundingError) << endl;
			}
			
			REPORT( DETAILED,  " The error is of the order " << mpfr_get_exp( roundingError) );
			realPrecOld     = realPrecCurrent;
			realPrecCurrent = mpfr_get_exp( roundingError);
			
			if (realPrecCurrent == realPrecOld){
				throw "impossible to reach precision with current coefficients";
			}
			if (realPrecCurrent > targetPrec){
				levelSignalList.clear();
				g +=   realPrecCurrent - targetPrec;
			}
		}
		/* print signals with charateristics */
		for (uint32_t j=0; j< levelSignalList.size(); j++)
			REPORT( DETAILED,  levelSignalList[j]->getName() << " size="<<levelSignalList[j]->getSize() << " weight=" << levelSignalList[j]->getWeight() );
		
		
//		vhdl << tab << "R <= " << join("level",2*degree) << ";" << endl;		 
	}

	PolynomialEvaluator::~PolynomialEvaluator() {
	}

}



//		outputLevel << join("level",0);
//		outputLevelSize  = coef[degree]->getSize();       /* the width of the level0 signal = width of (C_n) */
//		outputLevelShift = coef[degree]->getWeight();

// 		/*TODO REMOVE: iterate and finally generate */
//		REPORT( DETAILED,  tab << declare(outputLevel.str(), coef[degree]->getSize()) << " <= " << join("C",degree) << ";" << endl;
//		REPORT( DETAILED,  " -- " << "Error starts " << mpfr_get_d(roundingError, GMP_RNDN) << endl; /* this is obviously 0 */




//even

//				vhdl << tab << declare( outputLevel.str(), outputLevelSize) << " <= " << inputLevel.str() << " * " << "X" << ";" << endl;
//				
//				
//				/* comment vhdl. TODO remove after this becomes stable */
//				vhdl << " -- " << "Input level size = " << inputLevelSize << " shift =" << inputLevelShift << endl;
//				vhdl << " -- " << "X           size = " << xWidth         << " shift =" << 1-xWidth << endl; 

//odd


//					int msb1 = inputSignal->getWeight() inputLevelShift;
//					int msb2 = coef[coef.size()-(i+3)/2]->getSize() +  coef[coef.size()-(i+3)/2]->getWeight();

//					int lsb1 = inputLevelShift;
//					int lsb2 = coef[coef.size()-(i+3)/2]->getWeight(); /* the coef shift */
//				
//					int upperBound = max (msb1, msb2);
//					int lowerBound;
//					if (lsb1 < lsb2)
//						lowerBound = lsb2-g;
//					else
//						lowerBound = lsb2;

//					outputLevelSize = upperBound - lowerBound + 1;
//					outputLevelShift = lowerBound;
//				
//					vhdl << tab << declare( outputLevel.str(), outputLevelSize) << " <= " << align(upperBound+1-msb1, inputLevel.str()+range(inputLevelSize-1, inputLevelSize-(outputLevelSize-(upperBound+1-msb1))), lsb1 - outputLevelShift) << " + " 
//						                                                                  << align(upperBound+1-msb2, join("C", coef.size()-(i+3)/2), lsb2 - outputLevelShift) <<  ";" << endl;
//					vhdl << " -- " << "Coef        size = " << coef[coef.size()-(i+1)/2]->getSize() << " shift =" << coef[coef.size()-(i+3)/2]->getWeight() << endl; 
//					vhdl << " -- " << "Input level size = " << inputLevelSize << " shift =" << inputLevelShift << endl;


//last add

//				}else{
////					/* here we NEED to round in order to loose just 1/2 ULPs */
////					int msb1 = inputLevelSize + inputLevelShift;
////					int msb2 = coef[coef.size()-(i+3)/2]->getSize() +  coef[coef.size()-(i+3)/2]->getWeight();

////					int lsb1 = inputLevelShift;
////					int lsb2 = coef[coef.size()-(i+3)/2]->getWeight(); /* the coef shift */
////				
////					int upperBound = max (msb1, msb2);
////					int lowerBound = min (lsb1, lsb2);

////					outputLevelSize = upperBound - lowerBound + 1;
////					outputLevelShift = lowerBound;

////					vhdl << tab << declare( outputLevel.str(), outputLevelSize) << " <= " << align(upperBound+1-msb1, inputLevel.str()+range(inputLevelSize-1, inputLevelSize-(outputLevelSize-(upperBound+1-msb1))), lsb1 - outputLevelShift) << " + " 
////						                                                                  << align(upperBound+1-msb2, join("C", coef.size()-(i+3)/2), lsb2 - outputLevelShift) <<  ";" << endl;
////					vhdl << " -- " << "Input level size = " << inputLevelSize << " shift =" << inputLevelShift << endl;
////					vhdl << " -- " << "Coef        size = " << coef[coef.size()-(i+1)/2]->getSize() << " shift =" << coef[coef.size()-(i+3)/2]->getWeight() << endl; 
////					addOutput("R", outputLevelSize);
//				}



