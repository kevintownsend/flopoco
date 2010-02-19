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
#include <iomanip>
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
		
		int g = 0;                       /* the number of guard bits */			
		mpfr_t roundingError;            /* keep track of the rounding error */
		mpfr_init2(roundingError, 100);  /* initialize the precision */

		/* ================== I/O declarations ==================*/
		addInput("X", y->getSize());
		for (uint32_t i=0; i <= unsigned(degree); i++)
			addInput(join("C",i), coef[i]->getSize() );

		/* initialize first output level */
		LevelSignal* outputSignal      = new LevelSignal(join("level",0), coef[degree]->getSize(), coef[degree]->getWeight());
		LevelSignal* firstOutputSignal = new LevelSignal(outputSignal);
		LevelSignal* inputSignal;
		
		int realPrecOld = 0;
		int realPrecCurrent = 0;
		
		vector<LevelSignal*> levelSignalList;
		
		while ( realPrecCurrent > targetPrec){
			levelSignalList.push_back(firstOutputSignal);
			mpfr_set_si(roundingError, 0 , GMP_RNDN);  /* rounding is initialized to 0 */
			
			REPORT( DETAILED,  " -------------------- g=" << setw(2) << g << " -------------------- " );
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
		for (uint32_t j=0; j< levelSignalList.size(); j++){
			REPORT( DETAILED,  levelSignalList[j]->getName() << " size="<<levelSignalList[j]->getSize() << " weight=" << levelSignalList[j]->getWeight() );
		}	
			
		/* Gappa Style */
		REPORT(DETAILED, "======================================================");
		REPORT(DETAILED, "======================================================");
		REPORT(DETAILED, "======================================================");

		vector<int> yGuard(100); // negative values => truncation
		//initialization of the yGuard bits; these bits should evolve
		for (uint32_t i=1; i<=unsigned(degree); i++)
			yGuard[i]= 0;
		
		vector<int> aGuard(100); // positive values refer to "real" guard bits
		for (uint32_t i=0; i<unsigned(degree); i++)
			aGuard[i] = 1;
		
		/* the abs of the maximal value of y */
		mpfr_t maxABSy;
		mpfr_init2 ( maxABSy, 1000);
		mpfr_set_ui( maxABSy, 2, GMP_RNDN);
		mpfr_pow_si( maxABSy, maxABSy, y->getSize(), GMP_RNDN);
		mpfr_add_si( maxABSy, maxABSy, -1, GMP_RNDN);
		mpfr_set_exp( maxABSy, mpfr_get_exp(maxABSy)+y->getWeight()-y->getSize());
		REPORT(DETAILED, "Abs max value of y is " << mpfr_get_d( maxABSy, GMP_RNDN)); 
		
			
		vector<mpfr_t*> ykT_y(100); //the yk tild. The max absolute value of ykT
		for (uint32_t i=1; i<= unsigned(degree); i++){
			mpfr_t *cykT;
			cykT = (mpfr_t*) malloc(sizeof(mpfr_t));
			mpfr_init2(*cykT, 1000);
			if ( yGuard[i] < 0){
				mpfr_set_ui(*cykT, 2, GMP_RNDN);
				mpfr_pow_si(*cykT, *cykT, y->getWeight()-(signed(y->getSize())+yGuard[i]), GMP_RNDN); 
			}else
				mpfr_set_si(*cykT, 0, GMP_RNDN);
			ykT_y[i] = cykT;
		}
		
		for (uint32_t i=1; i<= unsigned(degree); i++)
			REPORT(DETAILED, "|y"<<i<<"T-y|="<< mpfr_get_d(*ykT_y[i],GMP_RNDN));


		vector<mpfr_t*> pikPT_pikP(100); //the pi k prime tild - p k prime ( trunc error)
		for (uint32_t i=1; i< unsigned(degree); i++){ //not needed of n
			mpfr_t *cpikPT_pikP;
			cpikPT_pikP = (mpfr_t*) malloc(sizeof(mpfr_t));
			mpfr_init2(*cpikPT_pikP, 1000);
	
			mpfr_set_ui(*cpikPT_pikP, 2, GMP_RNDN);
			mpfr_pow_si(*cpikPT_pikP, *cpikPT_pikP, coef[degree-i]->getWeight()-(signed(coef[degree-i]->getSize())+aGuard[degree-i]), GMP_RNDN); 
			pikPT_pikP[i] = cpikPT_pikP; //.push_back(cpikPT_pikP);
		}
		
		for (uint32_t i=1; i< unsigned(degree); i++)
			REPORT(DETAILED, "|pi"<<i<<"TP - pi"<<i<<"T|="<< mpfr_get_exp(*pikPT_pikP[i]));
		
		
		vector<mpfr_t*> sigmakP_sigmak(100);
		//initialize sigma0P_sigma0
		mpfr_t *sigmanP_sigman;
		sigmanP_sigman =(mpfr_t*) malloc( sizeof( mpfr_t));
		mpfr_init2 ( *sigmanP_sigman, 1000);
		mpfr_set_ui( *sigmanP_sigman, 0, GMP_RNDN);
		sigmakP_sigmak[0] = sigmanP_sigman;		

		vector<mpfr_t*> pikP_pik(100);
		vector<mpfr_t*> sigmakP(100);

		vector<mpfr_t*> a(100);
		for (uint32_t i=0; i<=unsigned(degree);i++){
			mpfr_t *ak;
			ak =(mpfr_t*) malloc( sizeof( mpfr_t));
			mpfr_init2 ( *ak, 1000);
			mpfr_set_ui( *ak, 2, GMP_RNDN);
			mpfr_pow_si( *ak, *ak, coef[i]->getSize(), GMP_RNDN);
			mpfr_add_si( *ak, *ak , -1, GMP_RNDN);
			mpfr_set_exp( *ak, (mpfr_get_d(*ak, GMP_RNDZ)!=0?mpfr_get_exp(*ak):0) + coef[i]->getWeight() - coef[i]->getSize());
			a[i]=ak;
		}
		for (uint32_t i=0; i<= unsigned(degree); i++)
			REPORT(DETAILED, "a"<<i<<"="<< mpfr_get_exp(*a[i]));
		
		sigmakP[0] = a[degree];	

		mpfr_t *yy;
		yy =(mpfr_t*) malloc( sizeof( mpfr_t));
		mpfr_init2 ( *yy, 1000);
		mpfr_set_ui( *yy, 2, GMP_RNDN);
		mpfr_pow_si( *yy, *yy, y->getSize(), GMP_RNDN);
		mpfr_add_si( *yy, *yy , -1, GMP_RNDN);
		mpfr_set_exp( *yy, ( mpfr_get_d(*yy, GMP_RNDZ)!=0?mpfr_get_exp(*yy):0) + y->getWeight() - y->getSize());

		REPORT(DETAILED, "y="<< mpfr_get_exp(*yy));

		vector<mpfr_t*> pikPT(100);

		vector<unsigned> sigmakPSize(100),   pikPTSize(100);
		vector<signed> sigmakPWeight(100), pikPTWeight(100);

		sigmakPSize[0]   = coef[degree]->getSize();
		sigmakPWeight[0] = coef[degree]->getWeight();

		for (uint32_t i=1; i<=unsigned(degree); i++){
			mpfr_t *t; //the pi
			t = (mpfr_t *) malloc( sizeof( mpfr_t ));			
			mpfr_init2( *t, 100);
			mpfr_mul ( *t, *ykT_y[i], *sigmakP[i-1], GMP_RNDN);

			REPORT(DEBUG, "(ykT-y)sigma(k-1)P exp" << mpfr_get_exp(*t) << " double=" << mpfr_get_d(*t, GMP_RNDN)); 

			mpfr_t *t2;
			t2 = (mpfr_t *) malloc (sizeof (mpfr_t));
			mpfr_init2( *t2, 1000);
			mpfr_mul ( *t2 , *yy , *sigmakP_sigmak[i-1], GMP_RNDN);

			REPORT(DEBUG, "y(sigma(k-1)P-sigma(k-1) exp" << mpfr_get_exp(*t2) << " double=" << mpfr_get_d(*t2, GMP_RNDN)); 
			
			mpfr_add (*t, *t , *t2, GMP_RNDN);
			pikP_pik[i] = t;

			REPORT( DEBUG, "----->|pikP"<<i<<"-pik| exp" << mpfr_get_exp( *pikP_pik[i]) << " double="<< mpfr_get_d( *pikP_pik[i], GMP_RNDN) );

			pikPTSize[i]   = (y->getSize()+yGuard[i]) + sigmakPSize[i-1];
			pikPTWeight[i] =  y->getWeight() + sigmakPWeight[i-1]; 			

			//compute value of pi k P T
			mpfr_t *h;
			h = (mpfr_t *) malloc( sizeof(mpfr_t));
			mpfr_init2(*h, 1000);
			mpfr_set_ui( *h, 2, GMP_RNDN);
			mpfr_pow_si( *h, *h, pikPTSize[i], GMP_RNDN);
			mpfr_add_si( *h, *h , -1, GMP_RNDN);
			mpfr_set_exp( *h, (mpfr_get_d(*h, GMP_RNDZ)!=0? mpfr_get_exp(*h):0) + pikPTWeight[i] - pikPTSize[i]);

			pikPT[i] = h;
			////////////////////////////////////////////////////////////////////

			mpfr_t *t3; //the sigma
			if (i < unsigned(degree)){
				t3 = (mpfr_t *) malloc( sizeof( mpfr_t ));
				mpfr_init2( *t3, 1000);
				mpfr_add( *t3, *pikPT_pikP[i], *pikP_pik[i],GMP_RNDN); 		   
			}else{
				t3 = (mpfr_t *) malloc( sizeof( mpfr_t ));
				mpfr_init2( *t3, 1000);
				mpfr_set( *t3, *pikP_pik[i], GMP_RNDN); 		   
			}
			sigmakP_sigmak[i] = t3;

			REPORT( DEBUG, "----->|sigmakP"<<i<<"-sigmak| exp" << mpfr_get_exp( *sigmakP_sigmak[i]) << " double=" << mpfr_get_d( *sigmakP_sigmak[i], GMP_RNDN));

		
			int maxMSB = max ( coef[degree-i]->getWeight(), pikPTWeight[i] ) + 1;
			sigmakPSize[i]   = maxMSB - (coef[degree-i]->getWeight() - coef[degree-i]->getSize() - aGuard[degree-i]) +1;
			sigmakPWeight[i] = maxMSB + 1; 

			mpfr_t *r;
			r = (mpfr_t *) malloc( sizeof(mpfr_t));
			mpfr_init2( *r, 1000);
			mpfr_add ( *r, *pikPT[i], *a[degree-i], GMP_RNDN);			

			sigmakP[i] = r; 
		}
		
		REPORT( DETAILED, "Error (order) for P(y)=" << mpfr_get_exp( *sigmakP_sigmak[degree])+1);
		
	}

	PolynomialEvaluator::~PolynomialEvaluator() {
	}

}



