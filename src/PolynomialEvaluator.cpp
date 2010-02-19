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

		degree_ = coef.size()-1;
		setName( join("PolynomialEvaluator_d",degree_) );
		
		REPORT(DETAILED, "Polynomial to evaluate: " << printPolynomial(coef, y));
		
		/* ================== I/O declarations ==================*/
		addInput("X", y_->getSize());
		for (uint32_t i=0; i <= unsigned(degree_); i++)
			addInput(join("C",i), coef_[i]->getSize() );

		/* Gappa Style */
		
		yGuard_.reserve(100);
		aGuard_.reserve(100);
		nYGuard_.reserve(100);
		
		maxBoundY = 5;
		maxBoundA = 10;
		
		/*init vectors */
		for (uint32_t i=1; i<=unsigned(degree_)+1; i++){
			yGuard_[i] = 0; //maxBoundY;
//			nYGuard[i] = 0;
		}
		
		for (uint32_t i=0; i<unsigned(degree_)+1; i++)
			aGuard_[i] = 0;
		
		/* the abs of the maximal value of y */
		
		mpfr_init2 ( maxABSy, 1000);
		mpfr_set_ui( maxABSy, 2, GMP_RNDN);
		mpfr_pow_si( maxABSy, maxABSy, y_->getSize(), GMP_RNDN);
		mpfr_add_si( maxABSy, maxABSy, -1, GMP_RNDN);
		mpfr_set_exp( maxABSy, mpfr_get_exp(maxABSy)+y_->getWeight()-y_->getSize());
		REPORT(DETAILED, "Abs max value of y is " << mpfr_get_d( maxABSy, GMP_RNDN)); 
		
		
	   sol = false;
		
		while (!sol){
			while ((nextStateA()) && (!sol)){
				while ( (nextStateY())&& (!sol) ){
					currentPrec = errorEstimator(yGuard_, aGuard_);
					if (currentPrec <= targetPrec ){
						sol = true;
					}
				} 
			}	
		}		
		
		REPORT(DETAILED, "yuppy " << errorEstimator(yGuard_, aGuard_));
		
	}
	
	
	int PolynomialEvaluator::errorEstimator(vector<int> &yGuard, vector<int> &aGuard){
//		cout << "asdasd";
		
		
		for (int j=0; j<=degree_; j++)
			cout << "ga["<<j<<"]="<<aGuard[j]<<" "; 
		cout << endl;				

		for (int j=1; j<=degree_; j++)
			cout << "gy["<<j<<"]="<<yGuard[j]<<" "; 
		cout << endl;	
	
	

		vector<mpfr_t*> ykT_y(100); //the yk tild. The max absolute value of ykT
		for (uint32_t i=1; i<= unsigned(degree_); i++){
			mpfr_t *cykT;
			cykT = (mpfr_t*) malloc(sizeof(mpfr_t));
			mpfr_init2(*cykT, 1000);
			if ( yGuard[i] < 0){
				mpfr_set_ui(*cykT, 2, GMP_RNDN);
				mpfr_pow_si(*cykT, *cykT, y_->getWeight()-(signed(y_->getSize())+yGuard[i]), GMP_RNDN); 
			}else
				mpfr_set_si(*cykT, 0, GMP_RNDN);
			ykT_y[i] = cykT;
		}
		
		for (uint32_t i=1; i<= unsigned(degree_); i++)
			REPORT(DETAILED, "|y"<<i<<"T-y|="<< mpfr_get_d(*ykT_y[i],GMP_RNDN));


		vector<mpfr_t*> pikPT_pikP(100); //the pi k prime tild - p k prime ( trunc error)
		for (uint32_t i=1; i< unsigned(degree_); i++){ //not needed of n
			mpfr_t *cpikPT_pikP;
			cpikPT_pikP = (mpfr_t*) malloc(sizeof(mpfr_t));
			mpfr_init2(*cpikPT_pikP, 1000);
	
			mpfr_set_ui(*cpikPT_pikP, 2, GMP_RNDN);
			mpfr_pow_si(*cpikPT_pikP, *cpikPT_pikP, coef_[degree_-i]->getWeight()-(signed(coef_[degree_-i]->getSize())+aGuard[degree_-i]), GMP_RNDN); 
			pikPT_pikP[i] = cpikPT_pikP; //.push_back(cpikPT_pikP);
		}
		
		for (uint32_t i=1; i< unsigned(degree_); i++)
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
		for (uint32_t i=0; i<=unsigned(degree_);i++){
			mpfr_t *ak;
			ak =(mpfr_t*) malloc( sizeof( mpfr_t));
			mpfr_init2 ( *ak, 1000);
			mpfr_set_ui( *ak, 2, GMP_RNDN);
			mpfr_pow_si( *ak, *ak, coef_[i]->getSize(), GMP_RNDN);
			mpfr_add_si( *ak, *ak , -1, GMP_RNDN);
			mpfr_set_exp( *ak, (mpfr_get_d(*ak, GMP_RNDZ)!=0?mpfr_get_exp(*ak):0) + coef_[i]->getWeight() - coef_[i]->getSize());
			a[i]=ak;
		}
		for (uint32_t i=0; i<= unsigned(degree_); i++)
			REPORT(DETAILED, "a"<<i<<"="<< mpfr_get_exp(*a[i]));
		
		sigmakP[0] = a[degree_];	

		mpfr_t *yy;
		yy =(mpfr_t*) malloc( sizeof( mpfr_t));
		mpfr_init2 ( *yy, 1000);
		mpfr_set_ui( *yy, 2, GMP_RNDN);
		mpfr_pow_si( *yy, *yy, y_->getSize(), GMP_RNDN);
		mpfr_add_si( *yy, *yy , -1, GMP_RNDN);
		mpfr_set_exp( *yy, ( mpfr_get_d(*yy, GMP_RNDZ)!=0?mpfr_get_exp(*yy):0) + y_->getWeight() - y_->getSize());

		REPORT(DETAILED, "y="<< mpfr_get_exp(*yy));

		vector<mpfr_t*> pikPT(100);

		vector<unsigned> sigmakPSize(100),   pikPTSize(100);
		vector<signed> sigmakPWeight(100), pikPTWeight(100);

		sigmakPSize[0]   = coef_[degree_]->getSize();
		sigmakPWeight[0] = coef_[degree_]->getWeight();

		for (uint32_t i=1; i<=unsigned(degree_); i++){
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

			pikPTSize[i]   = (y_->getSize()+yGuard[i]) + sigmakPSize[i-1];
			pikPTWeight[i] =  y_->getWeight() + sigmakPWeight[i-1]; 			

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
			if (i < unsigned(degree_)){
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

		
			int maxMSB = max ( coef_[degree_-i]->getWeight(), pikPTWeight[i] ) + 1;
			sigmakPSize[i]   = maxMSB - (coef_[degree_-i]->getWeight() - coef_[degree_-i]->getSize() - aGuard[degree_-i]) +1;
			sigmakPWeight[i] = maxMSB + 1; 

			mpfr_t *r;
			r = (mpfr_t *) malloc( sizeof(mpfr_t));
			mpfr_init2( *r, 1000);
			mpfr_add ( *r, *pikPT[i], *a[degree_-i], GMP_RNDN);			

			sigmakP[i] = r; 
		}
		
		REPORT( DETAILED, "Error (order) for P(y)=" << mpfr_get_exp( *sigmakP_sigmak[degree_])+1);
		return mpfr_get_exp( *sigmakP_sigmak[degree_])+1;
	} 



	PolynomialEvaluator::~PolynomialEvaluator() {
	}

}



