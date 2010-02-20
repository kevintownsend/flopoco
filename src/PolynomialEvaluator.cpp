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

	extern vector<Operator*> oplist;
	
		PolynomialEvaluator::PolynomialEvaluator(Target* target, vector<FixedPointCoefficient*> coef, YVar* y, int targetPrec):
		Operator(target), coef_(coef), y_(y), targetPrec_(targetPrec){
		srcFileName = "PolynomialEvaluator";

		degree_ = coef.size()-1;
		setName( join("PolynomialEvaluator_d",degree_) );
		
		REPORT(DETAILED, "Polynomial to evaluate: " << printPolynomial(coef, y));
		
		/* ================== I/O declarations ==================*/
		addInput("Y", y_->getSize());
		for (uint32_t i=0; i <= unsigned(degree_); i++)
			addInput(join("a",i), coef_[i]->getSize()+1);

		/* Gappa Style */
		
		yGuard_.reserve(100);
		aGuard_.reserve(100);
		nYGuard_.reserve(100);
		
		maxBoundY = 5;
		maxBoundA = 5;
		
		/*init vectors */
		for (uint32_t i=1; i<=unsigned(degree_)+1; i++){
			yGuard_[i] = 0; //maxBoundY;
			nYGuard_[i] = 0;
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
		

		/* design space exploration */				
		sol = false;
		while (!sol){
			while (((!sol) && nextStateA()) ){
				while ((!sol) && (nextStateY())  ){
					currentPrec = errorEstimator(yGuard_, aGuard_);
					if (currentPrec <= targetPrec ){
						sol = true;
					}
				} 
			}	
		}		
		
		REPORT(DETAILED, "yuppy " << errorEstimator(yGuard_, aGuard_));
		
		for (uint32_t i=0; i<=unsigned(degree_); i++){
			if (i==0)
				vhdl << tab << declare( join("sigmaP",i), 1+coef_[degree_-i]->getSize()) << " <= a"<<degree_<<";"<<endl; 
			else{
				vhdl << tab << declare( join("yT",i) , 1+y_->getSize()+yGuard_[i]) << " <= \"0\" & Y"<<range(y_->getSize()-1, -yGuard_[i]) << ";" << endl;
				SignedIntMultiplier* sm = new SignedIntMultiplier ( target,     1+y_->getSize()+yGuard_[i], sigmakPSize[i-1]+1);
				oplist.push_back(sm);
				
				inPortMap ( sm, "X", join("yT",i));
				inPortMap ( sm, "Y", join("sigmaP",i-1));
				outPortMap (sm, "R", join("piP",i));
				
				vhdl << instance ( sm, join("Product_",i) );
				syncCycleFromSignal(join("piP",i)); 
				nextCycle();///////////////////////////////////////////////////
				
				if (i<unsigned(degree_)){
					vhdl << tab << declare( join("piPT",i), pikPTSize[i]+1 ) << " <= " << join("piP",i)<<range(pikPSize[i], pikPSize[i] - pikPTSize[i] ) << ";" <<endl; // coef_[i]->getSize()+1+y_->getSize()+yGuard_[i]-1, coef_[i]->getSize()+1+y_->getSize()+yGuard_[i]-1 - pikPTSize[i]) << ";" << endl;   
					IntAdder* sa = new IntAdder (target, sigmakPSize[i]+1);
					oplist.push_back(sa);
				
				
					vhdl << tab << declare( join("op1_",i), sigmakPSize[i]+1 ) << " <= " << "(" << rangeAssign(sigmakPWeight[i] - coef_[degree_-i]->getWeight()-1,0, join("a",degree_-i)+of(coef_[degree_-i]->getSize()))
						                                                               << " & " << join("a",degree_-i) << " & "<< zg(aGuard_[degree_-i],0) << ");"<<endl;
						                                                               
					vhdl << tab << declare( join("op2_",i), sigmakPSize[i]+1 ) << " <= " << "(" << rangeAssign(sigmakPWeight[i]-pikPTWeight[i]-1,0, join("piPT",i)+of(pikPTSize[i])) 
						                                                               << " & " << join("piPT",i) << range(pikPTSize[i], pikPTSize[i] - pikPTWeight[i] - (coef_[degree_-i]->getSize()-coef_[degree_-i]->getWeight() + aGuard_[degree_ -i]))
						                                                               << ");" << endl;

					inPortMap ( sa, "X", join("op1_",i) );
					inPortMap ( sa, "Y", join("op2_",i) );
					inPortMapCst ( sa, "Cin", "'0'");
					outPortMap( sa, "R", join("sigmaP",i));
				
					vhdl << instance ( sa, join("Sum",i));
					syncCycleFromSignal( join("sigmaP",i) );
					nextCycle();///////////////////////////////////////////////////                                                                   
				                                                                   
				}else{
					IntAdder* sa = new IntAdder (target, sigmakPSize[i]+1);
					oplist.push_back(sa);

					vhdl << tab << declare( join("op1_",i), sigmakPSize[i]+1 ) << " <= "<< "(" << rangeAssign(sigmakPWeight[i]-pikPWeight[i],0, join("piP",i)+of(pikPSize[i]+1)) 
						                                                               << " & " << join("piP",i)<<range(pikPSize[i],0)  << ");" << endl;

					vhdl << tab << declare( join("op2_",i), sigmakPSize[i]+1 ) << " <= " << "(" << rangeAssign(sigmakPWeight[i]-coef_[degree_-i]->getWeight()-1,0, join("a",degree_-i)+of(coef_[degree_-i]->getSize()))
						                                                               << " & " << join("a",degree_-i) << " & "<< zg(sigmakPSize[i]-sigmakPWeight[i]-(coef[degree_-i]->getSize()-coef[degree_-i]->getWeight()) ,0) << ");"<<endl;
				
					inPortMap ( sa, "X", join("op1_",i) );
					inPortMap ( sa, "Y", join("op2_",i) );
					inPortMapCst ( sa, "Cin", "'0'");
					outPortMap( sa, "R", join("sigmaP",i));
		
					vhdl << instance ( sa, join("Sum",i));
					syncCycleFromSignal( join("sigmaP",i) );


					addOutput("R", sigmakPSize[i]+1);
					vhdl << tab << "R <= " << join("sigmaP",i) << ";" << endl;
				}			
			}		
		}
	}
	
	
	int PolynomialEvaluator::errorEstimator(vector<int> &yGuard, vector<int> &aGuard){
		ostringstream s1, s2;		
		
		for (int j=0; j<=degree_; j++)
			s1 << "aG["<<j<<"]="<<aGuard[j]<<" "; 

		for (int j=1; j<=degree_; j++)
			s2 << "yG["<<j<<"]="<<yGuard[j]<<" "; 

		REPORT(DETAILED, "------------------------------------------------------------");
		REPORT(DETAILED, s1.str());
		REPORT(DETAILED, s2.str());
		REPORT(DETAILED, "---------- HEADER ----------");
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
			REPORT(DETAILED, "degree of |y"<<i<<"T - y|="<< (mpfr_get_d(*ykT_y[i],GMP_RNDN)!=0?mpfr_get_exp(*ykT_y[i]):0));


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
			REPORT(DETAILED, "|pi"<<i<<"TP - pi"<<i<<"T|="<< (mpfr_get_d(*pikPT_pikP[i],GMP_RNDN)!=0? mpfr_get_exp(*pikPT_pikP[i]):0));
		
		
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
			REPORT(DETAILED, "weight of a"<<i<<"="<< mpfr_get_exp(*a[i]));
		
		sigmakP[0] = a[degree_];	

		mpfr_t *yy;
		yy =(mpfr_t*) malloc( sizeof( mpfr_t));
		mpfr_init2 ( *yy, 1000);
		mpfr_set_ui( *yy, 2, GMP_RNDN);
		mpfr_pow_si( *yy, *yy, y_->getSize(), GMP_RNDN);
		mpfr_add_si( *yy, *yy , -1, GMP_RNDN);
		mpfr_set_exp( *yy, ( mpfr_get_d(*yy, GMP_RNDZ)!=0?mpfr_get_exp(*yy):0) + y_->getWeight() - y_->getSize());

		REPORT(DETAILED, "weight of y="<< mpfr_get_exp(*yy));

		vector<mpfr_t*> pikPT(100);

		sigmakPSize.clear();
		pikPTSize.clear();
		pikPSize.clear();
		sigmakPWeight.clear();
		pikPTWeight.clear();
		pikPWeight.clear();
		
		sigmakPSize.reserve(100);
		pikPTSize.reserve(100);
		pikPSize.reserve(100);
		sigmakPWeight.reserve(100);
		pikPTWeight.reserve(100);
		pikPWeight.reserve(100);


		sigmakPSize[0]   = coef_[degree_]->getSize();
		sigmakPWeight[0] = coef_[degree_]->getWeight();
		
		REPORT(DETAILED, "----------   END  ----------");
		for (uint32_t i=1; i<=unsigned(degree_); i++){
			mpfr_t *t; //the pi
			t = (mpfr_t *) malloc( sizeof( mpfr_t ));			
			mpfr_init2( *t, 100);
			mpfr_mul ( *t, *ykT_y[i], *sigmakP[i-1], GMP_RNDN);

			REPORT(DETAILED, "|(y"<<i<<"T - y)*sigma"<<i-1<<"P| weight=" << (mpfr_get_d(*t,GMP_RNDN)!=0?mpfr_get_exp(*t):0));

			pikPWeight[i]  = y_->getWeight() + sigmakPWeight[i-1]; 
			pikPSize[i]   = (y_->getSize()+yGuard[i]) + sigmakPSize[i-1];

			mpfr_t *t2;
			t2 = (mpfr_t *) malloc (sizeof (mpfr_t));
			mpfr_init2( *t2, 1000);
			mpfr_mul ( *t2 , *yy , *sigmakP_sigmak[i-1], GMP_RNDN);

			REPORT(DETAILED, "|y*(sigma"<<i-1<<"P - sigma"<<i-1<<")| weight=" <<  (mpfr_get_d(*t2,GMP_RNDN)!=0?mpfr_get_exp(*t2):0));
			
			mpfr_add (*t, *t , *t2, GMP_RNDN);
			pikP_pik[i] = t;

			REPORT( DETAILED, "|pi"<<i<<"P - pi"<<i<<"| weight=" <<  (mpfr_get_d(*pikP_pik[i],GMP_RNDN)!=0?mpfr_get_exp(*pikP_pik[i]):0)); 

			
			pikPTWeight[i] =  y_->getWeight() + sigmakPWeight[i-1]; 
			pikPTSize[i] =    y_->getWeight() + sigmakPWeight[i-1]  + (coef_[degree_-i]->getSize()+aGuard[i]-coef_[degree_-i]->getWeight());
	
	
			REPORT( DEBUG, " pikPTSize="<< pikPTSize[i] << " pikPTWeight[i]="<<pikPTWeight[i]);
			REPORT( DEBUG, "pikPTWeight[i]"<<pikPTWeight[i]<<" coef_[degree_-i]->getSize() "<< coef_[degree_-i]->getSize() << " coef_[degree_-i]->getWeight() "<<coef_[degree_-i]->getWeight() << " aGuard =" << aGuard[i] ); 

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

			REPORT( DETAILED, "|sigma"<<i<<"P - sigma"<<i<<"| weight=" << (mpfr_get_d(*sigmakP_sigmak[i],GMP_RNDN)!=0?mpfr_get_exp(*sigmakP_sigmak[i]):0)); 

			if (i!=unsigned(degree_)){
				int maxMSB = max ( coef_[degree_-i]->getWeight(), pikPTWeight[i]);
				sigmakPSize[i]   = maxMSB + 1 - (coef_[degree_-i]->getWeight() - coef_[degree_-i]->getSize() - aGuard[degree_-i]);
				sigmakPWeight[i] = maxMSB + 1;
			
				REPORT( DEBUG, " sigmakPSize="<< sigmakPSize[i] << " sigmakPWeight[i]="<<sigmakPWeight[i]);
			}else{
				int maxMSB = max ( coef_[degree_-i]->getWeight(), pikPWeight[i]);
				int minLSB = min ( coef_[degree_-i]->getWeight()-coef_[degree_-i]->getSize(),pikPWeight[i]-pikPSize[i]); 
				sigmakPSize[i]   = 1 + maxMSB - minLSB + 1;//maxMSB + 1 - (coef_[degree_-i]->getWeight() - coef_[degree_-i]->getSize() - aGuard[degree_-i]);
				sigmakPWeight[i] = 1 + maxMSB;
			
				REPORT( DEBUG, " sigmakPSize="<< sigmakPSize[i] << " sigmakPWeight[i]="<<sigmakPWeight[i]);
			}
			mpfr_t *r;
			r = (mpfr_t *) malloc( sizeof(mpfr_t));
			mpfr_init2( *r, 1000);
			mpfr_add ( *r, *pikPT[i], *a[degree_-i], GMP_RNDN);			

			sigmakP[i] = r;
			REPORT( DETAILED, "|sigma"<<i<<"P| weight=" << (mpfr_get_d(*r,GMP_RNDN)!=0?mpfr_get_exp(*r):0)); 
		}
		
		REPORT( DETAILED, "Error (order) for P(y)=" << mpfr_get_exp( *sigmakP_sigmak[degree_])+1);
		return mpfr_get_exp( *sigmakP_sigmak[degree_])+1;
	} 



	PolynomialEvaluator::~PolynomialEvaluator() {
	}

}



