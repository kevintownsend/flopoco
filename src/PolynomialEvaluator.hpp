#ifndef PolynomialEvaluator_HPP
#define PolynomialEvaluator_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"

#include "SignedIntMultiplier.hpp"
#include "IntAdder.hpp"
#include "Operator.hpp"


namespace flopoco{
#define MINF -16384

	extern map<string, double> emptyDelayMap;
	
	/**
	 * The FixedPointCoefficinet class provides the objects for representing 
	 * fixed-point coefficinets. It is used by the polynomial evaluator
	 */
	class FixedPointCoefficient
	{
		public:
		
			FixedPointCoefficient(FixedPointCoefficient *x){
				size_= x->getSize();
				weight_ = x->getWeight();
			}

			/**
			 * Constructor. Initializes the attributes.
			 */
			FixedPointCoefficient(unsigned size, int weight){
				size_   = size;
				weight_  = weight;
			}

      		FixedPointCoefficient(unsigned size, int weight, mpfr_t valueMpfr){
				size_   = size;
				weight_  = weight;
				valueMpfr_=(mpfr_t*)malloc(sizeof(mpfr_t));
				mpfr_init2(*valueMpfr_, mpfr_get_prec(valueMpfr));
				mpfr_set(*valueMpfr_, valueMpfr, GMP_RNDN);
			}
			
			/** Destructor */
			~FixedPointCoefficient();

			/** 
			 * Fetch the weight of the MSB
			 * @return weight of the MSB
			 */
			int getWeight(){
				return weight_;
			} 

			/** 
			 * Fetch the size (width) of the coefficient
			 * @return coefficient width
			 */
			unsigned getSize(){
				return size_;
			}

			void setSize(unsigned s){
				size_ = s;
			}

			
			mpfr_t* getValueMpfr(){
				return valueMpfr_;
			}


		protected:
			int size_;        /**< The width (in bits) of the coefficient. In theory this equals to the
			                        minimum number of bits required to represent abs(value) */
			int weight_;      /**< The shift value. Coefficient is represented as value * 2 (shift) */
			
			mpfr_t* valueMpfr_;
	};

	/**
	 * The YVar class
	 */
	class YVar: public FixedPointCoefficient
	{
		public:
			/** constructor */	
			YVar( unsigned size, int weight, const int sign = 0 ):
				FixedPointCoefficient(size, weight), sign_(sign){
			}

			/** Destructor */
			~YVar(){};

			/**
			 * Fetch the variable size (if known). 0 = unknown
			 * @return sign
			 */ 
			int getSign(){
				return sign_;
			}
			
		protected:
			int sign_;        /**< The variable sign in known */         
	};
	
	/** class LevelSignal. Models The level signals */
	class LevelSignal{
		public: 
			LevelSignal(string name, unsigned size, int shift):
				name_(name), size_(size), shift_(shift){
			}
			
			LevelSignal(LevelSignal* l){
				name_  = l->getName();
				size_  = l->getSize();
				shift_ = l->getWeight();
			}
			
			/* Destructor */
			~LevelSignal(){};

			string getName(){
				return name_;
			}
			
			unsigned getSize(){
				return size_;
			}
			
			int getWeight(){
				return shift_;
			}
			
		protected:
			string   name_; /**TODO Comment */
			unsigned size_;
			int      shift_;
	}; 
	
	
	

	/** 
	 * The PolynomialEvaluator class. Constructs the oprator for evaluating
	 * a polynomial shiftth fix-point coefficients.
	 */
	class PolynomialEvaluator : public Operator
	{
		public:

			/**
			 * The polynomial evaluator class. FIXME the parameters are subhect to change
			 * TODO document them
			 */
			PolynomialEvaluator(Target* target, vector<FixedPointCoefficient*> coef, YVar* y, int targetPrec);

			/** Destructor */
			~PolynomialEvaluator();

			/** print the polynomial */
			string printPolynomial( vector<FixedPointCoefficient*> coef, YVar* y){
				ostringstream polyS;
				for(uint32_t i=0; i <= coef.size()-1; i++)
					polyS <<  "a"<<i<< "*2^("<< coef[i]->getWeight() <<")*"<< join("X^",i) << "*2^"<< y->getWeight();
				return polyS.str();
			}


			/** the absolute error caused by a truncation where the resulting
			    LSB is found at this index (trunkIndex)*/
			mpfr_t* getTruncError( int trunkIndex) {
				mpfr_t* tmp;
				mpfr_init2 ( *tmp, 500);
				mpfr_set_si( *tmp, 2 , GMP_RNDN);
				mpfr_pow_si( *tmp, *tmp, trunkIndex, GMP_RNDN);
				return tmp;				
			}


			int errorEstimator(vector<int> &yGuard, vector<int> &aGuard);

			bool nextStateY(){
				if (! sol){
					int carry = 1;
					for (int i=1; i<=degree_+1;i++){
						if ((nYGuard_[i] == maxBoundY) && ( carry==1)){
							nYGuard_[i] = 0;
							carry = 1;
						}else{
							nYGuard_[i]+=carry;
							carry=0;
						}
					}
				
					for (int i=1; i<=degree_; i++)
						yGuard_[i] = -(maxBoundY - nYGuard_[i]);
				
					if ( nYGuard_[degree_+1] != 0)
						return true;
					else
						return false;
				}else
					return false;
			}

			bool nextStateA(){
				if (! sol){
					int carry = 1;
					for (int i=0; i<=degree_;i++){
						if ((aGuard_[i] == maxBoundA) && ( carry==1)){
							aGuard_[i] = 0;
							carry = 1;
						}else{
							aGuard_[i]+=carry;
							carry=0;
						}
					}
					if ( aGuard_[degree_] != 0)
						return true;
					else
						return false;
				}
				else
					return false;
			}
			
			unsigned getOutputSize(){
				return wR;
			}
			
			vector<FixedPointCoefficient*> getCoeffParamVector(){
				return coef_;
			}

		protected:
			unsigned wR;
			
			vector<FixedPointCoefficient*> coef_;
			YVar* y_;
			int targetPrec_;
			int degree_;
			mpfr_t maxABSy;
			
			
			
			vector<int> yGuard_; // negative values => truncation
			vector<int> nYGuard_; 


			vector<int> aGuard_; // positive values refer to "real" guard bits

			int maxBoundY;
			int maxBoundA;
			
			int currentPrec;
			bool sol;
			
			
			vector<signed> sigmakPSize,  pikPSize, pikPTSize;
		    vector<signed> sigmakPWeight, pikPWeight, pikPTWeight;


	};
}

#endif
