#ifndef PolynomialEvaluator_HPP
#define PolynomialEvaluator_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"

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
			/**
			 * Constructor. Initializes the attributes.
			 */
			FixedPointCoefficient(unsigned size, int weight){
				size_   = size;
				weight_  = weight;
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


		protected:
			int size_;        /**< The width (in bits) of the coefficient. In theory this equals to the
			                        minimum number of bits required to represent abs(value) */
			int weight_;      /**< The shift value. Coefficient is represented as value * 2 (shift) */
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

		protected:
			vector<FixedPointCoefficient*> coef_;
			YVar* y_;
			int targetPrec_;
	};
}

#endif
