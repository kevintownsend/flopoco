#ifndef TableGenerator_HPP
#define TableGenerator_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "Operator.hpp"
#include "Table.hpp"
#include "PolynomialEvaluator.hpp"
#include "HOTBM/sollya.h"	// Do NOT use libsollya from user's environment

#include "HOTBM/Function.hh"
#include "HOTBM/MPPolynomial.hh"
#include "UtilSollya.hh"
#include "PiecewiseFunction.hh"

namespace flopoco{

	/** The TableGenerator class.  */
	class TableGenerator : public Table {

		public:
       TableGenerator(Target* target, string func, int wInX, int wOutX, int n);
			 /* TODO: Doxygen parameters*/ 
			TableGenerator(Target* target, string func, int wInX, int wOutX, int n,double xmin, double xmax, double scale);

			/**
			 * TableGenerator destructor
			 */
			~TableGenerator();
			
			MPPolynomial* getMPPolynomial(sollya_node_t t);
			vector<FixedPointCoefficient*> getPolynomialCoefficients(sollya_node_t t, sollya_chain_t c);
			vector<FixedPointCoefficient*> getPolynomialCoefficients(sollya_node_t t, int* sizeList);
      vector<vector<FixedPointCoefficient*> > getPolynomialCoefficientsVector();
      void printPolynomialCoefficientsVector();
      void updateMinWeightParam(int i, FixedPointCoefficient* zz);
      vector<FixedPointCoefficient*> getCoeffParamVector();
      void printCoeffParamVector();
      mpfr_t *getMaxApproxError();
      void generateDebug();
      void generateDebugPwf();
      sollya_chain_t makeIntPtrChainCustomized(int m, int n, int precshift, int msize);
      /************************************************/
      /********Virtual methoods from class Table*******/
      mpz_class function(int x);
			int    double2input(double x);
			double input2double(int x);
			mpz_class double2output(double x);
			double output2double(mpz_class x);
	    /************************************************/
      
      
		protected:
			int wInX_;   /**< TODO: Description*/ 
			int wOutX_;  /**< TODO: Description*/
			Function *f;
			vector< vector<FixedPointCoefficient*> > polyCoeffVector;
			vector<FixedPointCoefficient*> coeffParamVector;
			mpfr_t *maxError;
      PiecewiseFunction *pwf;
	};
}
#endif
