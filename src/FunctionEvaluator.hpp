#ifndef FunctionEvaluator_HPP
#define FunctionEvaluator_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "TableGenerator.hpp"
#include "FunctionEvaluator.hpp"

namespace flopoco{

	/** The FunctionEvaluator class */
	class FunctionEvaluator : public Operator
	{
	public:
		/**
		 * The FunctionEvaluator constructor
		 */
		FunctionEvaluator(Target* target, string func, int wInX, int wOutX, int n,double xmin, double xmax, double scale);

		/**
		 * FunctionEvaluator destructor
		 */
		~FunctionEvaluator();

	private:
		
		TableGenerator *tg;
		PolynomialEvaluator *pe;
		
	};

}

#endif
