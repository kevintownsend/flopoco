#ifndef _GENERICEVALUATOR_HH_
#define _GENERICEVALUATOR_HH_

#include "../Operator.hpp"

namespace flopoco {

enum EvaluationMethod
{
	Hotbm,
	Polynomial,
};

Operator * NewEvaluator(Target * target, char const * func, std::string const & uniqueName,
	int wI, int wO, int o, double xmin, double xmax, double scale,
	EvaluationMethod method = Hotbm);
	
}

#endif
