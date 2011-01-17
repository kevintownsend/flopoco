#ifndef ATANPOW_HPP
#define ATANPOW_HPP

#include "../Operator.hpp"
#include "GenericEvaluator.hpp"

namespace flopoco{

	struct AtanPow : Operator
	{
		AtanPow(Target * target, int wE, int wF, int o, EvaluationMethod method = Hotbm);
		virtual ~AtanPow();

		//virtual void outputVHDL(std::ostream& o, std::string name);

		int wE;
		int wF;
		int order;
	private:
		Operator * NewInterpolator(int wI, int wO, int o, double xmin, double xmax, double scale);
	
		EvaluationMethod method_;
		Operator * T[3];
	};

}

#endif
