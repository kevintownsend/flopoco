#ifndef LNSMUL_HPP
#define LNSMUL_HPP

#include "../Operator.hpp"

namespace flopoco{
	struct LNSMul : Operator
	{
		LNSMul(Target * target, int wE, int wF);
		virtual ~LNSMul();

		virtual void outputVHDL(std::ostream& o, std::string name);
		virtual void setOperatorName();	
	
	private:
		int wE;
		int wF;
	};
}
#endif
