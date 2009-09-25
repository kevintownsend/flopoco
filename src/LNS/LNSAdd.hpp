#ifndef LNSADD_HPP
#define LNSADD_HPP

#include "../Operator.hpp"


namespace flopoco{

	// LNSAdd: implements f+(z) = log2(1 + 2^z)
	// z in fixed-point, wE integral bits, wF fractional bits
	struct LNSAdd : Operator
	{
		LNSAdd(Target * target, int wE, int wF, int o);
		virtual ~LNSAdd();

		virtual void outputVHDL(std::ostream& o, std::string name);

		//virtual void fillTestCase(mpz_class a[]);

		int wE;
		int wF;
		int order;
	private:
	
		Operator * t[3];
	};

}
#endif
