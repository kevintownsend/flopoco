#ifndef ATANPOW_HPP
#define ATANPOW_HPP

#include "../Operator.hpp"

struct AtanPow : Operator
{
	AtanPow(Target * target, int wE, int wF, int o);
	virtual ~AtanPow();

	virtual void outputVHDL(std::ostream& o, std::string name);

	int wE;
	int wF;
	int order;
private:
	
	Operator * T[3];
};

#endif
