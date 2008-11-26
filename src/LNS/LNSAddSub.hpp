#ifndef LNSADDSUB_HPP
#define LNSADDSUB_HPP

#include "../Operator.hpp"
#include "LNSAdd.hpp"
#include "CotranHybrid.hpp"

struct LNSAddSub : Operator
{
	LNSAddSub(Target * target, int wE, int wF);
	virtual ~LNSAddSub();

	virtual void outputVHDL(std::ostream& o, std::string name);
	virtual void setOperatorName();	
	
private:
	int wE;
	int wF;
	int j;
	
	CotranHybrid * addsub;
};

#endif
