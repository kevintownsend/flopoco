#ifndef COTRAN_HPP
#define COTRAN_HPP

#include "../Operator.hpp"
#include <memory>
#include "CotranTables.hpp"
#include "LNSAdd.hpp"

struct Cotran : Operator
{
	Cotran(Target * target, int wE, int wF, int j = -1, int wECotran = -1, int o = 1);
	virtual ~Cotran();

	virtual void outputVHDL(std::ostream& o, std::string name);

	//virtual TestIOMap getTestIOMap();
	//virtual void fillTestCase(mpz_class a[]);

private:
	int wE;
	int wF;
	int j;
	int wECotran;
	
	void select_j();
	
	std::auto_ptr<CotranF1Table> f1;
	std::auto_ptr<CotranF2Table> f2;
	std::auto_ptr<CotranF3Table> f3;
	std::string f1_name;
	std::string f2_name;
	std::string f3_name;
	
	LNSAdd * sb;
	
	int wEssZero;
	int DBMaxInput;
};


#endif

