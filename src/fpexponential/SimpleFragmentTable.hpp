#ifndef SIMPLEFRAGMENTTABLE_HPP
#define SIMPLEFRAGMENTTABLE_HPP

#include "../Table.hpp"
#include "Fragment.hpp"

class SimpleFragmentTable : public Table
{
public:
	SimpleFragmentTable(Target* target, Fragment* fragment);

	~SimpleFragmentTable();

	mpz_class function(int x);

	Fragment * fragment; /**< holds all the useful information about size, accuracy, etc */
	
#if 0
	int    double2input(double x);
	
	double input2double(int x);
	
	mpz_class double2output(double x);
	
	double output2double(mpz_class x);
#endif
};






#endif //STDFRAGMENTEXPTABLE_HPP

