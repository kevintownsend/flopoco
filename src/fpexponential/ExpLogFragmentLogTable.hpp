#ifndef EXPLOGFRAGMENTLOGTABLE_HPP
#define EXPLOGFRAGMENTLOGTABLE_HPP

#include "../Table.hpp"
#include "Fragment.hpp"

class ExpLogFragmentLogTable : public Table
{
public:
	ExpLogFragmentLogTable(Target* target, Fragment* fragment);
	
	~ExpLogFragmentLogTable();

	mpz_class function(int x);

	Fragment * fragment; /**< holds all the useful information about size, accuracy, etc */
	
};






#endif //LOGFRAGMENTLOGTABLE_HPP

