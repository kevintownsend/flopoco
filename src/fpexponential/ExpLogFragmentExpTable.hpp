#ifndef EXPLOGFRAGMENTEXPTABLE_HPP
#define EXPLOGFRAGMENTEXPTABLE_HPP

#include "../Table.hpp"
#include "Fragment.hpp"

class ExpLogFragmentExpTable : public Table
{
public:
	ExpLogFragmentExpTable(Target* target, Fragment* fragment);
	
~ExpLogFragmentExpTable();

	mpz_class function(int x);

	Fragment * fragment; /**< holds all the useful information about size, accuracy, etc */
	
};






#endif //LOGFRAGMENTEXPTABLE_HPP

