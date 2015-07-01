#ifndef SELFUNCTIONTABLE_H
#define SELFUNCTIONTABLE_H

#include "Table.hpp"
#include <gmpxx.h>

namespace flopoco
{
	class SelFunctionTable : public Table
	{
		public:
			SelFunctionTable(Target* target);
			virtual ~SelFunctionTable();
			mpz_class function(int x);
		protected:
		private:
	};
}
#endif // SELFUNCTIONTABLE_H
