#ifndef SELFUNCTIONTABLE_H
#define SELFUNCTIONTABLE_H

#include "Table.hpp"
#include <gmpxx.h>

namespace flopoco
{
	class SelFunctionTable : public Table
	{
		public:
			SelFunctionTable(Target* target, float dmin, float dmax, int nbd, int nbw, int digit, int base, int wIn, int wOut);
			virtual ~SelFunctionTable();
			mpz_class function(int x);
		protected:
			float Dmin;
			float Dmax;
			float ro;
			int nbBitD;
			int nbBitW;
			int digitSet; /**Must not exceed 10!!!**/
			int radix;
		private:
	};
}
#endif // SELFUNCTIONTABLE_H
