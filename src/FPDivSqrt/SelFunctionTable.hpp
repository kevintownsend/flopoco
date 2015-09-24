#ifndef SELFUNCTIONTABLE_H
#define SELFUNCTIONTABLE_H

#include "Table.hpp"
#include <gmpxx.h>

namespace flopoco
{

	// dmin and dmax are 0.5 and 1 if there is no prescaling; In case of prescaling they can be larger.
	// nbd and nbw are the number of bits of d and w to input to the selection table
	// digit is alpha (the max absolute value of the digit)
	class SelFunctionTable : public Table
	{
		public:
			SelFunctionTable(Target* target, double dmin, double dmax, int nbd, int nbw, int digit, int base, int wIn, int wOut);
			virtual ~SelFunctionTable();
			mpz_class function(int x);
		protected:
			double Dmin;
			double Dmax;
			double ro;
			int nbBitD;
			int nbBitW;
			int digitSet; /**Must not exceed 10!!!**/
			int radix;
		private:
	};
}
#endif // SELFUNCTIONTABLE_H
