#ifndef FirstKCMTable_HPP
#define FirstKCMTable_HPP

#include "../Table.hpp"

namespace flopoco{

	class FirstKCMTable : public Table
	{
	public:
		FirstKCMTable(Target* target, int wIn, int wOut, mpz_class C);
	
		~FirstKCMTable();


		mpz_class function(int x);

		mpz_class C_; //the constant

		int wOut_; /**< the width of the output */
	
	};

}
#endif //FirstKCMTable_HPP

