#ifndef __FPEXP_HPP
#define __FPEXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "DualTable.hpp"

class Fragment;
class FPNumber;


namespace flopoco{


	class FPExp : public Operator
	{
	public:

		class magicExpTable: public DualTable {
		public:
			static const int sizeH = 28;
			static const int sizeL = 8;

			magicExpTable(Target* target) : 
				DualTable(target, 9, 36, 0, 511) {
				ostringstream name; 
				name <<"MagicSPExpTable";
				setName(name.str());
			};

			mpz_class function(int x){
				mpz_class h, l;
				mpfr_t a, y, yh,yl, one;
				mpfr_init2(a, wIn);
				mpfr_set_ui(a, x, GMP_RNDN);
				mpfr_div_2si(a, a, wIn, GMP_RNDN); // now a in [O1[
				mpfr_init2(yh, 27); 
				mpfr_exp(yh, a, GMP_RNDN);
				mpfr_init2(one, 16);
				mpfr_set_d(one, 1.0, GMP_RNDN);
				mpfr_sub(yh, yh, one, GMP_RNDN); // e^x -1 in [0,2[

				mpfr_mul_2si(yh, yh, 26, GMP_RNDN);
				mpfr_get_z(h.get_mpz_t(), yh,  GMP_RNDN);

				mpfr_set_ui(a, x, GMP_RNDN);
				mpfr_div_2si(a, a, 2*wIn, GMP_RNDN); // now a in [O1[. 2^-9

				mpfr_init2(y, 600); 
				mpfr_exp(y, a, GMP_RNDN); // e^(2^-9 x)
				mpfr_sub(y, y, a, GMP_RNDN); // e^(2^-9 x) -x 
				mpfr_sub(y, y, one, GMP_RNDN); // e^(2^-9 x) -x -1
				// Now add the round bit which has value 2^{-26}. This will turn the truncation of e^z-1 to 17 bits into a rounding
				mpfr_div_2si(one, one, 26, GMP_RNDN); //  2^-26
				mpfr_add(y, y, one, GMP_RNDN); // e^(2^-9 x) -x -1 + roundbit

				//now scale so that first bit 
				mpfr_set_d(one, 1.0, GMP_RNDN);
				mpfr_mul_2si(one, one, 18+9, GMP_RNDN); //  2^-26
 				mpfr_init2(yl, 9); 
				mpfr_mul(yl, y, one, GMP_RNDN); // y scaled up to  [0..511] (actually [0..258])
				mpfr_get_z(l.get_mpz_t(), yl,  GMP_RNDN);

				// debug
				if((h>=(1<<27)) || l>=512 || h<0 || l<0)
					cout << "!!!!!" << h << l <<endl;

				//cout << x << "\t" << h << "\t" << l <<endl;
				mpfr_clears(y, yh, yl, a, one, NULL);
						
				return l + (h<<9);
			};
		};

		FPExp(Target* target, int wE, int wF);
		~FPExp();
		
		// Overloading the virtual functions of Operator
		// void outputVHDL(std::ostream& o, std::string name);
		
		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);
		
	private:
		int wE, wF;
		int k; // Size of the address bits for the first table
		int result_length, g;
		// Fragment *f;
		double area, max_error;
	};

}
#endif
