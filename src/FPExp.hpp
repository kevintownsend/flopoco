#ifndef __FPEXP_HPP
#define __FPEXP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "Table.hpp"
#include "DualTable.hpp"

class Fragment;
class FPNumber;


namespace flopoco{


	class FPExp : public Operator
	{
	public:

		class magicExpTable: public DualTable {
		public:
			static const int sizeH = 27;
			static const int sizeL = 9;

			magicExpTable(Target* target) :

				DualTable(target, 9, 36, 0, 511) {
				ostringstream name; 
				srcFileName="FPExp::MagicSPExpTable";
				name <<"MagicSPExpTable";
				setName(name.str());
			};

			mpz_class function(int x){
				mpz_class h, l;
				mpfr_t a, y, yh,yl, one;

				// convert x to 2's compliment
				int xs=x;
				if(xs>>(wIn-1))
					xs-=(1<<wIn);

				mpfr_init2(a, wIn);

				mpfr_set_si(a, xs, GMP_RNDN);
				mpfr_div_2si(a, a, wIn, GMP_RNDN); // now a in [-1/2, 1/2[
				mpfr_init2(yh, 27); 
				mpfr_exp(yh, a, GMP_RNDN); // e^x in [0.6,1.7[
				mpfr_init2(one, 16);
				mpfr_set_d(one, 1.0, GMP_RNDN);

				mpfr_mul_2si(yh, yh, 26, GMP_RNDN);
				mpfr_get_z(h.get_mpz_t(), yh,  GMP_RNDN);  // Should be a 27-bit number

				mpfr_set_ui(a, x, GMP_RNDN);
				mpfr_div_2si(a, a, 2*wIn, GMP_RNDN); // now a in [0,1[. 2^-9

				mpfr_init2(y, 600); 
				mpfr_exp(y, a, GMP_RNDN); // e^(2^-9 x)
				mpfr_sub(y, y, a, GMP_RNDN); // e^(2^-9 x) -x 
				mpfr_sub(y, y, one, GMP_RNDN); // e^(2^-9 x) -x -1

				//now scale so that first bit 
				mpfr_set_d(one, 1.0, GMP_RNDN);
				mpfr_mul_2si(one, one, 26, GMP_RNDN); //  2^-26
				mpfr_init2(yl, 9); 
				mpfr_mul(yl, y, one, GMP_RNDN); // y scaled up to  [0..511] (actually [0..258])
				mpfr_get_z(l.get_mpz_t(), yl,  GMP_RNDN);

				// debug
				if((h>=(1<<27)) || l>=512 || h<0 || l<0)
					REPORT(0, "Ouch!!!!!" <<"x=" << x << " " << xs << "    " << h << " " << l );

				//cout << x << "\t" << h << "\t" << l <<endl;
				mpfr_clears(y, yh, yl, a, one, NULL);
						
				return l + (h<<9);
			};

		};


		class firstExpTable: public Table {
		public:

			firstExpTable(Target* target, int wIn, int wOut) : 
				Table(target, wIn, wOut) {
				ostringstream name; 
				srcFileName="FPExp::firstExpTable";
				name <<"firstExpTable_" << wIn << "_" << wOut;
				setName(name.str());
			};

			mpz_class function(int x){
				mpz_class h;
				mpfr_t a, y;

				// convert x to 2's compliment
				int xs=x;
				if(xs>>(wIn-1))
					xs-=(1<<wIn);

				mpfr_init2(a, wIn);
				mpfr_set_si(a, xs, GMP_RNDN);
				mpfr_div_2si(a, a, wIn, GMP_RNDN); // now a in [-1/2, 1/2[
				mpfr_init2(y, wOut); 
				mpfr_exp(y, a, GMP_RNDN); // in [0.6, 1.7], MSB is 1

				mpfr_mul_2si(y, y, wOut-1, GMP_RNDN);
				mpfr_get_z(h.get_mpz_t(), y,  GMP_RNDN);

				// debug
				if((h>=(mpz_class(1)<<wOut)) || h<0)
					REPORT(0, "Ouch!!!!!" << h);

				//cout << x << "\t" << h << "\t" << l <<endl;
				mpfr_clears(y, a, NULL);
						
				return h;
			};
		};

		/** Table used for small precisions instead of the polynomial approximator */
		class LowerExpTable: public Table {
		public:
			int k;
			int lsb;
			
			LowerExpTable(Target* target, int k_, int wIn_, int lsb_) :
				// input size is k, output size lsb-(2*k)+1
				Table(target, wIn_, lsb_-(2*k_)+1),   k(k_),   lsb(lsb_) {
				ostringstream name; 
				srcFileName="FPExp::LowerExpTable";
				name <<"LowerExpTable_" << k << "_" << lsb;
				setName(name.str());
			};

			mpz_class function(int x){
				mpz_class l;
				mpfr_t a, y,yl, one;
				mpfr_init2(one, 16);
				mpfr_set_d(one, 1.0, GMP_RNDN);

				mpfr_init2(a, wIn);
				mpfr_set_ui(a, x, GMP_RNDN);
				mpfr_div_2si(a, a, k+wIn, GMP_RNDN); // now a in [0,1[. 2^-k

				mpfr_init2(y, 600);       
				mpfr_exp(y, a, GMP_RNDN); // e^a
				mpfr_sub(y, y, a, GMP_RNDN); // e^a - a 

				mpfr_sub(y, y, one, GMP_RNDN); // e^a -a -1

				//now scale to an int. 
				// a in  [0,1[. 2^-k  so   msb(e^a-a-1 = a^2/2 +...) is 2^-2k 
				// to send it to weight wOut-1 we need to multiply by 2^(2k + wOut-1)  
				mpfr_set_d(one, 1.0, GMP_RNDN);
				mpfr_mul_2si(one, one, 2*k+wOut-1, GMP_RNDN); // 
				mpfr_init2(yl, wOut); 
				mpfr_mul(yl, y, one, GMP_RNDN); // y scaled up)
				mpfr_get_z(l.get_mpz_t(), yl,  GMP_RNDN);

				// // debug -- shouldnt this be in Table?
				// if( l>=(1<<wOut) || l<0)
				// 	REPORT(0, "Ouch!!!!!" <<"x=" << x << " " << xs << "    " << h << " " << l );

				//cout << x << "\t" << h << "\t" << l <<endl;
				mpfr_clears(y, yl, a, one, NULL);
						
				return l;

			};
		};

		/** The constructor with manual control of all options
		    * @param wE exponent size
		    * @param wF fraction size
		    * @param k size of the input to the first table 
		    * @param d  degree of the polynomial approximation (if k=d=0, the constructor tries to compute sensible values)
		    * @param guardBits number of gard bits, defaults to 3, which is enough for faithful rounding 
		    * @param fullInput boolean, if true input mantissa is of size wE+wF+1, so that input shift doesn't padd it with 0s (useful for FPPow)
		    */

		FPExp(Target* target, int wE, int wF, int k, int d, int guardBits=3, bool fullInput=false);
		~FPExp();
		
		// Overloading the virtual functions of Operator
		// void outputVHDL(std::ostream& o, std::string name);
		
		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);
		TestCase* buildRandomTestCase(int i);

	private:
		int wE; /**< Exponent size */
		int wF; /**< Fraction size */
		int k;  /**< Size of the address bits for the first table  */
		int d;  /**< Degree of the polynomial approximation */
		int g;  /**< Number of guard bits */
	};

}
#endif
