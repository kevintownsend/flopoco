#include <iostream>
#include <sstream>
#include "FixConstant.hpp"
#include "utils.hpp"

using namespace std;

namespace flopoco{


	FixConstant::FixConstant(const int MSB_, const int LSB_, const bool signed_, const mpfr_t val_) : 
		MSB(MSB_), LSB(LSB_), width(MSB_-LSB_+1), isSignedFormat(signed_) {
		mpfr_init2(fpValue, width);
		mpfr_set(fpValue, val_, GMP_RNDN); // TODO check no error?
	}




#if 0 // Warning the following code is unfinished and untested
	FixConstant::FixConstant(const bool signed_, const mpfr_t val_) : 
		isSignedFormat(signed_) {

		bool sign;
		mpz_class zval;
		LSB =  mpfr_get_z_exp (zval.get_mpz_t(),val_);

		// conversion to (sign, absval)
		if(isSignedFormat) {
			if(zval<0) {
				sign=true; 
				zval=-zval;
			}
			else 
				sign=false;
		}
		else {// unsigned
			if(zval<0) {
				throw("Signed input to FixConstant(unsigned)");
			}
		}

		// Normalisation?
		cout << "Normalisation";
		while(zval & mpz_class(1) == 0) {
			cout << ".";
			zval = zval>>1;
			LSB++;
		}
		MSB=0;
		while(zval!=0) {
			cout << "*";
			zval = zval>>1;
			MSB++;			
		}

		if(isSignedFormat)
			MSB++;
		/* -1 s'écrit */

		mpfr_init2(fpValue, width);
		mpfr_set(fpValue, val_, GMP_RNDN); // TODO check no error?
	}
#endif



	FixConstant::~FixConstant(){
		mpfr_clear(fpValue);
	} 
	

	bool FixConstant::isSigned() {
		return isSignedFormat;
	}

	int FixConstant::getMSB() {
		return MSB;
	}
	int FixConstant::getLSB() {
		return LSB;
	}


	std::string FixConstant::getBitVector(int margins){
		// use the function in utils.hpp
		if(isSignedFormat)
			return signedFixPointNumber(fpValue, MSB, LSB);
		else
			return unsignedFixPointNumber(fpValue, MSB, LSB);
	}
}
