#ifndef LONGACC2FP_HPP
#define LONGACC2FP_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Operator.hpp"
#include "Shifters.hpp"
#include "LZOC.hpp"

class LongAcc2FP : public Operator
{
public:
	LongAcc2FP(Target* target, int MaxMSBX, int LSBA, int MSBA, int wE_out, int wF_out);
	~LongAcc2FP();

	int MaxMSBX;
	int LSBA; 
	int MSBA; 
	int wE_out; 
	int wF_out; 
	

	// Overloading the virtual functions of Operator
	void output_vhdl(ostream& o, string name);

private:
	LZOC* leadZOCounter;
	Shifter* leftShifter;
	
	int sizeAcc;         //  = MSBA-LSBA+1;
	int expBias;//the exponent bias value
	int wOutLZOC; //the number of bits that the leading zero/one conunter outputs the result on

};


#endif
