#ifndef FPADDER_HPP
#define FPADDER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"


/** The FPAdder class */
class FPAdder : public Operator
{
public:
	FPAdder(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR);
	~FPAdder();

	/** The width of the exponent for the input X */
	int wEX; 
	/** The width of the fraction for the input X */
	int wFX; 
	/** The width of the exponent for the input Y */
	int wEY; 
	/** The width of the fraction for the input Y */
	int wFY; 
	/** The width of the exponent for the output R */
	int wER; 
	/** The width of the fraction for the output R */
	int wFR;
	/** Signal if the output of the operator is to be or not normalized*/


	
	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);


	
private:
	int wF;
	int wE;
	int nBlock;
	int wN;
  int wLastBlock;
  int wRanges;
	int wTree;  
  int rangeLen;
	int rangeIdx;
	int wAddr;
};

#endif
