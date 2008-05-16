#ifndef FPADDER_HPP
#define FPADDER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "LZOC.hpp"
#include "Shifters.hpp"
#include "FloFP.hpp"

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


	/** Method which outputs a string of zeros*/
	string zero_generator(int n, int margins);
		
	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);

	/** overloaded method */
	TestCaseList generateStandardTestCases(int n);
	/** overloaded method */
	TestCaseList generateRandomTestCases(int n);
	
private:
	LZOC* leadingZeroCounter;
	Shifter* leftShifter;
	Shifter* rightShifter;

	int wF;
	int wE;
	int wOutLZC;
	int sizeRightShift;
};

#endif
