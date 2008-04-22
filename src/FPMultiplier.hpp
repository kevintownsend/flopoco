#ifndef FPMULTIPLIERS_HPP
#define FPMULTIPLIERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntMultiplier.hpp"

/** The FPMultiplier class. Left and right multipliers are perfectly
    symmetrical, so both are instances of the FPMultiplier class. Only the
    name of the VHDL instance changes */


class FPMultiplier : public Operator
{
public:
	FPMultiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm);
	~FPMultiplier();

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
	bool normalized;
	/** The pipeline depth of the integer multiplier used to multiply significands*/
	int IntMultPipelineDepth;

	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);

	/** Method which outputs a string of zeros*/
	string zero_generator(int n, int margins);

	/** overloaded method */
	TestCaseList generateStandardTestCases(int n);
	/** overloaded method */
	TestCaseList generateRandomTestCases(int n);

private:
	/** The integer multiplier object */
	IntMultiplier* intmult;
	/** The width of a signal which reunites a leading 0, and a part of the significands' product */
	int reunion_signal_width;
	/** The width of the chunks for the addition used for rounding */
	int addition_chunk_width;
	/** The number of parts that the reunion_signal will be split in in order to breake the addition in small pieces and pipeline it */
	int reunion_signal_parts;
	/** The width of the last chunk of the reunion signal. This width is usually smaller than the addition chunk width*/
	int addition_last_chunk_width;
  
};

#endif
