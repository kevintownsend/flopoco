#ifndef INTMULTIPLIER_HPP
#define INTMULTIPLIER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntAdder.hpp"

/** The Integer Multiplier class. Receives at input two numbers of 
wInX and wInY widths and outputs a result having the width wInX+wInY */


class IntMultiplier : public Operator
{
public:
	IntMultiplier(Target* target, int wInX, int wInY);
	~IntMultiplier();

	int wInX;
	int wInY;
	int wOut;

	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);
	string zero_generator(int n, int margins);

private:
	int IntAddPipelineDepth;
	int partsX;
	int partsY; 
	int number_of_zerosX;
	int number_of_zerosY;
	int multiplier_width_X;
	int multiplier_width_Y;
	int multiplier_width_avg;
	bool reverse; 
	IntAdder *intadd;
};
#endif
