#ifndef KARATSUBA_HPP
#define KARATSUBA_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

/** 
 * The Integer Multiplier class. Receives at input two numbers of 
 * wInX and wInY widths and outputs a result having the width wOut=wInX+wInY 
 **/
class Karatsuba : public Operator
{
public:
	Karatsuba(Target* target, int wInX, int wInY);
	~Karatsuba();
	
	/** the width (in bits) of the input X  */
	int wInX; 
	/** the width (in bits) of the input Y  */
	int wInY; 
	/** the width (in bits) of the output R  */
	int wOut; 

	void BuildCombinationalKaratsuba(std::ostream& o, int L, int R, string lName, string rName, string resultName, int depth, string branch);


	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);
	/** A method which generates strings of zeros */ 
	string zero_generator(int n, int margins); 
	
	TestIOMap getTestIOMap();
	void fillTestCase(mpz_class a[]);

private:
	int multiplierXWidth;
	
};
#endif
