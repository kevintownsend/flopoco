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
#include "IntAdder.hpp"

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

	TestIOMap getTestIOMap();
	void fillTestCase(mpz_class a[]);
	
private:
	/** The integer adder object */
	IntAdder *intadd1; 
	/** The integer adder object */
	IntAdder *intadd2; 
	/** The integer adder object */
	IntAdder *intaddClose1; 
	/** The integer adder object */
	IntAdder *intaddClose2; 
	/** The integer adder object */
	IntAdder *intaddClose3; 

	/** The integer adder object */
	IntAdder *intaddFar1; 
	/** The integer adder object */
	IntAdder *intaddFar2; 
	/** The integer adder object */
	IntAdder *intaddFar3; 

	LZOC* leadingZeroCounter;
	Shifter* leftShifter;
	Shifter* rightShifter;
	
	IntAdder* adderList[8];

	int wF;
	int wE;
	int wOutLZC;
	int sizeRightShift;
	
	int swapDifferencePipelineDepth;
	int closePathDepth;
	int farPathDepth;
	int maxPathDepth;
	
};

#endif
