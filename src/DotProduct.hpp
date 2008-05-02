#ifndef DOTPRODUCT_HPP
#define DOTPRODUCT_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FPMultiplier.hpp"
#include "LongAcc.hpp"

/** The DotProduct class.  */
class DotProduct : public Operator
{
public:
	DotProduct(Target* target, int wE, int wFX, int wFY, int MaxMSBX, int LSBA, int MSBA );
	~DotProduct();

	/** The width of the exponent for the inputs X and Y*/
	int wE; 
	/** The width of the fraction for the input X */
	int wFX;
	/** The width of the fraction for the input Y */
	int wFY; 
	/** Maximum expected weight of the MSB of the summand */
	int MaxMSBX; 
	/** The weight of the LSB of the accumulator; determines the final accuracy of the result.*/	
	int LSBA;
	/** The weight of the MSB of the accumulator; has to greater than that of the maximal expected result*/
	int MSBA;

	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);


private:
/** instance of a FPMultiplier */
FPMultiplier* fpMultiplier;
/** instance of a LongAcc */
LongAcc* longAcc;  
};

#endif
