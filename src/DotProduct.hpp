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

	/**
	 * The DotProduct constructor
	 * @param[in]		target	the target device
	 * @param[in]		wE			the width of the exponent for the inputs X and Y
	 * @param[in]		wFX     the width of the fraction for the input X
	 * @param[in]		wFY     the width of the fraction for the input Y
	 * @param[in]		MaxMSBX	maximum expected weight of the MSB of the summand
	 * @param[in]		LSBA    The weight of the LSB of the accumulator; determines the final accuracy of the result
	 * @param[in]		MSBA    The weight of the MSB of the accumulator; has to greater than that of the maximal expected result
	**/ 
	DotProduct(Target* target, int wE, int wFX, int wFY, int MaxMSBX, int LSBA, int MSBA );

	/**
	 * DotProduct destructor
	 */
	~DotProduct();

	/**
	 * Method belonging to the Operator class overloaded by the DotProduct class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);

	/** 
	 * Sets the default name of this operator
	 */
	void setOperatorName(); 
	
protected:
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

private:
/** instance of a FPMultiplier */
FPMultiplier* fpMultiplier;
/** instance of a LongAcc */
LongAcc* longAcc;  
};

#endif
