#ifndef FPMULTIPLIERS_HPP
#define FPMULTIPLIERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntMultiplier.hpp"

/** The FPMultiplier class */
class FPMultiplier : public Operator
{
public:
	
	/**
	 * The FPMutliplier constructor
	 * @param[in]		target		the target device
	 * @param[in]		wEX			the the with of the exponent for the f-p number X
	 * @param[in]		wFX			the the with of the fraction for the f-p number X
	 * @param[in]		wEY			the the with of the exponent for the f-p number Y
	 * @param[in]		wFY			the the with of the fraction for the f-p number Y
	 * @param[in]		wER			the the with of the exponent for the multiplication result
	 * @param[in]		wFR			the the with of the fraction for the multiplication result
	 **/
	FPMultiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm);

	/**
	 * FPMultiplier destructor
	 */
	~FPMultiplier();

	/**
	 * Method belonging to the Operator class overloaded by the FPMultiplier class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);

	/**
	 * Gets the signals which are interesting for TestCases.
	 * @see TestIOMap
	 */
	TestIOMap getTestIOMap();

	/**
	 * Gets the correct value associated to one or more inputs.
	 * @param a the array which contains both already filled inputs and
	 *          to be filled outputs in the order specified in getTestIOMap.
	 */
	void fillTestCase(mpz_class a[]);

	/** 
	 * Sets the default name of this operator
	 */
	void setOperatorName(); 

protected:
	
	int  wEX_;                  /**< The width of the exponent for the input X */
	int  wFX_;                  /**< The width of the fraction for the input X */
	int  wEY_;                  /**< The width of the exponent for the input Y */
	int  wFY_;                  /**< The width of the fraction for the input Y */
	int  wER_;                  /**< The width of the exponent for the output R */
	int  wFR_;                  /**< The width of the fraction for the output R */
	bool normalized_;	       /**< Signal if the output of the operator is to be or not normalized*/
	int  intMultPipelineDepth_; /**< The pipeline depth of the integer multiplier used to multiply significands*/


private:

	IntMultiplier* intmult_;     /**< The integer multiplier object */
	int reunionSignalWidth_;     /**< The width of a signal which reunites a leading 0, and a part of the significands' product */
	int additionChunkWidth_;     /**< The width of the chunks for the addition used for rounding */
	int reunionSignalParts_; 	 /**< The number of parts that the reunion_signal will be split in in order to breake the addition in small pieces and pipeline it */
	int additionLastChunkWidth_; /**< The width of the last chunk of the reunion signal. This width is usually smaller than the addition chunk width*/
 
};
#endif
