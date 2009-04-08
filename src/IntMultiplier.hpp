#ifndef INTMULTIPLIER_HPP
#define INTMULTIPLIER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntAdder.hpp"

/** 
 * The Integer Multiplier class. Receives at input two numbers of 
 * wInX and wInY widths and outputs a result having the width wOut=wInX+wInY 
 **/
class IntMultiplier : public Operator
{
public:
	/** 
	 * The constructor of the IntMultiplier class
	 * @param target argument of type Target containing the data for which this operator will be optimized
	 * @param wInX integer argument representing the width in bits of the input X 
	 * @param wInY integer argument representing the width in bits of the input Y
	 **/
	IntMultiplier(Target* target, int wInX, int wInY);
	
	/** IntMultiplier destructor */
	~IntMultiplier();
	
	/**
	 * Method belonging to the Operator class overloaded by the IntMultiplier class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);


	void emulate(TestCase* tc);

protected:

	int wInX_; /**< the width (in bits) of the input X  */
	int wInY_; /**< the width (in bits) of the input Y  */
	int wOut_; /**< the width (in bits) of the output R  */

private:

	IntAdder *intAdd_;            /**< The integer adder object */
	int      IntAddPipelineDepth; /**< The pipeline depth of the adder */
	int      partsX_; 	          /**< The number of parts that the input X will be split in */
	int      partsY_; 	          /**< The number of parts that the input Y will be split in */
	int      numberOfZerosX_; 	  /**< The number of zeros that the input X will be padded with so that it's length reaches a multiple of the suggested multiplier width */
	int      numberOfZerosY_; 	  /**< The number of zeros that the input Y will be padded with so that it's length reaches a multiple of the suggested multiplier width */
	int      multiplierWidthX_;   /**< The X width of the multiplier */
	int      multiplierWidthY_;   /**< The Y width of the multiplier */
	bool     reverse_; 	          /**< Signals if we are doing the multiplication X*Y of Y*X */
	IntAdder *intAdd1_; 
	IntAdder *intAdd2_; 

    //Specific IntMultiplier methods
    /**
	 * A method which pads the inputs with 0. 
	 * @param[in,out] o - the stream to which the new line will be added
	 **/
	void padInputs(std::ostream& o);

	/**
	 * A method which splits the input given by name into [parts] parts, each part having the width part_width
	 * @param[in,out]	o 			- the stream to which the code will be inputed
	 * @param[in]		name		- the name of the input string
	 * @param[in]		parts		- the number of parts the input will be split in
	 * @param[in]		part_width	- the width of the part	
	 **/
	void split(std::ostream& o, std::string name, int parts, int part_width);			  

	/**
	 * A does the multiplications of the parts X and Y were splitted in
	 * @param[in,out] o - the stream to which the multiplications are added
	 **/
	void doMultiplications(std::ostream& o);

	/**
	 * Decompose the multiplication of two chunks into high and low parts
	 * @param[in,out] o - the stream to which the decomposition is written
	 **/
	void decomposeLowHigh(std::ostream& o);

	/**
	 * Regroup the decomposition result into vectors low and high vectors
	 * @param[in,out] o - the stream to which the regrouping is written
	 **/
	void regroupLowHigh(std::ostream& o);

	/**
	 * Adder mapping
	 * @param[in,out] 	o 			- the stream to which the mapping is written
	 * @param[in]		left_term	- the name of the left term of the addition as a string
	 * @param[in]		right_term	- the name of the right term of the addition as a string
	 * @param[in]		result		- the name of the result of the addition
	 **/
	void mapAdder(std::ostream& o,std::string left_term, std::string right_term, std::string result );

	/**
	 * Connect the low part of the multiplication to the low part of the adder structure
	 * @param[in,out] 	o 			- the stream
	 **/
	void connectLow(std::ostream& o);

	/**
	 * Connect the high part of the multiplication to a buffer
	 * @param[in,out] 	o 			- the stream
	 **/
	void highToBuffer(std::ostream& o);

	/**
	 * Connect the buffer to the high part of the adder structure
	 * @param[in,out] 	o 			- the stream
	 **/
	void connectBuffers(std::ostream& o);

	/**
	 * Link together the components of the low part of the adder structure
	 * @param[in,out] 	o 			- the stream
	 **/
	void linkLowAdderStructure(std::ostream& o);

	/**
	 * Link together the components of the high part of the adder structure
	 * @param[in,out] 	o 			- the stream
	 **/
	void linkHighAdderStructure(std::ostream& o);

	/**
	 * The partial bits structure outputs the low part of the result by using short additions
	 * @param[in,out] 	o 			- the stream
	 **/
	void connectPartialBits(std::ostream& o);

	/**
	 * The addition which gives the high part of the product between X and Y is pipelined to obtain better performance
	 * @param[in,out] 	o 			- the stream
	 **/
	void pipelineAddition(std::ostream& o);

	/**
	 * Delay the low bits of the result with as many clock counts as the high bits addition pipeline take
	 * @param[in,out] 	o 			- the stream
	 **/
	void delayPartialBits(std::ostream& o);
};
#endif
