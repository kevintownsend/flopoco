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
	IntMultiplier(Target* target, int wInX, int wInY);
	~IntMultiplier();
	
	/** the width (in bits) of the input X  */
	int wInX; 
	/** the width (in bits) of the input Y  */
	int wInY; 
	/** the width (in bits) of the output R  */
	int wOut; 


	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);
	/** A method which generates strings of zeros */ 
	string zero_generator(int n, int margins); 
	
	/** Overloaded method */
	virtual TestCaseList generateStandardTestCases(int n);
	/** Overloaded method */ 
	virtual TestCaseList generateRandomTestCases(int n); 

	TestIOMap getTestIOMap();
	void fillTestCase(mpz_class a[]);

private:
    //Specific IntMultiplier methods
	void pad_inputs(std::ostream& o);
	void split(std::ostream& o, std::string name, int parts, int part_width);			  
	void do_multiplications(std::ostream& o);
	void decompose_low_high(std::ostream& o);
	void regroup_low_high(std::ostream& o);
	void map_adder(std::ostream& o,std::string left_term, std::string right_term, std::string result );
	void connect_low(std::ostream& o);
	void high_to_buffer(std::ostream& o);
	void connect_buffers(std::ostream& o);
	void link_low_adder_structure(std::ostream& o);
	void link_high_adder_structure(std::ostream& o);
	void connect_partial_bits(std::ostream& o);
	void pipeline_addition(std::ostream& o);
	void delay_partial_bits(std::ostream& o);

	// Specific attributes
	/** The integer adder object */
	IntAdder *intadd; 
	/** The pipeline depth of the adder */
	int IntAddPipelineDepth; 
	/** The number of parts that the input X will be split in */
	int partsX; 
	/** The number of parts that the input Y will be split in */
	int partsY; 
	/** The number of zeros that the input X will be padded with so that it's length reaches a multiple of the suggested multiplier width */
	int number_of_zerosX; 
	/** The number of zeros that the input Y will be padded with so that it's length reaches a multiple of the suggested multiplier width */
	int number_of_zerosY; 
	/** The X width of the multiplier */
	int multiplier_width_X; 
	/** The Y width of the multiplier */
	int multiplier_width_Y; 
	/** Signals if we are doing the multiplication X*Y of Y*X */
	bool reverse; 
	/** The recommended addition chunk width so that the design reaches the desired frequency */
	int addition_chunk_width; 
	/** The number of pipeline levels that the last addition of the algorithm is split in so that we reach the desired frequency */
	int pipe_levels; 
	
};
#endif
