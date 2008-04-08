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

	int wInX;
	int wInY;
	int wOut;

	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);
	//a zero generator methodgi
	string zero_generator(int n, int margins);

private:
	//====================== Specific IntMultiplier methods
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

	//====================== Specific attributes
	int IntAddPipelineDepth;
	int partsX;
	int partsY; 
	int number_of_zerosX;
	int number_of_zerosY;
	int multiplier_width_X;
	int multiplier_width_Y;
	int multiplier_width_avg;
	bool reverse; 
	int addition_chunk_width;
	int pipe_levels;
	IntAdder *intadd;
};
#endif
