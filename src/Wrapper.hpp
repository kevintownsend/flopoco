#ifndef WRAPPER_HPP
#define WRAPPER_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Wrapper.hpp"

/** 
 * A wrapper is a VHDL entity that places registers before and after
 * an operator, so that you can synthesize it and get delay and area,
 * without the synthesis tools optimizing out your design because it
 * is connected to nothing.  
 **/
class Wrapper : public Operator
{
public:
	// The operator to wrap
	Operator* op;
	Wrapper(Target* target, Operator* op);
	~Wrapper();
	void output_vhdl(ostream& o, string name);
};


#endif
