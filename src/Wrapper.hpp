#ifndef WRAPPER_HPP
#define WRAPPER_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Wrapper.hpp"

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
