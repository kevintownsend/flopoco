#ifndef MUX_HPP
#define MUX_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Operator.hpp"

/** Constructs a multiplexer  */
class Mux : public Operator
{
public:
	Mux(Target* target, int wIn, int n);
	~Mux();

	/** The width the multiplexer inputs*/
	int wIn; 
	/** the number of multiplexer inputs */
	int n;
	/** the number of bits of the selector*/
	int wAddr; 

	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);

private:

};
#endif
