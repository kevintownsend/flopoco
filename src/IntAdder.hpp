#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>

#include "Operator.hpp"


/** The IntAdder class for experimenting with adders.

Not useful to users, but useful to evaluate numerical values of
_fastcarry_delay for a new target  */
class IntAdder : public Operator
{
public:
	IntAdder(Target* target, int wIn, const int p=0);
	~IntAdder();
	/** the width for X, Y and R*/
	int wIn;


	/** Overloading the virtual functions of Operator */
	void output_vhdl(std::ostream& o, std::string name);

	TestIOMap getTestIOMap();
	void fillTestCase(mpz_class a[]);

private:
	/** the maximum chunk size for an addition so that the requested frequency can theoretically be used*/
	int chunk_size;
	/** the last chunk size - the last one is slightly smaller*/
	int last_chunk_size; 
	/** the number of pieline levels of the sequential operator*/ 
	int pipe_levels;
	
	//force pipeline
	int forcePipeline;
};

#endif
