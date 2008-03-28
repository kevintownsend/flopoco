#ifndef INTADDERS_HPP
#define INTADDERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
//#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"


/** The IntAdder class for experimenting with adders.

Not useful to users, but useful to evaluate numerical values of
_fastcarry_delay for a new target  */

class IntAdder : public Operator
{
public:
  IntAdder(Target* target, int wIn);
  ~IntAdder();

  int wIn;


  // Overloading the virtual functions of Operator
  void output_vhdl(std::ostream& o, std::string name);


private:
  int chunk_size;
  int last_chunk_size; // the last one is slightly smaller
  int pipe_levels;
};



#endif
