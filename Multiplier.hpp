#ifndef MULTIPLIERS_HPP
#define MULTIPLIERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntMultiplier.hpp"

/** The Multiplier class. Left and right multipliers are perfectly
    symmetrical, so both are instances of the Multiplier class. Only the
    name of the VHDL instance changes */


class Multiplier : public Operator
{
public:
  Multiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm);
  ~Multiplier();


	int wEX; 
	int wFX; 
	int wEY; 
	int wFY; 
	int wER; 
	int wFR;
	//bool pipelined;
	bool normalized;
	int IntMultPipelineDepth;
  // Overloading the virtual functions of Operator
  void output_vhdl(std::ostream& o, std::string name);

  //void setup_pipeline();

private:
IntMultiplier* intmult;
  /* if boolean true, the corresponding level signal is registered*/ 
  //bool level_registered [128];
};



#endif
