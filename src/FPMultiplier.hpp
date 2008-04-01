#ifndef FPMULTIPLIERS_HPP
#define FPMULTIPLIERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntMultiplier.hpp"

/** The FPMultiplier class. Left and right multipliers are perfectly
    symmetrical, so both are instances of the FPMultiplier class. Only the
    name of the VHDL instance changes */


class FPMultiplier : public Operator
{
public:
  FPMultiplier(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, int norm);
  ~FPMultiplier();


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

   
   string zero_generator(int n, int margins);
private:
IntMultiplier* intmult;
  
};



#endif
