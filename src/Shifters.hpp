#ifndef SHIFTERS_HPP
#define SHIFTERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
//#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"


/** The Shifter class. Left and right shifters are perfectly
		symmetrical, so both are instances of the Shifter class. Only the
		name of the VHDL instance changes */

typedef enum {Left, Right} ShiftDirection;

class Shifter : public Operator
{
public:
	Shifter(Target* target, int wIn, int maxShift, ShiftDirection dir);
	~Shifter();

	int wIn;
	int maxShift;

	int wOut;
	int wShiftIn;

	// Overloading the virtual functions of Operator
	void output_vhdl(std::ostream& o, std::string name);

	void setup_pipeline();

private:
	ShiftDirection direction;
	/* if boolean true, the corresponding level signal is registered*/ 
	bool level_registered [128];
};



#endif
