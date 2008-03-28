#ifndef INTMULTIPLIER_HPP
#define INTMULTIPLIER_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"


/** The Integer Multiplier class. Receives at input two numbers of 
wInX and wInY widths and outputs a result having the width wInX+wInY */


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


private:
 int partsX;
 int partsY; 
 int number_of_zerosX;
 int number_of_zerosY;
 int multiplier_width;
 bool reverse; 

};



#endif
