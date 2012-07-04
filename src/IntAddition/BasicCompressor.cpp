#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"
#include <gmpxx.h>

#include "BasicCompressor.hpp"

using namespace std;


// personalized parameter
//string BasicCompressor::operatorInfo = "UserDefinedInfo list parameter;



BasicCompressor::BasicCompressor(Target * target, vector<int> height)
	:Operator(target)
{
	ostringstream name;
	
	int c1=height[1];
	int c0=height[0];
	
	int w = c1+c0;
	int wOut = intlog2(c0 + c1*2);
	
        for (int i = 0; i < height.size(); i++)
	{
	  
	}
	name << "Compressor_"<< c0 << c1 << "_" << wOut;
	setName(name.str());
	setCopyrightString("Illyes K. & Popa B. 2012");
        
	addInput("X0", c0);
	addOutput("R", wOut);
	
	if (c1!= 0 )
	{
	addInput("X1", c1);
	
	vhdl << tab << declare("X", c0 + c1) << " <= X1 & X0;\n";
	}
	else
	{
	  vhdl << tab << declare("X", c0 + c1) << " <= X0;\n";
	}
	vhdl << "with X select R <= \n";
	for (mpz_class i = 0; i < (1 << w); i++) {
	  if ( i < (1 << c0))
	  {
		vhdl << tab << "\"" << unsignedBinary(popcnt(i),wOut) << "\" when \""
		     << unsignedBinary(i,w) << "\", \n";
	  }
	  else
	  {
		vhdl << tab << "\"" << unsignedBinary(popcnt(i>>c0)+ popcnt(i),wOut) << "\" when \""
		     << unsignedBinary(i,w) << "\", \n";
	  }
	}
	vhdl << tab << "\"" << std::string (wOut, '-') << "\" when others;\n" << endl;
};


void BasicCompressor::emulate(TestCase * tc)
{
	mpz_class sx = tc->getInputValue("X");
	tc->addExpectedOutput("R", popcnt(sx));
}


