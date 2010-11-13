/*
  Floating-point pipeline generator for FloPoCo
 
  Author : Florent de Dinechin
 
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

  All rights reserved
 */

#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <cstdlib>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "FPPipeline.hpp"

using namespace std;

namespace flopoco{

  extern vector<Operator*> oplist;



  FPPipeline::FPPipeline(Target* target, string func, int wE_, int wF_): 
    Operator(target), wE(wE_), wF(wF_) {
	

    /* Convert the input string into a sollya evaluation tree */
    expTree = parseString(func.c_str());

    // Name HAS to be unique!
    // will cause weird bugs otherwise
    ostringstream complete_name;
    complete_name << ""; 

    /* If conversion did not succeed (i.e. parse error)
     * throw an exception */
    if (expTree == 0)
      throw "Unable to parse input function.";
		
  }

  FPPipeline::~FPPipeline() {
  }



}
#endif //HAVE_SOLLYA


