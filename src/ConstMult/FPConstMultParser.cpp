/*
 * A correctly-rounded multiplier by an arbitrary real constant for FloPoCo

  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License, 2008-2010.
*/

#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../sollya.h" // TODO : fix upstream Sollya, or fix in FloPoCo
#include "../utils.hpp"
#include "../Operator.hpp"
#include "FPConstMult.hpp"
#include "FPConstMultParser.hpp"
#include "../FPNumber.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;





	FPConstMultParser::~FPConstMultParser() {
		// TODO but who cares really
	}


}
#endif //HAVE_SOLLYA
