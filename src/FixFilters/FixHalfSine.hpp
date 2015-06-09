#ifndef FIXHALFSINE_HPP
#define FIXHALFSINE_HPP

#include "Operator.hpp"
#include "FixFIR.hpp"
#include "utils.hpp"

#include <iostream>
#include <sstream>
#include <vector>

#include <string.h>


namespace flopoco{

	class FixHalfSine : public FixFIR
	{
	public:

		FixHalfSine(Target* target, int lsb_, int N_);

		virtual ~FixHalfSine();

	private: 
		int N; /* FixFIR::n = 2*N */
	};

}


#endif
