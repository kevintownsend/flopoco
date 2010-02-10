#ifndef TableGenerator_HPP
#define TableGenerator_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "Operator.hpp"
#include "HOTBM/sollya.h"	// Do NOT use libsollya from user's environment

namespace flopoco{

	/** The TableGenerator class.  */
	class TableGenerator : public Operator{

		public:

			 /* TODO: Doxygen parameters*/ 
			TableGenerator(Target* target, int wIn, int wOut );

			/**
			 * TableGenerator destructor
			 */
			~TableGenerator();

		protected:
			int wIn_;   /**< TODO: Description*/ 
			int wOut_;  /**< TODO: Description*/
	};
}
#endif
