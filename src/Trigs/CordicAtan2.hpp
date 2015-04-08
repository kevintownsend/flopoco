#ifndef CORDICATAN2_HPP
#define CORDICATAN2_HPP

#include "FixAtan2.hpp"
#include "utils.hpp"


#include <vector>

namespace flopoco{ 

	
	class CordicAtan2 : public FixAtan2 {
	  
	  public:

		// Possible TODO: add an option to obtain angles in radians or whatever

		/** Constructor: w is the input and output size, all signed fixed-point number. 
		 Angle is output as a signed number between 00...00 and 11...1111, for 2pi(1-2^-w)
		      pi is 0100..00, etc.
		Actual position of the fixed point in the inputs doesn't matter as long as it is the same for x and y

		*/
		CordicAtan2(Target* target, int wIn, int wOut, int method=0, map<string, double> inputDelays = emptyDelayMap);

		// destructor
		~CordicAtan2();
		


	private:
		int w;                     /**< input and output size (two's complement each, including a sign bit) */
		int	maxIterations;         /**< index at which iterations stop */
		int gXY;                   /**< number of guard bits on the (X,Y) datapath */
		int gA;                    /**< number of guard bits on the Angle datapath */
		bool negateByComplement;   /**< An architecture parameter: we negate negative values to obtain the first octant */
		vector<mpfr_t> atani;      /**< */

		void computeGuardBitsForCORDIC();
		
	};

}

#endif
