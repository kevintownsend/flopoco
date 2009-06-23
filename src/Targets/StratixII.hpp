/*
 * A model of Stratix II FPGA 
 *
 * Author : Florent de Dinechin
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#ifndef STRATIXII_HPP
#define  STRATIXII_HPP
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../Target.hpp"

/** Class for representing an StratixII target */
class StratixII : public Target
{
public:

	/** The default constructor. */  
	StratixII() : Target()	{
		fastcarryDelay_ = 3.5e-11; 	// aproximately right    
		elemWireDelay_  = 0.3e-11; 	// ???
		lut2lutDelay_   = 1.5e-10; 	// ???
		lutDelay_       = 0.378e-9; // 378 ps  
		ffDelay_        = 0.127e-9; // 127 ps LE register clock-to-output max delay for -3 speed grade
		multXInputs_    = 36;
		multYInputs_    = 36;
		// all these values are set precisely to match the Stratix 2
		lut2_ 				= 0.162e-9; // obtained from Handbook
		lut3_				= 0.280e-9; // obtained from Handbook
		lut4_				= 0.378e-9; // obtained from Handbook
		innerLABcarryDelay_	= 0.146e-9; // obtained from Quartus 2 Chip Planner
		interLABcarryDelay_	= 0.245e-9; // obtained from Quartus 2 Chip Planner
		shareOutToCarryOut_	= 0.172e-9; // obtained from Quartus 2 Chip Planner
		muxStoO_			= 0.189e-9; // obtained from Quartus 2 Chip Planner by subtraction
		fdCtoQ_				= 0.352e-9; // obtained from Quartus 2 Chip Planner by subtraction
		carryInToSumOut_	= 0.125e-9;	// obtained from Quartus 2 Chip Planner
		
	}
	
	/** The destructor */
	virtual ~StratixII() {}

	/** overloading the virtual functions of Target
	 * @see the target class for more details 
	 */
	double carryPropagateDelay();
	double adderDelay(int size);
	double localWireDelay();
	double lutDelay();
	double ffDelay();
	double distantWireDelay(int n);
	bool   suggestSubmultSize(int &x, int &y, int wInX, int wInY);
	bool   suggestSubaddSize(int &x, int wIn);
	bool   suggestSlackSubaddSize(int &x, int wIn, double slack);

private:

	double fastcarryDelay_; /**< The delay of the fast carry chain */
	double lut2lutDelay_;   /**< The delay between two LUTs */
	double ffDelay_;   		/**< The delay between two flipflops (not including elemWireDelay_) */
	double elemWireDelay_;  /**< The elementary wire dealy (for computing the distant wire delay) */
	double lutDelay_;       /**< The LUT delay (in seconds)*/
	
	// Added by Sebi
	double lut2_;           	/**< The LUT delay for 2 inputs */
	double lut3_;           	/**< The LUT delay for 3 inputs */
	double lut4_;           	/**< The LUT delay for 4 inputs */
	double innerLABcarryDelay_;	/**< The wire delay between the upper and lower parts of a LAB --> R4 & C4 interconnects */	
	double interLABcarryDelay_;	/**< The approximate wire between two LABs --> R24 & C16 interconnects */	
	double shareOutToCarryOut_;	/**< The delay between the shared arithmetic out of one LAB and the carry out of the following LAB */	
	double muxStoO_;			/**< The delay of the MUX right after the 3-LUT of a LAB */	
	double fdCtoQ_;				/**< The delay of the FlipFlop. Also contains an approximate Net Delay experimentally determined */	
	double carryInToSumOut_;	/**< The delay between the carry in and the adder outup of one LAB */
};
#endif
