/*
 * A model of Stratix IV FPGA optimized for (EP4SE230F29C2 speed grade 2)
 *
 * Author : Florent de Dinechin, Sebastian Banescu
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

#ifndef STRATIXIV_HPP
#define  STRATIXIV_HPP
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../Target.hpp"


namespace flopoco{

	/** Class for representing an StratixIV target */
	class StratixIV : public Target
	{
	public:

		/** The default constructor. */  
		StratixIV() : Target()	{
			id_                 = "StratixIV";
			sizeOfBlock_ 		= 4608; 	// the size of a primitive block is 2^9 * 9
			fastcarryDelay_ 	= 0.9e-11; 	// *aproximately right    
			elemWireDelay_  	= 0.105e-11;	// *obtained from Quartus 2 Chip Planner
			lut2lutDelay_   	= 1.5e-10;	// ???
			lutDelay_       	= 0.355e-9; 	// 355 ps  
			ffDelay_        	= 0.069e-9; 	// *obtained from Quartus 2 Chip Planner
			multXInputs_    	= 36;
			multYInputs_    	= 36;
			lutInputs_		= 6;
			almsPerLab_		= 10;		// there are 10 ALMs per LAB
			// all these values are set precisely to match the Stratix 2
			lut2_ 			= 0.184e-9; 	// *obtained from Quartus 2 Chip Planner
			lut3_			= 0.228e-9; 	// *obtained from Quartus 2 Chip Planner
			lut4_			= 0.235e-9; 	// *obtained from Handbook
			innerLABcarryDelay_	= 0.067e-9; 	// *obtained from Quartus 2 Chip Planner
			interLABcarryDelay_	= 0.129e-9; 	// *obtained from Quartus 2 Chip Planner
			shareOutToCarryOut_	= 0.172e-9; 	// obtained from Quartus 2 Chip Planner
			muxStoO_		= 0.189e-9; 	// obtained from Quartus 2 Chip Planner by subtraction
			fdCtoQ_			= 0.214e-9; 	// *obtained from Quartus 2 Chip Planner by subtraction
			carryInToSumOut_	= 0.291e-9;	// *obtained from Quartus 2 Chip Planner
			// DSP parameters
			totalDSPs_ 		= 161;		
			nrConfigs_ 		= 4;		// StratixIV has 9, 12, 18, 36 bit multipliers by default
		
			multiplierWidth_[0] 	= 9;
			multiplierWidth_[1] 	= 12;
			multiplierWidth_[2] 	= 18;
			multiplierWidth_[3] 	= 36;
		
			multiplierDelay_[0] 	= 3.156e-9; 	// *obtained experimentaly from Quartus 2
			multiplierDelay_[1] 	= 3.069e-9; 	// *obtained experimentaly from Quartus 2
			multiplierDelay_[2] 	= 2.744e-9; 	// *obtained experimentaly from Quartus 2
			multiplierDelay_[3] 	= 3.604e-9; 	// *obtained experimentaly from Quartus 2
		}
	
		/** The destructor */
		virtual ~StratixIV() {}

		/** overloading the virtual functions of Target
		 * @see the target class for more details 
		 */
		double carryPropagateDelay();
		double adderDelay(int size);
		void   getAdderParameters(double &k1, double &k2, int size);
		double localWireDelay();
		double lutDelay();
		double ffDelay();
		double distantWireDelay(int n);
		bool   suggestSubmultSize(int &x, int &y, int wInX, int wInY);
		bool   suggestSubaddSize(int &x, int wIn);
		bool   suggestSlackSubaddSize(int &x, int wIn, double slack);
		int    getIntMultiplierCost(int wInX, int wInY);
		long   sizeOfMemoryBlock();
		DSP*   createDSP(); 
		int    getEquivalenceSliceDSP();
		int    getNumberOfDSPs();
		void   getDSPWidths(int &x, int &y);
		int    getIntNAdderCost(int wIn, int n);
	
	private:

		double fastcarryDelay_; 	/**< The delay of the fast carry chain */
		double lut2lutDelay_;   	/**< The delay between two LUTs */
		double ffDelay_;   		/**< The delay between two flipflops (not including elemWireDelay_) */
		double elemWireDelay_;  	/**< The elementary wire dealy (for computing the distant wire delay) */
		double lutDelay_;      	 	/**< The LUT delay (in seconds)*/
	
		// Added by Sebi
		double lut2_;           	/**< The LUT delay for 2 inputs */
		double lut3_;           	/**< The LUT delay for 3 inputs */
		double lut4_;           	/**< The LUT delay for 4 inputs */
		double innerLABcarryDelay_;	/**< The wire delay between the upper and lower parts of a LAB --> R4 & C4 interconnects */	
		double interLABcarryDelay_;	/**< The approximate wire between two LABs --> R24 & C16 interconnects */	
		double shareOutToCarryOut_;	/**< The delay between the shared arithmetic out of one LAB and the carry out of the following LAB */	
		double muxStoO_;		/**< The delay of the MUX right after the 3-LUT of a LAB */	
		double fdCtoQ_;			/**< The delay of the FlipFlop. Also contains an approximate Net Delay experimentally determined */	
		double carryInToSumOut_;	/**< The delay between the carry in and the adder outup of one LAB */
		int    almsPerLab_;		/**< The number of ALMs contained by a LAB */
	
		// DSP parameters
		int	totalDSPs_;		/**< The total number of DSP blocks available on this target */	
		int	nrConfigs_;		/**< The number of distinct predefinded multiplier widths */
		int 	multiplierWidth_[4];	/**< The multiplier width available */
		double 	multiplierDelay_[4];	/**< The corresponding delay for each multiplier width available */
		double	inputRegDelay_[4];	/**< The input register delay to DSP block for each multiplier width available */
		double	pipe2OutReg2Add; 	/**< The DPS block pipeline register to output register delay in two-multipliers adder mode */
		double  pipe2OutReg4Add; 	/**< The DPS block pipeline register to output register delay in four-multipliers adder mode */
	
	};
}
#endif
