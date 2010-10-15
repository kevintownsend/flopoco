/*
 * A model of FPGA that works well enough for Virtex-5 chips  (LX, speedgrade -3, e.g. xc5vlx30-3) 
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

#ifndef VIRTEX5_HPP
#define VIRTEX5_HPP
#include "../Target.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>


namespace flopoco{

	/** Class for representing an Virtex5 target */
	class Virtex5 : public Target
	{
	public:
		/** The default constructor. */  
		Virtex5() : Target()	{
			id_             = "Virtex5";
			vendor_         = "Xilinx";
			sizeOfBlock_ 	= 18432;	// the size of a primitive block is 2^11 * 9
			maxFrequencyMHz_= 500;
			// all these values are set more or less randomly, to match  virtex 5 more or less
			fastcarryDelay_ = 0.023e-9; //s   
			elemWireDelay_  = 0.436e-9;
			lutDelay_       = 0.086e-9; 
			multXInputs_    = 25;
			multYInputs_    = 18;
			// all these values are set precisely to match the Virtex5
			fdCtoQ_         = 0.272e-9; //the deterministic delay + an approximate NET delay
			lut2_           = 0.086e-9;
			lut3_           = 0.086e-9; //TODO
			lut4_           = 0.086e-9; //TODO
			muxcyStoO_      = 0.305e-9;
			muxcyCINtoO_    = 0.023e-9;
			ffd_            = 0.022e-9;
			muxf5_          = 0.291e-9;
			slice2sliceDelay_       = 0.436e-9;
			xorcyCintoO_    = 0.300e-9;
			lutInputs_ 		= 6;
			nrDSPs_ 		= 160; // XC5VLX30 has 1 column of 32 DSPs
			dspFixedShift_ 	= 17; 
			
			DSPMultiplierDelay_       = 2.387e-9;
			DSPAdderDelay_            = 1.889e-9;
			DSPCascadingWireDelay_        = 0.235e-9;
			DSPToLogicWireDelay_ = 0.436e-9;

			RAMDelay_ = 1.647e-9; //TODO
			RAMToLogicWireDelay_ = 0.222e-9; 
		}

		/** The destructor */
		virtual ~Virtex5() {}

		/** overloading the virtual functions of Target
		 * @see the target class for more details 
		 */
		double carryPropagateDelay();
		double adderDelay(int size);
		double comparatorDelay(int size);
		double comparatorConstDelay(int size);
		
		double DSPMultiplierDelay(){ return DSPMultiplierDelay_;}
		double DSPAdderDelay(){ return DSPAdderDelay_;}
		double DSPCascadingWireDelay(){ return DSPCascadingWireDelay_;}
		double DSPToLogicWireDelay(){ return DSPToLogicWireDelay_;}
		double LogicToDSPWireDelay(){ return DSPToLogicWireDelay_;}
		
		double RAMDelay() { return RAMDelay_; }
		double RAMToLogicWireDelay() { return RAMToLogicWireDelay_; }
		double LogicToRAMWireDelay() { return RAMToLogicWireDelay_; }
		
		void   getAdderParameters(double &k1, double &k2, int size);
		double localWireDelay(int fanout = 1);
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
		void   getDSPWidths(int &x, int &y, bool sign = false);
		int    getIntNAdderCost(int wIn, int n);	
	
	private:
		double fastcarryDelay_; /**< The delay of the fast carry chain */
		double lutDelay_;       /**< The delay between two LUTs */
		double elemWireDelay_;  /**< The elementary wire dealy (for computing the distant wire delay) */
	
		double fdCtoQ_;         /**< The delay of the FlipFlop. Also contains an approximate Net Delay experimentally determined */
		double lut2_;           /**< The LUT delay for 2 inputs */
		double lut3_;           /**< The LUT delay for 3 inputs */
		double lut4_;           /**< The LUT delay for 4 inputs */
		double muxcyStoO_;      /**< The delay of the carry propagation MUX, from Source to Out*/
		double muxcyCINtoO_;    /**< The delay of the carry propagation MUX, from CarryIn to Out*/
		double ffd_;            /**< The Flip-Flop D delay*/
		double muxf5_;          /**< The delay of the almighty mux F5*/
		double slice2sliceDelay_;       /**< This is approximate. It approximates the wire delays between Slices */
		double xorcyCintoO_;    /**< the S to O delay of the xor gate */
		int nrDSPs_;			/**< Number of available DSPs on this target */
		int dspFixedShift_;		/**< The amount by which the DSP block can shift an input to the ALU */
		
		double DSPMultiplierDelay_;
		double DSPAdderDelay_;
		double DSPCascadingWireDelay_;
		double DSPToLogicWireDelay_;
		
		double RAMDelay_;
		double RAMToLogicWireDelay_;

	};

}
#endif
