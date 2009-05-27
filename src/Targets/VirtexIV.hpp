/*
 * A model of FPGA that works well enough for Virtex-4 chips. 
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

#ifndef VIRTEXIV_HPP
#define VIRTEXIV_HPP
#include "../Target.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

/** Class for representing an VirtexIV target */
class VirtexIV : public Target
{
public:
	/** The default constructor. */  
	VirtexIV() : Target()	{
		// all these values are set more or less randomly, to match  virtex 4 more or less
		fastcarryDelay_ = 3.4e-11; //s   
		elemWireDelay_  = 0.3e-11;
		lut2lutDelay_   = 1.5e-10;
		lutDelay_       = 1.5e-9; 
		multXInputs_    = 18;
		multYInputs_    = 18;
		// all these values are set precisely to match the Virtex4
		fdCtoQ_         = 0.272e-9; //the deterministic delay + an approximate NET delay
		lut2_           = 0.147e-9;
		lut3_           = 0.147e-9; //TODO
		lut4_           = 0.147e-9; //TODO
		muxcyStoO_      = 0.278e-9;
		muxcyCINtoO_    = 0.034e-9;
		ffd_            = 0.017e-9;
		muxf5_          = 0.291e-9;
		netDelay_       = 0.436e-9;
		xorcyCintoO_    = 0.273e-9;
	}

	/** The destructor */
	virtual ~VirtexIV() {}

	/** overloading the virtual functions of Target
	 * @see the target class for more details 
	 */
	double carryPropagateDelay();
	double adderDelay(int size);
	double localWireDelay();
	double lutDelay();
	double distantWireDelay(int n);
	bool   suggestSubmultSize(int &x, int &y, int wInX, int wInY);
	bool   suggestSubaddSize(int &x, int wIn);
	bool   suggestSlackSubaddSize(int &x, int wIn, double slack);

private:
	double fastcarryDelay_; /**< The delay of the fast carry chain */
	double lut2lutDelay_;   /**< The delay between two LUTs */
	double elemWireDelay_;  /**< The elementary wire dealy (for computing the distant wire delay) */
	double lutDelay_;       /**< The LUT delay (in seconds)*/
	
	// Added by Bogdan
	double fdCtoQ_;         /**< The delay of the FlipFlop. Also contains an approximate Net Delay experimentally determined */
	double lut2_;           /**< The LUT delay for 2 inputs */
	double lut3_;           /**< The LUT delay for 3 inputs */
	double lut4_;           /**< The LUT delay for 4 inputs */
	double muxcyStoO_;      /**< The delay of the carry propagation MUX, from Source to Out*/
	double muxcyCINtoO_;    /**< The delay of the carry propagation MUX, from CarryIn to Out*/
	double ffd_;            /**< The Flip-Flop D delay*/
	double muxf5_;          /**< The delay of the almighty mux F5*/
	double netDelay_;       /**< This is approximate. It approximates the wire delays between Slices */
	double xorcyCintoO_;    /**< the S to O delay of the xor gate */
};

#endif
