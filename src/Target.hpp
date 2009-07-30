/*
 * The abstract class that models different chips for delay and area. 
 * Classes for real chips inherit from this one. They should be in subdirectory Targets
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

#ifndef TARGET_HPP
#define TARGET_HPP

#include "DSP.hpp"

/** Abstract target Class. All classes which model real chips inherit from this class */
class Target
{
 public:
 	/** The default constructor. Creates a pipelined target, with 4 inputs/LUT, 
 	 * with a desired frequency of 400MHz and which is allowed to use hardware multipliers
 	 */ 
	Target()   {
		pipeline_          = true;
		lutInputs_         = 4;
		frequency_         = 400000000.;
		useHardMultipliers_= true;
	}
	
	/** The destructor */
	virtual ~Target() {}

	// Architecture-related methods
	/** Returns the number of inputs that the LUTs have on the specifi device
	 * @return the number of inputs for the look-up tables (LUTs) of the devie
	 */
	int lutInputs();
	
	/** Function for determining the submultiplication sizes so that the design is able 
	 * to function at a desired frequency ( multiplicatio is considerd X * Y )
	 * @param[in,out] x the size of the submultiplication for x
	 * @param[in,out] y the size of the submultiplication for y
	 * @param[in] wInX the width of X
	 * @param[in] wInY the width of Y
	 */
	virtual bool suggestSubmultSize(int &x, int &y, int wInX, int wInY)=0;

	/** Function for determining the subadition sizes so that the design is able 
	 * to function at a desired frequency ( addition is considerd X + Y )
	 * @param[in,out] x the size of the subaddition for x and y
	 * @param[in] wIn the widths of X and Y
	 */
	virtual bool suggestSubaddSize(int &x, int wIn)=0; 	

	/** Function for determining the subadition sizes so that the design is able 
	 * to function at a desired frequency ( addition is considerd X + Y )
	 * @param[in,out] x the size of the subaddition for x and y
	 * @param[in] wIn the widths of X and Y
	 * @param[in] slack the time delay consumed out of the input period 
	 */
	virtual bool   suggestSlackSubaddSize(int &x, int wIn, double slack)=0;

	// Delay-related methods
	/** Function which returns the lutDelay for this target
	 * @return the LUT delay
	 */
	virtual double lutDelay() =0;
	
	/** Function which returns the flip-flop Delay for this target
		 (not including any net delay)
	 * @return the flip-flop delay
	 */
	virtual double ffDelay() =0;
	
	/** Function which returns the carry propagate delay
	 * @return the carry propagate dealy
	 */
	virtual double carryPropagateDelay() =0;

	/** Function which returns addition delay for an n bit addition
	 * @param n the number of bits of the addition (n-bit + n-bit )
	 * @return the addition delay for an n bit addition
	 */
	virtual double adderDelay(int n) =0;

	/** Function which returns the local wire delay (local routing)
	 * @return local wire delay 
	 */
	virtual double localWireDelay() =0;
	
	/** Function which returns the size of a primitive memory block,which could be recognized by the sintetyzer as a dual block.
	* @return the size of the primitive memory block
	*/	
	
	virtual long sizeOfMemoryBlock() = 0 ;

	/** Function which returns the distant wire delay.
	 * @return distant wire delay 
	 */
	virtual double distantWireDelay(int n) =0;

	// Methods related to target behaviour and performance
	/** Sets the target to pipelined */
	void setPipelined();                
	
	/**< Sets the target to combinatorial */    
	void setNotPipelined();                 
	
	/** Returns true if the target is pipelined, otherwise false
	 * @return if the target is pipelined
	 */
	bool isPipelined();
	
	/** Returns the desired frequency for this target in Hz
	 * @return the frequency
	 */
	double frequency();

	/** Returns the desired frequency for this target in MHz
	 * @return the frequency
	 */
	double frequencyMHz();
	
	/** Sets the desired frequency for this target
	 * @param f the desired frequency
	 */	
	void setFrequency(double f);

	/** Sets the use of hardware multipliers 
	 * @param v use or not harware multipliers
	 */
	void setUseHardMultipliers(bool v);

	/** Returns true if the operator for this target is allowed to use hardware multipliers
	 * @return the status of hardware multipliers usage
	 */
	bool getUseHardMultipliers(); 

	/** Function which returns the number of LUTs needed to implement
	 *	 a multiplier of the given width
	 * @param	wInX the width (in bits) of the first operand
	 * @param	wInY the width (in bits) of the second operand
	 */
	virtual int multiplierLUTCost(int wInX, int wInY) =0;
	
	/** Constructs a specific DSP to each target */
	
	virtual DSP* createDSP() = 0;
	

protected:
	int    lutInputs_;          /**< The number of inputs for the LUTs */
	bool   pipeline_;           /**< True if the target is pipelined/ false otherwise */
	double frequency_;          /**< The desired frequency for the operator in Hz */
	int    multXInputs_;        /**< The size for the X dimmension of the hardware multipliers */
	int    multYInputs_;        /**< The size for the Y dimmension of the hardware multipliers */
	bool   useHardMultipliers_; /**< If true the operator design can use the hardware multipliers */
	long sizeOfBlock;		/**<The size of a primitive memory block> */
	
};


#endif
