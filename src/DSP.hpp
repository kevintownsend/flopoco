/*
 * The class that models different digital signal processors. 
 * Classes for real chips instantiate this one giving it apropriate parameter values.
 *
 * Author : Sebastian Banescu, Radu Tudoran
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

#ifndef DSP_HPP
#define DSP_HPP

class DSP
{
 public:
 	/** The default constructor. Creates a DSP configuration valid for all targets: 
	 * 18x18 multipliers, no shifter, using one adder inside the block and is
	 * not configured to do multiply-accumulate
 	 */ 
	DSP()   {
			multiplierWidth_ = 18;
			fixedShift_ = 0;
			nrAdders_ = 1;
			multAccumulate_ = false;
	}
	
	/** The destructor */
	virtual ~DSP() {}

	
protected:
	int    multiplierWidth_;	/**< The recommended with of the multiplier used by this DSP */
	int    fixedShift_;         /**< The shift capacity of one input of the DSP block */
	int    nrAdders_;          	/**< The number of adders used by this DSP reflects the number of multipliers used and the number of outputs */
	bool   multAccumulate_;     /**< If true the DSP block will be configured as Multiply-Accumulate */
	
	// attributes used for tilling algorithm
	int 	positionX;		/**< The X coordinates for the 2 points that mark the position of the DSP block inside the tiling */
	int		positionY;		/**< The Y coordinates for the 2 points that mark the position of the DSP block inside the tiling */
	DSP*   	shiftIn;        	/**< The DSP block from which this block obtains a shifted operand (Virtex4) */
	DSP*	shiftOut;			/**< The DSP block to which this block provides a shifted operand (Virtex4) */
	DSP*	addOperands;		/**< The DSP blocks whose multiplication results may be added to this one depending on the number of adders used (StratixII)*/
};


#endif
