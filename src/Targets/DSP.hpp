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


namespace flopoco{

	class DSP
	{
	public:
		/** The default constructor. Creates a DSP configuration valid for all targets: 
		 * 18x18 multipliers, no shifter, using one adder inside the block and is
		 * not configured to do multiply-accumulate
		 */ 
		DSP()   {
			maxMultiplierWidth_ = 18;
			maxMultiplierHeight_ = 18;
			multiplierWidth_ = 0;
			multiplierHeight_ = 0;
			fixedShift_ = 0;
			nrAdders_ = 1;
			multAccumulate_ = false;
		}
	
		DSP(int Shift,int maxMultWidth,int maxMultHeight);
	
		/** The destructor */
		virtual ~DSP() {}

		/** Returns the with of the multiplier that this DSP block is using 
		 * @return the with of the multiplier that this DSP block is using */
		int getMultiplierWidth();
	
		/** Assigns the designated width to the multiplier of the DSP block
		 * @param w width assigned to the multiplier of the DSP block.
		 */
		void setMultiplierWidth(int w);
		
		/** Returns the height of the multiplier that this DSP block is using 
		 * @return the height of the multiplier that this DSP block is using */
		int getMultiplierHeight();
	
		/** Assigns the designated height to the multiplier of the DSP block
		 * @param w height assigned to the multiplier of the DSP block.
		 */
		void setMultiplierHeight(int w);
	
		/** Returns the amount by which an input can be shifted inside the DSP block 
		 * @return the amount by which an input can be shifted inside the DSP block */
		int getShiftAmount();
	
		/** Returns the maximum with of the multiplier that this DSP block is using 
		 * @return the maximum with of the multiplier that this DSP block is using */
		int getMaxMultiplierWidth();
	
		
		/** Returns the maximum height of the multiplier that this DSP block is using 
		 * @return the maximum height of the multiplier that this DSP block is using */
		int getMaxMultiplierHeight();
	
	
		/** Returns the number of adders this DSP block is using
		 * @return the number of adders this DSP block is using */
		int getNumberOfAdders();
	
		/** Assign the number of adders used by this DSP block in order to know the number of addition operands.
		 * @param o	number of addition operands
		 */
		void setNumberOfAdders(int o);
	
		/** Returns TRUE if the DSP block is configured to perform multiply-accumulate
		 * @return TRUE if the DSP block is configured to perform multiply-accumulate */
		bool isMultiplyAccumulate();
	
		/** Sets a flag that will configure the DSP block to multiply-accumulate when TRUE.
		 * @param m value assigned to multiply-accumulate flag
		 */
		void setMultiplyAccumulate(bool m);
	
		/** Procedure that return by the value of its two parameters the coordinates
		 * of the top-right corner of the DSP within a tiling.
		 * @param x the horizontal coordinate of the top-right corner
		 * @param y the vertical coordinate of the top-right corner
		 */
		void getTopRightCorner(int &x, int &y);
	
		/** Assigns the x and y coordinates of the top-right corner of the DSP within the tiling.
		 * @param x the horizontal coordinate of the top-right corner
		 * @param y the vertical coordinate of the top-right corner
		 */
		void setTopRightCorner(int x, int y);
	
		/** Procedure that return by the value of its two parameters the coordinates
		 * of the bottom-left corner of the DSP within a tiling.
		 * @param x the horizontal coordinate of the bottom-left corner
		 * @param y the vertical coordinate of the bottom-left corner
		 */
		void getBottomLeftCorner(int &x, int &y);
	
		/** Assigns the x and y coordinates of the bottom-left corner of the DSP within the tiling.
		 * @param x the horizontal coordinate of the bottom-left corner
		 * @param y the vertical coordinate of the bottom-left corner
		 */
		void setBottomLeftCorner(int x, int y);
	
		/** Returns a reference to the DSP object whose output is the shifted input
		 * for this DSP block. This method will be mainly used for Virtex architecture. 
		 * @return a reference to the DSP object whose output is the shifted input for this DSP block*/
		DSP* getShiftIn();
	
		/** Assign the DSP object whose output is the shifted input for this DSP block.
		 *  This method will be mainly used for Virtex architecture.
		 * @param d the DSP object whose output is the shifted input for this DSP block.
		 */
		void setShiftIn(DSP* d);
	
		/** Returns a reference to the DSP object whose shifted input is the output
		 * of this DSP block. This method will be mainly used for Virtex architecture 
		 * @return a reference to the DSP object whose shifted input is the output of this DSP block*/
		DSP* getShiftOut();
	
		/** Assign the DSP object whose shifted input is the output of this DSP block.
		 *  This method will be mainly used for Virtex architecture.
		 * @param d the DSP object whose shifted input is the output of this DSP block.
		 */
		void setShiftOut(DSP* d);
	
		/** Returns references to all DSP objects that must be added with this DSP.
		 * @return  references to all DSP objects that must be added with this DSP */
		DSP** getAdditionOperands();
	
		/** Assigns an array containing the operands of the addition that this DSP is part of.
		 * @param o the array containing the operands of the addition that this DSP is part of.
		 */
		void setAdditionOperands(DSP** o);
	
		/** Swaps the width and height of the DSP block. This is used in case of asymetrical DSP
		 * blocks when we want to do a vertical/horizontal replacement of the block and we need
		 * to rotate the block by 90 degrees.
		 */
		void rotate();
		
		void allocatePositions(int dimension);
		void push(int X,int Y);
		int pop();
		
		
		int *Xpositions;
		int *Ypositions;
		
	
	protected:
		int    multiplierWidth_;	/**< The recommended with of the multiplier used by this DSP */
		int    maxMultiplierWidth_;	/**< The recommended with of the multiplier used by this DSP */
		int    multiplierHeight_;	/**< The recommended with of the multiplier used by this DSP */
		int    maxMultiplierHeight_;	/**< The recommended with of the multiplier used by this DSP */

		int    fixedShift_;         /**< The shift capacity of one input of the DSP block */
		int    nrAdders_;          	/**< The number of adders used by this DSP reflects the number of multipliers used and the number of outputs */
		bool   multAccumulate_;     /**< If true the DSP block will be configured as Multiply-Accumulate */
	
		// attributes used for tilling algorithm
		int 	positionX_[2];		/**< The X coordinates for the 2 points that mark the position of the DSP block inside the tiling */
		int		positionY_[2];		/**< The Y coordinates for the 2 points that mark the position of the DSP block inside the tiling */
		DSP*   	shiftIn_;        	/**< The DSP block from which this block obtains a shifted operand (Virtex4) */
		DSP*	shiftOut_;			/**< The DSP block to which this block provides a shifted operand (Virtex4) */
		DSP**	addOperands_;		/**< The DSP blocks whose multiplication results may be added to this one depending on the number of adders used (StratixII)*/
	private:
		/** Assigns the shift amount that is used within this DSP block.
		 * @param s shift amount used within this DSP block.
		 */
		void setShiftAmount(int s);
	
		
	
		int pos;
		int max_pos;

	};

}

#endif
