/*
 * A generic class for tables of values
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

#ifndef TABLE_HPP
#define TABLE_HPP
#include <gmpxx.h>

#include "Operator.hpp"

/** A basic hardware look-up table for FloPoCo. 

 If the input to your table are negative, etc, or if you want to
 define errors, or... then derive a class from this one.

*/



class Table : public Operator
{
 public:

	/** Input width (in bits)*/
	int wIn;

	/** Output width (in bits)*/
	int wOut;

	/** minimal input value (default 0) */
	int minIn; 

	/** maximal input value (default 2^wIn-1) */
	int maxIn; 
	
	/**
	 * The Table constructor
	 * @param[in] target the target device
	 * @param[in] wIn    the with of the input in bits
	 * @param[in] wOut   the with of the output in bits  
	 **/
	Table(Target* target, int _wIn, int _wOut, int _minIn=0, int _maxIn=-1);

	virtual ~Table() {};



	/** The function that will define the values contained in the table
	 * @param[in] x  input to the table, an integer value between minIn and maxIn
	 * @return    an mpz integer  between 0 and 2^wOut-1 
	 */
	virtual mpz_class function(int x) =0;


	void outputVHDL(ostream& o, string name);

	/** A function that translates an real value into an integer input
		 This function should be overridden by an implementation of Table
	*/
	virtual int    double2input(double x);

	/** A function that translates an integer input value into its real semantics
		 This function should be overridden by an implementation of Table
	*/
	virtual double input2double(int x);

	/** A function that translates an real value into an integer output
		 This function should be overridden by an implementation of Table
	*/
	virtual  mpz_class double2output(double x);

	/** A function that translates an integer output value into its real semantics
		 This function should be overridden by an implementation of Table
	*/
	virtual double output2double(mpz_class x);
	
	/** A function that returns an estimation of the size of the table in LUTs. Your mileage may vary thanks to boolean optimization */
	int size_in_LUTs();
 private:
	bool full; /**< true if there is no "don't care" inputs, i.e. minIn=0 and maxIn=2^wIn-1 */
};


#endif
