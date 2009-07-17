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

#ifndef DUALTABLE_HPP
#define DUALTABLE_HPP
#include <gmpxx.h>

#include "Operator.hpp"
#include <math.h>
#include <sstream>
#include <cstdlib>

/** A basic hardware look-up table for FloPoCo. 

 If the input to your table are negative, etc, or if you want to
 define errors, or... then derive a class from this one.

*/

using namespace std;

class DualTable : public Operator
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
		
	Target *target;
 
	/**
	 * The Table constructor
	 * @param[in] target the target device
	 * @param[in] wIn    the with of the input in bits
	 * @param[in] wOut   the with of the output in bits  
	 **/
	DualTable(Target* target, int _wIn, int _wOut, int _minIn=0, int _maxIn=-1);

	virtual ~DualTable() {};

	
	/** The function that will define the values contained in the table
	 * @param[in] x  input to the table, an integer value between minIn and maxIn
	 * @return    an mpz integer  between 0 and 2^wOut-1 
	 */
	
	virtual mpz_class function(int x) =0;

	/** Overloading the method of Operator */
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
	
	/**
	 The inner class of the DualTable which is the fundamental entity which is recognized as a dual port memory.
	*/
	class primitiveDualMemory: public Operator{
		
		public:
		
		/** Input width (in bits)*/
		int wIn;

		/** Output width (in bits)*/
		int wOut;

		/** minimal input value (default 0) */
		int minIn; 

		/** maximal input value (default 2^wIn-1) */
		int maxIn; 
	
		stringstream inputData;
		
		void setInputData(std::ostringstream&  data)
		{
			inputData.str("");
			inputData<<data.str();
		}
		
				
	primitiveDualMemory(Target* target, int _wIn, int _wOut, int _minIn, int _maxIn) : 
	Operator(target),
	wIn(_wIn), wOut(_wOut), minIn(_minIn), maxIn(_maxIn) 
	{
	setCopyrightString("Radu Tudoran (2009)");
	
	// Set up the IO signals

	addInput ("X1"  , wIn);
	addOutput ("Y1"  , wOut);
	addInput ("X2"  , wIn);
	addOutput ("Y2"  , wOut);			
		
	}

	~primitiveDualMemory() {}
		
	
	void outputVHDL(std::ostream& o, std::string name) {
	
	int i,x;
	mpz_class y;
		
	vhdl<<tab<<" process (clk)"<<endl;
	vhdl<<tab<<" variable memVar: ROM_Array := ( "<<endl;
	
	int count=0;
	for (x = minIn; x <= maxIn; x++) {
		inputData>>y;
		vhdl 	<< tab<<tab << "\"" << unsignedBinary(y, wOut)  ;
		if(x!=maxIn)
			vhdl<< "\", ";
		else
			vhdl<<"\" ";
		count++;
		if(count==4)
		{vhdl<<endl;
			count=0;
		}
	}
	vhdl<<");"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if clk'event and clk = '1' then "<<endl;
	vhdl<<tab<<tab<<tab<<" Y1 <= memVar(conv_integer(X1)); "<<endl;
	vhdl<<tab<<tab<<tab<<" Y2 <= memVar(conv_integer(X2)); "<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
	
	licence(o);
	stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);

	o<<"type ROM_Array is array (0 to "<< intpow2(wIn)-1<< ") of std_logic_vector("<<wOut-1 <<" downto 0);"<<endl;
	
	o << buildVHDLComponentDeclarations();	
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);		
	
	o<<buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);
		
	}
		
		private:
			
		
	};


	bool full; /**< true if there is no "don't care" inputs, i.e. minIn=0 and maxIn=2^wIn-1 */
	int nrOfMemBlocks;
	int limitSingleMemory ;
	
	primitiveDualMemory * primitiveBlock;
	primitiveDualMemory ** primitiveBlocks ;
	
};


#endif
