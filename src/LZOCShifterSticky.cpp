/*
 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
 *
 * Author : Florent de Dinechin, Bogdan Pasca
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

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "LZOCShifterSticky.hpp"

using namespace std;

/** The LZOCShifterSticky constructor
*@param[in] target the target device for this operator
*@param[in] wIn the width of the mantissa input
*@param[in] wOut the width of the mantissa output

The width of the lzo output will be floor(log2(wIn+1)). It should be accessed as operator->wCount.
*/

// TODO to be optimal for FPAdder, we have to provide a way to disable the sticky computation.

LZOCShifterSticky::LZOCShifterSticky(Target* target, int wIn, int wOut, bool compute_sticky, const int countType) :
	Operator(target), wIn(wIn), wOut(wOut), compute_sticky(compute_sticky), zoc(countType), countType(countType) {

	ostringstream name; 
	name <<"LZOCShifterSticky_"<<wIn<<"_"<<wOut<<"_"<<compute_sticky;

	unique_name=name.str();

	if (target->is_pipelined())
		set_sequential();
	else
		set_combinatorial();	

	//Set up the internal architecture signals

	// Terminology: level i is the signal on which we will test wether the 2^(i-1) leading bits are zeroes.
	// i goes from wCount to 0

	// Input will be copied to level wCount, possibly padded to the next power of two with zeroes
	int i, p2i;


	// from the end of the pipeline to the beginning
	i=0;
	size[0] = wOut; // size of the result
	while(size[i]-wOut < wIn){
		i++;
		p2i = 1<<(i-1);
		size[i] = size[i-1] + p2i;
		cout<<"size "<<i<<" is="<<size[i]<<endl;
		// Invariant: size[i] = wOut + 2^i -1
	}
	
	// the attribute that gives the number of bits of the LZO count
	wCount = i;
	// should be identical to : wCount = intlog2(wIn+1); // +1 for the case all zeroes
	//Set up the IO signals
	add_input ("I", wIn);
	if (zoc==-1) add_input ("OZb");  
	add_output("Count", wCount);
	add_output("O", wOut);
	if(compute_sticky)
		add_output("Sticky");
	
	// the first stage doesn't need to register zeroes
	size[wCount] = 1<<wCount;


	if(verbose){
		cout << "  wCount=" << wCount << "    Level sizes: ";
		for(i=0; i<=wCount; i++) cout << size[i]<<" ";
		cout <<endl;
	}

	for (int i=wCount; i>=0; i--){
		ostringstream levelName, leveldName, stickyName;
		levelName << "level"  << i;
		stickyName << "sticky" << i;
		level[i] = levelName.str();
		if (is_sequential() && i!=wCount) {
			add_registered_signal(level[i], size[i]);
			leveld[i] = level[i] + "_d" ;
			if(compute_sticky)
				add_registered_signal(stickyName.str());
		}
		else {
			add_signal(           level[i], size[i]);
			leveld[i] = level[i];
			if(compute_sticky)
				add_signal(           stickyName.str());
		}
		// The signal that holds the leading zero count
		if(i<wCount) {
			ostringstream name;
			name << "count" << i;
			if (is_sequential()) 
				add_delay_signal_no_reset(name.str(),1, i+1);				
			else 
				add_signal(name.str());
		}
	}

	if (is_sequential()){
		set_pipeline_depth(wCount); //was wCount-1
		if (zoc==-1)
		add_delay_signal_no_reset("sozb",1, wCount-1);
	}
}

/** The LZOCShifterSticky destructor
*/
LZOCShifterSticky::~LZOCShifterSticky() {}



/**
 * Method belonging to the Operator class overloaded by the LZOCShifterSticky class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/

void LZOCShifterSticky::output_vhdl(std::ostream& o, std::string name) {
	Licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	new_architecture(o,name);

	output_vhdl_signal_declarations(o);	

	// connect input to level wCount, possibly with 0 padding

	begin_architecture(o);
	o << tab << "level" << wCount<< " <= " ;
	if(wIn == size[wCount]) // no padding needed
		o << "I;" <<endl;
	else 
		o << "I & (" << size[wCount]-wIn-1 << " downto 0 => '0');" << endl ;

	if (zoc==-1)
	o<<tab<<"sozb <= OZb;"<<endl;
	
	if(compute_sticky)
		o<<tab<<"sticky" << wCount << " <= '1';"<<endl;

	for  (int i = wCount; i>=1; i--){
		int p2i = 1 << i;
		int p2io2 = 1 << (i-1);

		//=====================================
		o << tab << "count" << i-1 << " <= '1' when " << leveld[i] << "(" << size[i]-1 << " downto " << size[i]-p2io2 << ") = (" << p2io2-1 << " downto 0 => "; 
		if (zoc==-1)
			o<<get_delay_signal_name("sozb",wCount-i);
		else
			o<<"'"<<zoc<<"'";
		o<<" )   else '0';" << endl; 
		//=====================================
		// REM in the following,  size[i]-size[i-1]-1 is almost always equal to p2io2, except in the first stage
		if (i==wCount){
			o << tab << level[i-1] << " <= " ;
			o << " " << "("<<leveld[i] << "(" << size[i]/2 -1 << " downto " << 0 << ") & "<<zero_generator(size[i-1]-size[i]/2, 0)<<" )  when count" << i-1 << "='1'" << endl;
			
			if (size[i-1]-size[i]>0)
				o << tab << tab << tab <<  "else (" << leveld[i] << "(" << size[i]-1 << " downto 0) &  "<<zero_generator(size[i-1]-size[i],0)<<") ;" << endl; 
			else
				o << tab << tab << tab <<  "else " << leveld[i] << "(" << size[i]-1 << " downto "<<size[i]-size[i-1]<<");" << endl; 
		}
		else
		{
			o << tab << level[i-1] << " <= " ;
			o << " " << leveld[i] << "(" << size[i]-1 << " downto " << size[i]-size[i-1] << ")   when count" << i-1 << "='0'" << endl;
			o << tab << tab << tab <<  "else " << leveld[i] << "(" << size[i-1]-1 << " downto 0);" << endl; 
		}
		
		if(compute_sticky){
			o << tab << "sticky" << i-1 << " <= '0' when (count" << i-1 << "='1' or " << leveld[i] << "(" << size[i]-size[i-1]-1 << " downto 0) = (" << size[i]-size[i-1]-1 << " downto 0 => '0'))  else sticky" << i;
			if(is_sequential() && i<wCount) o << "_d";
			o << ";" << endl;
		}
			
		o << endl;
	}

	if(compute_sticky)
		o << tab << "sticky <= sticky0_d;" << endl;
	o << tab << "O      <= level0_d;" << endl;
	o << tab << "Count  <= ";
	for (int i=wCount-1; i >=0; i--){
		ostringstream name;
		name << "count" << i;
		if (is_sequential())
			o << get_delay_signal_name(name.str(), i+1);
		else
			o << name.str();
		if (i>0) o << " & ";
	}
	
	o << ";" << endl;


	if (is_sequential()){		
		output_vhdl_registers(o); o<<endl;
	}

	end_architecture(o);
	
}

TestIOMap LZOCShifterSticky::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("I"));
	//tim.add(*get_signal_by_name("OZb"));
	tim.add(*get_signal_by_name("Count"));
	tim.add(*get_signal_by_name("O"));
//	if (compute_sticky)
//		tim.add(*get_signal_by_name("Sticky"));
	
	return tim;
}

void LZOCShifterSticky::fillTestCase(mpz_class a[])
{

	if (countType>=0)
	{
		mpz_class& si     = a[0];
		mpz_class& scount = a[1];
		mpz_class& so     = a[2];
	
		/* Count the leading zero/one s */
		int j=(wIn-1);//the index of the MSB of the input
		int bit = (countType == 0) ? 0 : 1; //what are we counting
		for (j = (wIn-1); j >= 0; j--)
				if (mpz_tstbit(si.get_mpz_t(), j) != bit)
				break;
		
		int icount =(wIn-1)-j;
		scount = icount; //the count result
	
		//compute max value on wOut bits
		maxValue=2;
		for (int i=2;i<=wOut;i++)
		maxValue=maxValue*2;
		maxValue--;
	
		mpz_class inputValue;
		inputValue=si;
	
		//if we are counting zeros
		if (countType==0) 
			while (!((inputValue<=maxValue)&&(2*inputValue>maxValue)))
				if (inputValue>maxValue)
					inputValue=inputValue/2;
				else
					inputValue=inputValue*2;
		else
		{
			int restOfBits = wIn - icount;
			if (icount>0)
			{
			//	cout<<"icount="<<icount<<endl;
				mpz_class ones=1;
				for (int i=1;i<=icount;i++)
					ones = ones*2;
			
				ones=ones-1;
			
		//		cout<<"ones="<<ones<<endl;
						
				for (int i=1;i<=restOfBits;i++)
					ones=ones*2;
				inputValue=inputValue-ones; // the input without the leading ones
			}

			if ((wIn<=wOut) || ((wIn>wOut) && (restOfBits<wOut) ))	//shift result in place	
				for (int i=1;i<=(wOut-restOfBits);i++)
					inputValue=inputValue*2;
			else
				for (int i=1;i<=restOfBits-wOut;i++)
					inputValue=inputValue/2;
		}
		
		
		
		
		
		
		so=inputValue;
	}
	
	
}



