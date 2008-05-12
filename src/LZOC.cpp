/*
 * A leading zero/one counter for FloPoCo
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

#include "LZOC.hpp"

using namespace std;

/** The LZOC constructor
*@param[in] target the target device for this operator
*@param[in] wIn the width of the input
*@param[in[ wOut the width of the output 
*/
LZOC::LZOC(Target* target, int wIn, int wOut) :
	Operator(target), wIn(wIn), wOut(wOut)  {

	p2wOut = 1<<wOut; // no need for GMP here
	ostringstream name; 
	name <<"LZOC_"<<wIn<<"_"<<wOut;

	unique_name=name.str();

	//Set up the IO signals
	add_input ("I", wIn);
	add_input ("OZB");  
	add_output("O", wOut);

	if (target->is_pipelined())
		set_sequential();
	else
		set_combinatorial();	

	//Set up the internal architecture signals
	if (is_sequential()){
			for (int i=wOut; i>0; i--){
			ostringstream signalName;
			signalName<<"level"<<i;
			add_delay_signal(signalName.str(), (1<<i),2);
			}
	}
	else{
		add_signal("tmpO", wOut);
	
		for (int i=wOut; i>0; i--){
			ostringstream signalName;
			signalName<<"level"<<i;
			add_signal(signalName.str(), (1<<i));
		}
	}
	 
	if (is_sequential())
	set_pipeline_depth(2*wOut); 
	 
	 
}


/** The LZOC destructor
*/
LZOC::~LZOC() {}



/**
 * Method belonging to the Operator class overloaded by the LZOC class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/

void LZOC::output_vhdl(std::ostream& o, std::string name) {
	Licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	new_architecture(o,name);
	
	output_vhdl_signal_declarations(o);	
	for (int i=wOut; i>0; i--){
		o << tab << "signal tmpO"<<i<<"  : std_logic_vector("<<wOut-1<<" downto "<<i-1<<");"<<endl;
		o << tab << "signal tmpO"<<i<<"_d: std_logic_vector("<<wOut-1<<" downto "<<i-1<<");"<<endl;				
	}
	
	
	begin_architecture(o);
	
	if (is_sequential()){
		
		output_vhdl_registers(o); o<<endl;
		
		o << tab << "process(clk, rst)"<<endl;
		o << tab << "  begin"<<endl;
		o << tab << "if clk'event and clk = '1' then"<<endl;
		o << tab << "     if rst = '1' then"<<endl;
		for (int i=wOut; i>0; i--)
		o << tab << "        tmpO"<<i<<"_d <= ("<<wOut-1<<" downto "<<i-1<<" => '0');"<<endl;				
		o << tab << "     else"<<endl;
		for (int i=wOut; i>0; i--)
		o << tab << "        tmpO"<<i<<"_d <= tmpO"<<i<<";"<<endl;				
		o << tab << "     end if;"<<endl;
		o << tab << "  end if;"<<endl;
		o << tab << "end process;"<<endl;
		
		
		
		// connect first stage to I
		cout<<"p2wOut is = "<<p2wOut<<endl;
		if (p2wOut==wIn)
			o << "  level"<<wOut<<" <=  I ;" << endl;
		else if (p2wOut>wIn) // pad input with zeroes/ones function of what we count. If LZC pad with 1, else pad with 
			o << "  level"<<wOut<<" <=  I & ("<<p2wOut-wIn-1<<" downto 0 => not(ozb)) ;" << endl;
		else if (p2wOut<wIn)
			o << "  level"<<wOut<<" <=  I("<<wIn-1<<" downto "<<wIn - p2wOut <<") ;" << endl;
		// recursive structure
		for (int i=wOut; i>1; i--) {
			if (i!=wOut)
			o << "  tmpO"<<i<<"("<<wOut-1<<" downto "<<i<<") <= tmpO"<<i+1<<"_d("<<wOut-1<<" downto "<<i<<");"<<endl;
			
			o << "  tmpO"<<i<<"("<<i-1<<") <= '1' when level"<<i<<"_d("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<") = ("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<" => ozb) else '0';"<< endl;
			o << "  level"<<i-1<<" <= level"<<i<<"_d_d("<<(1<<(i-1))-1<<" downto 0) when tmpO"<<i<<"_d("<<i-1<<") = '1'"<< endl
				<< "               else level"<<i<<"_d_d("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<");" << endl;
		}
		o << "  tmpO1(0) <= '1' when level1_d(1) =  ozb else '0';"<< endl;
		o << "  tmpO"<<1<<"("<<wOut-1<<" downto "<<1<<") <= tmpO"<<2<<"_d("<<wOut-1<<" downto "<<1<<");"<<endl;
		
		o << " O <= tmpO1_d;"<<endl;
	}
	else{
		// connect first stage to I
		cout<<"p2wOut is = "<<p2wOut<<endl;
		if (p2wOut==wIn)
			o << "  level"<<wOut<<" <=  I ;" << endl;
		else if (p2wOut>wIn) // pad input with zeroes/ones function of what we count. If LZC pad with 1, else pad with 
			o << "  level"<<wOut<<" <=  I & ("<<p2wOut-wIn-1<<" downto 0 => not(ozb)) ;" << endl;
		else if (p2wOut<wIn)
			o << "  level"<<wOut<<" <=  I("<<wIn-1<<" downto "<<wIn - p2wOut <<") ;" << endl;
		// recursive structure
		for (int i=wOut; i>1; i--) {
			o << "  tmpO("<<i-1<<") <= '1' when level"<<i<<"("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<") = ("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<" => ozb) else '0';"<< endl;
			o << "  level"<<i-1<<" <= level"<<i<<"("<<(1<<(i-1))-1<<" downto 0) when tmpO("<<i-1<<") = '1'"<< endl
				<< "               else level"<<i<<"("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<");" << endl;
		}
		o << "  tmpO(0) <= '1' when level1(1) =  ozb else '0';"<< endl;
		o << " O <= tmpO;"<<endl;
	}
	
	end_architecture(o);
}


