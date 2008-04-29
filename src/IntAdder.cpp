/*
 * An integer adder for FloPoCo

 It may be pipelined to arbitrary frequency.
 Also useful to derive the carry-propagate delays for the subclasses of Target

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

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "IntAdder.hpp"

using namespace std;









IntAdder::IntAdder(Target* target, int wIn) :
	Operator(target), wIn(wIn) {

	ostringstream name;
	name <<"IntAdder_"<<wIn;
	unique_name=name.str();

	// Set up the IO signals
	add_input ("X", wIn);
	add_input ("Y", wIn);
	add_output ("R", wIn);

	if (target->is_pipelined())
		set_sequential();
	else
		set_combinatorial();
 
	chunk_size = (int)floor( (1./target->frequency() - target->lut_delay()) / target->carry_propagate_delay()); // 1 if no need for pipeline
	pipe_levels = wIn/chunk_size + 1; // =1 when no pipeline
	if(pipe_levels==1){
		add_registered_signal("iR", wIn);
		set_pipeline_depth(1); 
	}

	else{
		set_pipeline_depth(pipe_levels); 
		//  if(c2_pipe_levels) c2_chunk_size=sizeSummand;
		cout << tab <<"Estimated delay will be " << target->adder_delay(wIn) <<endl; 
		cout << tab << "chunk="<<chunk_size << " freq=" << 1e-6/target->adder_delay(chunk_size) <<"  levels="<<pipe_levels <<endl;
		if(chunk_size > (wIn/pipe_levels)+2)
			chunk_size = (wIn/pipe_levels)+2;
		cout << tab << "after chunk balancing, chunk=="<<chunk_size << " freq=" << 1e-6/target->adder_delay(chunk_size) <<"  levels="<<pipe_levels <<endl;
		last_chunk_size = wIn - (pipe_levels-1)*chunk_size;
		cout << tab << "last chunk=="<<last_chunk_size <<endl;
		for(int i=0; i<pipe_levels; i++){
			ostringstream snamex, snamey, snamer;
			snamex <<"ix_"<<i;
			snamey <<"iy_"<<i;
			snamer <<"ir_"<<i;
			int size;
			if(i<pipe_levels-1)
				size = chunk_size+1;
			else
				size = last_chunk_size;
			add_delay_signal(snamex.str(), size, i);
			add_delay_signal(snamey.str(), size, i);
			add_delay_signal(snamer.str(), size, pipe_levels-i);
		}
	}

}


IntAdder::~IntAdder() {
}







void IntAdder::output_vhdl(std::ostream& o, std::string name) {
	ostringstream signame;
	Licence(o,"Florent de Dinechin (2007)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);

	o << "architecture arch of " << name  << " is" << endl;
	
	output_vhdl_signal_declarations(o);

	o << "begin" << endl;
	if(is_sequential()){
		if(pipe_levels==1){
			o << tab << "iR <= X+Y;" <<endl;
			o << tab << "R <= iR_d;" <<endl;      
		}
		else{
			// Initialize the chunks
			for(int i=0; i<pipe_levels; i++){
				int maxIndex;
				if(i==pipe_levels-1)
					maxIndex = wIn-1;
				else 
					maxIndex = i*chunk_size+chunk_size-1;
				ostringstream snamex, snamey, snamer;
				snamex <<"ix_"<<i;
				snamey <<"iy_"<<i;
				snamer <<"ir_"<<i;
				o << tab << snamex.str() << " <= ";
				if(i < pipe_levels-1)
					o << "\"0\" & ";
				o << "X(" << maxIndex<< " downto " << i*chunk_size << ");" << endl;
				o << tab << snamey.str() << " <= ";
				if(i < pipe_levels-1)
					o << "\"0\" & ";
				o << "Y(" << maxIndex<< " downto " << i*chunk_size << ");" << endl;
			}
			// then the pipe_level adders of size chunk_size
			for(int i=0; i<pipe_levels; i++){
				int size;
				if(i==pipe_levels-1) 
					size=last_chunk_size -1 ;
				else  
					size = chunk_size;
				ostringstream snamex, snamey, snamer;
				snamex <<"ix_"<<i;
				snamey <<"iy_"<<i;
				snamer <<"ir_"<<i;
				o << tab << snamer.str() << " <= "
					<< get_delay_signal_name(snamex.str(), i)  
					<< " + "
					<< get_delay_signal_name(snamey.str(), i);
				// add the carry in
				if(i>0) {
					ostringstream carry;
					carry <<"ir_"<<i-1;
					o  << " + ((" << size  << " downto 1 => '0') & " << get_delay_signal_name(carry.str(), 1) << "(" << chunk_size << "))";
				}
				o << ";" << endl;
			}
			// Then the output to R 
			for(int i=0; i<pipe_levels; i++){
				int maxIndex, size;
				if(i==pipe_levels-1) {
					maxIndex = wIn-1;
					size=last_chunk_size;
				}
				else  {
					maxIndex = i*chunk_size+chunk_size-1;
					size = chunk_size;
				}
				ostringstream snamer;
				snamer <<"ir_"<<i;
				o << tab << "R(" << maxIndex << " downto " << i*chunk_size << ")  <=  "
					<<  get_delay_signal_name(snamer.str(), pipe_levels -i) << "(" << size-1 << " downto 0);" << endl;
			}
		}
		output_vhdl_registers(o);
	}
	else{
		o << tab << "R <= X+Y;" <<endl;
	}
	o << "end architecture;" << endl << endl;
}



