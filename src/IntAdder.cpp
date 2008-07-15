/*
 * An integer adder for FloPoCo
 *
 * It may be pipelined to arbitrary frequency.
 * Also useful to derive the carry-propagate delays for the subclasses of Target
 *
 * Authors : Florent de Dinechin, Bogdan Pasca
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

/**
 * The IntAdder constructor
 * @param[in]		target		the target device
 * @param[in]		wIn			  the with of the inputs and output
 **/
IntAdder::IntAdder(Target* target, int wIn, const int p) :
	Operator(target), wIn(wIn), forcePipeline(p) {

	ostringstream name;
	name <<"IntAdder_"<<wIn;
	unique_name=name.str();

	// Set up the IO signals
	add_input ("X"  , wIn);
	add_input ("Y"  , wIn);
	add_input ("Cin", 1  );
	add_output("R"  , wIn);

	if (target->is_pipelined())
		set_sequential();
	else
		set_combinatorial();

	if (is_sequential()){
		//the maximum chunk size so that (optimisically) the target freqency can be reached. 
		chunk_size = (int)floor( (1./target->frequency() - target->lut_delay()) / target->carry_propagate_delay()); 
		//the pipe_levels gives th number of chunks that the addition must be split in so that the frequency will ve reached
		
		//set the pipeline depth 
		if (wIn<=chunk_size)
			pipe_levels=0;
		else
			pipe_levels=wIn/chunk_size;
			
		set_pipeline_depth(pipe_levels+p); 
		if(verbose) {
			cout << tab <<"Estimated delay will be " << target->adder_delay(wIn) <<endl; 
			cout << tab << "chunk="<<chunk_size << " freq=" << 1e-6/target->adder_delay(chunk_size) <<" levels="<<pipe_levels <<endl;
		}
		if(chunk_size > (wIn/(pipe_levels+1))+2)
			chunk_size = (wIn/(pipe_levels+1))+2;
		if(verbose)
			cout << tab << "after chunk balancing, chunk=="<<chunk_size << " freq=" << 1e-6/target->adder_delay(chunk_size) <<"  levels="<<pipe_levels <<endl;
		last_chunk_size = wIn - (pipe_levels)*chunk_size;
		if(verbose)
			cout << tab << "last chunk=="<<last_chunk_size <<endl;

		//if no pipelining is required - then these signals are not present
		if (pipe_levels>0)
			for(int i=0; i<=pipe_levels; i++){
				ostringstream snamex, snamey, snamer, cOutR, iSum;
				int size;
				
				snamex <<"ix_"<<i;
				snamey <<"iy_"<<i;
				iSum<<"iSum"<<i;	//isum will contain the 1 bit extended sum of ix and iy
				snamer <<"ir_"<<i;
				cOutR<<"cOutR"<<i;
				
				if(i<pipe_levels)
					size = chunk_size;
				else
					size = last_chunk_size;
							
				add_delay_signal_bus(snamex.str(), size, i);
				add_delay_signal_bus(snamey.str(), size, i);
				
				if (i<pipe_levels)	
					add_signal(iSum.str(),size+1); //is larger because inputs are padded with one 0 each in MSB position
				else
					add_signal(iSum.str(),size);
				
				add_delay_signal_bus(snamer.str(), size, pipe_levels-i);
				
				if (i<pipe_levels)
					add_registered_signal_with_sync_reset(cOutR.str(), 1);
			}	
		
	}
	if (p==1)
	add_registered_signal_with_sync_reset("r1reg",wIn);
	
}



/**
 *  destructor
 */
IntAdder::~IntAdder() {
}


/**
 * Method belonging to the Operator class overloaded by the IntAdder class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void IntAdder::output_vhdl(std::ostream& o, std::string name) {
	ostringstream signame;
	Licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);

	o << "architecture arch of " << name  << " is" << endl;
	
	output_vhdl_signal_declarations(o);

	o << "begin" << endl;
	if(is_sequential()){
		if(pipe_levels==0){
			//addition simplification
			if (forcePipeline==0)
			o << tab << "R <= X + Y + Cin;" <<endl;      
			else
			{
			o<< tab <<" r1reg <= X + Y + Cin;"<<endl;
			o<< tab <<" R<=r1reg_d;"<<endl;
			}
			
		}
		else{
			// Initialize the chunks
			for(int i=0; i<=pipe_levels; i++){
				if (!((i==pipe_levels)&&(last_chunk_size==0)))
				{
					int maxIndex;
					if(i==pipe_levels)
						maxIndex = wIn-1;
					else 
						maxIndex = i*chunk_size+chunk_size-1;
						
					ostringstream snamex, snamey, snamer;
					snamex <<"ix_"<<i;
					snamey <<"iy_"<<i;
					snamer <<"ir_"<<i;
				
					o << tab << snamex.str() << " <= ";
					if(i < pipe_levels)
						o << "";
						o << "X(" << maxIndex<< " downto " << i*chunk_size << ");" << endl;
				
					o << tab << snamey.str() << " <= ";
					if(i < pipe_levels)
						o << "";
						o << "Y(" << maxIndex<< " downto " << i*chunk_size << ");" << endl;
				}
			}
			// then the pipe_level adders of size chunk_size
			for(int i=0; i<=pipe_levels; i++){
				if (!((i==pipe_levels)&&(last_chunk_size==0)))
				{
					int size;
					if(i==pipe_levels) 
						size= last_chunk_size -1 ;
					else  
						size = chunk_size;
					
					ostringstream snamex, snamey, snamer, iSum, cOutR, cOutRM;
					snamex <<"ix_"<<i;
					snamey <<"iy_"<<i;
					snamer <<"ir_"<<i;
					iSum << "iSum"<<i;
					cOutR << "cOutR"<<i;
					cOutRM<< "cOutR"<<i-1;
					
					o<<endl<<endl;
					
					o << tab << iSum.str() << " <= ";
					if(i!=pipe_levels)
						o<< "(\"0\" & ";
						o<< get_delay_signal_name(snamex.str(), i) ; 
					if(i!=pipe_levels)	
						o<< ")";
						o<< " + ";
					if(i!=pipe_levels)	
						o<< "(\"0\" & ";
						o<< get_delay_signal_name(snamey.str(), i);
					if(i!=pipe_levels)	
						o<< ")";
					
					// add the carry in
					if(i>0) {
						ostringstream carry;
						carry <<"ir_"<<i-1;
						o  << " + " << get_delay_signal_name(cOutRM.str() , 1);
					}else
						o  << " + Cin";
					
					o << ";" << endl;
					
					if(i<pipe_levels) 
						o << tab<<snamer.str() <<" <= "<< iSum.str()<<"("<<size-1<<" downto 0);"<<endl;
					else 
						o << tab<<snamer.str() <<" <= "<< iSum.str()<<"("<<size<<" downto 0);"<<endl;
						
					if(i<pipe_levels) 
						o << tab<< cOutR.str() <<" <= "<< iSum.str()<<"("<<size<<");"<<endl;
				}
			}
			// Then the output to R 
			for(int i=0; i<=pipe_levels; i++){
				if (!((i==pipe_levels)&&(last_chunk_size==0)))
				{
					int maxIndex, size;
					if(i==pipe_levels) {
						maxIndex = wIn-1;
						size=last_chunk_size;
					}
					else  {
						maxIndex = i*chunk_size+chunk_size-1;
						size = chunk_size;
					}
					if (forcePipeline==0)
					{
						ostringstream snamer;
						snamer <<"ir_"<<i;
						o << tab << "R(" << maxIndex << " downto " << i*chunk_size << ")  <=  "
							<<  get_delay_signal_name(snamer.str(), pipe_levels -i) << "(" << size-1 << " downto 0);" << endl; 
					}
					else
					{
					ostringstream snamer;
						snamer <<"ir_"<<i;
						o << tab << "r1reg_d(" << maxIndex << " downto " << i*chunk_size << ")  <=  "
							<<  get_delay_signal_name(snamer.str(), pipe_levels -i) << "(" << size-1 << " downto 0);" << endl; 
					o<<tab<<"R<=r1reg_d;"<<endl;
					}
					
				}
			}
		}
		output_vhdl_registers(o);
	}
	else{
		o << tab << "R <= X + Y + Cin;" <<endl;
	}
	o << "end architecture;" << endl << endl;
}



TestIOMap IntAdder::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("X"));
	tim.add(*get_signal_by_name("Y"));
	tim.add(*get_signal_by_name("Cin"));
	tim.add(*get_signal_by_name("R"));
	return tim;
}

void IntAdder::fillTestCase(mpz_class a[])
{
	mpz_class& svX = a[0];
	mpz_class& svY = a[1];
	mpz_class& svC = a[2];
	mpz_class& svR = a[3];

	svR = svX + svY + svC;
	// Don't allow overflow
	mpz_clrbit(svR.get_mpz_t(),wIn); 
	
}


