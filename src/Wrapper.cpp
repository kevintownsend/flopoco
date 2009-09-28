/*
 * A generic wrapper generator for FloPoCo. 
 *
 * A wrapper is a VHDL entity that places registers before and after
 * an operator, so that you can synthesize it and get delay and area,
 * without the synthesis tools optimizing out your design because it
 * is connected to nothing.
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
#include "Operator.hpp"
#include "Wrapper.hpp"

namespace flopoco{

	Wrapper::Wrapper(Target* target, Operator *op):
		Operator(target), op_(op)
	{
		setCopyrightString("Florent de Dinechin (2007)");
		/* the name of the Wrapped operator consists of the name of the operator to be 
			wrapped followd by _Wrapper */
		setName(op_->getName() + "_Wrapper");
		
		//this operator is a sequential one	even if Target is unpipelined
		setSequential();	
	

		// Copy the signals of the wrapped operator 
		// This replaces addInputs and addOutputs
		for(int i=0; i < op->getIOListSize(); i++)	{
			Signal* s = op->getIOListSignal(i);
			if(s->type() == Signal::in) 
				addInput(s->getName(), s->width());
			if(s->type() == Signal::out) 
				addOutput(s->getName(), s->width());
		}

		
 		string idext;

		// copy inputs
		for(int i=0; i < op->getIOListSize(); i++){
			Signal* s = op->getIOListSignal(i);
			 if(s->type() == Signal::in) {
				 idext = "i_"+s->getName();
				 vhdl << tab << declare(idext, s->width()) << " <= " << s->getName() << ";" << endl;
			}
		}		

		// register inputs
		setCycle(1);

		for(int i=0; i < op->getIOListSize(); i++){
			Signal* s = op->getIOListSignal(i);
			 if(s->type() == Signal::in) {
				 idext = "i_"+s->getName();
				 inPortMap (op, s->getName(), idext);
			 }
		}		


		// port map the outputs
		for(int i=0; i < op->getIOListSize(); i++){
			Signal* s = op->getIOListSignal(i);
			if(s->type() == Signal::out) {
				idext = "o_" + s->getName();
				outPortMap (op, s->getName(), idext);
			}
		}

		// The VHDL for the instance
		vhdl << instance(op, "test");

		// Advance cycle to the cycle of the outputs
		syncCycleFromSignal(idext, false); // this is the last output
		nextCycle();

		// copy the outputs
		for(int i=0; i < op->getIOListSize(); i++){
			Signal* s = op->getIOListSignal(i);
			if(s->type() == Signal::out) {
				string idext = "o_" + s->getName();
				vhdl << tab << s->getName() << " <= " << use(idext) << ";" <<endl;
			}
		}
		
	}

	Wrapper::~Wrapper() {
	}


}
