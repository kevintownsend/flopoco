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
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Wrapper.hpp"


/** 
 * The Wrapper constructor
 * @param[in] target the target device
 * @param[in] op the operator to be wrapped
 **/
Wrapper::Wrapper(Target* target, Operator *op):
	Operator(target), op(op)
{
  /* the name of the Wrapped operator consists of the name of the operator to be 
  wrapped followd by _Wrapper */
  set_operator_name("","Wrapper");
		
	if (!target->is_pipelined()) 	
		set_sequential();	
	
	// Copy the signals of the wrapped operator
	for(int i=0; i < op->getIOListSize(); i++)	{
		ioList.push_back( new Signal ( * op->getIOListSignal(i) ) );
	}
		
	// declare internal registered signals
	for(int i=0; i < op->getIOListSize(); i++){
		string idext = "i_" + (op->getIOListSignal(i))->id();
		//the clock is not registred
		if ( op->getIOListSignal(i)->id()!="clk")
			add_registered_signal(idext, op->getIOListSignal(i)->width());
	}

	//set pipeline parameters
	set_pipeline_depth(2 + op->pipeline_depth());
}

/**
 * The destructor
 **/
Wrapper::~Wrapper() {
}

/** Overloaded method which sets the operator name 
 *@param[in] prefix the prefix to be put in front of the default name
 *@param[in] postfix the postfix to be put at the end of the default name
 */
void Wrapper::set_operator_name(std::string prefix, std::string postfix){
	unique_name = prefix + "_" + op->getOperatorName() + "_" + postfix;
}

/**
 * Method belonging to the Operator class overloaded by the Wrapper class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void Wrapper::output_vhdl(ostream& o, string name) {

	Licence(o,"Florent de Dinechin (2007)");
	Operator::StdLibs(o);

	output_vhdl_entity(o);
	o << "architecture arch of " << name  << " is" << endl;

	// the operator to wrap
	op->Operator::output_vhdl_component(o);
	// The local signals
	output_vhdl_signal_declarations(o);

	o << "begin\n";
	o << "--wrapper operator"<<endl;
	// connect inputs
	for(int i=0; i<op->getIOListSize(); i++){
		string idext = "i_" + op->getIOListSignal(i)->id() ;
		if ((op->getIOListSignal(i)->type() == Signal::in) && (op->getIOListSignal(i)->id()!="clk"))
			o << tab << idext << " <=  " << op->getIOListSignal(i)->id() << ";" << endl;
	}

	// the instance
	o << tab << "test:" << op->getOperatorName() << "\n"
		<< tab << tab << "port map ( ";
	for(int i=0; i < op->getIOListSize(); i++) {
		Signal s = *op->getIOListSignal(i);
		if(i>0) 
			o << tab << tab << "           ";
		string idext = "i_" + op->getIOListSignal(i)->id() ;
		if (op->getIOListSignal(i)->type() == Signal::in)
			if (op->getIOListSignal(i)->id()!="clk")
				o << op->getIOListSignal(i)->id()  << " =>  " << idext << "_d";
			else
				o << op->getIOListSignal(i)->id()  << " =>  clk";
		else
			o << op->getIOListSignal(i)->id()  << " =>  " << idext;
		if (i < op->getIOListSize()-1) 
			o << "," << endl;
	}
	o << ");" <<endl;

	// the delays 
	output_vhdl_registers(o);

	// connect outputs
	for(int i=0; i<op->getIOListSize(); i++){
		string idext = "i_" + op->getIOListSignal(i)->id() ;
		if(op->getIOListSignal(i)->type() == Signal::out)
			o << tab << op->getIOListSignal(i)->id() << " <=  " << idext << "_d;" << endl;
	}
	
	o << "end architecture;" << endl << endl;
		
}

