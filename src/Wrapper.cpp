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

Wrapper::Wrapper(Target* target, Operator *op):
	Operator(target), op_(op)
{
	/* the name of the Wrapped operator consists of the name of the operator to be 
	wrapped followd by _Wrapper */
	setOperatorName();
		
	if (!target->isPipelined()) 	
		setSequential();	
	
	// Copy the signals of the wrapped operator
	for(int i=0; i < op->getIOListSize(); i++)	{
		ioList_.push_back( new Signal ( * op->getIOListSignal(i) ) );
	}
		
	// declare internal registered signals
	for(int i=0; i < op->getIOListSize(); i++){
		string idext = "i_" + (op->getIOListSignal(i))->getSignalName();
		//the clock is not registred
		if ( op->getIOListSignal(i)->getSignalName()!="clk")
			addDelaySignal(idext, op->getIOListSignal(i)->width());
	}

	//set pipeline parameters
	setPipelineDepth(2 + op->getPipelineDepth());
}

Wrapper::~Wrapper() {
}

void Wrapper::setOperatorName(){
	uniqueName_ = op_->getOperatorName() + "_Wrapper";
}

void Wrapper::outputVHDL(ostream& o, string name) {

	licence(o,"Florent de Dinechin (2007)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	op_->Operator::outputVHDLComponent(o);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	o << "--wrapper operator"<<endl;
	// connect inputs
	for(int i=0; i < op_->getIOListSize(); i++){
		string idext = "i_" + op_->getIOListSignal(i)->getSignalName() ;
		if ((op_->getIOListSignal(i)->type() == Signal::in) && (op_->getIOListSignal(i)->getSignalName()!="clk"))
			o << tab << idext << " <=  " << op_->getIOListSignal(i)->getSignalName() << ";" << endl;
	}

	// the instance
	o << tab << "test:" << op_->getOperatorName() << "\n"
		<< tab << tab << "port map ( ";
	for(int i=0; i < op_->getIOListSize(); i++) {
		Signal s = *op_->getIOListSignal(i);
		if(i>0) 
			o << tab << tab << "           ";
		string idext = "i_" + op_->getIOListSignal(i)->getSignalName() ;
		if (op_->getIOListSignal(i)->type() == Signal::in)
			if (op_->getIOListSignal(i)->getSignalName()!="clk")
				o << op_->getIOListSignal(i)->getSignalName()  << " =>  " << idext << "_d";
			else
				o << op_->getIOListSignal(i)->getSignalName()  << " =>  clk";
		else
			o << op_->getIOListSignal(i)->getSignalName()  << " =>  " << idext;
		if (i < op_->getIOListSize()-1) 
			o << "," << endl;
	}
	o << ");" <<endl;

	// the delays 
	outputVHDLRegisters(o);

	// connect outputs
	for(int i=0; i<op_->getIOListSize(); i++){
		string idext = "i_" + op_->getIOListSignal(i)->getSignalName() ;
		if(op_->getIOListSignal(i)->type() == Signal::out)
			o << tab << op_->getIOListSignal(i)->getSignalName() << " <=  " << idext << "_d;" << endl;
	}
	
	endArchitecture(o);
}

