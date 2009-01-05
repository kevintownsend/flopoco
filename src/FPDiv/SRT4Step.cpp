/*
 * A step of SRT4 division for  FloPoCo
 *
 * Author : Jeremie Detrey, Florent de Dinechin
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
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "../utils.hpp"

#include "SRT4Step.hpp"

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


SRT4Step::SRT4Step(Target* target, int wF) :
	Operator(target), wF(wF) {

	int i, j;
	ostringstream name;

	name<<"SRT4Step_"<<wF; 
	uniqueName_ = name.str(); 

	setCombinatorial();
		
	addInput ("x", wF+3);
	addInput ("d", wF+1);
	addInput ("dtimes3", wF+3);
	addOutput ("q", 3);
	addOutput ("w", wF+3);
	


	// -------- Pipeline setup--------------------

	if(isSequential())	// TODO	
		setPipelineDepth(0);
	else
		setPipelineDepth(0);

	


	// Signals-------------------------------
	addSignal("sel", 5);
	addSignal("qi", 3);
	addSignal("qTimesD", wF+4);
	addSignal("x0", wF+4);
	addSignal("w0", wF+4); // Beware was wF+4 downto 1
}

SRT4Step::~SRT4Step() {
}


void SRT4Step::outputVHDL(std::ostream& o, std::string name) {
  
	ostringstream signame, synch1, synch2, xname,zeros, zeros1, zeros2, str1, str2;
	int i;

	licence(o,"Jeremie Detrey, Florent de Dinechin (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);	
			
	outputVHDLSignalDeclarations(o);	  
	beginArchitecture(o);
	outputVHDLRegisters(o); o<<endl;

	o << tab << "sel <= x(" << wF+2 << " downto " << wF-1 << ") & d(" << wF-1 << ");" << endl; 
  	o << tab << "with sel select" << endl;
   o << tab << "qi <= " << endl;
	o << tab << tab << "\"001\" when \"00010\" | \"00011\"," << endl;
	o << tab << tab << "\"010\" when \"00100\" | \"00101\" | \"00111\"," << endl;
	o << tab << tab << "\"011\" when \"00110\" | \"01000\" | \"01001\" | \"01010\" | \"01011\" | \"01101\" | \"01111\"," << endl;
	o << tab << tab << "\"101\" when \"11000\" | \"10110\" | \"10111\" | \"10100\" | \"10101\" | \"10011\" | \"10001\"," << endl;
	o << tab << tab << "\"110\" when \"11010\" | \"11011\" | \"11001\"," << endl;
	o << tab << tab << "\"111\" when \"11100\" | \"11101\"," << endl;
#if 0 
	// FPLibrary's version, leaves more freedom to the optimizer 
	// No noticeable difference in synthesis results, and it adds a lot of warnings to the simulation
	o << tab << tab << "\"000\" when \"00000\" | \"00001\" | \"11110\" | \"11111\"," << endl;
	o << tab << tab << "\"---\" when others;" << endl;
#else
	o << tab << tab << "\"000\" when others;" << endl;
#endif
	o << endl;
	o << tab << "with qi select" << endl;
   o << tab << tab << "qTimesD <= "<< endl ;
	o << tab << tab << tab << "\"000\" & d            when \"001\" | \"111\"," << endl;
	o << tab << tab << tab << "\"00\" & d & \"0\"     when \"010\" | \"110\"," << endl;
	o << tab << tab << tab << "\"0\" & dTimes3             when \"011\" | \"101\"," << endl;
#if 0 // FPLibrary's version, leaves more freedom to the optimizer
	o << tab << tab << tab << "(" << wF+3 << " downto 0 => '0') when \"000\"," << endl;
	o << tab << tab << tab << "(" << wF+3 << " downto 0 => '-') when others;" << endl;
#else
	o << tab << tab << tab << "(" << wF+3 << " downto 0 => '0') when others;" << endl;
#endif
	o << endl;
	o << tab << "x0 <= x & \"0\";" << endl;
	o << tab << "with qi(2) select" << endl;
   o << tab << "w0 <= x0 - qTimesD when '0'," << endl;
	o << tab << "      x0 + qTimesD when others;" << endl;
	o << endl;
	o << tab << "q <= qi;" << endl;
  	o << tab << "w <= w0(" << wF+1 << " downto 0) & \"0\";" << endl;
	endArchitecture(o);
}

