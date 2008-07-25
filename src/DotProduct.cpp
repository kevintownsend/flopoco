/*
 * Dot Product unit for FloPoCo
 *
 * Author : Bogdan Pasca
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
#include "utils.hpp"
#include "Operator.hpp"

#include "DotProduct.hpp"

using namespace std;
extern vector<Operator*> oplist;

DotProduct::DotProduct(Target* target, int wE, int wFX, int wFY, int MaxMSBX, int LSBA, int MSBA ):
	Operator(target), wE(wE), wFX(wFX), wFY(wFY), MaxMSBX(MaxMSBX), LSBA(LSBA), MSBA(MSBA)  {
	
	setOperatorName();
	
	/* Set up the I/O signals of of the entity */
	addInput("X", 2 + 1 + wE + wFX);
	addInput("Y", 2 + 1 + wE + wFY);
	addOutput("A", MSBA-LSBA+1); //the width of the output represents the accumulator size


  /* This operator is sequential.*/
	setSequential();

  /* Instantiate one FPMultiplier used to multiply the two inputs fp numbers, X and Y */
  int fpMultiplierResultExponentWidth = wE;
  int fpMultiplierResultFractionWidth = wFX + wFY + 1;
  int fpMultiplierResultNormalization = 1;
	fpMultiplier = new FPMultiplier(target, wE, wFX, wE, wFY, fpMultiplierResultExponentWidth, fpMultiplierResultFractionWidth, fpMultiplierResultNormalization);
	oplist.push_back(fpMultiplier);
  
  /* Instantiate one LongAcc used accumulate the products produced by fpMultiplier */
	longAcc = new LongAcc(target, wE, fpMultiplierResultFractionWidth, MaxMSBX, LSBA, MSBA);
	oplist.push_back(longAcc);
  
  /* Signal declaration */
  addSignal("fpMultiplierResultExponent",wE);
  addSignal("fpMultiplierResultSignificand", fpMultiplierResultFractionWidth + 1);
  addSignal("fpMultiplierResultException", 2);
  addSignal("fpMultiplierResultSign",1);
  addSignal("fpMultiplierResult",2 + 1 + wE + fpMultiplierResultFractionWidth);
   
  /* Set the pipeline depth of this operator */
  setPipelineDepth(fpMultiplier->getPipelineDepth() + longAcc->getPipelineDepth() );
}

DotProduct::~DotProduct() {
}

void DotProduct::setOperatorName(){
	ostringstream name;
	/* Set up the name of the entity */
	name <<"DotProduct_"<<wE<<"_"<<wFX<<"_"<<wFY<<"_"
			 <<(MaxMSBX>=0?"":"M")<<abs(MaxMSBX)<<"_"
			 <<(LSBA>=0?"":"M")<<abs(LSBA)<<"_"
			 <<(MSBA>=0?"":"M")<<abs(MSBA) ;
	uniqueName_=name.str();
}

void DotProduct::outputVHDL(std::ostream& o, std::string name) {
	licence(o,"Bogdan Pasca (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o, name);
	fpMultiplier->outputVHDLComponent(o); //output the code of the fpMultiplier component
	longAcc->outputVHDLComponent(o);      //output the code of the longAcc component
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);

	/* Map the signals on the fp multiplier */
	o << tab << "fp_multiplier: " << fpMultiplier->getOperatorName() << endl;
	o << tab << "    port map ( X => X, " << endl;
	o << tab << "               Y => Y, " << endl;
	o << tab << "               ResultExponent => fpMultiplierResultExponent, " << endl;
	o << tab << "               ResultSignificand => fpMultiplierResultSignificand, " << endl;
	o << tab << "               ResultException => fpMultiplierResultException, " << endl;
	o << tab << "               ResultSign => fpMultiplierResultSign, " << endl;
	o << tab << "               clk => clk," << endl;
	o << tab << "               rst => rst " << endl;
	o << tab << "                         );" << endl; 
	o << endl;

	/*rebuild the output under the form of a fp number */
	o << tab << "fpMultiplierResult <= fpMultiplierResultException & "
	                                <<"fpMultiplierResultSign & "
	                                <<"fpMultiplierResultExponent & "
	                                <<"fpMultiplierResultSignificand("<<wFX + wFY<<" downto "<< 0 <<");"<<endl;

	/* Map the signals on the long accumulator */
	o << tab << "long_acc: " << longAcc->getOperatorName() << endl;
	o << tab << "    port map ( X => fpMultiplierResult, " << endl;
	o << tab << "               A => A, " << endl;
	o << tab << "               clk => clk," << endl;
	o << tab << "               rst => rst " << endl;
	o << tab << "                         );" << endl; 
	o << endl;

	endArchitecture(o);
} 


