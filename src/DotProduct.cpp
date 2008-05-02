/*
 * Floating Point Multiplier for FloPoCo
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


/**
 * The DotProduct constructor
 * @param[in]		target	the target device
 * @param[in]		wE			the width of the exponent for the inputs X and Y
 * @param[in]		wFX     the width of the fraction for the input X
 * @param[in]		wFY     the width of the fraction for the input Y
 * @param[in]		MaxMSBX	maximum expected weight of the MSB of the summand
 * @param[in]		LSBA    The weight of the LSB of the accumulator; determines the final accuracy of the result
 * @param[in]		MSBA    The weight of the MSB of the accumulator; has to greater than that of the maximal expected result
**/ 
DotProduct::DotProduct(Target* target, int wE, int wFX, int wFY, int MaxMSBX, int LSBA, int MSBA ):
	Operator(target), wE(wE), wFX(wFX), wFY(wFY), MaxMSBX(MaxMSBX), LSBA(LSBA), MSBA(MSBA)  {
	ostringstream name;

	/* Set up the name of the entity */
	name <<"DotProduct_"<<wE<<"_"<<wFX<<"_"<<wFY<<"_"
			 <<(MaxMSBX>=0?"":"M")<<abs(MaxMSBX)<<"_"
			 <<(LSBA>=0?"":"M")<<abs(LSBA)<<"_"
			 <<(MSBA>=0?"":"M")<<abs(MSBA) ;
	unique_name=name.str();
	
	/* Set up the I/O signals of of the entity */
	add_input("X", 2 + 1 + wE + wFX);
	add_input("Y", 2 + 1 + wE + wFY);
	add_output("A", MSBA-LSBA+1); //the width of the output represents the accumulator size
	

  /* This operator is sequential.*/
	set_sequential();

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
  add_signal("fpMultiplierResultExponent",wE);
  add_signal("fpMultiplierResultSignificand", fpMultiplierResultFractionWidth + 1);
  add_signal("fpMultiplierResultException", 2);
  add_signal("fpMultiplierResultSign",1);
  add_signal("fpMultiplierResult",2 + 1 + wE + fpMultiplierResultFractionWidth);
   
  /* Set the pipeline depth of this operator */
  set_pipeline_depth(fpMultiplier->pipeline_depth() + longAcc->pipeline_depth() );
}



/**
 * DotProduct destructor
 */
DotProduct::~DotProduct() {
}



/**
 * Method belonging to the Operator class overloaded by the DotProduct class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void DotProduct::output_vhdl(std::ostream& o, std::string name) {
  Licence(o,"Bogdan Pasca (2008)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
  new_architecture(o, name);
  fpMultiplier->output_vhdl_component(o); //output the code of the fpMultiplier component
	longAcc->output_vhdl_component(o);      //output the code of the longAcc component
	output_vhdl_signal_declarations(o);
  begin_architecture(o);
  
  /* Map the signals on the fp multiplier */
	o << tab << "fp_multiplier: " << fpMultiplier->unique_name << endl;
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
  o << tab << "long_acc: " << longAcc->unique_name << endl;
	o << tab << "    port map ( X => fpMultiplierResult, " << endl;
	o << tab << "               A => A, " << endl;
	o << tab << "               clk => clk," << endl;
	o << tab << "               rst => rst " << endl;
	o << tab << "                         );" << endl; 
	o << endl;
  
	end_architecture(o);
} 


