/*
 * A barrel shifter for FloPoCo
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

#include "Shifters.hpp"

using namespace std;





  // TODO there is a small inefficiency here, as most bits of s don't need to be copied all the way down





Shifter::Shifter(Target* target, int wIn, int maxShift, ShiftDirection direction) :
  Operator(target), wIn(wIn), maxShift(maxShift), direction(direction) {

  wOut = wIn + maxShift;
  wShiftIn = intlog2(maxShift);
  ostringstream name;
  if(direction==Left)
    name <<"LeftShifter_";
  else
    name <<"RightShifter_";

  name<<wIn<<"_by_max_"<<maxShift;

  unique_name=name.str();

  // Set up the IO signals
  add_input ("X", wIn);
  add_input ("S", wShiftIn);  
  add_output("R", wOut);

 if(target->is_pipelined()) 
   set_sequential();
 else
   set_combinatorial();
 
  // evaluate the pipeline

  double critical_path = 0.0;
  for (int i=0; i<wShiftIn; i++) 
    level_registered[i] = false;

  if(is_sequential()) {
    for (int i=0; i<wShiftIn; i++) {
      /* approximate delay of this stage */
      double stage_delay = target->lut_delay() + target->distant_wire_delay((1<<i));
      if (critical_path + stage_delay > 1/target->frequency()) {
	critical_path=stage_delay; 	// reset critical path
	level_registered[i+1]= true;
	increment_pipeline_depth();
      }
      else{ 
	critical_path += stage_delay; // augment critical path
	level_registered[i+1] = false;
      }
    }
    // register the last level
    level_registered[wShiftIn] = true;
    increment_pipeline_depth();
  }

  // Set up the intermediate signals 
  add_signal("level0", wIn);
  //  o << "  signal level0: std_logic_vector("<< wIn <<"-1 downto 0);" <<endl;
  for (int i=1; i<=wShiftIn; i++) {
    ostringstream sname;
    sname << "level"<<i;
    if (level_registered[i])
      add_registered_signal(sname.str(), wIn + (1<<i) -1 );
    else
      add_signal(sname.str(), wIn + (1<<i) -1);
    //o << "  signal level"<<i<<": std_logic_vector("<< wIn <<"+"<<intpow2(i)-1<<"-1 downto 0);" <<endl;
  }

  // The shift input has to be delayed as well now
  if(pipeline_depth()>=1) 
    add_delay_signal("ps", wShiftIn, pipeline_depth()-1); 
  else
    add_signal("ps", wShiftIn);

}


Shifter::~Shifter() {
}







void Shifter::output_vhdl(std::ostream& o, std::string name) {
  ostringstream signame;
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);
  output_vhdl_entity(o);

  o << "architecture arch of " << name  << " is" << endl;
  
  output_vhdl_signal_declarations(o);

  o << "begin" << endl;
  int stage=0;
  o << "   level0 <=  X ;" << endl;
  o << "   ps <=  s;" <<endl;
  ostringstream psname;
  psname << "ps";
  for (int i=0; i<wShiftIn; i++) {
    ostringstream lname;
    lname << "level"<<i;
    if (level_registered[i]) { // use the registered signal instead
      lname << "_d";
      // and use next stage of ps
      psname << "_d",
      // add a synchronisation barrier here
      o <<"  ----- synchro barrier ------- " <<endl;
      stage++;
    }
			
    o << "   level"<<i+1<<" <=  ";
    o << "("<<intpow2(i)-1<<" downto 0 => '0') & "<<lname.str()
      <<"  when "<<psname.str()<<"("<<i<<") = '"<<(direction==Right?1:0)<<"'   else  ";
    o << lname.str()<<" & ("<<intpow2(i)-1<<" downto 0 => '0');" << endl;
    
  }
  if(level_registered[wShiftIn])
    o <<"  ----- synchro barrier ------- " <<endl;
 
 o << "   R <=  level"<<wShiftIn;
  if(level_registered[wShiftIn]) 
    o << "_d";
  o << "("<< wOut-1<<" downto 0);" << endl << endl;


  if(is_sequential())
    output_vhdl_registers(o);
  o << "end architecture;" << endl << endl;
}



