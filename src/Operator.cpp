/*
 * The base Operator class, every operator should inherit it
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
#include <fstream>
#include <string>
#include <sstream>

#include "Operator.hpp"




void  Operator::add_input(const std::string name, const int width) {
  ioList.push_back(new Signal(name, Signal::in, width));
  number_of_inputs ++;
}

void  Operator::add_output(const std::string name, const int width) {
  ioList.push_back(new Signal(name, Signal::out, width));
  number_of_outputs ++;
}

void  Operator::add_FP_input(const std::string name, const int wE, const int wF) {
  ioList.push_back(new Signal(name, Signal::in, wE, wF));
  number_of_inputs ++;
}

void  Operator::add_FP_output(const std::string name, const int wE, const int wF) {
  ioList.push_back(new Signal(name, Signal::out, wE, wF));
  number_of_outputs ++;
}

void  Operator::add_signal(const std::string name, const int width) {
  signalList.push_back(new Signal(name, Signal::wire, width));
}

void  Operator::add_registered_signal(const std::string name, const int width) {
  signalList.push_back(new Signal(name, Signal::registered, width));
  ostringstream o;
  o << name <<"_d";
  signalList.push_back(new Signal(o.str(), Signal::wire, width));
  has_registers=true;
}

void  Operator::add_registered_signal_with_reset(const std::string name, const int width) {
  signalList.push_back(new Signal(name, Signal::registered_with_asynch_reset, width));
  ostringstream o;
  o << name <<"_d";
  signalList.push_back(new Signal(o.str(), Signal::wire, width));
  has_registers_with_asynch_reset=true;
}

void  Operator::add_registered_signal_with_synch_reset(const std::string name, const int width) {
  signalList.push_back(new Signal(name, Signal::registered_with_synch_reset, width));
  ostringstream o;
  o << name <<"_d";
  signalList.push_back(new Signal(o.str(), Signal::wire, width));
  has_registers_with_synch_reset=true;
}



string  Operator::add_delay_signal(const string name, const int width, const int delay) {
  ostringstream o;
  o << name;
  // if the delay is zero it is equivalent to add_signal
  if(delay>0) {
    for (int i=0; i<delay; i++){
      signalList.push_back(new Signal(o.str(), Signal::registered, width));    
      o  <<"_d";
    }
    has_registers=true;
  }
  signalList.push_back(new Signal(o.str(), Signal::wire, width));    
  
  return o.str();
}

string  Operator::get_delay_signal_name(const string name, const int delay) {
  ostringstream o;
  o << name;
  for (int i=0; i<delay; i++){
    o  <<"_d";
  }
  return o.str();
}


      
void  Operator::output_vhdl_signal_declarations(std::ostream& o) {
  for (int i=0; i<this->signalList.size(); i++){
    Signal* s = this->signalList[i];
    o<<tab<<  s->toVHDL() << ";" << endl;
  }
//   // Then the registered signals
//   for (int i=0; i<this->registerList.size(); i++){
//     Signal* s = this->registerList[i];
//     o<<tab<<  s->toVHDL() << ";" << endl;
//   }
//   // Then the registered signals with reset
//   for (int i=0; i<this->register_with_resetList.size(); i++){
//     Signal* s = this->register_with_resetList[i];
//     o<<tab<<  s->toVHDL() << ";" << endl;    
//   }
}



void  Operator::output_vhdl_registers(std::ostream& o) {


  // First registers without a reset
  if (has_registers) {
    o << tab << "process(clk)  begin\n"
      << tab << tab << "if clk'event and clk = '1' then\n";
    for(int i=0; i<signalList.size(); i++) {
      Operator::Signal *s = signalList[i];
      if(s->type()==Signal::registered) 
	o << tab <<tab << tab << s->id() <<"_d" << " <=  " << s->id() <<";\n";
    }
    o << tab << tab << "end if;\n";
    o << tab << "end process;\n"; 
  }
  
  // then registers with a reset
  if (has_registers_with_asynch_reset) {
    o << "  process(clk, rst)" << endl;
    o << "    begin" << endl;
    o << "      if rst = '1' then" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Operator::Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_asynch_reset)
         if (s->width()>1) 
	         o << tab <<tab << tab << s->id() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
         else
            o << tab <<tab << tab << s->id() <<"_d" << " <=  '0';\n";
    }
    o << "      elsif clk'event and clk = '1' then" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Operator::Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_asynch_reset) 
	o << tab <<tab << tab << s->id() <<"_d" << " <=  " << s->id() <<";\n";
    }
    o << "      end if;" << endl;
    o << "    end process;" << endl;
  }

  // then registers with synchronous reset
  if (has_registers_with_synch_reset) {
    o << "  process(clk, rst)" << endl;
    o << "    begin" << endl;
    o<<  "    if clk'event and clk = '1' then" << endl;
    o << "      if rst = '1' then" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Operator::Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_synch_reset)
         if (s->width()>1) 
	         o << tab <<tab << tab << s->id() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
         else
            o << tab <<tab << tab << s->id() <<"_d" << " <=  '0';\n";
    }
    o << "      else" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Operator::Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_synch_reset) 
	o << tab <<tab << tab << s->id() <<"_d" << " <=  " << s->id() <<";\n";
    }
    o << "      end if;" << endl;
    o << "    end if;" << endl;
    o << "    end process;" << endl;
  }
   


}

void Operator::output_vhdl_component(std::ostream& o, std::string name) {
  o << tab << "component " << name << " is" << endl
    << tab << "   port ( ";
  for (int i=0; i<this->ioList.size(); i++){
    Signal* s = this->ioList[i];
    if (i>0) // align signal names 
      o<<tab<<"          ";
    o<<  s->toVHDL();
    if(i < this->ioList.size()-1)  o<<";" << endl;
  }
  if (is_sequential()) {
    o <<";"<<endl<<"          clk: in std_logic";
    o <<";"<<endl<<"          rst: in std_logic";
  }
  o << "  );"<<endl;
  o << tab << "end component;" << endl;
}

void Operator::output_vhdl_component(std::ostream& o) {
  this->output_vhdl_component(o,  this->unique_name); 
}


void Operator::output_vhdl_entity(std::ostream& o) {
  o << "entity " << unique_name << " is" << endl
    << "   port ( ";
  for (int i=0; i<this->ioList.size(); i++){
    Signal* s = this->ioList[i];
    if (i>0) // align signal names 
      o<<"          ";
    o<<  s->toVHDL();
    if(i < this->ioList.size()-1)  o<<";" << endl;
  }
  if (is_sequential()) {
    o <<";"<<endl<<"          clk: in std_logic";
    o <<";"<<endl<<"          rst: in std_logic";
  }
  
  o << "  );"<<endl;
  o << "end entity;" << endl << endl;
}




void Operator::Licence(std::ostream& o, std::string authorsyears){
  o<<"--------------------------------------------------------------------------------"<<endl;
  // centering the unique name
  int s;
  if(unique_name.size()<76) s = (76-unique_name.size())/2; else s=0;
  o<<"--"; for(int i=0; i<s; i++) o<<" "; o <<unique_name <<endl; 
  o<<"-- This operator is part of the Infinite Virtual Library FloPoCoLib"<<endl
   <<"-- and is distributed under the terms of the GNU Lesser General Public Licence."<<endl
   <<"-- Authors: " << authorsyears <<endl
   <<"--------------------------------------------------------------------------------"<<endl;
}

void Operator::output_vhdl(std::ostream& o) {
  this->output_vhdl(o,  this->unique_name); 
}




void Operator::add_test_case(Operator::TestCase t){
  // Check that there are the right number of inputs
  if (t.input.size() != number_of_inputs) {
    cerr << "Error in add_test_case: expected "<<number_of_inputs << " inputs, got " << t.input.size() <<endl;
  }
  // Check that there are the right number of outputs
  if (t.expected_output.size() != number_of_outputs) {
    cerr << "Error in add_test_case: expected "<<number_of_outputs << " outputs, got " << t.expected_output.size() <<endl;
  }
  // Check that the sizes match
  // TODO
  // Insert into test case list
  // TODO

}


void Operator::output_vhdl_test(std::ostream& o) {
}



bool Operator::is_sequential() {
  return _is_sequential; 
}
void Operator::set_sequential() {
  _is_sequential=true; 
}
void Operator::set_combinatorial() {
  _is_sequential=false; 
}
int Operator::pipeline_depth() {
  return _pipeline_depth; 
}
void Operator::increment_pipeline_depth() {
  _pipeline_depth++; 
}
void Operator::set_pipeline_depth(int d) {
  _pipeline_depth=d; 
}

