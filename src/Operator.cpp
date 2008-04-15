/* vim: set tabstop=8 softtabstop=2 shiftwidth=2: */
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
  if(_signal_map.find(name) != _signal_map.end()) {
    std::ostringstream o;
    o << "ERROR in add_input, signal " << name<< " seems to already exist";
    throw o.str();
  }
  Signal *s = new Signal(name, Signal::in, width) ;
  ioList.push_back(s);
  _signal_map[name] = s ;
  number_of_inputs ++;
}

void  Operator::add_output(const std::string name, const int width) {
  if(_signal_map.find(name) != _signal_map.end()) {
    std::ostringstream o;
    o << "ERROR in add_input, signal " << name << " seems to already exist";
    throw o.str();
  }
  Signal *s = new Signal(name, Signal::out, width) ;
  ioList.push_back(s);
  _signal_map[name] = s ;
    number_of_outputs ++;
  }

void  Operator::add_FP_input(const std::string name, const int wE, const int wF) {
  if(_signal_map.find(name) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  Signal *s = new Signal(name, Signal::in, wE, wF) ;
  ioList.push_back(s);
  _signal_map[name] = s ;
  number_of_inputs ++;
}

void  Operator::add_FP_output(const std::string name, const int wE, const int wF) {
  if(_signal_map.find(name) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  Signal *s = new Signal(name, Signal::out, wE, wF) ;
  ioList.push_back(s);
  _signal_map[name] = s ;
  number_of_outputs ++;
}


void  Operator::add_signal(const std::string name, const int width) {
  if(_signal_map.find(name) != _signal_map.end()) {
    std::ostringstream o;
    o << "ERROR in add_signal, signal " << name << " seems to already exist";
    throw o.str();  
  }
  Signal *s = new Signal(name, Signal::wire, width);
  _signal_map[name] = s ;
  signalList.push_back(s);
}



void  Operator::add_registered_signal(const std::string name, const int width) {
  if(_signal_map.find(name) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  Signal *s;
  s = new Signal(name, Signal::registered, width);
  signalList.push_back(s);
  _signal_map[name] = s ;

  string o =  name + "_d";
  if(_signal_map.find(o) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  s = new Signal(o, Signal::wire, width);
  signalList.push_back(s);
  _signal_map[name] = s ;

  has_registers=true;
}



void  Operator::add_registered_signal_with_async_reset(const std::string name, const int width) {
  Signal *s;

  if(_signal_map.find(name) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  s = new Signal(name, Signal::registered_with_async_reset, width);
  signalList.push_back(s);
  _signal_map[name] = s ;

  string o = name + "_d";
  if(_signal_map.find(o) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  s = new Signal(o, Signal::wire, width);
  signalList.push_back(s);
  _signal_map[name] = s ;

  has_registers_with_async_reset=true;
}




void  Operator::add_registered_signal_with_sync_reset(const std::string name, const int width) {
  Signal *s;

  if(_signal_map.find(name) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  s = new Signal(name, Signal::registered_with_sync_reset, width);
  signalList.push_back(s);
  _signal_map[name] = s ;

  string o = name + "_d";
  if(_signal_map.find(o) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  s = new Signal(o, Signal::wire, width);
  signalList.push_back(s);
  _signal_map[name] = s ;

  has_registers_with_sync_reset=true;
}





string  Operator::add_delay_signal(const string name, const int width, const int delay) {
  ostringstream o;
  Signal *s;
  o << name;
  // if the delay is zero it is equivalent to add_signal
  if(delay>0) {
    for (int i=0; i<delay; i++){
      if(_signal_map.find(o.str()) != _signal_map.end()) {
	cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
	exit(EXIT_FAILURE);
      }
      s = new Signal(o.str(), Signal::registered, width);
      signalList.push_back(s);    
      _signal_map[name] = s ;
      o  <<"_d";
    }
    has_registers=true;
  }

  if(_signal_map.find(o.str()) != _signal_map.end()) {
    cerr << "ERROR in add_input , signal " << name<< " seems to already exist" << endl;
    exit(EXIT_FAILURE);
  }
  s = new Signal(o.str(), Signal::wire, width);
  signalList.push_back(s);    
  _signal_map[name] = s ;
  
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


Signal * Operator::get_signal_by_name(string name) {
  return _signal_map[name];
}


/*
string  Operator::add_level_delay_signal(const string name, const int start_level, const int width, const int delay) {
  ostringstream o;
  int l=start_level;
  //o << name;
  // if the delay is zero it is equivalent to add_signal
  if(delay>0) {
    for (int i=0; i<delay; i++){
      o.str("");
      o<<name<<"_level"<<l;
      signalList.push_back(new Signal(o.str(), Signal::registered_with_synch_reset, width));    
      l++;
    }
    has_registers_with_synch_reset=true;
  }
  signalList.push_back(new Signal(o.str(), Signal::wire, width));    
  
  return o.str();
}
*/

      
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
      Signal *s = signalList[i];
      if(s->type()==Signal::registered) 
	o << tab <<tab << tab << s->id() <<"_d" << " <=  " << s->id() <<";\n";
    }
    o << tab << tab << "end if;\n";
    o << tab << "end process;\n"; 
  }
  
  // then registers with a reset
  if (has_registers_with_async_reset) {
    o << tab << "process(clk, rst)" << endl;
    o << tab << tab << "begin" << endl;
    o << tab << tab << tab << "if rst = '1' then" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_async_reset)
         if (s->width()>1) 
	         o << tab <<tab << tab << s->id() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
         else
            o << tab <<tab << tab << s->id() <<"_d" << " <=  '0';\n";
    }
    o << tab << tab << tab << "elsif clk'event and clk = '1' then" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_async_reset) 
	o << tab <<tab << tab << s->id() <<"_d" << " <=  " << s->id() <<";\n";
    }
    o << tab << tab << tab << "end if;" << endl;
    o << tab << tab << "end process;" << endl;
  }

  // then registers with synchronous reset
  if (has_registers_with_sync_reset) {
    o << tab << "process(clk, rst)" << endl;
    o << tab << tab << "begin" << endl;
    o<<  "    if clk'event and clk = '1' then" << endl;
    o << tab << tab << tab << "if rst = '1' then" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_sync_reset)
         if (s->width()>1) 
	         o << tab <<tab << tab << s->id() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
         else
            o << tab <<tab << tab << s->id() <<"_d" << " <=  '0';\n";
    }
    o << tab << tab << tab << "else" << endl;
    for(int i=0; i<signalList.size(); i++) {
      Signal *s = signalList[i];
      if(s->type()==Signal::registered_with_sync_reset) 
	o << tab <<tab << tab << s->id() <<"_d" << " <=  " << s->id() <<";\n";
    }
    o << tab << tab << tab << "end if;" << endl;
    o << tab << tab << "end if;" << endl;
    o << tab << tab << "end process;" << endl;
  }
}

void Operator::output_vhdl_component(std::ostream& o, std::string name) {
  o << tab << "component " << name << " is" << endl;
  if (ioList.size() > 0)
  {
    o << tab << tab << "port ( ";
    for (int i=0; i<this->ioList.size(); i++){
      Signal* s = this->ioList[i];
      if (i>0) // align signal names 
	o<<tab<<"          ";
      o<<  s->toVHDL();
      if(i < this->ioList.size()-1)  o<<";" << endl;
    }
    o << tab << ");"<<endl;
  }
  o << tab << "end component;" << endl;
}

void Operator::output_vhdl_component(std::ostream& o) {
  this->output_vhdl_component(o,  this->unique_name); 
}


void Operator::output_vhdl_entity(std::ostream& o) {
  o << "entity " << unique_name << " is" << endl;
  if (ioList.size() > 0)
  {
    o << tab << "port ( ";

    for (int i=0; i<this->ioList.size(); i++){
      Signal* s = this->ioList[i];
      if (i>0) // align signal names 
	o<<"          ";
      o<<  s->toVHDL();
      if(i < this->ioList.size()-1)  o<<";" << endl;
    }
  
    o << tab << ");"<<endl;
  }
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


/**
 * A new line macro 
 * @param[in,out] o - the stream to which the new line will be added
 **/
void Operator::new_line(std::ostream& o)
{
o<<endl;
}

/**
 * A new architecture macro 
 * @param[in,out] o 	- the stream to which the new architecture line will be added
 * @param[in]     name	- the name of the entity corresponding to this architecture
 **/
void Operator::new_architecture(std::ostream& o, std::string name)
{
o << "architecture arch of " << name  << " is" << endl;
}


/**
 * A begin architecture macro 
 * @param[in,out] o 	- the stream to which the begin line will be added
 **/
void Operator::begin_architecture(std::ostream& o)
{
o << "begin" << endl;
}

/**
 * A end architecture macro 
 * @param[in,out] o 	- the stream to which the begin line will be added
 **/
void Operator::end_architecture(std::ostream& o)
{
o << "end architecture;" << endl << endl;
}


void Operator::output_vhdl(std::ostream& o) {
  this->output_vhdl(o,  this->unique_name); 
}

bool Operator::is_sequential() {
  return _is_sequential; 
}
void Operator::set_sequential() {
  _is_sequential=true; 
  add_input("clk");
  add_input("rst");
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

void Operator::output_final_report() {
  cout << "Entity " << unique_name <<":"<< endl;
  if(this->pipeline_depth()!=0)
    cout << tab << "Pipeline depth = " << pipeline_depth() << endl;
  else
    cout << tab << "Not pipelined"<< endl;
}

