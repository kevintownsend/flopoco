/*
 * The abstract class that models different chips for delay and area. 

 * Classes for real chips inherit from this one. They should be in subdirectory Targets
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

#ifndef TARGET_HPP
#define TARGET_HPP



class Target
{
 public:
    
  virtual double lut_delay() =0;
  virtual double carry_propagate_delay() =0;
  virtual double adder_delay(int n) =0;
  virtual double local_wire_delay() =0;
  virtual double distant_wire_delay(int n) =0;



  void set_pipelined();
  void set_not_pipelined();
  bool is_pipelined();
  void set_lut_inputs(int n);
  int lut_inputs();
  int mult_x_inputs(); // may be 18 
  int mult_y_inputs(); // may be 24 
  double frequency();
  void set_frequency(double f);

  Target()   {
    _pipeline=true;
    _lut_inputs=4;
    _frequency = 400000000.;
  }

  virtual ~Target() {}

protected:
  int _lut_inputs;
  //std::string name;
  bool _pipeline;
  double _frequency; // Hz
  int _mult_x_inputs; 
  int _mult_y_inputs;
};


#endif
