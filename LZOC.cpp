/*
 * A leading zero/one counter for FloPoCo
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

#include "LZOC.hpp"

using namespace std;








LZOC::LZOC(Target* target, int wIn, int wOut) :
  Operator(target), wIn(wIn), wOut(wOut)  {


  p2wOut = 1<<wOut; // no need for GMP here
  ostringstream name; 
  name <<"LZOC_"<<wIn<<"_"<<wOut;

  unique_name=name.str();

  // Set up the IO signals
  add_input ("I", wIn);
  add_input ("OZB");  
  add_output("O", wOut);
 }

LZOC::~LZOC() {}








void LZOC::output_vhdl(std::ostream& o, std::string name) {
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);
  output_vhdl_entity(o);
  o << "architecture arch of " << name  << " is" << endl;
  for (int i=wOut; i>0; i--) 
    o << "  signal level"<<i<<": std_logic_vector("<<(1<<i)-1<<" downto 0);" <<endl;
  
  o << "begin" << endl;
  // connect first stage to X
  if (p2wOut==wIn)
    o << "  level"<<wOut<<" <=  X ;" << endl;
  else if (p2wOut>wIn) // pad input with zeroes
    o << "  level"<<wOut<<" <=  X & ("<<p2wOut-wIn-1<<" downto 0 => '0') ;" << endl;
  else if (p2wOut<wIn)
    o << "  level"<<wOut<<" <=  X("<<wIn-1<<" downto "<<wIn - p2wOut <<") ;" << endl;
  // recursive structure
  for (int i=wOut; i>1; i--) {
    o << "  O("<<i-1<<") <= '1' when level"<<i<<" = ("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<" => ozb) else '0';"<< endl;
    o << "  level"<<i-1<<" <= level"<<i<<"("<<(1<<(i-1))-1<<" downto 0) when O("<<i-1<<") = '1'"<< endl
      << "               else level"<<i<<"("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<");" << endl;
  }
  o << "  O(0) <= '1' when level1(1) =  ozb else '0';"<< endl;
  o << "end architecture;" << endl << endl;
}


