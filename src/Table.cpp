/*
 * A generic class for tables of values
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
#include "utils.hpp"
#include "Table.hpp"

using namespace std;




int Table::double2input(double x){
  cerr << "Error, double2input is being used and has not been overriden";
  return 1;
}

double Table::input2double(int x) {
  cerr << "Error, input2double is being used and has not been overriden";
  return 1;
}

mpz_class Table::double2output(double x){
  cerr << "Error, double2output is being used and has not been overriden";
  return 0;
}

double Table::output2double(mpz_class x) {
  cerr << "Error, output2double is being used and has not been overriden";
  return 1;
}







void Table::outputComponent(ostream& o, string name)
{
  o << "  component " << name << " is" << endl
    << "    port ( x : in  std_logic_vector(" << wIn-1  << " downto 0);" << endl
    << "           y : out std_logic_vector(" << wOut - 1 << " downto 0) );" << endl
    << "  end component;" << endl << endl;

}

void Table::outputComponent(ostream& o, string name, int component_id)
{
  ostringstream n;
  n << name << "_" << component_id;
  outputComponent(o,  n.str()); 
}








void Table::output(ostream& o, string name, int component_id)
{
  ostringstream n;
  n << name << "_" << component_id;
  output(o,  n.str()); 
}


void Table::output(ostream& o, string name)
{
  if (wIn > 10)
    cerr << "Avertissement : production d'une table avec " << wIn << " bits en entree" << endl;

  o<<"library ieee;\nuse ieee.std_logic_1164.all;\nuse ieee.std_logic_arith.all;"<<endl
   <<"use ieee.std_logic_unsigned.all;\nlibrary work;\nuse work.pkg_fp_log.all;"<<endl<<endl;
  o << "entity  " << name << " is" << endl
    << "    port ( x : in  std_logic_vector(" << wIn-1  << " downto 0);" << endl
    << "           y : out std_logic_vector(" << wOut - 1 << " downto 0) );" << endl
    << "end entity;" << endl << endl;

  // retrait du code (sauf pour la première ligne) */
  const int spaces = 9;
  int i,x;
  mpz_class y;
  char* margin = new char[spaces + 1];

  o << "architecture arch of " << name  << " is" << endl
    << "begin" << endl
    << "  with x select" << endl
    << "    y <= ";

  // mise en retrait du code
  for (i = 0; i < spaces; i++) margin[i] = ' ';
  margin[i] = '\0';

  for (x = minIn; x <= maxIn; x++) {
    y=this->function(x);
    //    cout << x <<"  "<< y << endl;
    o << "\"";
      printBinPosNumGMP(o, y, wOut);
    o << "\" when \"";
    printBinNum(o, x, wIn);
    o << "\"," << endl
      << margin;
  }
  o << "\"";
  for (i = 0; i < wOut; i++) o << "-";
  o << "\" when others;" << endl;
  delete [] margin;
  o << "end architecture;" << endl << endl;
}


int Table::size_in_LUTs() {
  return wOut*int(intpow2(wIn-4));
}
