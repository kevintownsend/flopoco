/*
 * Utility for converting a FP number into its binary representation,
 * for testing etc
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
#include <fstream>
#include <string>
#include <vector>
#include <mpfr.h>

#include"../utils.hpp"
using namespace std;


static void usage(char *name){
  cerr << "\nUsage: "<<name<<" wE wF x\n" ;
  cerr << "  x is only a double at the moment, sorry";
  exit (EXIT_FAILURE);
}




int check_strictly_positive(char* s, char* cmd) {
  int n=atoi(s);
  if (n<=0){
    cerr<<"ERROR: got "<<s<<", expected strictly positive number."<<endl;
    usage(cmd);
  }
  return n;
}

int main(int argc, char* argv[] )
{
  if(argc != 4) usage(argv[0]);
  int wE = check_strictly_positive(argv[1], argv[0]);
  int wF = check_strictly_positive(argv[2], argv[0]);
  double x = atof(argv[3]);

  int exponent = 0;
  uint64_t biased_exponent;
  int sign = (x<0?1:0);
  double significand = (x<0?-x:x);
  while(significand < 1.0) {
    significand *= 2.0;
    exponent --;
  }
  while(significand >= 2.0) {
    significand *= 0.5;
    exponent ++;
  }

  // add exponent bias
  biased_exponent = exponent + (1<<(wE-1))-1;

  //TODO check exponent is within range
  // TODO truncated, not rounded ! 
    cerr<<"Warning, input was truncated, not rounded\n";

  //exn bits
  if(x==0)
    cout << "00";
  else
    cout << "01";

  // sign bit
  cout << sign;

  // exponent
  printBinNum(cout, biased_exponent, wE);
//   for (int i=0; i<wE; i++) {
//     cout << (exponent -((exponent>>1)<<1));
//     exponent = exponent >>1;
//   }
    
  // significand
  significand -= 1.0;
  for (int i=0; i<wF; i++) {
    if (significand >= 0.5) {
      cout << "1";
      significand = (significand - 0.5)* 2.0;
    }
    else
      {
      cout << "0";
      significand *= 2.0;
    }
  }
  cout<<endl;

  if(significand!=0)
    cerr<<"Warning, input was truncated, not rounded\n";

  return 0;
}
