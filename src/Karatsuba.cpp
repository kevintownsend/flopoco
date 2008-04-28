/*
 * An integer multiplier for FloPoCo
 *
 * Authors : Bogdan Pasca
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
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"
#include "Karatsuba.hpp"

using namespace std;
extern vector<Operator*> oplist;

/** 
 * The constructor of the Karatsuba class
 * @param target argument of type Target containing the data for which this operator will be optimized
 * @param wInX integer argument representing the width in bits of the input X 
 * @param wInY integer argument representing the width in bits of the input Y
 **/
Karatsuba:: Karatsuba(Target* target, int wInX, int wInY) :
  Operator(target), wInX(wInX), wInY(wInY), wOut(wInX + wInY){
 
  int depth, i, j;
  ostringstream name, nameH, nameL;

  //gets information if the multiplier is or not to be pipelined
  if(target->is_pipelined())
    set_sequential();
  else
    set_combinatorial();  

  /* Name Setup procedure
   *  The name has the format: IntMultiplier_wInX_wInY
   *  wInX = width of the X input
   *  wInY = width of the Y input
   */  
  name.str("");;
  name <<"Karatsuba_"<<wInX<<"_"<<wInY;
  unique_name = name.str(); 
  
  /* Set up the IO signals
   * X and Y have wInX and wInY bits respectively 
   * R has wOut bits where wOut = (wInX + WInY) bits
   */
  add_input ("X", wInX);
  add_input ("Y", wInY);
  add_output("R", wOut);
  
  //set parameters
  int multiplierWidthDiscarded;
  bool status = target->suggest_submult_size(multiplierWidth, multiplierWidthDiscarded, wInX, wInY); 
  if (!status)
    cout<<"Warning: the requested frequency cannot be reached. "<<endl;
  
  cout<<"The suggested width for the multiplier is: "<<  multiplierWidth<<endl;

  if (is_sequential()){}
  else{
    
  }
  
  

}



/**
 * Karatsuba destructor
 */
Karatsuba::~Karatsuba() {
}





/**
 * Method belonging to the Operator class overloaded by the Karatsuba class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void Karatsuba::output_vhdl(std::ostream& o, std::string name) {
  ostringstream tempO;
  
  Licence(o,"Bogdan Pasca (2008)");
  Operator::StdLibs(o);  
  output_vhdl_entity(o);

 
  if (is_sequential()){}
  else
    BuildCombinationalKaratsuba(tempO,wInX-1,0,"X","Y","R", 0,  "mb");

  new_architecture(o,name);
  output_vhdl_signal_declarations(o);
  begin_architecture(o);
  o<<tempO.str();
  
  
  end_architecture(o);
}


void Karatsuba::BuildCombinationalKaratsuba(std::ostream& o, int L, int R, string lName, string rName, string resultName, int depth, string branch)
{
int termWidth;
int basisBits;
ostringstream t1LeftTerm, t1RightTerm, t1Result, t0Result, t2Result, t0Plust2, t1Minus_t0Plust2,t0ShiftPlust2;

//cout<<"Enter recursively with L="<<L<<" R="<<R<<" Left"<<lName<<" right="<<rName<<" ResultName="<<resultName<<" depth="<<depth<<" branch="<<branch<<endl;
// std::cin.get();


  if (L-R+1 <= multiplierWidth){
   
    o << resultName << " <= "<< lName <<"("<<L<<" downto "<<R<<") * "<<
                                rName <<"("<<L<<" downto "<<R<<");"<<endl;
  //cout<<"Solved with L="<<L<<"R="<<R<<endl;
  }
  else{ 
    /* parameters setup */
    termWidth = L - R + 1; 
    if (termWidth % 2 == 0)
      basisBits = termWidth / 2;
    else
      basisBits = termWidth / 2 + 1;
    
  //  cout<<"The width in bits of the basis: "<<basisBits<<endl;
           
    t1LeftTerm       <<"t1LeftTerm_"       <<depth<<"_"<<branch;
    t1RightTerm      <<"t1RightTerm_"      <<depth<<"_"<<branch;
    t1Result         <<"t1Result_"         <<depth<<"_"<<branch;
    t0Result         <<"t0Result_"         <<depth<<"_"<<branch;
    t2Result         <<"t2Result_"         <<depth<<"_"<<branch;
    t0Plust2         <<"t0Plust2_"         <<depth<<"_"<<branch;
    t1Minus_t0Plust2 <<"t1Minus_t0Plust2_" <<depth<<"_"<<branch;
    t0ShiftPlust2    <<"t0ShiftPlust2_"    <<depth<<"_"<<branch;
  
    add_signal( t1LeftTerm.str(), basisBits+1); 
    add_signal(t1RightTerm.str(), basisBits+1);
    o<<  t1LeftTerm.str() <<" <= (\"0"<<  ((termWidth % 2 == 0)?"":"0") <<"\" & "<<lName<<"("<< L            <<" downto "<< L - basisBits +((termWidth % 2 == 0)?1:2)<<")) + "
                 <<"(\"0\" & "<<lName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;    
    o<< t1RightTerm.str() <<" <= (\"0"<<  ((termWidth % 2 == 0)?"":"0") <<"\" & "<<rName<<"("<< L            <<" downto "<< L - basisBits +((termWidth % 2 == 0)?1:2)<<")) + "
                 <<"(\"0\" & "<<rName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;    
    
   //                           <<"(\"0\" & "                                      <<rName<<"("<< basisBits -1 <<" downto "<< 0                <<"));"<<endl;
        //                 <<"(\"0\" & "<<lName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;    
    
    add_signal(t1Result.str(),2*(basisBits+1));
    add_signal(t2Result.str(),2*basisBits);
    add_signal(t0Result.str(),2* ((termWidth % 2 == 0)?basisBits:(basisBits-1))  );
    
    BuildCombinationalKaratsuba(o, basisBits,   0            , t1LeftTerm.str(), t1RightTerm.str(), t1Result.str(), depth+1, branch + "_cb");
    //cout<<"Right Brench, depth = "<<depth <<endl;
 //   BuildCombinationalKaratsuba(o, basisBits-1, 0            , lName           , rName            , t2Result.str(), depth+1, branch + "_rb");
    BuildCombinationalKaratsuba(o, L-basisBits+((termWidth % 2 == 0)?1:2)-1, L-2*basisBits+((termWidth % 2 == 0)?1:2), lName , rName, t2Result.str(), depth+1, branch + "_rb");
    //cout<<"Left Brench, depth="<<depth<<endl;
    BuildCombinationalKaratsuba(o, L          , L-basisBits + ((termWidth % 2 == 0)?1:2) , lName           , rName            , t0Result.str(), depth+1, branch + "_lb");
    
    
    add_signal(t0Plust2.str()        ,2* basisBits+1 );
    add_signal(t1Minus_t0Plust2.str(),2*(basisBits+1));
    add_signal(t0ShiftPlust2.str()   ,2*(L-R+1)          ); 
     
    o<< t0Plust2.str()         <<" <= "<<"( \"0\" & "<<t0Result.str() << ") + ( \"0\" & "<< t2Result.str()<<");"<<endl;
    o<< t1Minus_t0Plust2.str() <<" <= "<< t1Result.str() << " - "<<"(\"0\" & "<< t0Plust2.str()<<")"<<";"<<endl;
    
    o<< t0ShiftPlust2.str()    <<" <= "<< "("<< t0Result.str()                            <<" & "<< zero_generator(2*basisBits,0)<<") + "
                                       << "("<< zero_generator( 2*(L-R+1) - 2*basisBits,0)<<" & "<< t2Result.str()               <<");"<<endl;
    o<< resultName             <<" <= "<< t0ShiftPlust2.str() << " + "
                                       << "("<<zero_generator(2*(L-R+1)-(3*basisBits+2),0)<<" & "<< t1Minus_t0Plust2.str()<<" & "<<zero_generator(basisBits,0)<<");"<<endl;   
  }
}                            




/**
 * A zero generator method which takes as input two arguments and returns a string of zeros with quotes as stated by the second argurment
 * @param[in] n        integer argument representing the number of zeros on the output string
 * @param[in] margins  integer argument determining the position of the quotes in the output string. The options are: -2= no quotes; -1=left quote; 0=both quotes 1=right quote
 * @return returns a string of zeros with the corresonding quotes given by margins
 **/
string Karatsuba::zero_generator(int n, int margins)
{
  ostringstream left,full, right, zeros;
  int i;

  //generate the zeros without quotes
  for (i=1; i<=n;i++)
    zeros<<"0";

  //generate the 3 strings with left quote, both and right quote
  left<<"\""<<zeros.str();
  full<<left.str()<<"\"";
  right<<zeros.str()<<"\"";

  //margins determines which of the strings are outputed
  switch(margins){
    case -2: return zeros.str(); break;
    case -1: return left.str();  break;
    case  0: return full.str();  break;
    case  1: return right.str(); break;
    default: return full.str();
  }
}

TestCaseList Karatsuba::generateRandomTestCases(int n)
{
	// TODO
	/* Signals */
	Signal sx = *get_signal_by_name("X");
	Signal sy = *get_signal_by_name("Y");
	Signal sr = *get_signal_by_name("R");

	TestCaseList tcl;	/* XXX: Just like Lyon's Transporation company. :D */
	mpz_class x, y, r;

	for (int i = 0; i < n; i++)	
	{
		x = getLargeRandom(sx.width());
		y = getLargeRandom(sy.width());
		r = x * y;

		TestCase tc;
		tc.addInput(sx, x);
		tc.addInput(sy, y);
		tc.addExpectedOutput(sr, r);
		tc.addComment(x.get_str() + " * " + y.get_str() + " = " + r.get_str());
		tcl.add(tc);
	}

	return tcl;
}

TestCaseList Karatsuba::generateStandardTestCases(int n)
{
	// TODO
	return TestCaseList();
}
