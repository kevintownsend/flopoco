/*
 * A multiplier by an integer constant for FloPoCo
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
#include "../utils.hpp"
#include "../Operator.hpp"

#include "ShiftAddOp.hpp"
#include "ShiftAddDag.hpp"

#include "IntConstMult.hpp"
//#include "rigo.h"

using namespace std;


// TODO Pipeline. fuse NegShifts. Rewriting in general







IntConstMult::IntConstMult(Target* _target, int _xsize, mpz_class _n) :
  Operator(_target), xsize(_xsize), n(_n){
  ostringstream name; 
  name <<"IntConstMult_"<<xsize<<"_"<<n;
  unique_name=name.str();

  implementation = new ShiftAddDag(this);
  nsize = intlog2(n);
  rsize = nsize+xsize;

  add_input("X", xsize);
  add_output("R", rsize);


  bits = new int[nsize];
  BoothCode = new int[nsize+1];
  // Build the binary representation
  mpz_class nn = n;
  int l = 0;
  while(nn!=0) {
    bits[l] = (nn.get_ui())%2;
    l++;
    nn = nn>>1;
  }
  if(verbose) {
    cout<<" "; 
    for (int i=nsize-1; i>=0; i--)    cout << ((int) bits[i]);   
    cout << endl;
  }
  recodeBooth();
  if(verbose) printBoothCode();
  if(verbose) cout << "   Non-zero digits:" << nonZeroInBoothCode <<endl;


#if 0
  if(false && n==mpz_class("254876276031724631276054471292942"))
    {
      const int PX=0;
      cerr<<"Optimization by rigo.c"<< endl;//                    x    s    y
      	/*
	*/
      implementation->computeVarSizes(); 
      implementation->result = implementation->sao.size()-1;
    }
  else 
  if(n==mpz_class("1768559438007110"))
    {
      const int PX=0;
      cerr<<"Optimization by rigo.c"<< endl;//                   
      implementation->addOp( new ShiftAddOp(implementation, NegShift,   PX, 0) );       // 1  mx = -u0
      implementation->addOp( new ShiftAddOp(implementation, Add,   0, 19,  0) );        // 2  u3 = u0 << 19 + u0
      implementation->addOp( new ShiftAddOp(implementation, Shift,   2, 20) );          // 3  u103 = u3 << 20
      implementation->addOp( new ShiftAddOp(implementation, Add,   3, 4,   3) );        // 4  u203 = u103 << 4  + u103
      implementation->addOp( new ShiftAddOp(implementation, Add,   0, 14,  1) );        // 5  u7 = u0 << 14 + mx
      implementation->addOp( new ShiftAddOp(implementation, Add,   5, 6,  0) );         // 6  u6 = u7 << 6 + u0
      implementation->addOp( new ShiftAddOp(implementation, Add,   6, 10,  0) );        // 7  u5 = u6 << 10 + u0
      implementation->addOp( new ShiftAddOp(implementation, Shift, 7,  16    ));         // 8  u1 = u5 << 16
      implementation->addOp( new ShiftAddOp(implementation, Add,   8, 0,   4) ) ;       // 9  u101 = u1 + u203
      implementation->addOp( new ShiftAddOp(implementation, Add,   0, 21,  1) );        // 10 u107 = u0 << 21 + mx
      implementation->addOp( new ShiftAddOp(implementation, Add,   10, 18,  0) );       // 11 u106 = u107 << 18 + u0
      implementation->addOp( new ShiftAddOp(implementation, Add,   11, 4,   1) );       // 12 u105 = u106 << 4 + mx
      implementation->addOp( new ShiftAddOp(implementation, Add,   12, 5,   0) );       // 13 u2 = u105 << 5 + u0
      implementation->addOp( new ShiftAddOp(implementation, Shift, 13, 1) );            // 14 u102 = u2 << 1
      implementation->addOp( new ShiftAddOp(implementation, NegShift,   14, 0) );       // 15 mu102 = - u102
      implementation->addOp( new ShiftAddOp(implementation, Add,   14, 2,   15) );       // 16 u202 = u102 << 2  + mu102
      implementation->addOp( new ShiftAddOp(implementation, Add,   9, 0,   16) );        // R = u101 + u202

      	/*
0  u0 = x
1  mx = -u0
2  u3 = u0 << 19 + u0
3  u103 = u3 << 20
4  u203 = u103 << 4  + u103
5  u7 = u0 << 14 + mx
6  u6 = u7 << 6 + u0
7  u5 = u6 << 10 + u0
8  u1 = u5 << 16
9  u101 = u1 + u203
10 u107 = u0 << 21 + mx
11 u106 = u107 << 18 + u0
12 u105 = u106 << 4 + mx
13 u2 = u105 << 5 + u0
14 u102 = u2 << 1
15 mu102 = - u102
16 u202 = u102 << 2  + mu102
R = u101 + u202
	*/
      implementation->computeVarSizes(); 
      implementation->result = implementation->saolist.size()-1;
    }
  else 

#endif

    buildMultBoothTree();
  if(verbose) showShiftAddDag();

}

IntConstMult::~IntConstMult() {
  delete [ ] bits;
  delete [ ] BoothCode;
  delete implementation;
}


void IntConstMult::printBoothCode() {
  for (int i=nsize; i>=0; i--) {
    if(BoothCode[i]==0)       cout << "0"; 
    else if(BoothCode[i]==1)  cout << "+" ;
    else if(BoothCode[i]==-1) cout << "-" ;   
  }
  cout << endl;
}

void IntConstMult::recodeBooth() {
  int i;
  int *b, *c, *d;
  nonZeroInBoothCode = 0;

  b = new int[nsize+1];
  for (i=0; i<nsize; i++)
    b[i] = bits[i];
  b[nsize]=0;

  c = new int[nsize+1];
	
  c[0] = 0;
  for (i=0; i<nsize; i++) {
    if (b[i] + b[i+1] + c[i] >= 2)
      c[i+1] = 1;
    else 
      c[i+1] = 0;
    BoothCode[i] = b[i] + c[i] -2*c[i+1];
  }
  BoothCode[nsize] = c[nsize];

  // This method generates +0- for 3 instead of 0++.
  // As the -1 will cost us more, we rewrite that

  if (nsize>=1) {
    for (i=nsize-1; i>=0; i--) {
      if ((BoothCode[i]==-1) && (BoothCode[i+1]==0) && (BoothCode[i+2]==1)) {
	BoothCode[i]=1;
	BoothCode[i+1]=1;
	BoothCode[i+2]=0;
      }
    }
  }
  for (i=0; i<=nsize; i++) {
    if (BoothCode[i] != 0) 
      nonZeroInBoothCode ++;
  }
  delete [] c; delete [] b;
}





// The simple, sequential version: the DAG is a rake
void  IntConstMult::buildMultBooth(){
  int v,k,i;
  ShiftAddOp *z;
  i=0;

  // build the opposite of the input
  ShiftAddOp* MX = new ShiftAddOp(implementation, NegShift, implementation->PX , 0);

  while(0==BoothCode[i]) i++;
  // first op is a shift, possibly of 0
  if (1==BoothCode[i]) {
    if(i==0) // no need to add a variable
      z = implementation->PX;
    else {
      z = new ShiftAddOp(implementation, Shift,   implementation->PX, i);
    }
  }
  else {
    if(i==0) 
      z = MX;
    else {
      z = new ShiftAddOp(implementation, Shift,   MX, i);
    }
  }
  i++;
  
  for (k=1; k<nonZeroInBoothCode; k++) {
    while(0==BoothCode[i]) i++;
    if (1==BoothCode[i]) 
      z = new ShiftAddOp(implementation, Add,   implementation->PX, i,  z);
    else
      z = new ShiftAddOp(implementation, Add,   MX, i,  z);
    i++;
  }
  // Which variable number holds the result?
  implementation->result = z;
}






// The same as the previous, but builds a balanced tree.
void  IntConstMult::buildMultBoothTree(){
  int v,k,i,j,nk;
  ShiftAddOp *z;
  ShiftAddOp**  level;
  int*     shifts;
  
  i=0;
  while(0==BoothCode[i]) i++;
  int globalshift=i;
  if (nonZeroInBoothCode==1) { // a power of two
    if(i==0) // no need to add a variable
      z = implementation->PX;
    else {
      z = new ShiftAddOp(implementation, Shift,   implementation->PX, i);
    }
    implementation->result = z;
  }
  else { // at least two non-zero bits

    // build the opposite of the input
    ShiftAddOp* MX = new ShiftAddOp(implementation, NegShift, implementation->PX , 0);

    // fill an initial array with Xs and MXs
    level = new ShiftAddOp*[nonZeroInBoothCode];
    shifts = new int[nonZeroInBoothCode];
    for (j=0; j<nonZeroInBoothCode-1; j++) {
      if (1==BoothCode[i]) 
	level[j] = implementation->PX;
      else
	level[j] = MX; 
      shifts[j] = i - globalshift;
      i++;
      while(0==BoothCode[i]) i++;
    }
    level[j] = implementation->PX;
    shifts[j] = i-globalshift;

    k=nonZeroInBoothCode;
    while(k!=1) {
      nk=k>>1;
      for (j=0; j<nk; j++) {
	
	level[j] = implementation->provideShiftAddOp(Add, level[2*j+1], (shifts[2*j+1]-shifts[2*j]),  level[2*j]);
	shifts[j] = shifts[2*j];
      }
      if(nk<<1 != k) {
	level[j] = level[2*j];
	shifts[j] = shifts[2*j];
	k=nk+1;
      }
      else 
	k=nk;
    }
    if(globalshift==0)
      implementation->result = level[0];
    else
      implementation->result = implementation->provideShiftAddOp(Shift, level[0], globalshift);
    
//     cout << "(";
//     for (j=0; j<nonZeroInBoothCode; j++) {
//       if (level[j]==MX) cout <<"-"; else cout<<"+";
//       cout << "<<"<<shifts[j];
//     }
//     cout<<")<<"<<globalshift;

  }

  delete level;
  delete shifts;
  //exit(0);
  cerr << "    Number of adders: "<<implementation->saolist.size() << endl;
}



void IntConstMult::showShiftAddDag(){
  cout<<"ShiftAddDag"<<endl;
  for (int i=0; i<implementation->saolist.size(); i++) {
    cout << *(implementation->saolist[i]) <<endl;
  }
};









void IntConstMult::output_vhdl(std::ostream& o, std::string name) {
  Licence(o,"Florent de Dinechin (2007)");
  Operator::StdLibs(o);
  output_vhdl_entity(o);

  o << "architecture arch of " << name  << " is" << endl;
  implementation->output_vhdl(o);
  o << "end architecture;" << endl << endl;
    
}


void optimizeLefevre(const vector<mpz_class>& constants) {
  

};

TestCaseList IntConstMult::generateRandomTestCases(int n_)
{
	/* Signals */
	Signal sx = *get_signal_by_name("X");
	Signal sr = *get_signal_by_name("R");

	TestCaseList tcl;	/* XXX: Just like Lyon's Transporation company. :D */
	mpz_class x, r;

	for (int i = 0; i < n_; i++)	
	{
		x = getLargeRandom(sx.width());
		r = x * n;

		TestCase tc;
		tc.addInput(sx, x);
		tc.addExpectedOutput(sr, r);
		tc.addComment(x.get_str() + " * " + n.get_str() + " = " + r.get_str());
		tcl.add(tc);
	}

	return tcl;
}

TestCaseList IntConstMult::generateStandardTestCases(int n_)
{
	// TODO
	return TestCaseList();
}
