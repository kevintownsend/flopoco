/*
 * Shift-and-add DAG for integer constant multiplication
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

using namespace std;

mpz_class ShiftAddDag::computeConstant(OpType op, ShiftAddOp* i, int s, ShiftAddOp* j) {
  switch(op) {
  case X:
    if(i!=NULL || j!=NULL) {
      cerr << "ERROR Unexpected non-null pointer in computeConstant(X,..). \n";
      exit(EXIT_FAILURE);
    }
    return 1;

  case Add: 
    if(i==NULL || j==NULL) {
      cerr << "ERROR Unexpected null pointer in computeConstant(Add, ...), exiting. \n";
      exit(EXIT_FAILURE);
    }
    // the constant by which this variable multiplies 
    return (i->n << s) +  j->n;

  case Shift:
    if(i==NULL || j!=NULL) {
      cerr << "ERROR In inputs to computeConstant(Shift,...), exiting. \n";
      exit(EXIT_FAILURE);
    }
    return i->n << s ;

  case NegShift:
    if(i==NULL || j!=NULL) {
      cerr << "ERROR In ShiftAddOp, cannot construct such a NegShift, exiting. \n";
      exit(EXIT_FAILURE);
    }
    return - i->n<<s;
  }

}



// This method looks up in the current Dag if the requireed op
// exists, and either returns a pointer to it, or creates the
// corresponding node.
ShiftAddOp* ShiftAddDag::provideShiftAddOp(OpType op, ShiftAddOp* i, int s, ShiftAddOp* j){
  mpz_class n=this->computeConstant(op, i, s, j);
  for(int ii=0; ii<this->saolist.size(); ii++) {
    if (n==saolist[ii]->n) 
      return saolist[ii];
  }
  return new ShiftAddOp(this, op, i, s, j);
}



void ShiftAddDag::output_vhdl(std::ostream& o) {
  int k,i,j;
  // Signal declarations
  for (i=0; i<saolist.size(); i++) 
    o << tab << "signal " <<saolist[i]->name<<" : std_logic_vector("<<saolist[i]->size-1<<" downto 0);"<<endl;
  
  // Architecture
  o << "begin" << endl;
  for (i=0; i<saolist.size(); i++) {
    ShiftAddOp *p = saolist[i];
    o<<tab<<"-- " << *p <<endl;

    switch(p->op) {
    case X: 
      cerr << "ERROR unexpected ShiftAddOp(X) in output_vhdl()"; exit(EXIT_FAILURE);
      break;

    case Add:
      if(p->s==0) {
	o << tab << p->name << " <= " ;
	// The i part
	if(p->size>p->i->size+1) // need to sign-extend x
	  o <<"( (" << p->size-1 << " downto " << p->i->size <<" => '" << (p->i->n >= 0 ? "0" : "1" ) << "') & " << p->i->name << ")";
	else 
	  o << p->i->name;
	o<<" + " ;
	// the y part
	if(p->size>p->j->size+1) // need to sign-extend y
	  o << "( (" << p->size-1 <<" downto " << p->j->size <<" => '" << (p->j->n >= 0 ? "0" : "1" ) << "') & " << p->j->name << ") ;";
	else 
	  o << p->j->name << ";";
      }

      else { // Add with actual shift
	if(p->s >= p->j->size) { // Simpler case when the two words to add are disjoint; size=p->i->size+s+1
	  //                   xxxxx
	  //                           yyyyy
	  // The lower bits of the sum are those of y, possibly sign-extended but otherwise untouched
	  o << tab << p->name << "("<< p->s - 1 <<" downto 0) <= "
	    << "(" <<  p->s-1 <<" downto " << p->j->size <<" => " << p->j->name << "(" << p->j->size-1 << ")) & "    // sign extend y
	  << p->j->name << ";   -- lower bits untouched"<<endl;
	  // The higher bits (size-1 downto s) of the result are those of x, possibly plus 11...1 if y was negative
	  o << tab << p->name << "("<<p->size-1<<" downto "<< p->s<<") <= " << p->i->name ;
	  if(p->j->n < 0) 
	    o<<" + (" << p->size-1 <<" downto " <<  p->s <<" => " << p->j->name << "(" << p->j->size-1 << ")) ";
	  o<< ";   -- sum of higher bits"<<endl;     
	}
	else{ // p->j->size>s.        Cases:      xxxxx              xxxxxx
	  //                                yyyyyyyyyy             yyyyyyyyyyyy
	  // so we may need to sign-extend Vx, or Vy, or even both.
	  // In both cases the lower bits of the result (s-1 downto 0) are untouched
	  o << tab << p->name << "("<<p->s-1<<" downto 0) <= " << p->j->name <<"("<< p->s-1<<" downto 0);   -- lower bits untouched"<<endl;
	  // The higher bits of the result are sum/diff
	  o << tab << p->name << "("<<p->size-1<<" downto " <<  p->s << ") <= "; // vz (size-1 downto s) <=
	  // The x part
	  o <<"(";
	  if(p->size >= p->i->size +  p->s +1) {// need to sign-extend vx. If the constant is positive padding with 0s is enough
	    o <<" (" << p->size-1 << " downto " << p->i->size +  p->s <<" => ";
	    if(p->i->n >= 0) // pad with 0s 
	      o<< "'0'";
	    else // sign extend
	      o<< p->i->name << "(" << p->i->size-1 << ")";
	    o << ") & ";
	  }
	  o << p->i->name << "("<< p->i->size -1 <<" downto 0)";
	  o << ")   +   (" ;
	  // the y part
	  if(p->size>=p->j->size+1) {// need to sign-extend vy. If the constant is positive padding with 0s is enough
	    o<<" (" << p->size-1 << " downto " << p->j->size <<" => ";
	    if(p->j->n >= 0) // pad with 0s 
	      o<< "'0'";
	    else // sign extend
	      o << p->j->name << "(" << p->j->size-1 << ")";
	    o << ") & ";
	  }
	  o << p->j->name << "("<< p->j->size -1 <<" downto " <<  p->s << ") ); " <<endl;
	}
      }
    
    break;

    case Shift:
    case NegShift:
      o << tab << p->name << " <= " ;
      if(p->op == NegShift)   
	o << "("<< p->size -1 <<" downto 0 => '0') - " ; 
      if (p->s == 0) 
	o << p->i->name <<";"<<endl; 
      else
	o  << p->i->name <<" & ("<< p->s - 1 <<" downto 0 => '0');"<<endl; 
      break;
    }     
  }
  // Sometimes the size of the result variable is one bit more than r.
  o << tab << "r <= " << result->name << "("<< icm->rsize-1 <<" downto 0);"<<endl;
}





