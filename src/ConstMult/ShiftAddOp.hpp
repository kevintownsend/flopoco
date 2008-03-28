#ifndef SHIFTADDOP_HPP
#define SHIFTADDOP_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

class ShiftAddDag;

/* The class IntConstMult has two nested classes:

 ShiftAddOp defines a shift-and-add operation.

 ShiftAddDag deals with the implementation of an IntConstMult as a
   vector of ShiftAddOp. It defines the intermediate variables with
   their bit sizes and provide methods for evaluating the cost of an
   implementation. Hopefully. Some day.

*/

typedef enum {X, Add, Shift, NegShift} OpType;

// The following class defines a shift-and-add operation, of one of 4 types:
//   Add(z,i, s, y)     Vz  <-   Vi<<s  + Vy 
//   Shift(i,s)         Vz  <-   Vi<<s 
//   NShift(i,s)        Vz  <-   (-Vi)<<s 
//   i, y and z are variable identifiers.
// This is single-assignment code, therefore 
//   1/ the ShiftAddOp constructor builds a new variable number for the destination.
//   2/ the ShiftAddOp object holds all the information related to its destination variable Vz. 



class ShiftAddOp {
public:
  ShiftAddDag* impl;
  OpType op;
  ShiftAddOp* j;
  ShiftAddOp* i;
  int s; // The shift is on i
  mpz_class n; //the constant by which this variable multiplies
  string name; // its string representation  
  int size; // the size of n
  bool isRegister; // if true, buffer the output in a pipeline implementation 
  ShiftAddOp(ShiftAddDag* impl, OpType op, ShiftAddOp* i=NULL, int s=0, ShiftAddOp* j=NULL);
  
  ~ShiftAddOp(){};
  
  friend std::ostream& operator<<(std::ostream& o, const ShiftAddOp& sao ) // output
  {    
    o << sao.name << " <-  ";
    switch(sao.op) {
    case X:        o << " X"; break;
    case Add:      o << sao.i->name << "<<" << sao.s << "  + " << sao.j->name;   break;
    case Shift:    o << " " << sao.i->name << "<<" << sao.s;                     break;
    case NegShift: o << "-" << sao.i->name << "<<" << sao.s;   break;
    }   
    return o;
  }
};

#endif
