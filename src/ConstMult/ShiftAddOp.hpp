#ifndef SHIFTADDOP_HPP
#define SHIFTADDOP_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

class ShiftAddDag;


typedef enum {X, Add, Shift, Neg} OpType;


/**
 Class ShiftAddOp defines a shift-and-add operation, of one of 3 types:
   Add(z,i, s, y)     Vz  <-   Vi<<s  + Vy 
   Shift(i,s)         Vz  <-   Vi<<s 
   Neg(i,s)        Vz  <-   (-Vi) 
   i, y and z are variable identifiers.
 This is single-assignment code, therefore 
   1/ the ShiftAddOp constructor builds a new variable number for the destination.
   2/ the ShiftAddOp object holds all the information related to its destination variable Vz. 
*/


class ShiftAddOp {
public:
	ShiftAddDag* impl;
	OpType op;
	ShiftAddOp* j;
	ShiftAddOp* i;

	/**  The shift on i*/
	int s; 

	/** the constant by which this variable multiplies */
	mpz_class n; 

	/** string representation of the constant */
	string name;  

	/** size of the constant */
	int size; 

	/** cost in term of full-adders */
	int cost_in_full_adders;

	// The following attributes are useful to pipeline a ShiftAddDag
	// Contrary to the previous, which are local, they only make sense
	// in the context of a ShiftAddDag. They will be set by the IntConstMult
	// constructor.

	/** is this operation registered?*/
	bool is_registered;

	/** does it need to be delayed (in a DAG where one branch has a
		 deeper pipeline than the other branch, the shallowest must be
		 delayed by the difference. 

		 delayed_by is the max of all the i_delayed_by or j_delayed_by for all the parents of this in the DAG

	*/
	int delayed_by;

	int i_delayed_by;

	int j_delayed_by;

	/** Constructor */
	ShiftAddOp(ShiftAddDag* impl, OpType op, ShiftAddOp* i=NULL, int s=0, ShiftAddOp* j=NULL);
	
	~ShiftAddOp(){};


	
	friend std::ostream& operator<<(std::ostream& o, const ShiftAddOp& sao ) // output
	{    
		o << sao.name << " <-  ";
		switch(sao.op) {
		case X:        o << " X"; break;
		case Add:      o << sao.i->name << "<<" << sao.s << "  + " << sao.j->name;   break;
		case Shift:    o << " " << sao.i->name << "<<" << sao.s;                     break;
		case Neg:      o << "-" << sao.i->name;   break;
		}   
		return o;
	}
};

#endif
