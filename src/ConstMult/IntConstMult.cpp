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
extern vector<Operator*> oplist;


void reset_visited(ShiftAddOp* sao) {
	if (sao!=NULL) {
		sao->already_visited=false;
		switch(sao->op) {
		case X:
			break;
		case Add:
		case Sub:
		case RSub:
			reset_visited(sao->i);
			reset_visited(sao->j);
			break;
		case Shift:
		case Neg:
			reset_visited(sao->i);
			break;
		}
	}
}


// No longer used

int compute_tree_depth(ShiftAddOp* sao) {
	int ipd, jpd;
	if (sao==NULL)
		return 0;
	else {
		switch(sao->op) {
		case X:
			return 0;
		case Add:
		case Sub:
		case RSub:
			ipd = compute_tree_depth(sao->i);
			jpd = compute_tree_depth(sao->j);
			return 1+max(ipd, jpd);
		case Shift:
			ipd = compute_tree_depth(sao->i);
			return ipd;
		case Neg:
			ipd = compute_tree_depth(sao->i);
			return 1+ipd;
		}
		return 0;
	}
}

// Do not forget to call reset_visited before calling this one.
int compute_total_cost_rec(ShiftAddOp* sao) {
	if (sao==NULL || sao->already_visited)
		return 0;
	else {
		sao->already_visited=true;
		switch(sao->op) {
		case X:
			return 0;
		case Add:
		case Sub:
		case RSub:
			return sao->cost_in_full_adders + compute_total_cost_rec(sao->i) + compute_total_cost_rec(sao->j);
		case Shift:
		case Neg:
			return sao->cost_in_full_adders + compute_total_cost_rec(sao->i);
		}
	}
	throw string("This exception in IntConstMult::compute_total_cost_rec should never happen");
}

int compute_total_cost(ShiftAddOp* sao) {
	reset_visited(sao);
	return compute_total_cost_rec(sao);
}



/**
 * Depth-first traversal of the DAG to build the pipeline.
 * @param partial_delay accumulates the delays of several stages

 Starting from the leaves, we accumulate partial delays until target_period is reached.
 Then pipeline level will be inserted.
 
 */

void IntConstMult::build_pipeline(ShiftAddOp* sao, double &partial_delay) {
	string iname, jname, isignal, jsignal;
	double idelay,jdelay, max_children_delay;
	int size, isize, jsize, shift, adder_size; 
	bool use_pipelined_adder;
	IntAdder* adder;
	
	if (sao==NULL)
		return;
	else {

		// First check that the sub-DAG has not been already visited

		bool already_visited=true;
		try { 
			getSignalByName(sao->name);
		} catch (std::string s) {
			already_visited=false;
		}
		if(already_visited)
			return;

		// A few variables to make the code below easier to read
		ShiftAddOpType op = sao->op;
		size = sao->size;
		shift = sao->s;
			
		switch(op) {
		case X:
			partial_delay=0;
			setCycle(0, false);
			return;

		case Add:
		case Sub:
		case RSub:
			// A few variables to make the code below easier to read
			isize = sao->i->size;
			jsize = sao->j->size;

			build_pipeline(sao->i, idelay);
			if(sao->i != sao->j) {
				build_pipeline(sao->j, jdelay);
			}
			iname = sao->i->name; 
			jname = sao->j->name; 

			adder_size = sao->cost_in_full_adders+1;
			vhdl << endl << tab << "-- " << *sao <<endl; // comment what we're doing
			setCycleFromSignal(iname, false);
			syncCycleFromSignal(jname, false);

			max_children_delay = max(idelay,jdelay);

			// Now decide what kind of adder we will use, and compute the remaining delay
			use_pipelined_adder=false;
			if (isSequential()) {
				// First case: using a plain adder fits within the current pipeline level
				double tentative_delay = max_children_delay + target_->adderDelay(adder_size) + target_->localWireDelay();
				if(tentative_delay <= 1./target_->frequency()) {
					use_pipelined_adder=false;
					partial_delay = tentative_delay;					
				}
				else { 
					// register the children 
					nextCycle();
					// Is a standard adder OK ?
					tentative_delay = target_->ffDelay() + target_->localWireDelay() + target_->adderDelay(adder_size);
					if(tentative_delay <= 1./target_->frequency()) {
						use_pipelined_adder=false;
						partial_delay = tentative_delay;					
					}
					else { // Need to instantiate an IntAdder
						use_pipelined_adder=true;
						adder = new IntAdder(target_, adder_size);
						adder->changeName(getName() + "_" + sao->name + "_adder");
						oplist.push_back(adder);

						partial_delay = target_->adderDelay(adder->getLastChunkSize());
					}
				}
			}
			// Now generate VHDL

			if(shift==0) { // Add with no shift -- this shouldn't happen with current DAGs so te following code is mostly untested
				if(op==Sub || op==RSub)
					throw string("In IntConstMult::build_pipeline, Sub and RSub with zero shift currently unimplemented"); // TODO
				isignal = sao->name + "_High_L";  
				jsignal = sao->name + "_High_R"; 

				// The i part
				vhdl << tab << declare(isignal, size) << " <= ";
				if(size>isize+1) // need to sign-extend x
					vhdl <<"(" << size-1 << " downto " << isize <<" => '" << (sao->i->n >= 0 ? "0" : "1" ) << "') & ";
				vhdl << use(iname) << ";" << endl;
				// the y part
				vhdl << tab << declare(jsignal, size) << " <= ";
				if(size>jsize) // need to sign-extend y
					vhdl << "(" << size-1 <<" downto " << jsize <<" => '" << (sao->j->n >= 0 ? "0" : "1" ) << "') & ";
				vhdl << use(jname) << ";" << endl;

				if(use_pipelined_adder) { // Need to use an IntAdder subcomponent
					inPortMap  (adder, "X", isignal);
					inPortMap  (adder, "Y", jsignal);
					inPortMapCst  (adder, "Cin", "'0'");
					outPortMap (adder, "R",sao->name);
					vhdl << instance(adder, sao->name + "_adder");
				}
				else
					vhdl << tab << declare(sao->name, size) << " <= " << use(isignal) << " + " << use(jsignal) << ";" << endl;
			}


			else { // Add with actual shift
				if(op == Add || op==RSub) {
					if(shift >= jsize) { // Simpler case when the two words to add are disjoint; size=isize+s+1
						//                        jjjjjj
						//             +/-  iiii
						// TODO perf: use an IntAdder here when needed
						// The lower bits of the sum are those of y, possibly sign-extended but otherwise untouched
						vhdl << tab << declare(sao->name, sao->size) << "("<< shift - 1 <<" downto 0) <= " ;
						if(shift>jsize) {
							vhdl << "(" <<  shift-1 <<" downto " << jsize << " => ";
							if(sao->j->n >= 0)  vhdl << "'0'"; // pad with 0s 
							else                vhdl << use(jname) << "(" << jsize-1 << ")";// sign extend
							vhdl << ") & ";
						}
						vhdl << use(jname) << ";   -- lower bits untouched"<<endl;

						if(op == Add) {
							// The higher bits (size-1 downto s) of the result are those of x, possibly plus 11...1 if y was negative
							vhdl << tab << sao->name << "("<<sao->size-1<<" downto "<< shift<<") <= " << use(iname) ;
							if(sao->j->n < 0) { 
								vhdl <<" + (" << sao->size-1 <<" downto " <<  shift <<" => " << use(jname) << "(" << jsize-1 << ")) "
									  << ";   -- sum of higher bits"<<endl;
							}
							else 
								vhdl << ";   -- higher bits also untouched"<<endl;
						}
						else {// op == RSub
							// The higher bits (size-1 downto s) of the result are those of -x, possibly plus 11...1 if y was negative
							vhdl << tab << sao->name << "("<<sao->size-1<<" downto "<< shift<<") <= " ;
							if(sao->j->n < 0) 
								vhdl <<"(" << sao->size-1 << " downto " <<  shift <<" => " << use(jname) << "(" << jsize-1 << ")) ";
							else
								vhdl <<"(" << sao->size-1 << " downto " <<  shift <<" => " << "'0') ";
							vhdl << " - " << use(iname) << ";   -- sum of higher bits"<<endl;
						}
					} // end if (shift >= jsize)
					else{ 
						// jsize>s.        Cases:      xxxxx              xxxxxx
						//                                yyyyyyyyyy             yyyyyyyyyyyy
						// so we may need to sign-extend Vx, or Vy, or even both.
						// The higher bits of the result are sum/diff
						isignal = sao->name + "_High_L";  
						jsignal = sao->name + "_High_R"; 
						// The x part
						vhdl << tab << declare(isignal,  size - shift) << " <= ";
						if(size >= isize +  shift +1) { // need to sign-extend vx. If the constant is positive, padding with 0s is enough
							vhdl <<" (" << size-1 << " downto " << isize +  shift <<" => ";
							if(sao->i->n >= 0)   vhdl << "'0'";// pad with 0s 
							else                 vhdl << use(iname) << "(" << isize-1 << ")"; // sign extend
							vhdl << ") & ";
						}
						vhdl << use(iname) << "("<< isize -1 <<" downto 0) ;" << endl;
						// the y part
						vhdl << tab << declare(jsignal,  size - shift) << " <= ";
						if(size >= jsize+1) {// need to sign-extend vy. If the constant is positive padding with 0s is enough
							vhdl <<" (" << size-1 << " downto " << jsize <<" => ";
							if(sao->j->n >= 0)  vhdl << "'0'"; // pad with 0s 
							else                vhdl << use(jname) << "(" << jsize-1 << ")";// sign extend
							vhdl << ") & ";
						}
						vhdl << use(jname) << "("<< jsize -1 <<" downto " <<  shift << "); " << endl;
					
						// do the sum
						if(use_pipelined_adder) {
							inPortMap  (adder, "X", jsignal);
							if(op==Add) {
								inPortMap  (adder, "Y", isignal);
								inPortMapCst  (adder, "Cin", "'0'");
							} 
							else { // RSub
								string isignalneg = isignal+"_neg"; 
								vhdl << declare(isignalneg, size - shift) << " <= not " << use(isignal);
								inPortMap  (adder, "Y", isignal);
								inPortMapCst  (adder, "Cin", "'1'");
							}
							string resname=sao->name+"_h";
							outPortMap (adder, "R",resname);
							vhdl << instance(adder, sao->name + "_adder");

							syncCycleFromSignal(resname, false);
							//nextCycle();
							vhdl << tab << declare(sao->name, sao->size) << "("<<size-1<<" downto " <<  shift << ") <= " << use(resname) + ";" << endl;
						}
						else
							vhdl << tab << declare(sao->name, sao->size) << "("<<size-1<<" downto " <<  shift << ") <= " // vz (size-1 downto s)
								  << use(jsignal) << (op==Add ? " + " : "-") << use(isignal) << ";   -- sum of higher bits" << endl; 
			
						// In both cases the lower bits of the result (s-1 downto 0) are untouched
						vhdl << tab << sao->name << "("<<shift-1<<" downto 0) <= " << use(jname) <<"("<< shift-1<<" downto 0);   -- lower bits untouched"<<endl;

					} // end if (shift >= jsize) else
				} // end if(op == Add || op == RSub) 
				else { // op=Sub 
					// Do a normal subtraction of size size
					isignal = sao->name + "_L";  
					jsignal = sao->name + "_R"; 
					vhdl << tab << declare(isignal,  size) << " <= ";
					if(size > isize+shift) {// need to sign-extend vx. If the constant is positive padding with 0s is enough
						vhdl <<" (" << size-1 << " downto " << isize+shift <<" => ";
						if(sao->i->n >= 0)   vhdl << "'0'";// pad with 0s 
						else                 vhdl << use(iname) << "(" << isize-1 << ")";// sign extend
						vhdl << ") & ";
					}
					vhdl << use(iname) << " & (" << shift-1 << " downto 0 => '0')" << ";" << endl;

					vhdl << tab << declare(jsignal,  size) << " <= ";
					vhdl <<" (" << size-1 << " downto " << jsize <<" => ";
					if(sao->j->n >= 0)   vhdl << "'0'";// pad with 0s 
					else                 vhdl << use(jname) << "(" << jsize-1 << ")";// sign extend
					vhdl << ") & ";
					vhdl << use(jname) << ";" << endl;
					
					// do the subtraction
					if(use_pipelined_adder) {
						string jsignalneg = jsignal+"_neg"; 
						vhdl << declare(jsignalneg, size) << " <= not " << use(jsignal);
						inPortMap  (adder, "X", isignal);
						inPortMap  (adder, "Y", jsignalneg);
						inPortMapCst  (adder, "Cin", "'1'");
						string resname=sao->name+"_h";
						outPortMap (adder, "R",resname);
						vhdl << instance(adder, sao->name + "_adder");

						syncCycleFromSignal(resname, false);
						//nextCycle();
						vhdl << tab << declare(sao->name, size) << " <=  " << use(resname) + ";" << endl;
						
					}
				}
			}
			
			return;

			// shift and neg almost identical
		case Shift:
		case Neg:
			isize = sao->i->size;

			double local_delay;
			if(op == Neg){   
				local_delay = target_->adderDelay(sao->cost_in_full_adders);
			}
			else 
				local_delay=0;

			build_pipeline(sao->i, idelay);

			iname = sao->i->name; 
			setCycleFromSignal(iname, false);

			if(isSequential() 
				&& idelay +  target_->localWireDelay() + local_delay > 1./target_->frequency()
				&& sao->i->op != X) {
				// This resets the partial delay to that of this ShiftAddOp
				nextCycle();
				partial_delay =  target_->ffDelay() + target_->adderDelay(sao->cost_in_full_adders);
			}
			else{ // this ShiftAddOp and its child will be in the same pipeline level
				partial_delay = idelay + target_->localWireDelay() + local_delay;
			}
			vhdl << tab << declare(sao->name, sao->size) << " <= " ;
			// TODO use a pipelined IntAdder when necessary
			if(op == Neg)   
				vhdl << "("<< sao->size -1 <<" downto 0 => '0') - " << use(iname) <<";"<<endl; 
			else { // Shift
				if (shift == 0) 
					vhdl << use(iname) <<";"<<endl; 
				else
					vhdl << use(iname) <<" & ("<< shift - 1 <<" downto 0 => '0');"<<endl;
			}
			break;

		}   
	}
}



IntConstMult::IntConstMult(Target* _target, int _xsize, mpz_class _n) :
	Operator(_target), n(_n), xsize(_xsize){
	ostringstream name; 

	setCopyrightString("Florent de Dinechin (2007)");

	//C++ wrapper for GMP does not work properly on win32, using mpz2string
	name <<"IntConstMult_"<<xsize<<"_"<<mpz2string(n);
	uniqueName_=name.str();

	implementation = new ShiftAddDag(this);
	nsize = intlog2(n);
	rsize = intlog2(n * ((mpz_class(1)<<xsize)-1));

	addInput("inX", xsize);
	addOutput("R", rsize);
	
	cout << "  Power of two" <<intlog2(n) << " " <<  (mpz_class(1) << intlog2(n)) << endl;

	if((mpz_class(1) << (intlog2(n)-1)) == n) { // power of two
		if(verbose) 
			cout << "  Power of two" << endl;
		vhdl << tab << "R <= inX & " << rangeAssign(intlog2(n)-2, 0, "'0'") << ";" << endl;
	}
	else {
		vhdl << tab << declare("X", xsize) << " <= inX;" << endl; // no longer useful but harmless
		bits = new int[nsize];
		BoothCode = new int[nsize+1];
		// Build the binary representation -- I'm sure there is a mpz method for it
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



		// Build in implementation a tree constant multiplier 
		buildMultBoothTree();
		if(verbose) showShiftAddDag();
		int cost=compute_total_cost(implementation->result);
		if(verbose) {
			cout << "Estimated bare cost (not counting pipeline overhead) : " << cost << " FA/LUT" << endl;
		}

		
		double delay=0.0;
		// recursively build the pipeline in the vhdl stream
		build_pipeline(implementation->result, delay);
		
		// copy the top of the DAG into variable R
		vhdl << endl << tab << "R <= " << implementation->result->name << "("<< rsize-1 <<" downto 0);"<<endl;
		
	}
}




	// One hand-coded Lefevre multiplier, for comparison purposes -- to be inserted somewhere in the constructor

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
			implementation->addOp( new ShiftAddOp(implementation, Neg,   PX) );       // 1  mx = -u0
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
			implementation->addOp( new ShiftAddOp(implementation, Neg,   14) );       // 15 mu102 = - u102
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
	//	else
#endif


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
	int *b, *c;
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
	int k,i;
	ShiftAddOp *z;
	i=0;

	// build the opposite of the input
	ShiftAddOp* MX = new ShiftAddOp(implementation, Neg, implementation->PX);

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
	int k,i,j,nk;
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
		ShiftAddOp* MX = new ShiftAddOp(implementation, Neg, implementation->PX);

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

	}

	delete level;
	delete shifts;

	if(verbose) cerr << "   Number of adders: "<<implementation->saolist.size() << endl;
}



void IntConstMult::showShiftAddDag(){
	cout<<"    ShiftAddDag"<<endl;
	for (uint32_t i=0; i<implementation->saolist.size(); i++) {
		cout << "     "<<*(implementation->saolist[i]) <<endl;
	}
};








void optimizeLefevre(const vector<mpz_class>& constants) {
	

};


void IntConstMult::emulate(TestCase *tc){
	mpz_class svX = tc->getInputValue("inX");
	mpz_class svR = svX * n;
	tc->addExpectedOutput("R", svR);
}



void IntConstMult::buildStandardTestCases(TestCaseList* tcl){
	TestCase *tc;

	tc = new TestCase(this); 
	tc->addInput("inX", mpz_class(0));
	emulate(tc);
	tc->addComment("Multiplication by 0");
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("inX", mpz_class(1));
	emulate(tc);
	tc->addComment("Multiplication by 1");
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("inX", mpz_class(2));
	emulate(tc);
	tc->addComment("Multiplication by 2");
	tcl->add(tc);

	tc = new TestCase(this); 
	tc->addInput("inX", (mpz_class(1) << xsize) -1);
	emulate(tc);
	tc->addComment("Multiplication by the max positive value");
	tcl->add(tc);

//	tc = new TestCase(this); 
//	tc->addInput("inX", (mpz_class(1) << (xsize) -1) + mpz_class(1));
//	emulate(tc);
//	tc->addComment("Multiplication by 10...01");
//	tcl->add(tc);


}

