#include "ConstMult/FixRealKCM.hpp"
#include "IntMultiplier.hpp"
#include "FixedPointFunctions/FunctionTable.hpp"
#include "IntConstDiv.hpp"

// works only with sollya
#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers */
#include "gmp.h"
#include "mpfr.h"
#include <assert.h>
// for debug
#include <signal.h>
#include "FixSinCos.hpp"

using namespace std;
using namespace flopoco;


#define SUBCYCLEPIPELINING 0
#define USEBITHEAP 1

FixSinCos::FixSinCos(Target * target, int w_):Operator(target), w(w_)
{
	srcFileName="FixSinCos";
	
	// definition of the name of the operator
	ostringstream name;
	name << "FixSinCos_" << w ;
	setName(name.str());

	setCopyrightString("Florent de Dinechin, Guillaume Sergent (2012)");

	// declaring inputs
	addInput("X", w+1);

	// declaring outputs
	addOutput("S", w+1);
	addOutput("C", w+1);

	/*********************************** RANGE REDUCTION **************************************/

	// the argument is reduced into (0,1/4) because one sin/cos
	// computation in this range can always compute the right sin/cos
	// 2's complement: 0 is always positive
	// to handle 2's complement efficiently, we consider values of X
	// in [-1,0[ as they were in [1,2[, just adding 2 and treating
	// the binary representation of X as if it was an unsigned fixed-point
	// in [0,2[
	vhdl << tab << declare ("X_sgn") << " <= X" << of (w) << ";" << endl;
	vhdl << tab << declare ("Q") << " <= X" << of (w-1) << ";" << endl;
	vhdl << tab << declare ("O") << " <= X" << of (w-2) << ";" << endl;
	vhdl << tab << declare ("Y",w-2) << " <= X " << range (w-3,0) << ";" << endl;
	// now X -> X_sgn + Q*.5 + O*.25 + Y where Q,O \in {0,1} and Y \in {0,.25}

	// Y_prime = .25 - y
	// perhaps do a logic ~ at a cost of 1 ulp ?
	//vhdl << declare ("Y_prime",w-2) << " <= "
	//     << "2**" << (w-2) << " - Y;" << endl;
	// if we do an arithmetic 1's complement it will do **** since
	// 1-0 doesn't fit
	manageCriticalPath(target->localWireDelay(w-2) + target->lutDelay());
	vhdl << tab << declare ("Y_prime", w-2) << " <= " << "not Y;" << endl;

	// we need to know the number of guard bits _now_ to have a good Y_in
	// some precision is lost with the optimized multipliers
	const int g=6;
	//const int g=4; //guard bits
	//const int g=3; //guard bits
	//we take 4 guard bits even if error < 8 ulp because rounding will take
	//another half final ulp
	int wIn = w-2+g;

	// Y_in = / Y_prime if O=1
	//        \ Y if O=0
	// and we extend the precision to make it look as if Y_prime was
	// 0.[1, wIn times] - Y: 1/2**g error reduction in Y_prime
	// (which should be arithmetic 1-Y)
	vhdl << tab << declare ("Y_in",wIn) << " <= Y_prime & "
	     << '"' << std::string (g, '1') << '"' << " when O='1' else Y & "
	     << '"' << std::string (g, '0') << '"' << ";"
	     << endl;

	// Exch = Q ^ O
	vhdl << tab << declare ("Exch") << " <= Q xor O;" << endl;

	// now we do manual polynomial computation to calculate sin (pi*y)
	// and cos (pi*y) using Taylor polynomials: with y'=pi*y
	// sin (4*pi*y) = y' - y'³/6
	// cos (4*pi*y) = 1 - y'²/2
	// this works if y' (or y) is small enough
	// to accomplish this we decompose x (==Y_in in the vhdl) to x = a + y
	// where a \in [0,1/4[ and {sin,cos} (pi*a) is tabulated
	// and y \in [0,1b-n[ is a small enough argument 
	// then we use the addition formulae (where Sin(x)=sin(pi*x)):
	// Sin (a+y) = Sin a Cos y + Cos a Sin y
	//           = Sin a - Sin a (1 - Cos y) + Cos a Sin y
	// Cos (a+y) = Cos a Cos y - Sin a Sin y
	//           = Cos a - Cos a (1 - Cos y) - Sin a Sin y
	// wA is the precision of A, beginning from the ,..[] bit
	// yes I know it's moronic to int->float->int
	// pi⁴/24=4.0587121 so (pi*y)⁴/24 is < (1+eps) ulp /this is not the thing to do actually/
	// iff 4 * (wA+2) - 2 >= w+g (2 as the log_2 of pi⁴/24)
	// the minimal a is therefore ceil((wIn)/4) - 2
	int wA = (int) ceil ((double) (w+g+2)/4.) - 2;
	if (wA <= 3)
		wA = 3;



	/*********************************** THE TABLES **************************************/

	// upper log2 of memory width needed for table
	int wg_log2 = (int) ceil ((double) log2 (w+g));
	int bram_words = target->sizeOfMemoryBlock() / (1u << wg_log2);
	int bram_words_log2 = (int) floor ((double) log2 (bram_words));
	// FIXME: assumes that both
	// -the memory width is a power of 2
	// -the bram can always be used as dual-port (one for sin, one for cos)
	//  with width >=w+g
	if (wA <= bram_words_log2-1)
		wA = bram_words_log2-1;
	int wY = wIn-wA, wZ = wY+2;
	// vhdl:split (Y_in -> A & Y_red)
	vhdl << tab << declare ("A",wA) << " <= Y_in" << range (wIn-1,wIn-wA) << ";" << endl;
	vhdl << tab << declare ("Y_red",wY) << " <= Y_in" << range (wIn-wA-1,0) << ';' << endl;
	// vhdl:lut (A -> A_cos_pi_tbl, SinPiA)
	FunctionTable *sin_table, *cos_table;
	{
		ostringstream omu; // one minus (guardless) ulp
		omu << "(1 - 1b-" << w << ")";
		// we'll include the half-ulp for final rounding
		// directly in the tables
		// TODO: a better check on wA to ensure there isn't too much
		// error using the with-rounding {Sin,Cos}PiA in multipliers
		ostringstream halfulp;
		halfulp << "1b-" << (w+1);
		sin_table = new FunctionTable (target, omu.str()
						     + " * sin(pi*x/4) + "
						     + halfulp.str()
					     , wA, -(w+g), -1);
		cos_table = new FunctionTable (target, omu.str()
						     + " * cos(pi*x/4) + "
						     + halfulp.str()
					     , wA, -(w+g), -1);
	}
	sin_table -> changeName(getName() + "_SinTable");
	cos_table -> changeName(getName() + "_CosTable");
	oplist.push_back (sin_table);
	oplist.push_back (cos_table);
	// SinPiA and CosPiA already contain the rounding constant
	outPortMap (sin_table, "Y", "SinPiA");
	inPortMap (sin_table, "X", "A");
	outPortMap (cos_table, "Y", "CosPiA");
	inPortMap (cos_table, "X", "A");
	vhdl << instance (sin_table, "sin_table");
	vhdl << instance (cos_table, "cos_table");

	// the results have precision w+g


	/*********************************** THE MULTIPLIER BY PI **************************************/

	// now, evaluate Sin Y_red and 1 - Cos Y_red
	// vhdl:cmul[pi] (Y_red -> Z)
	map<string, double> pi_mult_inputDelays;
	pi_mult_inputDelays["X"] = getCriticalPath();


#if 1 
	// For w=32: latency=18 (1 for the const mult), 781 + 1170	
	FixRealKCM *pi_mult;
	// the 2 extra bits to Z are added by the KCM multiplier
	pi_mult = new FixRealKCM (target, 0, wY-1, false, 0, "pi", 1.0, pi_mult_inputDelays); 
	oplist.push_back (pi_mult);
	inPortMap (pi_mult, "X", "Y_red");
	outPortMap (pi_mult, "R", "Z");
	vhdl << instance (pi_mult, "pi_mult");

	syncCycleFromSignal("Z", pi_mult->getOutputDelay("R"));
#else
	// For w=32: latency=17 (4 for the const mult), 863+1334
	// So this one seems to loose
	// TODO design a FixReal shift-and-add
	mpfr_t myPi;
	mpfr_init2 (myPi, wZ);
	mpfr_const_pi (myPi, GMP_RNDN);
	mpfr_mul_2si(myPi, myPi, wZ-2, GMP_RNDN);
	mpz_class mpzPi;
	mpfr_get_z(mpzPi.get_mpz_t(), myPi, GMP_RNDN);

	IntConstMult* pi_mult;
	pi_mult = new IntConstMult(target, wZ-2, mpzPi); 
	oplist.push_back (pi_mult);
	inPortMap (pi_mult, "X", "Y_red");
	outPortMap (pi_mult, "R", "Zfull");
	vhdl << instance (pi_mult, "pi_mult");
	int wZfull = getSignalByName("Zfull")->width();
	vhdl << tab << declare("Z",wZ) << " <= Zfull" << range (wZfull-1, wZfull-wZ) << ";" << endl; 
#endif


	/*********************************** THE SQUARER **************************************/


	map<string, double> sqr_z_inputDelays;
#if SUBCYCLEPIPELINING
	sqr_z_inputDelays["X"] = pi_mult->getOutputDelay("R");
	sqr_z_inputDelays["Y"] = pi_mult->getOutputDelay("R");
#else
	nextCycle();
#endif

	// vhdl:sqr (Z -> Z2o2)
	// we have no truncated squarer as of now
	/*IntSquarer *sqr_z;
	sqr_z = new IntSquarer (target, wZ);
	oplist.push_back (sqr_z);
	inPortMap (sqr_z, "X", "Z");
	outPortMap (sqr_z, "R", "Z2o2_ext");
	vhdl << instance (sqr_z, "sqr_z");
	// so now we truncate unnecessarily calculated bits of Z2o2_ext
	int wZ2o2 = 2*wZ - (w+g);
	vhdl << declare ("Z2o2",wZ2o2) << " <= Z2o2_ext"
	     << range (wZ-1,wZ-wZ2o2) << ";" << endl;*/
	// so we use a truncated multiplier instead
	IntMultiplier *sqr_z;
	int wZ2o2 = 2*wZ - (w+g)-1;
	if (wZ2o2 < 2)
		wZ2o2 = 2; //for sanity
	vhdl << tab << "-- First truncate the inputs of the multiplier to the precision of the output" << endl;
	vhdl << tab << declare("Z_truncToZ2", wZ2o2) << " <= Z" << range(wZ-1, wZ-wZ2o2) << ";" << endl;
	sqr_z = new IntMultiplier (target, wZ2o2, wZ2o2, wZ2o2, false, 1.0, sqr_z_inputDelays);
	oplist.push_back (sqr_z);
	inPortMap (sqr_z, "Y", "Z_truncToZ2");
	inPortMap (sqr_z, "X", "Z_truncToZ2");
	outPortMap (sqr_z, "R", "Z2o2");
	vhdl << instance (sqr_z, "sqr_z");


	/*********************************** THE CUBING UNIT **************************************/


	// remember the cycle and the critical path here, because we'll get back in time here
	map<string, double> Z3_inputDelays;
#if SUBCYCLEPIPELINING
	setSignalDelay("Z2o2", sqr_z->getOutputDelay("R")) ; // TODO should be done automatically by outPortMap 
	syncCycleFromSignal("Z2o2", sqr_z->getOutputDelay("R")); // TODO Then this would not need to pass a CP
	Z3_inputDelays["X"] = sqr_z->getOutputDelay("R");
	Z3_inputDelays["Y"] = sqr_z->getOutputDelay("R");
#else
	syncCycleFromSignal("Z2o2"); // TODO Then this would not need to pass a CP
	nextCycle();
#endif

	// vhdl:mul (Z, Z2o2 -> Z3)


	int wZ3 = 3*wZ - 2*(w+g) -1; // -1 for the div by 2
	int wZ3o6 = 3*wZ - 2*(w+g) -2;
	if (wZ3 < 2)
		wZ3 = 2;
	if (wZ3o6 < 2)
		wZ3o6 = 2; //using 1 will generate bad vhdl

#if 1
	if(wZ3<=11) {
		vhdl << tab << "-- First truncate Z" << endl;
		vhdl << tab << declare("Z_truncToZ3o6", wZ3o6) << " <= Z" << range(wZ-1, wZ-wZ3o6) << ";" << endl;
		FunctionTable *z3o6Table;
		z3o6Table = new FunctionTable (target, "x^3/6)", wZ3o6, -wZ3o6-2, -3);
		z3o6Table -> changeName(getName() + "_Z3o6Table");
		oplist.push_back (z3o6Table);
		inPortMap (z3o6Table, "X", "Z_truncToZ3o6");
		outPortMap (z3o6Table, "Y", "Z3o6");
		vhdl << instance (z3o6Table, "z3o6Table");
		
	}
	else {
		// TODO: replace all this with an ad-hoc unit
		vhdl << tab << "-- First truncate the inputs of the multiplier to the precision of the output" << endl;
		vhdl << tab << declare("Z2o2_truncToZ3", wZ3) << " <= Z2o2" << range(wZ2o2-1, wZ2o2-wZ3) << ";" << endl;
		vhdl << tab << declare("Z_truncToZ3", wZ3) << " <= Z" << range(wZ-1, wZ-wZ3) << ";" << endl;
		IntMultiplier *Z3;
		Z3 = new IntMultiplier (target, wZ3, wZ3, wZ3,  
		                              0.5, false);
		oplist.push_back (Z3);
		outPortMap (Z3, "R", "Z3o2");
		inPortMap (Z3, "Y", "Z2o2_truncToZ3");
		inPortMap (Z3, "X", "Z_truncToZ3");
		vhdl << instance (Z3, "z3_compute");
		syncCycleFromSignal("Z3o2", Z3->getOutputDelay("R"));
		

		map<string, double> cdiv_3_inputDelays;
#if SUBCYCLEPIPELINING
		cdiv_3_inputDelays["X"] = Z3->getOutputDelay("R");
#else
		nextCycle();
#endif
		
		IntConstDiv *cdiv_3;
		cdiv_3 = new IntConstDiv (target, wZ3, 3, -1, cdiv_3_inputDelays);
		oplist.push_back (cdiv_3);
		outPortMap (cdiv_3, "Q", "Z3o6");
		inPortMap (cdiv_3, "X", "Z3o2");
		vhdl << instance (cdiv_3, "cdiv_3");
		syncCycleFromSignal("Z3o6", cdiv_3->getOutputDelay("Q"));

	}
#else
	vhdl << tab << declare("Z_truncToZ3", wZ3) << " <= Z" << range(wZ-1, wZ-wZ3) << ";" << endl;
	ProductIR z3o6 = ProductIR::identity(wZ3).toPow(3).simplifyInPlace().div(6).quo.setMSB(-3);
	// target wZ3o6, not wZ3 when output is concerned
	// 1 ulp of output trunc as of now
	ProductIR z3o6_t = z3o6.truncate (mpz_class(1) << (wZ3*3-wZ3o6));
	size_t wZ3o6_ext = z3o6_t.data.size();
	GenericBinaryPolynomial *z3o6gbp;
	z3o6gbp = new GenericBinaryPolynomial (target, z3o6_t, Z3_inputDelays);
	oplist.push_back (z3o6gbp);
	outPortMap (z3o6gbp, "R", "Z3o6_before_trunc");
	inPortMap (z3o6gbp, "X", "Z_truncToZ3");
	vhdl << instance (z3o6gbp, "z3o6gbp");
	syncCycleFromSignal ("Z3o6_before_trunc", z3o6gbp->getOutputDelay("R"));
	if (wZ3o6_ext >= wZ3o6) { 
		vhdl << tab << declare("Z3o6",wZ3o6) << " <= Z3o6_before_trunc"
		     << range(wZ3o6_ext-1, wZ3o6_ext-wZ3o6) << ";\n";
	} else {
		vhdl << tab << declare("Z3o6",wZ3o6)
		     << " <= Z3o6_before_trunc & \"" 
		     << std::string (wZ3o6-wZ3o6_ext,'0')
		     << "\";" << endl;
	}
#endif


	/*********************************** Z-Z3o6 **************************************/


#if SUBCYCLEPIPELINING
#else
		nextCycle();
#endif


	// vhdl:sub (Z, Z3o6 -> sinZ)
	manageCriticalPath(target->adderDelay(wZ));
	vhdl << tab << declare ("SinZ", wZ)
	     << " <= Z - Z3o6;" << endl;
	setSignalDelay("SinZ", getCriticalPath());



	// and now, evaluate Sin Y_in and Cos Y_in
	// Cos Y_in:
	// vhdl:slr (Z2o2 -> Z2o2)



	/*********************************** Reconstruction of Sine **************************************/

	// First get back to the cycle of Z2
#if SUBCYCLEPIPELINING
	setCycleFromSignal("Z2o2", getSignalDelay("Z2o2"));
#else
	setCycleFromSignal("Z2o2");
	nextCycle();
#endif

	// // vhdl:id (CosPiA -> C_out_1)
	// vhdl:mul (Z2o2, CosPiA -> Z2o2CosPiA)
	vhdl << tab << "-- First truncate the larger input of the multiplier to the precision of the output" << endl;
	vhdl << tab << declare("CosPiA_truncToZ2o2", wZ2o2) << " <= CosPiA" << range(w+g-1, w+g-wZ2o2) << ";" << endl;
	IntMultiplier *c_out_2;
	c_out_2 = new IntMultiplier (target, wZ2o2, wZ2o2, wZ2o2, false, 1.0);
	oplist.push_back (c_out_2);
	inPortMap (c_out_2, "X", "Z2o2");
	inPortMap (c_out_2, "Y", "CosPiA_truncToZ2o2");
	outPortMap (c_out_2, "R", "Z2o2CosPiA");
	vhdl << instance (c_out_2, "c_out_2_compute");


#if SUBCYCLEPIPELINING
	syncCycleFromSignal("Z2o2CosPiA", c_out_2->getOutputDelay("R"));
#else
	syncCycleFromSignal("Z2o2CosPiA");
	nextCycle();
#endif

	// get back to the cycle of sinZ, certainly later than the SinPiA
#if SUBCYCLEPIPELINING
	setCycleFromSignal("SinZ", getSignalDelay("SinZ")); 
#else
	setCycleFromSignal("SinZ"); 
	nextCycle();
#endif

	// vhdl:mul (SinZ, SinPiA -> SinZSinPiA)
	IntMultiplier *c_out_3;

	vhdl << tab << "-- First truncate the larger input of the multiplier to the precision of the output" << endl;
	vhdl << tab << declare("SinPiA_truncToZ", wZ) << " <= SinPiA" << range(w+g-1, w+g-wZ) << ";" << endl;
	c_out_3 = new IntMultiplier (target, wZ, wZ, wZ, false);
	oplist.push_back (c_out_3);
	inPortMap (c_out_3, "Y", "SinPiA_truncToZ");
	inPortMap (c_out_3, "X", "SinZ");
	outPortMap (c_out_3, "R", "SinZSinPiA");
	vhdl << instance (c_out_3, "c_out_3_compute");

	setCycleFromSignal ("Z2o2CosPiA");
	nextCycle();

	// TODO: critical path supposed suboptimal (but don't know how to fix)
	manageCriticalPath(target->localWireDelay() + target->adderDelay(w+g));
	vhdl << tab << declare ("CosZCosPiA_plus_rnd", w+g)
	     << " <= CosPiA - Z2o2CosPiA;" << endl;

	syncCycleFromSignal ("SinZSinPiA");
	nextCycle();
	
	manageCriticalPath(target->localWireDelay() + target->adderDelay(w+g));
	vhdl << tab << declare ("C_out_rnd_aux", w+g)
	     << " <= CosZCosPiA_plus_rnd - SinZSinPiA;" << endl;



	/*********************************** Reconstruction of Sine **************************************/

	// Sin Y_in:

	// First get back to the cycle of Z2 (same as Cos_y_red):
	// it is certainly later than SinPiA

#if SUBCYCLEPIPELINING
	// TODO
#else
	setCycleFromSignal("Z2o2");
#endif

	// // vhdl:id (SinPiA -> S_out_1)
	// vhdl:mul (Z2o2, SinPiA -> Z2o2SinPiA)
	vhdl << tab << "-- First truncate the larger input of the multiplier to the precision of the output" << endl;
	vhdl << tab << declare("SinPiA_truncToZ2o2", wZ2o2) << " <= SinPiA" << range(w+g-1, w+g-wZ2o2) << ";" << endl;
	IntMultiplier *s_out_2;
	s_out_2 = new IntMultiplier (target, wZ2o2, wZ2o2, wZ2o2, false);
	oplist.push_back (s_out_2);
	inPortMap (s_out_2, "X", "Z2o2");
	inPortMap (s_out_2, "Y", "SinPiA_truncToZ2o2");
	outPortMap (s_out_2, "R", "Z2o2SinPiA");
	vhdl << instance (s_out_2, "s_out_2_compute");
	syncCycleFromSignal("Z2o2SinPiA");


	// get back to the cycle of SinZ, certainly later than the CosPiA
	setCycleFromSignal("SinZ", getSignalDelay("SinZ"));

	// vhdl:mul (SinZ, CosPiA -> SinZCosPiA)
	nextCycle();

	vhdl << tab << "-- First truncate the larger input of the multiplier to the precision of the output" << endl;
	vhdl << tab << declare("CosPiA_truncToSinZ", wZ) << " <= CosPiA" << range(w+g-1, w+g-wZ) << ";" << endl;

	IntMultiplier *s_out_3;
	s_out_3 = new IntMultiplier (target, wZ, wZ, wZ, false);
	oplist.push_back (s_out_3);
	inPortMap (s_out_3, "X", "SinZ");
	inPortMap (s_out_3, "Y", "CosPiA_truncToSinZ");
	outPortMap (s_out_3, "R", "SinZCosPiA");
	vhdl << instance (s_out_3, "s_out_3_compute");
	syncCycleFromSignal("SinZCosPiA");

	setCycleFromSignal ("Z2o2SinPiA");
	nextCycle();

	manageCriticalPath(target->localWireDelay() + target->adderDelay(w+g));
	vhdl << tab << declare ("CosZSinPiA_plus_rnd", w+g)
	     << " <= SinPiA - Z2o2SinPiA;" << endl;

	syncCycleFromSignal ("SinZCosPiA");
	nextCycle();

	manageCriticalPath(target->localWireDelay() + target->adderDelay(w+g));
	vhdl << tab << declare ("S_out_rnd_aux", w+g)
	     << " <= CosZSinPiA_plus_rnd + SinZCosPiA;" << endl;


	//Final synchronization
	syncCycleFromSignal("C_out_rnd_aux");

	vhdl << tab << declare ("C_out", w)
	     << " <= C_out_rnd_aux" << range (w+g-1, g) << ';' << endl;
	vhdl << tab << declare ("S_out", w)
	     << " <= S_out_rnd_aux" << range (w+g-1, g) << ';' << endl;

	// now just add the signs to the signals and convert them to
	// 2's complement
	vhdl << tab << declare ("S_wo_sgn", w)
	     << " <= C_out when Exch = '1' else S_out;" << endl;
	vhdl << tab << declare ("C_wo_sgn", w)
	     << " <= S_out when Exch = '1' else C_out;" << endl;
	vhdl << tab << declare ("S_wo_sgn_ext", w+1)
	     << " <= '0' & S_wo_sgn;" << endl
	     << tab << declare ("C_wo_sgn_ext", w+1)
	     << " <= '0' & C_wo_sgn;" << endl;
	vhdl << tab << declare ("S_wo_sgn_neg", w+1)
	     << " <= (not S_wo_sgn_ext) + 1;" << endl;
	vhdl << tab << declare ("C_wo_sgn_neg", w+1)
	     << " <= (not C_wo_sgn_ext) + 1;" << endl;
	//vhdl << tab << declare("S_sgn")
	//     << " <= X_sgn;" << endl;
	vhdl << tab << declare("C_sgn")
	     << " <= Q xor X_sgn;" << endl;
	vhdl << tab << "S <= S_wo_sgn_ext when X_sgn = '0'"
	     << " else S_wo_sgn_neg;" << endl
	     << tab << "C <= C_wo_sgn_ext when C_sgn = '0'"
	     << " else C_wo_sgn_neg;" << endl;

	REPORT(INFO, " wA=" << wA <<" wZ=" << wZ <<" wZ2=" << wZ2o2 <<" wZ3=" << wZ3 );

	// For LateX in the paper
	// cout << "     " << w <<  "   &   " << g <<"   &   " << wA << "   &   " << wZ << "   &   " << wZ2o2 << "   &   " << wZ3 << "   \\\\ \n \\hline" <<  endl;
};





void FixSinCos::emulate(TestCase * tc)
{
	mpz_class sx = tc->getInputValue ("X");
	mpfr_t x, sind, cosd, sinu, cosu, pixd, pixu, one_minus_ulp;
	mpz_t sind_z, cosd_z, sinu_z, cosu_z;
	//one extra bit of precision because we temporarily hold the sign
	//in the mantissa
	mpfr_init2 (x, 2+w);
	mpz_init2 (sind_z, 1+w);
	mpz_init2 (cosd_z, 1+w);
	mpz_init2 (sinu_z, 1+w);
	mpz_init2 (cosu_z, 1+w);
	mpfr_inits (sind, cosd, sinu, cosu, pixd, pixu, (mpfr_ptr) 0);
	mpfr_init2 (one_minus_ulp, 1+w);
	//these 3 following roundings are exact
	mpfr_set_ui (one_minus_ulp, 1UL, GMP_RNDD);
	mpfr_div_2ui (one_minus_ulp, one_minus_ulp, w, GMP_RNDD);
	mpfr_ui_sub (one_minus_ulp, 1UL, one_minus_ulp, GMP_RNDD);

	mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); // this rounding is exact
	mpfr_div_2si (x, x, w, GMP_RNDD); // this rounding is acually exact
	// now handle the sign of x
	if (mpfr_cmp_ui (x, 1UL) >= 0) { // if (x >= 1.f)
		// 2's complement: just sub 2
		mpfr_add_si (x, x, -2, GMP_RNDD); //this rnd is exact
	}
	int i=0, ep; // ep: extra precision
	do {
		ep = 1 << i;
		mpfr_set_prec (sind, 2+w+ep);
		mpfr_set_prec (cosd, 2+w+ep);
		mpfr_set_prec (sinu, 2+w+ep);
		mpfr_set_prec (cosu, 2+w+ep);
		mpfr_set_prec (pixd, 2+w+ep);
		mpfr_set_prec (pixu, 2+w+ep);
		mpfr_const_pi (pixd, GMP_RNDD);
		mpfr_const_pi (pixu, GMP_RNDU);
		mpfr_mul (pixd, pixd, x, GMP_RNDD);
		mpfr_mul (pixu, pixu, x, GMP_RNDU);
		if (mpfr_cmp_ui_2exp (x, 1UL, -1) < 0
		 && mpfr_cmp_si_2exp (x, -1L, -1) > 0) { // if (|x| < 0.5f)
			// then sin is increasing near x
			mpfr_sin (sind, pixd, GMP_RNDD);
			mpfr_sin (sinu, pixu, GMP_RNDU);
		} else {
			// then sin is decreasing near x
			mpfr_sin (sind, pixu, GMP_RNDD); // use the upper x for the lower sin
			mpfr_sin (sinu, pixd, GMP_RNDU);
		}
		if (mpfr_cmp_ui (x, 0UL) > 0) {
			mpfr_cos (cosd, pixu, GMP_RNDD); // cos decreases from 0 to pi
			mpfr_cos (cosu, pixd, GMP_RNDU);
		} else {
			mpfr_cos (cosd, pixd, GMP_RNDD); // cos increases from -pi to 0
			mpfr_cos (cosu, pixu, GMP_RNDU);
		}	
		// multiply by the (1-ulp) factor now
		mpfr_mul (sind, sind, one_minus_ulp, GMP_RNDD);
		mpfr_mul (cosd, cosd, one_minus_ulp, GMP_RNDD);
		mpfr_mul (sinu, sinu, one_minus_ulp, GMP_RNDU);
		mpfr_mul (cosu, cosu, one_minus_ulp, GMP_RNDU);
		// if the cosine is <0 we must add 2 to it
		// (2's complement)
		if (mpfr_cmp_ui (cosd, 0) < 0) {
			mpfr_add_ui (cosd, cosd, 2UL, GMP_RNDD); // exact rnd
		}
		// same as before
		if (mpfr_cmp_ui (cosu, 0) < 0) {
			mpfr_add_ui (cosu, cosu, 2UL, GMP_RNDU);
		}
		// and we have to do the same for sin
		if (mpfr_cmp_ui (sind, 0) < 0) {
			mpfr_add_ui (sind, sind, 2UL, GMP_RNDD); // exact rnd
		}
		if (mpfr_cmp_ui (sinu, 0) < 0) {
			mpfr_add_ui (sinu, sinu, 2UL, GMP_RNDU);
		}
		mpfr_mul_2si (sind, sind, w, GMP_RNDD); // exact rnd here
		mpfr_mul_2si (cosd, cosd, w, GMP_RNDD);
		mpfr_mul_2si (sinu, sinu, w, GMP_RNDU);
		mpfr_mul_2si (cosu, cosu, w, GMP_RNDU);
		mpfr_get_z (sind_z, sind, GMP_RNDD); // there can be a real rounding here
		mpfr_get_z (cosd_z, cosd, GMP_RNDD);
		mpfr_get_z (sinu_z, sinu, GMP_RNDU); // there can be a real rounding here
		mpfr_get_z (cosu_z, cosu, GMP_RNDU);
		// now we test if the upper results are the lower ones + 1
		// as we want them to
		mpz_sub_ui (sinu_z, sinu_z, 1UL);
		mpz_sub_ui (cosu_z, cosu_z, 1UL);
		if (mpz_cmp (sind_z, sinu_z) == 0 &&
		    mpz_cmp (cosd_z, cosu_z) == 0) {
			// the rounding are really what we want
			// TODO: detect wraparounds, or do the sign conversion
			// later than there
			break;
		} else {
			i++;
		}
	} while (i<16); // for sanity
	if (i==16)
		REPORT(INFO,"Computation failure with 65536 extra"
		            "precision bits, something got wrong\n");
	mpz_class sind_zc (sind_z), cosd_zc (cosd_z);
	// and since the operator does only faithful rounding
	// we add also as expected results the upper roundings
	mpz_class neg_zero (1); neg_zero <<= w;
	mpz_class ones (neg_zero - 1);
	mpz_class all_ones ((neg_zero << 1) - 1);
	// ones is actually the greatest representable signal
	tc->addExpectedOutput ("S", sind_zc);
	tc->addExpectedOutput ("C", cosd_zc);
	if (sind_zc != ones) {
		tc->addExpectedOutput ("S", (sind_zc+1) & all_ones);
	}
	if (cosd_zc != ones) {
		tc->addExpectedOutput ("C", (cosd_zc+1) & all_ones);
	}
	mpfr_clears (x, sind, cosd, sinu, cosu,
	             pixd, pixu, one_minus_ulp, NULL);
	mpz_clear (sind_z);
	mpz_clear (cosd_z);
	mpz_clear (sinu_z);
	mpz_clear (cosu_z);
}


void FixSinCos::buildStandardTestCases(TestCaseList * tcl)
{
	TestCase* tc;
	tc = new TestCase (this);
	tc -> addInput ("X",mpz_class(0));
	emulate(tc);
	tcl->add(tc);
	tc = new TestCase (this);
	tc -> addInput ("X",mpz_class(1));
	emulate(tc);
	tcl->add(tc);
	mpz_class ones(1);
	ones <<= w;
	ones -= 1;
	tc = new TestCase (this);
	tc -> addInput ("X",ones);
	emulate(tc);
	tcl->add(tc);
}

#endif // SOLLYA

