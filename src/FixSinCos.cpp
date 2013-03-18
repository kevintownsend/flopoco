#include "ConstMult/FixRealKCM.hpp"
#include "IntMultiplier.hpp"
#include "FixedPointFunctions/FunctionTable.hpp"
#include "IntConstDiv.hpp"
#include "BitHeap.hpp"
#include "IntMultipliers/FixSinPoly.hpp"

// TODOs 
// Compare not-ing Y and negating it properly

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

FixSinCos::FixSinCos(Target * target, int w_, float ratio):Operator(target), w(w_)
{
	srcFileName="FixSinCos";
	
	// definition of the name of the operator
	ostringstream name;
	name << "FixSinCos_" << w ;
	setName(name.str());

	setCopyrightString("Florent de Dinechin, Guillaume Sergent (2012)");




	// everybody needs many digits of Pi
	mpfr_init2(constPi, 10*w);
	mpfr_const_pi( constPi, GMP_RNDN);
	
	//compute the scale factor		
	mpfr_init2(scale, w+1);
	mpfr_set_d(scale, -1.0, GMP_RNDN);           // exact
	mpfr_mul_2si(scale, scale, -w, GMP_RNDN); // exact
	mpfr_add_d(scale, scale, 1.0, GMP_RNDN);     // exact
	REPORT(DEBUG, "scale=" << printMPFR(scale, 15));
	

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

	// notY = .25 - y
	// do a logic ~ at a cost of 1 ulp ?
	// It saves one carry-propagate latency (one cycle)
	// but enlarges the constant multiplier input by g bits
	// and adds one ulp of error
	manageCriticalPath(target->localWireDelay(w-2) + target->lutDelay());
	vhdl << tab << declare ("notY", w-2) << " <= " << "not Y;" << endl;

	// we need to know the number of guard bits _now_ to have a good Y_in
	// some precision is lost with the optimized multipliers
	const int g=6;
	//const int g=4; //guard bits
	//const int g=3; //guard bits
	//we take 4 guard bits even if error < 8 ulp because rounding will take
	//another half final ulp
	int wIn = w-2+g;

	// Y_in = / notY if O=1
	//        \ Y if O=0
	// and we extend the precision to make it look as if notY was
	// 0.[1, wIn times] - Y: 1/2**g error reduction in notY
	// (which should be arithmetic 1-Y)
	vhdl << tab << declare ("Y_in",wIn) << " <= notY & "
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
	int wY = wIn-wA;
	int wZ = wY+2;
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
		sin_table = new FunctionTable (target, 
		                               omu.str() + " * sin(pi*x/4) + " + halfulp.str(),
		                               wA, 
		                               -(w+g), 
		                               -1);
		cos_table = new FunctionTable (target, 
		                               omu.str() + " * cos(pi*x/4) + " + halfulp.str(),
		                               wA, 
		                               -(w+g), 
		                               -1);
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
	sqr_z = new IntMultiplier (target, wZ2o2, wZ2o2, wZ2o2, false, ratio, sqr_z_inputDelays);
	oplist.push_back (sqr_z);
	inPortMap (sqr_z, "Y", "Z_truncToZ2");
	inPortMap (sqr_z, "X", "Z_truncToZ2");
	outPortMap (sqr_z, "R", "Z2o2");
	vhdl << instance (sqr_z, "sqr_z");


	/*********************************** Z-Z^3/6 **************************************/

	int wZ3 = 3*wZ - 2*(w+g) -1; // -1 for the div by 2
	int wZ3o6 = 3*wZ - 2*(w+g) -2;
	if (wZ3 < 2)
		wZ3 = 2;
	if (wZ3o6 < 2)
		wZ3o6 = 2; //using 1 will generate bad vhdl
	REPORT(DETAILED, "wZ3 = " << wZ3 << "            wZ3o6 = " << wZ3o6 );
	
	if(wZ3<=10) {
		vhdl << tab << "-- First truncate Z" << endl;
		vhdl << tab << declare("Z_truncToZ3o6", wZ3o6) << " <= Z" << range(wZ-1, wZ-wZ3o6) << ";" << endl;
		FunctionTable *z3o6Table;
		z3o6Table = new FunctionTable (target, "x^3/6", wZ3o6, -wZ3o6-2, -3);
		z3o6Table -> changeName(getName() + "_Z3o6Table");
		oplist.push_back (z3o6Table);
		inPortMap (z3o6Table, "X", "Z_truncToZ3o6");
		outPortMap (z3o6Table, "Y", "Z3o6");
		vhdl << instance (z3o6Table, "z3o6Table");
		
		manageCriticalPath(target->adderDelay(wZ));
		vhdl << tab << declare ("SinZ", wZ)
		     << " <= Z - Z3o6;" << endl;
		setSignalDelay("SinZ", getCriticalPath());
		
	}
	else {
		// TODO: replace all this with an ad-hoc unit
		FixSinPoly *fsp =new FixSinPoly(target, 
		                                -wA-1, //msbin
		                                -w-g, // lsbin
		                                true, // truncated
		                                -wA-1, // msbOut_ = 0,
		                                -w-g, // lsbout
		                                false);
		oplist.push_back (fsp);
		inPortMap (fsp, "X", "Z");
		outPortMap(fsp, "R", "SinZ");
		vhdl << instance (fsp, "ZminusZ3o6");
	}


	// vhdl:sub (Z, Z3o6 -> SinZ)



	// and now, evaluate Sin Y_in and Cos Y_in
	// Cos Y_in:
	// vhdl:slr (Z2o2 -> Z2o2)



	/*********************************** Reconstruction of cosine **************************************/

	// First get back to the cycle of Z2
#if SUBCYCLEPIPELINING
	setCycleFromSignal("Z2o2", getSignalDelay("Z2o2"));
#else
	setCycleFromSignal("Z2o2");
	nextCycle();
#endif







#if 1 // No bit heap
	// // vhdl:id (CosPiA -> C_out_1)
	// vhdl:mul (Z2o2, CosPiA -> Z2o2CosPiA)

	vhdl << tab << "--  truncate the larger input of each multiplier to the precision of its output" << endl;
	vhdl << tab << declare("CosPiA_truncToZ2o2", wZ2o2) << " <= CosPiA" << range(w+g-1, w+g-wZ2o2) << ";" << endl;
	IntMultiplier *c_out_2;
	c_out_2 = new IntMultiplier (target, wZ2o2, wZ2o2, wZ2o2, false, ratio);
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

	// get back to the cycle of SinZ, certainly later than the SinPiA
#if SUBCYCLEPIPELINING
	setCycleFromSignal("SinZ", getSignalDelay("SinZ")); 
#else
	setCycleFromSignal("SinZ"); 
	nextCycle();
#endif

	vhdl << tab << "--  truncate the larger input of each multiplier to the precision of its output" << endl;
	vhdl << tab << declare("SinPiA_truncToZ", wZ) << " <= SinPiA" << range(w+g-1, w+g-wZ) << ";" << endl;


	// vhdl:mul (SinZ, SinPiA -> SinZSinPiA)
	IntMultiplier *c_out_3;

	c_out_3 = new IntMultiplier (target, wZ, wZ, wZ, false, ratio);
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

	vhdl << tab << declare ("C_out", w)
	     << " <= C_out_rnd_aux" << range (w+g-1, g) << ';' << endl;




	// ---------------------------------------------
#else //Bit heap computing Cos Z  ~   CosPiA - Z2o2*cosPiA - sinZ*SinPiA
	//
	// cosPiA   xxxxxxxxxxxxxxxxxxggggg
	//

#define TINKERCOS 1
#if TINKERCOS
	int gMult=0;
#else
	int g1 = IntMultiplier::neededGuardBits(wZ, wZ, wZ);
	int g2 = IntMultiplier::neededGuardBits(wZ2o2, wZ2o2, wZ2o2);
	int gMult=max(g1,g2);
#endif

	REPORT(0, "wZ2o2=" << wZ2o2 << "    wZ=" << wZ << "    g=" << g << "    gMult=" << gMult);
	BitHeap* bitHeapCos = new BitHeap(this, w+g+gMult, "Sin"); 

	// Add CosPiA to the bit heap
	bitHeapCos -> addUnsignedBitVector(gMult, "CosPiA", w+g);

	vhdl << tab << "--  truncate the larger input of each multiplier to the precision of its output" << endl;
	vhdl << tab << declare("CosPiA_truncToZ2o2", wZ2o2) << " <= CosPiA" << range(w+g-1, w+g-wZ2o2) << ";" << endl;
	vhdl << tab << "--  truncate the larger input of each multiplier to the precision of its output" << endl;
	vhdl << tab << declare("SinPiA_truncToZ", wZ) << " <= SinPiA" << range(w+g-1, w+g-wZ) << ";" << endl;


#if TINKERCOS
	setCycleFromSignal ("Z2o2");
	nextCycle();
	vhdl <<  tab << declare("Z2o2CosPiA", 2*wZ2o2) << " <= Z2o2 * CosPiA_truncToZ2o2" << ";" << endl;
	nextCycle();
	// add it to the bit heap
	
	for (int i=0; i<wZ2o2-1; i++)
		bitHeapCos->addBit(i, "not Z2o2CosPiA"+of(i+wZ2o2));
	bitHeapCos->addBit(wZ2o2-1, "Z2o2CosPiA"+of(wZ2o2+wZ2o2-1));
	for (int i=wZ2o2-1; i<w+g+gMult; i++)
		bitHeapCos->addConstantOneBit(i);
	bitHeapCos->addConstantOneBit(0);


	setCycleFromSignal ("SinZ");
	nextCycle();
	vhdl <<  tab << declare("SinZSinPiA", 2*wZ) << " <=  SinZ *SinPiA_truncToZ" << ";" << endl;
	nextCycle();
	
	// add it to the bit heap	
	for (int i=0; i<wZ-1; i++)
		bitHeapCos->addBit(i, "not SinZSinPiA"+of(i+wZ));
	bitHeapCos->addBit(wZ-1, "SinZSinPiA"+of(wZ+wZ-1));
	for (int i=wZ-1; i<w+g+gMult; i++)
		bitHeapCos->addConstantOneBit(i);
	bitHeapCos->addConstantOneBit(0);

	//	vhdl <<  tab << declare("Z2o2CosPiA", 2*wZ) << " <= " << range(w+g-1, w+g-wZ) << ";" << endl;
#else
	// First virtual multiplier
	new IntMultiplier (this,
	                   bitHeapCos,
	                   getSignalByName("Z2o2"),
	                   getSignalByName("CosPiA_truncToZ2o2"),
	                   wZ2o2, wZ2o2, wZ2o2,
	                   gMult,
	                   true, // negate
	                   false, // signed inputs
	                   ratio);

	     

	// Second virtual multiplier
	new IntMultiplier (this,
	                   bitHeapCos,
	                   getSignalByName("SinZ"),
	                   getSignalByName("SinPiA_truncToZ"),
	                   wZ, wZ, wZ,
	                   gMult,
	                   true, // negate
	                   false, // signed inputs
	                   ratio
	                   );
#endif
	     
	// The round bit is in the table already
	bitHeapCos -> generateCompressorVHDL();	
	vhdl << tab << declare ("C_out", w) << " <= " << bitHeapCos -> getSumName() << range (w+g+gMult-1, g+gMult) << ';' << endl;

#endif

	/*********************************** Reconstruction of sine **************************************/
 //Bit heap computing   SinPiA - Z2o2*sinPiA + sinZ*CosPiA
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
	s_out_2 = new IntMultiplier (target, wZ2o2, wZ2o2, wZ2o2, false, ratio);
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
	s_out_3 = new IntMultiplier (target, wZ, wZ, wZ, false, ratio);
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
	syncCycleFromSignal("C_out");

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




FixSinCos::~FixSinCos(){
		mpfr_clears (scale, constPi, NULL);		
	 };






void FixSinCos::emulate(TestCase * tc)
{
	// TODO eventually have only one shared FixSinCos emulate code

	mpfr_t z, rsin, rcos;
		mpz_class sin_z, cos_z;
		mpfr_init2(z, 10*w);
		mpfr_init2(rsin, 10*w); 
		mpfr_init2(rcos, 10*w); 

		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		
		/* Compute correct value */
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); //  exact
		mpfr_div_2si (z, z, w, GMP_RNDN); // exact
	
		// No need to manage sign bit etc: modulo 2pi is the same as modulo 2 in the initial format
		mpfr_mul(z, z, constPi, GMP_RNDN);

		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		mpfr_mul(rsin, rsin, scale, GMP_RNDN);
		mpfr_mul(rcos, rcos, scale, GMP_RNDN);

		mpfr_add_d(rsin, rsin, 6.0, GMP_RNDN); // exact rnd here
		mpfr_add_d(rcos, rcos, 6.0, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDN); // exact rnd here

		// Rounding down
		mpfr_get_z (sin_z.get_mpz_t(), rsin, GMP_RNDD); // there can be a real rounding here
		mpfr_get_z (cos_z.get_mpz_t(), rcos, GMP_RNDD); // there can be a real rounding here
		sin_z -= mpz_class(6)<<(w);
		cos_z -= mpz_class(6)<<(w);

		tc->addExpectedOutput ("S", sin_z);
		tc->addExpectedOutput ("C", cos_z);

		// Rounding up
		mpfr_get_z (sin_z.get_mpz_t(), rsin, GMP_RNDU); // there can be a real rounding here
		mpfr_get_z (cos_z.get_mpz_t(), rcos, GMP_RNDU); // there can be a real rounding here
		sin_z -= mpz_class(6)<<(w);
		cos_z -= mpz_class(6)<<(w);

		tc->addExpectedOutput ("S", sin_z);
		tc->addExpectedOutput ("C", cos_z);
		
		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
}


void FixSinCos::buildStandardTestCases(TestCaseList * tcl)
{
		mpfr_t z;
		mpz_class zz;
		TestCase* tc;
				
		mpfr_init2(z, 10*w);
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X", mpz_class(0));
		emulate(tc);
		tcl->add(tc);
					
		tc = new TestCase (this);
		tc->addComment("Pi/4-eps");
		mpfr_set_d (z, 0.24, GMP_RNDD); 
		mpfr_mul_2si (z, z, w-1, GMP_RNDD); 
		mpfr_get_z (zz.get_mpz_t(), z, GMP_RNDD);  
		tc -> addInput ("X", zz);
		emulate(tc);
		tcl->add(tc);
				
		tc = new TestCase (this);
		tc->addComment("Pi/6");
		mpfr_set_d (z, 0.166666666666666666666666666666666, GMP_RNDD); 
		mpfr_mul_2si (z, z, w-1, GMP_RNDD); 
		mpfr_get_z (zz.get_mpz_t(), z, GMP_RNDD);  
		tc -> addInput ("X", zz);
		emulate(tc);
		tcl->add(tc);
				
		tc = new TestCase (this);
		tc->addComment("Pi/3");
		mpfr_set_d (z, 0.333333333333333333333333333333, GMP_RNDD); 
		mpfr_mul_2si (z, z, w-1, GMP_RNDD); 
		mpfr_get_z (zz.get_mpz_t(), z, GMP_RNDD);  
		tc -> addInput ("X", zz);
		emulate(tc);
		tcl->add(tc);
				
		
		mpfr_clears (z, NULL);

}

#endif // SOLLYA

