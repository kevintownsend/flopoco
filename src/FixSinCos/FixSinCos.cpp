
// TODOs 
// Move FixSinPoly here!
// Compare not-ing Y and negating it properly

// works only with sollya
#ifdef HAVE_SOLLYA
#include "../ConstMult/FixRealKCM.hpp"
#include "../IntMultiplier.hpp"
#include "../FixFunctions/FunctionTable.hpp"
#include "../IntConstDiv.hpp"
#include "../BitHeap.hpp"
#include "../IntMultipliers/FixSinPoly.hpp"

#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers */
#include "mpfr.h"
#include "FixSinCos.hpp"
#include "Table.hpp"


using namespace std;
using namespace flopoco;


#define SUBCYCLEPIPELINING 0
#define USEBITHEAP 1
#define LARGE_PREC 1000 // 1000 bits should be enough for everybody
#define LOW_PREC 0


////////////////////////////////////// SinCosTable ///////////////////
// argRedCase=1 for full table, 4 for quadrant
FixSinCos::SinCosTable::SinCosTable(Target* target, int wIn, int wOut, int argRedCase_, FixSinCos* parent_):
	Table(target, wIn, wOut),argRedCase(argRedCase_), parent(parent_){
	ostringstream name;
	srcFileName = "FixSinCos::SinCosTable";
	name << "SinCosTable_" << wIn << "_" << wOut;
	setName(name.str());
	//	outDelayMap["Y"]=target->RMADelay();
}

FixSinCos::SinCosTable::~SinCosTable(){
};



mpz_class FixSinCos::SinCosTable::function (int x){
	mpz_class sin,cos,testSin, testCos;
	mpfr_t a, c, s;

	mpfr_init2(c, LARGE_PREC); // cosine
	mpfr_init2(s, LARGE_PREC); // sine

	// 2's complement convertion of x into xs
	int xs=x;
	if ( (xs>>(wIn-1)) && (argRedCase==1) )
		xs-=(1<<wIn);

	// a will be used for the value in xs
	mpfr_init2(a,10*wIn-1); 
	mpfr_set_si(a, xs, GMP_RNDN);

	//REPORT(0,"Evaluating function on point x="<<x<<" positive value xs="<<xs<<" converted value a="<<printMPFR(a, 10));

	//divide by 2^w then we get a true fixpoint number between -1 and 1.
	if (argRedCase==1){
		if (xs>>(wIn-1))
			xs-=(1<<wIn);
		mpfr_div_2si(a,a,wIn-1, GMP_RNDN);
		//REPORT(0,"a divided by 2^"<<wIn<<" a="<<printMPFR(a,10));
	}

	else if (argRedCase==4)	{ //quadrant
		mpfr_div_2si(a,a,wIn+1, GMP_RNDN);
		//REPORT(0,"a divided by 2^"<<wIn<<" a="<<printMPFR(a,10));
	}

	else	{
		mpfr_div_2si(a,a,wIn+2, GMP_RNDN);
	}
	mpfr_mul(a, a, parent->constPi, GMP_RNDN);

	mpfr_sin_cos(s, c, a, GMP_RNDN); //function evaluation on point x

	mpfr_mul(s, s, parent->scale, GMP_RNDN); //rescale the sine
	mpfr_mul(c, c, parent->scale, GMP_RNDN); //rescale the cosine

	//REPORT(0," s="<<printMPFR(s,10)<<"; c="<<printMPFR(c,10));

	mpfr_mul_2si(s,s,wOut/2-1, GMP_RNDN); //get the true value
	mpfr_mul_2si(c,c,wOut/2-1, GMP_RNDN); //get the true value

	mpfr_get_z(sin.get_mpz_t(),s, GMP_RNDN); //rounding into mpz
	mpfr_get_z(cos.get_mpz_t(),c, GMP_RNDN);

	//REPORT(0,"Calculated values before 2's complement test: sin="<<sin.get_mpz_t()<<"; cos="<<cos.get_mpz_t());

	// no more need intermediates a, c, and s
	mpfr_clears(a, c, s, NULL);

	// prepare 2's complement: testSin/testCos set to the heaviest bit of sin or cos	
	mpz_div_2exp(testSin.get_mpz_t(),sin.get_mpz_t(), (wOut/2)-1);
	mpz_div_2exp(testCos.get_mpz_t(),cos.get_mpz_t(), (wOut/2)-1);

	// check if negative, then 2's complement
	if(mpz_cmp_si(sin.get_mpz_t(),0)<0){
		//		mpz_com(sin.get_mpz_t(),sin.get_mpz_t());
		sin+=mpz_class(1)<<wOut/2; 
	}

	if (mpz_cmp_si(cos.get_mpz_t(),0)<0){
		//mpz_com(cos.get_mpz_t(),cos.get_mpz_t());
		cos+=mpz_class(1)<<wOut/2; 
	}

	// REPORT(0," function() returns. Value: "<<(sin+(cos<<wIn))<<" ( sin=" << sin<<" , cos="<<cos<<  " )");
	return ( cos  + ( sin << wOut/2 ) ); 

}







////////////////////////////////////// FixSinCos ///////////////////

FixSinCos::FixSinCos(Target * target, int w_, float ratio):Operator(target), w(w_){
	int g=-42; // silly value from the start, because so many different paths may assign it (or forget to do so) 
	srcFileName="FixSinCos";

	// definition of the name of the operator
	ostringstream name;
	name << "FixSinCos_" << w ;
	setName(name.str());

	setCopyrightString("Florent de Dinechin, Antoine Martinet, Guillaume Sergent, (2013)");

	// everybody needs many digits of Pi
	mpfr_init2(constPi, 10*w);
	mpfr_const_pi( constPi, GMP_RNDN);

	//compute the scale factor		
	mpfr_init2(scale, w+1);
	mpfr_set_d(scale, -1.0, GMP_RNDN);           // exact
	mpfr_mul_2si(scale, scale, -w, GMP_RNDN); // exact
	mpfr_add_d(scale, scale, 1.0, GMP_RNDN);     // exact
	//REPORT(DEBUG, "scale=" << printMPFR(scale, 15));

	// declaring inputs
	addInput("X", w+1);

	// declaring outputs
	addOutput("S", w+1);
	addOutput("C", w+1);

	// These are borders between small-precision cases for which we generate simpler architectures 

	// plain tabulation fits LUTs
	bool wSmallerThanBorder1 = ( w <= target->lutInputs() );    

	//  table with quadrant reduction fits LUTs
	bool wSmallerThanBorder2 = ( w-2 <= target->lutInputs() );   

	// plain tabulation fits BlockRAM
	bool wSmallerThanBorder3 = ( (w+1)*(mpz_class(2)<<(w+1)) <= target->sizeOfMemoryBlock() );

	//  table with quadrant reduction fits BlockRAM
	bool wSmallerThanBorder4 = ( (w+1)*(mpz_class(2)<<(w+1-3)) <= target->sizeOfMemoryBlock() );


	// For all the other methods we need to address a table of sincos with wA bits of the input
	// Let us compute wA such that these bits fit in a blockRAM
	int wAtemp=3;
	int gMax=4; // max of g for all paths
	while((mpz_class(2)<<wAtemp)*(w+1+gMax) < target->sizeOfMemoryBlock()) wAtemp++; 
	int wA = wAtemp--;
	REPORT(DETAILED, "Bits addressing the table for this size : " << wA);  


	// order-1 architecture: sinZ \approx Z, cosZ \approx 1
	int gOrder1Arch=3; // TODO check
	// Y will be smaller than 1/4 => Z will be smaller than pi*2^(-wA-2), or Z<2^-wA 
	// for sinZ we neglect something in the order of Z^2/2 : we want 2^(-2wA-1) < 2^(-w-g)    2wA+1>w+g
	bool wSmallerThanBorderFirstOrderTaylor = w < 2*wA-1-gOrder1Arch;

	// order-2 architecture: sinZ \approx Z, cosZ \approx 1-Z^2/2
	// now we neglect something in the order of Z^3/6 (smaller than Z^3/4): we want 2^(-3wA-2) < 2^(-w-g)    2wA+1>w+g
	int gOrder2Arch=3; // TODO check
	bool wSmallerThanBorderSecondOrderTaylor = w < 3*wA-2-gOrder2Arch;

	REPORT(0*DEBUG, "Boundaries on the various cases for w="<<w << ": " << wSmallerThanBorder1 << wSmallerThanBorder2 << wSmallerThanBorder3 << wSmallerThanBorder4 << wSmallerThanBorderFirstOrderTaylor<< wSmallerThanBorderSecondOrderTaylor);

	if (wSmallerThanBorder1 || (!wSmallerThanBorder2 && wSmallerThanBorder3))
	{
		REPORT(DETAILED, "Simpler architecture: Using plain table" );
		scT= new SinCosTable(target, w+1, 2*(w+1), 1, this);

		addSubComponent(scT); // adding the table

		//int sinCosSize = 2*(w_+1); // size of output (sine plus cosine in a same number, sine on high weight bits)
		vhdl << tab << declare("sinCosTabIn", w+1) << " <= X;" << endl;// signal declaration

		inPortMap(scT, "X", "sinCosTabIn");
		outPortMap(scT, "Y", "SC"); // ports mapping

		vhdl << instance(scT, "sinCosTable" ); //instanciation

		vhdl << tab << declare("Sine", w+1) << " <= SC" << range(2*(w+1)-1, w+1) << ";" << endl;// signal declaration
		vhdl << tab << declare("Cosine", w+1) << " <= SC" << range(w, 0) << ";" << endl;// signal declaration
		//delete (scT);
		vhdl << tab<< "S <= Sine;"<<endl;
		vhdl << tab<< "C <= Cosine;"<<endl;

		//REPORT(0, " scT deleted successfully");
	}

	else if (wSmallerThanBorder2 || (!wSmallerThanBorder3 && wSmallerThanBorder4)) 	{
		REPORT(DETAILED, "Simpler architecture: Using plain table with quadrant reduction");
		/*********************************** RANGE REDUCTION **************************************/
		// the argument is reduced into (0,1/2) because one sin/cos
		// computation in this range can always compute the right sin/cos
		// 2's complement: 0 is always positive
		vhdl << tab << declare ("X_sgn") << " <= X" << of (w) << ";" << endl;
		vhdl << tab << declare ("Q") << " <= X" << of (w-1) << ";" << endl;
		vhdl << tab << declare ("Y_remain",w-1) << " <= X " << range (w-2,0) << ";" << endl;

		vhdl << tab << declare ("Exch") << " <= Q;" << endl;

		vhdl << tab << declare("C_sgn")
			<< " <= Q xor X_sgn;" << endl; //sign of cosin

		vhdl << tab << declare ("Y_in",w-1) << " <= Y_remain;" << endl;

		/*********************************** REDUCED TABLE **************************************/
		scT= new SinCosTable(target, w-1, 2*(w+1), 4, this);

		addSubComponent(scT); // adding the table to vhdl component list

		//int sinCosSize = 2*(w-2); // size of output (sine plus cosine in a same number, sine on high weight bits)
		vhdl << tab << declare("sinCosRedTabIn", w-1) << " <=  Y_in;" << endl;// signal declaration

		inPortMap(scT, "X", "sinCosRedTabIn");
		outPortMap(scT, "Y", "SC_red"); // ports mapping

		vhdl << instance(scT, "sinCosRedTable" ); //instanciation
		//vhdl << instance(scT, "cosTable" );

		vhdl << tab << declare("S_out", w+1) << " <= SC_red " << range( 2*(w+1)-1 , w+1 ) << ";" << endl;// signal declaration
		vhdl << tab << declare("C_out", w+1) << " <= SC_red " << range( w, 0 ) << ";" << endl;// signal declaration
		//delete (scT);

		/*********************************** Reconstruction of both sine and cosine **************************************/

		vhdl << tab << declare ("S_wo_sgn", w+1)
			<< " <= C_out when Exch = '1' else S_out;" << endl; //swap sine and cosine if q xor o.
		vhdl << tab << declare ("C_wo_sgn", w+1)
			<< " <= S_out when Exch = '1' else C_out;" << endl;


		vhdl << tab << declare ("S_wo_sgn_ext", w+1)
			<< " <= S_wo_sgn;" << endl
			<< tab << declare ("C_wo_sgn_ext", w+1)
			<< " <= C_wo_sgn;" << endl; //S_wo_sgn_ext and C_wo_sgn are the positive versions of the sine and cosine

		vhdl << tab << declare ("S_wo_sgn_neg", w+1)
			<< " <= (not S_wo_sgn_ext) + 1;" << endl;
		vhdl << tab << declare ("C_wo_sgn_neg", w+1)
			<< " <= (not C_wo_sgn_ext) + 1;" << endl; //2's complement. We have now the negative version of the sine and cosine results


		vhdl << tab << "S <= S_wo_sgn_ext when X_sgn = '0'"
			<< " else S_wo_sgn_neg;" << endl
			<< tab << "C <= C_wo_sgn_ext when C_sgn = '0'"
			<< " else C_wo_sgn_neg;" << endl; //set the correspondant value to C and S according to input sign


	}




	else { // From now on we will have an argument reduction
		/*********************************** RANGE REDUCTION **************************************/ 
		addComment("The argument is reduced into (0,1/4)");
		vhdl << tab << declare ("X_sgn") << " <= X" << of (w) << ";  -- sign" << endl;
		vhdl << tab << declare ("Q") << " <= X" << of (w-1) << ";  -- quadrant" << endl;
		vhdl << tab << declare ("O") << " <= X" << of (w-2) << ";  -- octant" << endl;
		vhdl << tab << declare ("Y",w-2) << " <= X " << range (w-3,0) << ";" << endl;

		// now X -> X_sgn + Q*.5 + O*.25 + Y where Q,O \in {0,1} and Y \in {0,.25}

		// notY = .25 - y
		//  ?
		// It saves one carry-propagate latency (one cycle)
		// but enlarges the constant multiplier input by g bits
		// and adds one ulp of error
		addComment("Computing .25-Y :  we do a logic NOT, at a cost of 1 ulp");
		manageCriticalPath(target->localWireDelay(w-2) + target->lutDelay());
		vhdl << tab << declare ("notY", w-2) << " <= " << "not Y;" << endl;





		if (wSmallerThanBorderFirstOrderTaylor) {
				REPORT(DETAILED,"Simpler architecture: Using only first order Taylor");

				// now we know which value for g
				g = gOrder1Arch; 
				int wYIn=w-2+g;
				
				vhdl << tab << declare ("Y_in", wYIn) << " <= notY & "
						 << '"' << std::string (g, '1') << '"' << " when O='1' else Y & "
						 << '"' << std::string (g, '0') << '"' << ";"
						 << endl;

				int wY = wYIn-wA; // size of Y_red

				vhdl << tab << declare ( "A", wA) << " <= Y_in " << range(wYIn-1, wYIn-wA) << ";" << endl;
				vhdl << tab << declare ("Y_red", wY) << " <= Y_in" << range (wYIn-wA-1,0) << ";" << endl;

				//------------------------------------SinCosTable building for A -------------------------------------
				scT= new SinCosTable(target, wA, 2*(w+1+g), 8, this);
				addSubComponent(scT); // adding the table to vhdl component list
				inPortMap(scT, "X", "A");
				outPortMap(scT, "Y", "SCA");
				vhdl << instance(scT, "sinCosATable" ); 

				vhdl << tab << declare("SA", w+1+g) << " <= SCA " << range( 2*(w+1+g)-1 , w+1+g ) << ";" << endl;
				vhdl << tab << declare("CA", w+1+g) << " <= SCA " << range( w+g, 0 ) << ";" << endl;

				//-------------------------------- MULTIPLIER BY PI ---------------------------------------

				map<string, double> pi_mult_inputDelays;
				pi_mult_inputDelays["X"] = getCriticalPath();
				int wZ=w-wA+g;

				// the 2 extra bits to Z are added by the KCM multiplier
				pi_mult = new FixRealKCM (target,-w-g,-2-wA-1, false, -w-g , "pi", 1.0, pi_mult_inputDelays); //lsbOut=-7 because we add 6 bits between entrance and exit. 
				oplist.push_back (pi_mult);
				inPortMap (pi_mult, "X", "Y_red");
				outPortMap (pi_mult, "R", "Z");
				int wZz=getSignalByName("Z")->width();
				REPORT(DEBUG, "wZ=" <<wZ<<";"<<" wZz="<<wZz<<";");
				vhdl << instance (pi_mult, "pi_mult");

				//----------------- 1-Y_red ( = cos(Y_red) ) -----------------

				//		vhdl << declare("Z_op", Y_red_sqr_size) << "<= '0' & "<< og(Y_red_sqr_size-1)<< " - ( '0' & Y_red_sqr" << range(Y_red_sqr_size-1,1) << " );" << endl;
				int m=4;
				//---------------------------- Sine computation ------------------------
				vhdl << tab <<  declare("SinACosY_red",w+1+g) << " <= SA;"<<endl;
				vhdl << tab << declare("CosASinY_red", 2*wZ+m ) << " <= CA" << range( w+g, w+1+g-wZ-m ) <<"*Z;" <<endl;
				vhdl << tab << declare("PreSinX", w+1+g) << " <= SinACosY_red + ( " << zg(wA-2+1) << " & CosASinY_red" << range( 2*wZ+m-1, wZ+m-1 ) << " );"<<endl;

				//---------------------------- Cosine computation -------------------------------
				vhdl << tab << declare("CosACosY_red", w+1+g ) << " <= CA;" << endl;
				vhdl << tab << declare("SinASinY_red", 2*wZ+m ) << " <= SA" << range( w+g, w+1+g-wZ-m ) << "*Z;" << endl;
				vhdl << tab << declare("PreCosX", w+1+g) << " <= CosACosY_red - ( " << zg(wA-2+1) << " & SinASinY_red" << range( 2*wZ+m-1, wZ+m-1 )<< " );" << endl;

			
				// Something suspiscious here, we do not use the last guard bit.
				vhdl << tab << declare ("C_out", w) << " <= PreCosX" << range (w+g, g+1) << ';' << endl;
				vhdl << tab << declare ("S_out", w) << " <= PreSinX" << range (w+g, g+1) << ';' << endl;
			}



		else if (wSmallerThanBorderSecondOrderTaylor) {

		REPORT(DETAILED,"Using first-order Taylor for sine and second-order for cosine");
		g=gOrder2Arch; // TODO check
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


		int wY = wIn-wA; // size of Y_red
		int wZ = wY+2; // size 

		vhdl << tab << declare ( "A", wA) << " <= Y_in " << range(wIn-1, wIn-wA) << ";" << endl;
		vhdl << tab << declare ("Y_red", wY) << " <= Y_in" << range (wIn-wA-1,0) << ';' << endl;

		//------------------------------------SinCosTable building for A -------------------------------------

		/*********************************** REDUCED TABLE **************************************/
		scT= new SinCosTable(target, wA, w+ g, 8, this);

		addSubComponent(scT); // adding the table to vhdl component list

		//int sinCosSize = 2*(wA-2); // size of output (sine plus cosine in a same number, sine on high weight bits)
		vhdl << tab << declare("sinCosATabIn", wA) << " <= A;" << endl;// signal declaration

		inPortMap(scT, "X", "sinCosATabIn");
		outPortMap(scT, "Y", "SCA"); // ports mapping

		vhdl << instance(scT, "sinCosATable" ); //instanciation
		//vhdl << instance(scT, "cosTable" );

		vhdl << tab << declare("SA", wA+1) << " <= SCA " << range( 2*(wA+1)-1 , wA+1 ) << ";" << endl;// signal declaration
		vhdl << tab << declare("CA", wA+1) << " <= SCA " << range( wA, 0 ) << ";" << endl;// signal declaration
		//delete (scT);

		//-------------------------------- MULTIPLIER BY PI ---------------------------------------

		// For w=32: latency=18 (1 for the const mult), 781 + 1170	
		map<string, double> pi_mult_inputDelays;
		pi_mult_inputDelays["X"] = getCriticalPath();


		// the 2 extra bits to Z are added by the KCM multiplier
		pi_mult = new FixRealKCM (target, 0, wY-1, false, 0, "pi", 1.0, pi_mult_inputDelays); 
		oplist.push_back (pi_mult);
		inPortMap (pi_mult, "X", "Y_red");
		outPortMap (pi_mult, "R", "Z");
		vhdl << instance (pi_mult, "pi_mult");


		//--------------------------- SQUARER --------------------------------
		const int Y_red_sqr_size = 2*(wY+2);
		vhdl << declare("Y_red_sqr", Y_red_sqr_size ) << "<= Z*Z;" << endl;;


		//----------------- 1-Y_red ( = cos(Y_red) ) -----------------

		vhdl << declare("Z_op", Y_red_sqr_size) << "<= '0' & "<< og(Y_red_sqr_size-1)<< " - ( '0' & Y_red_sqr" << range(Y_red_sqr_size-1,1) << " );" << endl;

		//---------------------------- Sine computation ------------------------
		vhdl << declare("SinACosY_red",Y_red_sqr_size+wA+1) << " <= SA*Z_op;"<<endl;
		vhdl << declare("CosASinY_red", wY+wA+1) << " <= CA*Y_red;" <<endl;
		vhdl << declare("PreSinX",Y_red_sqr_size+wA+1) << " <= SinACosY_red + CosASinY_red;"<<endl;

		//---------------------------- Cosine computation ------------------------------------
		vhdl << declare("CosACosY_red", Y_red_sqr_size+wA+1 ) << " <= CA*Z_op;" << endl;
		vhdl << declare("SinASinY_red", wY+wA+1) << " <= SA*Y_red;" << endl;
		vhdl << declare("PreCosX", Y_red_sqr_size+wA+1) << " <= CosACosY_red - SinASinY_red;" << endl;

		// TODO sizes probably won't match
#if 0
		//--------------------------- Reconstruction of both sine and cosine -----------------------------

		vhdl << tab << declare ("S_wo_sgn", Y_red_sqr_size+wA+1)
			<< " <= PreCosX when Exch = '1' else PreSinX;" << endl; //swap sine and cosine if q xor o.
		vhdl << tab << declare ("C_wo_sgn",Y_red_sqr_size+wA+1)
			<< " <= PreSinX when Exch = '1' else PreCosX;" << endl;


		vhdl << tab << declare ("S_wo_sgn_ext", Y_red_sqr_size+wA+1)
			<< " <= S_wo_sgn;" << endl
			<< tab << declare ("C_wo_sgn_ext", Y_red_sqr_size+wA+1)
			<< " <= C_wo_sgn;" << endl; //S_wo_sgn_ext and C_wo_sgn are the positive versions of the sine and cosine

		vhdl << tab << declare ("S_wo_sgn_neg", Y_red_sqr_size+wA+1)
			<< " <= (not S_wo_sgn_ext) + 1;" << endl;
		vhdl << tab << declare ("C_wo_sgn_neg", Y_red_sqr_size+wA+1)
			<< " <= (not C_wo_sgn_ext) + 1;" << endl; //2's complement. We have now the negative version of the sine and cosine results


		vhdl << tab << "S <= S_wo_sgn_ext"<< range(Y_red_sqr_size+wA,Y_red_sqr_size+wA-w) << " when X_sgn = '0'"
			<< " else S_wo_sgn_neg"<< range(Y_red_sqr_size+wA,Y_red_sqr_size+wA-w) << ";" << endl
			<< tab << "C <= C_wo_sgn_ext "<< range(Y_red_sqr_size+wA,Y_red_sqr_size+wA-w) << " when C_sgn = '0'"
			<< " else C_wo_sgn_neg"<< range(Y_red_sqr_size+wA,Y_red_sqr_size+wA-w) << ";" << endl; //set the correspondant value to C and S according to input sign
#endif
		}




		else	{
			REPORT(DETAILED,"Generic case: Using third-order Taylor");

			// we need to know the number of guard bits _now_ to have a good Y_in
			// some precision is lost with the optimized multipliers
			g=4;
			//const int g=4; //guard bits
			//const int g=3; //guard bits
			//we take 4 guard bits even if error < 8 ulp because rounding will take
			//another half final ulp
			int wIn = w-2+g;

			vhdl << tab << declare ("Y_in",wIn) << " <= notY & "
					 << '"' << std::string (g, '1') << '"' << " when O='1' else Y & "
					 << '"' << std::string (g, '0') << '"' << ";"
					 << endl;


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
		  wA = (int) ceil ((double) (w+g+2)/4.) - 2;
			if (wA <= 3)
				wA = 3;



			/*********************************** THE TABLES **************************************/
			
#if 0
			// FIXME: assumes that both
			// -the memory width is a power of 2
			// -the bram can always be used as dual-port (one for sin, one for cos)
			//  with width >=w+g
			// upper log2 of memory width needed for table
			int wg_log2 = (int) ceil ((double) log2 (w+g));
			int bram_words = target->sizeOfMemoryBlock() / (1u << wg_log2);
			int bram_words_log2 = (int) floor ((double) log2 (bram_words));
			if (wA <= bram_words_log2-1)
				wA = bram_words_log2-1;
#endif
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

			//	int wZ3 = 3*wZ - 2*(w+g) -1; // -1 for the div by 2
			int wZ3o6 = 3*wZ - 2*(w+g) -2;
			if (wZ3o6 < 2)
				wZ3o6 = 2; //using 1 will generate bad vhdl
			
			if(wZ3o6<=12) {
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






		
#if 1 		// No bit heap
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
			nextCycle();
			
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
			
		REPORT(INFO, " wA=" << wA <<" wZ=" << wZ <<" wZ2=" << wZ2o2 <<" wZ3o6=" << wZ3o6 );

		// For LateX in the paper
		//	cout << "     " << w <<  "   &   "  << wA << "   &   " << wZ << "   &   " << wZ2o2 << "   &   " << wZ3o6 << "   \\\\ \n \\hline" <<  endl;

		} // closes if for generic case



		// When we arrive here we should have two signals C_out and S_out, each of size w 
		addComment("--- Final reconstruction of both sine and cosine ---");

		vhdl << tab << declare("C_sgn") << " <= X_sgn xor Q;" << endl;
		vhdl << tab << declare ("Exch") << " <= Q xor O;" << endl;

		vhdl << tab << declare ("S_wo_sgn", w)
				 << " <= C_out when Exch = '1' else S_out;" << endl; //swap sine and cosine if q xor o
		vhdl << tab << declare ("C_wo_sgn", w)
				 << " <= S_out when Exch = '1' else C_out;" << endl;
		
		
		vhdl << tab << declare ("S_wo_sgn_ext", w+1)
				 << " <= '0' & S_wo_sgn;" << endl
				 << tab << declare ("C_wo_sgn_ext", w+1)
				 << " <= '0' & C_wo_sgn;" << endl; //S_wo_sgn_ext and C_wo_sgn are the positive versions of the sine and cosine
		
		vhdl << tab << declare ("S_wo_sgn_neg", w+1)
				 << " <= (not S_wo_sgn_ext) + 1;" << endl;
		vhdl << tab << declare ("C_wo_sgn_neg", w+1)
				 << " <= (not C_wo_sgn_ext) + 1;" << endl; //2's complement. We have now the negative version of the sine and cosine results
		
		vhdl << tab << "S <= S_wo_sgn_ext when X_sgn = '0'"
				 << " else S_wo_sgn_neg;" << endl;

		vhdl << tab << "C <= C_wo_sgn_ext when C_sgn = '0'"
				 << " else C_wo_sgn_neg;" << endl; //set the correspondant value to C and S according to input sign
	}

};




FixSinCos::~FixSinCos(){
	if(scT) free(scT);
	if(pi_mult) free(pi_mult);
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

	//REPORT(DEBUG, "Emulate scale=" << printMPFR(scale, 15));
	
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

