
#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicAtan2.hpp"
#include "IntMult/IntMultiplier.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "ShiftersEtc/LZOC.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "FixFunctions/FixFunctionByPiecewisePoly.hpp"
#include "FixFunctions/BipartiteTable.hpp"


using namespace std;
// TODO reduce the CORDIC datapaths (Y with leading 0s and X with leading ones)
// Or write an exact version of the XY datapath

/* TODO Debugging:
There are still a few last-bit errors with the current setup in
./flopoco -verbose=3 -pipeline=no CordicAtan2 8 TestBenchFile -2


They are all for vectors of small norm 

One is unavoidable, it is atan2 (0,0)
   TODO add a don't care to the test framework?
*/
namespace flopoco{

	// TODO Options:
	// An option to output the angle in the (atan(2^-i)) basis
	// an an option for outputting the norm of the vector as well (scaled or not)

#define CORDIC 8
#define CORDIC_SCALING 9
#define INVMULTATAN 0















	CordicAtan2::CordicAtan2(Target* target, int w_, int method, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
	
		int stage;
		srcFileName="CordicAtan2";
		setCopyrightString ( "Matei Istoan, Florent de Dinechin (2012-...)" );
		useNumericStd_Unsigned();

		ostringstream name;
		name << "CordicAtan2_"<< w_ << "_uid" << getNewUId();
		setNameWithFreq( name.str() );
	
		// everybody needs many digits of Pi (used by emulate etc)
		mpfr_init2(constPi, 10*w);
		mpfr_const_pi( constPi, GMP_RNDN);
	
		mpfr_t  zatan;	// used only by CORDIC but who cares
		mpfr_init2(zatan, 10*w);
	
	
		// declaring inputs. 
		// man atan2 says "atan2(y,x) is atan(y/x)" so we will provide the inputs in the same order...
		addInput  ( "Y"  , w, true );
		addInput  ( "X"  , w, true );
	
		// declaring output
		addOutput  ( "A"  , w, 2 );
	
		setCriticalPath(getMaxInputDelays(inputDelays));
		//		manageCriticalPath( target->lutDelay());
	
		int degree=method & 7;
		method=method & 0xF8;
		REPORT(DEBUG, "method=" << method << "  degree=" << degree);
	
		//Defining the various parameters according to method

		// When pipelining I noticed that we perform a subtraction for the comparison X<Y in parallel to the negation anyway
		// so negateByComplement should always be false, 
		//except if some day we do an ASIC target where it will cost less hardware.
		int sizeXYInRR;
		bool doScaling;
		switch(method) {
		case CORDIC:
			negateByComplement=false;
			sizeXYInRR = w;
			computeGuardBitsForCORDIC();
			doScaling=false;
			break;
		case CORDIC_SCALING:
			negateByComplement=false; // it would entail a larger LZC and larger shifter.
			computeGuardBitsForCORDIC();
			sizeXYInRR = w;
			doScaling=true;
			break;
		case INVMULTATAN:
			negateByComplement=false;
			gXY=0;
			gA=0;
			sizeXYInRR = w;
			doScaling=true;
			break;
		}

		///////////// VHDL GENERATION
	
		/////////////////////////////////////////////////////////////////////////////
		// 
		//    First range reduction
		//
		/////////////////////////////////////////////////////////////////////////////
		setCriticalPath(getMaxInputDelays(inputDelays));

		vhdl << tab << declare("sgnX") << " <= X" << of(w-1) << ";" << endl;
		vhdl << tab << declare("sgnY") << " <= Y" << of(w-1) << ";" << endl;
	

		manageCriticalPath( target->adderDelay(w));
		if (negateByComplement)	{
			vhdl << tab << declare("pX", sizeXYInRR) << " <=      X  & " << zg(gXY)<< ";" << endl;
			vhdl << tab << declare("pY", sizeXYInRR) << " <=      Y  & " << zg(gXY)<< ";" << endl;
			vhdl << tab << declare("mX", sizeXYInRR) << " <= (not X) & " << og(gXY)<< ";  -- negation by not, implies one ulp error." << endl;
			vhdl << tab << declare("mY", sizeXYInRR) << " <= (not Y) & " << og(gXY)<< ";  -- negation by not, implies one ulp error. " << endl;
		}else {
			vhdl << tab << declare("pX", sizeXYInRR) << " <= X;" << endl;
			vhdl << tab << declare("pY", sizeXYInRR) << " <= Y;" << endl;
			vhdl << tab << declare("mX", sizeXYInRR) << " <= (" << zg(w) << " - X);" << endl;
			vhdl << tab << declare("mY", sizeXYInRR) << " <= (" << zg(w) << " - Y);" << endl;
		}

		// TODO: replace the following with LUT-based comparators
		// and first maybe experiment with synthesis tools
		vhdl << tab << declare("XmY", w+1) << " <= (sgnX & X)-(sgnY & Y);" << endl;
		vhdl << tab << declare("XpY", w+1) << " <= (sgnX & X)+(sgnY & Y);" << endl;
		vhdl << tab << declare("XltY") << " <= XmY" << of(w) <<";" << endl;
		vhdl << tab << declare("mYltX") << " <= not XpY" << of(w) <<";" << endl;
		// Range reduction: we define 4 quadrants, each centered on one axis (these are not just the sign quadrants)
		// Then each quadrant is decomposed in its positive and its negative octant.

		manageCriticalPath( target->lutDelay() + target->localWireDelay());

		vhdl << tab << "-- quadrant will also be the angle to add at the end" <<endl;
		vhdl << tab << declare("quadrant", 2) << " <= " << endl;
		vhdl << tab << tab << "\"00\"  when (not sgnX and not XltY and     mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"01\"  when (not sgnY and     XltY and     mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"10\"  when (    sgnX and     XltY and not mYltX)='1' else"    << endl;
		vhdl << tab << tab << "\"11\";"    << endl;
	


		int sizeXYR = sizeXYInRR-1; // no need for sign bit any longer
		vhdl << tab << declare("XR", sizeXYR) << " <= " << endl;
		vhdl << tab << tab << "pX" << range(sizeXYR-1, 0) << " when quadrant=\"00\"   else " << endl;
		vhdl << tab << tab << "pY" << range(sizeXYR-1, 0) << " when quadrant=\"01\"   else " << endl;
		vhdl << tab << tab << "mX" << range(sizeXYR-1, 0) << " when quadrant=\"10\"   else " << endl;
		vhdl << tab << tab << "mY" << range(sizeXYR-1, 0) << ";"    << endl;
	
		vhdl << tab << declare("YR", sizeXYR) << " <= " << endl;
		vhdl << tab << tab << "pY" << range(sizeXYR-1, 0) << " when quadrant=\"00\" and sgnY='0'  else " << endl;
		vhdl << tab << tab << "mY" << range(sizeXYR-1, 0) << " when quadrant=\"00\" and sgnY='1'  else " << endl;
		vhdl << tab << tab << "pX" << range(sizeXYR-1, 0) << " when quadrant=\"01\" and sgnX='0'  else " << endl;
		vhdl << tab << tab << "mX" << range(sizeXYR-1, 0) << " when quadrant=\"01\" and sgnX='1'  else " << endl;
		vhdl << tab << tab << "pY" << range(sizeXYR-1, 0) << " when quadrant=\"10\" and sgnY='0'  else " << endl;
		vhdl << tab << tab << "mY" << range(sizeXYR-1, 0) << " when quadrant=\"10\" and sgnY='1'  else " << endl;
		vhdl << tab << tab << "pX" << range(sizeXYR-1, 0) << " when quadrant=\"11\" and sgnX='0'  else "    << endl;
		vhdl << tab << tab << "mX" << range(sizeXYR-1, 0) << " ;"    << endl;

		vhdl << tab << declare("finalAdd") << " <= " << endl;
		vhdl << tab << tab << "'1' when (quadrant=\"00\" and sgnY='0') or(quadrant=\"01\" and sgnX='1') or (quadrant=\"10\" and sgnY='1') or (quadrant=\"11\" and sgnX='0')" << endl;
		vhdl << tab << tab << " else '0';  -- this information is sent to the end of the pipeline, better compute it here as one bit"    << endl; 


		////////////////////////////////////////////////////////////////////////////
		// 
		//    Second range reduction, scaling 
		//
		////////////////////////////////////////////////////////////////////////////

		if(doScaling) {
			manageCriticalPath( target->lutDelay() + target->localWireDelay());
			vhdl << tab << declare("XorY", sizeXYR-1) << " <= XR" << range(sizeXYR-1,1) << " or YR" << range(sizeXYR-1,1) << ";" << endl;
			// The LZC
			LZOC* lzc = new	LZOC(target, sizeXYR-1);
			addSubComponent(lzc);
			
			inPortMap(lzc, "I", "XorY");
			inPortMapCst(lzc, "OZB", "'0'");
			outPortMap(lzc, "O", "S"); 
			vhdl << instance(lzc, "lzc");
			syncCycleFromSignal("S");
			nextCycle();

			//setCycleFromSignal("lzo");
			//		setCriticalPath( lzc->getOutputDelay("O") );
			
			// The two shifters are two instance of the same component
			Shifter* lshift = new Shifter(target, sizeXYR, sizeXYR-1, Shifter::Left);   
			addSubComponent(lshift);
			
			inPortMap(lshift, "S", "S");
			
			inPortMap(lshift, "X", "XR");
			outPortMap(lshift, "R", "XRSfull");
			vhdl << instance(lshift, "Xshift");
			vhdl << tab << declare("XRS", sizeXYR) << " <=  XRSfull " << range(sizeXYR-1,0) << ";" << endl;

			inPortMap(lshift, "X", "YR");
			outPortMap(lshift, "R", "YRSfull");
			vhdl << instance(lshift, "Yshift");
			vhdl << tab << declare("YRS", sizeXYR) << " <=  YRSfull " << range(sizeXYR-1,0) << ";" << endl;

			syncCycleFromSignal("YRSfull");
			nextCycle();
			
			//syncCycleFromSignal("small_absZ0_normd_full");
			//setCriticalPath( getOutputDelay("R") );
			//double cpsmall_absZ0_normd = getCriticalPath();
			
			//int  = getSignalByName("small_absZ0_normd_full")->width() - (wF-pfinal+2);
		}
		else{ //scalingRR
			vhdl << tab << declare("XRS", sizeXYR) << " <=  XR;" << endl;
			vhdl << tab << declare("YRS", sizeXYR) << " <=  YR;" << endl;
		}







		int sizeZ=w-2+gA; // w-2 because two bits come from arg red 

		string finalZ; // to hold the atan before reconstruction


		if(method==CORDIC || method==CORDIC_SCALING) {

			////////////////////////////////////////////////////////////////////////////
			// 
			//   	 CORDIC iterations	
			//
			////////////////////////////////////////////////////////////////////////////

			// Fixed-point considerations:
			// Y -> 0 and X -> K.sqrt(x1^2+y1^2)
			// Max value attained by X is sqrt(2)*K which is smaller than 2

			int zMSB=-1;      // -1 because these two bits have weight 0 and -1, but we must keep the sign
			int zLSB = zMSB-sizeZ+1;
			int sizeX = w+gXY;
			int sizeY = sizeX;

			//		REPORT(DEBUG, "sizeXY=" << sizeXY << "   sizeXYR=" << sizeXYR);

			vhdl << tab << declare("X1", sizeX) << " <= '0' & XRS & " << zg(sizeX-sizeXYR-1) << ";" <<endl;			
			vhdl << tab << declare("Y1", sizeY) << " <= '0' & YRS & " << zg(sizeY-sizeXYR-1) << ";" <<endl;			
			stage=1; 
 
			manageCriticalPath( target->adderDelay(sizeX));
			vhdl << tab << "--- Iteration " << stage << " : sign is known positive ---" << endl;
			vhdl << tab << declare(join("YShift", stage), sizeX) << " <= " << rangeAssign(sizeX-1, sizeX-stage, "'0'") << " & Y" << stage << range(sizeX-1, stage) << ";" << endl;
			vhdl << tab << declare(join("X", stage+1), sizeX) << " <= " 
					 << join("X", stage) << " + " << join("YShift", stage) << " ;" << endl;
		
			vhdl << tab << declare(join("XShift", stage), sizeY) << " <= " << zg(stage) << " & X" << stage << range(sizeY-1, stage) << ";" <<endl;			
			vhdl << tab << declare(join("Y", stage+1), sizeY) << " <= " 
					 << join("Y", stage) << " - " << join("XShift", stage) << " ;" << endl;
			

			//create the constant signal for the arctan
			mpfr_set_d(zatan, 1.0, GMP_RNDN);
			mpfr_div_2si(zatan, zatan, stage, GMP_RNDN);
			mpfr_atan(zatan, zatan, GMP_RNDN);
			mpfr_div(zatan, zatan, constPi, GMP_RNDN);
			mpfr_t roundbit;
			mpfr_init2(roundbit, 30); // should be enough for anybody
			mpfr_set_d(roundbit, 1.0, GMP_RNDN);
			mpfr_div_2si(roundbit, roundbit, w, GMP_RNDN); // roundbit is in position 2^-w
		
			REPORT(DEBUG, "stage=" << stage << "  atancst=" << printMPFR(zatan));

			mpfr_add(zatan, zatan, roundbit, GMP_RNDN);
			vhdl << tab << declare("Z2", sizeZ) << " <= " << unsignedFixPointNumber(zatan, zMSB, zLSB) << "; -- initial atan, plus round bit" <<endl;


			for(stage=2; stage<=maxIterations;    stage++, sizeY--){
				// Invariant: sizeX-sizeY = stage-2
				vhdl << tab << "--- Iteration " << stage << " ---" << endl;

				manageCriticalPath( target->localWireDelay(sizeX+1) + target->adderDelay(max(sizeX,sizeZ)) );
				vhdl << tab << declare(join("sgnY", stage))  << " <= " <<  join("Y", stage)  <<  of(sizeY-1) << ";" << endl;
			
				if(-2*stage+1 >= -w+1-gXY) { 
					vhdl << tab << declare(join("YShift", stage), sizeX) 
							 << " <= " << rangeAssign(sizeX-1, sizeX -(sizeX-sizeY+stage), join("sgnY", stage))   
							 << " & Y" << stage << range(sizeY-1, stage) << ";" << endl;
				
					vhdl << tab << declare(join("X", stage+1), sizeX) << " <= " 
							 << join("X", stage) << " - " << join("YShift", stage) << " when " << join("sgnY", stage) << "=\'1\'     else "
							 << join("X", stage) << " + " << join("YShift", stage) << " ;" << endl;
				}
				else {	// autant pisser dans un violon
					vhdl << tab << declare(join("X", stage+1), sizeX) << " <= " << join("X", stage) << " ;" << endl;
				}

				vhdl << tab << declare(join("XShift", stage), sizeY) << " <= " << zg(2) << " & X" << stage << range(sizeX-1, sizeX - sizeY + 2) << ";" <<endl;			
				vhdl << tab << declare(join("YY", stage+1), sizeY) << " <= " 
						 << join("Y", stage) << " + " << join("XShift", stage) << " when " << join("sgnY", stage) << "=\'1\'     else "
						 << join("Y", stage) << " - " << join("XShift", stage) << " ;" << endl;
				vhdl << tab << declare(join("Y", stage+1), sizeY-1) << " <= " << join("YY", stage+1) << range(sizeY-2, 0) << ";" <<endl;
			

				//create the constant signal for the arctan
				mpfr_set_d(zatan, 1.0, GMP_RNDN);
				mpfr_div_2si(zatan, zatan, stage, GMP_RNDN);
				mpfr_atan(zatan, zatan, GMP_RNDN);
				mpfr_div(zatan, zatan, constPi, GMP_RNDN);
				REPORT(DEBUG, "stage=" << stage << "  atancst=" << printMPFR(zatan));		
				// rounding here in unsignedFixPointNumber()
				vhdl << tab << declare(join("atan2PowStage", stage), sizeZ) << " <= " << unsignedFixPointNumber(zatan, zMSB, zLSB) << ";" <<endl;
				vhdl << tab << declare(join("Z", stage+1), sizeZ) << " <= " 
						 << join("Z", stage) << " + " << join("atan2PowStage", stage) << " when " << join("sgnY", stage) << "=\'0\'      else "
						 << join("Z", stage) << " - " << join("atan2PowStage", stage) << " ;" << endl;

			} //end for loop
		
			// Give the time to finish the last rotation
			// manageCriticalPath( target->localWireDelay(w+1) + target->adderDelay(w+1) // actual CP delay
			//                     - (target->localWireDelay(sizeZ+1) + target->adderDelay(sizeZ+1))); // CP delay that was already added

			manageCriticalPath( target->localWireDelay(w+1) + target->adderDelay(2) );

			vhdl << tab << declare("finalZ", w) << " <= Z" << stage << of(sizeZ-1) << " & Z" << stage << range(sizeZ-1, sizeZ-w+1) << "; -- sign-extended and rounded" << endl;
		}





		else if (method==INVMULTATAN) {
			//FixFunctionByPiecewisePoly* recipTable;
			Operator* recipTable;
			//FixFunctionByPiecewisePoly* atanTable;
			Operator* atanTable;
			int msbRecip, lsbRecip, msbProduct, lsbProduct, msbAtan, lsbAtan;
			msbAtan = -2; // bits 0 and -1 come from the range reduction
			lsbAtan = -w+1;
			msbRecip = 1; // 1/x <= 2
			msbProduct = -1 ; // y/x between 0 and 1 but the faithful product may overflow a bit. 
			if(degree==0){ // both tables are correctly rounded
				lsbRecip = -w+1; // see error analysis in the paper
				lsbProduct = -w+1;
			}
			else{ // both tables are faithful
				lsbRecip = -w; // see error analysis in the paper
				lsbProduct = -w;				
			}
			// recip table computes 2/(1+x) once we have killed the MSB of XRS, which is always 1.
			vhdl << tab << declare("XRm1", w-2) << " <= XRS" << range(w-3,0)  << "; -- removing the MSB which is constantly 1" << endl;
			ostringstream invfun;
			invfun << "2/(1+x)-1b"<<lsbRecip;
			if(degree==1) {
				//BipartiteTable *deltaX2Table
				recipTable = new BipartiteTable(target,
													invfun.str(),
													-w+2,  // XRS was between 1/2 and 1. XRm1 is between 0 and 1/2
													msbRecip + 1, // +1 because internally uses signed arithmetic and we want an unsigned result
													lsbRecip);
			}
			else {
				recipTable = new FixFunctionByPiecewisePoly(target, 
																invfun.str(),
																-w+2,  // XRS was between 1/2 and 1. XRm1 is between 0 and 1/2
																msbRecip + 1, // +1 because internally uses signed arithmetic and we want an unsigned result
																lsbRecip,
																degree,
																true /*finalRounding*/
																);
			}
			recipTable->changeName(join("reciprocal_uid", getNewUId()));
			addSubComponent(recipTable);			
			inPortMap(recipTable, "X", "XRm1");
			outPortMap(recipTable, "Y", "R0"); 
			vhdl << instance(recipTable, "recipTable");
			vhdl << tab << declareFixPoint("R", false, msbRecip, lsbRecip) << " <= unsigned(R0" << range(msbRecip-lsbRecip  , 0) << "); -- removing the sign  bit" << endl;
			vhdl << tab << declareFixPoint("YRU", false, -1, -w+1) << " <= unsigned(YRS);" << endl;

			if(target->plainStupidVHDL()) { // generate a "*"
				vhdl << tab << declareFixPoint("P", false, msbRecip -1 +1, lsbRecip-w+1) << " <= R*YRU;" << endl;
				resizeFixPoint("PtruncU", "P", msbProduct, lsbProduct);
				vhdl << tab << declare("P_slv", msbProduct-lsbProduct+1)  << " <=  std_logic_vector(PTruncU);" << endl;
			}
			else{ // generate an IntMultiplier
				IntMultiplier::newComponentAndInstance(this,
																							 "divMult",     // instance name
																							 "R",  // x
																							 "YRU", // y
																							 "P",       // p
																							 msbProduct, lsbProduct
																							 );
			}


			ostringstream atanfun;
			atanfun << "atan(x)/pi";
			if(degree==1) {
				atanTable = new BipartiteTable(target,
													atanfun.str(),
													lsbProduct,
													msbAtan,
													lsbAtan);
			}
			else {
				atanTable  = new FixFunctionByPiecewisePoly(target,
																atanfun.str(),
																lsbProduct,
																msbAtan,
																lsbAtan,
																degree,
																true /*finalRounding*/
																);
			}
			atanTable->changeName(join("atan_uid", getNewUId()));
			addSubComponent(atanTable);			
			inPortMap(atanTable, "X", "P_slv");
			outPortMap(atanTable, "Y", "atanTableOut"); 
			vhdl << instance(atanTable, "atanTable");

			vhdl << tab << declare("finalZ", w) << " <= \"00\" & atanTableOut;" << endl;
			
			//			vhdl << tab << declare(finalZ, sizeZ) << " <= '0' & finalZu; -- adding back a sign  bit for the reconstruction" << endl;
			 
		}




		////////////////////////////////////////////////////////////////////////////
		// 
		//                            reconstruction
		//
		////////////////////////////////////////////////////////////////////////////

		vhdl << tab << declare("qangle", w) << " <= (quadrant & " << zg(w-2) << ");" << endl;
		vhdl << tab << "A <= "
				 << tab << tab << "     qangle + finalZ  when finalAdd='1'" << endl
				 << tab << tab << "else qangle - finalZ;" << endl;
	};


	CordicAtan2::~CordicAtan2(){
		mpfr_clears (kfactor, constPi, NULL);		
	};




	void CordicAtan2::computeGuardBitsForCORDIC(){	
#define ROUNDED_ROTATION 0 // 0:trunc 
	
#if ROUNDED_ROTATION
			REPORT(DEBUG, "Using rounded rotation trick");
#endif

			// ulp = weight of the LSB of the result is 2^(-w+1)
			// half-ulp is 2^-w
			// 1/pi atan(2^-w) < 1/2. 2^-w therefore after w-1 interations,  method error will be bounded by 2^-w 
			maxIterations = w-1;
			//error analysis for the (x,y) datapath
			double eps;  //error in ulp

			if(negateByComplement)
				eps=1; // initial neg-by-not
			else
				eps=0;

			double shift=0.5;
			for(int stage=1; stage<=maxIterations; stage++){
#if ROUNDED_ROTATION
				eps = eps + eps*shift + 0.5; // 0.5 assume rounding in the rotation.
#else
				eps = eps + eps*shift + 1.0; // 1.0 assume truncation in the rotation.
#endif
				shift *=0.5;
			}
			// guard bits depend only on the number of iterations
			gXY = (int) ceil(log2(eps));  
			//error analysis for the A datapath
			eps = maxIterations*0.5; // only the rounding error in the atan constant
			gA = 1 + (int) ceil(log2(eps)); // +1 for the final rounding 
			REPORT(DEBUG, "Error analysis computes eps=" << eps << " ulps on the XY datapath, hence  gXY=" << gXY);
			REPORT(DEBUG, "Error analysis computes eps=" << eps <<  " ulps on the A datapath, hence  gA=" << gA );
	} 








	void CordicAtan2::emulate(TestCase * tc) 
	{
		mpfr_t x,y,a;
		mpfr_init2(x, 10*w);
		mpfr_init2(y, 10*w);
		mpfr_init2(a, 10*w);

		mpz_class az;

		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		// interpret as signed two'ss complement
		if (1==(svX >> (w-1))) // sign bit
			svX -= (mpz_class(1)<<w);
		if (1==(svY >> (w-1))) // sign bit
			svY -= (mpz_class(1)<<w);
		/* Compute correct value */
		
		mpfr_set_z (x, svX.get_mpz_t(), GMP_RNDN); //  exact
		mpfr_set_z (y, svY.get_mpz_t(), GMP_RNDN); //  exact
	
		mpfr_atan2(a, y, x, GMP_RNDN); // a between -pi and pi
		mpfr_div(a, a, constPi, GMP_RNDN); // a between -1 and 1

		// Now convert a to fix point
		// Align to fix point by adding 6 -- if we just add 4 there is a 1-bit shift in case a<0
		mpfr_add_d(a, a, 6.0, GMP_RNDN);
		mpfr_mul_2si (a, a, w-1, GMP_RNDN); // exact scaling 

		mpz_class mask = (mpz_class(1)<<w) -1; 

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDD); // there can be a real rounding here
		az -= mpz_class(6)<<(w-1);
		az &= mask; 
 		tc->addExpectedOutput ("A", az);

		mpfr_get_z (az.get_mpz_t(), a, GMP_RNDU); // there can be a real rounding here
		az -= mpz_class(6)<<(w-1);
		az &= mask; 
 		tc->addExpectedOutput ("A", az);
		
		// clean up
		mpfr_clears (x,y,a, NULL);		
	}






	void CordicAtan2::buildStandardTestCases(TestCaseList * tcl) 
	{
		mpfr_t z;
		mpz_class zz;
		TestCase* tc;

		// 0
		tc = new TestCase (this);
		tc -> addInput ("X", mpz_class(1)<< (w-2));
		tc -> addInput ("Y", mpz_class(0));
		emulate(tc);
		tcl->add(tc);


		//pi/4
		tc = new TestCase (this);
		tc -> addInput ("X", mpz_class(1)<< (w-2));
		tc -> addInput ("Y", mpz_class(1)<< (w-2));
		emulate(tc);
		tcl->add(tc);

		// pi/2
		tc = new TestCase (this);
		tc -> addInput ("X", mpz_class(0));
		tc -> addInput ("Y", mpz_class(1)<< (w-2));
		emulate(tc);
		tcl->add(tc);


		// 3pi/4
		tc = new TestCase (this);
		tc -> addInput ("X", (mpz_class(1)<< (w)) -  (mpz_class(1)<< (w-2)) );
		tc -> addInput ("Y", mpz_class(1)<< (w-2));
		emulate(tc);
		tcl->add(tc);


	}
	

}

