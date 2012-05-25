#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCos.hpp"

using namespace std;

namespace flopoco{



	CordicSinCos::CordicSinCos(Target* target, int wIn_, int wOut_, int reducedIterations_, map<string, double> inputDelays) 
		: Operator(target), wIn(wIn_), wOut(wOut_), reducedIterations(reducedIterations_)
	{

		int stage;
		srcFileName="CordicSinCos";
		setCopyrightString ( "Matei Istoan, Florent de Dinechin (2012-...)" );

		ostringstream name;
		name << "CordicSinCos_" << (reducedIterations==1?"reducedIterations":"") << wIn << "_" << wOut;
		if(target->isPipelined())
			name  <<"_f" << target->frequencyMHz();
		else 
			name << "_comb";
		name << "_uid" << getNewUId();
		setName( name.str() );

		if(wIn<12){
			REPORT(INFO, "wIn is small, are you sure you don't want to tabulate this operator in a ROM?");
		}

		//TODO check all the following
		if (reducedIterations == 1)
			maxIterations=((wOut+1)>>1);
		else
			maxIterations = wOut+1;
		
#define ROUNDED_ROTATION 1 // 0:trunc 

#if ROUNDED_ROTATION
		REPORT(DEBUG, "Using rounded rotation trick");
#endif

		//error analysis
		double eps;  //error in ulp
		eps=0.5; //initial rounding of kfactor
		double shift=0.5;
		for(stage=1; stage<=maxIterations; stage++){
#if ROUNDED_ROTATION
			eps = eps + eps*shift + 0.5; // 0.5 assume rounding in the rotation.
#else
			eps = eps + eps*shift + 1.0; // 1.0 assume truncation in the rotation.
#endif
			shift *=0.5;
		}

		eps+=1; // the final neg-by-not
		REPORT(DEBUG, "Error analysis computes eps=" << eps << " ulps (before final rounding)");

		// guard bits depend only on the number of iterations
		g = 1+(int) ceil(log2(eps)); // +1 for the final rounding 

		
		// *********    internal precision and fixed-point alignment **************

 		// The input is as follows:
		// s has weight 2^0
		// q has weight 2^-1
		// o has weight 2^-2
		// Purpose: have the LSB of z, cosine, sine at weight -w
		// This means that the Cos and Sin datapath will have w+1 bits
		//   (sign at weight zero)
		// while the Z datapath starts on w-1 bits (sign bit at weight  -2,
		//  and decreasing, so the invariant is: sign bit at weight -stage-1)

		w = wOut-1 + g; // -1 because of sign

		REPORT(DEBUG, "wIn=" << wIn << " wOut=" << wOut 
		       << "   MaxIterations: " << maxIterations 
		       << "  Guard bits g=" << g << "  Neg. weight of LSBs w=" << w );


		// everybody needs many digits of Pi
		mpfr_init2(constPi, 10*w);
		mpfr_const_pi( constPi, GMP_RNDN);

		//compute the scale factor		
		mpfr_init2(scale, wOut+2);
		mpfr_set_d(scale, -1.0, GMP_RNDN);           // exact
		mpfr_mul_2si(scale, scale, -wOut+1, GMP_RNDN); // exact
		mpfr_add_d(scale, scale, 1.0, GMP_RNDN);     // exact
		REPORT(DEBUG, "scale=" << printMPFR(scale, 15));
		

		// declaring inputs
		addInput  ( "X"  , wIn, true );

		// declaring output
		addOutput  ( "C"  , wOut, 2 );
		addOutput  ( "S"  , wOut, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		manageCriticalPath( target->lutDelay());
		
		//reduce the argument X to [0, 1/2)
		vhdl << tab << declare("sgn") << " <= X(" << wIn-1 << ");  -- sign" << endl;
		vhdl << tab << declare("q") << " <= X(" << wIn-2 << ");  -- quadrant" << endl;
		vhdl << tab << declare("o") << " <= X(" << wIn-3 << ");  -- octant" << endl;
		vhdl << tab << declare("sqo", 3) << " <= sgn & q & o;  -- sign, quadrant, octant" << endl;
		//vhdl << tab << declare("sq", 2) << " <= sgn & q;  -- sign, quadrant" << endl;

		vhdl << tab << declare("qrot0", 3) << " <= sqo +  \"001\"; -- rotate by an octant" << endl; 
		vhdl << tab << declare("qrot", 2) << " <= qrot0(2 downto 1); -- new quadrant: 00 is around the origin" << endl; 
		vhdl << tab << declare("sred") << " <= o xor qrot(0); -- reduced sign" << endl; 
		// y is built out of the remaining wIn-2 bits



		int zMSB=-2;   // -2 'cause we live in a quadrant, initial angle is in -0.25pi, 0.25pi
		int zLSB=-w-1; // better have a bit more accuracy here, it seems 
		// extract the bits below the octant bit, the new sign will be at the weight of the octant bit
		int sizeZ =  zMSB - zLSB +1;
 		int sizeY=wIn-2;
		vhdl << tab << declare("Yp", sizeZ ) << "<= " ;
		if(sizeZ >= sizeY) {
			vhdl << "X" << range(sizeY-1,0); // sizeY-1 = wIn-4,  i.e. MSB = -3
			if(sizeZ > sizeY)
				vhdl << " & " << zg(sizeZ-sizeY) << ";" << endl;
		}
		else // sizeZ < sizeY
			vhdl << "X" << range(sizeY-1, sizeY-sizeZ) <<";" << endl;

		vhdl << tab << "--  This Yp is in -pi/4, pi/4. Now start CORDIC with angle atan(1/2)" << endl;

		//create the C1, S1, X1 and D1 signals for the first stage
		mpfr_t temp, zatan;

 		mpfr_init2(kfactor, 10*w);
		mpfr_init2(temp, 10*w);
		mpfr_set_d(kfactor, 1.0, GMP_RNDN);
		for(int i=1; i<=maxIterations; i++){
			mpfr_set_d(temp, 1.0, GMP_RNDN);
			mpfr_div_2si(temp, temp, 2*i, GMP_RNDN);
			mpfr_add_d(temp, temp, 1.0, GMP_RNDN);
			
			mpfr_mul(kfactor, kfactor, temp, GMP_RNDN);
		}
		mpfr_sqrt(kfactor, kfactor, GMP_RNDN);
		mpfr_d_div(kfactor, 1.0, kfactor, GMP_RNDN);

		mpfr_mul(kfactor, kfactor, scale, GMP_RNDN);
		
		REPORT(DEBUG, "kfactor=" << printMPFR(kfactor, 15));
		mpfr_clear(temp);
		

		vhdl << tab << declare("Cos1", w+1) << " <= " <<  unsignedFixPointNumber(kfactor, 0, -w) << ";" 
		     << "-- scale factor, about " << printMPFR(kfactor, 15) << endl;
		vhdl << tab << declare("Sin1", w+1) << " <= " << zg(w+1) << ";" << endl;
		vhdl << tab << declare("Z1", sizeZ) << "<= Yp;" << endl;
		vhdl << tab << declare("D1") << "<= Yp" << of(sizeZ-1) << ";" << endl;
		
				
		//create the stages of micro-rotations

		//build the cordic stages
		for(stage=1; stage<=maxIterations; stage++){
			double initCyclePathLength = getCriticalPath();
			
			manageCriticalPath(target->localWireDelay(w) + target->lutDelay());
			syncCycleFromSignal(join("D", stage));
			
			//shift Xin and Yin with 2^n positions to the right
			// Cosine is always positive, but sine may be negative and thus need sign extend
			vhdl << tab << declare(join("CosShift", stage), w+1) << " <= " << zg(stage) << " & Cos" << stage << range(w, stage) << ";" <<endl;
			
			vhdl << tab << declare(join("sgnSin", stage))  << " <= " <<  join("Sin", stage)  <<  of(w) << ";" << endl; 
			vhdl << tab << declare(join("SinShift", stage), w+1) 
			     << " <= " << rangeAssign(w, w+1-stage, join("sgnSin", stage))   
			     << " & Sin" << stage << range(w, stage) << ";" << endl;
			
			setCycleFromSignal(join("CosShift", stage));
			syncCycleFromSignal(join("SinShift", stage));
			manageCriticalPath(target->localWireDelay(w) + target->adderDelay(w) + target->lutDelay());
			
#if ROUNDED_ROTATION  // rounding of the shifted operand, should save 1 bit in each addition
			
#if 1
			vhdl << tab << declare(join("CosShiftRoundBit", stage)) << " <= " << join("Cos", stage)  << of(stage-1) << ";" << endl;
			vhdl << tab << declare(join("SinShiftRoundBit", stage)) << " <= " << join("Sin", stage) << of(stage-1) << ";" <<endl;
#else
			vhdl << tab << declare(join("CosShiftRoundBit", stage)) << " <= '0';" << endl;
			vhdl << tab << declare(join("SinShiftRoundBit", stage)) << " <= '0';" <<endl;
#endif
			vhdl << tab << declare(join("CosShiftNeg", stage), w+1) << " <= " << rangeAssign(w, 0, join("D", stage)) << " xor " << join("CosShift", stage)   << " ;" << endl;
			vhdl << tab << declare(join("SinShiftNeg", stage), w+1) << " <= (not " << rangeAssign(w, 0, join("D", stage)) << ") xor " << join("SinShift", stage)   << " ;" << endl;

			vhdl << tab << declare(join("Cos", stage+1), w+1) << " <= " 
			     << join("Cos", stage) << " + " << join("SinShiftNeg", stage) << " +  not (" << join("D", stage) << " xor " << join("SinShiftRoundBit", stage) << ") ;" << endl;

			vhdl << tab << declare(join("Sin", stage+1), w+1) << " <= " 
			     << join("Sin", stage) << " + " << join("CosShiftNeg", stage) << " + (" << join("D", stage) << " xor " << join("CosShiftRoundBit", stage) << ") ;" << endl;

#else
// truncation of the shifted operand
			vhdl << tab << declare(join("Cos", stage+1), w+1) << " <= " 
			     << join("Cos", stage) << " - " << join("SinShift", stage) << " when " << join("D", stage) << "=\'0\' else "
			     << join("Cos", stage) << " + " << join("SinShift", stage) << " ;" << endl;

			vhdl << tab << declare(join("Sin", stage+1), w+1) << " <= " 
			     << join("Sin", stage) << " + " << join("CosShift", stage) << " when " << join("D", stage) << "=\'0\' else "
			     << join("Sin", stage) << " - " << join("CosShift", stage) << " ;" << endl;

#endif
			double csPath = getCriticalPath();
			
			
			//create the constant signal for the arctan
			mpfr_init2(zatan, 10*w);
			
			mpfr_set_d(zatan, 1.0, GMP_RNDN);
			mpfr_div_2si(zatan, zatan, stage, GMP_RNDN);
			mpfr_atan(zatan, zatan, GMP_RNDN);
			mpfr_div(zatan, zatan, constPi, GMP_RNDN);
			REPORT(DEBUG, "stage=" << stage << "  atancst=" << printMPFR(zatan, 15));		
			//create the arctangent factor to be added to Zin
			
			setCycleFromSignal(join("D", stage));
			setCriticalPath(initCyclePathLength);
			
			manageCriticalPath(target->localWireDelay(sizeZ) + target->adderDelay(sizeZ) + target->lutDelay());
			
			REPORT(DEBUG, "  sizeZ=" << sizeZ << "   zMSB="<<zMSB );

			if(stage<maxIterations) {
				// LSB is always -w 
				vhdl << tab << declare(join("atan2PowStage", stage), sizeZ) << " <= " << unsignedFixPointNumber(zatan, zMSB, zLSB) << ";" <<endl;
				
				vhdl << tab << declare(join("fullZ", stage+1), sizeZ) << " <= " 
				     << join("Z", stage) << " + " << join("atan2PowStage", stage) << " when " << join("D", stage) << "=\'1\' else "
				     << join("Z", stage) << " - " << join("atan2PowStage", stage) << " ;" << endl;
				vhdl << tab << declare(join("Z", stage+1), sizeZ-1) << " <= "<< join("fullZ", stage+1) << range(sizeZ-2, 0) << ";" << endl; 
				vhdl << tab << declare(join("D", (stage+1))) << " <= fullZ" << stage+1 << "(" << sizeZ-1 <<");" <<endl;
			}
			//decrement the size of Z
			sizeZ--;
			zMSB--;
		}
			
		if(reducedIterations == 0){ //regular struture; all that remains is to assign the outputs correctly
			
			syncCycleFromSignal(join("D", stage));
			
			//assign output
			manageCriticalPath(target->localWireDelay(w) + target->adderDelay(w));
			
			vhdl << tab << declare("redCos", w+1) << "<= " << join("Cos", stage) << ";" << endl;
			vhdl << tab << declare("redSin", w+1) << "<= " << join("Sin", stage) << ";" << endl;
			

		}
		else{	//reduced iterations structure; rotate by the remaining angle and then assign the angles
			
			manageCriticalPath(target->localWireDelay(sizeZ) + target->DSPMultiplierDelay());
		
			//multiply X by Pi
			FixRealKCM* piMultiplier = new FixRealKCM(target, -(sizeZ+stage-1), -stage-1, 1, -(sizeZ+stage-1), "pi");
			oplist.push_back(piMultiplier);
			
			inPortMap(piMultiplier, "X", join("X", stage));
			outPortMap(piMultiplier, "R", "xMulPi");
			vhdl << instance(piMultiplier, "piMultiplier") << endl;
		
			syncCycleFromSignal("xMulPi");
			setCriticalPath(piMultiplier->getOutputDelay("R"));


			sizeZ += 2;	// TODO: check 
						// NOTE: the first factor in the sum corresponds to the compensation due to Pi multiplication
			
			vhdl << tab << declare("truncC", sizeZ) << " <= " << join("C", stage) << range(w-1, w-sizeZ) << ";" << endl;
			vhdl << tab << declare("truncS", sizeZ) << " <= " << join("S", stage) << range(w-1, w-sizeZ) << ";" << endl;
			
			//multiply with the angle X to obtain the actual values for sine and cosine
			IntMultiplier* zmultiplier = new IntMultiplier(target, sizeZ, sizeZ, true, inDelayMap("X",getCriticalPath()), 1.0);
			oplist.push_back(zmultiplier);
			
			inPortMap(zmultiplier, "X", "truncC");
			inPortMap(zmultiplier, "Y", "xMulPi");
			outPortMap(zmultiplier, "R", "CX");
			vhdl << instance(zmultiplier, "xMultiplierC") << endl;
			
			inPortMap(zmultiplier, "X", "truncS");
			inPortMap(zmultiplier, "Y", "xMulPi");
			outPortMap(zmultiplier, "R", "SX");
			vhdl << instance(zmultiplier, "xMultiplierS") << endl;
			
			setCycleFromSignal("CX");
			syncCycleFromSignal("SX");
			nextCycle(); // TODO

			//setCriticalPath(zmultiplier->getOutputDelay("R"));
			 
			
			manageCriticalPath(target->localWireDelay(w) + target->adderDelay(w));
			
			vhdl << tab << declare("shortCX", sizeZ) << "<= CX(" << 2*sizeZ-1 << " downto " << sizeZ << ");" << endl;
			vhdl << tab << declare("shortSX", sizeZ) << "<= SX(" << 2*sizeZ-1 << " downto " << sizeZ << ");" << endl;
			
			vhdl << tab << declare("paddedShortCX", w) << " <= (" << w-1 << " downto " << sizeZ << " => CX(" << 2*sizeZ-1 << ")) & shortCX;"  << endl;
			vhdl << tab << declare("paddedShortSX", w) << " <= (" << w-1 << " downto " << sizeZ << " => SX(" << 2*sizeZ-1 << ")) & shortSX;"  << endl;
			
			vhdl << tab << declare("CsubShortSX", w) << " <= " << join("C", stage) << " - paddedShortSX;" << endl;
			vhdl << tab << declare("SaddShortCX", w) << " <= " << join("S", stage) << " + paddedShortCX;" << endl;
			
			//assign output
			manageCriticalPath(target->localWireDelay(w) + target->adderDelay(w));
			
			vhdl << tab << declare("reducedC", w) << "<= CsubShortSX;" << endl;
			vhdl << tab << declare("reducedS", w) << "<= SaddShortCX;" << endl;
			
		}


		
		vhdl << tab << "---- final reconstruction " << endl;
		
		vhdl << tab << declare("redCosNeg", w+1) << " <= (not redCos); -- negate by NOT, 1 ulp error"<< endl;
		vhdl << tab << declare("redSinNeg", w+1) << " <= (not redSin); -- negate by NOT, 1 ulp error"<< endl;
												   
		manageCriticalPath(target->localWireDelay(w) + target->lutDelay());
		
		vhdl << tab << "with qrot select" << endl
		     << tab << tab << declare("CosX0", w+1) << " <= " << endl;
		vhdl << tab << tab << tab << " redCos    when \"00\"," << endl;
		vhdl << tab << tab << tab << " redSinNeg when \"01\"," << endl;
		vhdl << tab << tab << tab << " redCosNeg when \"10\"," << endl;
		vhdl << tab << tab << tab << " redSin    when others;" << endl;

		vhdl << tab << "with qrot select" << endl
		      << tab << tab << declare("SinX0", w+1) << " <= " << endl;
		vhdl << tab << tab << tab << " redSin    when \"00\"," << endl;
		vhdl << tab << tab << tab << " redCos    when \"01\"," << endl;
		vhdl << tab << tab << tab << " redSinNeg when \"10\"," << endl;
		vhdl << tab << tab << tab << " redCosNeg when others;" << endl;
		
		
		manageCriticalPath(target->localWireDelay(w) +  target->adderDelay(1+w+1));

		vhdl << tab << declare("roundedCosX", wOut+1) << " <= CosX0" << range(w, w-wOut) << " + " << " (" << zg(wOut) << " & \'1\');" << endl;
		vhdl << tab << declare("roundedSinX", wOut+1) << " <= SinX0" << range(w, w-wOut) << " + " << " (" << zg(wOut) << " & \'1\');" << endl;
														   
		vhdl << tab << "C <= roundedCosX" << range(wOut, 1) << ";" << endl;
		vhdl << tab << "S <= roundedSinX" << range(wOut, 1) << ";" << endl;


		mpfr_clears (zatan, NULL);		
	};


	CordicSinCos::~CordicSinCos(){
		mpfr_clears (scale, kfactor, constPi, NULL);		
	 };


	void CordicSinCos::emulate(TestCase * tc) 
	{
		mpfr_t z, rsin, rcos;
		mpz_class sin_z, cos_z;
		mpfr_init2(z, 10*w);
		mpfr_init2(rsin, 10*w); 
		mpfr_init2(rcos, 10*w); 
		
														  

		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		
		/* Compute correct value */
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); //  exact
		mpfr_div_2si (z, z, wIn-1, GMP_RNDN); // exact
	
		// No need to manage sign bit etc: modulo 2pi is the same as modulo 2 in the initial format
		mpfr_mul(z, z, constPi, GMP_RNDN);

		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		mpfr_mul(rsin, rsin, scale, GMP_RNDN);
		mpfr_mul(rcos, rcos, scale, GMP_RNDN);

		mpfr_add_d(rsin, rsin, 6.0, GMP_RNDN); // exact rnd here
		mpfr_add_d(rcos, rcos, 6.0, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rsin, rsin, wOut-1, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rcos, rcos, wOut-1, GMP_RNDN); // exact rnd here

		// Rounding down
		mpfr_get_z (sin_z.get_mpz_t(), rsin, GMP_RNDD); // there can be a real rounding here
		mpfr_get_z (cos_z.get_mpz_t(), rcos, GMP_RNDD); // there can be a real rounding here
		sin_z -= mpz_class(6)<<(wOut-1);
		cos_z -= mpz_class(6)<<(wOut-1);

		tc->addExpectedOutput ("S", sin_z);
		tc->addExpectedOutput ("C", cos_z);

		// Rounding up
		mpfr_get_z (sin_z.get_mpz_t(), rsin, GMP_RNDU); // there can be a real rounding here
		mpfr_get_z (cos_z.get_mpz_t(), rcos, GMP_RNDU); // there can be a real rounding here
		sin_z -= mpz_class(6)<<(wOut-1);
		cos_z -= mpz_class(6)<<(wOut-1);

		tc->addExpectedOutput ("S", sin_z);
		tc->addExpectedOutput ("C", cos_z);
		
		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
	}






	void CordicSinCos::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpfr_t z;
		mpz_class zz;
		
		
		mpfr_init2(z, 10*w);
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X", mpz_class(0));
		emulate(tc);
		tcl->add(tc);
					
		tc = new TestCase (this);
		tc->addComment("Pi/4-eps");
		mpfr_set_d (z, 0.24, GMP_RNDD); 
		mpfr_mul_2si (z, z, wIn-1, GMP_RNDD); 
		mpfr_get_z (zz.get_mpz_t(), z, GMP_RNDD);  
		tc -> addInput ("X", zz);
		emulate(tc);
		tcl->add(tc);
				
		tc = new TestCase (this);
		tc->addComment("Pi/6");
		mpfr_set_d (z, 0.166666666666666666666666666666666, GMP_RNDD); 
		mpfr_mul_2si (z, z, wIn-1, GMP_RNDD); 
		mpfr_get_z (zz.get_mpz_t(), z, GMP_RNDD);  
		tc -> addInput ("X", zz);
		emulate(tc);
		tcl->add(tc);
				
		tc = new TestCase (this);
		tc->addComment("Pi/3");
		mpfr_set_d (z, 0.333333333333333333333333333333, GMP_RNDD); 
		mpfr_mul_2si (z, z, wIn-1, GMP_RNDD); 
		mpfr_get_z (zz.get_mpz_t(), z, GMP_RNDD);  
		tc -> addInput ("X", zz);
		emulate(tc);
		tcl->add(tc);
				
		
		mpfr_clears (z, NULL);
	}


	
	mpz_class CordicSinCos::fp2fix(mpfr_t x, int wI, int wF){
		mpz_class h;
		
		mpfr_mul_2si(x, x, wF, GMP_RNDN);
		mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDN);  
		
		return h;
	}

}




#if 0 
		// Rough simulation, for debug purpose
		double c,cc,s,ss,z, z0, dd,p,scale;
		int d;
		const double pi=3.14159265358979323846264338327950288419716939937508;
		c=mpfr_get_d(kfactor, GMP_RNDN);
		s=0.0;
		z=0.15625; // 1/6
		z0=z;
		p=0.5;
		scale=(double) (1<<19);
		for(stage=1; stage<=maxIterations; stage++){
			if(z>=0) d=0; else d=1;
			if(d==0) dd=1.0; else dd=-1.0;
			cc = c - dd*p*s;
			ss = s + dd*p*c;
			cout << stage << "\t atan=" << atan(p)/pi<< "  \t d=" << d << "\t z=" << z << "\t c=" << c << "\t s=" << s;
			cout  << "      \t z=" << (int)(z*scale) << "  \t c=" << (int)(c*scale*4) << "\t s=" << (int)(s*scale*4) << endl;
			z = z - dd*atan(p)/pi;
			c=cc;
			s=ss;
			p=p*0.5;
		}		
		cout  << "Should be  \t\t\t\t\t\t c=" << cos(pi*z0) << "  \t s=" << sin(pi*z0) << endl;
		//		manageCriticalPath(target->localWireDelay(wcs + g) + target->lutDelay());
#endif		
