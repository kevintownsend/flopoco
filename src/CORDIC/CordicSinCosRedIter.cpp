#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCosRedIter.hpp"

using namespace std;

namespace flopoco{

	CordicSinCosRedIter::CordicSinCosRedIter(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
		int wcs = 1+w, wx = 1+w;
		//TODO: verify the validity of the necessary guard bits
		
		//initialize testing random number
		guard = ceil(log2(wcs));
		
		mpfr_init2(scale, 10*w);
		mpfr_set_d(scale, -1.0, GMP_RNDN);
		mpfr_mul_2si(scale, scale, -w, GMP_RNDN);
		mpfr_add_d(scale, scale, 1.0, GMP_RNDN);
		
		srcFileName="CordicSinCosRedIter";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCosRedIter_" << 1+w <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCosRedIter_" << 1+w << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		addInput  ( "X"  , wcs, true );

		// declaring output
		addOutput  ( "C"  , wcs, 2 );
		addOutput  ( "S"  , wcs, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wcs + guard) + target->adderDelay(1+w) + target->lutDelay());
		
		//reduce the argument X to [0, 1/2)
		vhdl << tab << declare("absX", 1+w) << "<= (X xor (" << w << " downto 0 => X(" << wx-1 << ")))"
												   << " + " 
												   << "(" << zg(w, 0) << " & X(" << wx-1 << "));"<< endl;
		
		vhdl << tab << declare("reducedX", wx  + guard) << "<= absX(" << wx-1 << ") & \'0\' & absX(" << wx-3 << " downto 0) & " << zg(guard, 0) << ";" << endl;
		vhdl << tab << declare("quadrantX", 2) << " <= X(" << wx-1 << " downto " << wx-2 << ");" << endl;
		
		syncCycleFromSignal("reducedX");
		
		//create the C0, S0, X0 and D0 signals for the first stage
		//compute the scale factor
		mpfr_t kfactor, temp;
			
		mpfr_init2(kfactor, 10*wcs);
		mpfr_set_d(kfactor, 1, GMP_RNDN);
		for(int i=0; i<=ceil((wcs+1+guard)/2.0); i++){
			mpfr_init2(temp, 10*wcs);
			mpfr_set_d(temp, 1, GMP_RNDN);
			mpfr_mul_2si(temp, temp, -2*i, GMP_RNDN);
			mpfr_add_d(temp, temp, 1.0, GMP_RNDN);
			
			mpfr_mul(kfactor, kfactor, temp, GMP_RNDN);
		}
		mpfr_sqrt(kfactor, kfactor, GMP_RNDN);
		mpfr_d_div(kfactor, 1.0, kfactor, GMP_RNDN);

		mpfr_mul(kfactor, kfactor, scale, GMP_RNDN);		
		
		manageCriticalPath(target->localWireDelay(wcs + guard) + target->lutDelay());
			
		vhdl << tab << declare("C0", wcs + guard) << "<= " << zg(1, 0) << " & \"" << generateFixPointNumber(kfactor, 0, wcs-1+guard) << "\";" << endl;
		vhdl << tab << declare("S0", wcs + guard) << "<= " << zg(1, 0) << " & " << zg(wcs-1+guard) << ";" << endl;
		vhdl << tab << declare("X0", wx  + guard) << "<= reducedX;" << endl;
		vhdl << tab << declare("D0") << "<= absX(" << w << ");" << endl;
				
		//create the stages of micro-rotations
		int stage;
		
		wcs += guard;
		wx += guard;
		
		//build the cordic stages
		for(stage=0; stage<=ceil((wcs+1)/2.0); stage++){
			double initCyclePathLength = getCriticalPath();
			
			manageCriticalPath(target->localWireDelay(wcs) + target->lutDelay());
			syncCycleFromSignal(join("D", stage));
			
			//shift Xin and Yin with 2^n positions to the right
			if(stage==0){
				vhdl << tab << declare(join("CShift", stage), wcs) << " <= C" << stage << ";" <<endl;
			}else{
				if(stage==1)
					vhdl << tab << declare(join("CsignExtend", stage), stage) << " <= C" << stage << "(" << wcs-1 << ");" <<endl;
				else
					vhdl << tab << declare(join("CsignExtend", stage), stage) << " <= (others => C" << stage << "(" << wcs-1 << "));" <<endl;
				vhdl << tab << declare(join("Cshifted", stage), wcs-stage) << " <= C" << stage << range(wcs-1, stage) << ";" <<endl;
				vhdl << tab << declare(join("CShift", stage), wcs) << " <= CsignExtend" << stage << " & Cshifted" << stage << ";" <<endl;
			}
			if(stage==0){
				vhdl << tab << declare(join("SShift", stage), wcs) << " <= S" << stage << ";" <<endl;
			}else{
				if(stage==1)
					vhdl << tab << declare(join("SsignExtend", stage), stage) << " <= S" << stage << "(" << wcs-1 << ");" <<endl;
				else
					vhdl << tab << declare(join("SsignExtend", stage), stage) << " <= (others => S" << stage << "(" << wcs-1 << "));" <<endl;
				vhdl << tab << declare(join("Sshifted", stage), wcs-stage) << " <= S" << stage << range(wcs-1, stage) << ";" <<endl;
				vhdl << tab << declare(join("SShift", stage), wcs) << " <= SsignExtend" << stage << " & Sshifted" << stage << ";" <<endl;
			}
			
			setCycleFromSignal(join("CShift", stage));
			syncCycleFromSignal(join("SShift", stage));
			manageCriticalPath(target->localWireDelay(wcs + guard) + target->adderDelay(wcs) + target->lutDelay());
			
			vhdl << tab << declare(join("intC", stage+1), wcs) << " <= " << join("C", stage) << " + " << join("SShift", stage) << " when " << join("D", stage) << "=\'1\' else "
				 << join("C", stage) << " - " << join("SShift", stage) << " ;" << endl;
			vhdl << tab << declare(join("intS", stage+1), wcs) << " <= " << join("S", stage) << " - " << join("CShift", stage) << " when " << join("D", stage) << "=\'1\' else "
				 << join("S", stage) << " + " << join("CShift", stage) << " ;" << endl;

			double csPath = getCriticalPath();
			
			//create the constant signal for the arctan
			mpfr_t zatan, zpow2, constPi;
			
			mpfr_init2(zatan, 10*wcs);
			mpfr_init2(zpow2, 10*wcs);
			mpfr_init2(constPi, 10*wcs);
			
			mpfr_set_d(zpow2, 1.0, GMP_RNDN);
			mpfr_mul_2si(zpow2, zpow2, (-1)*stage, GMP_RNDN);
			mpfr_atan(zatan, zpow2, GMP_RNDN);
			mpfr_const_pi( constPi, GMP_RNDN);
			mpfr_div(zatan, zatan, constPi, GMP_RNDN);
					
			//create the arctangent factor to be added to Zin
			bool negAtan = false;
			std::string strConverted;
			mpz_class fixConverted;
			
			fixConverted = fp2fix(zatan, 0, wcs-1);
			
			if(fixConverted<0){
				fixConverted = fixConverted * (-1);
				negAtan = true;
			}
			strConverted = unsignedBinary(fixConverted, wx-1);
			
			mpfr_free_cache();
			
			setCycleFromSignal(join("D", stage));
			setCriticalPath(initCyclePathLength);
			
			manageCriticalPath(target->localWireDelay(wx + guard) + target->adderDelay(wx) + target->lutDelay());
			
			vhdl << tab << declare(join("atan2PowStage", stage), wx) << " <= \'0\' & \"" << strConverted << "\";" <<endl;
			
			if(negAtan){
				vhdl << tab << declare(join("intX", stage+1), wx) << " <= " << join("X", stage) << " - " << join("atan2PowStage", stage) << " when " << join("D", stage) << "=\'1\' else "
				 << join("X", stage) << " + " << join("atan2PowStage", stage) << " ;" << endl;
			}
			else{
				vhdl << tab << declare(join("intX", stage+1), wx) << " <= " << join("X", stage) << " + " << join("atan2PowStage", stage) << " when " << join("D", stage) << "=\'1\' else "
				 << join("X", stage) << " - " << join("atan2PowStage", stage) << " ;" << endl;
			}
			
 			//create Dout as the result of the comparison of intZout with 0
			vhdl << tab << declare(join("intD", (stage+1))) << " <= intX" << stage+1 << "(" << wx-1 <<");" <<endl;
			
			double xPath = getCriticalPath();
			
			if (getCycleFromSignal(join("intC", (stage+1))) == getCycleFromSignal(join("intD", (stage+1)))){
				setCriticalPath(max(xPath, csPath));
				syncCycleFromSignal(join("intD", stage+1));
			}else
				if (syncCycleFromSignal(join("intC", (stage+1)))){
					setCriticalPath(csPath);
					syncCycleFromSignal(join("intC", stage));
				}else{
					setCriticalPath(xPath);
					syncCycleFromSignal(join("intD", stage+1));
				}
					
			//create the outputs
			vhdl << tab << declare(join("C", stage+1), wcs) << " <= intC" << stage+1 << ";" <<endl;
			vhdl << tab << declare(join("S", stage+1), wcs) << " <= intS" << stage+1 << ";" <<endl;
			if(wx > 2)
				vhdl << tab << declare(join("X", stage+1), wx-1) << " <= intX" << stage+1 << "(" << wx-2 << " downto 0);" <<endl;
			else
				vhdl << tab << declare(join("X", stage+1)) << " <= intX" << stage+1 << "(0);" <<endl;
			vhdl << tab << declare(join("D", stage+1)) << " <= intD" << stage+1 << ";" <<endl;
			
			//decrement the size of Z
			wx--;
		}
		
		//multiply X by Pi
		FixRealKCM* piMultiplier = new FixRealKCM(target, -(wx+stage-1), -stage-1, 1, -(wx+stage-1), "pi");
		oplist.push_back(piMultiplier);
		
		inPortMap(piMultiplier, "X", join("X", stage));
		outPortMap(piMultiplier, "R", "xMulPi");
		vhdl << instance(piMultiplier, "piMultiplier") << endl;
		
		wx += 2;	// TODO: check 
					// NOTE: the first factor in the sum corresponds to the compensation due to Pi multiplication
		
		vhdl << tab << declare("truncC", wx) << " <= " << join("C", stage) << range(wcs-1, wcs-wx) << ";" << endl;
		vhdl << tab << declare("truncS", wx) << " <= " << join("S", stage) << range(wcs-1, wcs-wx) << ";" << endl;
		
		//multiply with the angle X to obtain the actual values for sine and cosine
		IntMultiplier* zmultiplier = new IntMultiplier(target, wx, wx, true, inDelayMap("X",getCriticalPath()), 1.0);
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
		setCriticalPath(zmultiplier->getOutputDelay("xMultiplierC"));
		 
		
		manageCriticalPath(target->localWireDelay(wcs) + target->adderDelay(wcs));
		
		vhdl << tab << declare("shortCX", wx) << "<= CX(" << 2*wx-1 << " downto " << wx << ");" << endl;
		vhdl << tab << declare("shortSX", wx) << "<= SX(" << 2*wx-1 << " downto " << wx << ");" << endl;
		
		vhdl << tab << declare("paddedShortCX", wcs) << " <= (" << wcs-1 << " downto " << wx << " => CX(" << 2*wx-1 << ")) & shortCX;"  << endl;
		vhdl << tab << declare("paddedShortSX", wcs) << " <= (" << wcs-1 << " downto " << wx << " => SX(" << 2*wx-1 << ")) & shortSX;"  << endl;
		
		vhdl << tab << declare("CsubShortSX", wcs) << " <= " << join("C", stage) << " - paddedShortSX;" << endl;
		vhdl << tab << declare("SaddShortCX", wcs) << " <= " << join("S", stage) << " + paddedShortCX;" << endl;
		
		//assign output
		manageCriticalPath(target->localWireDelay(1+w) + target->adderDelay(1+w));
		
		vhdl << tab << declare("reducedC", wcs) << "<= CsubShortSX;" << endl;
		vhdl << tab << declare("reducedS", wcs) << "<= SaddShortCX;" << endl;
		
		vhdl << tab << declare("negReducedC", wcs) << "<= (reducedC xor (" << wcs-1 << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(wcs-1, 0) << " & \'1\');"<< endl;
		vhdl << tab << declare("negReducedS", wcs) << "<= (reducedS xor (" << wcs-1 << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(wcs-1, 0) << " & \'1\');"<< endl;
												   
		manageCriticalPath(target->localWireDelay(1+w) + target->lutDelay() + target->adderDelay(1+w+1));
		
		declare("effectiveC", wcs);
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "effectiveC <= reducedC when \"00\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"10\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedC when others;" << endl;	//edit to signal error
		
		declare("effectiveS", wcs);
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "effectiveS <= reducedS when \"00\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedC when \"10\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedS when others;" << endl;	//edit to signal error
		
		vhdl << tab << declare("roundedEffectiveC", 1+w+1) << " <= effectiveC(" << wcs-1 << " downto " << guard-1 << ") "
														   << "+"
														   << " (" << zg(1+w, 0) << " & \'1\');" << endl;
		vhdl << tab << declare("roundedEffectiveS", 1+w+1) << " <= effectiveS(" << wcs-1 << " downto " << guard-1 << ") "
														   << "+"
														   << " (" << zg(1+w, 0) << " & \'1\');" << endl;
														   
		vhdl << tab << "C <= roundedEffectiveC(" << 1+w << " downto 1);" << endl;
		vhdl << tab << "S <= roundedEffectiveS(" << 1+w << " downto 1);" << endl;
	};


	void CordicSinCosRedIter::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpfr_t z, constPi, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		
		/* Compute correct value */
		mpfr_init2(z, 10*w);
		
		mpfr_init2(constPi, 10*w);
		
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, w, GMP_RNDN); // this rounding is acually exact
		
		mpfr_const_pi( constPi, GMP_RNDN);
		mpfr_mul(z, z, constPi, GMP_RNDN);
		
		mpfr_init2(rsin, 10*w); 
		mpfr_init2(rcos, 10*w); 


		mpz_init2 (rsin_z, 2*w);
		mpz_init2 (rcos_z, 2*w);
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		mpfr_mul(rsin, rsin, scale, GMP_RNDN);
		mpfr_mul(rcos, rcos, scale, GMP_RNDN);

		mpfr_add_d(rsin, rsin, 6.0, GMP_RNDN);
		mpfr_add_d(rcos, rcos, 6.0, GMP_RNDN);
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDN); // exact rnd here

		// Rounding down
		mpfr_get_z (rsin_z, rsin, GMP_RNDD); // there can be a real rounding here
		mpfr_get_z (rcos_z, rcos, GMP_RNDD); // there can be a real rounding here
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		sin_zc -= mpz_class(6)<<w;
		cos_zc -= mpz_class(6)<<w;

		tc->addExpectedOutput ("S", sin_zc);
		tc->addExpectedOutput ("C", cos_zc);
		
		// Rounding up
		mpfr_get_z (rsin_z, rsin, GMP_RNDU); // there can be a real rounding here
		mpfr_get_z (rcos_z, rcos, GMP_RNDU); // there can be a real rounding here
		mpz_class sin_zcu (rsin_z), cos_zcu (rcos_z);
		sin_zcu -= mpz_class(6)<<w;
		cos_zcu -= mpz_class(6)<<w;

		tc->addExpectedOutput ("S", sin_zcu);
		tc->addExpectedOutput ("C", cos_zcu);
		


		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
		mpfr_free_cache();
	}



	void CordicSinCosRedIter::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpf_t zinit;
		mpfr_t z;
		mpz_t z_z;
		
		//mpf_set_default_prec (1+wI+wF+guard);
		
		mpfr_init2(z, 1+w+ceil(log2(1 + w)));
		mpz_init2 (z_z, 1+w+ceil(log2(1 + w)));
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "1.5707963267949e0", 10);
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "0.5235987755983e0", 10);
		mpf_set_str (zinit, "0.16666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "0.78539816339745e0", 10);
		mpf_set_str (zinit, "0.25e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "1.0471975511966e0", 10);
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD);
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		mpfr_clears (z, NULL);
	}


	std::string CordicSinCosRedIter::generateFixPointNumber(mpfr_t x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		mpz_class h;
		
		if(x<0){
			mpfr_neg (x, x,  GMP_RNDN);
		}
		
		mpfr_mul_2si(x, x, wF, GMP_RNDN);
		
		mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDN);         
		result = unsignedBinary(h, size);
		
		return result;
	}
	
	mpz_class CordicSinCosRedIter::fp2fix(mpfr_t x, int wI, int wF){
		mpz_class h;
		
		mpfr_mul_2si(x, x, wF, GMP_RNDN);
		mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDN);  
		
		return h;
	}

}












