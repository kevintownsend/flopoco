#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCos.hpp"

using namespace std;

namespace flopoco{

	CordicSinCos::CordicSinCos(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
		int wcs = 1+w, wx = 1+w;
		//TODO: verify the validity of the necessary guard bits
		int guardx, guardcs;
		
		//initialize testing random number
		gmp_randinit_mt (state);
		
		guardcs = ceil(log2(wcs));
		guardx  = guardcs;
		
		srcFileName="CordicSinCos";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCos_" << 1+w <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCos_" << 1+w << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		addInput  ( "X"  , wcs, true );

		// declaring output
		addOutput  ( "C"  , wcs, 2 );
		addOutput  ( "S"  , wcs, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wcs + guardcs) + target->adderDelay(1+w) + target->lutDelay());
		
		//reduce the argument X to [0, 1/2)
		vhdl << tab << declare("absX", 1+w) << "<= (X xor (" << w << " downto 0 => X(" << wx-1 << ")))"
												   << " + " 
												   << "(" << zg(w, 0) << " & X(" << wx-1 << "));"<< endl;
		
		vhdl << tab << declare("reducedX", wx  + guardcs) << "<= absX(" << wx-1 << ") & \'0\' & absX(" << wx-3 << " downto 0) & " << zg(guardcs, 0) << ";" << endl;
		vhdl << tab << declare("quadrantX", 2) << " <= X(" << wx-1 << " downto " << wx-2 << ");" << endl;
		
		syncCycleFromSignal("reducedX");
		
		//create the C0, S0, X0 and D0 signals for the first stage
		mpf_t xinit;
		
		//compute the scale factor
		mpfr_t kfactor, temp;
		
		mpfr_init2(kfactor, 10*(wcs+guardcs));
		mpfr_set_d(kfactor, 1, GMP_RNDN);
		for(int i=0; i<(w+guardcs); i++){
			mpfr_init2(temp, 10*(wcs+guardcs));
			mpfr_set_d(temp, 1, GMP_RNDN);
			mpfr_mul_2si(temp, temp, (-2)*i, GMP_RNDN);
			mpfr_add_d(temp, temp, 1, GMP_RNDN);
			
			mpfr_mul(kfactor, kfactor, temp, GMP_RNDN);
		}
		mpfr_sqrt(kfactor, kfactor, GMP_RNDN);
		mpfr_d_div(kfactor, 1, kfactor, GMP_RNDN);
		
		mpf_init2(xinit, 10*(wcs+guardcs));
		mpfr_get_f(xinit, kfactor, GMP_RNDN);		
		
		manageCriticalPath(target->localWireDelay(wcs + guardcs) + target->lutDelay());
		
		vhdl << tab << declare("C0", wcs + guardcs) << "<= " << zg(1, 0) << " & \"" << generateFixPointNumber(xinit, 0, wcs-1+guardcs) << "\";" << endl;
		vhdl << tab << declare("S0", wcs + guardcs) << "<= " << zg(1, 0) << " & \"" << generateFixPointNumber(0.0, 0, wcs-1+guardcs) << "\";" << endl;
		vhdl << tab << declare("X0", wx  + guardcs) << "<= reducedX;" << endl;
		vhdl << tab << declare("D0") << "<= absX(" << w << ");" << endl;
		
		mpf_clear (xinit);
				
		//create the stages of micro-rotations
		int stage;
		
		wcs += guardcs;
		wx += guardx;
		
		//build the cordic stages
		for(stage=0; stage<wcs-1; stage++){
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
			manageCriticalPath(target->localWireDelay(wcs + guardcs) + target->adderDelay(wcs) + target->lutDelay());
			
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
			
			mpfr_set_d(zpow2, 1, GMP_RNDD);
			mpfr_mul_2si(zpow2, zpow2, (-1)*stage, GMP_RNDD);
			mpfr_atan(zatan, zpow2, GMP_RNDD);
			mpfr_const_pi( constPi, GMP_RNDD);
			mpfr_div(zatan, zatan, constPi, GMP_RNDD);
					
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
			
			manageCriticalPath(target->localWireDelay(wx + guardx) + target->adderDelay(wx) + target->lutDelay());
			
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
		
		manageCriticalPath(target->localWireDelay(1+w+1) + target->adderDelay(1+w+1));
		syncCycleFromSignal(join("D", stage));
		
		vhdl << tab << declare("preRoundedIntCout", 1+w+1) << "<= " << join("C", stage) << "(" << wcs-1 << " downto " << guardcs-1 << ");" << endl;
		vhdl << tab << declare("preRoundedIntSout", 1+w+1) << "<= " << join("S", stage) << "(" << wcs-1 << " downto " << guardcs-1 << ");" << endl;
		
		vhdl << tab << declare("roundedIntCout", 1+w+1) << "<= preRoundedIntCout "
															<< "+"
															<< " (" << zg(1+w, 0) << " & \'1\');" << endl;
		vhdl << tab << declare("roundedIntSout", 1+w+1) << "<= preRoundedIntSout "
															<< "+"
															<< " (" << zg(1+w, 0) << " & \'1\');" << endl;
															
		//assign output
		manageCriticalPath(target->localWireDelay(1+w) + target->adderDelay(1+w));
		
		vhdl << tab << declare("reducedC", 1+w) << "<= roundedIntCout(" << 1+w << " downto 1);" << endl;
		vhdl << tab << declare("reducedS", 1+w) << "<= roundedIntSout(" << 1+w << " downto 1);" << endl;
		
		vhdl << tab << declare("negReducedC", 1+w) << "<= (reducedC xor (" << w << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(w, 0) << " & \'1\');"<< endl;
		vhdl << tab << declare("negReducedS", 1+w) << "<= (reducedS xor (" << w << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(w, 0) << " & \'1\');"<< endl;
												   
		manageCriticalPath(target->localWireDelay(1+w) + target->lutDelay());
		
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "C <= reducedC when \"00\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"10\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedC when others;" << endl;	//edit to signal error
		
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "S <= reducedS when \"00\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedC when \"10\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedS when others;" << endl;	//edit to signal error
	};


	void CordicSinCos::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpfr_t z, constPi, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		int g = ceil(log2(1 + w));
		
		/* Compute correct value */
		mpfr_init2(z, 1+w+g);
		
		mpfr_init2(constPi, 1+w+g);
		
		mpfr_init2(rsin, 1+w+g); 
		mpfr_init2(rcos, 1+w+g); 
		mpz_init2 (rsin_z, 1+w+g);
		mpz_init2 (rcos_z, 1+w+g);
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDD); // this rounding is exact
		mpfr_div_2si (z, z, w, GMP_RNDD); // this rounding is acually exact
		
		mpfr_const_pi( constPi, GMP_RNDD);
		mpfr_mul(z, z, constPi, GMP_RNDD);
		
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDD); // exact rnd here
		mpfr_get_z (rsin_z, rsin, GMP_RNDN); // there can be a real rounding here
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDD); // exact rnd here
		mpfr_get_z (rcos_z, rcos, GMP_RNDN); // there can be a real rounding here

		// Set outputs 
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		tc->addExpectedOutput ("C", cos_zc);
		tc->addExpectedOutput ("S", sin_zc);

		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
		mpfr_free_cache();
	}


	void CordicSinCos::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpf_t zinit;
		mpfr_t z;
		mpz_t z_z;
		
		//mpf_set_default_prec (1+wI+wF+guardcs);
		
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


	std::string CordicSinCos::generateFixPointNumber(float x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		float xcopy = x;
		mpfr_t mx;
		mpz_class h;
		
		mpfr_init2 (mx, 1+wI+wF);
		
		if(xcopy<0){
			xcopy = xcopy * (-1);
		}
		
		mpfr_set_d(mx, xcopy, GMP_RNDD);
		mpfr_mul_2si(mx, mx, wF, GMP_RNDD);
		
		mpfr_get_z(h.get_mpz_t(), mx,  GMP_RNDD); 
        
        result = unsignedBinary(h, size);
        
        return result;
	}
	
	std::string CordicSinCos::generateFixPointNumber(mpf_t x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		mpfr_t mx;
		mpz_class h;
		
		mpfr_init2 (mx, 10*(1+wI+wF));
		
		if(x<0){
			mpf_neg (x, x);
		}
		
		mpfr_set_f(mx, x, GMP_RNDD);
		mpfr_mul_2si(mx, mx, wF, GMP_RNDD);
		
		mpfr_get_z(h.get_mpz_t(), mx,  GMP_RNDD);         
        result = unsignedBinary(h, size);
        
        return result;
	}
	
	mpz_class CordicSinCos::fp2fix(mpfr_t x, int wI, int wF){
		mpz_class h;
		
		mpfr_mul_2si(x, x, wF, GMP_RNDD);
        mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDD);  
		
		return h;
	}

}











