#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCosClassic.hpp"

using namespace std;

namespace flopoco{

	CordicSinCosClassic::CordicSinCosClassic(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
		int wcs = 1+w, wx = 1+w;
		//TODO: verify the validity of the necessary guard bits
		
		//initialize testing random number
		gmp_randinit_mt (state);
		
		guard = ceil(log2(wcs));
		
		srcFileName="CordicSinCosClassic";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCosClassic_" << 1+w <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCosClassic_" << 1+w << "_uid" << getNewUId();
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
		mpf_t xinit;
		
		//compute the scale factor
		mpfr_t kfactor, temp;
		
		mpfr_init2(kfactor, 10*(wcs+guard));
		mpfr_set_d(kfactor, 1, GMP_RNDN);
		for(int i=0; i<(w+guard); i++){
			mpfr_init2(temp, 10*(wcs+guard));
			mpfr_set_d(temp, 1, GMP_RNDN);
			mpfr_mul_2si(temp, temp, (-2)*i, GMP_RNDN);
			mpfr_add_d(temp, temp, 1, GMP_RNDN);
			
			mpfr_mul(kfactor, kfactor, temp, GMP_RNDN);
		}
		mpfr_sqrt(kfactor, kfactor, GMP_RNDN);
		mpfr_d_div(kfactor, 1, kfactor, GMP_RNDN);
		
		mpf_init2(xinit, 10*(wcs+guard));
		mpfr_get_f(xinit, kfactor, GMP_RNDN);		
		
		manageCriticalPath(target->localWireDelay(wcs + guard) + target->lutDelay());
		
		vhdl << tab << declare("C0", wcs + guard) << "<= " << zg(1, 0) << " & \"" << generateFixPointNumber(xinit, 0, wcs-1+guard) << "\";" << endl;
		vhdl << tab << declare("S0", wcs + guard) << "<= " << zg(1, 0) << " & \"" << generateFixPointNumber(0.0, 0, wcs-1+guard) << "\";" << endl;
		vhdl << tab << declare("X0", wx  + guard) << "<= reducedX;" << endl;
		vhdl << tab << declare("D0") << "<= absX(" << w << ");" << endl;
		
		mpf_clear (xinit);
				
		//create the stages of micro-rotations
		int stage;
		
		wcs += guard;
		wx += guard;
		
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
		
		manageCriticalPath(target->localWireDelay(1+w+1) + target->adderDelay(1+w+1));
		syncCycleFromSignal(join("D", stage));
		
		vhdl << tab << declare("preRoundedIntCout", 1+w+1) << "<= " << join("C", stage) << "(" << wcs-1 << " downto " << guard-1 << ");" << endl;
		vhdl << tab << declare("preRoundedIntSout", 1+w+1) << "<= " << join("S", stage) << "(" << wcs-1 << " downto " << guard-1 << ");" << endl;
		
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



	std::string CordicSinCosClassic::generateFixPointNumber(float x, int wI, int wF)
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
	
	std::string CordicSinCosClassic::generateFixPointNumber(mpf_t x, int wI, int wF)
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
	
	mpz_class CordicSinCosClassic::fp2fix(mpfr_t x, int wI, int wF){
		mpz_class h;
		
		mpfr_mul_2si(x, x, wF, GMP_RNDD);
        mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDD);  
		
		return h;
	}

}











