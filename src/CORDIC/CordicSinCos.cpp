#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCos.hpp"

using namespace std;

namespace flopoco{

	CordicSinCos::CordicSinCos(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_), wI(0), wF(w_)
	{
		int wIxy = wI+2, wFxy = wF, wIz = wI+2, wFz = wF;
		int wInxy = wIxy + wFxy + 1;
		int wInz = wIz + wFz + 1;
		//TODO: verify the validity of the necessary guard bits
		int guardz, guardxy;
		
		//initialize testing random number
		gmp_randinit_mt (state);
		
		guardxy = ceil(log2(1 + wI + wF))+3;
		guardz  = guardxy;
		
		srcFileName="CordicSinCos";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCos_" << 1+wI+wF <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCos_" << 1+wI+wF << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		addInput  ( "X"  , 1+wI+wF, true );

		// declaring output
		addOutput  ( "C"  , 1+wI+wF, 2 );
		addOutput  ( "S"  , 1+wI+wF, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wInxy + ceil(log2(1 + wI + wF))));
		
		//create the Z0, Y0 and D0 signals for the first stage
		mpf_t xinit;
		
		mpf_set_default_prec (1+wI+wF+guardxy);
		mpf_init2   (xinit, 1+wI+wF+guardxy);
		mpf_set_str (xinit, "0.607252935008881256169446834473e0", 10);
		
		vhdl << tab << declare("C0", wInxy + guardxy) << "<= " << zg(3, 0) << " & \"" << generateFixPointNumber(xinit, wIxy-2, wFxy+guardxy) << "\";" << endl;
		vhdl << tab << declare("S0", wInxy + guardxy) << "<= " << zg(3, 0) << " & \"" << generateFixPointNumber(0.0, wIxy-2, wFxy+guardxy) << "\";" << endl;
		vhdl << tab << declare("X0", wInz  + guardxy) << "<= X(" << wI+wF << ") & X(" << wI+wF << ") & X & " << zg(guardxy, 0) << ";" << endl;
		vhdl << tab << declare("D0") << "<= X(" << wIz+wFz-2 << ");" << endl;
		
		mpf_clear (xinit);
		
		setCycleFromSignal("D0");
				
		//create the stages of micro-rotations
		int stage;
		
		wFxy += guardxy;
		wFz += guardz;
		
		REPORT(DEBUG, 1+wI+2+wF+guardxy-2 << " stages");

		for(stage=0; stage<1+wI+2+wF+guardxy-2; stage++){
			
			manageCriticalPath(target->localWireDelay(1+wIxy+wFxy) + target->lutDelay());
			
			//shift Xin and Yin with 2^n positions to the right
			if(stage==0){
				vhdl << tab << declare(join("CShift", stage), 1+wIxy+wFxy) << " <= C" << stage << ";" <<endl;
			}else{
				if(stage==1)
					vhdl << tab << declare(join("CsignExtend", stage), stage) << " <= C" << stage << "(" << 1+wIxy+wFxy-1 << ");" <<endl;
				else
					vhdl << tab << declare(join("CsignExtend", stage), stage) << " <= (others => C" << stage << "(" << 1+wIxy+wFxy-1 << "));" <<endl;
				vhdl << tab << declare(join("Cshifted", stage), 1+wIxy+wFxy-stage) << " <= C" << stage << range(1+wIxy+wFxy-1, stage) << ";" <<endl;
				vhdl << tab << declare(join("CShift", stage), 1+wIxy+wFxy) << " <= CsignExtend" << stage << " & Cshifted" << stage << ";" <<endl;
			}
			if(stage==0){
				vhdl << tab << declare(join("SShift", stage), 1+wIxy+wFxy) << " <= S" << stage << ";" <<endl;
			}else{
				if(stage==1)
					vhdl << tab << declare(join("SsignExtend", stage), stage) << " <= S" << stage << "(" << 1+wIxy+wFxy-1 << ");" <<endl;
				else
					vhdl << tab << declare(join("SsignExtend", stage), stage) << " <= (others => S" << stage << "(" << 1+wIxy+wFxy-1 << "));" <<endl;
				vhdl << tab << declare(join("Sshifted", stage), 1+wIxy+wFxy-stage) << " <= S" << stage << range(1+wIxy+wFxy-1, stage) << ";" <<endl;
				vhdl << tab << declare(join("SShift", stage), 1+wIxy+wFxy) << " <= SsignExtend" << stage << " & Sshifted" << stage << ";" <<endl;
			}
			vhdl << tab << declare(join("intC", (stage+1)), 1+wIxy+wFxy) << " <= " << join("C", stage) <<  " + " <<  join("SShift", stage) <<  " when D" << stage <<"='1'" <<  endl
			     << tab << tab << " else "  << join("C", stage) <<  " - " <<  join("SShift", stage) << ";" << endl;
			vhdl << tab << declare(join("intS", (stage+1)), 1+wIxy+wFxy) << " <= " << join("S", stage) <<  " + " <<  join("CShift", stage) <<  " when D" << stage <<" ='0'" <<  endl
			     << tab << tab << " else "  << join("S", stage) <<  " - " <<  join("CShift", stage) << ";" << endl;

			
			//create the constant signal for the arctan
			mpfr_t zatan, zpow2, constPi;
			
			mpfr_init2(zatan, 200);
			mpfr_init2(zpow2, 200);
			mpfr_init2(constPi, 200);
			
			mpfr_set_d(zpow2, 1, GMP_RNDD);
			mpfr_mul_2si(zpow2, zpow2, (-1)*stage, GMP_RNDD);
			mpfr_atan(zatan, zpow2, GMP_RNDD);
			mpfr_mul_2si(zatan, zatan, 1, GMP_RNDD);
			mpfr_const_pi( constPi, GMP_RNDD);
			mpfr_div(zatan, zatan, constPi, GMP_RNDD);
					
			//create the arctangent factor to be added to Zin
			bool negAtan = false;
			std::string strConverted;
			mpz_class fixConverted;
			
			fixConverted = fp2fix(zatan, wIxy, wFxy);			
			
			if(fixConverted<0){
				fixConverted = fixConverted * (-1);
				negAtan = true;
			}
			strConverted = unsignedBinary(fixConverted, 1+wIz+wFz-1);
			
			mpfr_free_cache();
			
			manageCriticalPath(target->localWireDelay(1+wIz+wFz) + target->lutDelay());
			
			// setCycleFromSignal(join("X",stage));
			
			vhdl << tab << declare(join("atan2PowStage", stage), 1+wIz+wFz) << " <= \'0\' & \"" << strConverted << "\";" <<endl;
			if(negAtan){
				vhdl << tab << declare(join("newAtan2PowStage", stage), 1+wIz+wFz) << " <= atan2PowStage" << stage << " xor (" << 1+wIz+wFz-1 << " downto 0 => D" << stage <<");" <<endl;
				vhdl << tab << declare(join("cInX", stage)) << "<= D" << stage << ";" <<endl;
			}
			else{
				vhdl << tab << declare(join("newAtan2PowStage", stage), 1+wIz+wFz) << " <= atan2PowStage" << stage << " xor (" << 1+wIz+wFz-1 << " downto 0 => (not D" << stage << "));" <<endl;
				vhdl << tab << declare(join("cInX", stage)) << "<= not D" << stage << ";" <<endl;
			}
			

			vhdl << tab << declare(join("intX", stage+1), 1+wIz+wFz) << " <= " << join("X", stage) <<  " + " <<  join("newAtan2PowStage", stage) <<  " + " <<  join("cInX", stage) << ";" << endl;
 			//create Dout as the result of the comparison of intZout with 0
			vhdl << tab << declare(join("intD", (stage+1))) << " <= intX" << stage+1 << "(" << 1+wIz+wFz-1 <<");" <<endl;

			
			//create the outputs
			vhdl << tab << declare(join("C", stage+1), 1+wIxy+wFxy) << " <= intC" << stage+1 << ";" <<endl;
			vhdl << tab << declare(join("S", stage+1), 1+wIxy+wFxy) << " <= intS" << stage+1 << ";" <<endl;
			vhdl << tab << declare(join("X", stage+1), 1+wIz-1+wFz) << " <= intX" << stage+1 << "(" << 1+wIz+wFz-2 << " downto 0);" <<endl;
			vhdl << tab << declare(join("D", stage+1)) << " <= intD" << stage+1 << ";" <<endl;
			
			//decrement the size of Z
			wIz--;
		}
		
		manageCriticalPath(target->localWireDelay(1+wIxy+wFxy));
		
		vhdl << tab << declare("preRoundedIntCout", 1+wI+wF+1) << "<= " << join("C", stage) << "(" << wIxy-2+wFxy << " downto " << guardxy-1 << ");" << endl;
		vhdl << tab << declare("preRoundedIntSout", 1+wI+wF+1) << "<= " << join("S", stage) << "(" << wIxy-2+wFxy << " downto " << guardxy-1 << ");" << endl;
		
		vhdl << tab << declare("roundedIntCout", 1+wI+wF+1) << "<= preRoundedIntCout "
															<< "+"
															<< " (" << zg(1+wI+wF, 0) << " & \'1\');" << endl;
		vhdl << tab << declare("roundedIntSout", 1+wI+wF+1) << "<= preRoundedIntSout "
															<< "+"
															<< " (" << zg(1+wI+wF, 0) << " & \'1\');" << endl;
															
		//		setCycleFromSignal("roundedIntSout");
		//   syncCycleFromSignal("roundedIntCout");
		
		//assign output
		// manageCriticalPath(target->localWireDelay(1+wI+wF));
		
		vhdl << tab << "C" << "<= roundedIntCout(" << 1+wI+wF << " downto 1);" << endl;
		vhdl << tab << "S" << "<= roundedIntSout(" << 1+wI+wF << " downto 1);" << endl;
	};


	void CordicSinCos::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpfr_t z, constPi, rsin, rcos;
		
		/* Compute correct value */
		mpfr_init2(z, 1+wI+wF);
		mpfr_init2(rsin, 1+wI+wF); 
		mpfr_init2(rcos, 1+wI+wF); 
		
		mpfr_init2(constPi, 4*(1+wI+wF));
		
		mpz_class rsin_z, rcos_z;
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDD); // this rounding is exact
		mpfr_div_2si (z, z, wF, GMP_RNDD); // this rounding is acually exact
		
		mpfr_mul_2si(z, z, -1, GMP_RNDD);

		mpfr_const_pi( constPi, GMP_RNDD);
		mpfr_mul(z, z, constPi, GMP_RNDD);
		

		// Rounding down
		mpfr_sin(rsin, z, GMP_RNDD); 
		mpfr_cos(rcos, z, GMP_RNDD);
		
		mpfr_mul_2si (rsin, rsin, wF, GMP_RNDD); // exact rnd here
		mpfr_get_z (rsin_z.get_mpz_t(), rsin, GMP_RNDD); // there can be a real rounding here
		mpfr_mul_2si (rcos, rcos, wF, GMP_RNDD); // exact rnd here
		mpfr_get_z (rcos_z.get_mpz_t(), rcos, GMP_RNDD); // there can be a real rounding here

		// Set outputs 
		tc->addExpectedOutput ("C", rcos_z);
		tc->addExpectedOutput ("S", rsin_z);

		// Rounding up
		mpfr_sin(rsin, z, GMP_RNDU); 
		mpfr_cos(rcos, z, GMP_RNDU);
		
		mpfr_mul_2si (rsin, rsin, wF, GMP_RNDD); // exact rnd here
		mpfr_get_z (rsin_z.get_mpz_t(), rsin, GMP_RNDU); // should be exact, too
		mpfr_mul_2si (rcos, rcos, wF, GMP_RNDD); // exact rnd here
		mpfr_get_z (rcos_z.get_mpz_t(), rcos, GMP_RNDU); // there can be a real rounding here

		// Set outputs 
		tc->addExpectedOutput ("C", rcos_z);
		tc->addExpectedOutput ("S", rsin_z);

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
		
		//mpf_set_default_prec (1+wI+wF+guardxy);
		
		mpfr_init2(z, 1+wI+wF+ceil(log2(1 + wI + wF))+3);
		mpz_init2 (z_z, 1+wI+wF+ceil(log2(1 + wI + wF))+3);
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+wI+wF+ceil(log2(1 + wI + wF))+3);
		//mpf_set_str (zinit, "1.5707963267949e0", 10);
		mpf_set_str (zinit, "1.0e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, wF-1, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+wI+wF+ceil(log2(1 + wI + wF))+3);
		//mpf_set_str (zinit, "0.5235987755983e0", 10);
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, wF-1, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+wI+wF+ceil(log2(1 + wI + wF))+3);
		//mpf_set_str (zinit, "0.78539816339745e0", 10);
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, wF-1, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+wI+2+wF+ceil(log2(1 + wI + wF))+3);
		//mpf_set_str (zinit, "1.0471975511966e0", 10);
		mpf_set_str (zinit, "0.66666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD);
		
		mpfr_mul_2si (z, z, wF-1, GMP_RNDD); 
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
		
		mpfr_init2 (mx, 1+wI+wF);
		
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












#if 0
			setCycleFromSignal(join("CShift", stage));
			syncCycleFromSignal(join("SShift", stage));
			manageCriticalPath(target->localWireDelay(1+wFxy) + target->lutDelay());
			
			//complement the shifted Xin and Yin if necessary
			vhdl << tab << declare(join("newCShift", stage), 1+wIxy+wFxy) << " <= CShift" << stage << " xor (" << 1+wIxy+wFxy-1 << " downto 0 => D" << stage << ");" <<endl;
			vhdl << tab << declare(join("newSShift", stage), 1+wIxy+wFxy) << " <= SShift" << stage << " xor (" << 1+wIxy+wFxy-1 << " downto 0 => (not D" << stage << "));" <<endl;
			
			//compute the carry-ins
			vhdl << tab << declare(join("cInNewC", stage)) << "<= D" << stage << ";" <<endl;
			vhdl << tab << declare(join("cInNewS", stage)) << "<= not D" << stage << ";" <<endl;
			
			#if 0
			nextCycle();
			
			//create Xout and Yout
			IntAdder *mainAdder = new IntAdder(target, 1+wIxy+wFxy, inDelayMap("X",getCriticalPath()));
			oplist.push_back(mainAdder);
			
			inPortMap(mainAdder, "X", join("C", stage));
			inPortMap(mainAdder, "Y", join("newSShift", stage));
			inPortMap(mainAdder, "Cin", join("cInNewS", stage));
			outPortMap (mainAdder, "R", join("intC", (stage+1)));
			vhdl << instance(mainAdder, join("cAdder", stage)) << endl;
			
			inPortMap(mainAdder, "X", join("S", stage));
			inPortMap(mainAdder, "Y", join("newCShift", stage));
			inPortMap(mainAdder, "Cin", join("cInNewC", stage));
			outPortMap (mainAdder, "R", join("intS", (stage+1)));
			vhdl << instance(mainAdder, join("sAdder", stage)) << endl;
			
			setCycleFromSignal(join("intC", (stage+1)));
			syncCycleFromSignal(join("intS", (stage+1)));
			setCriticalPath(mainAdder->getOutputDelay("R"));
			double xyPath = getCriticalPath();
			
#else
			manageCriticalPath(target->localWireDelay() + target->adderDelay(1+wIxy+wFxy));
			vhdl << tab << declare(join("intC", (stage+1)), 1+wIxy+wFxy) << " <= " << join("C", stage) <<  " + " <<  join("newSShift", stage) <<  " + " <<  join("cInNewC", stage) << ";" << endl;
			vhdl << tab << declare(join("intS", (stage+1)), 1+wIxy+wFxy) << " <= " << join("S", stage) <<  " + " <<  join("newCShift", stage) <<  " + " <<  join("cInNewS", stage) << ";" << endl;
#endif
#endif





#if 0
			syncCycleFromSignal(join("newAtan2PowStage", stage));
			nextCycle();
			
			//create Zout
			mainAdder = new IntAdder(target, 1+wIz+wFz, inDelayMap("X",getCriticalPath()));
			oplist.push_back(mainAdder);
			
			inPortMap(mainAdder, "X", join("X", stage));
			inPortMap(mainAdder, "Y", join("newAtan2PowStage", stage));
			inPortMap(mainAdder, "Cin", join("cInX", stage));
			outPortMap (mainAdder, "R", join("intX", (stage+1)));
			vhdl << instance(mainAdder, join("xAdder", stage)) << endl;
			
			setCycleFromSignal(join("intX", (stage+1)));
			setCriticalPath(mainAdder->getOutputDelay("R"));
			manageCriticalPath(target->localWireDelay(1+wIz+wFz));
			
 			//create Dout as the result of the comparison of intZout with 0
			vhdl << tab << declare(join("intD", (stage+1))) << " <= intX" << stage+1 << "(" << 1+wIz+wFz-1 <<");" <<endl;
			
			double zPath = getCriticalPath();
			
			if (getCycleFromSignal(join("intC", (stage+1))) == getCycleFromSignal(join("intD", (stage+1))))
				setCriticalPath(max(zPath, getCriticalPath()));
			else
				if (syncCycleFromSignal(join("intC", (stage+1))))
					setCriticalPath(xyPath);
#endif
