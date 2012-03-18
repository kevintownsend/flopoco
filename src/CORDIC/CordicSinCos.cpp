#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCos.hpp"
#include "../ConstMult/FixRealKCM.hpp"

using namespace std;

namespace flopoco{

	CordicSinCos::CordicSinCos(Target* target, int wI_, int wF_, map<string, double> inputDelays) 
		: Operator(target), wI(wI_), wF(wF_)
	{
		int wIn = wI + 2*wF + 1;
		
		srcFileName="CordicSinCos";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCos_" << wIn <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCos_" << wIn << "_uid" << getNewUId();
		setName( name.str() );

		// declaring inputs
		addInput  ( "Z"  , 1+wI+wF, true );

		// declaring output
		addOutput  ( "Xout"  , 1+wI+wF, true );
		addOutput  ( "Yout"  , 1+wI+wF, true );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wIn) + target->lutDelay());
		
		//create the Z0, Y0 and D0 signals for the first stage
		vhdl << tab << declare("X0", wIn) << "<= \'0\' & \"" << generateFixPointNumber(1, wI, 2*wF) << "\";" << endl;
		vhdl << tab << declare("Y0", wIn) << "<= \'0\' & \"" << generateFixPointNumber(0, wI, 2*wF) << "\";" << endl;
		vhdl << tab << declare("Z0", wIn) << "<= Z & " << zg(wF, 0) << ";" << endl;
		vhdl << tab << declare("D0") << "<= Z(" << wI+wF << ");" << endl;
		
		setCycleFromSignal("D0");
				
		//create the wIn-1 stages of micro-rotations
		int stage;
		FixMicroRotation *microRotation;
		
		for(stage=0; stage<wIn-2; stage++){
			manageCriticalPath(target->localWireDelay(wIn));
			
			microRotation = new FixMicroRotation(target, wI, 2*wF, stage, inDelayMap("Xin",getCriticalPath()));
			oplist.push_back(microRotation);
			
			inPortMap(microRotation, "Xin", getParamName("X", stage));
			inPortMap(microRotation, "Yin", getParamName("Y", stage));
			inPortMap(microRotation, "Zin", getParamName("Z", stage));
			inPortMap(microRotation, "Din", getParamName("D", stage));
			outPortMap(microRotation, "Xout", getParamName("X", stage+1));
			outPortMap(microRotation, "Yout", getParamName("Y", stage+1));
			outPortMap(microRotation, "Zout", getParamName("Z", stage+1));
			outPortMap(microRotation, "Dout", getParamName("D", stage+1));
			vhdl << instance(microRotation, getParamName("microRotation", stage)) << endl;
			
			setCycleFromSignal(getParamName("X", stage+1));
			syncCycleFromSignal(getParamName("Y", stage+1));
			syncCycleFromSignal(getParamName("Z", stage+1));
			syncCycleFromSignal(getParamName("D", stage+1));
			setCriticalPath(microRotation->getOutputDelay("Dout"));
		}
		
		//divide by the constant to obtain the actual value of the sine and cosine
		FixRealKCM *constMultiplier = new FixRealKCM(target, -(2*wF), wI, 1, -(2*wF), string("0.607252935008881256169446834473"), 1.0, inDelayMap("X",getCriticalPath()));
		oplist.push_back(constMultiplier);
		
		//multiply x
		manageCriticalPath(target->localWireDelay(wIn+1) + target->lutDelay());
		
		vhdl << tab << declare("preXout", wIn+1) << "<= " << getParamName("X", stage) << " & '0';" << endl;
		inPortMap (constMultiplier, "X", "preXout");
		outPortMap(constMultiplier, "R", "intXout");
		vhdl << instance(constMultiplier, "constMultiplierX") << endl;
		
		//multiply y
		vhdl << tab << declare("preYout", wIn+1) << "<= " << getParamName("Y", stage) << " & '0';" << endl;
		inPortMap (constMultiplier, "X", "preYout");
		outPortMap(constMultiplier, "R", "intYout");
		vhdl << instance(constMultiplier, "constMultiplierY") << endl;
		
		setCriticalPath(constMultiplier->getOutputDelay("R"));
		setCycleFromSignal("intYout");
		syncCycleFromSignal("intXout");
		
		//perform rounding of the final result
		manageCriticalPath(target->localWireDelay(1+wI+wF+1) + target->lutDelay() + target->adderDelay(1+wI+wF+1));
		
		vhdl << tab << declare("roundedIntXout", 1+wI+wF+1) << "<= intXout(" << wIn << " downto " << wF << ") "
															<< "+"
															<< " (" << zg(1+wI+wF, 0) << " & \'1\');" << endl;
		vhdl << tab << declare("roundedIntYout", 1+wI+wF+1) << "<= intYout(" << wIn << " downto " << wF << ") "
															<< "+"
															<< " (" << zg(1+wI+wF, 0) << " & \'1\');" << endl;
															
		setCycleFromSignal("roundedIntXout");
		syncCycleFromSignal("roundedIntYout");
		
		//assign output
		manageCriticalPath(target->localWireDelay(1+wI+wF));
		
		vhdl << tab << "Xout" << "<= roundedIntXout(" << 1+wI+wF << " downto 1);" << endl;
		vhdl << tab << "Yout" << "<= roundedIntYout(" << 1+wI+wF << " downto 1);" << endl;
	};


	void CordicSinCos::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("Z");
		mpfr_t z, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		
		/* Compute correct value */
		mpfr_init2(z, 1+wI+wF);
		mpfr_init2(rsin, 1+wI+wF); 
		mpfr_init2(rcos, 1+wI+wF); 
		mpz_init2 (rsin_z, 1+wI+wF);
		mpz_init2 (rcos_z, 1+wI+wF);
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, wF, GMP_RNDN); // this rounding is acually exact
		
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		
		mpfr_mul_2si (rsin, rsin, wF, GMP_RNDN); // exact rnd here
		mpfr_get_z (rsin_z, rsin, GMP_RNDN); // there can be a real rounding here
		mpfr_mul_2si (rcos, rcos, wF, GMP_RNDN); // exact rnd here
		mpfr_get_z (rcos_z, rcos, GMP_RNDN); // there can be a real rounding here

		// Set outputs 
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		tc->addExpectedOutput ("Xout", cos_zc);
		tc->addExpectedOutput ("Yout", sin_zc);

		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
	}


	void CordicSinCos::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpfr_t z;
		mpz_t z_z;
		
		mpfr_init2(z, 1+wI+wF);
		mpz_init2 (z_z, 1+wI+wF);
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("Z",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2
		tc = new TestCase (this);
		mpfr_set_d (z, 1.5707963267949, GMP_RNDN); 
		mpfr_mul_2si (z, z, wF, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("Z",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6
		tc = new TestCase (this);
		mpfr_set_d (z, 0.5235987755983, GMP_RNDN); 
		mpfr_mul_2si (z, z, wF, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("Z",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4
		tc = new TestCase (this);
		mpfr_set_d (z, 0.78539816339745, GMP_RNDN); 
		mpfr_mul_2si (z, z, wF, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("Z",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3
		tc = new TestCase (this);
		mpfr_set_d (z, 1.0471975511966, GMP_RNDN); 
		mpfr_mul_2si (z, z, wF, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("Z",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		mpfr_clears (z, NULL);
	}

	//not very random for the moment
	//returns the case z=0
	TestCase* CordicSinCos::buildRandomTestCase(int i) 
	{
		TestCase* tc = new TestCase(this);
		
		tc = new TestCase (this);
		tc -> addInput ("Z",mpz_class(0));
		emulate(tc);
		
		return tc;
	}
	
	std::string CordicSinCos::generateFixPointNumber(float x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		float xcopy = x;
		mpfr_t mx;
		mpz_class h;
		
		mpfr_init(mx);
		
		if(xcopy<0){
			xcopy = xcopy * (-1);
		}
		
		mpfr_set_d(mx, xcopy, GMP_RNDN);
		mpfr_mul_2si(mx, mx, wF, GMP_RNDN);
		
		mpfr_get_z(h.get_mpz_t(), mx,  GMP_RNDN); 
        
        result = unsignedBinary(h, size);
        
        return result;
	}
	
	std::string CordicSinCos::getParamName(std::string s, int number)
	{
		std::stringstream aux;
		
		aux << s << number;
		
		return aux.str();
		
	}

}












