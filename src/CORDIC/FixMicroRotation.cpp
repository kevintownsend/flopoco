#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixMicroRotation.hpp"

using namespace std;

namespace flopoco{

	FixMicroRotation::FixMicroRotation(Target* target, int wI_, int wF_, int stage_, map<string, double> inputDelays) 
		: Operator(target), wI(wI_), wF(wF_), stage(stage_)
	{
		int wIn = wI + wF + 1;
		
		srcFileName="FixMicroRotation";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "FixMicroRotation_" << wIn <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId() << "_stage" << stage;
		else
			name << "FixMicroRotation_" << wIn << "_uid" << getNewUId() << "_stage" << stage;
		setName( name.str() );

		// declaring inputs
		addInput  ( "Xin"  , wIn, true );
		addInput  ( "Yin"  , wIn, true );
		addInput  ( "Zin"  , wIn, true );
		addInput  ( "Din" 			   );

		// declaring output
		addOutput  ( "Xout"  , wIn, true );
		addOutput  ( "Yout"  , wIn, true );
		addOutput  ( "Zout"  , wIn, true );
		addOutput  ( "Dout"  			 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wIn) + 2*target->lutDelay());
		
		//shift Xin and Yin with 2^n positions to the right
		if(stage==0){
			vhdl << tab << declare("XinShift", wIn) << " <= Xin;" <<endl;
		}else{
			if(stage==1)
				vhdl << tab << declare("Xin_signExtend", stage) << " <= Xin(" << wIn-1 << ");" <<endl;
			else
				vhdl << tab << declare("Xin_signExtend", stage) << " <= (others => Xin(" << wIn-1 << "));" <<endl;
			vhdl << tab << declare("Xin_shifted", wIn-stage) << " <= Xin" << range(wIn-1, stage) << ";" <<endl;
			vhdl << tab << declare("XinShift", wIn) << " <= Xin_signExtend & Xin_shifted;" <<endl;
		}
		if(stage==0){
			vhdl << tab << declare("YinShift", wIn) << " <= Yin;" <<endl;
		}else{
			if(stage==1)
				vhdl << tab << declare("Yin_signExtend", stage) << " <= Yin(" << wIn-1 << ");" <<endl;
			else
				vhdl << tab << declare("Yin_signExtend", stage) << " <= (others => Yin(" << wIn-1 << "));" <<endl;
			vhdl << tab << declare("Yin_shifted", wIn-stage) << " <= Yin" << range(wIn-1, stage) << ";" <<endl;
			vhdl << tab << declare("YinShift", wIn) << " <= Yin_signExtend & Yin_shifted;" <<endl;
		}
		
		setCycleFromSignal("XinShift");
		syncCycleFromSignal("YinShift");
		manageCriticalPath(target->localWireDelay(wIn) + target->lutDelay());
		
		//complement the shifted Xin and Yin if necessary
		vhdl << tab << declare("newXinShift", wIn) << " <= XinShift xor (" << wIn-1 << " downto 0 => Din);" <<endl;
		vhdl << tab << declare("newYinShift", wIn) << " <= YinShift xor (" << wIn-1 << " downto 0 => (not Din));" <<endl;
		
		//compute the carry-ins
		vhdl << tab << declare("cInNewX") << "<= Din;" <<endl;
		vhdl << tab << declare("cInNewY") << "<= not Din;" <<endl;
		
		nextCycle();
		
		//create Xout and Yout
		IntAdder *mainAdder = new IntAdder(target, wIn, inDelayMap("X",getCriticalPath()));
		oplist.push_back(mainAdder);
		
		inPortMap(mainAdder, "X", "Xin");
		inPortMap(mainAdder, "Y", "newYinShift");
		inPortMap(mainAdder, "Cin", "cInNewY");
		outPortMap (mainAdder, "R", "intXout");
		vhdl << instance(mainAdder, "xAdder") << endl;
		
		inPortMap(mainAdder, "X", "Yin");
		inPortMap(mainAdder, "Y", "newXinShift");
		inPortMap(mainAdder, "Cin", "cInNewX");
		outPortMap (mainAdder, "R", "intYout");
		vhdl << instance(mainAdder, "yAdder") << endl;
		
		setCycleFromSignal("intXout");
		syncCycleFromSignal("intYout");
		setCriticalPath(mainAdder->getOutputDelay("R"));
		
		double xyPath = getCriticalPath();
		
		//create the constant signal for the arctan
		mpfr_t zatan, zpow2;
		
		mpfr_init(zatan);
		mpfr_init(zpow2);
		
		mpfr_set_d(zpow2, 1, GMP_RNDN);
		mpfr_mul_2si(zpow2, zpow2, (-1)*stage, GMP_RNDN);
		mpfr_atan(zatan, zpow2, GMP_RNDN);
				
		//create the arctangent factor to be added to Zin
		bool negAtan = false;
		std::string strConverted;
		mpz_class fixConverted = fp2fix(zatan, wI, wF);
		
		if(fixConverted<0){
			fixConverted = fixConverted * (-1);
			negAtan = true;
		}
		strConverted = unsignedBinary(fixConverted, wIn-1);
		
		manageCriticalPath(target->localWireDelay(wIn) + 2*target->lutDelay());
		
		setCycleFromSignal("Zin");
		
		vhdl << tab << declare("atan2PowStage", wIn) << " <= \'0\' & \"" << strConverted << "\";" <<endl;
		if(negAtan){
			vhdl << tab << declare("newAtan2PowStage", wIn) << " <= atan2PowStage xor (" << wIn-1 << " downto 0 => Din);" <<endl;
			vhdl << tab << declare("cInZ") << "<= Din;" <<endl;
		}
		else{
			vhdl << tab << declare("newAtan2PowStage", wIn) << " <= atan2PowStage xor (" << wIn-1 << " downto 0 => (not Din));" <<endl;
			vhdl << tab << declare("cInZ") << "<= not Din;" <<endl;
		}
		
		syncCycleFromSignal("newAtan2PowStage");
		nextCycle();
		
		//create Zout
		inPortMap(mainAdder, "X", "Zin");
		inPortMap(mainAdder, "Y", "newAtan2PowStage");
		inPortMap(mainAdder, "Cin", "cInZ");
		outPortMap (mainAdder, "R", "intZout");
		vhdl << instance(mainAdder, "zAdder") << endl;
		
		setCycleFromSignal("intZout");
		setCriticalPath(mainAdder->getOutputDelay("R"));
		manageCriticalPath(target->localWireDelay(wIn) + target->lutDelay());
		
		//create Dout as the result of the comparison of intZout with 0
		vhdl << tab << declare("intDout") << " <= intZout(" << wIn-1 <<");" <<endl;
		
		double zPath = getCriticalPath();
		
		if (getCycleFromSignal("intXout") == getCycleFromSignal("intDout"))
			setCriticalPath(max(zPath, getCriticalPath()));
		else
			if (syncCycleFromSignal("intXout"))
				setCriticalPath(xyPath);
		
		//create the outputs
		vhdl << tab << "Xout" << " <= intXout;" <<endl;
		vhdl << tab << "Yout" << " <= intYout;" <<endl;
		vhdl << tab << "Zout" << " <= intZout;" <<endl;
		vhdl << tab << "Dout" << " <= intDout;" <<endl;
	};


	void FixMicroRotation::emulate(TestCase * tc) 
	{
		
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("Xin");
		mpz_class svY = tc->getInputValue("Yin");
		mpz_class svZ = tc->getInputValue("Zin");
		mpz_class svD = tc->getInputValue("Din");
		
		mpfr_t x, y, z, d, xint, yint, zint, rx, ry, rz, rd;
		mpfr_init2(x, 1+wI+wF);
		mpfr_init2(y, 1+wI+wF);
		mpfr_init2(z, 1+wI+wF);
		mpfr_init2(xint, 1+wI+wF);
		mpfr_init2(yint, 1+wI+wF);
		mpfr_init2(zint, 1+wI+wF);
		mpfr_init2(d, 1+wI+wF);
		mpfr_init2(rx, 1+wI+wF);
		mpfr_init2(ry, 1+wI+wF);
		mpfr_init2(rz, 1+wI+wF);
		mpfr_init2(rd, 1+wI+wF);
		
		mpz_t rx_z, ry_z, rz_z, rd_z;
		mpz_init2 (rx_z, 1+wI+wF);
		mpz_init2 (ry_z, 1+wI+wF);
		mpz_init2 (rz_z, 1+wI+wF);
		mpz_init2 (rd_z, 1+wI+wF);
	
		/* Compute correct value */
		mpfr_set_z (x, svX.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (x, x, wF, GMP_RNDN); // this rounding is acually exact
		
		mpfr_set_z (y, svY.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (y, y, wF, GMP_RNDN); // this rounding is acually exact
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, wF, GMP_RNDN); // this rounding is acually exact
		
		if(svD==0)
			mpfr_set_d (d, 1, GMP_RNDN); 
		else
			mpfr_set_d (d, -1, GMP_RNDN);
		mpfr_div_2si (d, d, wF, GMP_RNDN); // this rounding is acually exact
		
		mpfr_mul_2si (xint, x, stage*(-1), GMP_RNDN);
		mpfr_mul_2si (yint, y, stage*(-1), GMP_RNDN);
		
		mpfr_set_d   (zint, 1, GMP_RNDN); 
		mpfr_mul_2si (zint, zint, stage*(-1), GMP_RNDN);
		mpfr_atan    (zint, zint, GMP_RNDN);
		
		mpfr_mul (xint, xint, d, GMP_RNDN);
		mpfr_mul (yint, yint, d, GMP_RNDN);
		mpfr_mul (zint, zint, d, GMP_RNDN);
		
		mpfr_sub (rx, x, xint, GMP_RNDN);
		mpfr_add (ry, y, yint, GMP_RNDN);
		mpfr_sub (rz, z, zint, GMP_RNDN);
		
		if(mpfr_cmp_ui(rz, 0)<0)
			mpfr_set_d (rd, 1.0, GMP_RNDN); 
		else
			mpfr_set_d (rd, 0.0, GMP_RNDN); 

		// Set outputs
		mpfr_mul_2si (rx, rx, wF, GMP_RNDN); // exact rnd here
		mpfr_get_z (rx_z, rx, GMP_RNDN); // there can be a real rounding here 
		mpz_class rx_zc(rx_z);
		tc->addExpectedOutput ("Xout", rx_zc);
		
		mpfr_mul_2si (ry, ry, wF, GMP_RNDN); // exact rnd here
		mpfr_get_z (ry_z, ry, GMP_RNDN); // there can be a real rounding here 
		mpz_class ry_zc(ry_z);
		tc->addExpectedOutput ("Yout", ry_zc);
		
		mpfr_mul_2si (rz, rz, wF, GMP_RNDN); // exact rnd here
		mpfr_get_z (rz_z, rz, GMP_RNDN); // there can be a real rounding here 
		mpz_class rz_zc(rz_z);
		tc->addExpectedOutput ("Zout", rz_zc);
		
		//mpfr_mul_2si (rd, rd, wF, GMP_RNDN); // exact rnd here
		mpfr_get_z (rd_z, rd, GMP_RNDN); // there can be a real rounding here 
		mpz_class rd_zc(rd_z);
		tc->addExpectedOutput ("Dout", rd_zc);
	}


	void FixMicroRotation::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpfr_t z;
		mpz_t z_z;
		
		mpfr_init2(z, 1+wI+wF);
		mpz_init2 (z_z, 1+wI+wF);
		
		
		tc = new TestCase (this);
		tc -> addInput ("Xin",mpz_class(1));
		tc -> addInput ("Yin",mpz_class(0));
		tc -> addInput ("Zin",mpz_class(0));
		tc -> addInput ("Din",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		
		tc = new TestCase (this);
		tc -> addInput ("Xin",mpz_class(1));
		tc -> addInput ("Yin",mpz_class(0));
		tc -> addInput ("Din",mpz_class(0));
		
		mpfr_set_d (z, 1.5707963267949, GMP_RNDN); 
		mpfr_mul_2si (z, z, wF, GMP_RNDN); 
		mpfr_get_z (z_z, z, GMP_RNDN);  
		tc -> addInput ("Zin",mpz_class(z_z));
		
		emulate(tc);
		tcl->add(tc);
		
		mpfr_clears (z, NULL);
	}

	TestCase* FixMicroRotation::buildRandomTestCase(int i) {
		TestCase* tc = new TestCase(this);
		
		tc = new TestCase (this);
		tc -> addInput ("Xin",mpz_class(1));
		tc -> addInput ("Yin",mpz_class(0));
		tc -> addInput ("Zin",mpz_class(0));
		tc -> addInput ("Din",mpz_class(1));
		emulate(tc);
		
		return tc;
	}
	
		
	mpz_class FixMicroRotation::fp2fix(mpfr_t x, int wI, int wF){
		mpz_class h;
		
		mpfr_mul_2si(x, x, wF, GMP_RNDN);
        mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDN);  
		
		return h;
	}

}












