#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixMicroRotation.hpp"

using namespace std;

namespace flopoco{

	FixMicroRotation::FixMicroRotation(Target* target, int wIx_, int wFx_, int wIy_, int wFy_, int wIz_, int wFz_, int stage_, bool hasLargerZout_, map<string, double> inputDelays) 
		: Operator(target), wIx(wIx_), wFx(wFx_), wIy(wIy_), wFy(wFy_), wIz(wIz_), wFz(wFz_), stage(stage_), hasLargerZout(hasLargerZout_)
	{
		int wInx = wIx + wFx + 1;
		int wIny = wIy + wFy + 1;
		int wInz = wIz + wFz + 1;
		
		srcFileName="FixMicroRotation";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "FixMicroRotation_" << wInx <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId() << "_stage" << stage;
		else
			name << "FixMicroRotation_" << wInx << "_uid" << getNewUId() << "_stage" << stage;
		setName( name.str() );

		// declaring inputs
		addInput  ( "Xin"  , wInx, true );
		addInput  ( "Yin"  , wIny, true );
		addInput  ( "Zin"  , wInz, true );
		addInput  ( "Din" 			    );

		// declaring output
		addOutput  ( "Xout"  , wInx+1, true );
		addOutput  ( "Yout"  , wIny+1, true );
		if(hasLargerZout)
			addOutput  ( "Zout"  , wInz+1, true );
		else
			addOutput  ( "Zout"  , wInz, true );
		addOutput  ( "Dout"  			    );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wInx+1) + 2*target->lutDelay());
		
		//shift Xin and Yin with 2^n positions to the right
		if(stage==0){
			vhdl << tab << declare("XinShift", wInx+1) << " <= Xin & \'0\';" <<endl;
		}else{
			if(stage==1)
				vhdl << tab << declare("Xin_signExtend", stage) << " <= Xin(" << wInx-1 << ");" <<endl;
			else
				vhdl << tab << declare("Xin_signExtend", stage) << " <= (others => Xin(" << wInx-1 << "));" <<endl;
			vhdl << tab << declare("Xin_shifted", wInx-stage+1) << " <= Xin" << range(wInx-1, stage-1) << ";" <<endl;
			vhdl << tab << declare("XinShift", wInx+1) << " <= Xin_signExtend & Xin_shifted;" <<endl;
		}
		if(stage==0){
			vhdl << tab << declare("YinShift", wIny+1) << " <= Yin & \'0\';" <<endl;
		}else{
			if(stage==1)
				vhdl << tab << declare("Yin_signExtend", stage) << " <= Yin(" << wIny-1 << ");" <<endl;
			else
				vhdl << tab << declare("Yin_signExtend", stage) << " <= (others => Yin(" << wIny-1 << "));" <<endl;
			vhdl << tab << declare("Yin_shifted", wIny-stage+1) << " <= Yin" << range(wIny-1, stage-1) << ";" <<endl;
			vhdl << tab << declare("YinShift", wIny+1) << " <= Yin_signExtend & Yin_shifted;" <<endl;
		}
		
		vhdl << tab << declare("XinExtended", wInx+1) << " <= Xin & \'0\';" <<endl;
		vhdl << tab << declare("YinExtended", wIny+1) << " <= Yin & \'0\';" <<endl;
		
		setCycleFromSignal("XinShift");
		syncCycleFromSignal("YinShift");
		manageCriticalPath(target->localWireDelay(wInx+1) + target->lutDelay());
		
		//complement the shifted Xin and Yin if necessary
		vhdl << tab << declare("newXinShift", wInx+1) << " <= XinShift xor (" << wInx << " downto 0 => Din);" <<endl;
		vhdl << tab << declare("newYinShift", wIny+1) << " <= YinShift xor (" << wIny << " downto 0 => (not Din));" <<endl;
		
		//compute the carry-ins
		vhdl << tab << declare("cInNewX") << "<= Din;" <<endl;
		vhdl << tab << declare("cInNewY") << "<= not Din;" <<endl;
		
		nextCycle();
		
		//create Xout and Yout
		// TODO - should create two adders, to be fair, but can save up space
		// for the time being, as x and y have the same precisions
		IntAdder *mainAdder = new IntAdder(target, wInx+1, inDelayMap("X",getCriticalPath()));
		oplist.push_back(mainAdder);
		
		inPortMap(mainAdder, "X", "XinExtended");
		inPortMap(mainAdder, "Y", "newYinShift");
		inPortMap(mainAdder, "Cin", "cInNewY");
		outPortMap (mainAdder, "R", "intXout");
		vhdl << instance(mainAdder, "xAdder") << endl;
		
		inPortMap(mainAdder, "X", "YinExtended");
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
		mpz_class fixConverted;
		
		if(hasLargerZout)
			fixConverted = fp2fix(zatan, wIz, wFz+1);
		else
			fixConverted = fp2fix(zatan, wIz, wFz);
		
		if(fixConverted<0){
			fixConverted = fixConverted * (-1);
			negAtan = true;
		}
		if(hasLargerZout)
			strConverted = unsignedBinary(fixConverted, wInz);
		else
			strConverted = unsignedBinary(fixConverted, wInz-1);
		
		manageCriticalPath(target->localWireDelay(wInz+1) + 2*target->lutDelay());
		
		setCycleFromSignal("Zin");
		
		if(hasLargerZout){
			vhdl << tab << declare("atan2PowStage", wInz+1) << " <= \'0\' & \"" << strConverted << "\";" <<endl;
			if(negAtan){
				vhdl << tab << declare("newAtan2PowStage", wInz+1) << " <= atan2PowStage xor (" << wInz << " downto 0 => Din);" <<endl;
				vhdl << tab << declare("cInZ") << "<= Din;" <<endl;
			}
			else{
				vhdl << tab << declare("newAtan2PowStage", wInz+1) << " <= atan2PowStage xor (" << wInz << " downto 0 => (not Din));" <<endl;
				vhdl << tab << declare("cInZ") << "<= not Din;" <<endl;
			}
		}else{
			vhdl << tab << declare("atan2PowStage", wInz) << " <= \'0\' & \"" << strConverted << "\";" <<endl;
			if(negAtan){
				vhdl << tab << declare("newAtan2PowStage", wInz) << " <= atan2PowStage xor (" << wInz-1 << " downto 0 => Din);" <<endl;
				vhdl << tab << declare("cInZ") << "<= Din;" <<endl;
			}
			else{
				vhdl << tab << declare("newAtan2PowStage", wInz) << " <= atan2PowStage xor (" << wInz-1 << " downto 0 => (not Din));" <<endl;
				vhdl << tab << declare("cInZ") << "<= not Din;" <<endl;
			}
		}
		
		if(hasLargerZout)
			vhdl << tab << declare("ZinExtended", wInz+1) << "<= Zin & \'0\';" <<endl;
		else
			vhdl << tab << declare("ZinExtended", wInz) << "<= Zin;" <<endl;
		
		syncCycleFromSignal("newAtan2PowStage");
		nextCycle();
		
		//create Zout
		if(hasLargerZout)
			mainAdder = new IntAdder(target, wInz+1, inDelayMap("X",getCriticalPath()));
		else
			mainAdder = new IntAdder(target, wInz, inDelayMap("X",getCriticalPath()));
		oplist.push_back(mainAdder);
		
		inPortMap(mainAdder, "X", "ZinExtended");
		inPortMap(mainAdder, "Y", "newAtan2PowStage");
		inPortMap(mainAdder, "Cin", "cInZ");
		outPortMap (mainAdder, "R", "intZout");
		vhdl << instance(mainAdder, "zAdder") << endl;
		
		setCycleFromSignal("intZout");
		setCriticalPath(mainAdder->getOutputDelay("R"));
		manageCriticalPath(target->localWireDelay(wInz) + target->lutDelay());
		
		//create Dout as the result of the comparison of intZout with 0
		if(hasLargerZout)
			vhdl << tab << declare("intDout") << " <= intZout(" << wInz <<");" <<endl;
		else
			vhdl << tab << declare("intDout") << " <= intZout(" << wInz-1 <<");" <<endl;
		
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
		mpfr_init2(x, 1+wIx+wFx);
		mpfr_init2(y, 1+wIy+wFy);
		mpfr_init2(z, 1+wIz+wFz);
		mpfr_init2(xint, 1+wIx+wFx);
		mpfr_init2(yint, 1+wIy+wFy);
		mpfr_init2(zint, 1+wIz+wFz);
		mpfr_init2(d, 1+wIx+wFx);
		mpfr_init2(rx, 1+wIx+wFx);
		mpfr_init2(ry, 1+wIy+wFy);
		mpfr_init2(rz, 1+wIz+wFz);
		mpfr_init2(rd, 1+wIx+wFx);
		
		mpz_t rx_z, ry_z, rz_z, rd_z;
		mpz_init2 (rx_z, 1+wIx+wFx);
		mpz_init2 (ry_z, 1+wIy+wFy);
		mpz_init2 (rz_z, 1+wIz+wFz);
		mpz_init2 (rd_z, 1+wIx+wFx);
	
		/* Compute correct value */
		mpfr_set_z (x, svX.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (x, x, wFx, GMP_RNDN); // this rounding is acually exact
		
		mpfr_set_z (y, svY.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (y, y, wFy, GMP_RNDN); // this rounding is acually exact
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, wFz, GMP_RNDN); // this rounding is acually exact
		
		if(svD==0)
			mpfr_set_d (d, 1, GMP_RNDN); 
		else
			mpfr_set_d (d, -1, GMP_RNDN);
		mpfr_div_2si (d, d, wFx, GMP_RNDN); // this rounding is acually exact
		
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
		mpfr_mul_2si (rx, rx, wFx, GMP_RNDN); // exact rnd here
		mpfr_get_z (rx_z, rx, GMP_RNDN); // there can be a real rounding here 
		mpz_class rx_zc(rx_z);
		tc->addExpectedOutput ("Xout", rx_zc);
		
		mpfr_mul_2si (ry, ry, wFy, GMP_RNDN); // exact rnd here
		mpfr_get_z (ry_z, ry, GMP_RNDN); // there can be a real rounding here 
		mpz_class ry_zc(ry_z);
		tc->addExpectedOutput ("Yout", ry_zc);
		
		mpfr_mul_2si (rz, rz, wFz, GMP_RNDN); // exact rnd here
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
		
		mpfr_init2(z, 1+wIz+wFz);
		mpz_init2 (z_z, 1+wIz+wFz);
		
		
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
		mpfr_mul_2si (z, z, wFz, GMP_RNDN); 
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












