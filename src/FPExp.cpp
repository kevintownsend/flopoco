#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	// For NaN

#include "FPExp.hpp"
//#include "fpexponential/SPExpDualTable.hpp"
#include "FPNumber.hpp"
#include "ConstMult/IntIntKCM.hpp"
#include "Shifters.hpp"
#include "FunctionEvaluator.hpp"
#include "utils.hpp"

//#include "fpexponential/Fragment.hpp"
//#include "fpexponential/explore.h"

using namespace std;


namespace flopoco{
	extern vector<Operator *> oplist;


FPExp::FPExp(Target* target, int wE, int wF)
	: Operator(target), wE(wE), wF(wF)
{
	/* Generate unique name */
	{
		std::ostringstream o;
		o << "FPExp_" << wE << "_" << wF;
		uniqueName_ = o.str();
	}

	setCopyrightString("F. de Dinechin (2008)");
	srcFileName="FPExp";
	

	addFPInput("X", wE, wF);
	addFPOutput("R", wE, wF, 2);  // 2 because faithfully rounded

	// TODO a bit more science here. Not sure the following code works for any other g
	g=3; 


	int d=2; // degree: TODO

	// The two constants
	mpz_class mpzLog2, mpzInvLog2;
	
	mpfr_t mp2, mp1, mplog2, mpinvlog2;

	mpfr_inits(mp1, mp2, NULL);
	mpfr_set_si(mp1, 1, GMP_RNDN);
	mpfr_set_si(mp2, 2, GMP_RNDN);

	// 1/log2 ~ 1.44, truncated, on wE+3 bits
	mpfr_init2(mplog2, 3*(wE+wF+g));	// way too much precision
	mpfr_log(mplog2, mp2, GMP_RNDN);
	mpfr_init2(mpinvlog2, wE+3);	
	mpfr_div(mpinvlog2, mp1, mplog2, GMP_RNDD);
	mpfr_mul_2si(mpinvlog2, mpinvlog2, wE+2, GMP_RNDN); //Exact
	mpfr_get_z(mpzInvLog2.get_mpz_t(), mpinvlog2, GMP_RNDN);
	mpfr_clear(mplog2);

	//Computing log2 ~ 0.69 on wE+wF+g bits, rounded down, too
	mpfr_init2(mplog2, wE+wF+g);
	mpfr_log(mplog2, mp2, GMP_RNDD);
	mpfr_mul_2si(mplog2, mplog2, wE+wF+g, GMP_RNDN); //Exact
	mpfr_get_z(mpzLog2.get_mpz_t(), mplog2, GMP_RNDN);

#if 0 // remove to use KCMs
	std::string  cstLog2;
	{
		std::ostringstream o;
		o.fill('0');
		o.width(wE+wF+g);
		o << mpzLog2.get_str(2);
		cstLog2 = o.str();
	}

	
	std::string cstInvLog2;
	{
 		std::ostringstream o;
 		o.fill('0');
 		o.width(wE+3);
 		o << mpzInvLog2.get_str(2);
 		cstInvLog2 = o.str();
 	}
	
	vhdl << tab << declare("cstInvLog2", wE+3) << " <= \"" << cstInvLog2 << "\";" << endl;
	vhdl << tab  << declare("cstLog2", wE+wF+g) << " <= \"" << cstLog2 << "\";" << endl;
#endif

	mpfr_clears(mp1, mp2, mplog2, mpinvlog2, NULL);

	// int sizeY=wF+g-1;
	// int sizeY2=wF+g-1 -p;
	// int sizepolyIn=wF+g-2*p-3; // As in the FPT2005 paper
	// int sizepolyOut=wF+g-2*p-2; // As in the FPT2005 paper

	addConstant("wE", "positive", wE);
	addConstant("wF", "positive", wF);
	addConstant("g", "positive", g);

	int bias = (1<<(wE-1))-1;
	if(bias < wF+g){
		ostringstream e;
		e << "ERROR in FPExp, unable to build architecture if wF+g > 2^(wE-1)-1." <<endl;
		e << "      Try increasing wE." << endl;
		e << "      If you really need FPExp to work with such values, please report this as a bug :)" << endl;
		throw e.str();
		}


	vhdl << tab  << declare("Xexn", 2) << " <= X(wE+wF+2 downto wE+wF+1);" << endl;
	vhdl << tab  << declare("XSign") << " <= X(wE+wF);" << endl;
	vhdl << tab  << declare("XexpField", wE) << " <= X(wE+wF-1 downto wF);" << endl;
 	vhdl << tab  << declare("Xfrac", wF) << " <= X(wF-1 downto 0);" << endl;
	int e0 = bias - (wF+g);
	vhdl << tab  << declare("e0", wE+2) << " <= conv_std_logic_vector(" << e0 << ", wE+2);  -- bias - (wF+g)" << endl;
	vhdl << tab  << declare("shiftVal", wE+2) << " <= (\"00\" & XexpField) - e0; -- for a left shift" << endl;

	vhdl << tab  << "-- underflow when input is shifted to zero (shiftval<0), in which case exp = 1" << endl;
	vhdl << tab  << declare("expIs1") << " <= shiftVal(wE+1);" << endl;
 
	// As we don't have a signed shifter, shift first, complement next. TODO? replace with a signed shifter
	vhdl << tab << "--  mantissa with implicit bit" << endl;
	vhdl << tab  << declare("mXu", wF+1) << " <= \"1\" & Xfrac;" << endl;

 	// left shift
	int maxshift=wE-1+ wF+g; // maxX < 2^(wE-1); 
	Shifter* lshift = new Shifter(target, wF+1, maxshift , Shifter::Left);   
	oplist.push_back(lshift);
	int shiftInSize = lshift->getShiftInWidth();
	vhdl << tab  << declare("shiftValIn", shiftInSize) << " <= shiftVal" << range(shiftInSize-1, 0) << ";" << endl;

	outPortMap(lshift, "R", "fixX0");
	inPortMap(lshift, "S", "shiftValIn");
	inPortMap(lshift, "X", "mXu");
	vhdl << instance(lshift, "mantissa_shift");
	syncCycleFromSignal("fixX0");
	nextCycle();

	vhdl << tab  << "-- Partial overflow/underflow detection" << endl;
	vhdl << tab  << declare("oufl0") << " <= not shiftVal(wE+1) when shiftVal(wE downto 0) >= conv_std_logic_vector(" << maxshift << ", wE+1) else '0';" << endl;

 	int sizeXfix = wE-1 + 1+wF+g +1; // +1 for the sign == wE+wF+g +1
	vhdl << tab << declare("fixX", sizeXfix) << " <= " << "'0' & fixX0" << range(wE-1 + wF+g + wF+1 -1, wF) << ";" << endl;
	vhdl << tab << "-- two's compliment version mantissa" << endl;
	vhdl << tab  << declare("fixXsigned", sizeXfix) << " <= " << endl
		  << tab << "(" << sizeXfix-1 << " downto 0 => '0') when expIs1='1'" << endl
		  << tab << "else ((" << sizeXfix-1 << " downto 0 => '0') - fixX) when XSign = '1'" << endl
		  << tab << "else fixX;" << endl;

 	// fixXsigned(), on wE+3 bits,  is an int . 2^-2
	// cstInvLog2,   on wE+3 bits,  is an int . 2^{-wE-2}
	// The product is an int . 2^{-wE-4}

	vhdl << tab <<	declare("xMulIn",wE+3) << " <=  fixXsigned "<< range(sizeXfix-1, sizeXfix-wE-3) << "; -- 2's complement, rounded down" << endl;

	if((wE+wF)*target->normalizedFrequency() > 14)
		nextCycle();

	IntIntKCM *mulInvLog2 = new IntIntKCM(target, wE+3, mpzInvLog2, true /* signed input */);
	oplist.push_back(mulInvLog2);
	outPortMap(mulInvLog2, "R", "xInvLog2ToZero");
	inPortMap(mulInvLog2, "X", "xMulIn");
	vhdl << instance(mulInvLog2, "mulInvLog2");
	syncCycleFromSignal("xInvLog2ToZero");

	if((wE+wF)*target->normalizedFrequency() > 8)
		nextCycle(); 
	
	vhdl << tab <<	declare("xInvLog2",2*(wE+3)) << " <=   xInvLog2ToZero when XSign='0' else xInvLog2ToZero + ("<< rangeAssign(2*(wE+3)-1, wE+3, "'1'") << " & xMulIn); -- if X<0 xMulIn should be multiplied by RU(1/log2) = 1+invlog2" << endl;

	vhdl << tab << declare("K", wE+1) << " <= xInvLog2 " << range(2*(wE+3)-2, wE+4) << "; -- rounding down in all cases" <<endl;

	if((wE+wF)*target->normalizedFrequency() > 4)
		nextCycle();

	IntIntKCM *mulLog2 = new IntIntKCM(target, wE+1, mpzLog2, true /* signed input */);
	oplist.push_back(mulLog2);
	outPortMap(mulLog2, "R", "KLog2");
	inPortMap(mulLog2, "X", "K");
	vhdl << instance(mulLog2, "mulLog2");
	syncCycleFromSignal("KLog2");
	nextCycle();

	vhdl << tab << declare("Ypadded", sizeXfix) << " <= fixXsigned - KLog2" << range(wE+1 + wE+wF+g -1, wE+1 + wE+wF+g - sizeXfix) << ";\n";

	int sizeY=wF+g; // This is also the weight of Y's LSB

	vhdl << tab << declare("Y", sizeY) << " <= Ypadded" << range(sizeXfix-wE-2, 0) << ";\n";
	vhdl << tab << declare("ShouldBeZero", wE+1) << " <= Ypadded" << range(sizeXfix-1, sizeXfix-wE-1) << "; -- for debug TODO remove\n";

	vhdl << tab << "-- Now compute the exp of this fixed-point value" <<endl;

	int rWidth;


 
	if(wE==8 && wF==23) {
		//--------------------------------Special case for single precision-----------------
		// TODO and smaller
		// For g=3 we enter with Y on 26 bits. Better first pad it
		vhdl << tab << declare("Y0", 27) << " <= Y & '0';\n";

		sizeY=27;  // (23+3+1) TODO this is a difference to std case
		k=9;
		int sizeZ=18; // should be wF+g-k; 
		// CHECK k+sizeZ == sizeY


		vhdl << tab << declare("Addr1", k) << " <= Y0" << range(sizeY-1, sizeY-k) << ";\n";
		vhdl << tab << declare("Z", sizeZ) << " <= Y0" << range(sizeZ-1, 0) << ";\n";
		vhdl << tab << declare("Addr2", k) << " <= Z" << range(sizeZ-1, sizeZ-k) << ";\n";
		magicExpTable* table;
		table = new magicExpTable(target);
		oplist.push_back(table);
		outPortMap(table, "Y2", "lowerTerm0");
		inPortMap(table, "X2", "Addr2");
		outPortMap(table, "Y1", "expA0");
		inPortMap(table, "X1", "Addr1");
		vhdl << instance(table, "table");
		syncCycleFromSignal("expA0");
		if((wE+wF)*target->normalizedFrequency() > 8)
			nextCycle();
		if((wE+wF)*target->normalizedFrequency() > 16)
			nextCycle(); // To get high-speed BRam speed
		
		// TODO if needed: since the leading bit is always 1, maybe we can win one bit of accuracy on the input of the product
		vhdl << tab << "-- Adding back the leading bit" << endl;
		
		vhdl << tab << declare("expAh", 2) << " <= '0' & expA0(35);" << endl;
		vhdl << tab << declare("expA", 28) << " <= (\"01\" + expAh) & expA0" << range(34, 34-25) << ";" << endl;
		
		vhdl << tab << "-- Computing Z + (exp(z)-1-Z)" << endl;
		vhdl << tab << declare("lowerTerm", k) << " <= lowerTerm0" << range(k-1, 0) << ";" << endl;
		vhdl << tab << declare("expZminus1", sizeZ+1) << " <= ('0' & Z) + (" << rangeAssign(sizeZ, sizeZ-k, "'0'") << " & lowerTerm) ;" << endl;
		// TODO This is nonsense. No need to store the extra accuracy then.
		vhdl << tab << "-- Magical rounding-by-truncation to 17 bits, thanks to an half-ulp stored in the table" << endl;
		vhdl << tab << declare("expZminus1rounded", 17) << " <= expZminus1" << range(sizeZ, 2) << ";" << endl;

		vhdl << tab << declare("expAtruncated", 17) << " <= expA" << range(27, 27-16) << ";" << endl;

		if((wE+wF)*target->normalizedFrequency() > 4)
			nextCycle();

		vhdl << tab << declare("lowerProduct", 34) << " <= expAtruncated * expZminus1rounded;" << endl;
		if((wE+wF)*target->normalizedFrequency() > 14)
			nextCycle();

		vhdl << tab << "-- Final addition -- the product ranges on bit weights -7 to -40" << endl;
		vhdl << tab << declare("expY", 28) << " <= expA + (" << rangeAssign(27, 27-7, "'0'") << " & lowerProduct" << range(33, 33-19) << ");" << endl;

		rWidth=28; // for the rounding
	}


	// Generic polynomial evaluation, up to double-precision
	else{
#ifdef HAVE_SOLLYA
		// nY is in [0, 1]
		// values for double-precision. For the other ones: TODO
		k=11;
		d=2; // TODO for DP 2 is better
		int sizeZ=wF+g-k; 

		vhdl << tab << declare("Addr1", k) << " <= Y" << range(sizeY-1, sizeY-k) << ";\n";
		vhdl << tab << declare("Z", sizeZ) << " <= Y" << range(sizeZ-1, 0) << ";\n";
		firstExpTable* table;
		table = new firstExpTable(target, k, wF+g+1); // e^A-1 has MSB weight 1
		oplist.push_back(table);
		outPortMap(table, "Y", "expA0");
		inPortMap(table, "X", "Addr1");
		vhdl << instance(table, "table");
		syncCycleFromSignal("expA0");
		if((wE+wF)*target->normalizedFrequency() > 8)
			nextCycle();
		if((wE+wF)*target->normalizedFrequency() > 16)
			nextCycle(); // To get high-speed BRam speed
		vhdl << tab << "-- Adding back the leading bit" << endl;
		vhdl << tab << declare("expAh", 2) << " <= '0' & expA0(wF+g);" << endl;
		int expAsize = wF+g+2; // e^A has MSB weight 2
		vhdl << tab << declare("expA", expAsize) << " <= (\"01\" + expAh) & expA0" << range(expAsize-3, 0) << ";" << endl;

		int sizeZhigh=wF+g-2*k;
		vhdl << tab << declare("Zhigh", sizeZhigh) << " <= Z" << range(sizeZ-1, sizeZ-sizeZhigh) << ";\n";
		REPORT(LIST, "Generating the polynomial approximation, this may take some time");
		
		// We want the LSB value to be  2^(wF+g)
		FunctionEvaluator *fe = new FunctionEvaluator(target, "exp(x)-x-1, 0,1,1", sizeZhigh, wF+g, d);
		oplist.push_back(fe);
		inPortMap(fe, "X", "Zhigh");
		outPortMap(fe, "R", "expZmZm1_0");
		vhdl << instance(fe, "poly");
		syncCycleFromSignal("expZmZm1_0");
		// here we have in expZmZm1 the exponential  	x
		// int rWeight=fe->getRWeight(); // TODO should be 2, 3 if signed, but I get 4
		//	int rWidth=fe->getRWidth(); // TODO bug, it is the one before rounding
		rWidth=getSignalByName("expZmZm1_0")->width(); 

		// Alignment of expZmZm10:  MSB has weight -2*k, LSB has weight -(wF+g).
		int expZmZm1size = wF+g - 2*k +1;
		// FunctionEvaluator, in its ignorance, may have added MSB bits 
		vhdl << tab << declare("ShouldBeZero2", (rWidth- expZmZm1size)) << " <= expZmZm1_0" << range(rWidth-1, expZmZm1size)  << "; -- for debug to check it is always 0" <<endl;
		vhdl << tab << declare("expZmZm1", expZmZm1size) << " <= expZmZm1_0" << range(expZmZm1size-1, 0)  << "; -- the clean one" <<endl;
		
		
		vhdl << tab << "-- Computing Z + (exp(Z)-1-Z)" << endl;
		vhdl << tab << declare("expZminus1", sizeZ+1) << " <= ('0' & Z) + (" << rangeAssign(sizeZ, sizeZ-k, "'0'") << " & expZmZm1) ;" << endl;

		vhdl << tab << "-- Truncating expA to the same accuracy as expZminus1" << endl;
		vhdl << tab << declare("expAtruncated", sizeZ+1) << " <= expA" << range(expAsize-1, expAsize-sizeZ-1) << ";" << endl;

		if((wE+wF)*target->normalizedFrequency() > 4)
			nextCycle();

		vhdl << tab << declare("lowerProduct", 2*(sizeZ+1)) << " <= expAtruncated * expZminus1;" << endl;
		if((wE+wF)*target->normalizedFrequency() > 14)
			nextCycle();

		vhdl << tab << "-- Final addition -- the product MSB bit weight is -k+2 = "<< -k+2 << endl;
		//TODO fix here
		vhdl << tab << declare("expY", expAsize) << " <= expA + (" << rangeAssign(expAsize-1, expAsize-k-2, "'0'") 
		     << " & lowerProduct" << range(2*(sizeZ+1)-1, 2*(sizeZ+1)- (expAsize-k-2)) << ");" << endl;

		rWidth=expAsize; // for the rounding
			       
				       
#else
		throw string("FPExp requires Sollya for this precision, sorry.");
#endif
	}


	// The following is generic normalization/rounding code if we have an approx of exp(y) of size rwidth in expY.
	// We start a cycle here
	nextCycle();

	vhdl << tab << declare("needNoNorm") << " <= expY(" << rWidth-1 << ");" << endl;
	vhdl << tab << declare("roundBit") << " <= expY(" << rWidth-2-wF << ")  when needNoNorm = '1'    else expY(" <<  rWidth-3-wF << ") ;" << endl;

	vhdl << tab << "-- Rounding: all this should consume one row of LUTs" << endl; 
	vhdl << tab << declare("preRoundBiasSig", wE+wF+2)
		  << " <= conv_std_logic_vector(" << bias+1 << ", wE+2)  & expY" << range(rWidth-2, rWidth-2-wF+1) << " when needNoNorm = '1'" << endl
		  << tab << tab << "else conv_std_logic_vector(" << bias << ", wE+2)  & expY" << range(rWidth-3, rWidth-3-wF+1) << " ;" << endl;
	vhdl << tab << declare("roundNormAddend", wE+wF+2) << " <= K(" << wE << ") & K & "<< rangeAssign(wF-1, 1, "'0'") << " & roundBit;" << endl;

	if((wE+wF)*target->normalizedFrequency() > 16)
		nextCycle();
	vhdl << tab << declare("roundedExpSig", wE+wF+2) << " <= preRoundBiasSig + roundNormAddend when Xexn=\"01\" else "
		  << " \"000\" & (wE-2 downto 0 => '1') & (wF-1 downto 0 => '0');" << endl;
	
	vhdl << tab << declare("ofl1") << " <= not XSign and oufl0 and (not Xexn(1) and Xexn(0)); -- input positive, normal,  very large" << endl;
	vhdl << tab << declare("ofl2") << " <= not XSign and (roundedExpSig(wE+wF) and not roundedExpSig(wE+wF+1)) and (not Xexn(1) and Xexn(0)); -- input positive, normal, overflowed" << endl;
 	vhdl << tab << declare("ofl3") << " <= not XSign and Xexn(1) and not Xexn(0);  -- input was -infty" << endl;
 	vhdl << tab << declare("ofl") << " <= ofl1 or ofl2 or ofl3;" << endl;

	vhdl << tab << declare("ufl1") << " <= (roundedExpSig(wE+wF) and roundedExpSig(wE+wF+1))  and (not Xexn(1) and Xexn(0)); -- input normal" << endl;
 	vhdl << tab << declare("ufl2") << " <= XSign and Xexn(1) and not Xexn(0);  -- input was -infty" << endl;
	vhdl << tab << declare("ufl3") << " <= XSign and oufl0  and (not Xexn(1) and Xexn(0)); -- input negative, normal,  very large" << endl;

 	vhdl << tab << declare("ufl") << " <= ufl1 or ufl2 or ufl3;" << endl;

	vhdl << tab << declare("Rexn", 2) << " <= \"11\" when Xexn = \"11\"" << endl
 		  << tab << tab << "else \"10\" when ofl='1'" << endl
		  << tab << tab << "else \"00\" when ufl='1'" << endl
		  << tab << tab << "else \"01\";" << endl;
		
	vhdl << tab << "R <= Rexn & '0' & roundedExpSig" << range(wE+wF-1, 0) << ";" << endl;
	


}	

FPExp::~FPExp()
{
}



void FPExp::emulate(TestCase * tc)
{
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("X");

	/* Compute correct value */
	FPNumber fpx(wE, wF);
	fpx = svX;

	mpfr_t x, ru,rd;
	mpfr_init2(x,  1+wF);
	mpfr_init2(ru, 1+wF);
	mpfr_init2(rd, 1+wF); 
	fpx.getMPFR(x);
	mpfr_exp(rd, x, GMP_RNDD);
	mpfr_exp(ru, x, GMP_RNDU);
	FPNumber  fprd(wE, wF, rd);
	FPNumber  fpru(wE, wF, ru);
	mpz_class svRD = fprd.getSignalValue();
	mpz_class svRU = fpru.getSignalValue();
	tc->addExpectedOutput("R", svRD);
	tc->addExpectedOutput("R", svRU);
	mpfr_clears(x, ru, rd, NULL);
}
 


void FPExp::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		mpfr_t x, y;
		FPNumber *fx, *fy;
		// double d;

		mpfr_init2(x, 1+wF);
		mpfr_init2(y, 1+wF);



		tc = new TestCase(this); 
		tc->addFPInput("X", log(2));
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::plusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::minusDirtyZero);
		emulate(tc);
		tcl->add(tc);



		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 2.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 1.5);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", -1.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", -2.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", -3.0);
		emulate(tc);
		tcl->add(tc);




		tc = new TestCase(this); 
		tc->addComment("The largest number whose exp is finite");
		fx = new FPNumber(wE, wF, FPNumber::largestPositive);
		fx->getMPFR(x);
		mpfr_log(y, x, GMP_RNDN);
		//		cout << "A " << fx->getSignalValue() << endl;
		//		 d = mpfr_get_d(x, GMP_RNDN);
		// cout << d << endl;
		// d = mpfr_get_d(y, GMP_RNDN);
		// cout << d << endl;
		fy = new FPNumber(wE, wF, y); 
		tc->addFPInput("X", fy);
		emulate(tc);
		tcl->add(tc);
		delete(fx); 

		tc = new TestCase(this); 
		tc->addComment("The first number whose exp is infinite");
		mpfr_nextabove(y);
		fy = new FPNumber(wE, wF, y); 
		tc->addFPInput("X", fy);
		emulate(tc);
		tcl->add(tc);
		delete(fy);




		tc = new TestCase(this); 
		tc->addComment("The last number whose exp is nonzero");
		fx = new FPNumber(wE, wF, FPNumber::smallestPositive);
		fx->getMPFR(x);
		mpfr_log(y, x, GMP_RNDU);

		// cout << "A " << fx->getSignalValue() << endl;
		// d = mpfr_get_d(x, GMP_RNDN);
		// cout << d << endl;
		// d = mpfr_get_d(y, GMP_RNDN);
		// cout << d << endl;

		fy = new FPNumber(wE, wF, y); 
		tc->addFPInput("X", fy);
		emulate(tc);
		tcl->add(tc);
		delete(fx); 

		tc = new TestCase(this); 
		tc->addComment("The first number whose exp flushes to zero");
		mpfr_nextbelow(y);
		fy = new FPNumber(wE, wF, y); 
		tc->addFPInput("X", fy);
		emulate(tc);
		tcl->add(tc);
		delete(fy);



	
		mpfr_clears(x, y, NULL);
	}



}
