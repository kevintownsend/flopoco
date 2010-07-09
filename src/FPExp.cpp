#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	// For NaN

#include "FPExp.hpp"
#include "FPNumber.hpp"
#include "ConstMult/IntIntKCM.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "Shifters.hpp"
#include "IntMultiplier.hpp"
#include "FunctionEvaluator.hpp"
#include "utils.hpp"
#include "IntAdder.hpp"


using namespace std;

// TODO fix-point KCM by a real constant

namespace flopoco{
	extern vector<Operator *> oplist;

  FPExp::FPExp(Target* target, int wE, int wF) 
	: Operator(target), wE(wE), wF(wF)
{
    cout << "TODO!" <<endl;
    int k=11;
    int d=2;
    FPExp(target, wE, wF, k, d);
  }

  FPExp::FPExp(Target* target, int wE, int wF, int k, int d)
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

	// The two constants
	mpz_class mpzLog2, mpzInvLog2;
	
	mpfr_t mp2, mp1, mplog2, mpinvlog2;

	mpfr_inits(mp1, mp2, NULL);
	mpfr_set_si(mp1, 1, GMP_RNDN);
	mpfr_set_si(mp2, 2, GMP_RNDN);









	// 1/log2 ~ 1.44, truncated, on sizeFirstKCM bits
	int sizeFirstKCM=wE+4;
	mpfr_init2(mplog2, 3*(wE+wF+g));	// way too much precision
	mpfr_log(mplog2, mp2, GMP_RNDN);
	mpfr_init2(mpinvlog2, sizeFirstKCM);	
	mpfr_div(mpinvlog2, mp1, mplog2, GMP_RNDN);
	mpfr_mul_2si(mpinvlog2, mpinvlog2, sizeFirstKCM-1, GMP_RNDN); //Exact
	mpfr_get_z(mpzInvLog2.get_mpz_t(), mpinvlog2, GMP_RNDN);
	mpfr_clear(mplog2);

	//Computing log2 ~ 0.69 on wE+wF+g bits, rounded down, too
	mpfr_init2(mplog2, wE+wF+g);
	mpfr_log(mplog2, mp2, GMP_RNDN);
	mpfr_mul_2si(mplog2, mplog2, wE+wF+g, GMP_RNDN); //Exact
	mpfr_get_z(mpzLog2.get_mpz_t(), mplog2, GMP_RNDN);

	mpfr_clears(mp1, mp2, mplog2, mpinvlog2, NULL);

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
	vhdl << tab << declare("addOp0",sizeXfix) << " <= (fixX xor "<<rangeAssign(sizeXfix-1,0,"XSign")<<") and "<<rangeAssign(sizeXfix-1,0,"not(expIs1)")<<";"<<endl;
		
	vhdl << tab << declare("addCarry") << "<= XSign and not(expIs1);"<<endl;
//	nextCycle();
	
	IntAdder *twosComplMantissa = new IntAdder(target, sizeXfix);
	oplist.push_back(twosComplMantissa);
	
	inPortMap(twosComplMantissa, "X", "addOp0");
	inPortMapCst(twosComplMantissa, "Y", zg(sizeXfix,0));
	inPortMap(twosComplMantissa, "Cin", "addCarry");
	outPortMap(twosComplMantissa, "R", "fixXsigned");
	
	vhdl << instance( twosComplMantissa, "The2cMantissa") << endl;
	syncCycleFromSignal("fixXsigned");
	nextCycle(); ////////????FIXME
	
//	vhdl << tab  << declare("fixXsigned", sizeXfix) << " <= " << endl
//		  << tab << "(" << sizeXfix-1 << " downto 0 => '0') when expIs1='1'" << endl
//		  << tab << "else ((" << sizeXfix-1 << " downto 0 => '0') - fixX) when XSign = '1'" << endl
//		  << tab << "else fixX;" << endl;

 	// fixXsigned(), on wE+3 bits,  is an int . 2^-2
	// cstInvLog2,   on wE+3 bits,  is an int . 2^{-wE-2}
	// The product is an int . 2^{-wE-4}

	vhdl << tab <<	declare("xMulIn",sizeFirstKCM) << " <=  fixXsigned "<< range(sizeXfix-1, sizeXfix-sizeFirstKCM) << "; -- 2's complement, rounded down" << endl;

	if((wE+wF)*target->normalizedFrequency() > 14)
		nextCycle();


	// TODO all the following should be implemented as a single specific KCM
	IntIntKCM *mulInvLog2 = new IntIntKCM(target, sizeFirstKCM, mpzInvLog2, true /* signed input */);
	oplist.push_back(mulInvLog2);
	outPortMap(mulInvLog2, "R", "xInvLog2");
	inPortMap(mulInvLog2, "X", "xMulIn");
	vhdl << instance(mulInvLog2, "mulInvLog2");
	syncCycleFromSignal("xInvLog2");

	if((wE+wF)*target->normalizedFrequency() > 8)
		nextCycle(); 
	int xInvlog2size=getSignalByName("xInvLog2")->width();
	vhdl << tab << declare("Kunrounded", wE+2) << " <= xInvLog2 " << range(xInvlog2size-2, xInvlog2size -wE-3) << " + ( " << rangeAssign(wE+1,1,"'0'") << "  & '1'); -- adding round bit " <<endl;
	vhdl << tab << declare("K", wE+1) << " <= Kunrounded " << range(wE+1,1) << "; -- rounding to nearest" <<endl;

	// end TODO


	if((wE+wF)*target->normalizedFrequency() > 4)
		nextCycle();

#if 0

	IntIntKCM *mulLog2 = new IntIntKCM(target, wE+1, mpzLog2, true /* signed input */);
	oplist.push_back(mulLog2);
	outPortMap(mulLog2, "R", "KLog2");
	inPortMap(mulLog2, "X", "K");
	vhdl << instance(mulLog2, "mulLog2");
	syncCycleFromSignal("KLog2");
	nextCycle();

	vhdl << tab << declare("neg2Op",sizeXfix) << " <= not(KLog2" << range(wE+1 + wE+wF+g -1, wE+1 + wE+wF+g - sizeXfix)<<");"<<endl;
#else
	FixRealKCM *mulLog2 = new FixRealKCM(target, 0, wE, true  /* signed input */, -wF-g, "log(2)" );
	oplist.push_back(mulLog2);
	outPortMap(mulLog2, "R", "KLog2");
	inPortMap(mulLog2, "X", "K");
	vhdl << instance(mulLog2, "mulLog2");
	syncCycleFromSignal("KLog2");
	nextCycle();

	vhdl << tab << declare("neg2Op",wF+g) << " <= not(KLog2" << range(wF+g-1, 0) << ");"<<endl;
	vhdl << tab << declare("fixXsignedLSB",wF+g) << " <= fixXsigned" << range(wF+g-1, 0) << ";"<<endl;

#endif
	
	IntAdder *yPaddedAdder = new IntAdder(target,wF+g); // we know the leading bits will cancel out
	oplist.push_back(yPaddedAdder);
	
	inPortMap( yPaddedAdder, "X", "fixXsignedLSB");
	inPortMapCst ( yPaddedAdder, "Y", "neg2Op");
	inPortMapCst ( yPaddedAdder, "Cin", "'1'");
	outPortMap( yPaddedAdder, "R", "Y");
	
	vhdl << instance(yPaddedAdder, "theYAdder") << endl;
	syncCycleFromSignal("Y"); 
	nextCycle(); 

//	vhdl << tab << declare("Ypadded", sizeXfix) << " <= fixXsigned - KLog2" << range(wE+1 + wE+wF+g -1, wE+1 + wE+wF+g - sizeXfix) << ";\n";

	int sizeY=wF+g; // This is also the weight of Y's LSB

	//	vhdl << tab << declare("Y", sizeY) << " <= Ypadded" << range(sizeXfix-wE-2, 0) << ";\n";

	vhdl << tab << "-- Now compute the exp of this fixed-point value" <<endl;

	int rWidth;

	// nY is in [-1/2, 1/2]
	int sizeZ = wF+g-k; 
	int sizeExpA = wF+g+1; // e^A has MSB weight 1
	int sizeZhigh=wF+g-2*k;
	int sizeZxpZmZm1 = wF+g - 2*k +1;
	int sizeExpZm1 = sizeZ+1; // 
	int sizeMultIn = sizeZ; // sacrificing accuracy where it costs
	REPORT(2, "sizeZ=" << sizeZ);
	REPORT(2, "sizeExpA=" << sizeExpA);
	REPORT(2, "sizeZhigh=" << sizeZhigh);
	REPORT(2, "sizeZxpZmZm1=" << sizeZxpZmZm1);
	REPORT(2, "sizeExpZm1=" << sizeExpZm1);
	//REPORT(2, "=" << );

	vhdl << tab << declare("Addr1", k) << " <= Y" << range(sizeY-1, sizeY-k) << ";\n";
	vhdl << tab << declare("Z", sizeZ) << " <= Y" << range(sizeZ-1, 0) << ";\n";

 
	if(0 && wF==23) {
		//--------------------------------Special case for single precision-----------------
		// TODO and smaller
		// For g=3 we enter with Y on 26 bits. Better first pad it
		vhdl << tab << declare("Y0", 27) << " <= Y & '0';\n";
		sizeY=27;  // (23+3+1) TODO this is a difference to std case
		k=9;
		int sizeZ=18; // should be wF+g-k; 
		// CHECK k+sizeZ == sizeY


		vhdl << tab << declare("Addr2", k) << " <= Z" << range(sizeZ-1, sizeZ-k) << ";\n";
		magicExpTable* table;
		table = new magicExpTable(target);
		oplist.push_back(table);
		outPortMap(table, "Y2", "lowerTerm0");
		inPortMap(table, "X2", "Addr2");
		outPortMap(table, "Y1", "expA0");
		inPortMap(table, "X1", "Addr1");
		vhdl << instance(table, "table");

		vhdl << tab << declare("expA", 28) << " <=  expA0" << range(35, 8) << ";" << endl;

		syncCycleFromSignal("expA");
		if((wE+wF)*target->normalizedFrequency() > 8)
			nextCycle();
		if((wE+wF)*target->normalizedFrequency() > 16)
			nextCycle(); // To get high-speed BRam speed
		
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


		vhdl << tab << declare("Zhigh", sizeZhigh) << " <= Z" << range(sizeZ-1, sizeZ-sizeZhigh) << ";\n";

		if(wF<=23) { // Magic exp table works up to single precision

#if 0 //// stupid ISE not able to pack both tables in a single dual-port one
			LowerExpTable* lowertable;
			lowertable = new LowerExpTable(target, k, sizeZhigh, wF+g); // last parameter is -LSB of the result
			oplist.push_back(lowertable);
			outPortMap(lowertable, "Y", "expZmZm1_0");
			inPortMap(lowertable, "X", "Zhigh");
			vhdl << instance(lowertable, "expZmZm1_table");
#else
		vhdl << tab << declare("Addr2", k) << " <= Z" << range(sizeZ-1, sizeZ-k) << ";\n";
		magicExpTable* table;
		table = new magicExpTable(target);
		oplist.push_back(table);
		outPortMap(table, "Y2", "lowerTerm0");
		inPortMap(table, "X2", "Addr2");
		outPortMap(table, "Y1", "expA0");
		inPortMap(table, "X1", "Addr1");
		vhdl << instance(table, "table");

		vhdl << tab << declare("expA", 27) << " <=  expA0" << range(35, 9) << ";" << endl;
		vhdl << tab << declare("expZmZm1_0", k) << " <= lowerTerm0" << range(k-1, 0) << ";" << endl;


#endif
		}
		else { // use a polynomial evaluator
			firstExpTable* table;
			table = new firstExpTable(target, k, sizeExpA); // e^A-1 has MSB weight 1
			oplist.push_back(table);
			outPortMap(table, "Y", "expA");
			inPortMap(table, "X", "Addr1");
			vhdl << instance(table, "table");

			REPORT(LIST, "Generating the polynomial approximation, this may take some time");
			// We want the LSB value to be  2^(wF+g)
			FunctionEvaluator *fe;
			ostringstream function;
			function << "exp(x*1b-" << k << ")-x*1b-" << k << "-1, 0,1,1";
			fe = new FunctionEvaluator(target, function.str(), sizeZhigh, wF+g, d);
			oplist.push_back(fe);
			inPortMap(fe, "X", "Zhigh");
			outPortMap(fe, "R", "expZmZm1_0");
			vhdl << instance(fe, "poly");
		}

		syncCycleFromSignal("expA");

		syncCycleFromSignal("expZmZm1_0");

		if((wE+wF)*target->normalizedFrequency() > 8)
			nextCycle();
		if((wE+wF)*target->normalizedFrequency() > 16)
			nextCycle(); // To get high-speed BRam speed


		

		// here we have in expZmZm1 e^Z-Z-1

		// Alignment of expZmZm10:  MSB has weight -2*k, LSB has weight -(wF+g).
		// FunctionEvaluator, in its ignorance, may have added MSB bits 
		//		vhdl << tab << declare("ShouldBeZero2", (rWidth- sizeZxpZmZm1)) << " <= expZmZm1_0" << range(rWidth-1, sizeZxpZmZm1)  << "; -- for debug to check it is always 0" <<endl;
		vhdl << tab << declare("expZmZm1", sizeZxpZmZm1) << " <= expZmZm1_0" << range(sizeZxpZmZm1-1, 0)  << "; " <<endl;
		
		vhdl << tab << "-- Computing Z + (exp(Z)-1-Z)" << endl;
		vhdl << tab << declare("expZminus1", sizeExpZm1) << " <= ('0' & Z) + (" << rangeAssign(sizeZ, sizeZ-k+1, "'0'") << " & expZmZm1) ;" << endl;

		vhdl << tab << "-- Truncating expA to the same accuracy as expZminus1" << endl;
		vhdl << tab << declare("expArounded0", sizeMultIn+1) << " <= (expA" << range(sizeExpA-1, sizeExpA-sizeMultIn-1) << ") + '1';  -- rounding" << endl;
		vhdl << tab << declare("expArounded", sizeMultIn) << " <= expArounded0" << range(sizeMultIn, 1) << ";" << endl;
		//		vhdl << tab << declare("expZm1truncated", sizeMultIn) << " <= expZminus1" << range(sizeExpZm1-1, sizeExpZm1-sizeMultIn) << ";" << endl;

//		if((wE+wF)*target->normalizedFrequency() > 4)
			nextCycle(); //in any case this will be absorbed by the DSP blocks

		IntMultiplier *lowProd = new IntMultiplier(target, sizeMultIn, sizeExpZm1);
		oplist.push_back(lowProd);
		
		inPortMap(lowProd, "X", "expArounded");
		inPortMap(lowProd, "Y", "expZminus1");
		outPortMap(lowProd, "R", "lowerProduct");
		
		vhdl << instance(lowProd, "TheLowerProduct")<<endl;
		 
//		vhdl << tab << declare("lowerProduct", 2*sizeExpZm1) << " <= expAtruncated * expZminus1;" << endl;

		syncCycleFromSignal("lowerProduct");

//		if((wE+wF)*target->normalizedFrequency() > 14) //TODO will need to buffer out before addition
//			nextCycle();

		vhdl << tab << "-- Final addition -- the product MSB bit weight is -k+2 = "<< -k+2 << endl;
		
		map<string,double> finalAdderInDelayMap;
		finalAdderInDelayMap["Y"] = (lowProd->getOutDelayMap())["R"];
#if 0 // Trying to round the product instead of truncating it		Why doesn't that work?
		IntAdder *finalAdder = new IntAdder(target, sizeExpA+1, finalAdderInDelayMap, 2, 1, -1);
		oplist.push_back(finalAdder);
		
		vhdl << tab << declare("extendedLowerProduct",sizeExpA+1) << " <= (" << rangeAssign(sizeExpA-1, sizeExpA-k+1, "'0'") 
		     << " & lowerProduct" << range(sizeMultIn+sizeExpZm1-1, sizeMultIn+sizeExpZm1 - (sizeExpA-k+1) -1) << ");" << endl;		
		vhdl << tab << declare("extendedExpA",sizeExpA+1) << " <= expA & '1'; -- rounding bit for the product" << endl;		
		     
		     
		inPortMap(finalAdder, "X", "extendedExpA");
		inPortMap(finalAdder, "Y", "extendedLowerProduct");
		inPortMapCst(finalAdder, "Cin", "'0'");
		outPortMap(finalAdder, "R", "expY0");
		
		vhdl << instance(finalAdder,"TheFinalAdder") << endl;
		syncCycleFromSignal("expY0");
		vhdl << tab << declare("expY",sizeExpA) << " <= expY0" << range(sizeExpA, 1) << " ; -- rounding" << endl;		
		     		
#else  //  truncating it		
		IntAdder *finalAdder = new IntAdder(target, sizeExpA, finalAdderInDelayMap, 2, 1, -1);
		oplist.push_back(finalAdder);
		
		vhdl << tab << declare("extendedLowerProduct",sizeExpA) << " <= (" << rangeAssign(sizeExpA-1, sizeExpA-k+1, "'0'") 
		     << " & lowerProduct" << range(sizeMultIn+sizeExpZm1-1, sizeMultIn+sizeExpZm1 - (sizeExpA-k+1)) << ");" << endl;		
		     		     
		inPortMap(finalAdder, "X", "expA");
		inPortMap(finalAdder, "Y", "extendedLowerProduct");
		inPortMapCst(finalAdder, "Cin", "'0'");
		outPortMap(finalAdder, "R", "expY");
		
		vhdl << instance(finalAdder,"TheFinalAdder") << endl;
		syncCycleFromSignal("expY");
		     		
#endif

//		nextCycle(); //thsi one should be controlled by some parameter
		
//		cout << ">>>>>>>>>>>>>>>>>>>>...The delay on the adder output R is" << finalAdder->getOutDelayMap()["R"] << endl;
			
//		vhdl << tab << declare("expY", sizeExpA) << " <= expA + (" << rangeAssign(sizeExpA-1, sizeExpA-k+1, "'0'") 
//		     << " & lowerProduct" << range(2*sizeExpZm1-1, 2*sizeExpZm1- (sizeExpA-k+1)) << ");" << endl;

		rWidth=sizeExpA; // for the rounding
			       
				       
#else
		throw string("FPExp requires Sollya for this precision, sorry.");
#endif
	}


	// The following is generic normalization/rounding code if we have an approx of exp(y) of size rwidth in expY
	// with MSB of weight 2^1
	// We start a cycle here
	nextCycle();

	vhdl << tab << declare("needNoNorm") << " <= expY(" << rWidth-1 << ");" << endl;
	vhdl << tab << declare("roundBit") << " <= expY(" << rWidth-2-wF << ")  when needNoNorm = '1'    else expY(" <<  rWidth-3-wF << ") ;" << endl;

	vhdl << tab << "-- Rounding: all this should consume one row of LUTs" << endl; 
	vhdl << tab << declare("preRoundBiasSig", wE+wF+2)
		  << " <= conv_std_logic_vector(" << bias << ", wE+2)  & expY" << range(rWidth-2, rWidth-2-wF+1) << " when needNoNorm = '1'" << endl
		  << tab << tab << "else conv_std_logic_vector(" << bias-1 << ", wE+2)  & expY" << range(rWidth-3, rWidth-3-wF+1) << " ;" << endl;
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

	// One test out of 8 fully random (tests NaNs etc)
	// All the remaining ones test numbers with exponents between -wF-3 and wE-2,
        // For numbers outside this range, exp over/underflows or flushes to 1. 
 
	TestCase* FPExp::buildRandomTestCase(int i){
		TestCase *tc;
		tc = new TestCase(this); 
		mpz_class x;
		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);
		mpz_class bias = ((1<<(wE-1))-1);
		/* Fill inputs */
		if ((i & 7) == 0) { //fully random
		  x = getLargeRandom(wE+wF+3);
		}
		else
		  {
		    mpz_class e = (getLargeRandom(wE+wF) % (wE+wF+2) ) -wF-3; // Should be between -wF-3 and wE-2
		    //cout << e << endl;
		    e = bias + e;
		    mpz_class sign = getLargeRandom(1);
		    x  = getLargeRandom(wF) + (e << wF) + (sign<<(wE+wF)) + normalExn;
		  }
		tc->addInput("X", x);
		/* Get correct outputs */
		emulate(tc);
		return tc;
	}


}
