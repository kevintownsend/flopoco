#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	// For NaN

#include "FPExp.hpp"
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

	g=3; // TODO
	p=9; // TODO
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

#if 1 // remove to use KCMs
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
	vhdl << tab  << declare("shiftValIn", shiftInSize) << " <= shiftval" << range(shiftInSize-1, 0) << ";" << endl;

	outPortMap(lshift, "R", "fixX0");
	inPortMap(lshift, "S", "shiftValIn");
	inPortMap(lshift, "X", "mXu");
	vhdl << instance(lshift, "mantissa_shift");
	syncCycleFromSignal("fixX0");

	vhdl << tab  << "-- Partial overflow/underflow detection" << endl;
	vhdl << tab  << declare("oufl0") << " <= not shiftVal(wE+1) when shiftVal(wE downto 0) >= conv_std_logic_vector(" << maxshift << ", wE+1) else '0';" << endl;

 	int sizeXfix = wE-1 + 1+wF+g +1; // +1 for the sign
	vhdl << tab << declare("fixX", sizeXfix) << " <= " << "'0' & fixX0" << range(wE-1 + wF+g + wF+1 -1, wF) << ";" << endl;
	vhdl << tab << "-- two's compliment version mantissa" << endl;
	vhdl << tab  << declare("fixXsigned", sizeXfix) << " <= ((" << sizeXfix-1 << " downto 0 => '0') - fixX) when XSign = '1' else  fixX;" << endl;

 	// fixXsigned(), on wE+3 bits,  is an int . 2^-2
	// cstInvLog2,   on wE+3 bits,  is an int . 2^{-wE-2}
	// The product is an int . 2^{-wE-4}

	vhdl << tab <<	declare("xMulIn",wE+3) << " <=  fixXsigned "<< range(sizeXfix-1, sizeXfix-wE-3) << "; -- 2's complement, rounded down" << endl;
	IntIntKCM *mulInvLog2 = new IntIntKCM(target, wE+3, mpzInvLog2, true /* signed input */);
	oplist.push_back(mulInvLog2);
	outPortMap(mulInvLog2, "R", "xInvLog2ToZero");
	inPortMap(mulInvLog2, "X", "xMulIn");
	vhdl << instance(mulInvLog2, "mulInvLog2");

	vhdl << tab <<	declare("xInvLog2",2*(wE+3)) << " <=   xInvLog2ToZero when XSign='0' else xInvLog2ToZero + ("<< rangeAssign(2*(wE+3)-1, wE+3, "'1'") << " & xMulIn); -- if X<0 xMulIn should be multiplied by RU(1/log2) = 1+invlog2" << endl;

	vhdl << tab << declare("K", wE+1) << " <= xInvLog2 " << range(2*(wE+3)-2, wE+4) << "; -- rouding down in all cases" <<endl;


	IntIntKCM *mulLog2 = new IntIntKCM(target, wE+1, mpzLog2, true /* signed input */);
	oplist.push_back(mulLog2);
	outPortMap(mulLog2, "R", "KLog2");
	inPortMap(mulLog2, "X", "K");
	vhdl << instance(mulLog2, "mulLog2");

	vhdl << tab << declare("Ypadded", sizeXfix) << " <= fixXsigned - KLog2" << range(wE+1 + wE+wF+g -1, wE+1 + wE+wF+g - sizeXfix) << ";\n";
	vhdl << tab << declare("Y", wF+g) << " <= Ypadded" << range(wF+g-1, 0) << ";\n";
	vhdl << tab << declare("ShouldBeZero", wE) << " <= Ypadded" << range(sizeXfix-1, sizeXfix-wE) << "; -- for debug TODO remove\n";
	if(wF+g<=32) {
	vhdl << tab << "-- Now compute the exp of this fixed-point value" <<endl;
	REPORT(LIST, "Generating the polynomial approximation, this may take some time");
	// nY is in [0, 1]
	FunctionEvaluator *fe = new FunctionEvaluator(target, "exp(x), 0,1,1", wF+g, wF+g, d);
	oplist.push_back(fe);
	inPortMap(fe, "X", "Y");
	outPortMap(fe, "R", "expY");
	vhdl << instance(fe, "poly");
	syncCycleFromSignal("expY");
	// here we have in expY the exponential  	
	// int rWeight=fe->getRWeight(); // TODO should be 2, 3 if signed, but I get 4
	//	int rWidth=fe->getRWidth(); // TODO bug, it is the one before rounding
	int rWidth=getSignalByName("expY")->width(); 

 	vhdl << tab << declare("ShouldBeZero2", 2) << " <= expY" << range(rWidth-1, rWidth-2)  << "; -- for debug TODO remove\n";
	vhdl << tab << declare("needNorm") << " <= expY(" << rWidth-3 << ");" << endl;
	vhdl << tab << declare("roundBit") << " <= expY(" << rWidth-4-wF << ")  when needNorm = '1'    else expY(" <<  rWidth-5-wF << ") ;" << endl;


	vhdl << tab << "-- Rounding: all this should consume one row of LUTs" << endl; 
	vhdl << tab << declare("preRoundBiasSig", wE+wF+2)
		  << " <= conv_std_logic_vector(" << bias+1 << ", wE+2)  & expY" << range(rWidth-4, rWidth-4-wF+1) << " when needNorm = '1'" << endl
		  << tab << tab << "else conv_std_logic_vector(" << bias << ", wE+2)  & expY" << range(rWidth-5, rWidth-5-wF+1) << " ;" << endl;
	vhdl << tab << declare("roundNormAddend", wE+wF+2) << " <= K(" << wE << ") & K & "<< rangeAssign(wF-1, 1, "'0'") << " & roundBit;" << endl;
	vhdl << tab << declare("roundedExpSig", wE+wF+2) << " <= preRoundBiasSig + roundNormAddend when Xexn=\"01\" else "
		  << " \"000\" & (wE-2 downto 0 => '1') & (wF-1 downto 0 => '0');" << endl;
	
	vhdl << tab << declare("ofl1") << " <= not XSign and oufl0 and (not Xexn(1) and Xexn(0)); -- input positive, normal,  very large" << endl;
	vhdl << tab << declare("ofl2") << " <= not XSign and (roundedExpSig(wE+wF) and not roundedExpSig(wE+wF+1)) and (not Xexn(1) and Xexn(0)); -- input positive, normal, overflowed" << endl;
 	vhdl << tab << declare("ofl3") << " <= not Xsign and Xexn(1) and not Xexn(0);  -- input was -infty" << endl;
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


#if 0
	vhdl << tab << declare("preRoundSignificand", wF) << " <= nZ(wF+g-1 downto g) when sign = '0'    else  nZ(wF+g-2 downto g-1);" << endl;
	vhdl << tab << declare("preRoundExponent",wE+1) << " <= nK + (\"00\" & (wE-2 downto 1 => '1') & (not sign));" << endl;
	vhdl << tab << declare("preRoundExponentSignificand", wE+wF+1) << " <= preRoundExponent & preRoundSignificand;" << endl;
	vhdl << tab << declare("postRoundExponentSignificand", wE+wF+1) 
		  << " <= preRoundExponentSignificand + nZ(g-1) when sign = '0'    else  preRoundExponentSignificand + nZ(g-2);" << endl;
	vhdl << tab << declare("exponent",wE+1) << " <= postRoundExponentSignificand(wE+wF downto wF);" << endl;
	vhdl << tab << declare("significand", wF) << " <= postRoundExponentSignificand(wF-1 downto 0);" << endl;

	vhdl << tab << declare("ofl1") << " <= ofl0 or exponent(wE);" << endl;
	vhdl << tab << declare("ufl1") << " <= '1' when X(wE+wF+2 downto wE+wF+1) = \"00\"      else  ufl0;" << endl;
	vhdl << tab << declare("ofl2") << " <= '1' when X(wE+wF+2 downto wE+wF+1) = \"10\"      else ofl1 and (not ufl1);" << endl;
	vhdl << tab << "R(wE+wF+2 downto wE+wF+1) <= \"11\"                   when X(wE+wF+2 downto wE+wF+1) = \"11\" else" << endl
		  << "                               (not X(wE+wF)) & \"0\" when ofl2 = '1'                         else" << endl
		  << "                               \"01\";" << endl;
	vhdl << tab << "R(wE+wF downto 0) <= \"00\" & (wE-2 downto 0 => '1') & (wF-1 downto 0 => '0') when ufl1 = '1' else" << endl
		  << "                         \"0\" & exponent(wE-1 downto 0) & significand;" << endl;

#endif
#if 0

	fp_exp <<
		"\n"
		"  signal fixX : std_logic_vector(wE+wF+g-1 downto 0); -- fixed point, 2's complement X\n"
		"  signal nK0 : std_logic_vector(wE+4+wE downto 0);\n"
		"  signal nK1 : std_logic_vector(wE+4+wE+1 downto 0);\n"
		"  signal nK  : std_logic_vector(wE downto 0);\n"
		"  \n"
		"  signal nKLog20 : std_logic_vector(wE+wE-1+wF+g-1 downto 0);\n"
		"  signal nKLog2  : std_logic_vector(wE+wE-1+wF+g downto 0);\n"
		"  signal nY  : std_logic_vector(wE+wF+g-1 downto 0);\n"
		"  signal sign : std_logic;\n"
		"  signal unsigned_input : std_logic_vector(wF+g-2 downto 0);\n"
		"  \n"
		"  signal nZ : std_logic_vector(wF+g-1 downto 0);\n"
		"  \n"
		"  signal preRoundSignificand : std_logic_vector(wF-1 downto 0);\n"
		"  signal preRoundExponent : std_logic_vector(wE downto 0); \n"
		"  signal preRoundExponentSignificand : std_logic_vector(wE+wF downto 0);\n"
		"  signal postRoundExponentSignificand : std_logic_vector(wE+wF downto 0);\n"
		"  signal significand : std_logic_vector(wF-1 downto 0);\n"
		"  signal exponent : std_logic_vector(wE downto 0);\n"
		"\n"
		"  signal sticky : std_logic;\n"
		"  signal round  : std_logic;\n"
		"\n"
		"  signal fR0 : std_logic_vector(wF+1 downto 0);\n"
		"  signal fR1 : std_logic_vector(wF downto 0);\n"
		"  signal fR  : std_logic_vector(wF-1 downto 0);\n"
		"\n"
		"  signal eR : std_logic_vector(wE downto 0);\n"
		"  \n"
		"  signal ofl0 : std_logic;\n"
		"  signal ofl1 : std_logic;\n"
		"  signal ofl2 : std_logic;\n"
		"  signal ufl0 : std_logic;\n"
		"  signal ufl1 : std_logic;\n"
		"  \n"
		"begin\n"
		"  -- conversion of X to fixed point\n"
		"  shift : " << uniqueName_ << "_shift\n"
		"    generic map ( wE => wE,\n"
		"                  wF => wF,\n"
		"                  g => g )\n"
		"    port map ( fpX => X(wE+wF downto 0),\n"
		"               fixX  => fixX,\n"
		"               ofl => ofl0,\n"
		"               ufl => ufl0 );\n"
		"\n"
		"  -- compute E\n"
		"  nK0 <= fixX(wE+wF+g-2 downto wF+g-4) * cstInvLog2;\n"
		"  nK1 <= (\"0\" & nK0) - (\"0\" & cstInvLog2 & (wE+4-2 downto 0 => '0')) when fixX(wE+wF+g-1) = '1' else\n"
		"         \"0\" & nK0;\n"
		"  -- round E\n"
		"  nK <= nK1(wE+4+wE+1 downto 4+wE+1) + ((wE downto 1 => '0') & nK1(4+wE));\n"
		"\n"
		"  -- compute Y\n"
		"  nKLog20 <= nK(wE-1 downto 0) * cstLog2;\n"
		"  nKLog2  <= (\"0\" & nKLog20) - (\"0\" & cstLog2 & (wE-1 downto 0 => '0')) when nK(wE) = '1' else\n"
		"             \"0\" & nKLog20;\n"
		"\n"
		"  nY <= fixX - nKLog2(wE+wE-1+wF+g-1 downto wE-1);\n"
		"  sign <= nY(wF+g-1);\n"
		"  unsigned_input <= nY(wF+g-2 downto 0) when sign = '0' else (wF+g-2 downto 0 => '0') - nY(wF+g-2 downto 0);\n"
		"  \n";

	fp_exp << 
		"  label2 : " << uniqueName_ << "_exp_" << result_length - 1 << endl <<
		"    port map (x    => unsigned_input,\n"
		"              y    => nZ,\n"
		"              sign => sign);\n"
		"\n"
		"  preRoundSignificand(wF-1 downto 0) <= nZ(wF+g-1 downto g) when sign = '0' else\n"
		"                 nZ(wF+g-2 downto g-1);\n"
		"  preRoundExponent <= nK + (\"00\" & (wE-2 downto 1 => '1') & (not sign));\n"
		"  preRoundExponentSignificand <= preRoundExponent & preRoundSignificand;\n"
		"  postRoundExponentSignificand <= preRoundExponentSignificand + nZ(g-1) when sign = '0' else\n"
		"                                  preRoundExponentSignificand + nZ(g-2);\n"
		"  exponent <= postRoundExponentSignificand(wE+wF downto wF);\n"
		"  significand <= postRoundExponentSignificand(wF-1 downto 0);\n"
		"\n"
		"  ofl1 <= ofl0 or exponent(wE);\n"
		"\n"
		"  ufl1 <= '1' when X(wE+wF+2 downto wE+wF+1) = \"00\" else\n"
		"          ufl0;\n"
		"\n"
		"  ofl2 <= '1' when X(wE+wF+2 downto wE+wF+1) = \"10\" else\n"
		"          ofl1 and (not ufl1);\n"
		"  \n"
		"  R(wE+wF+2 downto wE+wF+1) <= \"11\"                   when X(wE+wF+2 downto wE+wF+1) = \"11\" else\n"
		"                                 (not X(wE+wF)) & \"0\" when ofl2 = '1'                         else\n"
		"                                 \"01\";\n"
		"\n"
		"  R(wE+wF downto 0) <= \"00\" & (wE-2 downto 0 => '1') & (wF-1 downto 0 => '0') when ufl1 = '1' else\n"
		"                         \"0\" & exponent(wE-1 downto 0) & significand;\n"
		"end architecture;\n";

#endif

}	

FPExp::~FPExp()
{
}

#if 0
// Overloading the virtual functions of Operator
void FPExp::outputVHDL(std::ostream& o, std::string name)
{
	stringstream fp_exp, fixp_exp, fixp_exp_tbl;

	f->generate(uniqueName_, fixp_exp, fixp_exp_tbl);

	std::string cstInvLog2, cstLog2;
	{
		mpz_class z;

		mpfr_t mp2, mp1, mp;
		mpfr_init2(mp, 2*(wE+wF+g));	// XXX: way too much precision
		mpfr_inits(mp1, mp2, 0, NULL);
		mpfr_set_si(mp1, 1, GMP_RNDN);
		mpfr_set_si(mp2, 2, GMP_RNDN);

		mpfr_log(mp, mp2, GMP_RNDN);
		mpfr_mul_2si(mp, mp, wE-1+wF+g, GMP_RNDN);
		mpfr_get_z(z.get_mpz_t(), mp, GMP_RNDN);
		{
			std::ostringstream o;
			o.fill('0');
			o.width(wE-1+wF+g);
			o << z.get_str(2);
			cstLog2 = o.str();
		}

		mpfr_mul_2si(mp, mp, -(wE-1+wF+g), GMP_RNDN);
		mpfr_div(mp, mp1, mp, GMP_RNDN);
		mpfr_mul_2si(mp, mp, wE+1, GMP_RNDN);
		mpfr_get_z(z.get_mpz_t(), mp, GMP_RNDN);
		{
			std::ostringstream o;
			o.fill('0');
			o.width(wE+1);
			o << z.get_str(2);
			cstInvLog2 = o.str();
		}

		mpfr_clears(mp1, mp2, mp, 0, NULL);
	}

	fp_exp <<
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"use ieee.std_logic_arith.all;\n"
		"use ieee.std_logic_unsigned.all;\n"
		"\n"
		"package pkg_" << uniqueName_ << " is\n"
		"  function min ( x, y : integer ) return integer;\n"
		"  function max ( x, y : integer ) return integer;\n"
		"  function fp_exp_shift_wx ( wE, wF, g : positive ) return positive;\n"
		"  function log2 ( x : positive ) return integer;\n"
		"  \n"
		"  component " << uniqueName_ << "_shift is\n"
		"    generic ( wE : positive;\n"
		"              wF : positive;\n"
		"              g : positive );\n"
		"    port ( fpX : in  std_logic_vector(wE+wF downto 0);\n"
		"           fixX  : out std_logic_vector(wE+wF+g-1 downto 0);\n"
		"           ofl : out std_logic;\n"
		"           ufl : out std_logic );\n"
		"  end component;\n"
		"  \n";

	fp_exp <<
		"  component " << uniqueName_ << " is" << endl <<
		"    generic ( wE : positive := " << wE << ";" << endl <<
		"              wF : positive := " << result_length - g << ";" << endl <<
		"              g : positive := " << g << " );" << endl;

	fp_exp <<
		"    port ( x : in  std_logic_vector(2+wE+wF downto 0);\n"
		"           r : out std_logic_vector(2+wE+wF downto 0) );\n"
		"  end component;\n"
		"end package;\n"
		"\n"
		"package body pkg_" << uniqueName_ << " is\n"
		"  function min ( x, y : integer ) return integer is\n"
		"  begin\n"
		"    if x <= y then\n"
		"      return x;\n"
		"    else\n"
		"      return y;\n"
		"    end if;\n"
		"  end function;\n"
		"\n"
		"  function max ( x, y : integer ) return integer is\n"
		"  begin\n"
		"    if x >= y then\n"
		"      return x;\n"
		"    else\n"
		"      return y;\n"
		"    end if;\n"
		"  end function;\n"
		"  \n"
		"  function fp_exp_shift_wx ( wE, wF, g : positive ) return positive is\n"
		"  begin\n"
		"    return min(wF+g, 2**(wE-1)-1);\n"
		"  end function;\n"
		"  \n"
		"  function log2 ( x : positive ) return integer is\n"
		"    variable n : natural := 0;\n"
		"  begin\n"
		"    while 2**(n+1) <= x loop\n"
		"      n := n+1;\n"
		"    end loop;\n"
		"    return n;\n"
		"  end function;\n"
		"\n"
		"end package body;\n"
		"\n"
		"-- Conversion de l'entrée en virgule fixe\n"
		"-- ======================================\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"use ieee.std_logic_arith.all;\n"
		"use ieee.std_logic_unsigned.all;\n"
		"library work;\n"
		"use work.pkg_" << uniqueName_ << ".all;\n"
		"\n"
		"entity " << uniqueName_ << "_shift is\n"
		"  generic ( wE : positive;\n"
		"            wF : positive;\n"
		"            g : positive );\n"
		"  port ( fpX : in  std_logic_vector(wE+wF downto 0);\n"
		"         fixX  : out std_logic_vector(wE+wF+g-1 downto 0);\n"
		"         ofl : out std_logic;\n"
		"         ufl : out std_logic );\n"
		"end entity;\n"
		"\n"
		"architecture arch of " << uniqueName_ << "_shift is\n"
		"  -- length of the fraction part of X\n"
		"  constant wX : integer  := fp_exp_shift_wx(wE, wF, g);\n"
		"  -- nombre d'étapes du décalage \n"
		"  -- (de l'ordre de log2(taille du nombre en virgule fixe))\n"
		"  constant n  : positive := log2(wX+wE-2)+1;\n"
		"\n"
		"  signal e0 : std_logic_vector(wE+1 downto 0);\n"
		"  signal eX : std_logic_vector(wE+1 downto 0);\n"
		"\n"
		"  signal mXu : std_logic_vector(wF downto 0);\n"
		"  signal mXs : std_logic_vector(wF+1 downto 0);\n"
		"\n"
		"  signal buf : std_logic_vector((n+1)*(wF+2**n+1)-1 downto 0);\n"
		"begin\n"
		"  -- évalue le décalage à effectuer\n"
		"  e0 <= conv_std_logic_vector(2**(wE-1)-1 - wX, wE+2);\n"
		"  eX <= (\"00\" & fpX(wE+wF-1 downto wF)) - e0;\n"
		"\n"
		"  -- underflow quand l'entrée est arrondie à zéro (donc au final exp(entrée) = 1) (?)\n"
		"  ufl <= eX(wE+1);\n"
		"  -- overflow (détection partielle en se basant uniquement sur l'exposant de l'entrée)\n"
		"  ofl <= not eX(wE+1) when eX(wE downto 0) > conv_std_logic_vector(wX+wE-2, wE+1) else\n"
		"         '0';\n"
		"\n"
		"  -- mantisse de l'entrée (rajoute le 1 implicite)\n"
		"  mXu <= \"1\" & fpX(wF-1 downto 0);\n"
		"  -- représentation signée de la mantisse\n"
		"  mXs <= (wF+1 downto 0 => '0') - (\"0\" & mXu) when fpX(wE+wF) = '1' else (\"0\" & mXu);\n"
		"\n"
		"  -- ajoute eX zéros à droite de la mantisse\n"
		"  buf(wF+1 downto 0) <= mXs;\n"
		"  shift : for i in 0 to n-1 generate\n"
		"    buf( (i+1)*(wF+2**n+1) + wF+2**(i+1) downto\n"
		"         (i+1)*(wF+2**n+1) )\n"
		"      <=   -- pas de décalage si eX(i) = 0\n"
		"           ( 2**i-1 downto 0 => buf(i*(wF+2**n+1) + wF+2**i) ) &\n"
		"           buf( i*(wF+2**n+1) + wF+2**i downto\n"
		"                i*(wF+2**n+1) )\n"
		"        when eX(i) = '0' else\n"
		"           -- décalage de 2 ^ i bits si eX(i) = 1\n"
		"           buf( i*(wF+2**n+1) + wF+2**i downto\n"
		"                i*(wF+2**n+1) ) &\n"
		"           ( 2**i-1 downto 0 => '0' );\n"
		"  end generate;\n"
		"\n"
		"  no_padding : if wX >= g generate\n"
		"    fixX <= buf(n*(wF+2**n+1)+wF+wE+wX-1 downto n*(wF+2**n+1)+wX-g);\n"
		"  end generate;\n"
		"\n"
		"  padding : if wX < g generate\n"
		"    fixX <= buf(n*(wF+2**n+1)+wF+wE+wX-1 downto n*(wF+2**n+1)) & (g-wX-1 downto 0 => '0');\n"
		"  end generate;\n"
		"\n"
		"end architecture;\n"
		"\n"
		"-- Exponentielle en virgule flottante\n"
		"-- ==================================\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"use ieee.std_logic_arith.all;\n"
		"use ieee.std_logic_unsigned.all;\n"
		"library work;\n"
		"use work.pkg_" << uniqueName_ << "_exp.all;\n"
		"use work.pkg_" << uniqueName_ << ".all;\n"
		"\n";

	fp_exp <<
		"entity " << uniqueName_ << " is" << endl <<
		"    generic ( wE : positive := " << wE << ";" << endl <<
		"              wF : positive := " << result_length - g << ";" << endl <<
		"              g : positive := " << g << " );" << endl <<
		"    port ( x : in  std_logic_vector(2+wE+wF downto 0);" << endl <<
		"           r : out std_logic_vector(2+wE+wF downto 0));" << endl <<
		"end entity;" << endl << endl <<
		"architecture arch of " << uniqueName_ << " is" << endl <<
		"  constant cstInvLog2 : std_logic_vector(wE+1 downto 0) := \"" << cstInvLog2 << "\";" << endl <<
		"  constant cstLog2 : std_logic_vector(wE-1+wF+g-1 downto 0) := \"" << cstLog2 << "\";" << endl;

	fp_exp <<
		"\n"
		"  signal fixX : std_logic_vector(wE+wF+g-1 downto 0); -- fixed point, 2's complement X\n"
		"  signal nK0 : std_logic_vector(wE+4+wE downto 0);\n"
		"  signal nK1 : std_logic_vector(wE+4+wE+1 downto 0);\n"
		"  signal nK  : std_logic_vector(wE downto 0);\n"
		"  \n"
		"  signal nKLog20 : std_logic_vector(wE+wE-1+wF+g-1 downto 0);\n"
		"  signal nKLog2  : std_logic_vector(wE+wE-1+wF+g downto 0);\n"
		"  signal nY  : std_logic_vector(wE+wF+g-1 downto 0);\n"
		"  signal sign : std_logic;\n"
		"  signal unsigned_input : std_logic_vector(wF+g-2 downto 0);\n"
		"  \n"
		"  signal nZ : std_logic_vector(wF+g-1 downto 0);\n"
		"  \n"
		"  signal preRoundSignificand : std_logic_vector(wF-1 downto 0);\n"
		"  signal preRoundExponent : std_logic_vector(wE downto 0); \n"
		"  signal preRoundExponentSignificand : std_logic_vector(wE+wF downto 0);\n"
		"  signal postRoundExponentSignificand : std_logic_vector(wE+wF downto 0);\n"
		"  signal significand : std_logic_vector(wF-1 downto 0);\n"
		"  signal exponent : std_logic_vector(wE downto 0);\n"
		"\n"
		"  signal sticky : std_logic;\n"
		"  signal round  : std_logic;\n"
		"\n"
		"  signal fR0 : std_logic_vector(wF+1 downto 0);\n"
		"  signal fR1 : std_logic_vector(wF downto 0);\n"
		"  signal fR  : std_logic_vector(wF-1 downto 0);\n"
		"\n"
		"  signal eR : std_logic_vector(wE downto 0);\n"
		"  \n"
		"  signal ofl0 : std_logic;\n"
		"  signal ofl1 : std_logic;\n"
		"  signal ofl2 : std_logic;\n"
		"  signal ufl0 : std_logic;\n"
		"  signal ufl1 : std_logic;\n"
		"  \n"
		"begin\n"
		"  -- conversion of X to fixed point\n"
		"  shift : " << uniqueName_ << "_shift\n"
		"    generic map ( wE => wE,\n"
		"                  wF => wF,\n"
		"                  g => g )\n"
		"    port map ( fpX => X(wE+wF downto 0),\n"
		"               fixX  => fixX,\n"
		"               ofl => ofl0,\n"
		"               ufl => ufl0 );\n"
		"\n"
		"  -- compute E\n"
		"  nK0 <= fixX(wE+wF+g-2 downto wF+g-4) * cstInvLog2;\n"
		"  nK1 <= (\"0\" & nK0) - (\"0\" & cstInvLog2 & (wE+4-2 downto 0 => '0')) when fixX(wE+wF+g-1) = '1' else\n"
		"         \"0\" & nK0;\n"
		"  -- round E\n"
		"  nK <= nK1(wE+4+wE+1 downto 4+wE+1) + ((wE downto 1 => '0') & nK1(4+wE));\n"
		"\n"
		"  -- compute Y\n"
		"  nKLog20 <= nK(wE-1 downto 0) * cstLog2;\n"
		"  nKLog2  <= (\"0\" & nKLog20) - (\"0\" & cstLog2 & (wE-1 downto 0 => '0')) when nK(wE) = '1' else\n"
		"             \"0\" & nKLog20;\n"
		"\n"
		"  nY <= fixX - nKLog2(wE+wE-1+wF+g-1 downto wE-1);\n"
		"  sign <= nY(wF+g-1);\n"
		"  unsigned_input <= nY(wF+g-2 downto 0) when sign = '0' else (wF+g-2 downto 0 => '0') - nY(wF+g-2 downto 0);\n"
		"  \n";

	fp_exp << 
		"  label2 : " << uniqueName_ << "_exp_" << result_length - 1 << endl <<
		"    port map (x    => unsigned_input,\n"
		"              y    => nZ,\n"
		"              sign => sign);\n"
		"\n"
		"  preRoundSignificand(wF-1 downto 0) <= nZ(wF+g-1 downto g) when sign = '0' else\n"
		"                 nZ(wF+g-2 downto g-1);\n"
		"  preRoundExponent <= nK + (\"00\" & (wE-2 downto 1 => '1') & (not sign));\n"
		"  preRoundExponentSignificand <= preRoundExponent & preRoundSignificand;\n"
		"  postRoundExponentSignificand <= preRoundExponentSignificand + nZ(g-1) when sign = '0' else\n"
		"                                  preRoundExponentSignificand + nZ(g-2);\n"
		"  exponent <= postRoundExponentSignificand(wE+wF downto wF);\n"
		"  significand <= postRoundExponentSignificand(wF-1 downto 0);\n"
		"\n"
		"  ofl1 <= ofl0 or exponent(wE);\n"
		"\n"
		"  ufl1 <= '1' when X(wE+wF+2 downto wE+wF+1) = \"00\" else\n"
		"          ufl0;\n"
		"\n"
		"  ofl2 <= '1' when X(wE+wF+2 downto wE+wF+1) = \"10\" else\n"
		"          ofl1 and (not ufl1);\n"
		"  \n"
		"  R(wE+wF+2 downto wE+wF+1) <= \"11\"                   when X(wE+wF+2 downto wE+wF+1) = \"11\" else\n"
		"                                 (not X(wE+wF)) & \"0\" when ofl2 = '1'                         else\n"
		"                                 \"01\";\n"
		"\n"
		"  R(wE+wF downto 0) <= \"00\" & (wE-2 downto 0 => '1') & (wF-1 downto 0 => '0') when ufl1 = '1' else\n"
		"                         \"0\" & exponent(wE-1 downto 0) & significand;\n"
		"end architecture;\n";

	o << fixp_exp_tbl.str() << fixp_exp.str() << fp_exp.str();
}
#endif


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
		tc->addFPInput("X", FPNumber::plusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::minusDirtyZero);
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

	
		mpfr_clears(x, y, NULL);
	}



}
