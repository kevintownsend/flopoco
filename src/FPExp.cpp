/*
  A floating-point exponential for FloPoCo
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	// For NaN

#include "FPExp.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"

#include "fpexponential/stdfragment.h"
#include "fpexponential/explore.h"

using namespace std;

namespace flopoco{

	FPExp::FPExp(Target* target, int wE, int wF)
		: wE(wE), wF(wF)
	{
		/* Generate unique name */
		{
			std::ostringstream o;
			o << "FPExp_" << wE << "_" << wF;
			uniqueName_ = o.str();
		}

		addFPInput("X", wE, wF);
		addFPOutput("R", wE, wF, 2);  // 2 because faithfully rounded

		int explore_size = wF;

		f = explore(explore_size);
		if (!f) throw std::string("FPExp::FPExp(): No fragment");

		result_length = f->prepare(area, max_error);

		g = intlog2(max_error) + 2;
		cout
			<< "    Estimated area (for fixed-point part): " << area << endl
			<< "    Maximum error: " << max_error << endl
			<< "    Internal precision: " << result_length << endl
			<< "    Precision: " << result_length - g << endl;
	}	

	FPExp::~FPExp()
	{
		delete f;
	}

	// Overloading the virtual functions of Operator
	void FPExp::outputVHDL(std::ostream& o, std::string name)
	{
		licence(o, "J. Detrey, F. de Dinechin, C. Klein, X. Pujol  (2008)");
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
			"           nX  : out std_logic_vector(wE+wF+g-1 downto 0);\n"
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
			"         nX  : out std_logic_vector(wE+wF+g-1 downto 0);\n"
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
			"    nX <= buf(n*(wF+2**n+1)+wF+wE+wX-1 downto n*(wF+2**n+1)+wX-g);\n"
			"  end generate;\n"
			"\n"
			"  padding : if wX < g generate\n"
			"    nX <= buf(n*(wF+2**n+1)+wF+wE+wX-1 downto n*(wF+2**n+1)) & (g-wX-1 downto 0 => '0');\n"
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
			"  signal nX : std_logic_vector(wE+wF+g-1 downto 0); -- fixed point, 2's complement X\n"
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
			"               nX  => nX,\n"
			"               ofl => ofl0,\n"
			"               ufl => ufl0 );\n"
			"\n"
			"  -- compute E\n"
			"  nK0 <= nX(wE+wF+g-2 downto wF+g-4) * cstInvLog2;\n"
			"  nK1 <= (\"0\" & nK0) - (\"0\" & cstInvLog2 & (wE+4-2 downto 0 => '0')) when nX(wE+wF+g-1) = '1' else\n"
			"         \"0\" & nK0;\n"
			"  -- round E\n"
			"  nK <= nK1(wE+4+wE+1 downto 4+wE+1) + ((wE downto 1 => '0') & nK1(4+wE));\n"
			"\n"
			"  -- compute Y\n"
			"  nKLog20 <= nK(wE-1 downto 0) * cstLog2;\n"
			"  nKLog2  <= (\"0\" & nKLog20) - (\"0\" & cstLog2 & (wE-1 downto 0 => '0')) when nK(wE) = '1' else\n"
			"             \"0\" & nKLog20;\n"
			"\n"
			"  nY <= nX - nKLog2(wE+wE-1+wF+g-1 downto wE-1);\n"
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
 
}
