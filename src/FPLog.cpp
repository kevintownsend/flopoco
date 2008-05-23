#include <fstream>
#include <sstream>
#include <math.h>	// for NaN

#include "FPLog.hpp"
#include "FloFP.hpp"
#include "utils.hpp"

#include "fplog/FirstInvTable.hpp"
#include "fplog/FirstLogTable.hpp"
#include "fplog/SecondInvTable.hpp"
#include "fplog/OtherLogTable.hpp"

#include "LZOC.hpp"
#include "Shifters.hpp"

using namespace std;

// XXX: Not well encapsulated
extern vector<Operator *> oplist;

// TODO: get rid of this
void mpfr_shift_left(mpfr_t& x, int s) {
  mpfr_t two;
  int i;

  mpfr_init2(two, 2);
  mpfr_set_ui(two, 2, GMP_RNDD);
  for(i=0; i<s; i++) {
    mpfr_mul(x, x, two, GMP_RNDD);
  }
  mpfr_clear(two);
}


FPLog::FPLog(Target* target, int wE, int wF)
	: Operator(target), wE(wE), wF(wF)
{
	/* Generate unique name */
	{
		std::ostringstream o;
		o << "FPLog_" << wE << "_" << wF;
		unique_name = o.str();
	}

	add_FP_input("x", wE, wF);
	add_FP_output("r", wE, wF);

	int i;

	gLog = 3; // TODO: almost randomly chosen
	target_prec = wF+((wF+1)>>1) +gLog;

	// First compute the precision of each iteration 

	// Stage 0
	p[0] = 0;
	a[0] = 5; 
	s[0] = wF+2;  
	sfullZ[0] = wF+2;

	p[1] = a[0]-1; 
	sbt[1] = wF+2 ;
	s[1] = wF+2;
	t[1] = 0;
	sfullZ[1] = wF+7;

	// Following stages -- this is true starting from stage 1, although
	// stage 1 needs a specific inverter table
	i=1;
	while(2*p[i] < wF+2){ // ensures 2*p[stages+1]>= wF+2, enough for faithful rounding
		a[i] = 4;
		p[i+1] = p[i] + a[i] - 1;

		// size before truncation
		sbt[i+1] = s[i] +  p[i] + 2;

		if(p[i+1]+sbt[i+1] <= target_prec) 
			{ // No truncation at all
				psize[i] = s[i];
				s[i+1] = sbt[i+1];
				t[i+1] = 0;
			}
		else
			{ // Truncate everybody to targetprec : 
				// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>eps[i]
				// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
				// We need only this many bits of Z[i] to compute it.
				psize[i] = target_prec - 2*p[i];
				if (psize[i]>s[i]) // in the first iterations
				  psize[i]=s[i];
				s[i+1] = target_prec - p[i+1];
				t[i+1] = sbt[i+1] - s[i+1];
			}

 
 
		sfullZ[i+1] =  sfullZ[i] + a[i] + p[i] + 1;
		i++;
	}  


	// Deduce the number of stages
	stages = i-1;

	// MSB of squarer input is p[stages+1];
	// LSB will be target_prec
	// size will be target_prec - p[stages+1]  

	cout<<"LogArch::LogArch: needs "<<stages<<" stages"<<endl;

	// Now allocate the various table objects

	it0 = new FirstInvTable(a[0], a[0]+1);
	lt0 = new FirstLogTable(a[0], target_prec, it0);
	it1 = new SecondInvTable(a[1], p[1]);
	for(i=1; i<=stages; i++) {
		lt[i] = new OtherLogTable(a[i], target_prec - p[i], p[i], i); 
	 }
	int computedG = intlog2(3*(stages+1));
}	

FPLog::~FPLog()
{
}

// Overloading the virtual functions of Operator
void FPLog::output_vhdl(std::ostream& o, std::string name)
{
	/* Reuse other FloPoCo stuff */
	int t_pipelined = target->is_pipelined();
	target->set_not_pipelined();

	LZOC *lzoc = new LZOC(target, wF, intlog2(wF));
	oplist.push_back(lzoc);

	if (t_pipelined) target->set_pipelined();

	Licence(o, "Cristian KLEIN (2008)");
	
	int i;
	mpfr_t two;
	mpfr_t log2;
	mpz_t zlog2;

	o <<
		"-------------------------------------------------------------------------------\n"
		"-- Left barrel shifter\n"
		"--\n"
		"-- Generics:\n"
		"--   - w : width of the input operand\n"
		"--   - n : number of shifting stages (width of the signal s)\n"
		"--\n"
		"-- Ports:\n"
		"--   - i [in]  : input signal\n"
		"--   - s [in]  : shift by s\n"
		"--   - o [out] :  output\n"
		"--\n"
		"-- Recursive structure with n stages.\n"
		"-------------------------------------------------------------------------------\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"\n"
		"entity " << unique_name << "_lshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"end entity;\n"
		"\n"
		"architecture arch of " << unique_name << "_lshift is\n"
		"  component " << unique_name << "_lshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"  end component;\n"
		"  signal o0 : std_logic_vector(w-1 downto 0);\n"
		"begin\n"
		"\n"
		"  check1: if 2**(n-1)>=w generate\n"
		"    o0 <=      i                      when s(n-1) = '0'\n"
		"          else (w-1 downto 0 => '0');\n"
		"  end generate;\n"
		"  check2: if 2**(n-1)<w generate\n"
		"    o0 <= i   when s(n-1) = '0' else\n"
		"          i(w-2**(n-1)-1 downto 0) & (2**(n-1)-1 downto 0 => '0');\n"
		"  end generate;\n"
		"  \n"
		"  ------------------------------------------------------------- Recursive stage\n"
		"\n"
		"  recursive : if n > 1 generate\n"
		"    shift0 : " << unique_name << "_lshift\n"
		"      generic map ( w => w,\n"
		"                    n => n-1 )\n"
		"      port map (  i => o0,\n"
		"                  s => s(n-2 downto 0),\n"
		"                  o => o   );\n"
		"  end generate;\n"
		"\n"
		"  ----------------------------------------------------------------- Final stage\n"
		"  single : if n = 1 generate\n"
		"    o <= o0;\n"
		"  end generate;\n"
		"\n"
		"end architecture; -------------------------------------------------------------\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"-------------------------------------------------------------------------------\n"
		"-- Right barrel shifter\n"
		"--\n"
		"-- Generics:\n"
		"--   - w : width of the input operand\n"
		"--   - n : number of shifting stages (width of the signal s)\n"
		"--\n"
		"-- Ports:\n"
		"--   - i [in]  : input signal\n"
		"--   - s [in]  : shift by s\n"
		"--   - o [out] :  output\n"
		"--\n"
		"-- Recursive structure with n stages.\n"
		"-------------------------------------------------------------------------------\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"\n"
		"entity " << unique_name << "_rshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"end entity;\n"
		"\n"
		"  architecture arch of " << unique_name << "_rshift is\n"
		"  component " << unique_name << "_rshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"  end component;\n"
		"  signal o0 : std_logic_vector(w-1 downto 0);\n"
		"begin\n"
		"\n"
		"  check1: if 2**(n-1)>=w generate\n"
		"    o0 <= i   when s(n-1) = '0' else\n"
		"     (w-1 downto 0 => '0');\n"
		"  end generate;\n"
		"  check2: if 2**(n-1)<w generate\n"
		"  o0 <= i   when s(n-1) = '0' else\n"
		"           (w-1 downto w-2**(n-1) => '0')  &  i(w-1 downto 2**(n-1));\n"
		"  end generate;\n"
		"  \n"
		"  ------------------------------------------------------------- Recursive stage\n"
		"\n"
		"  recursive : if n > 1 generate\n"
		"    shift0 : " << unique_name << "_rshift\n"
		"      generic map ( w => w,\n"
		"                    n => n-1 )\n"
		"      port map (  i => o0,\n"
		"                  s => s(n-2 downto 0),\n"
		"                  o => o   );\n"
		"  end generate;\n"
		"\n"
		"  ----------------------------------------------------------------- Final stage\n"
		"  single : if n = 1 generate\n"
		"    o <= o0;\n"
		"  end generate;\n"
		"\n"
		"end architecture; -------------------------------------------------------------\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"-------------------------------------------------------------------------------\n"
		"-- Leading-zero counter and normalization\n"
		"--\n"
		"-- Generics:\n"
		"--   - w : width of the input operand\n"
		"--   - n : number of LZC and shifting stages (width of the signal z)\n"
		"--\n"
		"-- Ports:\n"
		"--   - i [in]  : input signal\n"
		"--   - z [out] : number of leading zeros\n"
		"--   - o [out] : normalized signal\n"
		"--\n"
		"-- Recursive structure with n stages. At most 2^n-1 leading zeros are counted.\n"
		"-------------------------------------------------------------------------------\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"\n"
		"entity " << unique_name << "_lzc_norm is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i : in  std_logic_vector(w-1 downto 0);\n"
		"           z : out std_logic_vector(n-1 downto 0);\n"
		"           o : out std_logic_vector(w-1 downto 0) );\n"
		"end entity;\n"
		"\n"
		"architecture arch of " << unique_name << "_lzc_norm is\n"
		"  component " << unique_name << "_lzc_norm is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i : in  std_logic_vector(w-1 downto 0);\n"
		"           z : out std_logic_vector(n-1 downto 0);\n"
		"           o : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  signal z0 : std_logic;\n"
		"  signal o0 : std_logic_vector(w-1 downto 0);\n"
		"begin ---------------------------------------------- Test 2^(n-1) leading zeros\n"
		"  z0 <= '1' when i(w-1 downto w-2**(n-1)) = (w-1 downto w-2**(n-1) => '0') else\n"
		"        '0';\n"
		"\n"
		"  o0 <= i                                                    when z0 = '0' else\n"
		"        i(w-2**(n-1)-1 downto 0) & (2**(n-1)-1 downto 0 => '0');\n"
		"   z(n-1) <= z0;\n"
		"    ----------------------------------------------------------------- Final stage\n"
		"  single : if n = 1 generate\n"
		"    o <= o0;\n"
		"  end generate;\n"
		"  ------------------------------------------------------------- Recursive stage\n"
		"  recursive : if n > 1 generate\n"
		"    lzc_norm0 : " << unique_name << "_lzc_norm\n"
		"      generic map ( w => w,\n"
		"                    n => n-1 )\n"
		"      port map ( i => o0,\n"
		"                 z => z(n-2 downto 0),\n"
		"                 o => o );\n"
		"  end generate;\n"
		"\n"
		"end architecture; -------------------------------------------------------------\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"use ieee.std_logic_arith.all;\n"
		"use ieee.std_logic_unsigned.all;\n"
		"\n"
		"-- Entity fp_log defined here\n"
		" \n";

	// TODO replace with constants
	o << "entity " << name << " is " << endl;
	o << "  generic ( wE : positive := " << wE  << ";" << endl;
	o << "            wF : positive := "  << wF << ");" << endl;
	o << "  port ( X : in  std_logic_vector(2+wE+wF downto 0);" << endl;
	o << "         R : out std_logic_vector(2+wE+wF downto 0)";
	if (is_sequential())
		o  << ";" << endl << "         clk : in  std_logic  );" << endl;
	else
		o  << "  );" << endl;
		
	o << "end entity;" << endl;
	o << "architecture arch of "  << name << " is " << endl;

	lzoc->output_vhdl_component(o);
	o <<
		"\n"
		"  component " << unique_name << "_lshift is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"           s  : in std_logic_vector(n-1 downto 0);\n"
		"           o  : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"	
		"  component " << unique_name << "_rshift is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"           s  : in std_logic_vector(n-1 downto 0);\n"
		"           o  : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  component " << unique_name << "_lzc_norm is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i : in  std_logic_vector(w-1 downto 0);\n"
		"           z : out std_logic_vector(n-1 downto 0);\n"
		"           o : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  -- def of component range_red, and many constants and signals,  here\n";

	o <<   "  component " << unique_name << "_range_red is port ("<<endl;
	o <<   "            Y0 : in  std_logic_vector("<<wF+1<<" downto 0);"<<endl;
	o <<   "            A  : in std_logic_vector("<< a[0]-1 <<" downto 0);"<<endl;
	if(is_sequential()) 
		o << "          clk  : in std_logic;"<<endl;
	o <<   "            Z  : out std_logic_vector("<< s[stages+1]-1 <<" downto 0);"<<endl;
	o <<   "    almostLog  : out std_logic_vector("<< lt0->wOut-1 <<" downto 0)  );"<<endl;
	o <<   "  end component;"<<endl;
	o <<   "  constant g   : positive := "<<gLog<<";"<<endl;
	o <<   "  constant a0 : positive := "<< a[0] <<";"<<endl;
	o <<   "  constant log2wF : positive := "<< intlog2(wF) <<";"<<endl;
	o <<   "  constant targetprec : positive := "<< target_prec <<";"<<endl;
	o <<   "  constant sfinal : positive := "<< s[stages+1] <<";"<<endl;
	o <<   "  constant pfinal : positive := "<< p[stages+1] <<";"<<endl;
	o <<   "  constant lzc_size : positive := "<< max(intlog2(wF), intlog2(wE+p[stages+1]+1)) << ";" << endl;

	// the maximum shift distance in the "small" path ?
	// o << "  constant shiftvalsize : positive := "<< intlog2(wF + 1 - p[stages+1]) <<";"<<endl;

	// The log2 constant
	mpfr_init2(two, 2);
	mpfr_set_d(two, 2.0, GMP_RNDN);
	mpfr_init2(log2, wF+gLog);
	mpfr_log(log2, two, GMP_RNDN);
	mpfr_shift_left(log2, wF+gLog);
	mpz_init2(zlog2, wF+gLog);
	mpfr_get_z(zlog2, log2, GMP_RNDN);
	o << "  signal log2 : std_logic_vector(wF+g-1 downto 0) := \"";
	printBinPosNumGMP(o, mpz_class(zlog2), wF+gLog);
	o << "\";"<<endl;
	o << "  signal E0offset : std_logic_vector(wE-1 downto 0) := \"";
	printBinPosNumGMP(o, (mpz_class(1)<<(wE-1)) -2 + wE , wE);
	o << "\"; -- E0 + wE "<<endl;
	o << "  signal pfinal_s : std_logic_vector(log2wF -1 downto 0) := \"";
	printBinPosNumGMP(o, mpz_class(p[stages+1]), intlog2(wF));
	o << "\";"<<endl;

	o <<
		"  signal FirstBit : std_logic;\n"
		"  signal Y0 : std_logic_vector(wF+1 downto 0);\n"
		"  signal E  : std_logic_vector(wE-1 downto 0);\n"
		"  signal absE  : std_logic_vector(wE-1 downto 0);\n"
		"  signal absELog2  : std_logic_vector(wF+wE+g-1 downto 0);\n"
		"  signal absELog2_pad, LogF_normal_pad, Log_normal, Log_normal_normd  : std_logic_vector(wE+targetprec-1 downto 0);\n"
		"  signal E_small,ER  : std_logic_vector(wE-1 downto 0);\n"
		"  signal E_normal  : std_logic_vector(lzc_size-1 downto 0);\n"
		"  signal Log_small_normd, Log_g  : std_logic_vector(wF+g-1 downto 0);\n"
		"  signal EFR  : std_logic_vector(wE+wF-1 downto 0);\n"
		"  signal lzo : std_logic_vector(log2wF-1 downto 0);\n"
		"  signal shiftval : std_logic_vector(log2wF downto 0);\n"
		"  signal absZ0, absZ0s : std_logic_vector(wF-pfinal+1 downto 0);\n"
		"  signal Zfinal, Log1p_normal : std_logic_vector(sfinal-1 downto 0);\n"
		"  signal Z2o2_full: std_logic_vector(2*(sfinal-pfinal) -1 downto 0);\n"
		"  signal squarerIn: std_logic_vector(sfinal-pfinal-1 downto 0);\n"
		"  signal Z2o2_small_s, Z2o2: std_logic_vector(sfinal-pfinal downto 0);\n"
		"  signal Log_small, Z_small, Z2o2_small: std_logic_vector(wF+g+1 downto 0);  \n"
		"  signal almostLog, logF_normal : std_logic_vector(targetprec-1 downto 0);\n"
		"  signal E0_sub : std_logic_vector(1 downto 0);\n"
		"  signal sR, small, doRR, ufl, sticky, round: std_logic;\n"
		"begin\n"
		"\n"
		"  FirstBit <=  X(wF-1);\n"
		"  Y0 <=      \"1\"  & X(wF-1 downto 0) & \"0\" when FirstBit = '0'\n"
		"        else \"01\" & X(wF-1 downto 0);\n"
		"\n"
		"  E  <= (X(wE+wF-1 downto wF)) - (\"0\" & (wE-2 downto 1 => '1') & (not FirstBit));\n"
		"\n"
		"  sR <= '0'   when    X(wE+wF-1 downto wF)   =   '0' &(wE-2 downto 0 => '1')  -- binade [1..2)\n"
		"        else not X(wE+wF-1);                -- MSB of exponent\n"
		"\n"
		"  absE <= ((wE-1 downto 0 => '0') - E)   when sR = '1'\n"
		"          else E;\n"
		"\n"
		"  absELog2 <= absE * log2;\n"
		"  \n"
		"  lzoc1 : " << lzoc->unique_name << "\n"
		"    port map (  i => Y0(wF downto 1), ozb => FirstBit,  o => lzo);\n"
		"\n"
		"  shiftval <= ('0' & lzo) - ('0' & pfinal_s); \n"
		"\n"
		"  doRR <= shiftval(log2wF);             -- sign of the result\n"
		"\n"
		"  small <= '1' when ((E=(wE-1 downto 0 => '0')) and (doRR='0'))\n"
		"          else '0';\n"
		"\n"
		"-- The range reduction instance\n"
		"  rr: " << unique_name << "_range_red\n"
		"     port map ( A => X(wF-1 downto wF-a0), Y0 => Y0,\n";

	if(is_sequential())
		o << "                clk=>clk,"<<endl;

	o <<
			"                Z => Zfinal, almostLog => almostLog);\n"
			"  absZ0 <=   Y0(wF-pfinal+1 downto 0)          when (sR='0') else\n"
			"             ((wF-pfinal+1 downto 0 => '0') - Y0(wF-pfinal+1 downto 0));\n"
			"\n"
			"--  absZ0 <=   Y0(wF-pfinal downto 0)   xor (wF-pfinal downto 0 => sR);\n"
			"\n"
			"  lshiftsmall: " << unique_name << "_lshift\n"
			"    generic map (w => wF-pfinal+2, n => log2wF)\n"
			"    port map (  i => absZ0, s => shiftval(log2wF-1 downto 0), o => absZ0s );\n"
			"\n"
			"  -- Z2o2 will be of size sfinal-pfinal, set squarer input size to that\n"
			"  sqintest: if sfinal > wf+2 generate\n"
			"    squarerIn <= Zfinal(sfinal-1 downto pfinal) when doRR='1'\n"
			"                 else (absZ0s &  (sfinal-wF-3 downto 0 => '0'));  \n"
			"  end generate sqintest;\n"
			"  sqintest2: if sfinal <= wf+2 generate\n"
			"    squarerIn <= Zfinal(sfinal-1 downto pfinal) when doRR='1'\n"
			"                 else absZ0s(wF-pfinal+1 downto wf+2-sfinal);  \n"
			"  end generate sqintest2;\n"
			"\n"
			"  -- Z2o2 will be of size sfinal - pfinal -1, set squarer input size to that\n"
			"--  sqintest: if sfinal >= wf+3 generate\n"
			"--    squarerIn <= Zfinal(sfinal-1 downto pfinal+1) when doRR='1'\n"
			"--                 else (absZ0s &  (sfinal-wF-4 downto 0 => '0'));  \n"
			"--  end generate sqintest;\n"
			"--  sqintest2: if sfinal < wf+3 generate\n"
			"--    squarerIn <= Zfinal(sfinal-1 downto pfinal+1) when doRR='1'\n"
			"--                 else absZ0s(wF-pfinal+1 downto wf+3-sfinal);  \n"
			"--  end generate sqintest2;\n"
			" \n"
			"  Z2o2_full <= (squarerIn * squarerIn);\n"
			"  Z2o2 <= Z2o2_full (2*(sfinal-pfinal)-1  downto sfinal-pfinal-1);\n"
			"\n"
			"  Log1p_normal  <=   Zfinal  -  ((sfinal-1 downto sfinal-pfinal-1  => '0') & (Z2o2(sfinal-pfinal downto 2)));\n"
			"\n"
			"  LogF_normal <=   almostLog + ((targetprec-1 downto sfinal => '0') & Log1p_normal);\n"
			"\n"
			"  absELog2_pad <=   absELog2 & (targetprec-wF-g-1 downto 0 => '0');       \n"
			"  LogF_normal_pad <= (wE-1  downto 0 => LogF_normal(targetprec-1))  & LogF_normal;\n"
			"  \n"
			"  Log_normal <=  absELog2_pad  + LogF_normal_pad when sR='0'  \n"
			"                else absELog2_pad - LogF_normal_pad;\n"
			"\n"
			"  lzc_norm_0 : " << unique_name << "_lzc_norm\n"
			"    generic map (w => wE+targetprec, n => lzc_size)\n"
			"    port map (i => Log_normal, z => E_normal, o => Log_normal_normd);\n"
			"\n"
			"\n"
			"  rshiftsmall: " << unique_name << "_rshift\n"
			"    generic map (w => sfinal-pfinal+1,  n => log2wF) \n"
			"    port map (i => Z2o2,\n"
			"              s => shiftval(log2wF-1 downto 0),\n"
			"              o => Z2o2_small_s);\n"
			"\n"
			"  -- send the MSB to position pfinal\n"
			"  Z2o2_small <=  (pfinal-1 downto 0  => '0') & Z2o2_small_s & (wF+g-sfinal downto 0  => '0') ;\n"
			"\n"
			"  -- mantissa will be either Y0-z^2/2  or  -Y0+z^2/2,  depending on sR  \n"
			"\n"
			"  Z_small <= (absZ0s & (pfinal+g-1 downto 0 => '0'));\n"
			"  Log_small  <=       Z_small -  Z2o2_small when (sR='0')\n"
			"                else  Z_small +  Z2o2_small;\n"
			"\n"
			"  -- Possibly subtract 1 or 2 to the exponent, depending on the LZC of Log_small\n"
			"  E0_sub <=      \"11\" when Log_small(wF+g+1) = '1'\n"
			"            else \"10\" when Log_small(wF+g+1 downto wF+g) = \"01\"\n"
			"            else \"01\" ;\n"
			"\n"
			"  E_small <=  \"0\" & (wE-2 downto 2 => '1') & E0_sub\n"
			"               - ((wE-1 downto log2wF => '0') & lzo) ;\n"
			"\n"
			"  Log_small_normd <= Log_small(wF+g+1 downto 2) when Log_small(wF+g+1)='1'\n"
			"             else Log_small(wF+g downto 1)  when Log_small(wF+g)='1'  -- remove the first zero\n"
			"             else Log_small(wF+g-1 downto 0)  ; -- remove two zeroes (extremely rare, 001000000 only)\n"
			"                                               \n"
			"  ER <= E_small when small='1'\n"
			"        else E0offset - ((wE-1 downto lzc_size => '0') & E_normal);\n"
			"  -- works only if wE > lzc_size approx log2wF, OK for usual exp/prec\n"
			"\n"
			"  Log_g  <=  Log_small_normd (wF+g-2 downto 0) & \"0\" when small='1'           -- remove implicit 1\n"
			"        else Log_normal_normd(wE+targetprec-2 downto wE+targetprec-wF-g-1 );  -- remove implicit 1\n"
			"\n"
			"  sticky <= '0' when Log_g(g-2 downto 0) = (g-2 downto 0 => '0') else\n"
			"            '1';\n"
			"  round <= Log_g(g-1) and (Log_g(g) or sticky);\n"
			"\n"
			"  -- use a trick: if round leads to a change of binade, the carry propagation\n"
			"  -- magically updates both mantissa and exponent\n"
			"  EFR <= (ER & Log_g(wF+g-1 downto g)) + ((wE+wF-1 downto 1 => '0') & round); \n"
			"\n"
			"\n"
			"  -- The smallest log will be log(1+2^{-wF}) \\approx 2^{-wF}\n"
			"  -- The smallest representable number is 2^{-2^(wE-1)} \n"
			"  -- Therefore, if \n"
			"--    underflow : if max(wE, log2(wF)+1) > wE generate\n"
			"--      ufl <=      '1' when (eR2(wE0-1) = '1') or (eR = (wE-1 downto 0 => '0'))\n"
			"--             else '0';\n"
			"--    end generate;\n"
			"\n"
			"--    no_underflow : if max(wE, log2(wE+wF)+2) = wE generate\n"
			"      ufl <= '0';\n"
			"--    end generate;\n"
			"\n"
			"  R(wE+wF+2 downto wE+wF) <= \"110\" when ((X(wE+wF+2) and (X(wE+wF+1) or X(wE+wF))) or (X(wE+wF+1) and X(wE+wF))) = '1' else\n"
			"                               \"101\" when X(wE+wF+2 downto wE+wF+1) = \"00\"                                                       else\n"
			"                               \"100\" when X(wE+wF+2 downto wE+wF+1) = \"10\"                                                       else\n"
			"                               \"00\" & sR when (((Log_normal_normd(wE+targetprec-1)='0') and (small='0')) or ( (Log_small_normd (wF+g-1)='0') and (small='1'))) or (ufl = '1') else\n"
			"                               \"01\" & sR;\n"
			"\n"
			"  R(wE+wF-1 downto 0) <=  EFR;\n"
			"\n"
			"end architecture;\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"-------------------------------------------------------------------------------\n"
			"-- Range Reduction box\n"
			"-------------------------------------------------------------------------------\n"
			"library ieee;\n"
			"use ieee.std_logic_1164.all;\n"
			"use ieee.std_logic_arith.all;\n"
			"use ieee.std_logic_unsigned.all;\n"
			"\n"
			"\n";

	//---------------Range reduction entity----------------------------------------------

	o <<   "entity " << unique_name << "_range_red is port ("<<endl;
	o <<   "          Y0 : in  std_logic_vector("<<wF+1<<" downto 0);"<<endl;
	o <<   "          A  : in std_logic_vector("<< a[0]-1 <<" downto 0);"<<endl;
	if(is_sequential()) 
		o << "          clk  : in std_logic;"<<endl;
	o <<   "          Z  : out std_logic_vector("<< s[stages+1]-1 <<" downto 0);"<<endl;
	o <<   "  almostLog  : out std_logic_vector("<< lt0->wOut-1 <<" downto 0)  );"<<endl;
	o <<   "end entity;"<<endl<<endl;


	o << "architecture arch of " << unique_name << "_range_red is " <<endl<<endl;
	// All the components for the tables
	ostringstream nameIT0;
	nameIT0 << unique_name << "_invtable0_" << wE << "_" <<wF;
	it0->outputComponent(o, nameIT0.str());

	ostringstream nameLT0;
	nameLT0 << unique_name << "_logtable0_" << wE << "_" <<wF;
	lt0->outputComponent(o, nameLT0.str());

	for(i=1; i<=stages; i++) {
		ostringstream name3;
		name3 << unique_name << "_logtable"<<i<<"_" << wE << "_" <<wF;
		lt[i]->outputComponent(o, name3.str()); 
	}


	for (i=0; i<= stages; i++) {
		o << "   signal       A"<<i<<":  std_logic_vector("<< a[i] - 1  <<" downto 0);"<<endl;
	}

	
	for (i=1; i<= stages; i++)
			o << "   signal       B"<<i<<":  std_logic_vector("<< s[i] - a[i] - 1  <<" downto 0);"<<endl;

	for (i=0; i<= stages+1; i++)
		o << "   signal Z"<<i<<", Z"<<i<<"_d:  std_logic_vector("<< s[i] - 1  <<" downto 0);"<<endl;

	for (i=1; i<= stages; i++) {
		o << "   signal    epsZ"<<i<<":  std_logic_vector("<< s[i]+p[i]+1  <<" downto 0);"<<endl;
	}
	// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>eps[i]
	// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
	// We need only this many bits of Z[i] to compute it.
	for (i=1; i<= stages; i++) {
		o << "   signal      ZM"<<i<<":  std_logic_vector("<< psize[i] - 1  <<" downto 0);"<<endl;
	}

	o << "   signal       P0:  std_logic_vector("<<  it0->wOut + s[0] - 1  <<" downto 0);"<<endl;

	for (i=1; i<= stages; i++) {
		o << "   signal       P"<<i<<":  std_logic_vector("<< psize[i]+a[i] - 1  <<" downto 0);"<<endl;
	}

	o << "   signal       L0:  std_logic_vector("<< lt0->wOut -1 <<" downto 0);"<<endl;

	for (i=1; i<= stages; i++)
			o << "   signal       L"<<i<<":  std_logic_vector("<< lt[i]->wOut -1 <<" downto 0);"<<endl;

	// Note that all the Si have the size of L0: absence of carry out is proven.
	for (i=1; i<= stages+1; i++)
		o << "   signal S"<<i<<", S"<<i<<"_d:  std_logic_vector("<< lt0->wOut-1 <<" downto 0);"<<endl;

	o << "   signal    InvA0:  std_logic_vector("<< a[0]  <<" downto 0);"<<endl;
	//  o << "   " <<endl;



	o << "begin" <<endl;
	//  o << "   A0 <= Fx (wF-1 downto wF-"<<a[0]<<");" <<endl;
	o << "   A0 <= A;"<<endl;
	o << "   it0:"<<nameIT0.str()<<" port map (x=>A0, y=>InvA0);" <<endl; 
	o << "   lt0:"<<nameLT0.str()<<" port map (x=>A0, y=>L0);"<<endl;
	o << "   P0 <= InvA0 * Y0;" <<endl <<endl;
	if(is_sequential()) 
	{
		o << "   -- Synchronization barrier 1 " <<endl;
		o << "   process(clk)  begin\n     if clk'event and clk='1' then"<<endl;
		o << "     Z1_d <= P0("<< s[1] -1<<" downto 0);"<<endl;
		o << "     S1_d <= L0;"<<endl;
		o << "     end if;\n   end process;"<<endl<<endl;
	}
	else
	{
		o << "   Z1_d <= P0("<< s[1] -1<<" downto 0);"<<endl;
		o << "   S1_d <= L0;"<<endl;
	}

	for (i=1; i<= stages; i++) {

		o <<endl;
			//computation
		o << "   A"<<i<<" <= Z"<<i<<"_d(" << s[i] - 1  <<" downto "<< s[i] - a[i]  << ");"<<endl;
		o << "   B"<<i<<" <= Z"<<i<<"_d(" << s[i] - a[i] - 1  <<" downto 0 );"<<endl;
		o << "   lt"<<i<<":" << unique_name << "_logtable"<<i<<"_"<< wE <<"_"<<wF<<" port map (x=>A"<<i<<", y=>L"<<i<<");"<<endl;
		if(psize[i] == s[i])
			o << "   ZM"<<i<<" <= Z"<<i<< "_d;"<<endl;   
		else
			o << "   ZM"<<i<<" <= Z"<<i<<"_d(" <<s[i]-1   <<" downto "<< s[i]-psize[i]  << ");"<<endl;   
		o << "   P"<<i<<" <= A"<<i<<"*ZM"<<i<<";"<<endl;

		if(i==1) // special case for the first iteration
			{
				o << "   epsZ"<<i<<" <= ("<<s[i]+p[i]+1<<" downto 0 => '0') "
				     << "     when  A1 = ("<<a[1]-1<<" downto 0 => '0')"<<endl
				     << "       else (\"01\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"_d )"
				     << "  when ((A1("<<a[1]-1<<")='0') and (A1("<<a[1]-2<<" downto 0) /= ("<<a[1]-2<<" downto 0 => '0')))"<<endl
				     << "       else "
				     << "(\"1\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"_d  & \"0\") "
				     << ";"<<endl;
			}
		else 
			{
				o << "   epsZ"<<i<<" <=  ("<< s[i]+p[i]+1<<" downto 0 => '0') "
				     << "     when  A"<<i<<" = ("<<a[i]-1<<" downto 0 => '0')"<<endl
				     << "     else    (\"01\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"_d);"<<endl;
			}

		o << "   Z"<<i+1<<" <=   (\"0\" & B"<<i;
		if (s[i+1] > 1+(s[i]-a[i]))  // need to padd Bi
			o << " & ("<<s[i+1] - 1-(s[i]-a[i]) -1<<" downto 0 => '0') ";    
		o <<")"<<endl 
				 << "         - ( ("<<p[i]-a[i]<<" downto 0 => '0') & P"<<i;
		// either pad, or truncate P
		if(p[i]-a[i]+1  + psize[i]+a[i]  < s[i+1]) // size of leading 0s + size of p 
			 o << " & ("<<s[i+1] - (p[i]-a[i]+1  + psize[i]+a[i]) - 1 <<" downto 0 => '0')";  // Pad
		if(p[i]-a[i]+1  + psize[i]+a[i]  > s[i+1]) 
			//truncate
			o <<"("<< psize[i]+a[i] - 1  <<" downto "<<  p[i]-a[i]+1  + psize[i]+a[i]  - s[i+1] << " )";
		o << "  )"<< endl;

		o << "         + epsZ"<<i << "("<<s[i]+p[i]+1<<" downto "<<s[i]+p[i] +2 - s[i+1]<<")"
				 << ";"<<endl;
			

		o << "   S"<<i+1<<" <=   S"<<i<<"_d + (("<<lt0->wOut-1<<" downto "<<lt[i]->wOut<<" =>'0') & L"<<i<<");"<<endl;

		o << endl;
		if (is_sequential()) 
			{
				o << "   -- Synchronization barrier "<<i+1 <<endl;
				o << "   process(clk)  begin\n     if clk'event and clk='1' then"<<endl;
				o << "     Z"<<i+1<<"_d <=   Z"<<i+1<< ";"<<endl;
				o << "     S"<<i+1<<"_d <=   S"<<i+1<< ";"<<endl;
				o << "     end if;\n   end process;"<<endl<<endl;
			}
		else
			{
				o << "   Z"<<i+1<<"_d <=   Z"<<i+1<< ";"<<endl;
				o << "   S"<<i+1<<"_d <=   S"<<i+1<< ";"<<endl<<endl;
			}
	}


	o << "   Z <= Z"<<stages+1<<"_d;"<<endl;  
	o << "   almostLog <= S"<<stages+1<<"_d;"<<endl;  

	o << "end architecture;" <<endl<<endl;

	// All the tables 
	ostringstream  name4;
	name4 << unique_name << "_invtable0_" << wE << "_" <<wF;
	it0->output(o, name4.str());

	ostringstream name5;
	name5 << unique_name << "_logtable0_" << wE << "_" <<wF;
	lt0->output(o, name5.str());

	for(i=1; i<=stages; i++) {
		ostringstream name6;
		name6 << unique_name << "_logtable"<<i<<"_" << wE << "_" <<wF;
		lt[i]->output(o, name6.str()); 
	}
}

void FPLog::addTestCase(TestCaseList& tcl, FloFP x)
{
	Signal& sx = *get_signal_by_name("x");
	Signal& sr = *get_signal_by_name("r");
	Signal  sr_exc = (*get_signal_by_name("r")).getException();
	Signal  sr_sgn = (*get_signal_by_name("r")).getSign();
	Signal  sr_exp = (*get_signal_by_name("r")).getExponent();
	Signal  sr_man = (*get_signal_by_name("r")).getMantissa();

	int wE, wF;
	x.getPrecision(wE, wF);
	FloFP r(wE, wF);
	r = x.log();

	TestCase tc;
	tc.addInput(sx, x.getSignalValue());
	tc.addExpectedOutput(sr_exc, r.getExceptionSignalValue());
	tc.addExpectedOutput(sr_sgn, r.getSignSignalValue());
	if (r.getExceptionSignalValue() == 1)
	{
		/* FPLog only returns faithful rounding */
		tc.addExpectedOutput(sr, r.getRoundedDownSignalValue());
		tc.addExpectedOutput(sr, r.getRoundedUpSignalValue());
	}
	tcl.add(tc);
}

TestCaseList FPLog::generateStandardTestCases(int n)
{
	TestCaseList tcl;	/* XXX: Just like Lyon's Transportion Company. :D */
	FloFP x(wE, wF), r(wE, wF), zero(wE, wF), one(wE,wF);
	int i;

	/* Generate testcases for 0 ± 10ulp */
	for (x = 0.0, i = -10, x += i; i <= 10; i++, x++)
		addTestCase(tcl, x);

	/* Generate testcases for 1 ± 10ulp */
	for (x = 1.0, i = -10, x += i; i <= 10; i++, x++)
		addTestCase(tcl, x);

	/* Generate testcases for NaN ± 10ulp */
	for (x = NAN, i = -10, x += i; i <= 10; i++, x++)
		addTestCase(tcl, x);

	return tcl;
}

TestCaseList FPLog::generateRandomTestCases(int n)
{
	TestCaseList tcl;	/* XXX: Just like Lyon's Transportion Company. :D */
	FloFP x(wE, wF), r(wE, wF);

	for (int i = 0; i < n; i++)
	{
		x = getLargeRandom(wE+wF+1) + (mpz_class(1) << (wE + wF + 1));
		addTestCase(tcl, x);
	}

	return tcl;
}

