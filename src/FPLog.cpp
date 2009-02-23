/*
 * An FP logarithm for FloPoCo
 *
 * Author : Florent de Dinechin
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/
#include <fstream>
#include <sstream>
#include <math.h>	// for NaN
#include "FPLog.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"


using namespace std;

extern vector<Operator *> oplist;



FPLog::FPLog(Target* target, int wE, int wF)
	: Operator(target), wE(wE), wF(wF)
{

	setCopyrightString("F. de Dinechin, C. Klein  (2008)");

	ostringstream o;
	o << "FPLog_" << wE << "_" << wF;
	setName(o.str());

	setOperatorType();

	if(isSequential())
		cerr << endl << "************ WARNING: pipeline does not work yet" << endl << endl;

	addFPInput("X", wE, wF);
	addFPOutput("R", wE, wF, 2); // 2 because faithfully rounded

	int i;

	gLog = 5; // TODO : compute or tabulate 
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

	// Various variables that make life more concise
	sfinal =  s[stages+1];
	pfinal = p[stages+1];




	// On we go

	// MSB of squarer input is p[stages+1];
	// LSB will be target_prec
	// size will be target_prec - p[stages+1]  

	if(verbose)
		cout<<"FPLog: needs "<<stages<<" range reduction stages"<<endl;

	//TODO move somewhere -- removed temporarily to suppress a warning
	// int computedG = intlog2(3*(stages+1));


	vhdl << tab << declare("XExnSgn", 3) << " <=  X(wE+wF+2 downto wE+wF);";
	vhdl << tab << declare("FirstBit") << " <=  X(wF-1);" << endl;
	vhdl << tab << 	declare("sR") << " <= '0'   when  X(wE+wF-1 downto wF)   =   '0' &(wE-2 downto 0 => '1')  -- binade [1..2)" << endl
		  << "        else not X(wE+wF-1);                -- MSB of exponent" << endl;
	vhdl << tab << declare("Y0", wF+2) << " <=      \"1\"  & X(wF-1 downto 0) & \"0\" when FirstBit = '0'" 
		  << 		"        else \"01\" & X(wF-1 downto 0);" << endl;
	if(isSequential()) vhdl << tab << "-- Rem: the Y0 input is registered inside the RangeRed box" << endl;

	vhdl << tab << declare("absZ0", wF-pfinal+2) << " <=   Y0(wF-pfinal+1 downto 0)          when (" << use("sR") << "='0') else" << endl
		  << "             ((wF-pfinal+1 downto 0 => '0') - Y0(wF-pfinal+1 downto 0));" << endl;


	vhdl << tab << declare("E", wE) << " <= (X(wE+wF-1 downto wF)) - (\"0\" & (wE-2 downto 1 => '1') & (not FirstBit));" << endl;

	nextCycle();//////////////////////////////////////
	vhdl << tab << declare("absE", wE) << " <= ((wE-1 downto 0 => '0') - "<< use("E") << ")   when " << use("sR") << " = '1'" << endl
		  << "          else "<< use("E") << ";" << endl;
	vhdl << tab << declare("EeqZero") << " <= '1' when "<< use("E") << "=(wE-1 downto 0 => '0') else '0';" << endl;
	nextCycle();//////////////////////////////////////

	mpfr_t two, log2;
	mpz_t zlog2;

	// The log2 constant
	mpfr_init2(two, 2);
	mpfr_set_d(two, 2.0, GMP_RNDN);
	mpfr_init2(log2, wF+gLog);
	mpfr_log(log2, two, GMP_RNDN);
	mpfr_mul_2si(log2, log2, wF+gLog, GMP_RNDN); // shift left
	mpz_init2(zlog2, wF+gLog);
	mpfr_get_z(zlog2, log2, GMP_RNDN);
	vhdl << tab << declare("log2", wF+gLog) << " <= \"" << unsignedBinary(mpz_class(zlog2), wF+gLog) << "\";"<<endl;
	// TODO replace with a KCM
	vhdl << tab << declare("absELog2", wF+wE+gLog) << " <= " << use("absE") << " * log2;" << endl;

	// Back to cycle 1, after the 1-bit shift
	setCycle(1);//////////////////////////////////////

	lzoc = new LZOC(target, wF); 
	oplist.push_back(lzoc);
	
	vhdl << tab << declare("Y0h", wF) << " <= Y0(wF downto 1);" << endl; 
	inPortMap(lzoc, "I", "Y0h");
	inPortMap(lzoc, "OZB", "FirstBit");
	outPortMap(lzoc, "O", "lzo"); 
	vhdl << instance(lzoc, "lzoc1");
	syncCycleFromSignal("lzo");
	nextCycle();//////////////////////////////////////

	vhdl << tab << declare("pfinal_s", intlog2(wF)) << " <= \"" 
		  << unsignedBinary(mpz_class(p[stages+1]), intlog2(wF)) << "\";"<<endl;

	vhdl << tab << declare("shiftval", intlog2(wF)+1) << " <= ('0' & " << use("lzo") << ") - ('0' & pfinal_s); " << endl;
	vhdl << tab << declare("shiftvalin", intlog2(wF-pfinal+2)) 
		  << " <= shiftval(" << intlog2(wF-pfinal+2)-1 << " downto 0);" << endl;
	vhdl << tab << declare("doRR") << " <= shiftval(log2wF);             -- sign of the result" << endl;

	nextCycle();//////////////////////////////////////
	vhdl << tab << declare("small") << " <= " << use("EeqZero") << " and not " << use("doRR") << ";" << endl;


	// ao stands for "almost one"
	ao_lshift = new Shifter(target, wF-pfinal+2,  wF-pfinal+2, Left);   
	oplist.push_back(ao_lshift);

	inPortMap(ao_lshift, "X", "absZ0");
	inPortMap(ao_lshift, "S", "shiftvalin");
	outPortMap(ao_lshift, "R", "absZ0s");
	vhdl << instance(ao_lshift, "ao_lshift");



	//////////////////////////////////////////////
	setCycle(0);
	vhdl << tab << "-- The range reduction instance" << endl;

	rrbox = new LogRangeRed(target, this);
	oplist.push_back(rrbox);
	vhdl << tab << declare("rrA", a[0]) << " <= X(wF-1 downto wF-a0);" << endl;
	inPortMap(rrbox, "A", "rrA");
	inPortMap(rrbox, "Y0", "Y0");
	outPortMap(rrbox, "Z", "Zfinal");
	outPortMap(rrbox, "almostLog", "almostLog");
	vhdl << instance(rrbox, "rr");

	// Synchro between RR box and case almost 1
	setCycleFromSignal("Zfinal", false);	
	syncCycleFromSignal("absZ0s", false);
	nextCycle(); ///////////////////// TODO check this one is useful

	int absZ0sSize = getSignalByName("absZ0s")->width();
	vhdl << tab << "-- Z2o2 will be of size sfinal-pfinal, set squarer input size to that" << endl;
	vhdl << tab << declare("squarerIn", sfinal-pfinal) << " <= " 
		  << use("Zfinal") << "(sfinal-1 downto pfinal) when " << use("doRR") << "='1'" << endl;
	if(sfinal-pfinal>absZ0sSize)
		vhdl << tab << "                 else (" << use("absZ0s") << " & " << rangeAssign(sfinal-pfinal-absZ0sSize-1, 0, "'0'") << ");  " << endl;
	else  // sfinal-pfinal <= absZ0sSize
		vhdl << tab << "                 else " << use("absZ0s") << "" << range(absZ0sSize-1, absZ0sSize - (sfinal-pfinal)) << ";  " << endl<< endl;
	vhdl << tab << "-- Z2o2 will be of size sfinal - pfinal -1, set squarer input size to that" << endl;
	nextCycle(); ///////////////////// 
	vhdl << tab << declare("Z2o2_full", 2*(sfinal-pfinal)) << " <= (" << use("squarerIn") << " * " << use("squarerIn") << ");" << endl;
	vhdl << tab << declare("Z2o2", sfinal-pfinal+1) << " <= Z2o2_full (2*(sfinal-pfinal)-1  downto sfinal-pfinal-1);" << endl;
	nextCycle(); ///////////////////// 
	vhdl << tab << declare("Log1p_normal", sfinal) << " <=   " << use("Zfinal") << "  -  ((sfinal-1 downto sfinal-pfinal-1  => '0') & (" << use("Z2o2") << "(sfinal-pfinal downto 2)));" << endl;
	nextCycle(); ///////////////////// 
	vhdl << tab << declare("LogF_normal", target_prec) << " <=   " << use("almostLog") << " + ((targetprec-1 downto sfinal => '0') & Log1p_normal);" << endl;
	nextCycle(); ///////////////////// 
	vhdl << tab << declare("absELog2_pad", wE+target_prec) << " <=   " << use("absELog2") << " & (targetprec-wF-g-1 downto 0 => '0');       " << endl;
	vhdl << tab << declare("LogF_normal_pad", wE+target_prec) << " <= (wE-1  downto 0 => " << use("LogF_normal") << "(targetprec-1))  & " << use("LogF_normal") << ";" << endl;
	vhdl << tab << declare("Log_normal", wE+target_prec) << " <=  absELog2_pad  + LogF_normal_pad when " << use("sR") << "='0'  " << endl
		  << "                else absELog2_pad - LogF_normal_pad;" << endl;
	nextCycle(); ///////////////////// 

	final_norm = new LZOCShifterSticky(target, wE+target_prec, wE+target_prec, false, 0);
	oplist.push_back(final_norm);
	inPortMap(final_norm, "I", "Log_normal");
	outPortMap(final_norm, "Count", "E_normal");
	outPortMap(final_norm, "O", "Log_normal_normd");
	vhdl << instance(final_norm, "final_norm");

	// back to the squarer output
	setCycleFromSignal("Z2o2", false);
	nextCycle(); ///////////////////// 

	ao_rshift = new Shifter(target, s[stages+1]-p[stages+1]+1, s[stages+1]-p[stages+1]+1, Right) ;
	oplist.push_back(ao_rshift);
	inPortMap(ao_rshift, "X", "Z2o2");
	inPortMap(ao_rshift, "S", "shiftvalin");
	outPortMap(ao_rshift, "R", "Z2o2_small_s");
	vhdl << instance(ao_rshift, "ao_rshift");

	setCycleFromSignal("Z2o2_small_s", false);
	nextCycle(); ///////////////////// 

	vhdl << tab << "  -- send the MSB to position pfinal" << endl;
	int Z2o2_small_sSize = getSignalByName("Z2o2_small_s")->width();
	vhdl << tab << declare("Z2o2_small", wF+gLog+2) << " <=  (pfinal-1 downto 0  => '0') & Z2o2_small_s"
		  << range(Z2o2_small_sSize-1, Z2o2_small_sSize- (wF+gLog+2) + pfinal) << ";" << endl;

	vhdl << tab << "-- mantissa will be either Y0-z^2/2  or  -Y0+z^2/2,  depending on sR  " << endl;
	vhdl << tab << declare("Z_small", wF+gLog+2) << " <= " << use("absZ0s") << " & " << rangeAssign((wF+gLog+2)-absZ0sSize-1, 0, "'0'") << ";" << endl;
	vhdl << tab << declare("Log_small", wF+gLog+2) << " <=       Z_small -  Z2o2_small when (" << use("sR") << "='0')" << endl
		  << "              else  Z_small +  Z2o2_small;" << endl;


	nextCycle(); ///////////////////// 
	vhdl << tab << "-- Possibly subtract 1 or 2 to the exponent, depending on the LZC of " << use("Log_small") << "" << endl;
	vhdl << tab << declare("E0_sub", 2) << " <=      \"11\" when " << use("Log_small") << "(wF+g+1) = '1'" << endl
		  << "          else \"10\" when " << use("Log_small") << "(wF+g+1 downto wF+g) = \"01\"" << endl
		  << "          else \"01\" ;" << endl;
	vhdl << tab << declare("E_small", wE) << " <=  \"0\" & (wE-2 downto 2 => '1') & E0_sub" << endl
		  << "             - ((wE-1 downto log2wF => '0') & " << use("lzo") << ") ;" << endl;
	vhdl << tab << declare("Log_small_normd", wF+gLog) << " <= " << use("Log_small") << "(wF+g+1 downto 2) when " << use("Log_small") << "(wF+g+1)='1'" << endl
		  << "           else " << use("Log_small") << "(wF+g downto 1)  when " << use("Log_small") << "(wF+g)='1'  -- remove the first zero" << endl
		  << "           else " << use("Log_small") << "(wF+g-1 downto 0)  ; -- remove two zeroes (extremely rare, 001000000 only)" << endl ;

	setCycleFromSignal("E_normal", false);
	syncCycleFromSignal("E_small", false);
	nextCycle(); ///////////////////// 

	int E_normalSize = getSignalByName("E_normal")->width(); 
	vhdl << declare("E0offset", wE) << " <= \"" << unsignedBinary((mpz_class(1)<<(wE-1)) -2 + wE , wE) << "\"; -- E0 + wE "<<endl;
	vhdl << tab << declare("ER", wE) << " <= " << use("E_small") << " when " << use("small") << "='1'" << endl
		  << "      else E0offset - (" << rangeAssign(wE-1,  E_normalSize, "'0'") << " & " << use("E_normal") << ");" << endl
		  << "-- works only if wE > E_normalSize approx log2wF, OK for usual exp/prec" << endl;


	vhdl << tab << declare("Log_g", wF+gLog) << " <=  " << use("Log_small_normd") << "(wF+g-2 downto 0) & \"0\" when " << use("small") << "='1'           -- remove implicit 1" << endl
		  << "      else " << use("Log_normal_normd") << "(wE+targetprec-2 downto wE+targetprec-wF-g-1 );  -- remove implicit 1" << endl ;
	vhdl << tab << declare("sticky") << " <= '0' when Log_g(g-2 downto 0) = (g-2 downto 0 => '0')    else '1';" << endl;
	vhdl << tab << declare("round") << " <= Log_g(g-1) and (Log_g(g) or sticky);" << endl;
	vhdl << tab << "-- if round leads to a change of binade, the carry propagation" << endl
		  << tab << "-- magically updates both mantissa and exponent" << endl;
	// TODO an IntAdder here ?
	vhdl << tab << declare("EFR", wE+wF) << " <= (ER & Log_g(wF+g-1 downto g)) + ((wE+wF-1 downto 1 => '0') & round); " << endl;

	nextCycle(); ///////////////////// 
	vhdl << tab <<	"-- The smallest log will be log(1+2^{-wF}) \\approx 2^{-wF}" << endl
		  << tab << "-- The smallest representable number is 2^{-2^(wE-1)} " << endl;
	vhdl << tab << declare("ufl") << " <= '0';" << endl;
	vhdl << tab << "R(wE+wF+2 downto wE+wF) <= \"110\" when ((" << use("XExnSgn") << "(2) and (" << use("XExnSgn") << "(1) or " << use("XExnSgn") << "(0))) or (" << use("XExnSgn") << "(1) and " << use("XExnSgn") << "(0))) = '1' else" << endl
		  << "                               \"101\" when " << use("XExnSgn") << "(2 downto 1) = \"00\"                                                       else" << endl
		  << "                               \"100\" when " << use("XExnSgn") << "(2 downto 1) = \"10\"                                                       else" << endl
		  << "                               \"00\" & " << use("sR") << " when (((" << use("Log_normal_normd") << "(wE+targetprec-1)='0') and (" << use("small") << "='0')) or ( (" << use("Log_small_normd") << " (wF+g-1)='0') and (" << use("small") << "='1'))) or (ufl = '1') else" << endl
		  << "                               \"01\" & " << use("sR") << ";" << endl;
	vhdl << tab << "R(wE+wF-1 downto 0) <=  "<< use("EFR") << ";" << endl;
}	

FPLog::~FPLog()
{
#if 0 // probably to be performed on oplist
	delete lzoc;
	delete ao_rshift;
	delete ao_lshift;
	delete rrbox;
	delete final_norm;
#endif
}



void FPLog::outputVHDL(std::ostream& o, std::string name)
{
	licence(o);
	stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLComponentDeclarations();	
	o << buildVHDLSignalDeclarations();
	// constants
	o <<   "  constant g   : positive := "<<gLog<<";"<<endl;
	o <<   "  constant a0 : positive := "<< a[0] <<";"<<endl;
	o <<   "  constant wE : positive := " << wE <<";"<<endl;
	o <<   "  constant wF : positive := " << wF <<";"<<endl;
	o <<   "  constant log2wF : positive := "<< intlog2(wF) <<";"<<endl;
	o <<   "  constant targetprec : positive := "<< target_prec <<";"<<endl;
	o <<   "  constant sfinal : positive := "<< s[stages+1] <<";"<<endl;
	o <<   "  constant pfinal : positive := "<< p[stages+1] <<";"<<endl;
	beginArchitecture(o);		
	o<<buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);
}





void FPLog::emulate(TestCase * tc)
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
	fpx.getMPFR(x, false); // No fake zero here
	mpfr_log(rd, x, GMP_RNDD);
	mpfr_log(ru, x, GMP_RNDU);
#if 0
	mpfr_out_str (stderr, 10, 30, x, GMP_RNDN); cerr << " ";
	mpfr_out_str (stderr, 10, 30, rd, GMP_RNDN); cerr << " ";
	mpfr_out_str (stderr, 10, 30, ru, GMP_RNDN); cerr << " ";
	cerr << endl;
#endif
	FPNumber  fprd(wE, wF, rd);
	FPNumber  fpru(wE, wF, ru);
	mpz_class svRD = fprd.getSignalValue();
	mpz_class svRU = fpru.getSignalValue();
	tc->addExpectedOutput("R", svRD);
	tc->addExpectedOutput("R", svRU);
	mpfr_clears(x, ru, rd, NULL);
}
 

// TEST FUNCTIONS


void FPLog::buildStandardTestCases(TestCaseList* tcl){
	TestCase *tc;

	tc = new TestCase(this); 
	tc->addInput("X", 1.0);
	emulate(tc);
	tcl->add(tc);


}



// One test out of 8 fully random (tests NaNs etc)
// All the remaining ones test positive numbers.
// with special treatment for exponents 0 and -1.
 
void FPLog::buildRandomTestCases(TestCaseList* tcl, int n){

	TestCase *tc;
	mpz_class a;

	for (int i = 0; i < n; i++) {
		tc = new TestCase(this); 
		/* Fill inputs */
		if ((i & 7) == 0)
			a = getLargeRandom(wE+wF+3);
		else if ((i & 7) == 1) // exponent of 1
			a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-1)) << wF); 
		else if ((i & 7) == 2) // exponent of 0.5
			a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-2)) << wF); 
		else
			a  = getLargeRandom(wE+wF) + (mpz_class(1)<<(wE+wF+1)); // 010xxxxxx
		
		tc->addInput("X", a);
		/* Get correct outputs */
		emulate(tc);

		//		cout << tc->getInputVHDL();
		//    cout << tc->getExpectedOutputVHDL();
		// add to the test case list
		tcl->add(tc);
	}
}

