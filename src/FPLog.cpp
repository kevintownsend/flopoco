/*
  An FP logarithm for FloPoCo
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

*/

// TODO List: 
//  * test cases for boundary cases pfinal etc 
//  * finetune pipeline
//  * Port back the Arith paper
#include <fstream>
#include <sstream>
#include <math.h>	// for NaN
#include "FPLog.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"
#include "IntSquarer.hpp"
#include "fplogarithm/FirstInvTable.hpp"
#include "fplogarithm/FirstLogTable.hpp"
#include "fplogarithm/OtherLogTable.hpp"
#include "ConstMult/IntIntKCM.hpp"


using namespace std;


namespace flopoco{

	extern vector<Operator *> oplist;



	FPLog::FPLog(Target* target, int wE, int wF, int inTableSize)
		: Operator(target), wE(wE), wF(wF)
	{

		setCopyrightString("F. de Dinechin, C. Klein  (2008)");

		ostringstream o;

		o << "FPLog_" << wE << "_" << wF << "_" << inTableSize << "_";
		if(target->isPipelined()) 
			o << target->frequencyMHz() ;
		else
			o << "comb";
		setName(o.str());

		addFPInput("X", wE, wF);
		addFPOutput("R", wE, wF, 2); // 2 because faithfully rounded

		int i, bitsPerStage;

		if (inTableSize==0) {
			if (wF>=20)
				bitsPerStage=12;
			else
				bitsPerStage=8; // somehow arbitrary
		}
		else
			bitsPerStage=inTableSize;

		if(bitsPerStage <6)
			throw string("FPLog error: tables need at least 6 input bits");



		// First compute the parameters of each iteration 

		// stage 0, needs a specific inverter table
		p[0] = 0;
		a[0] = bitsPerStage;   
		// How many bits have been zeroed?  a0=5 is a lucky case but it messes up the following
		p[1] = a[0]-2; // if a0==5, would be a0-1

		// The number of guard bit
		// For each stage: 0.5 ulp for rounding the log table
		// For each stage after the first one:
		//                 1 ulp for truncating Z
		//                 1 ulp for truncating P
		//                 1 ulp for truncating EiY
		//        total 3.5 bit/stage
		// Plus for the rest of the computation:
		//                 1 ulp for the approximation error of neglecting the 3rd-order term
		//                 1 ulp for truncating Z before input to the square
		//                 1 ulp for truncating Z^2

		// the iterations i=0 and i=1 lead to no truncation

		i=1;
		gLog=4;
		while(3*p[i]+1 <= p[i]+wF+gLog){ // while the third-order term of log(1+Zi) is not negligible
			if(i==1) { 	// stage 1 cannot be more accurate than 2p1-1
				a[1] = p[1];  
				p[2] = 2*p[1]-1;
			}
			else {
				a[i] = bitsPerStage;
				p[i+1] = p[i] + a[i] - 1; // and we zero out a[i]-1 bits
			}
			i++;
			gLog=max(4, intlog2(3+0.5+0.5+3*i-1));
		}  
	
		// The number of stages, not counting stage 0
		stages = i-1;
		gLog=max(4, intlog2(3+0.5+0.5+3*stages));
	
		if(verbose>=2) {
			cerr << "> FPLog\t Initial parameters:" << endl;
			for(i=0; i<=stages; i++) {
				cerr << "> FPLog\t";
				cerr<<"\tp"<<i<<"=" << p[i];
				cerr<<"\ta"<<i<<"=" << a[i];
				cerr <<endl;
			}
		}
		// Now we probably have reduced too far
		pfinal = p[stages+1];
		int extraBits = pfinal - ((wF+gLog-2)>>1);
		int extraBitsperstage =  floor(((double)extraBits) / ((double)(stages+1)));
		if(verbose)
			cerr << "> FPLog\t before optimization:  pfinal=" << pfinal << "--- extraBits=" << extraBits << "---extraBitsperstage=" << extraBitsperstage << endl;
		if(bitsPerStage>6) {
			for (i=0; i<= stages; i++)
				a[i] -= extraBitsperstage;
		}
		else  {
			for (i=2; i<= stages; i++)
				a[i] -= extraBitsperstage;
		}
		// Recompute the pi
		p[1] = a[0]-2;
	
		for(i=1; i<=stages; i++){ // for faithful rounding
			if(i==1)   	// stage 1 cannot be more accurate than 2p1-1
				p[2] = p[1] + a[1] -1; // 2p[1]-1
			else
				p[i+1] = p[i] + a[i] - 1;
		}
		pfinal = p[stages+1];
		extraBits = pfinal -  ((wF+gLog-2)>>1);
		if(stages>=2) { // remove extra bits from stages 2 to stages
			extraBitsperstage =  floor(((double)extraBits) / ((double)(stages-1)));
			int extraBitsStage2 = extraBits - (stages-1)*extraBitsperstage;
			a[2] -= extraBitsperstage + extraBitsStage2;
			p[2+1] = p[2] + a[2] - 1;
			for(i=3; i<= stages; i++) {
				a[i] -= extraBitsperstage;
				p[i+1] = p[i] + a[i] - 1;
			} 
			pfinal = p[stages+1];
		}

		extraBits = pfinal -  ((wF+gLog-2)>>1);

		if(verbose>=2)
			cerr << "> FPLog\t after optimization:   pfinal=" << pfinal << "--- extraBits=" << extraBits << endl;
	
 
		if(verbose)
			cerr << "> FPLog"<<tab<<"Guard bits: " << gLog << endl;


		target_prec = wF + pfinal +gLog;
		if(verbose==2)
			cerr << "> FPLog"<<tab<<"Target precision: " << target_prec << endl;

		s[0] = wF+2;  
		psize[0] = s[0] + a[0]+1;
		//	sfullZ[0] = wF+2;
		sbt[1] = psize[0] - p[1] -2; // -2 because P0 begins with 01   was: wF+2 ;
		s[1] = sbt[1];
		t[1] = 0;
		//	sfullZ[1] = sfullZ[0]      +    a[0] + 1;  
		//         size of Y0          size of approx inv of A0

		for(i=1; i<=stages; i++) { 
			// size of Z before truncation; Zi = 0Bi - 0AiZi + EiYi ; 
			// AiZi has weight -2*p[i], and size a[i]+s[i]
			// EiYi has weight -2*p[i]-1 and size 1+p[i]+s[i]
			// except i=1: same weight but size may be 1+p[i]+s[i]+1
			if(i==1)
				sbt[i+1] =    max( a[i]+s[i]+1, 1+p[i]+sbt[i]+1);
			else
				sbt[i+1] =    max( a[i]+s[i]+1, 1+p[i]+sbt[i]);

			if(p[i+1]+sbt[i+1] <= target_prec) 
				{ // No truncation at all
					psize[i] = s[i];
					s[i+1] = sbt[i+1];
					t[i+1] = 0;
				}
			else
				{ // Truncate everybody to targetprec : 
					// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>Ei[i]
					// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
					// We need only this many bits of Z[i] to compute it.
					psize[i] = target_prec - 2*p[i];
					if (psize[i]>s[i]) // in the first iterations
						psize[i]=s[i];
					s[i+1] = target_prec - p[i+1];
					t[i+1] = sbt[i+1] - s[i+1];
				}
		}

		sfinal =  s[stages+1];

		// MSB of squarer input is p[stages+1];
		// LSB will be target_prec
		// size will be target_prec - p[stages+1]  


		if(verbose)
			cerr<<"> FPLog\t needs 1+"<<stages<<" range reduction stages"<<endl;
		if(verbose>=2) {
			for(i=0; i<=stages; i++) {
				cerr << "> FPLog\t";
				cerr<<"\tp"<<i<<"=" << p[i];
				cerr<<"\ta"<<i<<"=" << a[i];
				cerr<<"\ts"<<i<<"=" << s[i];
				cerr<<"\tpsize"<<i<<"=" << psize[i];
				cerr <<endl;
			}
			cerr << "> FPLog\t\tsfinal=" << sfinal << "\tpfinal=" << pfinal << endl;
		
		}


		// On we go with the vhdl


		addConstant("g",  "positive",          gLog);
		addConstant("wE", "positive",          wE);
		addConstant("wF", "positive",          wF);
		addConstant("log2wF", "positive",     intlog2(wF));
		addConstant("targetprec", "positive", target_prec);
		addConstant("sfinal", "positive",     s[stages+1]);
		addConstant("pfinal", "positive",     p[stages+1]);


		vhdl << tab << declare("XExnSgn", 3) << " <=  X(wE+wF+2 downto wE+wF);" << endl;
		vhdl << tab << declare("FirstBit") << " <=  X(wF-1);" << endl;
		vhdl << tab << 	declare("sR") << " <= '0'   when  X(wE+wF-1 downto wF)   =   '0' &(wE-2 downto 0 => '1')  -- binade [1..2)" << endl
			  << "        else not X(wE+wF-1);                -- MSB of exponent" << endl;
		vhdl << tab << declare("Y0", wF+2) << " <=      \"1\"  & X(wF-1 downto 0) & \"0\" when FirstBit = '0'" 
			  << 		"        else \"01\" & X(wF-1 downto 0);" << endl;
		if(isSequential()) vhdl << tab << "-- Rem: the Y0 input is registered inside the RangeRed box" << endl;

		vhdl << tab << declare("absZ0", wF-pfinal+2) << " <=   Y0(wF-pfinal+1 downto 0)          when (sR='0') else" << endl
			  << "             ((wF-pfinal+1 downto 0 => '0') - Y0(wF-pfinal+1 downto 0));" << endl;


		vhdl << tab << declare("E", wE) << " <= (X(wE+wF-1 downto wF)) - (\"0\" & (wE-2 downto 1 => '1') & (not FirstBit));" << endl;

		nextCycle();//////////////////////////////////////
		vhdl << tab << declare("absE", wE) << " <= ((wE-1 downto 0 => '0') - E)   when sR = '1'" << endl
			  << "          else E;" << endl;
		vhdl << tab << declare("EeqZero") << " <= '1' when E=(wE-1 downto 0 => '0') else '0';" << endl;
		nextCycle();//////////////////////////////////////

		// Back to cycle 1, after the 1-bit shift
		setCycle(1);//////////////////////////////////////

		lzoc = new LZOC(target, wF); 
		oplist.push_back(lzoc);
	
		vhdl << tab << declare("Y0h", wF) << " <= Y0(wF downto 1);" << endl; 
		inPortMap(lzoc, "I", "Y0h");
		inPortMap(lzoc, "OZB", "FirstBit");
		outPortMap(lzoc, "O", "lzo"); 
		vhdl << instance(lzoc, "lzoc1");
		setCycleFromSignal("lzo", false);
		nextCycle();//////////////////////////////////////

		vhdl << tab << declare("pfinal_s", intlog2(wF)) << " <= \"" 
			  << unsignedBinary(mpz_class(pfinal), intlog2(wF)) << "\";"<<endl;

		vhdl << tab << declare("shiftval", intlog2(wF)+1) << " <= ('0' & lzo) - ('0' & pfinal_s); " << endl;
		vhdl << tab << declare("shiftvalinL", intlog2(wF-pfinal+2)) 
			  << " <= shiftval(" << intlog2(wF-pfinal+2)-1 << " downto 0);" << endl;
		vhdl << tab << declare("shiftvalinR", intlog2(sfinal-pfinal+1)) 
			  << " <= shiftval(" << intlog2(sfinal-pfinal+1)-1 << " downto 0);" << endl;

		vhdl << tab << declare("doRR") << " <= shiftval(log2wF);             -- sign of the result" << endl;

		nextCycle();//////////////////////////////////////
		vhdl << tab << declare("small") << " <= EeqZero and not doRR;" << endl;


		// ao stands for "almost one"
		vhdl << tab << "-- The left shifter for the 'small' case" <<endl; 
		ao_lshift = new Shifter(target, wF-pfinal+2,  wF-pfinal+2, Shifter::Left);   
		oplist.push_back(ao_lshift);

		inPortMap(ao_lshift, "X", "absZ0");
		inPortMap(ao_lshift, "S", "shiftvalinL");
		outPortMap(ao_lshift, "R", "small_absZ0_normd_full");
		vhdl << instance(ao_lshift, "small_lshift");

		syncCycleFromSignal("small_absZ0_normd_full");

		int small_absZ0_normd_size = getSignalByName("small_absZ0_normd_full")->width() - (wF-pfinal+2);
		vhdl << tab << declare("small_absZ0_normd", small_absZ0_normd_size) << " <= small_absZ0_normd_full" << range(small_absZ0_normd_size -1, 0) << "; -- get rid of leading zeroes" << endl; 


		//////////////////////////////////////////////
		setCycle(0);
		vhdl << tab << "---------------- The range reduction box ---------------" << endl;



		vhdl << tab << declare("A0", a[0]) << " <= X" << range(wF-1,  wF-a[0]) << ";" << endl;

		vhdl << tab << "-- First inv table" << endl;
		FirstInvTable* it0 = new FirstInvTable(target, a[0], a[0]+1);
		oplist.push_back(it0);
		inPortMap       (it0, "X", "A0");
		outPortMap      (it0, "Y", "InvA0");
		vhdl << instance(it0, "itO");

		nextCycle();
		// TODO: somehow arbitrary
		if(target->isPipelined() && a[0]+1>=9) {  
			IntMultiplier* p0 = new IntMultiplier(target, a[0]+1, wF+2);
			oplist.push_back(p0);
			inPortMap  (p0, "X", "InvA0");
			inPortMap  (p0, "Y", "Y0");
			outPortMap (p0, "R", "P0");
			vhdl << instance(p0, "p0_mult");
		}
		else {
			if(verbose) cerr << "> FPLog: unpipelined multiplier for P0, implemented as * in VHDL" << endl;  
			vhdl << tab << declare("P0",  psize[0]) << " <= InvA0*Y0;" <<endl <<endl;
		}	

		setCycleFromSignal("P0", true);
		vhdl << tab << declare("Z1", s[1]) << " <= P0" << range (psize[0] - p[1]-3,  0) << ";" << endl;


		for (i=1; i<= stages; i++) {

			vhdl <<endl;
			//computation
			vhdl << tab << declare(join("A",i), a[i]) <<     " <= " << use(join("Z",i)) << range(s[i]-1,      s[i]-a[i]) << ";"<<endl;
			vhdl << tab << declare(join("B",i), s[i]-a[i]) << " <= " << use(join("Z",i)) << range(s[i]-a[i]-1, 0        ) << ";"<<endl;

			vhdl << tab << declare(join("ZM",i), psize[i]) << " <= " << use(join("Z",i)) ;
			if(psize[i] == s[i])
				vhdl << ";"<<endl;   
			else
				vhdl << range(s[i]-1, s[i]-psize[i])  << ";" << endl;   

			nextCycle();

			if(target->isPipelined() && a[0]+1>=9) {  
				IntMultiplier* pi = new IntMultiplier(target, a[i], psize[i]);
				oplist.push_back(pi);
				inPortMap  (pi, "X", join("A",i));
				inPortMap  (pi, "Y", join("ZM",i));
				outPortMap (pi, "R", join("P",i));
				vhdl << instance(pi, join("p",i)+"_mult");
			}
			else {
				if(verbose) cerr << "> FPLog: unpipelined multiplier for P"<<i<<", implemented as * in VHDL" << endl;  
				vhdl << tab << declare(join("P",i),  psize[i] + a[i]) << " <= " << use(join("A",i)) << "*" << use(join("ZM",i)) << ";" << endl;
			}


			int yisize = s[i]+p[i]+1; // +1 because implicit 1

			// While the multiplication takes place, we may prepare the other term 
			string Yi= join("Y",i);
			string Zi= join("Z",i);
			vhdl << tab << declare(Yi, yisize) << " <= " 
				  << " \"1\" & " << rangeAssign(p[i]-1,  0,  "'0'") << " & " << use(Zi) <<";"<<endl;

			// We truncate EiY to a position well above target_prec
			vhdl << tab << declare(join("EiY",i), s[i+1]) << " <= " ;
			if(i==1) { 
				vhdl <<  use(Yi) ;
				// now we may need to truncate or to pad Yi
				if(yisize > s[i+1]) // truncate Yi
					vhdl << range(yisize-1, yisize-s[i+1]);
				else if (yisize < s[i+1]) // pad Yi
					vhdl << " & " << rangeAssign(s[i+1] - yisize -1, 0, "'0'");
				vhdl  << "  when " << use(join("A",i)) << of(a[i]-1) << " = '1'" << endl
						<< tab << "  else  \"0\" & " << use(Yi);
				if(yisize > s[i+1]-1) // truncate Yi
					vhdl << range(yisize-1, yisize-s[i+1]+1);
				vhdl 	  << ";" << endl ;
			} 
			else { // i>1, general case
				vhdl << rangeAssign(2*p[i] -p[i+1] -2,  0,  "'0'")  << " & " << use(Yi) << range(yisize-1,  yisize - (s[i+1] - (2*p[i]-p[i+1]) )-1) << ";" << endl ;
			}
			nextCycle();
			// IntAdder here? TODO
			vhdl << tab << declare(join("EiYPB",i), s[i+1]) << " <= " 
				  << " (\"0\" & " << use(join("B",i));
			if (s[i+1] > 1+(s[i]-a[i]))  // need to padd Bi
				vhdl << " & " << rangeAssign(s[i+1] - 1-(s[i]-a[i]) -1,  0 , "'0'");    
			vhdl <<") + " << use(join("EiY",i)) << ";" << endl; 

			syncCycleFromSignal(join("P",i), false);
			nextCycle();


			vhdl << tab << declare(join("Pp", i), s[i+1])  << " <= " 
				  << rangeAssign(p[i]-a[i],  0,  "'0'") << " & " << use(join("P", i));
			// either pad, or truncate P
			if(p[i]-a[i]+1  + psize[i]+a[i]  < s[i+1]) // size of leading 0s + size of p 
				vhdl << " & "<< rangeAssign(s[i+1] - (p[i]-a[i]+1  + psize[i]+a[i]) - 1,    0,  "'0'");  // Pad
			if(p[i]-a[i]+1  + psize[i]+a[i]  > s[i+1]) 
				//truncate
				vhdl << range(psize[i]+a[i]-1,    p[i]-a[i]+1  + psize[i]+a[i] - s[i+1]);
			vhdl << ";"<<endl;

			vhdl << tab << declare(join("Z", i+1), s[i+1])  << " <= " 
				  << use(join("EiYPB",i)) << " - " << use(join("Pp", i)) << ";"<<endl;
		
		}

		vhdl << tab << declare("Zfinal", s[stages+1]) << " <= " << use(join("Z", stages+1)) << ";" << endl;  



		vhdl << tab << "--  Synchro between RR box and case almost 1" << endl;  
		setCycleFromSignal("Zfinal", false);	
		syncCycleFromSignal("small_absZ0_normd", false);
		nextCycle(); ///////////////////// TODO check this one is useful

		// In the small path we need Z2O2 accurate to  (wF+gLog+2) - pfinal
		// In the RR path we need Z2O2 accurate to sfinal-pfinal
		// Take the max. This is useful  for small precs only
		int squarerInSize;
		if((wF+gLog+2) - pfinal > sfinal-pfinal) 
			squarerInSize = (wF+gLog+2) - pfinal;
		else 
			squarerInSize = sfinal-pfinal;
				
		vhdl << tab << declare("squarerIn", squarerInSize) << " <= " 
			 << "Zfinal(sfinal-1 downto sfinal-"<< squarerInSize << ") when doRR='1'" << endl;
		if(squarerInSize>small_absZ0_normd_size)
			vhdl << tab << "                 else (small_absZ0_normd & " << rangeAssign(squarerInSize-small_absZ0_normd_size-1, 0, "'0'") << ");  " << endl;
		else  // sfinal-pfinal <= small_absZ0_normd_size
			vhdl << tab << "                 else small_absZ0_normd" << range(small_absZ0_normd_size-1, small_absZ0_normd_size - squarerInSize) << ";  " << endl<< endl;
		nextCycle(); /////////////////////
		if(squarerInSize<=17)
			vhdl << tab << declare("Z2o2_full", 2*squarerInSize) << " <= (squarerIn * squarerIn);" << endl;
		else{
			IntSquarer* sq = new IntSquarer(target, squarerInSize);
			oplist.push_back(sq);
			inPortMap  (sq, "X", "squarerIn");
			outPortMap (sq, "R", "Z2o2_full");
			vhdl << instance(sq, "squarer");
			setCycleFromSignal("Z2o2_full", true);
		}
		vhdl << tab << declare("Z2o2_normal", sfinal-pfinal-1) << " <= Z2o2_full ("<< 2*squarerInSize-1 << "  downto " << 2*squarerInSize - (sfinal-pfinal-1) << ");" << endl;
		nextCycle(); ///////////////////// 
		vhdl << tab << declare("Log1p_normal", sfinal) << " <=   Zfinal  -  ((pfinal downto 0  => '0') & " << use("Z2o2_normal")		<< ");" << endl;

		vhdl << endl << tab << "-- Now the log tables, as late as possible" << endl;
		setCycle(getCurrentCycle() - stages -1 , true);

		vhdl << tab << "-- First log table" << endl;
		FirstLogTable* lt0 = new FirstLogTable(target, a[0], target_prec, it0, this);
		oplist.push_back(lt0);
		inPortMap       (lt0, "X", "A0");
		outPortMap      (lt0, "Y", "L0");
		vhdl << instance(lt0, "ltO");
		vhdl << tab << declare("S1", lt0->wOut) << " <= L0;"<<endl;
		nextCycle();

		for (i=1; i<= stages; i++) {

			// TODO better pipeline the small input as late as possible than pipeline the large output
			OtherLogTable* lti = new OtherLogTable(target, a[i], target_prec - p[i], i, a[i], p[i]); 
			oplist.push_back(lti);
			inPortMap       (lti, "X", join("A", i));
			outPortMap      (lti, "Y", join("L", i));
			vhdl << instance(lti, join("lt",i));
			nextCycle();
			vhdl << tab << declare(join("S",i+1),  lt0->wOut) << " <= " 
				  <<  use(join("S",i))  << " + (" << rangeAssign(lt0->wOut-1, lti->wOut,  "'0'") << " & " << use(join("L",i)) <<");"<<endl;
		}

		vhdl << tab << declare("almostLog", lt0->wOut) << " <= " << use(join("S",stages+1)) << ";" << endl;  


		nextCycle(); ///////////////////// 
		vhdl << tab << declare("LogF_normal", target_prec) << " <=   almostLog + ((targetprec-1 downto sfinal => '0') & Log1p_normal);" << endl;
		// The log2 constant
		mpfr_t two, log2;
		mpz_t zlog2;
		mpfr_init2(two, 2);
		mpfr_set_d(two, 2.0, GMP_RNDN);
		mpfr_init2(log2, wF+gLog);
		mpfr_log(log2, two, GMP_RNDN);
		mpfr_mul_2si(log2, log2, wF+gLog, GMP_RNDN); // shift left
		mpz_init2(zlog2, wF+gLog);
		mpfr_get_z(zlog2, log2, GMP_RNDN);
#if 0
		// TODO replace with a KCM
		vhdl << tab << declare("log2", wF+gLog) << " <= \"" << unsignedBinary(mpz_class(zlog2), wF+gLog) << "\";"<<endl;
		vhdl << tab << declare("absELog2", wF+wE+gLog) << " <= absE * log2;" << endl;
		nextCycle(); ///////////////////// 
		vhdl << tab << declare("absELog2_pad", wE+target_prec) << " <=   absELog2 & (targetprec-wF-g-1 downto 0 => '0');       " << endl;
#else
		IntIntKCM* kcm=new IntIntKCM(target, wE, mpz_class(zlog2));
		oplist.push_back(kcm);
		// get back enough cycles to synchronize it with LogF_normal
		setCycle(getCurrentCycle() - kcm->getPipelineDepth());
		inPortMap       (kcm, "X", "absE");
		outPortMap      (kcm, "R", "absELog2");
		vhdl << instance(kcm, "Log2KCM");
		setCycle(getCurrentCycle() + kcm->getPipelineDepth(), false);
		nextCycle(); ///////////////////// 
		vhdl << tab << declare("absELog2_pad", wE+target_prec) << " <=   " 
			 << "absELog2 & (targetprec-wF-g-1 downto 0 => '0');       " << endl;		
#endif

		vhdl << tab << declare("LogF_normal_pad", wE+target_prec) << " <= (wE-1  downto 0 => LogF_normal(targetprec-1))  & LogF_normal;" << endl;
		vhdl << tab << declare("Log_normal", wE+target_prec) << " <=  absELog2_pad  + LogF_normal_pad when sR='0'  " << endl
			  << "                else absELog2_pad - LogF_normal_pad;" << endl;
		nextCycle(); ///////////////////// 

		final_norm = new LZOCShifterSticky(target, wE+target_prec, target_prec, intlog2(wE+(wF>>1))+1, false, 0); 
		oplist.push_back(final_norm);
		inPortMap(final_norm, "I", "Log_normal");
		outPortMap(final_norm, "Count", "E_normal");
		outPortMap(final_norm, "O", "Log_normal_normd");
		vhdl << instance(final_norm, "final_norm");

		// back to the squarer output
		setCycleFromSignal("Z2o2_full", false);
		nextCycle(); ///////////////////// 
		int Z2o2_small_size=(wF+gLog+2) - pfinal; // we need   (wF+gLog+2) - pfinal bits of Z2O2

		vhdl << tab << declare("Z2o2_small_bs", Z2o2_small_size)  << " <= Z2o2_full" << range(2*squarerInSize -1, 2*squarerInSize -Z2o2_small_size) << ";" << endl;
		ao_rshift = new Shifter(target, Z2o2_small_size, sfinal-pfinal+1, Shifter::Right) ;
		oplist.push_back(ao_rshift);
		inPortMap(ao_rshift, "X", "Z2o2_small_bs");
		inPortMap(ao_rshift, "S", "shiftvalinR");
		outPortMap(ao_rshift, "R", "Z2o2_small_s");
		vhdl << instance(ao_rshift, "ao_rshift");

		setCycleFromSignal("Z2o2_small_s", false);
		nextCycle(); ///////////////////// 

		vhdl << tab << "  -- send the MSB to position pfinal" << endl;
		int Z2o2_small_sSize = getSignalByName("Z2o2_small_s")->width();
		vhdl << tab << declare("Z2o2_small", wF+gLog+2) << " <=  (pfinal-1 downto 0  => '0') & " << use("Z2o2_small_s")
			  << range(Z2o2_small_sSize-1,  Z2o2_small_sSize - (wF+gLog+2) + pfinal) << ";" << endl;
	
		vhdl << tab << "-- mantissa will be either Y0-z^2/2  or  -Y0+z^2/2,  depending on sR  " << endl;
		vhdl << tab << declare("Z_small", wF+gLog+2) << " <= small_absZ0_normd & " << rangeAssign((wF+gLog+2)-small_absZ0_normd_size-1, 0, "'0'") << ";" << endl;
		vhdl << tab << declare("Log_small", wF+gLog+2) << " <=       Z_small -  Z2o2_small when (sR='0')" << endl
			  << "                else  Z_small +  Z2o2_small;" << endl;


		nextCycle(); ///////////////////// 
		vhdl << tab << "-- Possibly subtract 1 or 2 to the exponent, depending on the LZC of Log_small" << endl;
		vhdl << tab << declare("E0_sub", 2) << " <=   \"11\" when Log_small(wF+g+1) = '1'" << endl
			  << "          else \"10\" when Log_small(wF+g+1 downto wF+g) = \"01\"" << endl
			  << "          else \"01\" ;" << endl;
		// Is underflow possible?
		vhdl << tab <<	"-- The smallest log will be log(1+2^{-wF}) \\approx 2^{-wF}  = 2^" << -wF <<  endl
			  << tab << "-- The smallest representable number is 2^{1-2^(wE-1)} = 2^" << 1 -(1<<(wE-1))<< endl;
		if(1 -(1<<(wE-1)) < -wF) {
			vhdl << tab <<	"-- No underflow possible" <<  endl;
			vhdl << tab << declare("ufl") << " <= '0';" << endl;
			vhdl << tab << declare("E_small", wE) << " <=  (\"0\" & (wE-2 downto 2 => '1') & E0_sub)  -  ";
			if(wE>getSignalByName("lzo")->width())
				vhdl << "((wE-1 downto " << getSignalByName("lzo")->width() << " => '0') & lzo) ;" << endl;
			else
				vhdl << "lzo;" << endl;
		}
		else{
			vhdl << tab <<	"-- Underflow may happen" <<  endl;
			vhdl << tab << declare("E_small", wE+1) << " <=  (\"00\" & (wE-2 downto 2 => '1') & E0_sub)  -  (";
			vhdl << "'0' & lzo);" << endl;
			vhdl << tab << declare("ufl") << " <= E_small(wE);" << endl;
		}	
		vhdl << tab << declare("Log_small_normd", wF+gLog) << " <= Log_small(wF+g+1 downto 2) when Log_small(wF+g+1)='1'" << endl
			  << "           else Log_small(wF+g downto 1)  when Log_small(wF+g)='1'  -- remove the first zero" << endl
			  << "           else Log_small(wF+g-1 downto 0)  ; -- remove two zeroes (extremely rare, 001000000 only)" << endl ;

		setCycleFromSignal("E_normal", false);
		syncCycleFromSignal("E_small", false);
		nextCycle(); ///////////////////// 

		int E_normalSize = getSignalByName("E_normal")->width(); 
		vhdl << tab << declare("E0offset", wE) << " <= \"" << unsignedBinary((mpz_class(1)<<(wE-1)) -2 + wE , wE) << "\"; -- E0 + wE "<<endl;
		vhdl << tab << declare("ER", wE) << " <= E_small" << range(wE-1,0) << " when small='1'" << endl;
		if(wE>E_normalSize)
			vhdl << "      else E0offset - (" << rangeAssign(wE-1,  E_normalSize, "'0'") << " & E_normal);" << endl;
		else
			vhdl << "      else E0offset - E_normal;" << endl;
		vhdl << tab << declare("Log_g", wF+gLog) << " <=  Log_small_normd(wF+g-2 downto 0) & \"0\" when small='1'           -- remove implicit 1" << endl
			  << "      else Log_normal_normd(targetprec-2 downto targetprec-wF-g-1 );  -- remove implicit 1" << endl ;
		// Sticky is always 1 for a transcendental function !
		// vhdl << tab << declare("sticky") << " <= '0' when Log_g(g-2 downto 0) = (g-2 downto 0 => '0')    else '1';" << endl;
		vhdl << tab << declare("round") << " <= Log_g(g-1) ; -- sticky is always 1 for a transcendental function " << endl;
		vhdl << tab << "-- if round leads to a change of binade, the carry propagation magically updates both mantissa and exponent" << endl;
		// TODO an IntAdder here ?
		vhdl << tab << declare("EFR", wE+wF) << " <= (ER & Log_g(wF+g-1 downto g)) + ((wE+wF-1 downto 1 => '0') & round); " << endl;

		vhdl << tab << "R(wE+wF+2 downto wE+wF) <= \"110\" when ((XExnSgn(2) and (XExnSgn(1) or XExnSgn(0))) or (XExnSgn(1) and XExnSgn(0))) = '1' else" << endl
			  << "                              \"101\" when XExnSgn(2 downto 1) = \"00\"  else" << endl
			  << "                              \"100\" when XExnSgn(2 downto 1) = \"10\"  else" << endl
			  << "                              \"00\" & sR when (((Log_normal_normd(targetprec-1)='0') and (small='0')) or ( (Log_small_normd (wF+g-1)='0') and (small='1'))) or (ufl = '1') else" << endl
			  << "                               \"01\" & sR;" << endl;
		vhdl << tab << "R(wE+wF-1 downto 0) <=  EFR;" << endl;
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
		fpx.getMPFR(x); 
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
		mpz_class x;


		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addComment("1.0");
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);
		tc->addComment("The worst case of the error analysis: max cancellation, and full range reduction");
		x = (mpz_class(1) << wF) - (mpz_class(1) << (wF-pfinal+2)) // mantissa
			+ (((mpz_class(1) << (wE-1)) -2) << wF)  // exponent
			+ (mpz_class(1) << (wE+wF+1))	; // exn=010
		tc->addInput("X", x);
		emulate(tc);
		tcl->add(tc);



	}



	// One test out of 8 fully random (tests NaNs etc)
	// All the remaining ones test positive numbers.
	// with special treatment for exponents 0 and -1, 
	// and for the range reduction worst case.
 
	void FPLog::buildRandomTestCases(TestCaseList* tcl, int n){

		TestCase *tc;
		mpz_class a;
		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);

		for (int i = 0; i < n; i++) {
			tc = new TestCase(this); 
			/* Fill inputs */
			if ((i & 7) == 0)
				a = getLargeRandom(wE+wF+3);
			else if ((i & 7) == 1) // exponent of 1
				a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-1)) << wF) + normalExn; 
			else if ((i & 7) == 2) // exponent of 0.5
				a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-2)) << wF) + normalExn; 
			else if ((i & 7) == 3) { // worst case for range reduction
				tc->addComment("The worst case of the error analysis: max cancellation, and full range reduction");
				a = (mpz_class(1) << wF) - (mpz_class(1) << (wF-pfinal+2)) + getLargeRandom(wF-pfinal+2) // mantissa
					+ (((mpz_class(1) << (wE-1)) -2) << wF)  // exponent
					+ (mpz_class(1) << (wE+wF+1))	; // exn=010
			}
			else
				a  = getLargeRandom(wE+wF)  + normalExn; // 010xxxxxx
		
			tc->addInput("X", a);
			/* Get correct outputs */
			emulate(tc);
			// add to the test case list
			tcl->add(tc);
		}
	}
	TestCase* FPLog::buildRandomTestCases(int i){

		TestCase *tc;
		mpz_class a;
		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);

		tc = new TestCase(this); 
		/* Fill inputs */
		if ((i & 7) == 0)
			a = getLargeRandom(wE+wF+3);
		else if ((i & 7) == 1) // exponent of 1
			a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-1)) << wF) + normalExn; 
		else if ((i & 7) == 2) // exponent of 0.5
			a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-2)) << wF) + normalExn; 
		else if ((i & 7) == 3) { // worst case for range reduction
			tc->addComment("The worst case of the error analysis: max cancellation, and full range reduction");
			a = (mpz_class(1) << wF) - (mpz_class(1) << (wF-pfinal+2)) + getLargeRandom(wF-pfinal+2) // mantissa
				+ (((mpz_class(1) << (wE-1)) -2) << wF)  // exponent
				+ (mpz_class(1) << (wE+wF+1))	; // exn=010
		}
		else
			a  = getLargeRandom(wE+wF)  + normalExn; // 010xxxxxx
		
		tc->addInput("X", a);
		/* Get correct outputs */
		emulate(tc);
		// add to the test case list
		return tc;
	}

}
