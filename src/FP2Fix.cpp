/*
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "FP2Fix.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <math.h>
#include <locale>

#include <stdio.h>
#include <mpfr.h>

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0


   FP2Fix::FP2Fix(Target* target, int _LSBO, int _MSBO, int _Signed,int _wEI, int _wFI, bool _trunc_p) :
         Operator(target), wEI(_wEI), wFI(_wFI), Signed(_Signed), LSBO(_LSBO), MSBO(_MSBO), trunc_p(_trunc_p) {

      int MSB=MSBO;
      int LSB=LSBO;
      
      
      ostringstream name;

      if ((MSB < LSB)){
         cerr << " FP2Fix: Input constraint LSB <= MSB not met."<<endl;
         exit (EXIT_FAILURE);
      }

      int wFO = MSB - LSB + 1;
      mpz_class maxExpWE = mpz_class(1)<<(wEI-1);
      mpz_class minExpWE = 1 - maxExpWE;

      int eMax = static_cast<int>(maxExpWE.get_si()) - 1;
      int wFO0;
      if(eMax+1 < MSB + 1)
         wFO0 = eMax + 1 - LSB;
      else
         wFO0 = MSB + 1 - LSB;

      if (( maxExpWE < MSB ) || ( minExpWE > LSB)){
         cerr << " The exponent is too small for full coverage. Try increasing the exponent !"<<endl;
         exit (EXIT_FAILURE);
      }
      int absMSB = MSB>=0?MSB:-MSB;
      int absLSB = LSB>=0?LSB:-LSB;
      name<<"FP2Fix_" << wEI << "_" << wFI << (LSB<0?"M":"") << "_" << absLSB << "_" << (MSB<0?"M":"") << absMSB <<"_"<< (Signed==1?"S":"US") << "_" << (trunc_p==1?"T":"NT");
      setName(name.str());

      setCopyrightString("Fabrizio Ferrandi (2011)");

      /* Set up the IO signals */

      addFPInput ("I", wEI,wFI);
      addOutput ("O", MSB-LSB+1);

      /*	VHDL code description	*/
      vhdl << tab << declare("neA0",wEI) << " <= not I" << range(wEI+wFI-1,wFI) << ";"<<endl;
      mpz_class tobeshiftedof;
      if(eMax < MSB)
         tobeshiftedof = eMax + eMax;
      else
         tobeshiftedof = eMax+MSB;
      vhdl << tab << declare("tobeshiftedof",wEI) << " <= conv_std_logic_vector(" << tobeshiftedof << ", "<< wEI<<");"<<endl;

      Exponent_difference = new IntAdder(target, wEI);
      Exponent_difference->changeName(getName()+"Exponent_difference");
      oplist.push_back(Exponent_difference);
      inPortMap  (Exponent_difference, "X", "tobeshiftedof");
      inPortMap  (Exponent_difference, "Y", "neA0");
      inPortMapCst(Exponent_difference, "Cin", "'1'");
      outPortMap (Exponent_difference, "R","eA1");
      vhdl << instance(Exponent_difference, "Exponent_difference");
      
      setCycleFromSignal("eA1");
      setCriticalPath(Exponent_difference->getOutputDelay("R"));

      if (wFI+1 < wFO0+2)
      {
         vhdl << tab << declare("fA0",wFO0+2) << " <= \"1\" & I" << range(wFI-1, 0) << " & " << rangeAssign(wFO0-wFI,0,"'0'")<<";"<<endl;
      }
      else
      {
         vhdl << tab << declare("fA0",wFO0+2) << range(wFO0+1,1) << " <= \"1\" & I" << range(wFI-1, wFI-wFO0) << ";" << endl;
         manageCriticalPath(target->localWireDelay() + target->lutDelay());
         vhdl << tab << "fA0" << of(0) << " <= '0' when I" << range(wFI-wFO0-1, 0) << " = " << rangeAssign(wFI-wFO0-1,0,"'0'") << " else '1';"<<endl;
      }

      //FXP shifter mappings
      FXP_shifter = new FXP_Shift(target, wEI, wFO0);
      oplist.push_back(FXP_shifter);

      inPortMap (FXP_shifter, "fA", "fA0");
      inPortMap (FXP_shifter, "n", "eA1");
      outPortMap (FXP_shifter, "fR", "fA1");
      vhdl << instance(FXP_shifter, "FXP_shifter");

      syncCycleFromSignal("fA1");
      setCriticalPath(FXP_shifter->getOutputDelay("R"));

      if(trunc_p)
      {
         vhdl << tab << declare("fA2",wFO+1) <<  "<= " << rangeAssign(wFO-wFO0,0,"'0'") << " & fA1" << range(wFO0+1, 2)<< ";"<<endl;
      }
      else
      {
         manageCriticalPath(target->localWireDelay() + target->lutDelay());
         vhdl << tab << declare("round") << " <= fA1(1) and (fA1(2) or fA1(0));"<<endl;

         vhdl << tab << declare("fA2a",wFO+1) <<  "<= " << rangeAssign(wFO-wFO0,0,"'0'") << " & fA1" << range(wFO0+1, 2)<< ";"<<endl;
         vhdl << tab << declare("fA2b",wFO+1) <<  "<= " << rangeAssign(wFO,1,"'0'") << " & round;"<<endl;
         MantSum = new IntAdder(target, wFO+1);
         MantSum->changeName(getName()+"MantSum");
         oplist.push_back(MantSum);
         inPortMap  (MantSum, "X", "fA2a");
         inPortMap  (MantSum, "Y", "fA2b");
         inPortMapCst(MantSum, "Cin", "'0'");
         outPortMap (MantSum, "R","fA2");
         vhdl << instance(MantSum, "MantSum");

         setCycleFromSignal("fA2");
         setCriticalPath(MantSum->getOutputDelay("R"));
      }
      if (eMax+1 > MSB+1)
      {
         manageCriticalPath(target->localWireDelay() + target->lutDelay());
         vhdl << tab << declare("overFl0") << "<= '1' when I" << range(wEI+wFI-1,wFI) << " > conv_std_logic_vector("<< eMax+MSB << "," << wEI << ") else I" << of(wEI+wFI+2)<<";"<<endl;
      }
      else
      {
         vhdl << tab << declare("overFl0") << "<= I" << of(wEI+wFI+2)<<";"<<endl;
      }
      
      manageCriticalPath(target->localWireDelay() + target->lutDelay());
      if(Signed == 0)
         vhdl << tab << declare("overFl1") << " <= fA2" << of(wFO) << ";"<<endl;
      else
         vhdl << tab << declare("overFl1") << " <= fA2" << of(wFO) << " or fA2" << of(wFO-1)<<";"<<endl;
      if (minExpWE < LSB)
         vhdl << tab << declare("undrFl0") << " <= '1' when I" << range(wEI+wFI-1,wFI) << " < conv_std_logic_vector(" << eMax+LSB-1 << "," << wEI <<
               ") else not (I" << of (wEI+wFI+2) << " or I" << of(wEI+wFI+1) << ");" << endl;
      else
         vhdl << tab << declare("undrFl0") << " <= not (I" << of (wEI+wFI+2) << " or I" << of(wEI+wFI+1) << ");" << endl;
      vhdl << tab << declare("fA3", wFO) << of(0) << " <= fA2" << of(0) << " xor I" << of(wEI+wFI) << ";" << endl;
      for(int i=1; i < wFO; ++i)
         vhdl << tab << "fA3" << of(i) << " <= fA2" << of(i) << " xor I" << of(wEI+wFI) << ";" << endl;
      vhdl << tab << declare("fA3b", wFO) << " <= " << rangeAssign(wFO-1,1,"'0'") << " & I(" << wEI+wFI << ");" << endl;

      MantSum2 = new IntAdder(target, wFO);
      MantSum2->changeName(getName()+"MantSum2");
      oplist.push_back(MantSum2);
      inPortMap  (MantSum2, "X", "fA3");
      inPortMap  (MantSum2, "Y", "fA3b");
      inPortMapCst(MantSum2, "Cin", "'0'");
      outPortMap (MantSum2, "R","fA4");
      vhdl << instance(MantSum2, "MantSum2");

      setCycleFromSignal("fA4");
      setCriticalPath(MantSum2->getOutputDelay("R"));

      vhdl << tab << declare("eTest",2) << " <= (overFl0 or overFl1) & undrFl0;" << endl;

      manageCriticalPath(target->localWireDelay() + target->lutDelay());
      vhdl << tab << "with eTest select" << endl;
      vhdl << tab << tab << "O <= " << rangeAssign(wFO-1,0,"'0'") << " when \"01\"," << endl;
      vhdl << tab << tab << "I" << of(wEI+wFI) << " & (" << wFO-2 << " downto 0 => not I" << of(wEI+wFI) << ") when \"10\","<<endl;
      vhdl << tab << tab << "fA4 when others;"<<endl;
   }


   FP2Fix::~FP2Fix() {
   }


   void FP2Fix::emulate(TestCase * tc)
   {
      /* Get I/O values */
      mpz_class svI = tc->getInputValue("I");
      FPNumber  fpi(wEI, wFI, svI);
      mpfr_t i;
      mpfr_init2(i, 1+wFI);
      fpi.getMPFR(i);
      std::cerr << "FP " << printMPFR(i, 100) << std::endl;
      mpz_class svO;

      if(trunc_p)
         mpfr_get_z(svO.get_mpz_t(), i, GMP_RNDZ);
      else
         mpfr_get_z(svO.get_mpz_t(), i, GMP_RNDN);
      svO = svO >> LSBO;
      if (Signed != 0) {
         mpz_class tmpCMP = (mpz_class(1)  << (MSBO-LSBO))-1;
         if (svO > tmpCMP){ //negative number 
            mpz_class tmpSUB = (mpz_class(1) << (MSBO-LSBO+1));
            svO = svO - tmpSUB;
	 }
      }
      svO = svO & ((mpz_class(1) << (MSBO-LSBO+1))-1);
      std::cerr << "FIX " << svO << std::endl;
      tc->addExpectedOutput("O", svO);
      // clean-up
      mpfr_clear(i);
  }
  
  TestCase* FP2Fix::buildRandomTestCase(int i)
  {
     TestCase *tc;
     mpz_class a;
     tc = new TestCase(this); 
     mpz_class e = (getLargeRandom(wEI+wFI) % (MSBO-LSBO+1)); // Should be between 0 and MSBO-LSBO+1
     mpz_class normalExn = mpz_class(1)<<(wEI+wFI+1);
     mpz_class bias = ((1<<(wEI-1))-1);
     e = bias + e;
     a = getLargeRandom(wFI) + (e << wFI) + normalExn;
     tc->addInput("I", a);
     /* Get correct outputs */
     emulate(tc);
     return tc;		
   }
   
   void FP2Fix::buildRandomTestCases(TestCaseList* tcl, int n)
   {
      TestCase *tc;
      mpz_class a;
      for (int i = 0; i < n; i++)
      {
         tc = new TestCase(this); 
         /* Fill inputs */
	 mpz_class e = (getLargeRandom(wEI+wFI) % (MSBO-LSBO+1)); // Should be between 0 and MSBO-LSBO+1
	 mpz_class normalExn = mpz_class(1)<<(wEI+wFI+1);
	 mpz_class bias = ((1<<(wEI-1))-1);
	 e = bias + e;
	 a = getLargeRandom(wFI) + (e << wFI) + normalExn;
         tc->addInput("I", a);
         /* Get correct outputs */
         emulate(tc);
         tcl->add(tc);
      }
   }


   void FP2Fix::buildStandardTestCases(TestCaseList* tcl){

   }

}
