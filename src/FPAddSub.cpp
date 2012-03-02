
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "../utils.hpp"

#include "FPAddSub.hpp"

using namespace std;

//TODO +- inf for exponent => update exception

namespace flopoco{

	extern vector<Operator*> oplist;

#define DEBUGVHDL 0


FPAddSub::FPAddSub(Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR, map<string, double> inputDelays) :
		Operator(target), wEX(wEX), wFX(wFX), wEY(wEY), wFY(wFY), wER(wER), wFR(wFR) {

		srcFileName="FPAddSub";
			
		//parameter set up. For now all wEX=wEY=wER and the same holds for fractions
		wF = wFX;
		wE = wEX;
		
		//set operator name
		ostringstream name;
		name<<"FPAddSub_"<<wE<<"_"<<wF<<"_uid"<<getNewUId(); 
		setName(name.str()); 

		//size of signal holding shift amount
		sizeRightShift = intlog2(wF+3);

		/* Set up the IO signals */
		/* Inputs: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction) */
		vhdl << "-- Set up the IO signals" << endl;
		vhdl << "-- The FloPoCo floating point format: 2b(Exception) + 1b(Sign) + wEX bits (Exponent) + wFX bits(Fraction)" << endl;
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		addFPInput ("X", wE, wF);
		addFPInput ("Y", wE, wF);
		addFPOutput("Radd", wE, wF);
		addFPOutput("Rsub", wE, wF);

		//=========================================================================|
		//                          Swap/Difference                                |
		// ========================================================================|
		vhdl << "-- Exponent difference and swap  --" << endl;

		vhdl << "-- Determine which of the two inputs is larger an swap so that X>Y" << endl;
		vhdl << tab << declare("excExpFracX",2+wE+wF) << " <= X"<<range(wE+wF+2, wE+wF+1) << " & X"<<range(wE+wF-1, 0)<<";"<<endl;
		vhdl << tab << declare("excExpFracY",2+wE+wF) << " <= Y"<<range(wE+wF+2, wE+wF+1) << " & Y"<<range(wE+wF-1, 0)<<";"<<endl;

		
		setCriticalPath(getMaxInputDelays(inputDelays));
		manageCriticalPath(target->localWireDelay() + target->adderDelay(wE+1));
		vhdl << "-- Compute exponent difference" << endl;
		vhdl<< tab << declare("eXmeY",wE+1) << " <= (\"0\" & X"<<range(wE+wF-1,wF)<<") - (\"0\" & Y"<<range(wE+wF-1,wF)<<");"<<endl;
		vhdl<< tab << declare("eYmeX",wE+1) << " <= (\"0\" & Y"<<range(wE+wF-1,wF)<<") - (\"0\" & X"<<range(wE+wF-1,wF)<<");"<<endl;
		double cpeXmeY = getCriticalPath();
		
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		vhdl << "-- Flag for swap" << endl;
		if (wF < 30){
			manageCriticalPath(target->localWireDelay() + target->adderDelay(wE)); //comparator delay implemented for now as adder
			vhdl << tab << declare("swap") << " <= '0' when excExpFracX >= excExpFracY else '1';" << endl;
		}else{
			IntAdder *cmpAdder = new IntAdder(target, wE+wF+2+1);
			oplist.push_back(cmpAdder);
			
			vhdl << tab << declare("addCmpOp1",wE+wF+2+1) << "<= " << zg(1,0) << " & excExpFracX;"<<endl;
			vhdl << tab << declare("addCmpOp2",wE+wF+2+1) << "<= " << og(1,0) << " & not(excExpFracY);"<<endl;
			
			inPortMap(cmpAdder, "X", "addCmpOp1");
			inPortMap(cmpAdder, "Y", "addCmpOp2");
			inPortMapCst(cmpAdder, "Cin", "'1'");
			outPortMap (cmpAdder, "R", "cmpRes");
			
			vhdl << instance(cmpAdder, "cmpAdder") << endl;
			syncCycleFromSignal("cmpRes");
			setCriticalPath( cmpAdder->getOutputDelay("R") );
			vhdl<< tab << declare("swap") << " <= cmpRes"<<of(wE+wF+2)<<";"<<endl;
		}
		
		double cpswap = getCriticalPath();
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		
		vhdl << "-- Depending on the value of the swap flag, assign the corresponding values to the newX and newY signals" << endl;
		// Depending on the value of swap, assign the corresponding values to the newX and newY signals
		 
		vhdl<<tab<<declare("newX",wE+wF+3) << " <= X when swap = '0' else Y;"<<endl;
		vhdl<<tab<<declare("newY",wE+wF+3) << " <= Y when swap = '0' else X;"<<endl;
		
		//Break down the signals
		vhdl << "-- Break down the signals" << endl;
		vhdl << tab << declare("expX",wE) << "<= newX"<<range(wE+wF-1,wF)<<";"<<endl;
		vhdl << tab << declare("excX",2)  << "<= newX"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("excY",2)  << "<= newY"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("signX")   << "<= newX"<<of(wE+wF)<<";"<<endl;
		vhdl << tab << declare("signY")   << "<= newY"<<of(wE+wF)<<";"<<endl;
		vhdl << tab << declare("EffSub")  << " <= signX xor signY;"<<endl;
		vhdl << tab << declare("sdsXsYExnXY",6) << " <= signX & signY & excX & excY;"<<endl; 
		vhdl << tab << declare("sdExnXY",4) << " <= excX & excY;"<<endl; 
		manageCriticalPath(target->localWireDelay()+ target->lutDelay());
		vhdl << tab << declare("fracY",wF+1) << " <= "<< zg(wF+1)<<" when excY=\"00\" else ('1' & newY("<<wF-1<<" downto 0));"<<endl;
		double cpfracY = getCriticalPath();

		
		
		//exception bits: need to be updated but for not FIXME
		vhdl << "-- Set exception bits" << endl;
		manageCriticalPath(target->localWireDelay()+2*target->lutDelay());
		vhdl <<tab<<"with sdsXsYExnXY select "<<endl;
		vhdl <<tab<<declare("excRt",2) << " <= \"00\" when \"000000\"|\"010000\"|\"100000\"|\"110000\","<<endl
		<<tab<<tab<<"\"01\" when \"000101\"|\"010101\"|\"100101\"|\"110101\"|\"000100\"|\"010100\"|\"100100\"|\"110100\"|\"000001\"|\"010001\"|\"100001\"|\"110001\","<<endl
		<<tab<<tab<<"\"10\" when \"111010\"|\"001010\"|\"001000\"|\"011000\"|\"101000\"|\"111000\"|\"000010\"|\"010010\"|\"100010\"|\"110010\"|\"001001\"|\"011001\"|\"101001\"|\"111001\"|\"000110\"|\"010110\"|\"100110\"|\"110110\", "<<endl
		<<tab<<tab<<"\"11\" when others;"<<endl;
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl <<tab<<declare("signRadd") << "<= '0' when (sdsXsYExnXY=\"100000\" or sdsXsYExnXY=\"010000\") else signX;"<<endl;
		vhdl <<tab<<declare("signRsub") << "<= '0' when (sdsXsYExnXY=\"100000\" or sdsXsYExnXY=\"010000\") else (((signX and not swap) or (signX and swap and not signY)) or (signY and swap));"<<endl;
		
		
		setCycleFromSignal("swap");;
		if ( getCycleFromSignal("eYmeX") == getCycleFromSignal("swap") )
			setCriticalPath(max(cpeXmeY, cpswap));
		else{
			if (syncCycleFromSignal("eYmeX"))
				setCriticalPath(cpeXmeY);
		}
		manageCriticalPath(target->localWireDelay() + target->lutDelay());		// Multiplexer
		vhdl << tab << declare("expDiff",wE+1) << " <= eXmeY when swap = '0' else eYmeX;" << endl; 
		manageCriticalPath(target->localWireDelay() + target->eqConstComparatorDelay(wE+1));
		vhdl << tab << declare("shiftedOut") << " <= '1' when (expDiff >= " << wF+2 << ") else '0';" << endl;
		
		
		// ShiftVal=the number of positions that fracY must be shifted to the right
		vhdl << "-- Set the number of positions that the fractionary part of Y must be shifted to the right" << endl;
		if (wE>sizeRightShift) {
			manageCriticalPath(target->localWireDelay() + target->lutDelay());
			vhdl << tab << declare("shiftVal",sizeRightShift) << " <= expDiff(" << sizeRightShift-1 << " downto 0)"
				<< " when shiftedOut='0' else CONV_STD_LOGIC_VECTOR(" <<wFX+3 << "," <<sizeRightShift << ") ;" << endl; 
		}		
		else if (wE==sizeRightShift) {
			vhdl << tab << declare("shiftVal",sizeRightShift) << " <= expDiff;" << endl ;
		}
		else 	{ //  wE< sizeRightShift
			vhdl << tab << declare("shiftVal",sizeRightShift) << " <= CONV_STD_LOGIC_VECTOR(0," << sizeRightShift-wE <<") & expDiff;" << endl;
		}
		
		if ( getCycleFromSignal("fracY") == getCycleFromSignal("shiftVal") )
			setCriticalPath( max(cpfracY, getCriticalPath()) );
		else{
			if (syncCycleFromSignal("fracY"))
				setCriticalPath(cpfracY);
		}		
		
		
		// Shift right the significand of new Y with as many positions as the exponent difference suggests (alignment) 
		vhdl << "-- Shift right the significand of new Y with as many positions as the exponent difference suggests (alignment)" << endl;
		rightShifter = new Shifter(target,wF+1,wF+3, Shifter::Right, inDelayMap("X",getCriticalPath()));
		rightShifter->changeName(getName()+"_RightShifter");
		oplist.push_back(rightShifter);
		inPortMap  (rightShifter, "X", "fracY");
		inPortMap  (rightShifter, "S", "shiftVal");
		outPortMap (rightShifter, "R","shiftedFracY");
		vhdl << instance(rightShifter, "RightShifterComponent");
		syncCycleFromSignal("shiftedFracY");
		setCriticalPath(rightShifter->getOutputDelay("R"));
		nextCycle();
		setCriticalPath(0.0);
		double cpshiftedFracY = getCriticalPath();
		
		
		//sticky compuation in parallel with addition, no need for manageCriticalPath
		//FIXME: compute inside shifter
		// Compute sticky bit as the OR of the shifted out bits during the alignment
		manageCriticalPath(target->localWireDelay() + target->eqConstComparatorDelay(wF+1));
		vhdl << "-- Compute sticky bit as the or of the shifted out bits during the alignment" << endl;
		vhdl<<tab<< declare("sticky") << " <= '0' when (shiftedFracY("<<wF<<" downto 0)=CONV_STD_LOGIC_VECTOR(0,"<<wF<<")) else '1';"<<endl;
		double cpsticky = getCriticalPath();
		
		
		setCycleFromSignal("shiftedFracY");
		nextCycle();         
		setCriticalPath(0.0);
		setCriticalPath(cpshiftedFracY);
		// Pad fraction of Y [overflow][shifted frac having inplicit 1][guard][round]
		vhdl << "-- Pad fraction of Y [overflow][shifted frac having inplicit 1][guard][round]" << endl;
		vhdl << tab << declare("fracYfar", wF+4) << " <= \"0\" & shiftedFracY(" << 2*wF+3 << " downto " << wF+1 << ");" << endl;	
		vhdl << tab << declare("nfracYfar", wF+4) << " <= (\"0\" & shiftedFracY(" << 2*wF+3 << " downto " << wF+1 << ")) xor " << og(wF+4,0) << ";" << endl;
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		// Pad fraction of X [overflow][inplicit 1][fracX][guard bits]
		vhdl << "-- Pad fraction of X [overflow][inplicit 1][fracX][guard bits]" << endl;
		vhdl << tab << declare("fracXfar", wF+4) << " <= \"01\" & (newX(" << wF-1 << " downto 0)) & \"00\";" << endl;
		
		if (getCycleFromSignal("sticky")==getCycleFromSignal("fracXfar"))
			setCriticalPath( max (cpsticky, getCriticalPath()) );
		else
			if (syncCycleFromSignal("sticky"))
				setCriticalPath(cpsticky);
		manageCriticalPath(target->localWireDelay()+ target->lutDelay());
		
		// Result is always positive (the numbers are ordered in descending order, so even the subtraction gives a positive result)
		vhdl << "-- Perform the addition of the two fractional parts" << endl;
		fracAddFar = new IntAdder(target,wF+4, inDelayMap("X", getCriticalPath()));
		oplist.push_back(fracAddFar);
		inPortMap 	 (fracAddFar,   "X", "fracXfar"	  	  );
		inPortMap 	 (fracAddFar,   "Y", "fracYfar" 	  );
		inPortMapCst (fracAddFar,   "Cin", "\'0\'"	  	  );
		outPortMap	 (fracAddFar,   "R", "fracAddResult"  );
		vhdl << instance(fracAddFar, "fracAdder");
		
		vhdl << "-- Perform the subtraction of the two fractional parts" << endl;
		inPortMap  	 (fracAddFar,   "X", 	"fracXfar"		);
		inPortMap  	 (fracAddFar,   "Y", 	"nfracYfar"		);
		inPortMapCst (fracAddFar, "Cin", 	"\'1\'"			);
		outPortMap 	 (fracAddFar,   "R",	"fracSubResult"	);
		vhdl << instance(fracAddFar, "fracSub");
		
		syncCycleFromSignal("fracAddResult");
		syncCycleFromSignal("fracSubResult");
		setCriticalPath(fracAddFar->getOutputDelay("R"));
		
		if (getCycleFromSignal("sticky")==getCycleFromSignal("fracAddResult"))
			setCriticalPath(max(cpsticky, getCriticalPath()));
		else{
			if (syncCycleFromSignal("sticky"))
				setCriticalPath(cpsticky);
		}	
		
		// Shift in place
		vhdl << "-- Shift the significands in place" << endl;
		vhdl << tab << declare("fracGRS",wF+5) << "<= fracSubResult & sticky;" << endl;			//mantisa for subtraction
		vhdl << tab << declare("fracGRSAlt",wF+5) << "<= fracAddResult & \'0\' when (fracAddResult(" << wF+3 << ") = \'0\') else \'0\' & fracAddResult;"<<endl;		//mantisa for addition
	
		// Incremente the exponent
		vhdl << "-- Incremente the exponent" << endl;
		vhdl << tab << declare("extendedExpInc",wE+2) << "<= (\"00\" & expX) + '1';"<<endl;		//exponent for subtraction
		vhdl << tab << declare("extendedExp",wE+2) << "<= (\"00\" & expX);"<<endl;				//exponent for addition
		
		
		// Count leading zeros for the result of the subtraction
		vhdl << "-- Count leading zeros for the result of the subtraction" << endl;
		lzocs = new LZOCShifterSticky(target, wF+5, wF+5, intlog2(wF+5), false, 0, inDelayMap("I",getCriticalPath()));
		oplist.push_back(lzocs);
		inPortMap  (lzocs, "I", "fracGRS");
		outPortMap (lzocs, "Count","nZerosNew");
		outPortMap (lzocs, "O","shiftedFrac");
		vhdl << instance(lzocs, "LZC_component");
		syncCycleFromSignal("shiftedFrac");
		setCriticalPath(lzocs->getOutputDelay("O"));
		double cpshiftedFrac = getCriticalPath();		
		
		
		// Need to decide how much to add to the exponent
		// Update exponent
		vhdl << "-- Decide how much to add to the exponent and update the exponent" << endl;
		
		manageCriticalPath(target->localWireDelay() + target->adderDelay(wE+2));
		vhdl << tab << declare("updatedExpSub",wE+2) << " <= extendedExpInc - (" << zg(wE+2-lzocs->getCountWidth(),0) <<" & nZerosNew);"<<endl;
		vhdl << tab << declare("extendedExpAdd",wE+2) << "<= extendedExpInc when (fracAddResult(" << wF+3 << ") = \'1\') else extendedExp;"<<endl;
		vhdl << tab << declare("eqdiffsign")<< " <= '1' when nZerosNew="<<og(lzocs->getCountWidth(),0)<<" else '0';"<<endl; 
		
		
		// Concatenate exponent with fraction to absorb the possible carry out
		vhdl << "-- Concatenate exponent with fraction to absorb the possible carry out" << endl;
		vhdl<<tab<<declare("expFracSub",wE+2+wF+1)<<"<= updatedExpSub & shiftedFrac"<<range(wF+3,3)<<";"<<endl;
		vhdl<<tab<<declare("expFracAdd",wE+2+wF+1)<<"<= extendedExpAdd & fracGRSAlt"<<range(wF+3,3)<<";"<<endl;
		double cpexpFrac = getCriticalPath();
		
		
		// At least in parallel with previous 2 statements
		setCycleFromSignal("shiftedFrac");
		setCriticalPath(cpshiftedFrac);
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl << "-- Create signals for rounding - for subtraction" << endl;
		vhdl<<tab<<declare("stkSub")<<"<= shiftedFrac"<<of(1)<<" or shiftedFrac"<<of(0)<<";"<<endl;
		vhdl<<tab<<declare("rndSub")<<"<= shiftedFrac"<<of(2)<<";"<<endl;
		vhdl<<tab<<declare("grdSub")<<"<= shiftedFrac"<<of(3)<<";"<<endl;
		vhdl<<tab<<declare("lsbSub")<<"<= shiftedFrac"<<of(4)<<";"<<endl;
		
		vhdl << "-- Create signals for rounding - for addition" << endl;
		vhdl<<tab<<declare("stkAdd")<<"<= fracGRSAlt"<<of(1)<<" or fracGRSAlt"<<of(0)<<";"<<endl;
		vhdl<<tab<<declare("rndAdd")<<"<= fracGRSAlt"<<of(2)<<";"<<endl;
		vhdl<<tab<<declare("grdAdd")<<"<= fracGRSAlt"<<of(3)<<";"<<endl;
		vhdl<<tab<<declare("lsbAdd")<<"<= fracGRSAlt"<<of(4)<<";"<<endl;
		
		// Decide what to add to the guard bit
		vhdl << "-- Decide what to add to the guard bit" << endl;
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl<<tab<<declare("addToRoundBitSub")<<"<= '0' when (lsbSub='0' and grdSub='1' and rndSub='0' and stkSub='0')  else '1';"<<endl;
		vhdl<<tab<<declare("addToRoundBitAdd")<<"<= '0' when (lsbAdd='0' and grdAdd='1' and rndAdd='0' and stkAdd='0')  else '1';"<<endl;
		// Round
		vhdl << "-- Perform rounding" << endl;
		
		if (getCycleFromSignal("expFracAdd") == getCycleFromSignal("addToRoundBitAdd"))
			setCriticalPath(max(cpexpFrac, getCriticalPath()));
		else
			if (syncCycleFromSignal("expFracAdd"))
				setCriticalPath(cpexpFrac);
			
		IntAdder *raSub = new IntAdder(target, wE+2+wF+1, inDelayMap("X", getCriticalPath() ) );
		oplist.push_back(raSub);
		
		inPortMap( 	  raSub, "X",   "expFracSub");
		inPortMapCst( raSub, "Y",   zg(wE+2+wF+1,0) );
		inPortMap( 	  raSub, "Cin", "addToRoundBitSub");
		outPortMap(   raSub,  "R",  "RoundedExpFracSub");
		vhdl << instance(raSub, "roundingAdderSub");
		setCycleFromSignal("RoundedExpFracSub");
		
		inPortMap( 	  raSub, "X",   "expFracAdd");
		inPortMapCst( raSub, "Y",   zg(wE+2+wF+1,0) );
		inPortMap( 	  raSub, "Cin", "addToRoundBitAdd");
		outPortMap(   raSub,  "R",  "RoundedExpFracAdd");
		vhdl << instance(raSub, "roundingAdderAdd");
		setCycleFromSignal("RoundedExpFracAdd");
		
		setCriticalPath(raSub->getOutputDelay("R"));

		// Possible update to exception bits
		vhdl << "-- Possible update to exception bits" << endl;
		vhdl << tab << declare("upExcAdd",2)<<" <= RoundedExpFracAdd"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		vhdl << tab << declare("upExcSub",2)<<" <= RoundedExpFracSub"<<range(wE+wF+2,wE+wF+1)<<";"<<endl;
		
		vhdl << "-- Recreate addition result" << endl;
		vhdl << tab << declare("fracRAdd",wF)<<" <= RoundedExpFracAdd"<<range(wF,1)<<";"<<endl;
		vhdl << tab << declare("expRAdd",wE) <<" <= RoundedExpFracAdd"<<range(wF+wE,wF+1)<<";"<<endl;
		
		vhdl << "-- Recreate subtraction result" << endl;
		vhdl << tab << declare("fracRSub",wF)<<" <= RoundedExpFracSub"<<range(wF,1)<<";"<<endl;
		vhdl << tab << declare("expRSub",wE) <<" <= RoundedExpFracSub"<<range(wF+wE,wF+1)<<";"<<endl;

		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl << tab << declare("exExpExcAdd",4) << " <= upExcAdd & excRt;"<<endl;
		vhdl << tab << declare("exExpExcSub",4) << " <= upExcSub & excRt;"<<endl;
		
		vhdl << tab << "with (exExpExcAdd) select "<<endl;
		vhdl << tab << declare("excRt2Add",2) << "<= \"00\" when \"0000\"|\"0100\"|\"1000\"|\"1100\"|\"1001\"|\"1101\","<<endl
		<<tab<<tab<<"\"01\" when \"0001\","<<endl
		<<tab<<tab<<"\"10\" when \"0010\"|\"0110\"|\"0101\","<<endl
		<<tab<<tab<<"\"11\" when others;"<<endl;
		
		vhdl << tab << "with (exExpExcSub) select "<<endl;
		vhdl << tab << declare("excRt2Sub",2) << "<= \"00\" when \"0000\"|\"0100\"|\"1000\"|\"1100\"|\"1001\"|\"1101\","<<endl
		<<tab<<tab<<"\"01\" when \"0001\","<<endl
		<<tab<<tab<<"\"10\" when \"0010\"|\"0110\"|\"0101\","<<endl
		<<tab<<tab<<"\"11\" when others;"<<endl;
		
		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl<<tab<<declare("excRAdd",2) << " <= \"00\" when (eqdiffsign='1' and EffSub='1') else excRt2Add;"<<endl;
		vhdl<<tab<<declare("excRSub",2) << " <= \"00\" when (eqdiffsign='1' and EffSub='1') else excRt2Sub;"<<endl;
		

		// Assign result 
		vhdl << "-- Assign result of addition and of subtraction" << endl;
		vhdl << tab << declare("computedRAdd",wE+wF+3) << " <= excRAdd & ((signRadd and not EffSub) or (signRsub and EffSub)) & expRAdd & fracRAdd;" << endl;
		vhdl << tab << declare("computedRSub",wE+wF+3) << " <= excRSub & ((signRadd and EffSub) or (signRsub and not EffSub)) & expRSub & fracRSub;" << endl;
		vhdl << tab << "Radd <= computedRAdd when (EffSub = \'0\') else computedRSub;"<<endl;
		vhdl << tab << "Rsub <= computedRSub when (EffSub = \'1\') else computedRAdd;"<<endl;
	}

	FPAddSub::~FPAddSub() {
	}


	void FPAddSub::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
	
		/* Compute correct value */
		FPNumber fpx(wEX, wFX), fpy(wEY, wFY);
		fpx = svX;
		fpy = svY;
		mpfr_t x, y, radd, rsub;
		mpfr_init2(x, 1+wFX);
		mpfr_init2(y, 1+wFY);
		mpfr_init2(radd, 1+wFR);
		mpfr_init2(rsub, 1+wFR);
		fpx.getMPFR(x);
		fpy.getMPFR(y);
		mpfr_add(radd, x, y, GMP_RNDN);
		mpfr_sub(rsub, x, y, GMP_RNDN);

		// Set outputs 
		FPNumber  fpradd(wER, wFR, radd);
		mpz_class svRadd = fpradd.getSignalValue();
		tc->addExpectedOutput("Radd", svRadd);
		
		FPNumber  fprsub(wER, wFR, rsub);
		mpz_class svRsub = fprsub.getSignalValue();
		tc->addExpectedOutput("Rsub", svRsub);

		// clean up
		mpfr_clears(x, y, radd, rsub, NULL);
	}





	void FPAddSub::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		// Regression tests 
		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", -1.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", FPNumber::plusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", FPNumber::minusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::plusInfty);
		tc->addFPInput("Y", FPNumber::minusInfty);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::plusInfty);
		tc->addFPInput("Y", FPNumber::plusInfty);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addFPInput("X", FPNumber::minusInfty);
		tc->addFPInput("Y", FPNumber::minusInfty);
		emulate(tc);
		tcl->add(tc);
	
	}



	TestCase* FPAddSub::buildRandomTestCase(int i){

		TestCase *tc;
		mpz_class x,y;
		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);
		mpz_class negative  = mpz_class(1)<<(wE+wF);

		tc = new TestCase(this); 
		/* Fill inputs */
		if ((i & 7) == 0) {// cancellation, same exponent
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
		}
		else if ((i & 7) == 1) {// cancellation, exp diff=1
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			e++; // may rarely lead to an overflow, who cares
			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
		}
		else if ((i & 7) == 2) {// cancellation, exp diff=1
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
			e++; // may rarely lead to an overflow, who cares
			y  = getLargeRandom(wF) + (e << wF) + normalExn;
		}
		else if ((i & 7) == 3) {// alignment within the mantissa sizes
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow, who cares
			y  = getLargeRandom(wF) + (e << wF) + normalExn;
		}
		else if ((i & 7) == 4) {// subtraction, alignment within the mantissa sizes
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow
			y  = getLargeRandom(wF) + (e << wF) + normalExn + negative;
		}
		else if ((i & 7) == 5 || (i & 7) == 6) {// addition, alignment within the mantissa sizes
			mpz_class e = getLargeRandom(wE);
			x  = getLargeRandom(wF) + (e << wF) + normalExn;
			e +=	getLargeRandom(intlog2(wF)); // may lead to an overflow
			y  = getLargeRandom(wF) + (e << wF) + normalExn;
		}
		else{ //fully random
			x = getLargeRandom(wE+wF+3);
			y = getLargeRandom(wE+wF+3);
		}
		// Random swap
		mpz_class swap = getLargeRandom(1);
		if (swap == mpz_class(0)) {
			tc->addInput("X", x);
			tc->addInput("Y", y);
		}
		else {
			tc->addInput("X", y);
			tc->addInput("Y", x);
		}
		/* Get correct outputs */
		emulate(tc);
		return tc;
	}

}
