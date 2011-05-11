/*
 * A post-normalization unit for the FloPoCo Long Accumulator 
 *
 * Author : Bogdan Pasca
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


#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "utils.hpp"
#include "Operator.hpp"
#include "LongAcc2FP.hpp"
#include "IntAdder.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	LongAcc2FP::LongAcc2FP(Target* target, int LSBA, int MSBA, int wEOut, int wFOut): 
		Operator(target), 
		LSBA_(LSBA), MSBA_(MSBA), wEOut_(wEOut), wFOut_(wFOut)
	{
	
		ownTarget_ = target;
		ostringstream name;
		setCopyrightString("Bogdan Pasca (2008-2009)");		
		name <<"LongAcc2FP_"
			  <<(LSBA_>=0?"":"M")<<abs(LSBA_)<<"_"
			  <<(MSBA_>=0?"":"M")<<abs(MSBA_)<<"_"
			  <<wEOut_<<"_"<<wFOut_;
		setName(name.str());
	
		sizeAcc_ = MSBA - LSBA + 1;
	
		//inputs and outputs
		addInput  ("A", sizeAcc_);
		addInput  ("C", sizeAcc_);
		addInput  ("AccOverflow",1);
		addOutput ("R", 3 + wEOut_ + wFOut_);

	
		if (target->isPipelined()) 
			setSequential();
		else
			setCombinatorial();
	
		// Set up various architectural parameters
		sizeAcc_ = MSBA_-LSBA_+1;	

		vhdl << tab <<declare("signA") << " <= A" << of(sizeAcc_-1)<<";"<<endl;
		vhdl << tab <<declare("signC") << " <= C" << of(sizeAcc_-1)<<";"<<endl;
		vhdl << tab <<declare("AccOverflowFlag") << " <= AccOverflow;"<<endl;	

		IntAdder *a = new IntAdder(target, sizeAcc_);
		oplist.push_back(a);
	
		inPortMap( a,   "X",   "A");
		inPortMap( a,   "Y",   "C");
		inPortMapCst(a, "Cin", "'0'");
		outPortMap( a,  "R",   "cpR");
		vhdl << instance(a,    "CarryPropagation");
	
		syncCycleFromSignal("cpR");
		setCriticalPath( a->getOutputDelay("R"));
		setSignalDelay( "cpR", getCriticalPath());

		vhdl << tab << declare("resSign") << " <= cpR" << of(sizeAcc_-1) << ";" << endl;

		manageCriticalPath(target->localWireDelay() + target->lutDelay());	
		//detect Addition overflow
		vhdl << tab << declare("signConcat",3) << " <= signA & signC & resSign;" << endl;
		vhdl << tab << "with signConcat select " << endl;
		vhdl << tab << declare("ovf") << " <= '0' when \"000\", "<<endl;
		vhdl << tab << "       '1' when \"001\","<<endl;
		vhdl << tab << "       '0' when \"010\","<<endl;
		vhdl << tab << "       '0' when \"011\","<<endl;
		vhdl << tab << "       '0' when \"100\","<<endl;
		vhdl << tab << "       '0' when \"101\","<<endl;
		vhdl << tab << "       '1' when \"110\","<<endl;
		vhdl << tab << "       '0' when \"111\","<<endl;
		vhdl << tab << "       '0' when others;"<<endl;

		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl << tab << declare("ovf_updated") << " <= ovf or AccOverflow;" << endl;
		
		//count the number of zeros/ones in order to determine the value of the exponent
		//Leading Zero/One counter 
		
		setCycleFromSignal("cpR");
		setCriticalPath( getSignalDelay("cpR") );
		
		lzocShifterSticky_ = new LZOCShifterSticky(target, sizeAcc_,  wFOut_ + 1, intlog2(sizeAcc_), false, -1, inDelayMap("I",target->localWireDelay()+getCriticalPath())); 
		oplist.push_back(lzocShifterSticky_);
		countWidth_ = lzocShifterSticky_->getCountWidth();
	
		inPortMap(lzocShifterSticky_, "I", "cpR");
		inPortMap(lzocShifterSticky_, "OZb", "resSign");
		outPortMap(lzocShifterSticky_, "Count", "nZO");
		outPortMap(lzocShifterSticky_, "O"    , "resFrac");
		vhdl << tab << instance(lzocShifterSticky_, "InputLZOCShifter");
		
		syncCycleFromSignal("resFrac");
		setCriticalPath( lzocShifterSticky_->getOutputDelay("O") );	
		setSignalDelay("resFrac",getCriticalPath());

//		// c2MnZO is  -nZO in 2'c complement with a twist. Usually for 2's complement all bits must be inverted
//		// and then a 1 must be added to the LSB. We postpone the extra 1 addition until later.
//		manageCriticalPath(target->localWireDelay()+target->lutDelay());
//		vhdl << tab <<declare("c2MnZO", countWidth_+1) << " <= \"1\" & not(nZO);"<<endl; //the bit inversion is done in O(1)
	
		// the 1 is added here, as a carry in bit for the substraction MSBA - nZO,
		// which is an addition in 2's complement
		manageCriticalPath(target->localWireDelay() + target->adderDelay(countWidth_+1));
		vhdl << tab <<declare("expAdj", countWidth_+1) << " <= CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<countWidth_+1<<") - (\"0\" & nZO);"<<endl;
		expBias_ =  intpow2(wEOut-1) - 1;
		if (countWidth_+1 < wEOut_){
			// this case is encountered most often.
			// by adding expBias to expAdj, we always get a positive quantity, so no need to converte back to SM
			manageCriticalPath(target->localWireDelay() + target->adderDelay(wEOut_));
			//the exponent bias	
			vhdl << tab << declare("expBias",wEOut_) << " <= CONV_STD_LOGIC_VECTOR("<<expBias_<<","<<wEOut_<<");"<<endl; //fixed value
			vhdl << tab <<declare("expRAdjusted",wEOut_) << " <= expBias + "<< // + sign extended expAdj
				"(" << rangeAssign(wEOut_-1,countWidth_+1, use("expAdj")+of(countWidth_)) << " & expAdj);"<<endl;
			vhdl << tab <<declare("excBits",2) << " <=\"01\";"<<endl;
		} else if (countWidth_+1 == wEOut_){
			manageCriticalPath(target->localWireDelay() + target->adderDelay(wEOut_));
			vhdl << tab << declare("expBias",wEOut_) << " <= CONV_STD_LOGIC_VECTOR("<<expBias_<<","<<wEOut_<<");"<<endl; //fixed value
			vhdl << tab <<declare("expRAdjusted",wEOut_) <<" <= expBias + expAdj;"<<endl;
			vhdl << tab <<declare("excBits",2) << " <=\"01\";"<<endl;
		}else{
			// case encountered when we wish to put the contents of the accumulator in a floating point number 
			// with a very small exponent
			manageCriticalPath(target->localWireDelay() + target->adderDelay(countWidth_+1));
			vhdl << tab << declare("expBias",wEOut_) << " <= CONV_STD_LOGIC_VECTOR("<<expBias_<<","<<wEOut_<<");"<<endl; //fixed value
			vhdl << tab <<declare("expAdjustedExt",countWidth_+1)<<" <= (" << rangeAssign(countWidth_, wEOut_, "'0'") << " & expBias) + expAdj;"<<endl;
	
			vhdl << tab <<declare("signExpAdjustedExt") << " <= expAdjustedExt" << of(countWidth_)<<";"<<endl;

			manageCriticalPath(target->localWireDelay() + target->adderDelay(countWidth_));//FIXME it might not be enough for correct inference
			vhdl << tab <<declare("modulusExpAdjusted",countWidth_) << "<= ("<<rangeAssign(countWidth_-1,0, "signExpAdjustedExt")<<" xor expAdjustedExt"<<range(countWidth_-1,0)<<") + signExpAdjustedExt"<<";"<<endl;
		
			vhdl << tab <<declare("maxExponent",countWidth_)<<" <= CONV_STD_LOGIC_VECTOR("<<(mpz_class(1)<<wEOut_)-1<<" , "<<countWidth_<<");"<<endl;
	
//			manageCriticalPath(target->localWireDelay() + target->lutDelay());
			vhdl << tab <<declare("expOverflow")<<" <= '1' when ((modulusExpAdjusted > maxExponent) and (signExpAdjustedExt='0')) else '0';"<<endl;
			vhdl << tab <<declare("expUnderflow")<<" <= '1' when ((modulusExpAdjusted > maxExponent) and (signExpAdjustedExt='1')) else '0';"<<endl;
	
			declare("excBits",2);
			vhdl << tab <<"excBits(1) <= expOverflow and  not(expUnderflow);"<<endl;
			vhdl << tab <<"excBits(0) <= not(expOverflow) and not(expUnderflow);"<<endl;
			
			vhdl << tab <<declare("expRAdjusted",wEOut_) << " <= expAdjustedExt" << range(wEOut_-1,0)<<";"<<endl;
		}			
	
//		manageCriticalPath(target->localWireDelay() + target->lutDelay());
		vhdl << tab <<declare("excRes",2) << " <= excBits when ovf_updated='0' else \"10\";"<<endl;

		//get back to the fraction part
		setCycleFromSignal("resFrac");
		setCriticalPath( getSignalDelay("resFrac"));
		
		manageCriticalPath( target->localWireDelay() + target->lutDelay()); 
		vhdl << tab << declare("notResFrac",wFOut_+1) << " <= resFrac xor "<<rangeAssign(wFOut_,0,"resSign")<<";"<<endl;

		//a possible carry propagation
		IntAdder *b = new IntAdder(target,wFOut_ + 1, inDelayMap("X", target->localWireDelay() + getCriticalPath()) );
		oplist.push_back(b);
	
		inPortMap(b, "X", "notResFrac");
		inPortMapCst(b, "Y", zg(wFOut+1,0));
		inPortMap(b, "Cin", "resSign");
		outPortMap(b, "R", "resultFraction");
		vhdl << tab << instance(b, "carryPropagator");
		syncCycleFromSignal("resultFraction");		
		setCriticalPath( b->getOutputDelay("R"));
	
		vhdl << tab << "R <= excRes & resSign & expRAdjusted" << range(wEOut_-1,0)<<" & resultFraction" << range(wFOut_-1,0) <<";"<<endl;
		outDelayMap["R"] = getCriticalPath();
	}

	LongAcc2FP::~LongAcc2FP() {
	}


	void LongAcc2FP::buildRandomTestCases(TestCaseList* tcl, int n){

		TestCase *tc;
		mpz_class A,C;
		
		C=mpz_class(0);
		int chunkSize, k;
		
		ownTarget_->suggestSubaddSize(chunkSize, MSBA_-LSBA_+1);
		if (chunkSize>=MSBA_-LSBA_+1)
			k=1;
		else
			k= int(ceil( double(MSBA_-LSBA_+1)/double(chunkSize)));

		for (int i = 0; i < n; i++) {
			C= mpz_class(0);
			tc = new TestCase(this);
			
			A = getLargeRandom(MSBA_-LSBA_+1);

			tc->addInput("A", A);
			
			if (k==1)
				tc->addInput("C", mpz_class(0));
			else{
				for (int j=0;j<k-1;j++){
					C = C + (  getLargeRandom(1)<< (chunkSize*(j+1)) );
					tc->addInput("C", C);
				}
			
			
			}
			tc->addInput("AccOverflow",mpz_class(0));
			
			/* Get correct outputs */
			emulate(tc);
			tcl->add(tc);
		}
	}
	TestCase* LongAcc2FP::buildRandomTestCases(int i){

		TestCase *tc;
		mpz_class A,C;
		
		C=mpz_class(0);
		int chunkSize, k;
		
		ownTarget_->suggestSubaddSize(chunkSize, MSBA_-LSBA_+1);
		if (chunkSize>=MSBA_-LSBA_+1)
			k=1;
		else
			k= int(ceil( double(MSBA_-LSBA_+1)/double(chunkSize)));

		C= mpz_class(0);
		tc = new TestCase(this);
		
		A = getLargeRandom(MSBA_-LSBA_+1);

		tc->addInput("A", A);
		
		if (k==1)
			tc->addInput("C", mpz_class(0));
		else{
			for (int j=0;j<k-1;j++){
				C = C + (  getLargeRandom(1)<< (chunkSize*(j+1)) );
				tc->addInput("C", C);
			}
		
		
		}
		tc->addInput("AccOverflow",mpz_class(0));
		
		/* Get correct outputs */
		emulate(tc);
		return tc;	
	}


	void LongAcc2FP::emulate(TestCase *tc)
	{
		/* Get I/O values */
		mpz_class svA           = tc->getInputValue("A");
		mpz_class svAccOverflow = tc->getInputValue("AccOverflow");
		mpz_class svC           = tc->getInputValue("C");
		
	
		mpz_class newAcc;
		newAcc = svA + svC;
				  	
		mpz_class tmpSUB = (mpz_class(1) << (MSBA_ - LSBA_+1))    ;
		mpz_class tmpCMP = (mpz_class(1) << (MSBA_ - LSBA_  )) - 1;

		if (newAcc > tmpCMP)
			newAcc = newAcc - tmpSUB;

		mpfr_t x;
		mpfr_init2(x, 10000); //init to infinite prec
		mpfr_set_z(x, newAcc.get_mpz_t(), GMP_RNDN);

		mpfr_t cst, tmp2;
		mpfr_init2(cst, 10000); //init to infinite prec
		mpfr_init2(tmp2, 10000); //init to infinite prec


		mpfr_set_ui(cst, 2 , GMP_RNDN);
		mpfr_set_si(tmp2, LSBA_ , GMP_RNDN);
		mpfr_pow(cst, cst, tmp2, GMP_RNDN);

		mpfr_mul(x, x, cst, GMP_RNDN);

		mpfr_t myFP;
		mpfr_init2(myFP, wFOut_+1);
		
		mpfr_set(myFP, x, GMP_RNDD);
		FPNumber  fpr(wEOut_, wFOut_, myFP);
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);

		mpfr_set(myFP, x, GMP_RNDU);
		FPNumber  fpr2(wEOut_, wFOut_, myFP);
		mpz_class svR2 = fpr2.getSignalValue();
		tc->addExpectedOutput("R", svR2);



		// clean-up
		mpfr_clears(x, myFP, NULL);
	}
	
}
