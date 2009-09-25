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

		vhdl << tab <<declare("signA") << " <= " << use("A")<<of(sizeAcc_-1)<<";"<<endl;
		vhdl << tab <<declare("signC") << " <= " << use("C")<<of(sizeAcc_-1)<<";"<<endl;
		vhdl << tab <<declare("AccOverflowFlag") << " <= " << use("AccOverflow")<<";"<<endl;	

		IntAdder *a = new IntAdder(target, sizeAcc_);
		oplist.push_back(a);
	
		inPortMap( a,   "X",   "A");
		inPortMap( a,   "Y",   "C");
		inPortMapCst(a, "Cin", "'0'");
		outPortMap( a,  "R",   "cpR");
		vhdl << instance(a,    "CarryPropagation");
	
		syncCycleFromSignal("cpR");
		nextCycle();
		vhdl << tab << declare("resSign") << " <= " << use("cpR")<<of(sizeAcc_-1) << ";" << endl;
	
		//detect Addition overflow
		vhdl << tab << declare("signConcat",3) << " <= " << use("signA") << " & " << use("signC") << " & " << use("resSign") << ";" << endl;
		vhdl << tab << "with " << use("signConcat") << " select " << endl;
		vhdl << tab << declare("ovf") << " <= '0' when \"000\", "<<endl;
		vhdl << tab << "       '1' when \"001\","<<endl;
		vhdl << tab << "       '0' when \"010\","<<endl;
		vhdl << tab << "       '0' when \"011\","<<endl;
		vhdl << tab << "       '0' when \"100\","<<endl;
		vhdl << tab << "       '0' when \"101\","<<endl;
		vhdl << tab << "       '1' when \"110\","<<endl;
		vhdl << tab << "       '0' when \"111\","<<endl;
		vhdl << tab << "       '0' when others;"<<endl;
	
		vhdl << tab << declare("ovf_updated") << " <= " << use("ovf") << " or " << use("AccOverflow") << ";" << endl;
		
		//count the number of zeros/ones in order to determine the value of the exponent
		//Leading Zero/One counter 
	
		lzocShifterSticky_ = new LZOCShifterSticky(target, sizeAcc_,  wFOut_ + 1, intlog2(sizeAcc_), false, -1); // target, inputNbOfBits, outputNbOfBits, wCount, computeSticky, countType 
		oplist.push_back(lzocShifterSticky_);
		countWidth_ = lzocShifterSticky_->getCountWidth();
	
		inPortMap(lzocShifterSticky_, "I", "cpR");
		inPortMap(lzocShifterSticky_, "OZb", "resSign");
		outPortMap(lzocShifterSticky_, "Count", "nZO");
		outPortMap(lzocShifterSticky_, "O"    , "resFrac");
		vhdl << tab << instance(lzocShifterSticky_, "InputLZOCShifter");
		
		syncCycleFromSignal("resFrac");	
		//	nextCycle();
		//the exponent bias	
		vhdl << tab << declare("expBias",wEOut_) << " <= CONV_STD_LOGIC_VECTOR("<<expBias_<<","<<wEOut_<<");"<<endl; //fixed value

		// c2MnZO is  -nZO in 2'c complement with a twist. Usually for 2's complement all bits must be inverted
		// and then a 1 must be added to the LSB. We postpone the extra 1 addition until later.
		vhdl << tab <<declare("c2MnZO", countWidth_+1) << " <= \"1\" & not(nZO);"<<endl; //the bit inversion is done in O(1)
	
		// the 1 is added here, as a carry in bit for the substraction MSBA - nZO,
		// which is an addition in 2's complement
		vhdl << tab <<declare("expAdj", countWidth_+1) << " <= CONV_STD_LOGIC_VECTOR("<<MSBA_<<","<<countWidth_+1<<") + c2MnZO + 1;"<<endl;
	
		if (countWidth_+1 < wEOut_){
			// this case is encountered most often.
			// by adding expBias to expAdj, we always get a positive quantity, so no need to converte back to SM
			vhdl << tab <<declare("expRAdjusted",wEOut_) << " <= "<< use("expBias") <<" + "<< // + sign extended expAdj
				"(" << rangeAssign(wEOut_-1,countWidth_+1, use("expAdj")+of(countWidth_)) << " & " << use("expAdj")<<");"<<endl;
			vhdl << tab <<declare("excBits",2) << " <=\"01\";"<<endl;
		}
		else
			if (countWidth_+1 == wEOut_){
				vhdl << tab <<declare("expRAdjusted",wEOut_) <<" <= " << use("expBias") << " + " << use("expAdj") << ";"<<endl;
				vhdl << tab <<declare("excBits",2) << " <=\"01\";"<<endl;
			}	
			else
				{
					// case encountered when we wish to put the contents of the accumulator in a floating point number 
					// with a very small exponent
					vhdl << tab <<declare("expAdjustedExt",countWidth_+1)<<" <= " << "(" << rangeAssign(countWidth_, wEOut_, "'0'") << " & " << use("expBias") << ") + "<<use("expAdj")<<";"<<endl;
			
					vhdl << tab <<declare("signExpAdjustedExt") << " <= "<< use("expAdjustedExt")<<of(countWidth_)<<";"<<endl;
					vhdl << tab <<declare("modulusExpAdjusted",countWidth_) << "<= ("<<rangeAssign(countWidth_-1,0, use("signExpAdjustedExt"))<<" xor "<<use("expAdjustedExt")<<range(countWidth_-1,0)<<") + "<<
						use("signExpAdjustedExt")<<";"<<endl;
				
					vhdl << tab <<declare("maxExponent",countWidth_)<<" <= CONV_STD_LOGIC_VECTOR("<<(mpz_class(1)<<wEOut_)-1<<" , "<<countWidth_<<");"<<endl;
			
					vhdl << tab <<declare("expOverflow")<<" <= '1' when (("<<use("modulusExpAdjusted")<< " > " << use("maxExponent") << ") and ("<<use("signExpAdjustedExt")<<"='0')) else '0';"<<endl;
					vhdl << tab <<declare("expUnderflow")<<" <= '1' when (("<<use("modulusExpAdjusted")<< " > " << use("maxExponent") << ") and ("<<use("signExpAdjustedExt")<<"='1')) else '0';"<<endl;
			
					declare("excBits",2);
					vhdl << tab <<"excBits(1) <= " << use("expOverflow") << " and  not(" << use("expUnderflow")<<");"<<endl;
					vhdl << tab <<"excBits(0) <= not("<<use("expOverflow")<<") and not("<<use("expUnderflow")<<");"<<endl;
					
					vhdl << tab <<declare("expRAdjusted",wEOut_) << " <= " << use("expAdjustedExt")<< range(wEOut_-1,0)<<";"<<endl;
				}			
	
		vhdl << tab <<declare("excRes",2) << " <= " << use("excBits") << " when " << use("ovf_updated")<<"='0' else \"10\";"<<endl;

		//c2 of the fraction
		//	nextCycle();
		vhdl << tab << declare("notResFrac",wFOut_+1) << " <= not("<<use("resFrac")<<");"<<endl;

		//a possible carry propagation
		IntAdder *b = new IntAdder(target,wFOut_ + 1);
		oplist.push_back(b);
	
		inPortMap(b, "X", "notResFrac");
		inPortMapCst(b, "Y", zg(wFOut+1,0));
		inPortMap(b, "Cin", "resSign");
		outPortMap(b, "R", "postResFrac");
		vhdl << tab << instance(b, "carryPropagator");
	
		syncCycleFromSignal("postResFrac");		
	
		vhdl << tab << declare("resultFraction",wFOut_+1) << " <= "<<use("resFrac") << " when ("<<use("resSign")<<"='0') else "<<use("postResFrac")<<";"<<endl;
	
		vhdl << tab << "R <= "<<use("excRes") << " & " <<use("resSign")<<" & "<< use("expRAdjusted")<<range(wEOut_-1,0)<<" & "<<use("resultFraction") <<range(wFOut_-1,0) <<";"<<endl;
	}

	LongAcc2FP::~LongAcc2FP() {
	}

}
