/*
 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
 *
 * Authors : Florent de Dinechin, Bogdan Pasca
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
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "LZOCShifterSticky.hpp"

using namespace std;

LZOCShifterSticky::LZOCShifterSticky(Target* target, int wIn, int wOut, bool computeSticky, const int countType) :
	Operator(target), wIn_(wIn), wOut_(wOut), computeSticky_(computeSticky), countType_(countType) {
	
	ostringstream currLevel, prevLevel, currSticky, currCount;
	
	// -------- Parameter set up -----------------
	setEntityType( (countType_==-1?gen:spec) );
	wCount_ = intlog2(wIn_);
	cout << " wCount is: "<<wCount_<<endl;

	setOperatorName();
	setOperatorType();

	addInput ("I", wIn_);
	if (entityType_==gen) addInput ("OZb"); /* if we generate a generic LZOC */
	addOutput("Count", wCount_);
	addOutput("O", wOut_);
	if (computeSticky_)   addOutput("Sticky"); /* if we require a sticky bit computation */
	
	currLevel << "level" <<wCount_;
	currSticky <<"sticky"<<wCount_; 
	
	vhdl << tab << declare(currLevel.str(), wIn_) << " <= I ;"  <<endl; 
	if (computeSticky_) vhdl << tab << declare(currLevel.str(), 1  ) << " <= '0' ;"<<endl; //init sticky 
	
	for (int i=wCount_-1; i>=0; i--){
		currLevel.str(""); currLevel << "level" << i;
		prevLevel.str(""); prevLevel << "level" << i+1;
		currCount.str(""); currCount << "count" << i;
	
		vhdl << tab << declare(currCount.str(),1) << "<= '1' when " <<use(prevLevel.str())<<"("<< wIn_-1 <<" downto "<< wIn_ - intpow2(i) <<") = "
		     << "("<< wIn-1 <<" downto "<< wIn - intpow2(i) <<"=>"<< (countType_==-1? use("sozb"): countType_==0?"'0'":"'1'")<<") else '0';"<<endl;
		vhdl << tab << declare(currLevel.str(),wIn_) << "<= " << use(prevLevel.str()) << " when " << use(currCount.str()) << "='0' else "
		     << use(prevLevel.str()) << "(" << wIn - intpow2(i)-1 << " downto 0) & (" << intpow2(i)-1 << " downto 0 => '0');"<<endl; //FIXME should we pad with zero ?
//		     << (countType_==-1? "not("+use("sozb")+")": countType_==0?"'1'":"'0'")<<");"<<endl;
		      
	}
	vhdl << tab << "O <= "<<use(currLevel.str())<<";"<<endl;
	vhdl << tab << "Count <= ";
	for (int i=wCount_-1; i>=0; i--){
		currCount.str(""); currCount << "count" << i;
		vhdl << use(currCount.str());
		if (i>0) vhdl << " & ";
		else     vhdl << ";"<<endl;
	}
}

LZOCShifterSticky::~LZOCShifterSticky() {
}

void LZOCShifterSticky::setEntityType(entityType_t eType){
	entityType_ = eType; 
}

int LZOCShifterSticky::getCountWidth() const{
	return wCount_;
}

void LZOCShifterSticky::setOperatorName(){
	ostringstream name; 
	ostringstream computationInterface;
	/* The computationInterface refers to the VHDL entity ports. 
	   If computation is generic, then an extra input port is 
	   available on the entity for specifying the count type */ 
	if (countType_ < 0) 
		computationInterface<<"gen";
	else 
		computationInterface<<"spec";
	
	name<<"LZOCShifterSticky_"<<wIn_<<"_"<<wOut_<<"_"<<computeSticky_<<"_" //semantic name
	    <<computationInterface.str();
	
	uniqueName_=name.str();
}

void LZOCShifterSticky::outputVHDL(std::ostream& o, std::string name) {
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::stdLibs(o);
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);
	o << buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);
}


void LZOCShifterSticky::emulate(TestCase* tc)
{
	mpz_class si   = tc->getInputValue("I");
	
	mpz_class sozb = 42; //dummy value
	if (countType_ == -1) 
		sozb = tc->getInputValue("OZb");
	
	int sticky=0;
	int j, icount;
	
	/* Count the leading zero/one s */
	mpz_class bit = (countType_ == -1) ? sozb : (countType_ == 0 ? 0 : 1); /* what are we counting in the specific case */
	/* from the MSB towards the LSB, check if current bit of input = to the bit we test against */
	for (j = wIn_-1; j >= 0; j--)
			if (mpz_tstbit(si.get_mpz_t(), j) != bit)
				break;
	
	/* the number of bits is then equal to:
	 the index of the MSB - the index where we stoped previously */
	icount = (wIn_-1) - j;
	tc->addExpectedOutput("Count", icount);

	/* compute the max value on wOut_ bits */
	maxValue_ = mpzpow2(wOut_)-1;
	mpz_class inputValue = si;
	
	//compute output value and sticky
	mpz_class stickyTest =  (countType_==-1) ? (sozb==0?1:0) : (countType_ == 0 ? 1 : 0) ;
	if ((countType_==0) || (sozb==0)) 
		if (inputValue > 0) 
			while (!((inputValue<=maxValue_)&&(2*inputValue>maxValue_)))
				if (inputValue>maxValue_){
					if(mpz_tstbit(inputValue.get_mpz_t(), 0)==stickyTest)
						sticky=1;
					inputValue=inputValue/2;
				}else
					inputValue=inputValue*2;
		else {}
	else /* if we are counting ones */
	{
		int restOfBits = wIn_ - icount;
		if (icount>0){
			mpz_class ones = mpzpow2(icount)-1;
			ones *= mpzpow2(restOfBits);
			
			inputValue-=ones; // the input without the leading ones
		}

		if ((wIn_<=wOut_) || ((wIn_>wOut_) && (restOfBits<wOut_) ))	//shift result in place	
			inputValue *=mpzpow2(wOut_-restOfBits);
		else
			for (int i=1;i<=restOfBits-wOut_;i++){
				if(mpz_tstbit(inputValue.get_mpz_t(), 0)==stickyTest) //FIXME What do we count out when we count ones, the one or the zero?
					sticky=1;
				inputValue=inputValue/2;
			}
	}
			
	tc->addExpectedOutput("O",inputValue);
			
	if (computeSticky_){
		tc->addExpectedOutput("Sticky",sticky);
	}
}



