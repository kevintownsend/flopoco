/*
 * A barrel shifter for FloPoCo
 *
 * Author : Florent de Dinechin, Bogdan Pasca
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

#include "Shifters.hpp"

using namespace std;

// TODO there is a small inefficiency here, as most bits of s don't need to be copied all the way down

Shifter::Shifter(Target* target, int wIn, int maxShift, ShiftDirection direction, map<string, double> inputDelays) :
	Operator(target), wIn_(wIn), maxShift_(maxShift), direction_(direction) {
	
	ostringstream name;
	
	name << "Shift"<<(direction==Left?"Left":"Right")<<"_"<<wIn_<<"_by_max_"<<maxShift_;
	uniqueName_ = name.str();


	// -------- Parameter set up -----------------
	wOut_         = wIn_ + maxShift_;
	wShiftIn_     = intlog2(maxShift_);
	maxInputDelay_ = getMaxInputDelays(inputDelays); 

	addInput ("X", wIn_);
	addInput ("S", wShiftIn_);  
	addOutput("R", wOut_);

	setOperatorType();

	vhdl << tab << declare("level0",wIn_  ) << "<= X;" <<endl;
	vhdl << tab << declare("ps", wShiftIn_) << "<= S;" <<endl;
	
	// local variables
	double period = (1.0)/target->frequency();
	double stageDelay = 0.0;
	int    lastRegLevel = -1;
	int    unregisteredLevels = 0;
	int    dep = 0;
		
	for (int currentLevel=0; currentLevel<wShiftIn_; currentLevel++){
		//compute current level delay
		//TODO REMEMBER WHY THIS DOES WORK
		unregisteredLevels = currentLevel - lastRegLevel;
		if ( intpow2(unregisteredLevels-1) > wIn_+currentLevel+1 ) 
			dep = wIn+currentLevel+1 + unregisteredLevels -1;
		else
			dep = intpow2(unregisteredLevels-1);
		
		if (verbose)
			cout<<"depth = "<<dep<<" at i="<<currentLevel<<endl;	

		stageDelay = intlog(mpz_class(target->lutInputs()), mpz_class(dep)) * target->lutDelay() + (intlog(mpz_class(target->lutInputs()),mpz_class(dep))-1) * target->localWireDelay();
		if (lastRegLevel == -1)
			stageDelay+= maxInputDelay_;
		
		if (stageDelay>period){
			lastRegLevel = currentLevel;
			nextCycle(); ////////////////////////////////////////////////
		}
			ostringstream currentLevelName, nextLevelName;
			currentLevelName << "level"<<currentLevel;
			nextLevelName << "level"<<currentLevel+1;
			if (direction==Right){
				vhdl << tab << declare(nextLevelName.str(),wIn+intpow2(currentLevel+1)-1 ) 
					        <<"<=  ("<<intpow2(currentLevel)-1 <<" downto 0 => '0') & "<<use(currentLevelName.str())<<" when "<<use("ps")<<"("<<currentLevel<<") = '1' else "
					 << tab << use(currentLevelName.str()) <<" & ("<<intpow2(currentLevel)-1<<" downto 0 => '0');"<<endl;
			}else{
				vhdl << tab << declare(nextLevelName.str(),wIn+intpow2(currentLevel+1)-1 ) 
					        << "<= " << use(currentLevelName.str()) << " & ("<<intpow2(currentLevel)-1 <<" downto 0 => '0') when "<<use("ps")<<"("<<currentLevel<<") = '1' else "
					 << tab <<" ("<<intpow2(currentLevel)-1<<" downto 0 => '0') & "<< use(currentLevelName.str()) <<";"<<endl;
			}
			
	}
	ostringstream lastLevelName;
	lastLevelName << "level"<<wShiftIn_;
	if (direction==Right)
		vhdl << tab << "R <= "<<use(lastLevelName.str())<<"("<< wIn + intpow2(wShiftIn_)-1-1 << " downto " << wIn_ + intpow2(wShiftIn_)-1 - wOut_ <<");"<<endl;
	else
		vhdl << tab << "R <= "<<use(lastLevelName.str())<<"("<< wOut_-1 << " downto 0);"<<endl;

}

Shifter::~Shifter() {
}

void Shifter::setOperatorName(){
	ostringstream name;
	if(direction_==Left) name <<"LeftShifter_";
	else                 name <<"RightShifter_";
	name<<wIn_<<"_by_max_"<<maxShift_;
	uniqueName_=name.str();
}

void Shifter::outputVHDL(std::ostream& o, std::string name) {
	ostringstream signame;
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007,2008,2009)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);
	o << buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);
}

void Shifter::emulate(TestCase* tc)
{
	mpz_class sx = tc->getInputValue("X");
	mpz_class ss = tc->getInputValue("S");
	mpz_class sr ;
	
	while (ss > maxShift_)
		ss = getLargeRandom(wShiftIn_);

	mpz_class shiftAmmount;
	if (direction_==Right)
		shiftAmmount=maxShift_-ss;
	else
		shiftAmmount=ss;
	
	mpz_class shiftedInput;
	shiftedInput=sx;

	int i;
	for (i=0;i<shiftAmmount;i++)
			shiftedInput=shiftedInput*2;
	
	sr=shiftedInput;

	tc->addExpectedOutput("R", sr);

}



