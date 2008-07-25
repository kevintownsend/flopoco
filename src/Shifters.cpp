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

Shifter::Shifter(Target* target, int wIn, int maxShift, ShiftDirection direction) :
	Operator(target), wIn_(wIn), maxShift_(maxShift), direction_(direction) {

	wOut_ = wIn_ + maxShift_;
	wShiftIn_ = intlog2(maxShift_);
	setOperatorName();
	setOperatorType();

	// Set up the IO signals
	addInput ("X", wIn_);
	addInput ("S", wShiftIn_);  
	addOutput("R", wOut_);
 
	// evaluate the pipeline
	double criticalPath = 0.0;
	for (int i=0; i<wShiftIn_; i++) 
		levelRegistered_[i] = false;

	if(isSequential()) {
		for (int i=0; i<wShiftIn_; i++) {
			/* approximate delay of this stage */
			double stageDelay = /*target->lut_delay() + */ 1.2 * target->localWireDelay() * (1<<i) ;
			if (criticalPath + stageDelay > 1/target->frequency()) {
				criticalPath=stageDelay; 	// reset critical path
				levelRegistered_[i+1]= true;
				incrementPipelineDepth();
			}
			else{ 
				criticalPath += stageDelay; // augment critical path
				levelRegistered_[i+1] = false;
			}
		}
		// register the last level anyway
		if(!levelRegistered_[wShiftIn_]) {
			levelRegistered_[wShiftIn_] = true;
			incrementPipelineDepth();
		}
	}

	// Set up the intermediate signals 
	if (isSequential()){
		addSignal("level0", wIn_);
		
		for (int i=1; i<=wShiftIn_; i++) {
			ostringstream sname;
			sname << "level"<<i;
			if (levelRegistered_[i])
				addRegisteredSignalWithSyncReset(sname.str(), wIn_ + (1<<i) -1 );
			else
				addSignal(sname.str(), wIn_ + (1<<i) -1);
		}

		// The shift input has to be delayed as well now
		if(getPipelineDepth()>=1) 
			addDelaySignal("ps", wShiftIn_, getPipelineDepth()-1); 
		else
			addSignal("ps", wShiftIn_);

	}
	else{
		for (int i=0; i<=wShiftIn_; i++) {
				ostringstream sname;
				sname << "level"<<i;
				addSignal(sname.str(), wIn_ + (1<<i) -1);
		}	
	
		addSignal("ps", wShiftIn_);
	}
	
	if (verbose){
		cout <<"wShiftIn="<<wShiftIn_<<endl;
	}
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
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	if (isSequential())
	{
		outputVHDLRegisters(o);
		
		int stage=0;
		o << "   level0 <=  X ;" << endl;
		o << "   ps <=  s;" <<endl;
		ostringstream psname;
		psname << "ps";
		for (int i=0; i<wShiftIn_; i++) {
			ostringstream lname;
			lname << "level"<<i;
			if (levelRegistered_[i]) { // use the registered signal instead
				lname << "_d";
				// and use next stage of ps
				psname << "_d",
				// add a synchronisation barrier here
				o <<"  ----- synchro barrier ------- " <<endl;
				stage++;
			}
							
			o << "   level"<<i+1<<" <=  ";
			o << "("<<intpow2(i)-1<<" downto 0 => '0') & "<<lname.str()
				<<"  when "<<psname.str()<<"("<<i<<") = '"<<(direction_==Right?1:0)<<"'   else  ";
			o << lname.str()<<" & ("<<intpow2(i)-1<<" downto 0 => '0');" << endl;
			
		}
		if(levelRegistered_[wShiftIn_])
			o <<"  ----- synchro barrier ------- " <<endl;

		o << "   R <=  level"<<wShiftIn_;
		if(levelRegistered_[wShiftIn_]) 
			o << "_d";
		if (direction_==Left)
			o << "("<< wOut_-1<<" downto 0);" << endl << endl;
		else
			o << "("<< wIn_ + (1<<wShiftIn_) -2<<" downto "<<wIn_ + (1<<wShiftIn_) -1 - wOut_ <<");"<<endl << endl;
			
	}
	else
	{ //combinatorial version
		o << "   level0 <=  X ;" << endl;
		o << "   ps <=  s;" <<endl;
		ostringstream psname;
		psname << "ps";
		
		for (int i=0; i<wShiftIn_; i++) {
			ostringstream lname;
			lname << "level"<<i;
									
			o << " level"<<i+1<<" <=  ";
			o << "("<<intpow2(i)-1<<" downto 0 => '0') & "<<lname.str()
				<<"  when "<<psname.str()<<"("<<i<<") = '"<<(direction_==Right?1:0)<<"'   else  ";
			o << lname.str()<<" & ("<<intpow2(i)-1<<" downto 0 => '0');" << endl;
		}
		if (direction_==Left)
			o << "   R <=  level"<<wShiftIn_<<"("<< wOut_-1<<" downto 0);" << endl << endl;		
		else
			o << "   R <=  level"<<wShiftIn_<<"("<< wIn_ + (1<<wShiftIn_) -2<<" downto "<<wIn_ + (1<<wShiftIn_) -1 - wOut_ <<");" << endl << endl;
	}
	
	endArchitecture(o);		
}

TestIOMap Shifter::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("X"));
	tim.add(*getSignalByName("S"));
	tim.add(*getSignalByName("R"));
	return tim;
}

void Shifter::fillTestCase(mpz_class a[])
{
	mpz_class& sx   = a[0];
	mpz_class& ss   = a[1];
	mpz_class& sr   = a[2];
				
	while (ss>maxShift_)
	ss=getLargeRandom(wShiftIn_);

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
}
