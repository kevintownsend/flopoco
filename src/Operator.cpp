/*
 * The base Operator class, every operator should inherit it
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


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "Operator.hpp"
#include "utils.hpp"


void Operator::addInput(const std::string name, const int width, const bool isBus) {
	if (signalMap_.find(name) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addInput, signal " << name<< " seems to already exist";
		throw o.str();
	}
	Signal *s = new Signal(name, Signal::in, width, isBus) ; // default TTL and cycle OK
	s->setCycle(0);
	ioList_.push_back(s);
	signalMap_[name] = s ;
	numberOfInputs_ ++;
}

void Operator::addOutput(const std::string name, const int width, const int numberOfPossibleOutputValues, const bool isBus) {
	if (signalMap_.find(name) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addInput, signal " << name << " seems to already exist";
		throw o.str();
	}
	Signal *s = new Signal(name, Signal::out, width, isBus) ;
	s -> setNumberOfPossibleValues(numberOfPossibleOutputValues);
	ioList_.push_back(s);
	for(int i=0; i<numberOfPossibleOutputValues; i++) 
		testCaseSignals_.push_back(s);
	signalMap_[name] = s ;
	numberOfOutputs_ ++;
}

void Operator::addFPInput(const std::string name, const int wE, const int wF) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s = new Signal(name, Signal::in, wE, wF);
	s->setCycle(0);
	ioList_.push_back(s);
	signalMap_[name] = s ;
	numberOfInputs_ ++;
}

void Operator::addFPOutput(const std::string name, const int wE, const int wF, const int numberOfPossibleOutputValues) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s = new Signal(name, Signal::out, wE, wF) ;
	s -> setNumberOfPossibleValues(numberOfPossibleOutputValues);
	ioList_.push_back(s);
	for(int i=0; i<numberOfPossibleOutputValues; i++) 
		testCaseSignals_.push_back(s);
	signalMap_[name] = s ;
	numberOfOutputs_ ++;
}




void Operator::addSignalGeneric(const string name, const int width, const int delay, Signal::SignalType regType, bool isbus) {
	ostringstream o;
	Signal *s;

	o << name;
	if (isSequential() && delay > 0) { 	// if delay<=0,  it is equivalent to addSignal
		for (int i=0; i<delay; i++){
			if(signalMap_.find(o.str()) != signalMap_.end()) {
				std::ostringstream o;
				o << "ERROR in addSignalGeneric, signal " << name<< " seems to already exist";
				throw o.str();
			}
			s = new Signal(o.str(), regType, width, isbus);
			signalList_.push_back(s);    
			signalMap_[o.str()] = s ;
			o  << "_d";
		}
		if(regType==Signal::registeredWithoutReset)
			hasRegistersWithoutReset_ = true;
		if(regType==Signal::registeredWithSyncReset)
			hasRegistersWithSyncReset_ = true;
		if(regType==Signal::registeredWithAsyncReset)
			hasRegistersWithAsyncReset_ = true;
	}

	if (signalMap_.find(o.str()) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addSignalGeneric, signal " << name<< " seems to already exist";
		throw o.str();
	}
	s = new Signal(o.str(), Signal::wire, width, isbus);
	signalList_.push_back(s);    
	signalMap_[o.str()] = s ;
}




void Operator::addSignal(const std::string name, const int width) {
	addSignalGeneric(name,  width, 0, Signal::wire, //unused
						  false);
}



void Operator::addSignalBus(const std::string name, const int width) {
	addSignalGeneric(name,  width, 0, Signal::wire, //unused
						  true);
}

void Operator::addDelaySignal(const string name, const int width, const int delay) {
	addSignalGeneric(name,  width, delay, Signal::registeredWithoutReset, false);
}

void Operator::addDelaySignalBus(const string name, const int width, const int delay) {
	addSignalGeneric(name,  width, delay, Signal::registeredWithoutReset, true);
}

void Operator::addDelaySignalSyncReset(const string name, const int width, const int delay) {	
	addSignalGeneric(name,  width, delay, Signal::registeredWithSyncReset, false);
}

void Operator::addDelaySignalBusSyncReset(const string name, const int width, const int delay) {
	addSignalGeneric(name,  width, delay, Signal::registeredWithSyncReset, true);
}


string Operator::delaySignal(const string name, const int delay) {
	ostringstream o;
	Signal* s;
	bool isDeclared=true;
	
	if(signalMap_.find(name) ==  signalMap_.end()) {
		cerr << "WARNING in delaySignal, signal " << name << " not declared through addSignal" << endl;
		isDeclared=false;
	}

	if (delay<=0 || isSequential()==false)
		return name;
	else {
		if(isDeclared) {
			s=getSignalByName(name);
		}
		o << name;
		for (int i=0; i<delay; i++){
			o  << "_d";
		}
		return o.str();
	}
}

Signal * Operator::getSignalByName(string name) {
	ostringstream e;
	if(signalMap_.find(name) ==  signalMap_.end()) {
		e << "ERROR in getSignalByName, signal " << name<< " not declared";
		throw e.str();
	}
	return signalMap_[name];
}


void Operator::setName(std::string prefix, std::string postfix){
		ostringstream pr, po;
		if (prefix.length()>0)
			pr << prefix << "_"; 
		else 
			pr << "";
		if (postfix.length()>0)
			po << "_"<<postfix;
		else
			po << "";
		uniqueName_ = pr.str() + uniqueName_ + po.str();
}

void Operator::setName(std::string operatorName){
	uniqueName_ = operatorName;
}


void  Operator::changeName(std::string operatorName){
	commentedName_ = uniqueName_;
	uniqueName_ = operatorName;
}


string Operator::getName() const{
  return uniqueName_;
}

int Operator::getIOListSize() const{
  return ioList_.size();
}

vector<Signal*> * Operator::getIOList(){
  return &ioList_; 
}

Signal * Operator::getIOListSignal(int i){
  return ioList_[i];
}
			
 

void  Operator::outputVHDLSignalDeclarations(std::ostream& o) {
	for (unsigned int i=0; i < this->signalList_.size(); i++){
		Signal* s = this->signalList_[i];
		o<<tab<<  s->toVHDL() << ";" << endl;
	}

}

void  Operator::outputVHDLRegisters(std::ostream& o) {
	unsigned int i;
	// execute only if the operator is sequential, otherwise output nothing
	if (isSequential()){
		// First registers without a reset
		if (hasRegistersWithoutReset_) {
			o << tab << "process(clk)  begin\n"
			<< tab << tab << "if clk'event and clk = '1' then\n";
			for(i=0; i<signalList_.size(); i++) {
				Signal *s = signalList_[i];
				if(s->type()==Signal::registeredWithoutReset) 
					o << tab <<tab << tab << s->getName() <<"_d" << " <=  " << s->getName() <<";\n";
			}
			o << tab << tab << "end if;\n";
			o << tab << "end process;\n"; 
		}
		
		// then registers with a reset
		if (hasRegistersWithAsyncReset_) {
			o << tab << "process(clk, rst)" << endl;
			o << tab << tab << "begin" << endl;
			o << tab << tab << tab << "if rst = '1' then" << endl;
			for(i=0; i<signalList_.size(); i++) {
				Signal *s = signalList_[i];
				if(s->type()==Signal::registeredWithAsyncReset) {
					if ((s->width()>1)||(s->isBus())) 
						o << tab <<tab << tab << s->getName() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
					else
						o << tab <<tab << tab << s->getName() <<"_d" << " <=  '0';\n";
				}
			}
			o << tab << tab << tab << "elsif clk'event and clk = '1' then" << endl;
			for(i=0; i<signalList_.size(); i++) {
				Signal *s = signalList_[i];
				if(s->type()==Signal::registeredWithAsyncReset) 
					o << tab <<tab << tab << s->getName() <<"_d" << " <=  " << s->getName() <<";\n";
			}
			o << tab << tab << tab << "end if;" << endl;
			o << tab << tab << "end process;" << endl;
		}
		
		// then registers with synchronous reset
		if (hasRegistersWithSyncReset_) {
			o << tab << "process(clk, rst)" << endl;
			o << tab << tab << "begin" << endl;
			o<<  "    if clk'event and clk = '1' then" << endl;
			o << tab << tab << tab << "if rst = '1' then" << endl;
			for(i=0; i<signalList_.size(); i++) {
				Signal *s = signalList_[i];
				if(s->type()==Signal::registeredWithSyncReset) {
					if ((s->width()>1)||(s->isBus())) 
						o << tab <<tab << tab << s->getName() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
					else
						o << tab <<tab << tab << s->getName() <<"_d" << " <=  '0';\n";
				}
			}
			o << tab << tab << tab << "else" << endl;
			for(i=0; i<signalList_.size(); i++) {
				Signal *s = signalList_[i];
				if(s->type()==Signal::registeredWithSyncReset) 
					o << tab <<tab << tab << s->getName() <<"_d" << " <=  " << s->getName() <<";\n";
			}
			o << tab << tab << tab << "end if;" << endl;
			o << tab << tab << "end if;" << endl;
			o << tab << tab << "end process;" << endl;
		}
	}
}

void Operator::outputVHDLComponent(std::ostream& o, std::string name) {
	unsigned int i;
	o << tab << "component " << name << " is" << endl;
	if (ioList_.size() > 0)
	{
		o << tab << tab << "port ( ";
		if(isSequential()) {
			// add clk and rst signals which are no longer member of iolist
			o << "clk, rst : in std_logic;" <<endl;
		}
		for (i=0; i<this->ioList_.size(); i++){
			Signal* s = this->ioList_[i];
			if (i>0 || isSequential()) // align signal names 
				o<<tab<<"          ";
			o<<  s->toVHDL();
			if(i < this->ioList_.size()-1)  o<<";" << endl;
		}
		o << tab << ");"<<endl;
	}
	o << tab << "end component;" << endl;
}

void Operator::outputVHDLComponent(std::ostream& o) {
	this->outputVHDLComponent(o,  this->uniqueName_); 
}


void Operator::outputVHDLEntity(std::ostream& o) {
	unsigned int i;
	o << "entity " << uniqueName_ << " is" << endl;
	if (ioList_.size() > 0)
	{
		o << tab << "port ( ";

		if(isSequential()) {
			// add clk and rst signals which are no longer member of iolist
			o << "clk, rst : in std_logic;" <<endl;
		}
		for (i=0; i<this->ioList_.size(); i++){
			Signal* s = this->ioList_[i];
			if (i>0 || isSequential()) // align signal names 
				o<<"          ";
			o<<  s->toVHDL();
			if(i < this->ioList_.size()-1)  o<<";" << endl;
		}
	
		o << tab << ");"<<endl;
	}
	o << "end entity;" << endl << endl;
}


void Operator::setCopyrightString(std::string authorsYears){
	copyrightString_ = authorsYears;
}


void Operator::licence(std::ostream& o){
	licence(o, copyrightString_);
}


void Operator::licence(std::ostream& o, std::string authorsyears){
	o<<"--------------------------------------------------------------------------------"<<endl;
	// centering the unique name
	int s, i;
	if(uniqueName_.size()<76) s = (76-uniqueName_.size())/2; else s=0;
	o<<"--"; for(i=0; i<s; i++) o<<" "; o  << uniqueName_ << endl; 

	// if this operator was renamed from the command line, show the original name
	if(commentedName_!="") {
		if(commentedName_.size()<74) s = (74-commentedName_.size())/2; else s=0;
		o<<"--"; for(i=0; i<s; i++) o<<" "; o << "(" << commentedName_ << ")" << endl; 
	}

	o<<"-- This operator is part of the Infinite Virtual Library FloPoCoLib"<<endl
	 <<"-- and is distributed under the terms of the GNU Lesser General Public Licence."<<endl
	 <<"-- Authors: " << authorsyears <<endl
	 <<"--------------------------------------------------------------------------------"<<endl;
}

void Operator::outputVHDL(std::ostream& o) {
	this->outputVHDL(o, this->uniqueName_); 
}

bool Operator::isSequential() {
	return isSequential_; 
}

void Operator::setSequential() {
	isSequential_=true; 
	//	addInput("clk");
	// addInput("rst");
}

void Operator::setCombinatorial() {
	isSequential_=false; 
}

int Operator::getPipelineDepth() {
	return pipelineDepth_; 
}

void Operator::incrementPipelineDepth() {
	pipelineDepth_++; 
}

void Operator::setPipelineDepth(int d) {
	pipelineDepth_=d; 
}

void Operator::outputFinalReport() {
	cout << "Entity " << uniqueName_ <<":"<< endl;
	if(this->getPipelineDepth()!=0)
		cout << tab << "Pipeline depth = " << getPipelineDepth() << endl;
	else
		cout << tab << "Not pipelined"<< endl;
}



void Operator::setCycle(int cycle, bool report) {
	if(isSequential()) {
		currentCycle_=cycle;
		if(report)
			vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl ;
		// automatically update pipeline depth of the operator 
		if (currentCycle_ > pipelineDepth_) 
			pipelineDepth_ = currentCycle_;
	}
}

void Operator::nextCycle(bool report) {
	if(isSequential()) {
		currentCycle_ ++; 
		if(report)
			vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl ;
		// automatically update pipeline depth of the operator 
		if (currentCycle_ > pipelineDepth_) 
			pipelineDepth_ = currentCycle_;
	}
}


void Operator::setCycleFromSignal(string name, bool report) {
	ostringstream e;
	e << "ERROR in syncCycleFromSignal, "; // just in case

	if(isSequential()) {
		Signal* s;
		try {
			s=getSignalByName(name);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}

		if( s->getCycle() < 0 ) {
			ostringstream o;
			o << "signal " << name<< " doesn't have (yet?) a valid cycle";
		throw o.str();
		} 
		currentCycle_ = s->getCycle();
		if(report)
			vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl ;
		// automatically update pipeline depth of the operator 
		if (currentCycle_ > pipelineDepth_) 
			pipelineDepth_ = currentCycle_;
	}
}


void Operator::syncCycleFromSignal(string name, bool report) {
	ostringstream e;
	e << "ERROR in syncCycleFromSignal, "; // just in case

	if(isSequential()) {
		Signal* s;
		try {
			s=getSignalByName(name);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}

		if( s->getCycle() < 0 ) {
			ostringstream o;
			o << "signal " << name << " doesn't have (yet?) a valid cycle";
		throw o.str();
		} 
		// advance cycle if needed
		if (s->getCycle()>currentCycle_)
			currentCycle_ = s->getCycle();

		if(report)
			vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl ;
		// automatically update pipeline depth of the operator 
		if (currentCycle_ > pipelineDepth_) 
			pipelineDepth_ = currentCycle_;
	}
}



string Operator::declare(string name, const int width, bool isbus) {
	Signal* s;
	ostringstream e;
	// check the signals doesn't already exist
	if(signalMap_.find(name) !=  signalMap_.end()) {
		e << "ERROR in declare(), signal " << name<< " already exists";
		throw e.str();
	}
	// construct the signal (lifeSpan and cycle are reset to 0 by the constructor)
	s = new Signal(name, Signal::wire, width, isbus);
	// define its cycle 
	if(isSequential())
		s->setCycle(this->currentCycle_);
	// add the signal to signalMap and signalList
	signalList_.push_back(s);    
	signalMap_[name] = s ;
	return name;
}



string Operator::use(string name) {
	ostringstream e;
	e << "ERROR in use(), "; // just in case
	
	if(isSequential()) {
		Signal *s;
		try {
			s=getSignalByName(name);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}
		if(s->getCycle() < 0) {
			e << "signal " << name<< " doesn't have (yet?) a valid cycle";
			throw e.str();
		} 
		if(s->getCycle() > currentCycle_) {
			ostringstream e;
			e << "active cycle of signal " << name<< " is later than current cycle, cannot delay it";
			throw e.str();
		} 
		// update the lifeSpan of s
		s->updateLifeSpan( currentCycle_ - s->getCycle() );
		return s->delayedName( currentCycle_ - s->getCycle() );
	}
	else
		return name;
}




void Operator::outPortMap(Operator* op, string componentPortName, string actualSignalName){
	Signal* formal;
	Signal* s;
	ostringstream e;
	e << "ERROR in outPortMap(), "; // just in case
	// check the signals doesn't already exist
	if(signalMap_.find(actualSignalName) !=  signalMap_.end()) {
		e << "signal " << actualSignalName << " already exists";
		throw e.str();
	}
	try {
		formal=op->getSignalByName(componentPortName);
	}
	catch (string e2) {
		e << endl << tab << e2;
		throw e.str();
	}
	if (formal->type()!=Signal::out){
		e << "signal " << componentPortName << " of component " << op->getName() 
		  << " doesn't seem to be an output port";
		throw e.str();
	}
	int width = formal -> width();
	bool isbus = formal -> isBus();
	// construct the signal (lifeSpan and cycle are reset to 0 by the constructor)
	s = new Signal(actualSignalName, Signal::wire, width, isbus);
	// define its cycle 
	if(isSequential())
		s->setCycle( this->currentCycle_ + op->getPipelineDepth() );
	// add the signal to signalMap and signalList
	signalList_.push_back(s);    
	signalMap_[actualSignalName] = s ;

	// add the mapping to the mapping list of Op
	op->portMap_[componentPortName] = actualSignalName;
}


void Operator::inPortMap(Operator* op, string componentPortName, string actualSignalName){
	Signal* formal;
	ostringstream e;
	string name;
	e << "ERROR in inPortMap(), "; // just in case
	
	if(isSequential()) {
		Signal *s;
		try {
			s=getSignalByName(actualSignalName);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}
		if(s->getCycle() < 0) {
			ostringstream e;
			e << "signal " << actualSignalName<< " doesn't have (yet?) a valid cycle";
			throw e.str();
		} 
		if(s->getCycle() > currentCycle_) {
			ostringstream e;
			e << "active cycle of signal " << actualSignalName<< " is later than current cycle, cannot delay it";
			throw e.str();
		} 
		// update the lifeSpan of s
		s->updateLifeSpan( currentCycle_ - s->getCycle() );
		name = s->delayedName( currentCycle_ - s->getCycle() );
	}
	else
		name = actualSignalName;

	try {
		formal=op->getSignalByName(componentPortName);
	}
	catch (string e2) {
		e << endl << tab << e2;
		throw e.str();
	}
	if (formal->type()!=Signal::in){
		e << "signal " << componentPortName << " of component " << op->getName() 
		  << " doesn't seem to be an input port";
		throw e.str();
	}

	// add the mapping to the mapping list of Op
	op->portMap_[componentPortName] = name;
}



void Operator::inPortMapCst(Operator* op, string componentPortName, string actualSignal){
	Signal* formal;
	ostringstream e;
	string name;
	e << "ERROR in inPortMap(), "; // just in case

	try {
		formal=op->getSignalByName(componentPortName);
	}
	catch (string e2) {
		e << endl << tab << e2;
		throw e.str();
	}
	if (formal->type()!=Signal::in){
		e << "signal " << componentPortName << " of component " << op->getName() 
		  << " doesn't seem to be an input port";
		throw e.str();
	}

	// add the mapping to the mapping list of Op
	op->portMap_[componentPortName] = actualSignal;
}


string Operator::instance(Operator* op, string instanceName){
	ostringstream o;
	// TODO add checks here? Check that all the signals are covered for instance
	
	o << tab << instanceName << ": " << op->getName();
	if (isSequential()) 
		o << "  -- pipelineDepth="<< op->getPipelineDepth();
	o << endl;
	o << tab << tab << "port map (";
	// build vhdl and erase portMap_
	map<string,string>::iterator it;
	if(isSequential()) {
		o <<            " clk  => clk, " << endl;
		o <<  tab <<tab << "           rst  => rst, " << endl;
	}
	it=op->portMap_.begin();
	if(isSequential()) 
		o << tab << tab << "           " ;
	else
		o <<  " " ;
	o<< (*it).first << " => "  << (*it).second;
	//op->portMap_.erase(it);
	it++;
	for (  ; it != op->portMap_.end(); it++ ) {
		o << "," << endl;
		o <<  tab << tab << "           " << (*it).first << " => "  << (*it).second;
		//op->portMap_.erase(it);
	}
	o << ");" << endl;

	// add the operator to the subcomponent list 
	subComponents_[op->getName()]  = op;
	return o.str();
}
	



string Operator::buildVHDLSignalDeclarations() {
	ostringstream o;
	for(unsigned int i=0; i<signalList_.size(); i++) {
		Signal *s = signalList_[i];
		o << s->toVHDLDeclaration() << endl;
	}
	//now the signals from the I/O List which have the cycle>0
	for (unsigned int i=0; i<ioList_.size(); i++) {
		Signal *s = ioList_[i];
		if (s->getLifeSpan()>0){
			o << s->toVHDLDeclaration() << endl;	
		}
		
	}
	
	return o.str();	
}



string Operator::buildVHDLComponentDeclarations() {
	ostringstream o;
	for(map<string, Operator*>::iterator it = subComponents_.begin(); it !=subComponents_.end(); it++) {
		Operator *op = it->second;
		op->outputVHDLComponent(o);
		o<< endl;
	}
	return o.str();	
}

string  Operator::buildVHDLRegisters() {
	ostringstream o;

	// execute only if the operator is sequential, otherwise output nothing
	if (isSequential()){
		o << tab << "process(clk)  begin\n"
		  << tab << tab << "if clk'event and clk = '1' then\n";
		for(unsigned int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			if(s->getLifeSpan() >0) {
				for(int j=1; j <= s->getLifeSpan(); j++)
					o << tab <<tab << tab << s->delayedName(j) << " <=  " << s->delayedName(j-1) <<";" << endl;
			}
		}
		for(unsigned int i=0; i<ioList_.size(); i++) {
			Signal *s = ioList_[i];
			if(s->getLifeSpan() >0) {
				for(int j=1; j <= s->getLifeSpan(); j++)
					o << tab <<tab << tab << s->delayedName(j) << " <=  " << s->delayedName(j-1) <<";" << endl;
			}
		}

		o << tab << tab << "end if;\n";
		o << tab << "end process;\n"; 
	}
	return o.str();
}


void Operator::buildStandardTestCases(TestCaseList* tcl) {
	// Each operator should overload this method. If not, it is mostly harmless but deserves a warning.
	cerr << "WARNING: No standard test cases implemented for this operator" << endl;
}




void Operator::buildRandomTestCases(TestCaseList* tcl, int n){

	TestCase *tc;
	/* Generate test cases using random input numbers */
	for (int i = 0; i < n; i++) {
		tc = new TestCase(this); // TODO free all this memory when exiting TestBench
		/* Fill inputs */
		for (unsigned int j = 0; j < ioList_.size(); j++) {
			Signal* s = ioList_[j]; 
			if (s->type() == Signal::in) {
				mpz_class a = getLargeRandom(s->width());
				tc->addInput(s->getName(), a);
			}
		}
		/* Get correct outputs */
		emulate(tc);

		//		cout << tc->getInputVHDL();
		//    cout << tc->getExpectedOutputVHDL();


		// add to the test case list
		tcl->add(tc);
	}
}


void Operator::outputVHDL(std::ostream& o, std::string name) {
  
	licence(o);
	stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLComponentDeclarations();	
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);		
	o<<buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);

}



void Operator::emulate(TestCase * tc) {
		throw std::string("emulate() not implemented for ") + uniqueName_;
}
