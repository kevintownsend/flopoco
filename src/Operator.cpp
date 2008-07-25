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
#include "Operator.hpp"


void Operator::addInput(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addInput, signal " << name<< " seems to already exist";
		throw o.str();
	}
	Signal *s = new Signal(name, Signal::in, width) ;
	ioList_.push_back(s);
	signalMap_[name] = s ;
	numberOfInputs_ ++;
}

void Operator::addOutput(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addInput, signal " << name << " seems to already exist";
		throw o.str();
	}
	Signal *s = new Signal(name, Signal::out, width) ;
	ioList_.push_back(s);
	signalMap_[name] = s ;
	numberOfOutputs_ ++;
}

void Operator::addFPInput(const std::string name, const int wE, const int wF) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s = new Signal(name, Signal::in, wE, wF) ;
	ioList_.push_back(s);
	signalMap_[name] = s ;
	numberOfInputs_ ++;
}

void Operator::addFPOutput(const std::string name, const int wE, const int wF) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s = new Signal(name, Signal::out, wE, wF) ;
	ioList_.push_back(s);
	signalMap_[name] = s ;
	numberOfOutputs_ ++;
}

void Operator::addSignal(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addSignal, signal " << name << " seems to already exist";
		throw o.str();  
	}
	Signal *s = new Signal(name, Signal::wire, width);
	signalMap_[name] = s ;
	signalList_.push_back(s);
}

void Operator::addSignalBus(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		std::ostringstream o;
		o << "ERROR in addSignal, signal " << name << " seems to already exist";
		throw o.str();  
	}
	Signal *s = new Signal(name, Signal::wire, width, true);
	signalMap_[name] = s ;
	signalList_.push_back(s);
}

void Operator::addRegisteredSignalWithoutReset(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s;
	s = new Signal(name, Signal::registeredWithoutReset, width);
	signalList_.push_back(s);
	signalMap_[name] = s ;

	string o =  name + "_d";
	if (signalMap_.find(o) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o, Signal::wire, width);
	signalList_.push_back(s);
	signalMap_[name] = s ;

	hasRegistersWithoutReset_ = true;
}

void Operator::addRegisteredSignalWithAsyncReset(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s;
	s = new Signal(name, Signal::registeredWithAsyncReset, width);
	signalList_.push_back(s);
	signalMap_[name] = s ;

	string o = name + "_d";
	if (signalMap_.find(o) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o, Signal::wire, width);
	signalList_.push_back(s);
	signalMap_[name] = s ;

	hasRegistersWithAsyncReset_ = true;
}

void Operator::addRegisteredSignalWithSyncReset(const std::string name, const int width) {
	if (signalMap_.find(name) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	Signal *s;
	s = new Signal(name, Signal::registeredWithSyncReset, width);
	signalList_.push_back(s);
	signalMap_[name] = s ;

	string o = name + "_d";
	if (signalMap_.find(o) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o, Signal::wire, width);
	signalList_.push_back(s);
	signalMap_[name] = s ;

	hasRegistersWithSyncReset_ = true;
}

string Operator::addDelaySignal(const string name, const int width, const int delay) {
	ostringstream o;
	Signal *s;
	o << name;
	// if the delay is zero it is equivalent to addSignal
	if (delay > 0) {
		for (int i=0; i<delay; i++){
			if (signalMap_.find(o.str()) != signalMap_.end()) {
				cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
				exit(EXIT_FAILURE);
			}
			s = new Signal(o.str(), Signal::registeredWithSyncReset, width);
			signalList_.push_back(s);    
			signalMap_[name] = s ;
			o  << "_d";
		}
		hasRegistersWithSyncReset_ = true;
	}

	if (signalMap_.find(o.str()) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o.str(), Signal::wire, width);
	signalList_.push_back(s);    
	signalMap_[name] = s ;
	
	return o.str();
}

string Operator::addDelaySignalNoReset(const string name, const int width, const int delay) {
	ostringstream o;
	Signal *s;
	o << name;
	// if the delay is zero it is equivalent to addSignal
	if (delay > 0) {
		for (int i=0; i<delay; i++){
			if (signalMap_.find(o.str()) != signalMap_.end()) {
				cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
				exit(EXIT_FAILURE);
			}
			s = new Signal(o.str(), Signal::registeredWithoutReset, width);
			signalList_.push_back(s);    
			signalMap_[name] = s ;
			o  << "_d";
		}
		hasRegistersWithoutReset_ = true;
	}

	if (signalMap_.find(o.str()) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o.str(), Signal::wire, width);
	signalList_.push_back(s);    
	signalMap_[name] = s ;
	
	return o.str();
}

string Operator::addDelaySignalBus(const string name, const int width, const int delay) {
	ostringstream o;
	Signal *s;
	o << name;
	// if the delay is zero it is equivalent to addSignal
	if( delay > 0) {
		for (int i=0; i<delay; i++){
			if (signalMap_.find(o.str()) != signalMap_.end()) {
				cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
				exit(EXIT_FAILURE);
			}
			s = new Signal(o.str(), Signal::registeredWithSyncReset, width, true);
			signalList_.push_back(s);    
			signalMap_[name] = s ;
			o  << "_d";
		}
		hasRegistersWithSyncReset_ = true;
	}

	if (signalMap_.find(o.str()) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o.str(), Signal::wire, width, true);
	signalList_.push_back(s);    
	signalMap_[name] = s ;
	
	return o.str();
}

string Operator::addDelaySignalBusNoReset(const string name, const int width, const int delay) {
	ostringstream o;
	Signal *s;
	o << name;
	// if the delay is zero it is equivalent to addSignal
	if (delay > 0) {
		for (int i=0; i<delay; i++){
			if(signalMap_.find(o.str()) != signalMap_.end()) {
				cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
				exit(EXIT_FAILURE);
			}
			s = new Signal(o.str(), Signal::registeredWithoutReset, width, true);
			signalList_.push_back(s);    
			signalMap_[name] = s ;
			o  << "_d";
		}
		hasRegistersWithoutReset_ = true;
	}

	if (signalMap_.find(o.str()) != signalMap_.end()) {
		cerr << "ERROR in addInput , signal " << name<< " seems to already exist" << endl;
		exit(EXIT_FAILURE);
	}
	s = new Signal(o.str(), Signal::wire, width, true);
	signalList_.push_back(s);    
	signalMap_[name] = s ;
	
	return o.str();
}

string Operator::getDelaySignalName(const string name, const int delay) {
	ostringstream o;
	if (delay==0 || isSequential()==false)
		return name;
	else {
		o << name;
		for (int i=0; i<delay; i++){
			o  << "_d";
		}
		return o.str();
	}
}

Signal * Operator::getSignalByName(string name) {
	return signalMap_[name];
}

void Operator::setOperatorName(std::string prefix, std::string postfix){
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

void Operator::setOperatorName(std::string operatorName){
	uniqueName_ = operatorName;
}

void Operator::setOperatorType(){
	if (target_->isPipelined())
		setSequential();
	else
		setCombinatorial();	
}

void Operator::setCommentedName(std::string name){
	commentedName_ = name;
}

string Operator::getOperatorName() const{
  return uniqueName_;
}

int Operator::getIOListSize() const{
  return ioList_.size();
}

vector<Signal*> * Operator::getIOList(){
  return &ioList_; 
}

const Signal * Operator::getIOListSignal(int i){
  return ioList_[i];
}
			
void  Operator::outputVHDLSignalDeclarations(std::ostream& o) {
	for (int i=0; i < this->signalList_.size(); i++){
		Signal* s = this->signalList_[i];
		o<<tab<<  s->toVHDL() << ";" << endl;
	}
}

void  Operator::outputVHDLRegisters(std::ostream& o) {

	// First registers without a reset
	if (hasRegistersWithoutReset_) {
		o << tab << "process(clk)  begin\n"
			<< tab << tab << "if clk'event and clk = '1' then\n";
		for(int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			if(s->type()==Signal::registeredWithoutReset) 
				o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  " << s->getSignalName() <<";\n";
		}
		o << tab << tab << "end if;\n";
		o << tab << "end process;\n"; 
	}
	
	// then registers with a reset
	if (hasRegistersWithAsyncReset_) {
		o << tab << "process(clk, rst)" << endl;
		o << tab << tab << "begin" << endl;
		o << tab << tab << tab << "if rst = '1' then" << endl;
		for(int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			if(s->type()==Signal::registeredWithAsyncReset)
				 if ((s->width()>1)||(s->isBus())) 
								 o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
				 else
						o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  '0';\n";
		}
		o << tab << tab << tab << "elsif clk'event and clk = '1' then" << endl;
		for(int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			if(s->type()==Signal::registeredWithAsyncReset) 
				o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  " << s->getSignalName() <<";\n";
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
		for(int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			if(s->type()==Signal::registeredWithSyncReset)
				 if ((s->width()>1)||(s->isBus())) 
								 o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  (" << s->width()-1 <<" downto 0 => '0');\n";
				 else
						o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  '0';\n";
		}
		o << tab << tab << tab << "else" << endl;
		for(int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			if(s->type()==Signal::registeredWithSyncReset) 
				o << tab <<tab << tab << s->getSignalName() <<"_d" << " <=  " << s->getSignalName() <<";\n";
		}
		o << tab << tab << tab << "end if;" << endl;
		o << tab << tab << "end if;" << endl;
		o << tab << tab << "end process;" << endl;
	}
}

void Operator::outputVHDLComponent(std::ostream& o, std::string name) {
	o << tab << "component " << name << " is" << endl;
	if (ioList_.size() > 0)
	{
		o << tab << tab << "port ( ";
		for (int i=0; i<this->ioList_.size(); i++){
			Signal* s = this->ioList_[i];
			if (i>0) // align signal names 
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
	o << "entity " << uniqueName_ << " is" << endl;
	if (ioList_.size() > 0)
	{
		o << tab << "port ( ";

		for (int i=0; i<this->ioList_.size(); i++){
			Signal* s = this->ioList_[i];
			if (i>0) // align signal names 
				o<<"          ";
			o<<  s->toVHDL();
			if(i < this->ioList_.size()-1)  o<<";" << endl;
		}
	
		o << tab << ");"<<endl;
	}
	o << "end entity;" << endl << endl;
}

void Operator::licence(std::ostream& o, std::string authorsyears){
	o<<"--------------------------------------------------------------------------------"<<endl;
	// centering the unique name
	int s;
	if(uniqueName_.size()<76) s = (76-uniqueName_.size())/2; else s=0;
	o<<"--"; for(int i=0; i<s; i++) o<<" "; o  << uniqueName_ << endl; 

	// if this operator was renamed from the command line, show the original name
	if(commentedName_!="") {
		if(commentedName_.size()<74) s = (74-commentedName_.size())/2; else s=0;
		o<<"--"; for(int i=0; i<s; i++) o<<" "; o << "(" << commentedName_ << ")" << endl; 
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
	addInput("clk");
	addInput("rst");
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

