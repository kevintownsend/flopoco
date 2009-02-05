#include <iostream>
#include <sstream>
#include "Signal.hpp"

using namespace std;

Signal::Signal(const string name, const Signal::SignalType type, const int width, const bool isBus) : 
	name_(name), type_(type), wE_(0), wF_(0), width_(width), low_(0), high_(width-1),
	isSubSignal_(false), isBus_(isBus),	isFP_(false), ttl_(0), cycle_(0), lifeSpan_(0) {
	updateSignalName();
}

Signal::Signal(const string name, const Signal::SignalType type, const int wE, const int wF) : 
	name_(name), type_(type), wE_(wE), wF_(wF), width_(wE+wF+3),
	low_(0), high_(width_-1), isSubSignal_(false),isBus_(false), isFP_(true), ttl_(0), cycle_(0), lifeSpan_(0)
{
	updateSignalName();
}

Signal::~Signal(){}

const string& Signal::getName() const { 
	return id_; 
}

void Signal::updateSignalName()	{
	if (isSubSignal_ == false)
		id_ = name_;
	else
		{
			stringstream o;
			if (width_ == 1)
				o << name_ << "(" << low_ << ")";
			else
				o << name_ << "(" << high_ << " downto " << low_ << ")";
			id_ = o.str();
		}
}

int Signal::width() const{return width_;}
	
int Signal::wE() const {return(wE_);}

int Signal::wF() const {return(wF_);}
	
bool Signal::isFP() const {return isFP_;}

bool Signal::isBus() const {return isBus_;}

Signal::SignalType Signal::type() const {return type_;}
	
string Signal::toVHDL() {
	ostringstream o; 
	if(type()==Signal::wire || type()==Signal::registeredWithoutReset || type()==Signal::registeredWithAsyncReset || type()==Signal::registeredWithSyncReset) 
		o << "signal ";
	o << getName();
	o << " : ";
	if (type()==Signal::in)
		o << "in ";
	if(type()==Signal::out)
		o << "out ";
	
	if ((1==width())&&(!isBus_)) 
		o << "std_logic" ;
	else 
		if(isFP_) 
			o << " std_logic_vector(" << wE() <<"+"<<wF() << "+2 downto 0)";
		else
			o << " std_logic_vector(" << width()-1 << " downto 0)";
	return o.str();
}



string Signal::delayedName(int delay){
	ostringstream o;
#if 0
	o << getName();
	if(delay>0) {
		for (int i=0; i<delay; i++){
			o  << "_d";
		}
	}
#else // someday we need to civilize pipe signal names
	o << getName();
	if(delay>0) 
		o << "_d" << delay;
#endif
		return o.str();
}


string Signal::toVHDLDeclaration() {
	ostringstream o; 
	o << "signal " << getName();
	for (int i=1; i<=lifeSpan_; i++) {
			o << ", " << getName() << "_d" << i;
	}
	o << " : ";
	if ((1==width())&&(!isBus_)) 
		o << "std_logic" ;
	else 
		if(isFP_) 
			o << " std_logic_vector(" << wE() <<"+"<<wF() << "+2 downto 0)";
		else
			o << " std_logic_vector(" << width()-1 << " downto 0)";
	o << ";";
	return o.str();
}

Signal Signal::getSubSignal(int low, int high)
{
	if (low < low_)
		throw string("Attempted to return subsignal with smaller low index.");
	if (high > high_)
		throw string("Attempted to return subsignal with bigger high index."); 
	
	Signal s(name_, type_, high-low+1);
	s.low_ = low;
	s.high_ = high;
	s.isSubSignal_ = true;
	
		return s;
}
	
Signal Signal::getException()
{
	if (!isFP_)
		throw string("Not a floating point signal.");
	return getSubSignal(1+wE_+wF_, 2+wE_+wF_);
}

Signal Signal::getSign()
{
	if (!isFP_)
			throw string("Not a floating point signal.");
	return getSubSignal(wE_+wF_, wE_+wF_);
}

Signal Signal::getExponent()
{
	if (!isFP_)
		throw string("Not a floating point signal.");
	return getSubSignal(wF_, wE_+wF_-1);
}

Signal Signal::getMantissa()
{
	if (!isFP_)
		throw string("Not a floating point signal.");
	return getSubSignal(0, wF_-1);
}

void Signal::setTTL(uint32_t ttl) {
	ttl_ = ttl;
}

uint32_t Signal::getTTL() {
	return ttl_;
}

void Signal::setCycle(uint32_t cycle) {
	cycle_ = cycle;
}

uint32_t Signal::getCycle() {
	return cycle_;
}

void Signal::updateLifeSpan(uint32_t delay) {
	if(delay>lifeSpan_)
		lifeSpan_=delay;
	//TODO add TTL check
}

uint32_t Signal::getLifeSpan() {
	return lifeSpan_;
}

