#include <iostream>
#include <sstream>
#include "Signal.hpp"


Signal::Signal(const std::string name, const Signal::SignalType type, const int width, const bool isBus, const int ttl) : 
	name_(name), type_(type), wE_(0), wF_(0), width_(width), low_(0), high_(width-1),
	isSubSignal_(false), isBus_(isBus),	isFP_(false), ttl_(ttl) {
	updateSignalName();
}

Signal::Signal(const std::string name, const Signal::SignalType type, const int wE, const int wF) : 
	name_(name), type_(type), wE_(wE), wF_(wF), width_(wE+wF+3),
	low_(0), high_(width_-1), isSubSignal_(false),isBus_(false), isFP_(true), ttl_(0)
{
	updateSignalName();
}

Signal::~Signal(){}

const std::string& Signal::getSignalName() const { 
	return id_; 
}

void Signal::updateSignalName()	{
	if (isSubSignal_ == false)
		id_ = name_;
	else
		{
			std::stringstream o;
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
	
std::string Signal::toVHDL() {
	std::ostringstream o; 
	if(type()==Signal::wire || type()==Signal::registeredWithoutReset || type()==Signal::registeredWithAsyncReset || type()==Signal::registeredWithSyncReset) 
		o << "signal ";
	o << getSignalName();
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

Signal Signal::getSubSignal(int low, int high)
{
	if (low < low_)
		throw std::string("Attempted to return subsignal with smaller low index.");
	if (high > high_)
		throw std::string("Attempted to return subsignal with bigger high index."); 
	
	Signal s(name_, type_, high-low+1);
	s.low_ = low;
	s.high_ = high;
	s.isSubSignal_ = true;
	
		return s;
}
	
Signal Signal::getException()
{
	if (!isFP_)
		throw std::string("Not a floating point signal.");
	return getSubSignal(1+wE_+wF_, 2+wE_+wF_);
}

Signal Signal::getSign()
{
	if (!isFP_)
			throw std::string("Not a floating point signal.");
	return getSubSignal(wE_+wF_, wE_+wF_);
}

Signal Signal::getExponent()
{
	if (!isFP_)
		throw std::string("Not a floating point signal.");
	return getSubSignal(wF_, wE_+wF_-1);
}

Signal Signal::getMantissa()
{
	if (!isFP_)
		throw std::string("Not a floating point signal.");
	return getSubSignal(0, wF_-1);
}

void Signal::setTTL(uint32_t ttl) {
	ttl_ = ttl;
}

uint32_t Signal::getTTL() {
	return ttl_;
}

