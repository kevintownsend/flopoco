#ifndef __SIGNAL_HPP
#define __SIGNAL_HPP

#include <iostream>
#include <sstream>

#ifdef _WIN32
  #include "pstdint.h"
#else
  #include <inttypes.h>
#endif


/**
 * A class representing a signal. This is the basic block that operators use
 */
class Signal
{
public:
	/** The possible types of a signal*/
	typedef enum {
	in,                        /**< if the signal is an input signal */
	out,                       /**< if the signal is an output signal */
	wire,                      /**< if the signal is a wire (not registered) */
	registeredWithoutReset,    /**< if the signal is registered, but does not have a reset */
	registeredWithAsyncReset,  /**< if the signal is registered, and has an asynchronous reset */
	registeredWithSyncReset    /**< if the signal is registered, and has an synchronous reset */
	} SignalType;

	/** Signal constructor.
	 * The standard constructor for signals which are not floating-point.
	 * @param name      the name of the signal
	 * @param type      the type of the signal, @see SignalType
	 * @param width     the width of the signal
	 * @param isBus     the flag which signals if the signal is a bus (std_logic_vector)
	 */
	Signal(const std::string name, const SignalType type, const int width = 1, const bool isBus = false) : 
		name_(name), type_(type), wE_(0), wF_(0), width_(width), low_(0), high_(width-1),
		isSubSignal_(false), isBus_(isBus),	isFP_(false) {
		
		updateSignalName();
	}

	/** Signal constructor.
	 * The standard constructor for signals which are not floating-point.
	 * @param name      the name of the signal
	 * @param type      the type of the signal, @see SignalType
	 * @param width     the width of the signal
	 * @param isBus     the flag which signals if the signal is a bus (std_logic_vector)
	 */
	Signal(const std::string name, const SignalType type, const int wE, const int wF) : 
		name_(name), type_(type), wE_(wE), wF_(wF), width_(wE+wF+3),
		low_(0), high_(width_-1), isSubSignal_(false),isBus_(false), isFP_(true)
	{
		updateSignalName();
	}

	/** Signal destructor.
	 */		
	~Signal(){}
	
	/** Returns the name of the signal
	 * @return the name of the signal
	 */	
	const std::string& getSignalName() const { 
		return id_; 
	}

	/** Updates the name of the signal.
	 * It takes into consideration the fact that we might have subsignals
	 */	
	void updateSignalName()	{
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
	
	/** Returns the width of the signal
	 * @return the width of the signal
	 */	
	int width() const{return width_;}
	
	/** Returns the exponent width of the signal
	 * @return the width of the exponent if signal is isFP_
	 */	
	int wE() const {return(wE_);}

	/** Returns the fraction width of the signal
	 * @return the width of the fraction if signal is isFP_
	 */	
	int wF() const {return(wF_);}
	
	/** Reports if the signal is a floating-point signal
	 * @return if the signal is a FP siglal
	 */	
	bool isFP() const {return isFP_;}

	/** Reports if the signal has the bus flag active
	 * @return true if the signal is of bus type (std_logic_vector)
	 */		
	bool isBus() const {return isBus_;}

	/** Returns the type of the signal
	 * @return type of signal, @see SignalType
	 */	
	SignalType type() const {return type_;}
	
	/** outputs the VHDL code for declaring this signal 
	 * @return the VHDL for this signal. 
	 */	
	std::string toVHDL() {
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

	/** Returns a subsignal of the this signal
	 * @param low the low index of subsignal
	 * @param high the high index of the subsignal
	 * @return the corresponding subsignal
	 */	
	Signal getSubSignal(int low, int high)
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
	
	/** Returns a subsignal containing the exception bits of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getException()
	{
		if (!isFP_)
			throw std::string("Not a floating point signal.");
		return getSubSignal(1+wE_+wF_, 2+wE_+wF_);
	}

	/** Returns a subsignal containing the sign bit of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getSign()
	{
		if (!isFP_)
			throw std::string("Not a floating point signal.");
		return getSubSignal(wE_+wF_, wE_+wF_);
	}

	/** Returns a subsignal containing the exponent bits of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getExponent()
	{
		if (!isFP_)
			throw std::string("Not a floating point signal.");
		return getSubSignal(wF_, wE_+wF_-1);
	}

	/** Returns a subsignal containing the fraction bits of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getMantissa()
	{
		if (!isFP_)
			throw std::string("Not a floating point signal.");
		return getSubSignal(0, wF_-1);
	}

private:
	std::string   name_;        /**< The name of the signal */
	std::string   id_;          /**< The id of the signal. It is the same as name_ for regular signals, and is name_(high_-1 downto low_) for subsignals */
	SignalType    type_;        /**< The type of the signal, see SignalType */
	uint32_t      width_;       /**< The width of the signal */
	
	bool          isFP_;        /**< If the signal is of floating-point type */  
	uint32_t      wE_;          /**< The width of the exponent. Used for FP signals */
	uint32_t      wF_;          /**< The width of the fraction. Used for FP signals */
	
	bool          isSubSignal_; /**< If the signal is a subsignal */
	uint32_t      low_;         /**< The low index of the signal */
	uint32_t      high_;        /**< The high index of the signal */
	
	bool          isBus_;       /**< True is the signal is a bus (std_logic_vector)*/
};

#endif

