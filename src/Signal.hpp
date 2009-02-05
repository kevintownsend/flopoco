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
	Signal(const std::string name, const Signal::SignalType type, const int width = 1, const bool isBus=false);
	/** Signal constructor.
	 * The standard constructor for signals which are floating-point.
	 * @param name      the name of the signal
	 * @param type      the type of the signal, @see SignalType
	 * @param wE        the exponent width
	 * @param wF        the significand width
	 */
	Signal(const std::string name, const SignalType type, const int wE, const int wF);

	/** Signal destructor.
	 */		
	~Signal();
	
	/** Returns the name of the signal
	 * @return the name of the signal
	 */	
	const std::string& getName() const;

	/** Updates the name of the signal.
	 * It takes into consideration the fact that we might have subsignals
	 */	
	void updateSignalName();


	/** Returns the width of the signal
	 * @return the width of the signal
	 */	
	int width() const;

	
	/** Returns the exponent width of the signal
	 * @return the width of the exponent if signal is isFP_
	 */	
	int wE() const;

	/** Returns the fraction width of the signal
	 * @return the width of the fraction if signal is isFP_
	 */	
	int wF() const;
	
	/** Reports if the signal is a floating-point signal
	 * @return if the signal is a FP siglal
	 */	
	bool isFP() const;

	/** Reports if the signal has the bus flag active
	 * @return true if the signal is of bus type (std_logic_vector)
	 */		
	bool isBus() const;

	/** Returns the type of the signal
	 * @return type of signal, @see SignalType
	 */	
	SignalType type() const;
	
	/** outputs the VHDL code for declaring this signal 
	 * @return the VHDL for this signal. 
	 */	
	std::string toVHDL();

	/** obtain the name of a signal delayed by delay 
	 * @param delay*/
	std::string delayedName(int delay);


	/** outputs the VHDL code for declaring a signal with all its delayed versions
	 * This is the 2.0 equivalent of toVHDL()
	 * @return the VHDL for the declarations. 
	 */	
	std::string toVHDLDeclaration();

	/** Returns a subsignal of the this signal
	 * @param low the low index of subsignal
	 * @param high the high index of the subsignal
	 * @return the corresponding subsignal
	 */	
	Signal getSubSignal(int low, int high);
	
	/** Returns a subsignal containing the exception bits of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getException();

	/** Returns a subsignal containing the sign bit of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getSign();

	/** Returns a subsignal containing the exponent bits of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getExponent();


	/** Returns a subsignal containing the fraction bits of this signal
	 * @return the corresponding subsignal
	 */	
	Signal getMantissa();


	/** sets the time to live value (consider using the constructor instead)
	 * @return the corresponding subsignal
	 */	
	void setTTL(uint32_t ttl);


	/** obtain the time to live value of this signal
	 * @return the TTL value
	 */	
	uint32_t getTTL();


	/** sets the cycle at which the signal is active
	 */	
	void setCycle(uint32_t cycle) ;


	/** obtain the declared cycle of this signal
	 * @return the cycle
	 */	
	uint32_t getCycle();


	/** Updates the max delay associated to a signal
	 */	
	void updateLifeSpan(uint32_t delay) ;


	/** obtain max delay that has been applied to this signal
	 * @return the max delay 
	 */	
	uint32_t getLifeSpan() ;

private:
	std::string   name_;        /**< The name of the signal */
	std::string   id_;          /**< The id of the signal. It is the same as name_ for regular signals, and is name_(high_-1 downto low_) for subsignals */
	SignalType    type_;        /**< The type of the signal, see SignalType */
	uint32_t      width_;       /**< The width of the signal */

	uint32_t      ttl_;         /**< Time To Live for a delayed signal; used for early error reporting in delaySignal() */
	uint32_t      lifeSpan_;    /**< The max delay that will be applied to this signal; set by delaySignal(), should be eventually equal to ttl_, otherwise a warning may be produced */
	uint32_t      cycle_;       /**<  the cycle at which this signal is active in a pipelined operator. 0 means synchronized with the inputs */
	
	bool          isFP_;        /**< If the signal is of floating-point type */  
	uint32_t      wE_;          /**< The width of the exponent. Used for FP signals */
	uint32_t      wF_;          /**< The width of the fraction. Used for FP signals */
	
	bool          isSubSignal_; /**< If the signal is a subsignal */
	uint32_t      low_;         /**< The low index of the signal */
	uint32_t      high_;        /**< The high index of the signal */
	
	bool          isBus_;       /**< True is the signal is a bus (std_logic_vector)*/
};

#endif

