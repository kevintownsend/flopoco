#ifndef OPERATOR_HPP
#define OPERATOR_HPP
#include <vector>
#include <map>
#include <gmpxx.h>
#include "Target.hpp"
#include "Signal.hpp"
#include "TestCase.hpp"
#include "TestIOMap.hpp"

using namespace std;

// variables set by the command-line interface in main.cpp

extern  int verbose;
const std::string tab = "   ";

/**
 * This is a top-level class representing an Operator.
 * This class is inherited by all classes which will output a VHDL entity.
 */
class Operator
{
public:
	/** Default constructor.
	 * Creates a default operator instance. No parameters are passed to this constructor.
	 */
	Operator(){
		numberOfInputs_             = 0;
		numberOfOutputs_            = 0;
		hasRegistersWithoutReset_   = false;
		hasRegistersWithAsyncReset_ = false;
		hasRegistersWithSyncReset_  = false;
		pipelineDepth_              = 0;
		cycle_                      = 0;
	}
	
	/** Operator Constructor.
	 * Creates an operator instance with an instantiated target for deployment.
	 * @param target_ The deployment target of the operator.
	 */
	Operator(Target* target)  {
		target_                     = target;
		numberOfInputs_             = 0;
		numberOfOutputs_            = 0;
		hasRegistersWithoutReset_   = false;
		hasRegistersWithAsyncReset_ = false;
		hasRegistersWithSyncReset_  = false;
		pipelineDepth_              = 0;
		cycle_                      = 0;
	}

	/** Operator Destructor.
	 */	
	virtual ~Operator() {}


	/** Adds an input signal to the operator.
	 * Adds a signal of type Signal::in to the the I/O signal list.
	 * @param name  the name of the signal
	 * @param width the number of bits of the signal.
	 */
	void addInput  (const std::string name, const int width=1);
	

	/** Adds an output signal to the operator.
	 * Adds a signal of type Signal::out to the the I/O signal list.
	 * @param name  the name of the signal
	 * @param width the number of bits of the signal.
	 */
	void addOutput(const std::string name, const int width=1);
	

	/** Adds a floating point input signal to the operator.
	 * Adds a signal of type Signal::in to the the I/O signal list, 
	 * having the FP flag set on true. The total width of this signal will
	 * be wE + wF + 3. (2 bits for exception, 1 for sign)
	 * @param name the name of the signal
	 * @param wE   the width of the exponent
	 * @param wF   the withh of the fraction
	 */
	void addFPInput(const std::string name, const int wE, const int wF);


	/** Adds a floating point output signal to the operator.
	 * Adds a signal of type Signal::out to the the I/O signal list, 
	 * having the FP flag set on true. The total width of this signal will
	 * be wE + wF + 3. (2 bits for exception, 1 for sign)
	 * @param name the name of the signal
	 * @param wE   the width of the exponent
	 * @param wF   the withh of the fraction
	 */	
	void addFPOutput(const std::string name, const int wE, const int wF);

	
	/** Generic function that adds a signal to the signal list.
	 * The functions add*Signal* are implemented as instances of this functions.
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 * @param delay the delay of the signal. The number of register levels that this signal needs to be delayed. If negative, treated as zero (useful for pipeline synchronization)
	 * @param regType the register type (registeredWithoutReset,registeredWithAsyncReset,registeredWithSyncReset, see Signal::SignalType)
	 * @param isBus if true, a signal with a width of 1 is declared as std_logic_vector; if false it is declared as std_logic.
	 */	
	void addSignalGeneric(const string name, const int width, const int delay, Signal::SignalType regType, bool isBus);



	/** Adds a signal to the signal list.
	 * Adds a signal of type Signal::wire to the the signal list.
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 */	
	void addSignal(const std::string name, const int width=1);

	/** Adds a bus signal to the signal list.
	 * Adds a signal of type Signal::wire to the the signal list. 
	 * The signal added by this method has a flag indicating that it is a bus. 
	 * Even if the signal will have width=1, it will be declared as a standard_logic_vector(0 downto 0)
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 */	
	void addSignalBus(const std::string name, const int width=1);
	

	/** Adds a delayed signal (without reset) to the signal list.
	 * Adds to the the signal list
	 * one signal of type Signal::wire
	 * and max(0,delay) signals of type Signal::registeredWithoutReset. 
	 * The signal names are: name, name_d, name_d_d .. and so on. 
	 * This method is equivalent to addSignal if delay<=0
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 * @param delay the delay of the signal. The number of register levels that this signal needs to be delayed. If negative, treated as zero (useful for pipeline synchronization)
	 */	
	void addDelaySignal(const std::string name, const int width, const int delay=1);


	/** Adds a delayed signal with synchronous reset to the signal list.
	 * Adds to the the signal list
	 * one signal of type Signal::wire
	 * and max(0,delay) signals of type Signal::registeredWithSyncReset. 
	 * The signal names are: name, name_d, name_d_d .. and so on. 
	 * This method is equivalent to addSignal if delay<=0
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 * @param delay the delay of the signal. The number of register levels that this signal needs to be delayed. If negative, treated as zero (useful for pipeline synchronization)
	 */	
	void addDelaySignalSyncReset(const std::string name, const int width, const int delay=1);


	/** Adds a delayed signal bus with synchronous reset to the signal list.
	 * Adds to the the signal list
	 * one signal of type Signal::wire
	 * and max(0,delay) signals of type Signal::registeredWithSyncReset. 
	 * The signal names are: name, name_d, name_d_d .. and so on. 
	 * Each signal added by this method has a flag indicating that it is a bus. 
	 * Even if the signal will have width=1, it will be declared as a standard_logic_vector(0 downto 0)
	 * This method is equivalent to addSignalBus if delay<=0
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 * @param delay the delay of the signal. The number of register levels that this signal needs to be delayed. If negative, treated as zero (useful for pipeline synchronization)
	 */	
	void addDelaySignalBus(const std::string name, const int width, const int delay=1);


	/** Adds a delayed signal bus with synchronous reset to the signal list.
	 * Adds to the the signal list
	 * one signal of type Signal::wire
	 * and max(0,delay) signals of type Signal::registeredWithSyncReset. 
	 * The signal names are: name, name_d, name_d_d .. and so on. 
	 * Each signal added by this method has a flag indicating that it is a bus. 
	 * Even if the signal will have width=1, it will be declared as a standard_logic_vector(0 downto 0)
	 * This method is equivalent to addSignalBus if delay<=0
	 * @param name  the name of the signal
	 * @param width the width of the signal
	 * @param delay the delay of the signal. The number of register levels that this signal needs to be delayed. If negative, treated as zero (useful for pipeline synchronization)
	 */	
	void addDelaySignalBusSyncReset(const std::string name, const int width, const int delay=1);
	
	
	/** Returns the name of a signal at a certain delay.
	 * Returns a string of the form name_d_d_d..._d where #(_d)=delay
	 * @param name  the name of the signal
	 * @param delay the delay of the signal. The number of register levels that this signal needs to be delayed with.
	 If delay<=0 then no delay is inserted. This way the following code synchronizes signals from two paths with different delays:
	 getDelaySignalName(s1,  s2_delay - s1_delay)  // will be delayed if s2_delay>s1_delay
    getDelaySignalName(s2,  s1_delay - s2_delay)  // will be delayed if s1_delay>s2_delay

	 If the operator is not sequential, the string returned is simply
	 name. In principle, the code for a sequential operator is thus
	 gracefully degraded into combinatorial code. See FPAdder for an example.
	 */	 
	string  delaySignal(const string name, const int delay=1);


	/** Checks that each delayed signal is indeed used.
	*/
	void checkDelays();

	/** Sets Operator name to default name.
	 * This method must be overridden by all classes which extend Operator
	 *  
	*/
	virtual void setOperatorName(){ uniqueName_ = "UnknownOperator";};
	
	/** Sets Operator name to givenName.
	 * Sets the name of the operator to operatorName.
	 * @param operatorName new name of the operator
	*/
	void setOperatorName(std::string operatorName);
	
	/** Sets Operator name to prefix_(uniqueName_)_postfix
	 * @param prefix the prefix string which will be palced in front of the operator name
	 *               formed with the operator internal parameters
	 * @param postfix the postfix string which will be palced at the end of the operator name
	 *                formed with the operator internal parameters
	*/
	void setOperatorName(std::string prefix, std::string postfix);
	
	/** Sets the type of the operator. 
	 * The information is retrived from the deployment target 
	 */
	void setOperatorType();
	
	/** Sets the commented name of an operator. 
	 * This is the name which will be shown in the comment of the operator.
	 * @param name the name which will apear in the comment of the operator.
	 */
	void setCommentedName(std::string name);
	
	/** Return the operator name. 
	 * Returns a string value representing the name of the operator. 
	 * @return operator name
	 */
	string getOperatorName() const;
	
	/** Return the number of input+output signals 
	 * @return the size of the IO list. The total number of input and output signals
	 *         of the architecture.
	 */
	int getIOListSize() const;
	
	/** Returns a pointer to the list containing the IO signals.
	 * @return pointer to ioList 
	 */
	vector<Signal*> * getIOList();
	
	/** Returns a pointer a signal from the ioList.
	 * @param the index of the signal in the list
	 * @return pointer to the i'th signal of ioList 
	 */
	const Signal * getIOListSignal(int i);
		
	/** Returns a pointer to the signal having the name s. Throws an exception if the signal is not yet declared.
	  * @param s then name of the signal we want to return
	  * @return the pointer to the signal having name s 
	  */
	Signal* getSignalByName(string s);

	/** Outputs component declaration 
	 * @param o the stream where the component is outputed
	 * @param name the name of the VHDL component we want to output to o
	 */
	virtual void outputVHDLComponent(std::ostream& o, std::string name);
	
	/** Outputs the VHDL component code of the current operator
	 * @param o the stream where the component is outputed
	 */
	void outputVHDLComponent(std::ostream& o);  
	
	/** Function which outputs the processes which declare the registers ( assign name_d <= name )
	 * @param o the stream where the component is outputed
	 */
	void outputVHDLRegisters(std::ostream& o);
	
	/** Output the licence
	 * @param o the stream where the licence is going to be outputted
	 * @param authorsYears the names of the authors and the years of their contributions
	 */
	void licence(std::ostream& o, std::string authorsYears);
		
	/** Output the standard library paperwork 
	 * @param o the stream where the libraries will be written to
	 */
	static void stdLibs(std::ostream& o){
		o<<"library ieee;\nuse ieee.std_logic_1164.all;"<<endl 
		 <<"use ieee.std_logic_arith.all;"<<endl
		 <<"use ieee.std_logic_unsigned.all;"<<endl 
		 <<"library work;"<<endl<<endl;
	};
		
	/** Output the VHDL entity of the current operator.
	 * @param o the stream where the entity will be outputted
	 */
	void outputVHDLEntity(std::ostream& o);
	
	/** output all the signal declarations 
	 * @param o the stream where the signal deca
	 */
	void outputVHDLSignalDeclarations(std::ostream& o);

	/**
	 * A new line inline function
	 * @param[in,out] o the stream to which the new line will be added
	 **/
	inline void newLine(std::ostream& o) {	o<<endl; }
	
	/**
	 * A new architecture inline function
	 * @param[in,out] o 	- the stream to which the new architecture line will be added
	 * @param[in]     name	- the name of the entity corresponding to this architecture
	 **/
	inline void newArchitecture(std::ostream& o, std::string name){
		o << "architecture arch of " << name  << " is" << endl;
	}
	
	/**
	 * A begin architecture inline function 
	 * @param[in,out] o 	- the stream to which the begin line will be added
	 **/
	inline void beginArchitecture(std::ostream& o){
		o << "begin" << endl;
	}

	/**
	 * A end architecture inline function 
	 * @param[in,out] o 	- the stream to which the begin line will be added
	 **/
	inline void endArchitecture(std::ostream& o){
		o << "end architecture;" << endl << endl;
	}

	/** the main function outputs the VHDL for the operator 
	 * @param o the stream where the entity will outputted
	 * @param name the name of the architecture
	 */
	virtual void outputVHDL(std::ostream& o, std::string name) =0 ;
	
	/** the main function outputs the VHDL for the operator 
	 * @param o the stream where the entity will outputted
	 */	
	void outputVHDL(std::ostream& o);   // calls the previous with name = uniqueName

	/** True if the operator needs a clock signal; 
	 * It will also get a rst but doesn't need to use it.
	 */	
	bool isSequential();  
	
	/** Set the operator to sequential
	 */	
	void setSequential(); 
	
	/** Set the operator to combinatorial
	 */	
	void setCombinatorial();
	
	/** Set the depth of the pipeline
	 * @param d the depth of the pipeline for the operator
	 */	
	void setPipelineDepth(int d);
	
	/** Gets the pipeline depth of this operator 
	 * @return the pipeline depth of the operator
	*/
	int getPipelineDepth();
	
	/** Increments the pipeline depth of the current operator 
	*/
	void incrementPipelineDepth();
	
	/**
	 * Gets the signals which are interesting for TestCases.
	 * @see TestIOMap
	 */
	virtual TestIOMap getTestIOMap() {
		throw std::string("getTestIOMap: not implemented for ") + uniqueName_;
	}

	/**
	 * Gets the correct value associated to one or more inputs.
	 * @param a the array which contains both already filled inputs and
	 *          to be filled outputs in the order specified in getTestIOMap.
	 */
	virtual void fillTestCase(mpz_class a[]) {
		throw std::string("fillTestCorrectValues: not implemented for ") + uniqueName_;
	}
		
	/** Final report function, prints to the terminal.  By default
	 * reports the pipeline depth, but feel free to overload if you have any
	 * thing useful to tell to the end user
	*/
	virtual void outputFinalReport();	
	
#if 0
	/** Emulate the operator, using mpfr. 
	 */	
	virtual void emulate();	
#endif 


	//////////////////From here on we have methods of FloPoCoPipelineFramework2.0
	// The main incompatibility with previous FloPoCo is that a delayed signal is now added only once to signalList
	// and all its delayed versions are generated by buildVHDLSignalDeclarations and buildVHDLRegisters

	
	/** Define the current cycle 
	 * @param the new value of the current cycle */
	void setCycle(int cycle) ;


	/** Define the current cycle 
	 * @param the new value of the current cycle */
	void nextCycle() ;


	/** Set the current cycle to that of a signal 
	 * @param name is the signal name. It must have been defined before */
	void syncCycleFromSignal(string name) ;


	/** Declares a signal implicitely by having it appearing on the Left Hand Side of a VHDL assignment
	 * @param name is the name of the signal
	 * @param width is the width of the signal (optional, default 1)
	 * @param isbus: a signal of width 1 is declared as std_logic when false, as std_logic_vector when true (optional, default false)
	 * @return name
	 */
	string declare(string name, const int width=1, bool isbus=false);

	/** use a signal on the Right 
	 * @param name is the name of the signal
	 * @return name
	 */
	string use(string name);

	
	/** Declare an output mapping for an instance of a sub-component
	 * Also declares the local signal implicitely, with width taken from the component 	
	 * @param op is a pointer to the subcomponent
	 * @param componentPortName is the name of the port on the component
	 * @param actualSignalName is the name of the signal in This mapped to this port
	 * @return name
	 */
	void outPortMap(Operator* op, string componentPortName, string actualSignalName);


	/** use a signal as input of a subcomponent
	 * @param componentPortName is the name of the port on the component
	 * @param actualSignalName is the name of the signal in This mapped to this port
	 * @return name
	 */
	void inPortMap(Operator* op, string componentPortName, string actualSignalName);


	/** returns the VHDL for an instance of a sub-component. 
	 * @param componentPortName is the name of the port on the component
	 * @param actualSignalName is the name of the signal in This mapped to this port
	 * @return name
	 */
	string instance(Operator* op, string instanceName);

	
	/** build all the signal declarations from signals implicitely declared by declare().
	 *  This is the 2.0 equivalent of outputVHDLSignalDeclarations
	 */
	string buildVHDLSignalDeclarations();


	/** build all the registers from signals implicitely delayed by declare() 
	 *	 This is the 2.0 equivalent of outputVHDLSignalRegisters
	 */
	string buildVHDLRegisters();



	// TODO: add methods that allow for signals with reset.

	//////////////////End of FloPoCoPipelineFramework2.0

protected:    
	Target*         target_;     /**< The target on which the operator will be deployed */
	string          uniqueName_; /**< By default, a name derived from the operator class and the parameters */
	vector<Signal*> ioList_;     /**< The list of I/O signals of the operator */
	vector<Signal*> signalList_; /**< The list of internal signals of the operator */
	map<string, string>  portMap_;/**< Port map for an instance of this operator */
	map<string, double>    outDelayMap;
	ostringstream     vhdl;      /**< The internal stream to which the constructor will build the VHDL code */
private:
	int                    numberOfInputs_;             /**< The number of inputs of the operator */
	int                    numberOfOutputs_;            /**< The number of outputs of the operator */
	bool                   isSequential_;               /**< True if the operator needs a clock signal*/
	int                    pipelineDepth_;              /**< The pipeline depth of the operator. 0 for combinatorial circuits */
	map<string, Signal*>   signalMap_;                  /**< A container of tuples for recovering the signal based on it's name */ 
	bool                   hasRegistersWithoutReset_;   /**< True if the operator has registers without a reset */
	bool                   hasRegistersWithAsyncReset_; /**< True if the operator has registers having an asynch reset */
	bool                   hasRegistersWithSyncReset_;  /**< True if the operator has registers having a synch reset */
	string                 commentedName_;              /**< Usually is the default name of the architecture.  */
	int                    cycle_;                      /**< The current cycle, when building a pipeline */

};
#endif
