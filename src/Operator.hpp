#ifndef OPERATOR_HPP
#define OPERATOR_HPP
#include <vector>
#include <map>
#include <gmpxx.h>
#include "Target.hpp"
#include "Signal.hpp"
#include "TestCase.hpp"
#include <float.h>
#include <utility>
#include <vector>
#include "FlopocoStream.hpp"


using namespace std;

// variables set by the command-line interface in main.cpp


namespace flopoco {

	// Reporting levels
#define LIST 0
#define INFO 1
#define DETAILED 2
#define DEBUG 3
#define FULL 4


#define INNER_SEPARATOR "................................................................................"
#define DEBUG_SEPARATOR "________________________________________________________________________________"
#define OUTER_SEPARATOR "################################################################################"
#define REPORT(level, stream) { if ((level)<=(verbose)) cerr << "> " << srcFileName << ": " << stream << endl;} 



extern int verbose;
//extern bool combinatorialOperator; //FIXME: BAD design practice working with global variables



const std::string tab = "   ";
extern map<string, double> emptyDelayMap;





/**
 * This is a top-level class representing an Operator.
 * This class is inherited by all classes which will output a VHDL entity.
 */
class Operator
{

	static int uid;                  /**< The counter holding a unique id */

public:

	static int getNewUId(){
		Operator::uid++;
		return Operator::uid;
	}

	/** Operator Constructor.
	 * Creates an operator instance with an instantiated target for deployment.
	 * @param target_ The deployment target of the operator.
	 */
	Operator(Target* target, map<string, double> inputDelays = emptyDelayMap)  {
		target_                     = target;
		numberOfInputs_             = 0;
		numberOfOutputs_            = 0;
		hasRegistersWithoutReset_   = false;
		hasRegistersWithAsyncReset_ = false;
		hasRegistersWithSyncReset_  = false;
		pipelineDepth_              = 0;
		currentCycle_               = 0;
		needRecirculationSignal_    = false;
		inputDelayMap               = inputDelays;

		if (target_->isPipelined())
			setSequential();
		else
			setCombinatorial();	
			
		vhdl.disableParsing(!target_->isPipelined());	
	}

	/** Operator Destructor.
	 */	
	virtual ~Operator() {}


 /*****************************************************************************/
 /*         Paperwork-related methods (for defining an operator entity)       */
 /*****************************************************************************/

	/** Adds an input signal to the operator.
	 * Adds a signal of type Signal::in to the the I/O signal list.
	 * @param name  the name of the signal
	 * @param width the number of bits of the signal.
	 * @param isBus describes if this signal is a bus, that is, an instance of std_logic_vector
	 */
	void addInput  (const std::string name, const int width=1, const bool isBus=false);
	

	/** Adds an output signal to the operator.
	 * Adds a signal of type Signal::out to the the I/O signal list.
	 * @param name  the name of the signal
	 * @param width the number of bits of the signal.
	 * @param numberOfPossibleOutputValues (optional, defaults to 1) set to 2 for a faithfully rounded operator for instance
	 * @param isBus describes if this signal is a bus, that is, an instance of std_logic_vector
	 */	
	void addOutput(const std::string name, const int width=1, const int numberOfPossibleOutputValues=1, const bool isBus=false);
	

	/** Adds a floating point (FloPoCo format) input signal to the operator.
	 * Adds a signal of type Signal::in to the the I/O signal list, 
	 * having the FP flag set on true. The total width of this signal will
	 * be wE + wF + 3. (2 bits for exception, 1 for sign)
	 * @param name the name of the signal
	 * @param wE   the width of the exponent
	 * @param wF   the withh of the fraction
	 */
	void addFPInput(const std::string name, const int wE, const int wF);


	/** Adds a floating point (FloPoCo format) output signal to the operator.
	 * Adds a signal of type Signal::out to the the I/O signal list, 
	 * having the FP flag set on true. The total width of this signal will
	 * be wE + wF + 3. (2 bits for exception, 1 for sign)
	 * @param name the name of the signal
	 * @param wE   the width of the exponent
	 * @param wF   the withh of the fraction
	 * @param numberOfPossibleOutputValues (optional, defaults to 1) set to 2 for a faithfully rounded operator for instance
	 */	
	void addFPOutput(const std::string name, const int wE, const int wF, const int numberOfPossibleOutputValues=1);

	/** Adds a IEEE floating point input signal to the operator.
	 * Adds a signal of type Signal::in to the the I/O signal list, 
	 * having the FP flag set on true. The total width of this signal will
	 * be wE + wF + 1.  (1 bit for sign)
	 * @param name the name of the signal
	 * @param wE   the width of the exponent
	 * @param wF   the withh of the fraction
	 */
	void addIEEEInput(const std::string name, const int wE, const int wF);


	/** Adds a floating point output signal to the operator.
	 * Adds a signal of type Signal::out to the the I/O signal list, 
	 * having the FP flag set on true. The total width of this signal will
	 * be wE + wF + 1. (1 bit for sign)
	 * @param name the name of the signal
	 * @param wE   the width of the exponent
	 * @param wF   the withh of the fraction
	 * @param numberOfPossibleOutputValues (optional, defaults to 1) set to 2 for a faithfully rounded operator for instance
	 */	
	void addIEEEOutput(const std::string name, const int wE, const int wF, const int numberOfPossibleOutputValues=1);

	

	/** sets the copyright string, should be authors + year
	 * @param authorsYears the names of the authors and the years of their contributions
	 */
	void setCopyrightString(std::string authorsYears);


	/** Sets Operator name to givenName.
	 * Sets the name of the operator to operatorName.
	 * @param operatorName new name of the operator
	*/
	void setName(std::string operatorName = "UnknownOperator");
	
	/** This method should be used by an operator to change the default name of a sub-component. The default name becomes the commented name.
	 * @param operatorName new name of the operator
	*/
	void changeName(std::string operatorName);
	
	/** Sets Operator name to prefix_(uniqueName_)_postfix
	 * @param prefix the prefix string which will be placed in front of the operator name
	 *               formed with the operator internal parameters
	 * @param postfix the postfix string which will be placed at the end of the operator name
	 *                formed with the operator internal parameters
	*/
	void setName(std::string prefix, std::string postfix);
		
	/** Return the operator name. 
	 * Returns a string value representing the name of the operator. 
	 * @return operator name
	 */
	string getName() const;




 /*****************************************************************************/
 /*        VHDL-related methods (for defining an operator architecture)       */
 /*****************************************************************************/



	/* Functions related to pipeline management */

	/** Define the current cycle 
	 * @param the new value of the current cycle */
	void setCycle(int cycle, bool report=true) ;

	/** Return the current cycle 
	 * @return the current cycle */
	int getCurrentCycle(); 

	/** Define the current cycle 
	 * @param the new value of the current cycle */
	void nextCycle(bool report=true) ;

	/** get the critical path of the current cycle so far */
	double getCriticalPath() ;

	/** Set or reset the critical path of the current cycle  */
	void setCriticalPath(double delay) ;

	/** Adds to the critical path of the current stage, and insert a pipeline stage if needed
	 * @param the delay to add to the critical path of current pipeline stage */
	void addToCriticalPath(double delay) ;

	bool manageCriticalPath(double delay=0.0);

	/** Adds to the critical path of the current stage, and insert a pipeline stage if needed
	 * @param the delay to add to the critical path of current pipeline stage */
	void nextCycleCond(double delay) ;

	/** get the critical path delay associated to a given output of the operator
	 * @param the name of the output */
	double getOutputDelay(string s); 

	/** Set the current cycle to that of a signal. It may increase or decrease current cycle. 
	 * @param name is the signal name. It must have been defined before 
	 * @param report is a boolean, if true it will report the cycle 
	 */
	void setCycleFromSignal(string name, bool report=false) ;

	/** advance the current cycle to that of a signal. It may only increase current cycle. To synchronize
		 two or more signals, first call setCycleFromSignal() on the
		 first, then syncCycleFromSignal() on the remaining ones. It
		 will synchronize to the latest of these signals.  
		 * @param name is the signal name. It must have been defined before 
		 * @param report is a boolean, if true it will report the cycle 
		 */
	bool syncCycleFromSignal(string name, bool report=true) ;




	/** Declares a signal implicitely by having it appearing on the Left Hand Side of a VHDL assignment
	 * @param name is the name of the signal
	 * @param width is the width of the signal (optional, default 1)
	 * @param isbus: a signal of width 1 is declared as std_logic when false, as std_logic_vector when true (optional, default false)
	 * @param regType: the registring type of this signal. See also the Signal Class for more info
	 * @return name
	 */
	string declare(string name, const int width=1, bool isbus=false, Signal::SignalType regType = Signal::wire );
	// TODO: add methods that allow for signals with reset (when rewriting LongAcc for the new framework)


#if 1
	/** use a signal on the Right 
	 * @param name is the name of the signal
	 * @return name
	 */
	string use(string name);

	string use(string name, int delay);
#endif

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
	 * @param actualSignalName is the name of the signal (of this) mapped to this port
	 */
	void inPortMap(Operator* op, string componentPortName, string actualSignalName);

	/** use a constant signal as input of a subcomponent. 
	 * @param componentPortName is the name of the port on the component
	 * @param actualSignal is the constant signal to be mapped to this port
	 */
	void inPortMapCst(Operator* op, string componentPortName, string actualSignal);

	/** returns the VHDL for an instance of a sub-component. 
	 * @param op represents the operator to be port mapped 
	 * @param instanceName is the name of the instance as a label
	 * @return name
	 */
	string instance(Operator* op, string instanceName);

	/** returns the VHDL for an instance of a sub-component. 
	 * @param op represents the operator to be port mapped 
	 * @param instanceName is the name of the instance as a label
	 * @param clkName the name of the input clock signal
 	 * @param clkName the name of the input rst signal
	 * @return name
	 */
	string instance(Operator* op, string instanceName, string clkName, string rstName);

	/** adds attributes to the generated VHDL so that the tools use embedded RAM blocks for an instance
	 * @ param t a pointer to this instance
	 */
	void useHardRAM(Operator* t); 
	void useSoftRAM(Operator* t); 


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




 /*****************************************************************************/
 /*        Testing-related methods (for defining an operator testbench)       */
 /*****************************************************************************/

	/**
	 * Gets the correct value associated to one or more inputs.
	 * @param tc the test case, filled with the input values, to be filled with the output values.
	 * @see FPAdder for an example implementation
	 */
	virtual void emulate(TestCase * tc);
		
	/**
	 * Append standard test cases to a test case list. Standard test
	 * cases are operator-dependent and should include any specific
	 * corner cases you may think of. Never mind removing a standard test case because you think it is no longer useful!
	 * @param tcl a TestCaseList
	 */
	virtual void buildStandardTestCases(TestCaseList* tcl);


	/**
	 * Append random test cases to a test case list. There is a default
	 * implementation using a uniform random generator, but most
	 * operators are not exercised efficiently using such a
	 * generator. For instance, in FPAdder, the random number generator
	 * should be biased to favor exponents which are relatively close
	 * so that an effective addition takes place.
	 * @param tcl a TestCaseList
	 * @param n the number of random test cases to add
	 */
	virtual void buildRandomTestCaseList(TestCaseList* tcl, int n);

	/**
	 * Generate Random Test case identified by an integer . There is a default
	 * implementation using a uniform random generator, but most
	 * operators are not exercised efficiently using such a
	 * generator. For instance, in FPAdder, the random number generator
	 * should be biased to favor exponents which are relatively close
	 * so that an effective addition takes place.
	 * This function create a new TestCase (to be free after use)
	 * @param i the identifier of the test case to be generated
	 * @return TestCase*
	 */
	virtual TestCase* buildRandomTestCase(int i);



	Target* target() {
		return target_;
	}



 /*****************************************************************************/
 /*     From this point, we have methods that are not needed in normal use    */
 /*****************************************************************************/



	

	/** build all the signal declarations from signals implicitely declared by declare().
	 */
	string buildVHDLSignalDeclarations();

	/** build all the component declarations from the list built by instance().
	 */
	string buildVHDLComponentDeclarations();

	/** build all the registers from signals implicitely delayed by declare() 
	 *	 This is the 2.0 equivalent of outputVHDLSignalRegisters
	 */
	string buildVHDLRegisters();

	/** build all the type declarations.
	 */
	string buildVHDLTypeDeclarations();

	/** output the VHDL constants. */
	string buildVHDLConstantDeclarations();

	/** output the VHDL constants. */
	string buildVHDLAttributes();





	/** the main function outputs the VHDL for the operator.
		 If you use the modern (post-0.10) framework you no longer need to overload this method,
		 the default will do.
	 * @param o the stream where the entity will outputted
	 * @param name the name of the architecture
	 */
	virtual void outputVHDL(std::ostream& o, std::string name);
	
	/** the main function outputs the VHDL for the operator 
	 * @param o the stream where the entity will outputted
	 */	
	void outputVHDL(std::ostream& o);   // calls the previous with name = uniqueName





	/** True if the operator needs a clock signal; 
	 * It will also get a rst but doesn't need to use it.
	 */	
	bool isSequential();  


        /** True if the operator need a recirculation signal 
         *  TODO : change name
         */
        bool isRecirculatory();
	
	/** Set the operator to sequential.
		 You shouldn't need to use this method for standard operators 
		 (Operator::Operator()  calls it according to Target)
	 */	
	void setSequential(); 
	
	/** Set the operator to combinatorial
		 You shouldn't need to use this method for standard operators 
		 (Operator::Operator()  calls it according to Target)
	 */	
	void setCombinatorial();



        /** Set the operator to need a recirculation signal in order to 
                  trigger the pipeline work
         */
        void setRecirculationSignal();
	





	
	

	/** Returns a pointer to the signal having the name s. Throws an exception if the signal is not yet declared.
	  * @param s then name of the signal we want to return
	  * @return the pointer to the signal having name s 
	  */
	Signal* getSignalByName(string s);






	/** DEPRECATED Outputs component declaration 
	 * @param o the stream where the component is outputed
	 * @param name the name of the VHDL component we want to output to o
	 */
	virtual void outputVHDLComponent(std::ostream& o, std::string name);
	
	/**  DEPRECATED Outputs the VHDL component code of the current operator
	 * @param o the stream where the component is outputed
	 */
	void outputVHDLComponent(std::ostream& o);  
	
	/**  DEPRECATED Function which outputs the processes which declare the registers ( assign name_d <= name )
	 * @param o the stream where the component is outputed
	 */
	void outputVHDLRegisters(std::ostream& o);
	
		

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
	Signal * getIOListSignal(int i);
		
	



	/** DEPRECATED, better use setCopyrightString
		 Output the licence
	 * @param o the stream where the licence is going to be outputted
	 * @param authorsYears the names of the authors and the years of their contributions
	 */
	void licence(std::ostream& o, std::string authorsYears);


	/**  Output the licence, using copyrightString_
	 * @param o the stream where the licence is going to be outputted
	 */
	void licence(std::ostream& o);


	void pipelineInfo(std::ostream& o, std::string authorsYears);


	void pipelineInfo(std::ostream& o);

	/** DEPRECATED  Output the standard library paperwork 
	 * @param o the stream where the libraries will be written to
	 */
	static void stdLibs(std::ostream& o){
		o<<"library ieee;\nuse ieee.std_logic_1164.all;"<<endl 
		 <<"use ieee.std_logic_arith.all;"<<endl
		 <<"use ieee.std_logic_unsigned.all;"<<endl 
		 <<"library std;" << endl
		 <<"use std.textio.all;"<< endl 
		 <<"library work;"<<endl<< endl;
	};
		
	/** DEPRECATED  Output the VHDL entity of the current operator.
	 * @param o the stream where the entity will be outputted
	 */
	void outputVHDLEntity(std::ostream& o);
	
	/** DEPRECATED  output all the signal declarations 
	 * @param o the stream where the signal deca
	 */
	void outputVHDLSignalDeclarations(std::ostream& o);


	/** Add a VHDL type declaration. */
 	void addType(std::string name, std::string def);

	/** Add a VHDL constant. This may make the code easier to read, but more difficult to debug. */
	void addConstant(std::string name, std::string ctype, int cvalue);

	void addConstant(std::string name, std::string ctype, mpz_class cvalue);
	
	void addConstant(std::string name, std::string ctype, string cvalue);
	

	/** Add attribute, declaring the attribute name if it is not done already.
	 */ 
	void addAttribute(std::string attributeName,  std::string attributeType,  std::string object, std::string value );

	/**
	 * A new line inline function
	 * @param[in,out] o the stream to which the new line will be added
	 **/
	inline void newLine(std::ostream& o) {	o<<endl; }



	/** Final report function, prints to the terminal.  By default
	 * reports the pipeline depth, but feel free to overload if you have any
	 * thing useful to tell to the end user
	*/
	virtual void outputFinalReport();	
	
	
	/** Gets the pipeline depth of this operator 
	 * @return the pipeline depth of the operator
	*/
	int getPipelineDepth();




	/**
	* @return the output map containing the signal -> delay associations 
	*/	
	map<string, double> getOutDelayMap();
	
	/**
	* @return the output map containing the signal -> declaration cycle 
	*/	
	map<string, int> getDeclareTable();

//	/**
//	* @return the output map containing the signal -> declaration cycle 
//	*/	
//	vector<pair<string, int> > getUseTable(){
//		return useTable;
//	}
	
	FlopocoStream* getFlopocoVHDLStream(){
		return &vhdl;
	}

	void parse2();

	Target* getTarget(){
		return target_;
	}
	
	bool hasComponent(string s);
	

	int level; //printing issues

protected:    
	Target*             target_;          /**< The target on which the operator will be deployed */
	string              uniqueName_;      /**< By default, a name derived from the operator class and the parameters */
	vector<Signal*>     ioList_;          /**< The list of I/O signals of the operator */
	vector<Signal*>     testCaseSignals_; /**< The list of pointers to the signals in a test case entry. Its size also gives the dimension of a test case */
	vector<Signal*>     signalList_;      /**< The list of internal signals of the operator */
	map<string, Operator*> subComponents_;/**< The list of sub-components */
	map<string, string> portMap_;         /**< Port map for an instance of this operator */
	map<string, double> outDelayMap;      /**< Slack delays on the outputs */
	map<string, double> inputDelayMap;       /**< Slack delays on the inputs */
	FlopocoStream       vhdl;             /**< The internal stream to which the constructor will build the VHDL code */
	string              srcFileName;      /**< Used to debug and report.  */
	map<string, int>    declareTable;     /**< Table containing the name and declaration cycle of the signal */


private:
	int                    numberOfInputs_;             /**< The number of inputs of the operator */
	int                    numberOfOutputs_;            /**< The number of outputs of the operator */
	bool                   isSequential_;               /**< True if the operator needs a clock signal*/
	int                    pipelineDepth_;              /**< The pipeline depth of the operator. 0 for combinatorial circuits */
	map<string, Signal*>   signalMap_;                  /**< A container of tuples for recovering the signal based on it's name */ 
	map<string, pair<string, string> > constants_;      /**< The list of constants of the operator: name, <type, value> */
	map<string, string>    attributes_;                  /**< The list of attribute declarations (name, type) */
	map<string, string>    types_;                      /**< The list of type declarations (name, type) */
	map<pair<string,string>, string >  attributesValues_;/**< attribute values <attribute name, object (component, signal, etc)> ,  value> */
	bool                   hasRegistersWithoutReset_;   /**< True if the operator has registers without a reset */
	bool                   hasRegistersWithAsyncReset_; /**< True if the operator has registers having an asynch reset */
	bool                   hasRegistersWithSyncReset_;  /**< True if the operator has registers having a synch reset */
	string                 commentedName_;              /**< Usually is the default name of the architecture.  */
	string                 copyrightString_;            /**< Authors and years.  */
	int                    currentCycle_;               /**< The current cycle, when building a pipeline */
	double                 criticalPath_;               /**< The current delay of the current pipeline stage */
	bool                   needRecirculationSignal_;    /**< True if the operator has registers having a recirculation signal  */

};

} //namespace
#endif
