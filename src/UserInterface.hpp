/**
A generic user interface class that manages the command line and the documentation for various operators
Includes an operator factory inspired by David Thomas

For typical use, see src/ExpLog/FPExp.*

Author : David Thomas, Florent de Dinechin

Initial software.
Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL,  
2015-.
  All rights reserved.

*/
#ifndef UserInterface_hpp
#define UserInterface_hpp

#include "Operator.hpp"
#include <memory>

// Operator Factory, based on the one by David Thomas, with a bit of clean up.
// For typical use, see src/ShiftersEtc/Shifter  or   src/ExpLog/FPExp

namespace flopoco
{
	
		// Note: not using boost::function here, as it's likely to scare people, and also drags in quite a few header dependencies
	typedef OperatorPtr (*parser_func_t)(Target *,const vector<string> &);	
	class OperatorFactory;
	typedef shared_ptr<OperatorFactory> OperatorFactoryPtr;






	
	/** This is the class that manages a list of OperatorFactories, and the overall command line and documentation.
			Each OperatorFactory is responsible for the command line and parsing for one Operator sub-class. */
	class UserInterface
	{
	public:
		static void registerFactory(OperatorFactoryPtr factory);
		/**a helper factory function*/ 
		static void add(
										string name,
										string description, /**< for the HTML doc and the detailed help */ 
										string categories,	/**<  semicolon-seperated list of categories */
										string parameterList, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
										parser_func_t parser	);
		
		static unsigned getFactoryCount();
		static OperatorFactoryPtr getFactoryByIndex(unsigned i);
		static OperatorFactoryPtr getFactoryByName(string operatorName);


		////////////////// Parsing-related ///////////////////////////////
		static void parseGlobalOptions(const vector<string> &args);

		static void parseAll(Target* target, int argc, char* argv[]);
		//		static param_map_t  parseArguments(string opName, const vector<string> &args);
		static int checkStrictlyPositiveInt(vector<string>, string);
		static int checkOptionalInt(vector<string>, string);

		/** Provide a string with the full documentation. TODO: an HTML version*/
		static string getFullDoc();
		


		/** add an operator to the global (first-level) list, which is stored in its Target (not really its place, sorry).
				This method should be called by 
				1/ the main / top-level, or  
				2/ for sub-components that are really basic operators, 
				expected to be used several times, *in a way that is independent of the context/timing*.
				Typical example is a table designed to fit in a LUT or parallel row of LUTs
		*/

		static void addToGlobalOpList(OperatorPtr op);

		

		/** generates the code for operators in globalOpList, and all their subcomponents */
		static void outputVHDLToFile(ofstream& file);
		
		/** generates a report for operators in globalOpList, and all their subcomponents */
		static void finalReport(ostream & s);
		
		/** generates the code for operators in oplist, and all their subcomponents */
		static void outputVHDLToFile(vector<OperatorPtr> &oplist, ofstream& file);

		
		
	public:
		static vector<OperatorPtr>  globalOpList;  /**< Level-0 operators. Each of these can have sub-operators */
		
	private:
		static vector<OperatorFactoryPtr> sm_factoriesByIndex;
		static map<string,OperatorFactoryPtr> sm_factoriesByName;
	};



	





	
	/** This is the abstract class that each operator factory will inherit. 
			Each OperatorFactory is responsible for the command line and parsing for one Operator sub-class.  */
	class OperatorFactory
	{

	private:
		
		string m_name;
		string m_description;
		vector<string> m_categories;
		vector<string> m_paramNames;
		map<string,string> m_paramType;
		map<string,string> m_paramDoc;
		map<string,string> m_paramDefault; /* If equal to "", this parameter is mandatory (no default)*/
		parser_func_t m_parser;

	public:
		
		/** Implements a no-frills factory
				\param name Name for the operator. The factory will only respond for this name (case sensitive)
				\param description The short documentation
				\param categories A semi-colon seperated list of categories that the operator belongs to. Used to sort the doc.
				\param parameters A semicolon-separated list of parameter description, each being name(type)[=default]:short_description
				\param parser A function that can parse a vector of string arguments into an Operator instance
				\param testParameters Zero or more sets of arguments that the operator can be tested with.
		**/
		OperatorFactory(
						 string name,
						 string description, /**<  for the HTML doc and the detailed help */ 
						 string categories,	/**<  semicolon-seperated list of categories */
						 string parameters, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
						 parser_func_t parser	);
		
		virtual const string &name() const
		{ return m_name; }
		
		virtual const vector<string> &categories() const
		{ return m_categories; }
		

		/** Provide a string with the full documentation. TODO: an HTML version*/
		string getFullDoc();
		
		/*! Consumes zero or more string arguments, and creates an operator
			\param args The offered arguments start at index 0 of the vector, and it is up to the
			factory to check the types and whether there are enough.
			\param consumed On exit, the factory indicates how many of the arguments are used up.
		*/
		virtual OperatorPtr parseArguments(Target* target,  const vector<string> &args	)const
		{
			return m_parser(target, args);
		}


	};
}; // namespace flopoco

#endif
