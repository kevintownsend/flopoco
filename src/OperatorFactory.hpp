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
#ifndef flopoco_operator_factory_hpp
#define flopoco_operator_factory_hpp

#include "Operator.hpp"
#include <memory>

// Operator Factory, based on the one by David Thomas, with a bit of clean up.
// For typical use, see src/ShiftersEtc/Shifter.*

namespace flopoco
{
	
		// Note: not using boost::function here, as it's likely to scare people, and also drags in quite a few header dependencies
	typedef OperatorPtr (*parser_func_t)(Target *,const vector<string> &,int &);	
	class OperatorFactory;
	typedef shared_ptr<OperatorFactory> OperatorFactoryPtr;



	
	/** This is the class that manages a list of OperatorFactories, and the overall command line and documentation.
			Each OperatorFactory is responsible for the command line and parsing for one Operator sub-class. */
	class UserInterface
	{
	public:
		typedef pair<string, map<string, string>> param_map_t;

		static void registerFactory(OperatorFactoryPtr factory);
		static void add(
										string name,
										string description, /**< for the HTML doc and the detailed help */ 
										string categories,	/**<  semicolon-seperated list of categories */
										string parameterList, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
										parser_func_t parser	);
		

		static param_map_t  parseArguments(string opName, const vector<string> &args, int &consumed);
		static int checkStrictlyPositiveInt(UserInterface::param_map_t, string);
		static int checkOptionalInt(UserInterface::param_map_t, string);

		/** Provide a string with the full documentation. TODO: an HTML version*/
		static string getFullDoc();
		
		static unsigned getFactoryCount();
		static OperatorFactoryPtr getFactoryByIndex(unsigned i);
		static OperatorFactoryPtr findFactory(string operatorName);
		

		
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
		virtual OperatorPtr parseCommandLine(
																			 Target *target,
																			 const vector<string> &args,
																			 int &consumed
																			 )const
		{
			return m_parser(target, args, consumed);
		}



	};
}; // namespace flopoco

#endif
