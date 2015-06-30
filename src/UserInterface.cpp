#include "UserInterface.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>



namespace flopoco
{

	// This should be obsoleted soon. It is there only because random_main needs it
	void addOperator(OperatorPtr op) {
#if 0 // TODO will be done somewhere else
		if(cl_name!="")	{
			cerr << "Updating entity name to: " << cl_name << endl;
			op->changeName(cl_name);
			cl_name="";
		}
#endif
		UserInterface::globalOpList.push_back(op);
	}



	void UserInterface::addToGlobalOpList(OperatorPtr op) {
		bool alreadyPresent=false;
		// We assume all the operators added to GlobalOpList are unpipelined.
		for (auto i: globalOpList){
			if( op->getName() == i->getName() ) {
					alreadyPresent=true;
					// REPORT(DEBUG,"Operator::addToGlobalOpList(): " << op->getName() <<" already present in globalOpList");
				}
			}
			if(!alreadyPresent)
				globalOpList.push_back(op);
	}
	

	void UserInterface::outputVHDLToFile(ofstream& file){
		outputVHDLToFile(globalOpList, file);
	}

	
	/* The recursive method */
	void UserInterface::outputVHDLToFile(vector<OperatorPtr> &oplist, ofstream& file){
		string srcFileName = "Operator.cpp"; // for REPORT
		for(auto i: oplist) {
			try {
				REPORT(FULL, "---------------OPERATOR: "<<i->getName() <<"-------------");
				REPORT(FULL, "  DECLARE LIST" << printMapContent(i->getDeclareTable()));
				REPORT(FULL, "  USE LIST" << printVectorContent(  (i->getFlopocoVHDLStream())->getUseTable()) );

				// check for subcomponents
				if (! i->getOpList().empty() ){
					//recursively call to print subcomponent
					outputVHDLToFile(i->getOpList(), file);
				}
				i->getFlopocoVHDLStream()->flush();

				/* second parse is only for sequential operators */
				if (i->isSequential()){
					REPORT (FULL, "  2nd PASS");
					i->parse2();
				}
				i->outputVHDL(file);

			} catch (std::string s) {
					cerr << "Exception while generating '" << i->getName() << "': " << s <<endl;
			}
		}
	}


	
	void UserInterface::finalReport(ostream& s){
		s << endl<<"Final report:"<<endl;
		for(auto i: globalOpList) {
			i->outputFinalReport(s, 0);
		}
	}

	
	// Global factory list TODO there should be only one.
	vector<OperatorFactoryPtr> UserInterface::sm_factoriesByIndex;
	map<string,OperatorFactoryPtr> UserInterface::sm_factoriesByName;

	vector<OperatorPtr>  UserInterface::globalOpList;  /**< Level-0 operators. Each of these can have sub-operators */

	
	void UserInterface::registerFactory(OperatorFactoryPtr factory)	{
		if(sm_factoriesByName.find(factory->name())!=sm_factoriesByName.end())
			throw string("OperatorFactory - Factory with name '"+factory->name()+" has already been registered.");
		
		sm_factoriesByIndex.push_back(factory);
		sm_factoriesByName.insert(make_pair(factory->name(), factory));
	}

	unsigned UserInterface::getFactoryCount() {
		return sm_factoriesByIndex.size();
	}

	OperatorFactoryPtr UserInterface::getFactoryByIndex(unsigned i) {
		return sm_factoriesByIndex.at(i);
	}

	OperatorFactoryPtr UserInterface::getFactoryByName(string operatorName)	{
		return sm_factoriesByName[operatorName];
	}


	// Currently very rudimentary
	string OperatorFactory::getFullDoc(){
		ostringstream s;
		s << name() << ": " << m_description << endl << "  Parameters:"<<endl;
		for (unsigned i=0; i<m_paramNames.size(); i++) {
			string pname = m_paramNames[i];
			s << "    " << pname << " (" << m_paramType[pname] << "): " << m_paramDoc[pname] << "  ";
			if("" != m_paramDefault[pname])
				s << "  (optional, default value is " << m_paramDefault[pname] <<")";
			s<< endl;			
		}
		return s.str();
	}

	

	OperatorFactory::OperatorFactory(
						 string name,
						 string description, /**< only for the HTML doc and the detailed help */ 
						 string categories,	/**<  semicolon-seperated list of categories */
						 string parameters, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
						 parser_func_t parser  )
		: m_name(name), m_description(description), m_parser(parser)
	{
		// Parse the categories
		int start=0;
		while(start<(int)categories.size()){
			int end=categories.find(';', start);
			string part;
			if(end==-1)
				part=categories.substr(start, end);
			else
				part=categories.substr(start, end-start);
			if(part.size()!=0)
				m_categories.push_back(part);
			if(end==-1)
				break;
			start=end+1;
		}

		// Parse the parameter description
		// The internet says: this will remove newlines
		parameters.erase (remove (parameters.begin(), parameters.end(), '\n'), parameters.end());
		start=0;
		while(start<(int)parameters.size()){
			int end=parameters.find(';', start);
			string part;
			if(end==-1)
				part=parameters.substr(start, end);
			else
				part=parameters.substr(start, end-start);
			if(part.size()!=0) {
				int nameEnd = part.find('(', 0);
				string name=part.substr(0, nameEnd);
				cout << "Found parameter: {" << name<<"}";
				m_paramNames.push_back(name);
				
				int typeEnd = part.find(')', 0);
				string type=part.substr(nameEnd+1, typeEnd-nameEnd-1);
				cout << " of type  {" << type <<"}";
				if(type=="bool" || type=="int" || type=="real" || type=="string")
					m_paramType[name] = type;
				else {
					ostringstream s;
					s << "OperatorFactory: Type (" << type << ")  not a supported type.";
					throw s.str();
				}
				int j = typeEnd+1;
				m_paramDefault[name]="";
				if (part[j]=='=') {
					//parse default value
					j++;
					int defaultValEnd = part.find(':', 0);
					string defaultVal=part.substr(j, defaultValEnd-j);
					cout << " default to {" <<  defaultVal << "}  ";
					m_paramDefault[name]=defaultVal;
					j=defaultValEnd;
				}
				if(part[j]==':') {
					// description
					j++;
					while (part[j]==' ') j++; // remove leading spaces
					string description = part.substr(j, -1);
					m_paramDoc[name] = description;
					cout << " :  {" << description <<"}" << endl;
				}
				else throw("OperatorFactory: Error parsing parameter description");
			}
			if(end==-1)
				break;
			start=end+1;
			while (parameters[start]==' ') start++;
		}
	}
	

	void UserInterface::parseGlobalOptions(const vector<string> &args) {
		cout << "parsing global options" << endl;
	}


	
	// parseAll does:  while (something left){ 1/consumes global options, 2/ detects an operator name 3/ call the factory argument parser for this name}
	void UserInterface::parseAll(Target* target, int argc, char* argv[]) {
		// First convert the input arg to a vector of strings, for convenience
		vector<string> args;
		for(int i=1;i<argc;i++){ // skip the executable name
			args.push_back(string(argv[i]));
		}
		// Now the parsing itself. All the sub-parsers erase the data they consume from the string vectors
		try {
			cout << "args.size=" << args.size() <<endl;
			while(args.size() > 0) { // This loop is over the Operators that are passed on the command line
				parseGlobalOptions(args); 
				string opName = args[0];  // operator Name
				vector<string> opParams;
				opParams.push_back(args[0]); // place the operator name in position 0
				args.erase(args.begin());
				// First build the tentative param list: it removes complexity from the operator-level parser
				while(args.size()>0                        // there remains something to parse
							&& args[0].find("=") !=string::npos  // and it is a pair key=value
							&& args[0].front() != '-'              // and it is not a global option
						 ) {
					opParams.push_back(args[0]);
					args.erase(args.begin());
				}
				// Now we have consumed the parameters and we are ready to start parsing next global option or operator.

				cout << "Passing the following vector to the factory:" << endl;
				for (unsigned i=0; i< opParams.size(); i++)
					cout << " {" << opParams[i] << "}";
				cout << endl;
				OperatorFactoryPtr fp = getFactoryByName(opName);
				OperatorPtr op = fp->parseArguments(target, opParams);
				if(op!=NULL)	{// Some factories don't actually create an operator
					cout << "Adding operator" << endl;
					addOperator(op);
					cout << "... done" << endl;
				}
			}
		}catch(std::string &s){
			std::cerr<<"Error : "<<s<<"\n";
			//factory->Usage(std::cerr);
			exit(1);
		}catch(std::exception &s){
			std::cerr<<"Exception : "<<s.what()<<"\n";
			//factory->Usage(std::cerr);
			exit(1);	
		}


		
	}



	

	// The following are helper functions to make implementation of factory parsers trivial
	// Beware, args[0] is the operator name, so that we may look up the doc in the factories etc.

	int UserInterface::checkStrictlyPositiveInt(vector<string> args, string keyArg){
		vector<string>::iterator i = args.begin();
		string opName=*i;
		i++;
		while (i != args.end()){
			size_t eqPos = i->find('=');
			if(string::npos==eqPos || 0==eqPos)
				throw ("This doesn't seem to be a key=value pair: " + *i);
			string key= i->substr(0,eqPos);
			if(key==keyArg) {
				string val= i->substr(eqPos+1, string::npos);
				size_t end;
				int intval=stoi(val, &end);
				 if (val.length() == 0 || val.length() != end)
					 throw ("Expecting an int for parameter " + key + ", got "+val);
				if(intval>0)
					return intval;
				else
					throw ("Expecting strictly positive value for " + key + ", got " + val );
			}
			i++;
		} 
		throw ("Key "+keyArg+" not found in arguments of "+opName );
	}

	int UserInterface::checkOptionalInt(vector<string> args, string key){
		return 0;
	}



	
	void UserInterface::add( string name,
													 string description, /**< for the HTML doc and the detailed help */ 
													 string categories,	/**< semicolon-seperated list of categories */
													 string parameterList, /**< semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
													 parser_func_t parser	 ) {
		OperatorFactoryPtr factory(new OperatorFactory(name, description, categories, parameterList, parser));
		UserInterface::registerFactory(factory);
	}



	string UserInterface::getFullDoc(){
		ostringstream s; 
		for(unsigned i = 0; i<getFactoryCount(); i++) {
			OperatorFactoryPtr f =  UserInterface::getFactoryByIndex(i);
			s << f -> getFullDoc();
		}
		return s.str();
	}


}; // flopoco
