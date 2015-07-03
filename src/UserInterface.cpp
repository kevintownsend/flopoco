#include "UserInterface.hpp"
#include "FloPoCo.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>


namespace flopoco
{


	// Allocation of the global objects
	string UserInterface::outputFileName;
	string UserInterface::entityName=""; // used for the -name option
	int UserInterface::verbose;
	string UserInterface::targetFPGA;
	double UserInterface::targetFrequency;
	bool UserInterface::pipeline;
	bool UserInterface::clockEnable;
	bool UserInterface::useHardMult;
	bool UserInterface::plainVHDL;
	bool UserInterface::generateFigures;
	double UserInterface::unusedHardMultThreshold;
	int UserInterface::resourceEstimation;
	bool UserInterface::floorplanning;
	bool UserInterface::reDebug;
	bool UserInterface::flpDebug;

	void UserInterface::parseGenericOptions(vector<string> &args) {
		parseString(args, "name", &entityName, true); // not sticky: will be used, and reset, after the operator parser
		parseString(args, "outputFile", &outputFileName, true); // not sticky: will be used, and reset, after the operator parser
		parsePositiveInt(args, "verbose", &verbose, true); // sticky option
		parseFloat(args, "frequency", &targetFrequency, true); // sticky option
		parseBoolean(args, "useHardMult", &useHardMult, true);
		parseBoolean(args, "plainVHDL", &plainVHDL, true);
		parseBoolean(args, "generateFigures", &generateFigures, true);
		parseBoolean(args, "floorplanning", &floorplanning, true);
		parseBoolean(args, "reDebug", &reDebug, true );
		parseBoolean(args, "pipeline", &pipeline, true );
		//	parseBoolean(args, "", &  );
	}


	
	// Global factory list TODO there should be only one.
	vector<OperatorFactoryPtr> UserInterface::sm_factoriesByIndex;
	map<string,OperatorFactoryPtr> UserInterface::sm_factoriesByName;

	vector<OperatorPtr>  UserInterface::globalOpList;  /**< Level-0 operators. Each of these can have sub-operators */

	
	// This should be obsoleted soon. It is there only because random_main needs it
	void addOperator(OperatorPtr op) {
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
		cerr << "Output file: " << outputFileName <<endl;

	}

	
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

	// TODO make this case-insensitive
	OperatorFactoryPtr UserInterface::getFactoryByName(string operatorName)	{
		return sm_factoriesByName[operatorName];
	}



	void UserInterface::initialize(){
		// Initialize all the command-line options
		verbose=1;
		outputFileName="flopoco.vhdl";
		targetFPGA="Virtex5";
		targetFrequency=400e6;
		useHardMult=true;
	}

		
	void UserInterface::parseAll(int argc, char* argv[]) {
		initialize();

		// manage trivial cases
		if(argc==1) {
			cerr << getFullDoc();
			exit(EXIT_SUCCESS);
		}
		if(argc==2 && string(argv[1])=="BuildHTMLDoc") {
			buildHTMLDoc();
			exit(EXIT_SUCCESS);
		}

		// First convert for convenience the input arg list into
		// 1/ a (possibly empty) vector of global args / initial options,
		// 2/ a vector of operator specification, each being itself a vector of strings 
		vector<string> initialOptions;
		vector<vector<string>> operatorSpecs;


		vector<string> args;
		// convert all the char* to strings
		for (int i=1; i<argc; i++) // start with 1 to skip executable name
			args.push_back(string(argv[i]));

		// Build the global option list
		initialOptions.push_back("$$initialOptions$$");
		while(args.size() > 0 // there remains something to parse
					&& args[0].find("=") !=string::npos) {// and it is an option
			initialOptions.push_back(args[0]);
			args.erase(args.begin());
		}

		// Now there should be at least one operator specification
		while(args.size() > 0) { // This loop is over the Operators that are passed on the command line
			vector<string> opSpec;
			opSpec.push_back(string(args[0]));  // operator Name
			args.erase(args.begin());
			while(args.size() > 0 // there remains something to parse
						&& args[0].find("=") !=string::npos) {// and it is an option
				opSpec.push_back(args[0]);
				args.erase(args.begin());
			}
			operatorSpecs.push_back(opSpec);
		}
	
		cerr<< "Options:";
		for (auto i : initialOptions) cerr << "{"<<i<<"} ";
		cerr << endl;
		cerr<< "Operators:";
		cerr << endl;
		for (auto i : operatorSpecs) {
			for (auto j :i) 
				cerr << "{"<<j<<"} ";
			cerr << endl;
		}
		cerr << endl;

		// Now we have organized our input: do the parsing itself. All the sub-parsers erase the data they consume from the string vectors
		try {
			parseGenericOptions(initialOptions);
			initialOptions.erase(initialOptions.begin());
			if(initialOptions.size()>0){
				ostringstream s;
				cerr<< "Don't know what to do with the following global option(s) :" <<endl ;
				for (auto i : initialOptions)
					s << "  "<<i<<" ";
				s << endl;
				throw s.str();
			}
			
			for (auto opParams: operatorSpecs) {
				cerr<< "parsing operator ";
				for (auto j :opParams) 
					cerr << "{"<<j<<"} ";
				cerr << endl;

				string opName = opParams[0];  // operator Name
				// remove the generic options
				cerr << endl << targetFrequency << endl; 
				parseGenericOptions(opParams);
				cerr << endl << targetFrequency << endl; 

				// build the Target for this operator
				Target* target;
				// make this option case-insensitive, too
				std::transform(targetFPGA.begin(), targetFPGA.end(), targetFPGA.begin(), ::tolower);

					// This could also be a factory but it is less critical
				if(targetFPGA=="virtex4") target=new Virtex4();
				else if (targetFPGA=="virtex5") target=new Virtex5();
				else if (targetFPGA=="virtex6") target=new Virtex6();
				else if (targetFPGA=="spartan3") target=new Spartan3();
				else if (targetFPGA=="stratixii" || targetFPGA=="stratix2") target=new StratixII();
				else if (targetFPGA=="stratixiii" || targetFPGA=="stratix3") target=new StratixIII();
				else if (targetFPGA=="stratixiv" || targetFPGA=="stratix4") target=new StratixIV();
				else if (targetFPGA=="stratixv" || targetFPGA=="stratix5") target=new StratixV();
				else if (targetFPGA=="cycloneii" || targetFPGA=="cyclone2") target=new CycloneII();
				else if (targetFPGA=="cycloneiii" || targetFPGA=="cyclone3") target=new CycloneIII();
				else if (targetFPGA=="cycloneiv" || targetFPGA=="cyclone4") target=new CycloneIV();
				else if (targetFPGA=="cyclonev" || targetFPGA=="cyclone5") target=new CycloneV();
				else {
					throw("ERROR: unknown target: " + targetFPGA);
					}
				target->setPipelined(pipeline);
				target->setFrequency(1e6*targetFrequency);
				target->setUseHardMultipliers(useHardMult);
				target->setPlainVHDL(plainVHDL);
				target->setGenerateFigures(generateFigures);
				// Now build the operator
				OperatorFactoryPtr fp = getFactoryByName(opName);
				OperatorPtr op = fp->parseArguments(target, opParams);
				if(op!=NULL)	{// Some factories don't actually create an operator
					if(entityName!="") {
						op->changeName(entityName);
						entityName="";
					}
					//cerr << "Adding operator" << endl;
					addOperator(op);
				}
			}
		}catch(std::string &s){
			std::cerr<<"Error : "<<s<<"\n";
			//factory->Usage(std::cerr);
			exit(EXIT_FAILURE);
		}catch(std::exception &s){
			std::cerr<<"Exception : "<<s.what()<<"\n";
			//factory->Usage(std::cerr);
			exit(EXIT_FAILURE);	
		}

		// Now output to file
		ofstream file; 
		file.open(outputFileName.c_str(), ios::out);
		outputVHDLToFile(file); 
		file.close();
	}





	// Get the value corresponding to a key, case-insensitive
	string getVal(vector<string>& args, string keyArg){
		// convert to lower case. not efficient to do this each time but hey, this is a user interface.
		std::transform(keyArg.begin(), keyArg.end(), keyArg.begin(), ::tolower);
		vector<string>::iterator i = args.begin();
		// string opName=*i;
		i++;
		while (i != args.end()){
			size_t eqPos = i->find('=');
			if(string::npos==eqPos || 0==eqPos)
				throw ("This doesn't seem to be a key=value pair: " + *i);
			string key= i->substr(0,eqPos);
			// convert to lower case
			std::transform(key.begin(), key.end(), key.begin(), ::tolower);
			if(key==keyArg) {
				cerr <<"  found " << key << endl;
				string val= i->substr(eqPos+1, string::npos);
				cerr <<"  val= " << val << endl;
				// now remove this parameter from the args
				args.erase(i);
				return val;
				cerr <<"  val= " << val << endl;
			}
			i++;
		}
		return ""; // not found
	}
	
	// The following are helper functions to make implementation of factory parsers trivial
	// The code is not efficient but who cares: it is simple to maintain.
	// Beware, args[0] is the operator name, so that we may look up the doc in the factories etc.


	void UserInterface::throwMissingArgError(string opname, string key){
				throw (opname +": argument " + key + " not provided, and there doesn't seem to be a default value."
							 +"\n" +  getFactoryByName(opname) -> getFullDoc());

	}

	void UserInterface::parseString(vector<string> &args, string key, string* variable, bool genericOption){
		string val=getVal(args, key);
		if(val=="") {
			if(genericOption)
				return; // do nothing
			// key not given, use default value
			val = getFactoryByName(args[0])->getDefaultParamVal(key);
			if (val=="")
				throwMissingArgError(args[0], key);
		}
		*variable = val;
	}

	void UserInterface::parseBoolean(vector<string>& args, string key, bool* variable, bool genericOption){
		string val=getVal(args, key);
		if(val=="") {
			if(genericOption)
				return; // do nothing
			// key not given, use default value
			val = getFactoryByName(args[0])->getDefaultParamVal(key);
			if (val=="")
				throwMissingArgError(args[0], key);
		}
		if(val=="1" || val=="yes" || val=="true" || val=="Yes" || val=="True")
			*variable= true;
		else if(val=="0" || val=="no" || val=="false" || val=="No" || val=="False")
			*variable= false;
		else
				throw (args[0] +": expected boolean for argument " + key + ", got " + val);
	}

	void UserInterface::parseFloat(vector<string>& args, string key, double* variable, bool genericOption){
		string val=getVal(args, key);
		if(val=="") {
			if(genericOption)
				return; // do nothing
			// key not given, use default value
			val = getFactoryByName(args[0])->getDefaultParamVal(key);
			if (val=="")
				throwMissingArgError(args[0], key);
		}
		size_t end;
		int intval=stod(val, &end);
		if (val.length() == 0 || val.length() != end)
			throw (args[0] +": expecting a float for parameter " + key + ", got "+val);
		*variable= intval;
	}


	
	void UserInterface::parseInt(vector<string>& args, string key, int* variable, bool genericOption){
		string val=getVal(args, key);
		if(val=="") {
			if(genericOption)
				return; // do nothing
			// key not given, use default value
			val = getFactoryByName(args[0])->getDefaultParamVal(key);
			if (val=="")
				throwMissingArgError(args[0], key);
		}
		size_t end;
		int intval=stoi(val, &end);
		if (val.length() == 0 || val.length() != end)
			throw (args[0] +": expecting an int for parameter " + key + ", got "+val);
		*variable= intval;
	}


	
	void UserInterface::parsePositiveInt(vector<string> &args, string key, int* variable, bool genericOption){
		string val=getVal(args, key);
		if(val=="") {
			if(genericOption) {
				return; // option not found, but it was an option, so do nothing
			}
			else {			// key not given, use default value
				val = getFactoryByName(args[0])->getDefaultParamVal(key);
				if (val=="")
					throwMissingArgError(args[0], key);
 			}
		}
		size_t end;
		
		int intval=stoi(val, &end);
		if (val.length() == 0 || val.length() != end)
			throw (args[0] +": expecting an int for parameter " + key + ", got "+val);
		if(intval>=0)
			*variable = intval;
		else
			throw (args[0] +": expecting strictly positive value for " + key + ", got " + val );
	
	}


	void UserInterface::parseStrictlyPositiveInt(vector<string> &args, string key, int* variable, bool genericOption){
		string val=getVal(args, key);
		if(val=="") {
			if(genericOption)
				return; // do nothing
			// key not given, use default value (except if it is an initial option)
			if(args[0] != "$$initialOptions$$") {
					val = getFactoryByName(args[0])->getDefaultParamVal(key);
					if (val=="")
						throwMissingArgError(args[0], key);
			}
		}
		size_t end;
		int intval=stoi(val, &end);
		if (val.length() == 0 || val.length() != end)
			throw (args[0] +": expecting an int for parameter " + key + ", got "+val);
		if(intval>0)
			*variable = intval;
		else
			throw (args[0] +": expecting strictly positive value for " + key + ", got " + val );
	}


	




	
	
	void UserInterface::add( string name,
													 string description, /**< for the HTML doc and the detailed help */ 
													 string categories,	/**< semicolon-seperated list of categories */
													 string parameterList, /**< semicolon-separated list of parameters, each being name(type)[=default]:short_description  */
													 string extraHTMLDoc, /**< Extra information to go to the HTML doc, for instance links to articles or details on the algorithms */ 
													 parser_func_t parser	 ) {
		OperatorFactoryPtr factory(new OperatorFactory(name, description, categories, parameterList, extraHTMLDoc, parser));
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


	void UserInterface::buildHTMLDoc(){
		ofstream file;
		file.open("doc/web/operators.html", ios::out);
		file << "<!DOCTYPE html>" << endl;
		file << "<html>" << endl;
		file << "<head>" << endl;
		file << "<link rel=\"stylesheet\" href=\"flopoco.css\">" << endl;
		file << "<meta charset=\"utf-8\"> " << endl;
		file << "<title>FloPoCo user manual</title>" << endl;
		file << "</head>" << endl;
		file << "<body>" << endl;

		for(unsigned i = 0; i<getFactoryCount(); i++) {
			OperatorFactoryPtr f =  UserInterface::getFactoryByIndex(i);
			file << f -> getHTMLDoc();
		}
		file << "</body>" << endl;
		file << "</html>" << endl;		
		file.close();
	}





	////////////////// Operator factory /////////////////////////
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


	string OperatorFactory::getHTMLDoc(){
		ostringstream s;
		s << "<dl>"<<endl;
		s << "<dt class=\"operatorname\">" <<  name() << "</dt>"<< endl
			<< "<dd class=\"operatordescription\">"<< m_description << "</dd>" << endl
			<< "<dd><em>Parameters:</em> <dl>" << endl;
		for (unsigned i=0; i<m_paramNames.size(); i++) {
			string pname = m_paramNames[i];
			if("" != m_paramDefault[pname])
				s << "<span class=\"optionalparam\"> " ;
			s << "<dt> <code class=\"parametername\">" << pname << "</code>  (<code class=\"parametertype\">" << m_paramType[pname] << "</code>) " ;
			if("" != m_paramDefault[pname])
				s << "  (optional, default value is " << m_paramDefault[pname] <<")";
			s << "</dt>";
			s << "<dd>" << m_paramDoc[pname] << "</dd>";
			if("" != m_paramDefault[pname])
				s << " </span>";
			s<< endl;			
		}
		s << "</dl></dd>"<<endl;
		if("" != m_extraHTMLDoc)
			s << "<dd>" << m_extraHTMLDoc << "</dd>"<<endl;
		s << "</dl>"<<endl;
		return s.str();
	}

	
	string OperatorFactory::getDefaultParamVal(string& key){
		return  m_paramDefault[key];
	}


		

	OperatorFactory::OperatorFactory(
						 string name,
						 string description, /* for the HTML doc and the detailed help */ 
						 string categories,	/*  semicolon-seperated list of categories */
						 string parameters, /*  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
						 string extraHTMLDoc, /* Extra information to go to the HTML doc, for instance links to articles or details on the algorithms */ 
						 parser_func_t parser  )
		: m_name(name), m_description(description), m_extraHTMLDoc(extraHTMLDoc), m_parser(parser)
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
				//cout << "Found parameter: {" << name<<"}";
				m_paramNames.push_back(name);
				
				int typeEnd = part.find(')', 0);
				string type=part.substr(nameEnd+1, typeEnd-nameEnd-1);
				//cout << " of type  {" << type <<"}";
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
					//cout << " default to {" <<  defaultVal << "}  ";
					m_paramDefault[name]=defaultVal;
					j=defaultValEnd;
				}
				if(part[j]==':') {
					// description
					j++;
					while (part[j]==' ') j++; // remove leading spaces
					string description = part.substr(j, -1);
					m_paramDoc[name] = description;
					//cout << " :  {" << description <<"}" << endl;
				}
				else throw("OperatorFactory: Error parsing parameter description");
			}
			if(end==-1)
				break;
			start=end+1;
			while (parameters[start]==' ') start++;
		}
	}
	


}; // flopoco
