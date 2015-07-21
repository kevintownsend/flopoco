#include "UserInterface.hpp"
#include "FloPoCo.hpp"
#include "FPDivSqrt/Tools/NbBitsMin.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>

// TODO check the hard mult threshold

namespace flopoco
{

		
	// Colors from	https://github.com/Uduse/Escape-Sequence-Color-Header/blob/master/src/Escape_Sequences_Colors.h
	const char COLOR_NORMAL[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
	const char COLOR_BOLD_BLUE_NORMAL[] = { 0x1b, '[', '1', ';', '3', '4', ';', '4', '9', 'm', 0 };
	const char COLOR_BOLD[] = { 0x1b, '[', '1', 'm', 0 };
	const char COLOR_RED_NORMAL[] = { 0x1b, '[', '3', '1', ';', '4', '9', 'm', 0 };
	const char COLOR_BLUE_NORMAL[] = { 0x1b, '[', '3', '4', ';', '4', '9', 'm', 0 };
	const char COLOR_BOLD_RED_NORMAL[] = { 0x1b, '[', '1', ';', '3', '1', ';', '4', '9', 'm', 0 };
	const char* defaultFPGA="Virtex5";


	// Allocation of the global objects
	string UserInterface::outputFileName;
	string UserInterface::entityName=""; // used for the -name option
	int    UserInterface::verbose;
	string UserInterface::targetFPGA;
	double UserInterface::targetFrequencyMHz;
	bool   UserInterface::pipeline;
	bool   UserInterface::clockEnable;
	bool   UserInterface::useHardMult;
	bool   UserInterface::plainVHDL;
	bool   UserInterface::generateFigures;
	double UserInterface::unusedHardMultThreshold;
<<<<<<< HEAD
	int UserInterface::resourceEstimation;
	bool UserInterface::floorplanning;
	bool UserInterface::reDebug;
	bool UserInterface::flpDebug;

=======
	int    UserInterface::resourceEstimation;
	bool   UserInterface::floorplanning;
	bool   UserInterface::reDebug;
	bool   UserInterface::flpDebug;


		
		
	void UserInterface::main(int argc, char* argv[]) {
		try {
			sollya_lib_init();
			initialize();
			buildAll(argc, argv);
			outputVHDL();
			finalReport(cerr); 
			sollya_lib_close();
		}
		catch (string e) {
			cerr << endl << e;
		}
	}
	





	void UserInterface::parseGenericOptions(vector<string> &args) {
		parseString(args, "name", &entityName, true); // not sticky: will be used, and reset, after the operator parser
		parseString(args, "outputFile", &outputFileName, true); // not sticky: will be used, and reset, after the operator parser
		parseString(args, "target", &targetFPGA, true); // not sticky: will be used, and reset, after the operator parser
		parsePositiveInt(args, "verbose", &verbose, true); // sticky option
		parseFloat(args, "frequency", &targetFrequencyMHz, true); // sticky option
		parseFloat(args, "hardMultThreshold", &unusedHardMultThreshold, true); // sticky option
		parseBoolean(args, "useHardMult", &useHardMult, true);
		parseBoolean(args, "plainVHDL", &plainVHDL, true);
		parseBoolean(args, "generateFigures", &generateFigures, true);
		parseBoolean(args, "floorplanning", &floorplanning, true);
		parseBoolean(args, "reDebug", &reDebug, true );
		parseBoolean(args, "pipeline", &pipeline, true );
		//	parseBoolean(args, "", &  );
	}


	
>>>>>>> origin/newCLI
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
		
		// Messages for testbenches. Only works if you have only one TestBench
		Operator* op = globalOpList.back();
		if(op->getSrcFileName() == "TestBench"){
			cerr << "To run the simulation using ModelSim, type the following in 'vsim -c':" <<endl;
			cerr << tab << "vdel -all -lib work" <<endl;
			cerr << tab << "vlib work" <<endl;
			cerr << tab << "vcom " << outputFileName <<endl;
			cerr << tab << "vsim " << op->getName() <<endl;
			cerr << tab << "add wave -r *" <<endl;
			cerr << tab << "run " << ((TestBench*)op)->getSimulationTime() << "ns" << endl;
			cerr << "To run the simulation using gHDL, type the following in a shell prompt:" <<endl;
			string simlibs;
#if 0
			if(op->getStdLibType()==0 || op->getStdLibType()==-1)
				simlibs="--ieee=synopsys ";
			if(op->getStdLibType()==1)
				simlibs="--ieee=standard ";
#else
				simlibs="--ieee=standard --ieee=synopsys ";
#endif
			cerr <<  "ghdl -a " << simlibs << "-fexplicit "<< outputFileName <<endl;
			cerr <<  "ghdl -e " << simlibs << "-fexplicit " << op->getName() <<endl;
			cerr <<  "ghdl -r " << simlibs << op->getName() << " --vcd=" << op->getName() << ".vcd --stop-time=" << ((TestBench*)op)->getSimulationTime() << "ns" <<endl;
			cerr <<  "gtkwave " << op->getName() << ".vcd" << endl;
		}
		
	}


	void UserInterface::registerFactory(OperatorFactoryPtr factory)	{
		if(sm_factoriesByName.find(factory->name())!=sm_factoriesByName.end())
			throw string("OperatorFactory - Factory with name '"+factory->name()+" has already been registered.");
<<<<<<< HEAD

		sm_factoriesByIndex.push_back(factory);
=======
		
>>>>>>> origin/newCLI
		sm_factoriesByName.insert(make_pair(factory->name(), factory));
		sm_factoriesByIndex.push_back(factory);
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

	string categoryString(UserInterface::DocumentationCategory c){
		switch(c) {
		case UserInterface::ShiftersLZOCs:
			return "Shifters, Leading Zero Counters, etc";
		case UserInterface::BasicInteger:
			return "Basic Integer operators (pipelined)";
		case UserInterface::BasicFixPoint:
			return "Basic Fixed-point Operators";
		case UserInterface::BasicFloatingPoint:
			return "Basic Floating-point Operators";
		case UserInterface::CompositeFloatingPoint:
			return "Composite Floating-point Operators";
		case UserInterface::ElementaryFunctions:
			return "Elementary Functions in Fixed- or Floating-Point";
		case UserInterface::FunctionApproximation:
			return "Arbitrary Function Approximators";
		case UserInterface::ComplexFixPoint:
			return "Complex Fixed-Point Arithmetic Operators";
		case UserInterface::ComplexFloatingPoint:
			return "Complex Floating-Point Arithmetic Operators";
		case UserInterface::LNS:
			return "Logarithm Number System Operators";
		case UserInterface::Conversions:
			return "Conversions Between Various Number Formats";
		case UserInterface::TestBenches:
			return "Test Benches";
		case UserInterface::Miscellanous:
		return "Miscellanous";
		default: return"";
		}
	}
	
	void UserInterface::initialize(){
		// Initialize all the command-line options
		verbose=1;
		outputFileName="flopoco.vhdl";
		targetFPGA=defaultFPGA;
		targetFrequencyMHz=400;
		useHardMult=true;
		unusedHardMultThreshold=0.7;
	}


<<<<<<< HEAD
	void UserInterface::parseAll(int argc, char* argv[]) {
		initialize();
		// First convert the input arg to a vector of strings, for convenience
=======

	void UserInterface::buildAll(int argc, char* argv[]) {

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


>>>>>>> origin/newCLI
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
<<<<<<< HEAD
		// Now the parsing itself. All the sub-parsers erase the data they consume from the string vectors
		try {
			if(args.size()==0) {
				cerr << getFullDoc();
				exit(EXIT_SUCCESS);
			}
			if(args.size()==1 && args[0]=="BuildHTMLDoc") {
				buildHTMLDoc();
				exit(EXIT_SUCCESS);
			}
			if(args.size()==3 && args[0]=="NbBitsMin") {
				computeNbBit(atoi(args[1].c_str()), atoi(args[2].c_str()));
				exit(EXIT_SUCCESS);
			}

			//cout << "args.size=" << args.size() <<endl;
			while(args.size() > 0) { // This loop is over the Operators that are passed on the command line
				parseGenericOptions(args);
				string opName = args[0];  // operator Name
				vector<string> opParams;
				opParams.push_back(args[0]); // place the operator name in position 0
=======

		// Now there should be at least one operator specification
		while(args.size() > 0) { // This loop is over the Operators that are passed on the command line
			vector<string> opSpec;
			opSpec.push_back(string(args[0]));  // operator Name
			args.erase(args.begin());
			while(args.size() > 0 // there remains something to parse
						&& args[0].find("=") !=string::npos) {// and it is an option
				opSpec.push_back(args[0]);
>>>>>>> origin/newCLI
				args.erase(args.begin());
			}
			operatorSpecs.push_back(opSpec);
		}
	

		// Now we have organized our input: do the parsing itself. All the sub-parsers erase the data they consume from the string vectors
		try {
			parseGenericOptions(initialOptions);
			initialOptions.erase(initialOptions.begin());
			if(initialOptions.size()>0){
				ostringstream s;
				s << "Don't know what to do with the following global option(s) :" <<endl ;
				for (auto i : initialOptions)
					s << "  "<<i<<" ";
				s << endl;
				throw s.str();
			}
			
			for (auto opParams: operatorSpecs) {

				string opName = opParams[0];  // operator Name
				// remove the generic options
				parseGenericOptions(opParams);

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
				target->setFrequency(1e6*targetFrequencyMHz);
				target->setUseHardMultipliers(useHardMult);
				target->setPlainVHDL(plainVHDL);
				target->setGenerateFigures(generateFigures);
				// Now build the operator
				OperatorFactoryPtr fp = getFactoryByName(opName);
				if (fp==NULL){
					throw( "Can't find the operator factory for " + opName) ;
				}
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
	}


<<<<<<< HEAD
		// Now output to file
		ofstream file;
=======
	void UserInterface::outputVHDL() {
		ofstream file; 
>>>>>>> origin/newCLI
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
				//cerr <<"  found " << key << endl;
				string val= i->substr(eqPos+1, string::npos);
				//cerr <<"  val= " << val << endl;
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

<<<<<<< HEAD

	bool UserInterface::checkBoolean(vector<string> args, string key){
=======
	}

	void UserInterface::parseString(vector<string> &args, string key, string* variable, bool genericOption){
>>>>>>> origin/newCLI
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
		double dval=stod(val, &end);
		if (val.length() == 0 || val.length() != end)
			throw (args[0] +": expecting a float for parameter " + key + ", got "+val);
		*variable= dval;
	}


<<<<<<< HEAD

	int UserInterface::checkInt(vector<string> args, string key){
=======
	
	void UserInterface::parseInt(vector<string>& args, string key, int* variable, bool genericOption){
>>>>>>> origin/newCLI
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





<<<<<<< HEAD
	void UserInterface::parseGenericOptions(vector<string> &args) {
		cout << "parsing generic options" << endl;
		entityName=getVal(args, "name"); // will be used, and reset, after the operator parser
	}

=======
>>>>>>> origin/newCLI




	void UserInterface::add( string name,
<<<<<<< HEAD
													 string description, /**< for the HTML doc and the detailed help */
													 string categories,	/**< semicolon-seperated list of categories */
=======
													 string description, /**< for the HTML doc and the detailed help */ 
													 DocumentationCategory category,
													 string seeAlso,
>>>>>>> origin/newCLI
													 string parameterList, /**< semicolon-separated list of parameters, each being name(type)[=default]:short_description  */
													 string extraHTMLDoc, /**< Extra information to go to the HTML doc, for instance links to articles or details on the algorithms */
													 parser_func_t parser	 ) {
		OperatorFactoryPtr factory(new OperatorFactory(name, description, category, seeAlso, parameterList, extraHTMLDoc, parser));
		UserInterface::registerFactory(factory);
	}


#if 0
	const int outputToHTML=1;
	const int outputToConsole=2;
	
	string colorParameter(string s, int techno, bool optional) {
		string o
		if (techno==outputToHTML)
			o = "<code class=\"parametername\">" + s + "</code>";
		else 	if (techno==outputToConsole)
			o = (optional?COLOR_BOLD_RED_NORMAL:COLOR_BOLD) + s + COLOR_NORMAL;
		return o;
	}
#endif

	
	string UserInterface::getFullDoc(){
		ostringstream s;
<<<<<<< HEAD
=======
		s << "Usage: " << COLOR_BOLD << "flopoco  [options]  OperatorName parameters  [OperatorName parameters]..." << COLOR_NORMAL << endl;
		s << "  Both options and parameters are lists of " << COLOR_BOLD << "name=value" << COLOR_NORMAL << " pairs (with case-insensitive name)" << endl;
		s << COLOR_BLUE_NORMAL<< "Example: " << COLOR_NORMAL << "flopoco  frequency=300 target=Virtex5   FPExp  wE=8 wF=23 name=SinglePrecisionFPExp" << endl;
		s << "Generic options include:" << endl;
		s << "  " << COLOR_BOLD << "name" << COLOR_NORMAL << "=<string>:        override the the default entity name "<<endl;
		s << "  " << COLOR_BOLD << "outputFile" << COLOR_NORMAL << "=<string>:  override the the default output file name " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL <<endl;
		s << "  " << COLOR_BOLD << "pipeline" << COLOR_NORMAL << "=<0|1>:       pipelined operator, or not " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL << endl;
		s << "  " << COLOR_BOLD << "target" << COLOR_NORMAL << "=<string>:      target FPGA (default " << defaultFPGA << ") " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL<<endl;
		s << "     Supported targets: Stratix2...5, Virtex2...6, Cyclone2...5,Spartan3"<<endl;
		s << "  " << COLOR_BOLD << "frequency" << COLOR_NORMAL << "=<float>:    target frequency in MHz (default 400) " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL<<endl;
		s << "  " << COLOR_BOLD << "plainVHDL" << COLOR_NORMAL << "=<0|1>:      use plain VHDL (default), or not " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL << endl;
		s << "  " << COLOR_BOLD << "hardMultThreshold" << COLOR_NORMAL << "=<float>: unused hard mult threshold (O..1, default 0.7) " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL<<endl;
		s << "  " << COLOR_BOLD << "generateFigures" << COLOR_NORMAL << "=<0|1>:generate SVG graphics (default off) " << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL << endl;
		s << "  " << COLOR_BOLD << "verbose" << COLOR_NORMAL << "=<int>:        verbosity level (0-4, default=1)" << COLOR_RED_NORMAL << "(sticky option)" << COLOR_NORMAL<<endl;
		s << "Sticky options apply to the rest of the command line, unless changed again" <<endl;
>>>>>>> origin/newCLI
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
		s <<COLOR_BOLD_BLUE_NORMAL << name() << COLOR_NORMAL <<": " << m_description << endl;
		for (unsigned i=0; i<m_paramNames.size(); i++) {
			string pname = m_paramNames[i];
			s << "  " << ("" != m_paramDefault[pname]?COLOR_BOLD_RED_NORMAL:COLOR_BOLD) << pname <<COLOR_NORMAL<< " (" << m_paramType[pname] << "): " << m_paramDoc[pname] << "  ";
			if("" != m_paramDefault[pname])
<<<<<<< HEAD
				s << "  (optional, default value is " << m_paramDefault[pname] <<")";
			s<< endl;
=======
				s << COLOR_RED_NORMAL << "  (optional, default value is " << m_paramDefault[pname] <<")"<< COLOR_NORMAL;
			s<< endl;			
>>>>>>> origin/newCLI
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
<<<<<<< HEAD
						 string description, /* for the HTML doc and the detailed help */
						 string categories,	/*  semicolon-seperated list of categories */
						 string parameters, /*  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */
						 string extraHTMLDoc, /* Extra information to go to the HTML doc, for instance links to articles or details on the algorithms */
=======
						 string description, /* for the HTML doc and the detailed help */ 
						 UserInterface::DocumentationCategory category,
						 string seeAlso,
						 string parameters, /*  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
						 string extraHTMLDoc, /* Extra information to go to the HTML doc, for instance links to articles or details on the algorithms */ 
>>>>>>> origin/newCLI
						 parser_func_t parser  )
		: m_name(name), m_description(description), m_category(category), m_seeAlso(seeAlso), m_extraHTMLDoc(extraHTMLDoc), m_parser(parser)
	{
		int start;
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
					while (part[j]==' ' || part[j]=='\t') j++; // remove leading spaces and tabs
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
