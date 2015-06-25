#include "OperatorFactory.hpp"
#include <algorithm>

void addOperator(flopoco::Operator *op);

namespace flopoco
{	

	// Global factory list
	std::vector<OperatorFactoryPtr> UserInterface::sm_factoriesByIndex;
	std::map<std::string,OperatorFactoryPtr> UserInterface::sm_factoriesByName;

	void UserInterface::registerFactory(OperatorFactoryPtr factory)	{
		if(sm_factoriesByName.find(factory->name())!=sm_factoriesByName.end())
			throw std::string("OperatorFactory - Factory with name '"+factory->name()+" has already been registered.");
		
		sm_factoriesByIndex.push_back(factory);
		sm_factoriesByName.insert(std::make_pair(factory->name(), factory));
	}

	unsigned UserInterface::getFactoryCount() {
		return sm_factoriesByIndex.size();
	}

	OperatorFactoryPtr UserInterface::getFactoryByIndex(unsigned i) {
		return sm_factoriesByIndex.at(i);
	}

	OperatorFactoryPtr UserInterface::findFactory(std::string operatorName)	{
		return sm_factoriesByName[operatorName];
	}



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
						 std::string name,
						 std::string description, /**< only for the HTML doc and the detailed help */ 
						 std::string categories,	/**<  semicolon-seperated list of categories */
						 std::string parameters, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
						 parser_func_t parser  )
		: m_name(name), m_description(description), m_parser(parser)
	{
		// Parse the categories
		int start=0;
		while(start<(int)categories.size()){
			int end=categories.find(';', start);
			std::string part;
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
		parameters.erase (std::remove (parameters.begin(), parameters.end(), '\n'), parameters.end());
		start=0;
		while(start<(int)parameters.size()){
			int end=parameters.find(';', start);
			std::string part;
			if(end==-1)
				part=parameters.substr(start, end);
			else
				part=parameters.substr(start, end-start);
			if(part.size()!=0) {
				int nameEnd = part.find('(', 0);
				std::string name=part.substr(0, nameEnd);
				cout << "Found parameter: {" << name<<"}";
				m_paramNames.push_back(name);
				
				int typeEnd = part.find(')', 0);
				std::string type=part.substr(nameEnd+1, typeEnd-nameEnd-1);
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
					std::string defaultVal=part.substr(j, defaultValEnd-j);
					cout << " default to {" <<  defaultVal << "}  ";
					m_paramDefault[name]=defaultVal;
					j=defaultValEnd;
				}
				if(part[j]==':') {
					// description
					j++;
					while (part[j]==' ') j++; // remove leading spaces
					std::string description = part.substr(j, -1);
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
	
	
	UserInterface::param_map_t UserInterface::parseArguments(string opName, const std::vector<std::string> &args, int &consumed){}
	int UserInterface::checkStrictlyPositiveInt(UserInterface::param_map_t paramMap, std::string){
		return 0;
	}
	int UserInterface::checkOptionalInt(UserInterface::param_map_t paramMap, std::string){
		return 0;
	}



	
	void UserInterface::add( std::string name,
													 std::string description, /**< for the HTML doc and the detailed help */ 
													 std::string categories,	/**< semicolon-seperated list of categories */
													 std::string parameterList, /**< semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
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
