#include "OperatorFactory.hpp"

void addOperator(flopoco::Operator *op);

namespace flopoco
{	
	// Global factory list
	std::vector<OperatorFactoryPtr> OperatorFactoryHolding::sm_factoriesByIndex;
	std::map<std::string,OperatorFactoryPtr> OperatorFactoryHolding::sm_factoriesByName;

void OperatorFactoryHolding::registerFactory(OperatorFactoryPtr factory)
{
	if(sm_factoriesByName.find(factory->name())!=sm_factoriesByName.end())
		throw std::string("OperatorFactory - Factory with name '"+factory->name()+" has already been registered.");

	sm_factoriesByIndex.push_back(factory);
	sm_factoriesByName.insert(std::make_pair(factory->name(), factory));
}

unsigned OperatorFactoryHolding::getFactoryCount()
{
	return sm_factoriesByIndex.size();
}

OperatorFactoryPtr OperatorFactoryHolding::getFactoryByIndex(unsigned i)
{
	return sm_factoriesByIndex.at(i);
}

OperatorFactoryPtr OperatorFactoryHolding::findFactory(std::string operatorName)
{
	return sm_factoriesByName[operatorName];
}




	OperatorFactory::OperatorFactory(
						 std::string name,
						 std::string categories,	/**<  semicolon-seperated list of categories */
						 std::string parameters, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
						 std::string additionalDetails, /**< only for the HTML doc and the detailed help */ 
						 parser_func_t parser
)
		: m_name(name), m_additionalDetails(additionalDetails), m_parser(parser)
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
			if(type=="bool")
				m_paramType[name] = OperatorFactory::Bool;
			else if(type=="int")
				m_paramType[name] = OperatorFactory::Int;
			else if(type=="real")
				m_paramType[name] = OperatorFactory::Real;
			else if(type=="string")
				m_paramType[name] = OperatorFactory::String;
			else {
				ostringstream s;
				s << "OperatorFactory: Type (" << type << ")  not a supported type.";
				throw s.str();
			}
			int j = typeEnd+1;
			if (part[j]=='=') {
				//parsec default value
			}
			else if(part[j]==':') {
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

			
void OperatorFactoryHolding::add(
																 std::string name,
																 std::string categories,	/**<  semicolon-seperated list of categories */
																 std::string parameterList, /**<  semicolon-separated list of parameters, each being name(type)[=default]:short_description  */ 
																 std::string additionalDetails, /**< only for the HTML doc and the detailed help */ 
																 parser_func_t parser	
																 )
{
	OperatorFactoryPtr factory(new OperatorFactory(name, categories, parameterList, additionalDetails, parser));
	OperatorFactoryHolding::registerFactory(factory);
}


}; // flopoco
