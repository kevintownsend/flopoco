#include "OperatorFactory.hpp"

void addOperator(flopoco::Operator *op);

namespace flopoco
{	

std::vector<OperatorFactoryPtr> OperatorFactory::sm_factoriesByIndex;
std::map<std::string,OperatorFactoryPtr> OperatorFactory:: sm_factoriesByName;

void OperatorFactory::registerFactory(OperatorFactoryPtr factory)
{
	if(sm_factoriesByName.find(factory->name())!=sm_factoriesByName.end())
		throw std::string("OperatorFactory - Factory with name '"+factory->name()+" has already been registered.");

	sm_factoriesByIndex.push_back(factory);
	sm_factoriesByName.insert(std::make_pair(factory->name(), factory));
}

unsigned OperatorFactory::getFactoryCount()
{
	return sm_factoriesByIndex.size();
}

OperatorFactoryPtr OperatorFactory::getFactoryByIndex(unsigned i)
{
	return sm_factoriesByIndex.at(i);
}

OperatorFactoryPtr OperatorFactory::findFactory(std::string operatorName)
{
	return sm_factoriesByName[operatorName];
}




	DefaultOperatorFactory::SimpleOperatorFactory::SimpleOperatorFactory(
		std::string name,			
		std::string categories,
		usage_func_t usage,
		parser_func_t parser
)
	: m_name(name)
	, m_usage(usage)
	, m_parser(parser)
{
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
}

			
void DefaultOperatorFactory::add(
	std::string name,			
	std::string categories,	// semi-colon seperated list of categories
	usage_func_t usage,
	parser_func_t parser
){
	OperatorFactoryPtr factory(new SimpleOperatorFactory(name, categories, usage, parser));
	OperatorFactory::registerFactory(factory);
}


}; // flopoco
