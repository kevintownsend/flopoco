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

			
void OperatorFactoryHolding::add(
	std::string name,			
	std::string categories,	// semi-colon seperated list of categories
	usage_func_t usage,
	parser_func_t parser
){
	OperatorFactoryPtr factory(new OperatorFactory(name, categories, usage, parser));
	OperatorFactoryHolding::registerFactory(factory);
}


}; // flopoco
