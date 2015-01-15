#include "OPVariable.hpp"

namespace OperatorPipeline
{

OPVariable::OPVariable(OPValue* value, type_enum type, int* type_param, OPExpression* expresion)
{
	value_ = value;
	type_ = type;
	type_param_ = type_param;
	expression_ = expression;
}

OPVariable::~OPVariable()
{
	delete value_;
	delete expression_;
}

}
