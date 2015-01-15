#ifndef OPERATORPIPELINE_VARIABLE_HPP_
#define OPERATORPIPELINE_VARIABLE_HPP_

#include "OPValue.hpp"
#include "OPExpression.hpp"

namespace OperatorPipeline
{

typedef enum type_enum
{
	UINT,
	SINT
} type_enum;

class OPVariable
{
public:
	virtual OPValue* GetValue()=0;

	OPVariable(OPValue* value, type_enum type, int* type_param, OPExpression* expresion);
	virtual ~OPVariable();
	
private:
	OPValue *value_;
	type_enum type_;
	int* type_param_;
	OPExpression* expression;
};

}

#endif
