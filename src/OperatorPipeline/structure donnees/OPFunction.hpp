#ifndef OPERATORPIPELINE_FUNCTION_HPP_
#define OPERATORPIPELINE_FUNCTION_HPP_

#include "OPValue.hpp"
#include "OPExpression.hpp"

namespace OperatorPipeline
{

typedef enum function_enum
{
	LOG,
	EXP
} function_enum;

class OPFunction
{
public:
	virtual OPValue* GetValue()=0;

	OPFunction(function_enum function, OPExpression *operande);
	~OPFunction();
	
private:
	OPExpression *operande_;
	function_enum function_;
};

}

#endif
