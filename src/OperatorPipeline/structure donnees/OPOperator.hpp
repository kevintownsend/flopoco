#ifndef OPERATORPIPELINE_OPERATOR_HPP_
#define OPERATORPIPELINE_OPERATOR_HPP_

#include "OPValue.hpp"
#include "OPExpression.hpp"

namespace OperatorPipeline
{

class OPOperator
{
public:
	OPValue* GetResult(OPExpression op_left, OPExpression op_right)=0;

	OPOperator() = delete;
	virtual ~OPOperator() = delete;
};

}

#endif
