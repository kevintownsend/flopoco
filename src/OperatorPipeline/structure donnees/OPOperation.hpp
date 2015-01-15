#ifndef OPERATORPIPELINE_OPERATION_HPP_
#define OPERATORPIPELINE_OPERATION_HPP_

#include "OPValue.hpp"
#include "OPExpression.hpp"
#include "OPOperator.hpp"

namespace OperatorPipeline
{

class OPOperation
{
public:
	OPValue* GetValue();

	OPOperation(OPExpression* op_left, OPExpression* op_right, OPOperator* op);
	virtual ~OPOperation();
	
private:
    OPExpression* op_left_;
	OPExpression* op_right_;
	OPOperator* operator_;
};

}

#endif
