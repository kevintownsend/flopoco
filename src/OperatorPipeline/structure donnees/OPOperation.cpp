#include "OPOperation.hpp"

namespace OperatorPipeline
{

OPOperation::GetValue()
{
	//\todo
}

OPOperation::OPOperation(OPExpression* op_left, OPExpression* op_right, OPOperator* op)
{
    op_left_ = op_left;
	op_right_ = op_right;
	operator_ = op;
}

OPOperation::~OPOperation()
{
	delete operator_;
	delete op_right_;
	delete op_left_;
}

}


