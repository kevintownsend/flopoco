#ifndef OPERATORPIPELINE_EXPRESSION_HPP_
#define OPERATORPIPELINE_EXPRESSION_HPP_

#include "OPValue.hpp"

namespace OperatorPipeline
{

class OPExpression
{
public:
	virtual OPValue* GetValue()=0;

	virtual ~OPExpression();
};

}

#endif
