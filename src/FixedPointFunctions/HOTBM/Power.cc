#include "Power.hh"



Power::Power(int d_, Param &p_)
	: d(d_), p(p_), errPow(new PWPolynomial[2])
{
}

Power::~Power()
{
}

PWPolynomial *Power::getErrPow()
{
	return errPow;
}
