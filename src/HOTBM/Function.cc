#include "Function.hh"

Function::Function(string name_)
  : name(name_)
{
	/* Convert the input string into a sollya evaluation tree */
	node = parseString(name_.c_str());

	/* If conversion did not succeed (i.e. parse error)
	 * throw an exception */
	if (node == 0)
		throw "Unable to parse input function.";
}

Function::~Function()
{
	free_memory(node);
}

string Function::getName() const
{
  return name;
}

void Function::eval(mpfr_t r, mpfr_t x) const
{
	evaluateFaithful(r, node, x, getToolPrecision());
}

double Function::eval(double x) const
{
	mpfr_t mpX, mpR;
	double r;

	mpfr_inits(mpX, mpR, NULL);
	mpfr_set_d(mpX, x, GMP_RNDN);
	evaluateFaithful(mpR, node, mpX, getToolPrecision());
	r = mpfr_get_d(mpR, GMP_RNDN);
	mpfr_clears(mpX, mpR, NULL);

	return r;
}

sollya_node_t Function::getSollyaNode() const
{
	return node;
}

