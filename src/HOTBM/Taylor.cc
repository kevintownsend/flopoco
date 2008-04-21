#include <mpfr.h>
#include "Taylor.hh"
#include "sollya.h"

/*
 * Uses libsollya to compute the taylor polynomial
 */
Taylor::Taylor(Function &f_, int d, double x0)
  : poly(NULL)
{
	/* Convert parameters to their required type */
	mpfr_t mpx0;
	mpfr_init(mpx0);
	mpfr_set_d(mpx0, x0, GMP_RNDN);
	sollya_node_t nx0 = makeConstant(mpx0);

	/* Create input functions */
	sollya_node_t f = parseString(f_.getName().c_str());

	/* Call taylor */
	sollya_node_t nTaylor = taylor(f, d, nx0, getToolPrecision());

	/* Extract coefficients */
	int degree;
	sollya_node_t *nCoef;
	double *coef;
	mpfr_t temp;

	getCoefficients(&degree, &nCoef, nTaylor);
	coef = new double[degree+1];
	mpfr_init(temp);
	int i;
	for (i = 0; i <= degree; i++)
	{
		evaluateConstantExpression(temp, nCoef[i], getToolPrecision());
		coef[i] = mpfr_get_d(temp, GMP_RNDN);
	}

	/* Create the long-awaited polynomial */
	poly = new Polynomial(degree, coef);

	/* Cleanup */
	mpfr_clear(temp);
	delete coef;
	free_memory(f);
	free_memory(nx0);
}


Taylor::~Taylor()
{
  if (poly)
    delete poly;
}

Polynomial &Taylor::getPoly() const
{
  return *poly;
}

