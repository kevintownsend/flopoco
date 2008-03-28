#include "Term.hh"



Term::Term(int d_, double *k_, Param &p_)
  : d(d_), k(k_), p(p_), errMethod(new PWPolynomial[2]), refCount(0)
{
}

Term::~Term()
{
  delete[] errMethod;
}

PWPolynomial *Term::getErrMethod()
{
  return errMethod;
}
