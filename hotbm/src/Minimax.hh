#ifndef _MINIMAX_HH_
#define _MINIMAX_HH_

#include "Function.hh"
#include "MPPolynomial.hh"

using namespace std;



class Minimax {
public:
  Minimax(Function &f, double ia, double ib, int d);
  ~Minimax();

  MPPolynomial &getMPP() const;
  void getMPErr(mpfr_t mpErr_) const;

private:
  mpfr_t *buildSystem(Function &f, int d, mpfr_t *mpX);
  mpfr_t *solveSystem(int n, mpfr_t *mpM);
  mpfr_t *solveDicho(Function &f, MPPolynomial &mpP, mpfr_t mpIa, mpfr_t mpIb, int n);
  void findDicho(mpfr_t mpX, Function &f, MPPolynomial &mpP, mpfr_t mpIa, mpfr_t mpIb, int n);

  MPPolynomial *mpP;
  mpfr_t mpErr;
};

#endif // _MINIMAX_HH_
