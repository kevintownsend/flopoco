#include "Minimax.hh"

#define EPSILON 1.0e-12



Minimax::Minimax(Function &f, double ia, double ib, int d)
  : mpP(NULL)
{
	{
		fprintf(stderr, "\n>>> Call to sollya: %s",
			f.getName().c_str());
		mpfr_t x;
		mpfr_init(x);
		mpfr_set_d(x, ia, GMP_RNDN);
		mpfr_out_str(stderr, 10, 0, x, GMP_RNDN);
		printf(" ");
		mpfr_set_d(x, ib, GMP_RNDN);
		mpfr_out_str(stderr, 10, 0, x, GMP_RNDN);
		mpfr_clear(x);
		fprintf(stderr, "\n");
	}


  mp_prec_t oldPrec = mpfr_get_default_prec();
  mpfr_set_default_prec(80);

  cerr << "*";

  bool done = false, loop = false;
  int nTry = 1;
  while (!done) {
    mpfr_t mpIa, mpDiff, *mpX;

    mpfr_inits(mpErr, mpIa, mpDiff, NULL);
    mpfr_set_d(mpIa, ia, GMP_RNDN);

    mpX = new mpfr_t[d+2];
    for (int i = 0; i < d+2; i++) {
      mpfr_init(mpX[i]);
      mpfr_set_d(mpX[i], ib-ia, GMP_RNDN);
      mpfr_mul_ui(mpX[i], mpX[i], i, GMP_RNDN);
      mpfr_div_ui(mpX[i], mpX[i], d+1, GMP_RNDN);
      mpfr_add(mpX[i], mpX[i], mpIa, GMP_RNDN);
    }

    try {
      while (!done) {
	mpfr_t *mpM = buildSystem(f, d, mpX);
	mpfr_t *mpK = solveSystem(d+2, mpM);

	for (int i = 0; i < d+2; i++)
	  for (int j = 0; j < d+3; j++)
	    mpfr_clear(mpM[i*(d+3)+j]);
	delete[] mpM;

	mpP = new MPPolynomial(d, mpK);
	mpfr_abs(mpErr, mpK[d+1], GMP_RNDN);

	for (int i = 0; i < d+2; i++)
	  mpfr_clear(mpK[i]);
	delete[] mpK;

	mpfr_t *mpNX;
	try {
	  mpNX = solveDicho(f, *mpP, mpX[0], mpX[d+1], 1);
	}
	catch (const char *s) {
	  delete mpP;
	  mpP = NULL;
	  throw;
	}

	mpfr_set_ui(mpX[0], 0, GMP_RNDN);
	for (int i = 1; i < d+1; i++) {
	  mpfr_sub(mpX[i], mpX[i], mpNX[i], GMP_RNDN);
	  mpfr_div(mpX[i], mpX[i], mpNX[i], GMP_RNDN);
	  mpfr_abs(mpX[i], mpX[i], GMP_RNDN);
	  if (mpfr_cmp(mpX[i], mpX[0]) > 0)
	    mpfr_set(mpX[0], mpX[i], GMP_RNDN);
	}
	
	loop = mpfr_cmp(mpX[0], mpDiff) > 0;
	done = mpfr_cmp_d(mpX[0], EPSILON) < 0;
	mpfr_set(mpDiff, mpX[0], GMP_RNDN);

	for (int i = 0; i < d+2; i++)
	  mpfr_clear(mpX[i]);
	delete[] mpX;
	mpX = mpNX;

	if (!done) {
	  delete mpP;
	  mpP = NULL;
	  if (loop)
	    throw "Minimax::Minimax: Cannot stabilize to required precision.";
	}
      }
    }
    catch (const char *s) {
      if (++nTry > 3) {
	cerr << endl;
	throw "Minimax::Minimax: Unable to find minimax approximation with 320-bit precision.";
      }
      cerr << "!";
      mpfr_set_default_prec(2 * mpfr_get_default_prec());
    }

    for (int i = 0; i < d+2; i++)
      mpfr_clear(mpX[i]);
    delete[] mpX;

    if (!done)
      mpfr_clear(mpErr);
    mpfr_clears(mpIa, mpDiff, NULL);
  }

  for (; nTry < 3; nTry++)
    cerr << " ";

  mpfr_set_default_prec(oldPrec);

	mpfr_t x;
	mpfr_init(x);
	fprintf(stderr, "Coef:");
	int i;
	for (i = 0; i <= d; i++)
	{
		mpP->getMPK(x, i);
		fprintf(stderr," ");
		mpfr_out_str(stderr, 10, 0, x, GMP_RNDZ);
	}
	fprintf(stderr, "\nErr: ");

	mpfr_out_str(stderr, 10, 0, mpErr, GMP_RNDZ);
	fprintf(stderr, "\n");

}

Minimax::~Minimax()
{
  if (mpP)
    delete mpP;
  mpfr_clear(mpErr);
  mpfr_free_cache();
}

MPPolynomial &Minimax::getMPP() const
{
  return *mpP;
}

void Minimax::getMPErr(mpfr_t mpErr_) const
{
  mpfr_set(mpErr_, mpErr, GMP_RNDN);
}

mpfr_t *Minimax::buildSystem(Function &f, int d, mpfr_t *mpX)
{
  mpfr_t *mpM = new mpfr_t[(d+2)*(d+3)];

  for (int i = 0; i < d+2; i++) {
    for (int j = 0; j < d+3; j++)
      mpfr_init(mpM[i*(d+3)+j]);
    mpfr_set_ui(mpM[i*(d+3)], 1, GMP_RNDN);
    for (int j = 1; j <= d; j++)
      mpfr_mul(mpM[i*(d+3)+j], mpM[i*(d+3)+j-1], mpX[i], GMP_RNDN);
    mpfr_set_si(mpM[i*(d+3)+d+1], i%2 ? 1 : -1, GMP_RNDN);
    f.mpEval(mpM[i*(d+3)+d+2], mpX[i]);
  }

  return mpM;
}

mpfr_t *Minimax::solveSystem(int n, mpfr_t *mpM)
{
  mpfr_t *mpX = new mpfr_t[n];
  mpfr_t mpTmp;

  for (int i = 0; i < n; i++)
    mpfr_init(mpX[i]);
  mpfr_init(mpTmp);

  if (n == 1)
    mpfr_div(mpX[0], mpM[1], mpM[0], GMP_RNDN);
  else {
    mpfr_t *mpSM = new mpfr_t[(n-1)*n];

    for (int i = 0; i < n-1; i++) {
      mpfr_div(mpTmp, mpM[(i+1)*(n+1)], mpM[0], GMP_RNDN);
      for (int j = 0; j < n; j++) {
	mpfr_init(mpSM[i*n+j]);
	mpfr_mul(mpSM[i*n+j], mpTmp, mpM[j+1], GMP_RNDN);
	mpfr_sub(mpSM[i*n+j], mpM[(i+1)*(n+1)+j+1], mpSM[i*n+j], GMP_RNDN);
      }
    }

    mpfr_t *mpSX = solveSystem(n-1, mpSM);

    mpfr_set(mpX[0], mpM[n], GMP_RNDN);
    for (int i = 0; i < n-1; i++) {
      mpfr_set(mpX[i+1], mpSX[i], GMP_RNDN);
      mpfr_mul(mpTmp, mpM[i+1], mpX[i+1], GMP_RNDN);
      mpfr_sub(mpX[0], mpX[0], mpTmp, GMP_RNDN);
    }
    mpfr_div(mpX[0], mpX[0], mpM[0], GMP_RNDN);

    for (int i = 0; i < n-1; i++) {
      for (int j = 0; j < n; j++)
	mpfr_clear(mpSM[i*n+j]);
      mpfr_clear(mpSX[i]);
    }
    delete[] mpSM;
    delete[] mpSX;
  }

  mpfr_clear(mpTmp);

  return mpX;
}

mpfr_t *Minimax::solveDicho(Function &f, MPPolynomial &mpP, mpfr_t mpIa, mpfr_t mpIb, int n)
{
  int d = mpP.getD();
  mpfr_t *mpX = new mpfr_t[d-n+3];

  for (int i = 0; i < d-n+3; i++)
    mpfr_init(mpX[i]);
  mpfr_set(mpX[0], mpIa, GMP_RNDN);
  mpfr_set(mpX[d-n+2], mpIb, GMP_RNDN);

  if (n <= d) {
    try {
      mpfr_t *mpSX = solveDicho(f, mpP, mpIa, mpIb, n+1);
      try {
	for (int i = 0; i < d-n+1; i++)
	  findDicho(mpX[i+1], f, mpP, mpSX[i], mpSX[i+1], n);
      }
      catch (const char *) {
	for (int i = 0; i < d-n+2; i++)
	  mpfr_clear(mpSX[i]);
	delete[] mpSX;
	throw;
      }

      for (int i = 0; i < d-n+2; i++)
	mpfr_clear(mpSX[i]);
      delete[] mpSX;
    }
    catch (const char *) {
      for (int i = 0; i < d-n+3; i++)
	mpfr_clear(mpX[i]);
      delete[] mpX;
      throw;
    }
  }

  return mpX;
}

void Minimax::findDicho(mpfr_t mpX, Function &f, MPPolynomial &mpP, mpfr_t mpIa, mpfr_t mpIb, int n)
{
  mpfr_t mpIm;

  mpfr_init(mpIm);
  mpfr_add(mpIm, mpIa, mpIb, GMP_RNDN);
  mpfr_div_ui(mpIm, mpIm, 2, GMP_RNDN);

  if (!mpfr_cmp(mpIm, mpIa) || !mpfr_cmp(mpIm, mpIb))
    mpfr_set(mpX, mpIm, GMP_RNDN);
  else {
    mpfr_t mpFIa, mpFIm, mpFIb, mpTmp;

    mpfr_inits(mpFIa, mpFIm, mpFIb, mpTmp, NULL);

    f.mpEval(mpFIa, mpIa, n);
    mpP.eval(mpTmp, mpIa, n);
    mpfr_sub(mpFIa, mpFIa, mpTmp, GMP_RNDN);

    f.mpEval(mpFIb, mpIb, n);
    mpP.eval(mpTmp, mpIb, n);
    mpfr_sub(mpFIb, mpFIb, mpTmp, GMP_RNDN);

    mpfr_mul(mpTmp, mpFIa, mpFIb, GMP_RNDN);
    if (mpfr_sgn(mpTmp) > 0) {
      mpfr_clears(mpFIa, mpFIm, mpFIb, mpTmp, mpIm, NULL);
      throw "Minimax::findDicho: Unable to find a root in the given interval.";
    }

    f.mpEval(mpFIm, mpIm, n);
    mpP.eval(mpTmp, mpIm, n);
    mpfr_sub(mpFIm, mpFIm, mpTmp, GMP_RNDN);

    mpfr_mul(mpTmp, mpFIa, mpFIm, GMP_RNDN);
    if (mpfr_sgn(mpTmp) <= 0)
      findDicho(mpX, f, mpP, mpIa, mpIm, n);
    else
      findDicho(mpX, f, mpP, mpIm, mpIb, n);

    mpfr_clears(mpFIa, mpFIm, mpFIb, mpTmp, NULL);
  }

  mpfr_clear(mpIm);
}
