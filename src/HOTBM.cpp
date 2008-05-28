#include "HOTBM.hpp"
#include "HOTBM/Function.hh"
#include "HOTBM/Param.hh"
#include "HOTBM/HOTBMInstance.hh"
#include "HOTBM/Exhaustive.hh"
#include "utils.hpp"

HOTBM::HOTBM(Target* target, string func, int wI, int wO, int n)
	: inst(0), f(*new Function(func)), wI(wI), wO(wO)
{
	try {
		Param p(wI, wO, n);
		Exhaustive ex(f, p);

		inst = 0;
		inst = ex.getInstance();
		inst->roundTables();
	} catch (const char *s) {
		throw std::string(s);
	}
	if (!inst)
		throw std::string("HOTBM cound not be generated.");

	// TODO: Better unique name
	{
		std::ostringstream o;
		o << "hotbm_" << wI << "_" << wO << "_" << n;
		unique_name = o.str();
	}
	set_combinatorial();
	set_pipeline_depth(0);
	add_input("x", wI);
	add_output("r", wO+1);
}

HOTBM::~HOTBM()
{
	if (inst)
		delete inst;
}

// Overloading the virtual functions of Operator
void HOTBM::output_vhdl(std::ostream& o, std::string name)
{
	if (!inst)
		return;

	inst->genVHDL(o, name);
}

TestCaseList HOTBM::generateStandardTestCases(int n)
{
	// TODO
	return TestCaseList();
}

TestCaseList HOTBM::generateRandomTestCases(int n)
{
	Signal& sx = *get_signal_by_name("x");
	Signal& sr = *get_signal_by_name("r");

	mpz_class x, r;

	TestCaseList tcl;
	for (int i = 0; i < n; i++)
	{
		int outSign = 0;
		x = getLargeRandom(sx.width());

		mpfr_t mpX, mpR;
		mpfr_inits(mpX, mpR, 0);

		/* Convert a random signal to an mpfr_t in [0,1[ */
		mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
		mpfr_div_2si(mpX, mpX, wI, GMP_RNDN);

		/* Compute the function */
		f.eval(mpR, mpX);

		/* Compute the signal value */
		if (mpfr_signbit(mpR))
		{
			outSign = 1;
			mpfr_abs(mpR, mpR, GMP_RNDN);
		}
		mpfr_mul_2si(mpR, mpR, wO, GMP_RNDN);

		/* NOT A TYPO. HOTBM only guarantees faithful
		 * rounding, so we will round down here,
		 * add both the upper and lower neighbor.
		 */
		mpfr_get_z(r.get_mpz_t(), mpR, GMP_RNDD);

		TestCase tc;
		tc.addInput(sx, x);
		tc.addExpectedOutput(sr, r);
		tc.addExpectedOutput(sr, r+1);
		tcl.add(tc);
	}

	return tcl;
}

TestIOMap HOTBM::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("x"));
	tim.add(*get_signal_by_name("r"), 2); // Faithful rounding
	return tim;
}

void HOTBM::fillTestCase(mpz_class a[])
{
	/* Get inputs / outputs */
	mpz_class &x  = a[0];
	mpz_class &rd = a[1]; // rounded down
	mpz_class &ru = a[2]; // rounded up

	int outSign = 0;

	mpfr_t mpX, mpR;
	mpfr_inits(mpX, mpR, 0);

	/* Convert a random signal to an mpfr_t in [0,1[ */
	mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
	mpfr_div_2si(mpX, mpX, wI, GMP_RNDN);

	/* Compute the function */
	f.eval(mpR, mpX);

	/* Compute the signal value */
	if (mpfr_signbit(mpR))
	{
		outSign = 1;
		mpfr_abs(mpR, mpR, GMP_RNDN);
	}
	mpfr_mul_2si(mpR, mpR, wO, GMP_RNDN);

	/* NOT A TYPO. HOTBM only guarantees faithful
	 * rounding, so we will round down here,
	 * add both the upper and lower neighbor.
	 */
	mpfr_get_z(rd.get_mpz_t(), mpR, GMP_RNDD);
	ru = rd + 1;
}

