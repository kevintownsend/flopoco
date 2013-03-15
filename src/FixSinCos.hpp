#ifndef _FIXSINCOS_H
#define _FIXSINCOS_H

// works only with Sollya
#ifdef HAVE_SOLLYA

#include "Operator.hpp"

#include "utils.hpp"

using namespace flopoco;

class FixSinCos: public Operator {
	public:
		int w;


	public:
		FixSinCos(Target * target, int w);

		~FixSinCos();


		void emulate(TestCase * tc);

		void buildStandardTestCases(TestCaseList * tcl);


	private:
		mpfr_t scale;              /**< 1-2^(wOut-1)*/
		mpfr_t constPi;

};

#endif // HAVE_SOLLYA

#endif // header guard

