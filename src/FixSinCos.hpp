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
		FixSinCos(Target * target, int w = 2//, int unused = 0);
		);

		~FixSinCos() {
		};


		void emulate(TestCase * tc);

		void buildStandardTestCases(TestCaseList * tcl);

		void buildRandomTestCases(TestCaseList * tcl, int n);

		TestCase *buildRandomTestCases(int i);
};

#endif // HAVE_SOLLYA

#endif // header guard

