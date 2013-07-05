#ifndef _FIXSINCOS_H
#define _FIXSINCOS_H

// works only with Sollya
#ifdef HAVE_SOLLYA

#include <gmpxx.h>

#include "Operator.hpp"

#include "utils.hpp"

using namespace flopoco;




class FixSinCos: public Operator {
public:
	class SinCosTable: public Table {
	public:
		SinCosTable(Target* target, int a, int g, int argRedCase_, FixSinCos* parent);
		~SinCosTable();
	private:
		mpz_class function(int x);
		int argRedCase;     /**< argRedCase=1 for full table, 4 for quadrant */
		FixSinCos* parent;
	};


public:
	
	FixSinCos(Target * target, int w, float ratio=0.5);
	
	~FixSinCos();
	
	
	void emulate(TestCase * tc);
	
	void buildStandardTestCases(TestCaseList * tcl);
	

private:
	int w;
	mpfr_t scale;              /**< 1-2^(wOut-1)*/
	mpfr_t constPi;
};

#endif // HAVE_SOLLYA

#endif // header guard

