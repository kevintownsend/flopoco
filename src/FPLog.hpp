#ifndef __FPLOG_HPP
#define __FPLOG_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

/* Opaque! We don't want them to be known through the interface */
class FirstInvTable;
class FirstLogTable;
class SecondInvTable;
class OtherLogTable;
class FPNumber;

class FPLog : public Operator
{
public:
	FPLog(Target* target, int wE, int wF);
	~FPLog();

	// Overloading the virtual functions of Operator
	void outputVHDL(std::ostream& o, std::string name);

	TestIOMap getTestIOMap();
	void fillTestCase(mpz_class a[]);

private:
	int wE, wF;
	// The input sizes to the successive tables
	int a[42]; 
	// The intermediate precision: at step i, the exp part is bounded by
	//    1 <= m < 1+2^-p[i]
	int p[42]; 

	// The size of the product, and hence of the subword of Z[i] input to the mult
	int psize[42]; 

	// The total size of non-zero bits, will be limited to wF+g
	int s[42]; 

	// The numbers of bits to truncate to limit the size to wF+g
	int t[42];

	// size before truncation, should be equal to s[i]+t[i]
	int sbt[42];

	// The max number of stages
	int stages;
	
	// The guard bits 
	int gLog;

	// The table objects
	FirstInvTable* it0;
	SecondInvTable* it1;
	FirstLogTable* lt0;
	OtherLogTable* lt[42];

	// This is the size of the virtual mantissa:
	// 1.000(p[i] zeroes)00000xxxxxxx
	// the xxx which will be actually computed have p[i] bits less 
	int sfullZ[42]; 

	// A global boolean flag disabling the simulation of the fullsize mantissa
	// as soon as it would be more than 64 bits
	int fullSizeSim;

	// The target precision: numbers may be truncated so that their LSB has weight -target_prec 
	int target_prec;  
};

#endif
