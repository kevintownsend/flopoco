#ifndef IntMultiplierS_HPP
#define IntMultiplierS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "../Table.hpp"

namespace flopoco {


		


/** The IntMultiplier class, getting rid of Bogdan's mess.
*/
class IntMultiplier : public Operator {
public:
		/** An elementary LUT-based multiplier, written as a table so that synthesis tools don't infer DSP blocks for it*/
		class SmallMultTable: public Table {
		public:
			int dx, dy, wO;
			bool signedIO;
			SmallMultTable(Target* target, int dx, int dy, int wO, bool signedIO=false );
			mpz_class function(int x);
		};

    /**
    * The IntMultiplier constructor
    * @param[in] target           the target device
    * @param[in] wX             X multiplier size (including sign bit if any)
    * @param[in] wY             Y multiplier size (including sign bit if any)
    * @param[in] wTruncated       return result faithful to wX+wY-wTruncated bits (0 means full multiplier)
    * @param[in] signedIO     false=unsigned, true=signed
    * @param[in] ratio            DSP block use ratio
    **/
	IntMultiplier(Target* target, int wX, int wY, int wOut=0, bool signedIO = false, float ratio = 1.0, map<string, double> inputDelays = emptyDelayMap);

    /**
    *  Destructor
    */
    ~IntMultiplier();

    /**
    * The emulate function.
    * @param[in] tc               a list of test-cases
    */
    void emulate ( TestCase* tc );

		void buildStandardTestCases(TestCaseList* tcl);



protected:

	string PP(int i, int j);
	string PPTbl(int i, int j);
	string XY(int i, int j);
	string heap( int i, int j);

	void buildLogicOnly();


	void buildTiling();
	
	int wXdecl;                     /**< the width for X as declared*/
	int wYdecl;                     /**< the width for Y  as declared*/
	int wX;                         /**< the width for X after possible swap such that wX>wY */
	int wY;                         /**< the width for Y after possible swap such that wX>wY */
	int wOut;
	int wFull;                      /**< size of the full product: wX+wY-1 if signed, wX+wY if unsigned */
	int wTruncated;
	int signedIO;
	double ratio;
	double maxError; /**< the max absolute value error of this multiplier, in ulps of the result. Should be 0 for untruncated, 1 or a bit less for truncated.*/  
	int g ; /**< the number of guard bits*/
	
private:
	bool useDSP;
};

}
#endif
