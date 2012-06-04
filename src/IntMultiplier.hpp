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
			int dx, dy; 			
			SmallMultTable(Target* target, int dx, int dy );
			mpz_class function(int x);
		};

    /**
    * The IntMultiplier constructor
    * @param[in] target           the target device
    * @param[in] wX             X multiplier size (including sign bit if any)
    * @param[in] wY             Y multiplier size (including sign bit if any)
    * @param[in] wTruncated       return result faithful to wX+wY-wTruncated bits (0 means full multiplier)
    * @param[in] signedInputs     false=unsigned, true=signed
    * @param[in] ratio            DSP block use ratio
    **/
	IntMultiplier(Target* target, int wX, int wY, int wOut=0, bool signedInputs = false, float ratio = 1.0, bool useDSP=true, map<string, double> inputDelays = emptyDelayMap);

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
	
	int wX;                         /**< the width for X*/
	int wY;                         /**< the width for Y*/
	int wOut;
	int wTruncated;
	int signedInputs;
	float ratio;
	int g ; /**< the number of guard bits*/
	
private:
	bool useDSP;
};

}
#endif
