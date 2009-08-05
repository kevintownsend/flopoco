#ifndef IntTillingMult_HPP
#define InTtillingMult_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

/** The InttilingMult class */
class IntTilingMult : public Operator
{
public:
	/**
	 * The FPAdder constructor
	 * @param[in]		target		the target device
	 * @param[in]		MSB			the MSB of the input number
	 * @param[in]		LSB			the LSB of the input number
	 * @param[in]		wER			the with of the exponent for the convertion result
	 * @param[in]		wFR			the with of the fraction for the convertion result
	 */
	IntTilingMult(Target* target, int wInX, int wInY,float ratio);

	/**
	 * IntTilingMult destructor
	 */
	~IntTilingMult();


	void emulate(TestCase * tc);
	void buildStandardTestCases(TestCaseList* tcl);

private:
	/**
	 * Verifies that the indicated DSP block does not overlap with any 
	 * other block in the given configuration.
	 * @param config the tiling configuration containing DSP objects
	 * @param index the index of the DSP object we want to check is not overlapping in the configuration
	 * @return TRUE if no DSP blocks overlap on the tiling configuration
	 */
	bool checkOverlap(DSP** config, int index);
	
	/**
	 * Displaces the DSP block in the 'next' valid configuration, i.e.
	 * The DSP block in the configuration, indicated by the index parameter
	 * is moved in downward and to the left.
	 * @param config the tiling configuration containing DSP objects
	 * @param index the index of the DSP object we want to move within the configuration
	 * @return FALSE if we have reached the bottom-left most corner of the tiling
	 */
	bool move(DSP** config, int index);
	
	/**
	 * Repositions the indicated DSP object within the configuration in the 
	 * top-right most position possible.
	 * @param config the tiling configuration containing DSP objects
	 * @param index the index of the DSP object we want to move within the configuration
	 */
	void replace(DSP** config, int index);
	
	/**
	 * Initializes the tiling grid with the DSP-s being positioned in the right-most
	 * and top-most position possible without any overlapping of the DSP blocks
	 * @param config the tiling configuration containing DSP objects
	 * @param dspCount the number of DSP blocks in the configuration
	 */
	void initTiling(DSP** config, int dspCount);
	
	/**
	* This function returns the amount of extra displacemnet on the tiling board that a dsp block can have on vertical axis
	* @return extra height
	*/
	int getExtraHeight();

	/**
	* This function returns the amount of extra displacemnet on the tiling board that a dsp block can have on horizonta axis
	* @return extra width
	*/
	int getExtraWidth();
	

	/** The width of first input*/
	int wInX; 
	/** The width of second input*/
	int wInY; 
	/** The width of output */
	int wOut;
	/** The ratio between the number of DSPs and slices */
	float ratio;
	/* The working configuration of the tiling algorithm on DSPs */
	DSP** globalConfig;
	/** The best cost obtained so far*/
	float bestCost;
	/* The best configuration of the after tiling DSPs */
	DSP** bestConfig;
	/* The target */
	Target * target;

	/* The number of estimated DSPs that will be used according to this parameter */
	int nrDSPs;

	/** This function estimates the maximum number of DSPs that can be used with respect to the preference of the user */
	int estimateDSPs();

	/** This function computes the cost of the configuration received as input parameter */
	
	float computeCost(DSP** config);

	/** This function compares the cost of the current configuration( globalConfig) with the best configuration obtained so far(bestConfig). If the current one is better then the best one, 
	then it will copy the current configuration into the best one.
	*/
	void compareCost();
	
	/** This function is the backtracking function for finding the best configuration with respect to the user preference
	*/
	void tilingAlgorithm(int i, int n, bool rot,bool repl);
	


};


#endif
