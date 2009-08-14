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
	void initTiling(DSP** &config, int dspCount);
	
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
	
	/** The vector which will mark if a dsp in the tiling algorithm was or not rotated */
	bool* rot;

	/* The number of estimated DSPs that will be used according to this parameter */
	int nrDSPs;

	/** This function estimates the maximum number of DSPs that can be used with respect to the preference of the user */
	int estimateDSPs();

	/** This function computes the cost of the configuration received as input parameter */
	
	float computeCost(DSP** &config);

	/** This function compares the cost of the current configuration( globalConfig) with the best configuration obtained so far(bestConfig). If the current one is better then the best one, 
	then it will copy the current configuration into the best one.
	*/
	void compareCost();
	
		/** This function is the backtracking function for finding the best configuration with respect to the user preference
	*/
	void tilingAlgorithm(int i, int n, bool repl);
	
	/** This function will take a configuration and will try to maximize the multiplications that are realized in slices 
	* It will return the number of slices that are required to ferform those multiplications, and through the parameter partitions it will return the number of such partitions
	*/
	int partitionOfGridSlices(DSP** config,int &partitions);
	 
	 /** This function will fill the input matrix between the limits received as parameters with the value received as parameter */
	 
	 void fillMatrix(int **&matrix,int lw,int lh,int topleftX,int topleftY,int botomrightX,int botomrightY,int value);
	
	/** This function resets all the conections that exists between DSPs of a configuration */
	
	void resetConnections(DSP** &config);
	
	/** This function will create connections between the existing DSPs in the configuration according to the policies of Altera or Virtex */
	
	int bindDSPs(DSP** &config);
	
	/** This function will create the shifts between DSPs for Virtex boards */
	
	int bindDSPs4Virtex(DSP** &config);
	
	/** This function will sort the DSPs blocks by the topright corner */
	
	void sortDSPs(DSP** &config);

};


#endif
