#ifndef IntTillingMult_HPP
#define IntTillingMult_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"

/** The IntTillingMult class */
class IntTillingMult : public Operator
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
	IntTillingMult(Target* target, int wInX, int wInY,float ratio);

	/**
	 * IntTillingMult destructor
	 */
	~IntTillingMult();


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
	* This function returns the amount of extra displacemnet on the tilling board that a dsp block can have on vertical axis
	* @return extra height
	*/
	int getExtraHeight();

	/**
	* This function returns the amount of extra displacemnet on the tilling board that a dsp block can have on horizonta axis
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
	/* The working configuration of the tilling algorithm on DSPs */
	DSP** globalConfig;
	/* The best configuration of the after tilling DSPs */
	DSP** bestConfig;
	/* The target */
	Target * target;

	/* The number of estimated DSPs that will be used according to this parameter */
	int nrDSPs;

	int estimateDSPs();


	void tillingAlgorithm();
	


};


#endif
