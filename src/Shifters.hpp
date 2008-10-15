#ifndef SHIFTERS_HPP
#define SHIFTERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"

#include "Operator.hpp"

extern map<string, double> emptyDelayMap;
/** The types of shifting */
typedef enum {
			Left, /**< Left Shifter */
			Right /**< Right Shifter */
			} ShiftDirection;

/** The Shifter class. Left and right shifters are perfectly
		symmetrical, so both are instances of the Shifter class. Only the
		name of the VHDL instance changes */
class Shifter : public Operator
{
public:
	/**
	 * The Shifter constructor
	 * @param[in]		target		the target device
	 * @param[in]		wIn			  the with of the input
	 * @param[in]		maxShift	the maximum shift ammount
	 * @param[in]		direction	can be either Left of Right. Determines the shift direction
	 **/
	Shifter(Target* target, int wIn, int maxShift, ShiftDirection dir, map<string, double> inputDelays = emptyDelayMap);

	/** Destructor */
	~Shifter();

	/**
	 * Method belonging to the Operator class overloaded by the Shifter class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);

	/** 
	 * Sets the default name of this operator
	 */
	void setOperatorName(); 

	/**
	 * Gets the signals which are interesting for TestCases.
	 * @see TestIOMap
	 */
	TestIOMap getTestIOMap();
	
	/**
	 * Gets the correct value associated to one or more inputs.
	 * @param a the array which contains both already filled inputs and
	 *          to be filled outputs in the order specified in getTestIOMap.
	 */
	void fillTestCase(mpz_class a[]);

	/** Returns the number of bits of the sift ammount */
	int getShiftAmmount(){
		return wShiftIn_;
	}
protected:
	int wIn_;          /**< the width of the input*/
	int maxShift_;     /**< the maximum shift ammount*/
	int wOut_;         /**< the width of the output */
	int wShiftIn_; 	   /**< the number of bits of the input which determines the shift ammount*/
	map<string, double> inputDelays_;

private:
	ShiftDirection direction_;            /**< determines the shift direction. can be Left or Right */
	bool           levelRegistered_[128]; /**< if boolean true, the corresponding level signal is registered*/ 
	double maxInputDelay;
};
#endif
