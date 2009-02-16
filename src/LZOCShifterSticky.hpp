#ifndef LZOCShifterSticky_HPP
#define LZOCShifterSticky_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
extern map<string, double> emptyDelayMap;
/** 
 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
 */ 
class LZOCShifterSticky : public Operator
{
public:
	/** The different entity types for this operator */
	typedef enum {
				gen, /**< generic entity type. The OZB: in std_logic is present in entity */
				spec /**< specific entity type. The OZB: in std_logic is NOT present in entity */
				} entityType_t;

	/** 
	 * The LZOCShifterSticky constructor
	 * @param[in] target the target device for this operator
	 * @param[in] wIn the width of the mantissa input
	 * @param[in] wOut the width of the mantissa output
	 */
	LZOCShifterSticky(Target* target, int wIn, int wOut, bool compute_sticky, const int countType=-1, map<string, double> inputDelays = emptyDelayMap);
	
	/** The LZOCShifterSticky destructor */
	~LZOCShifterSticky();

	/** Sets the type of the entity for this operator
	 * @param eType the type of the entity for this operator. Can be
	 *              generic or specific
	 */
	void setEntityType(entityType_t eType);
	
	/** Returns the number of bits of the count
	 * @return the number of bits of the count
	 */
	int getCountWidth() const;
	
	
	/**
	 * Method belonging to the Operator class overloaded by the LZOCShifterSticky class
	 * @param[in,out] o     the stream where the current architecture will be outputed to
	 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
	 **/
	void outputVHDL(std::ostream& o, std::string name);
	
	/** Method for setting the operator name
	*/
	void setOperatorName();
	
	
	void emulate(TestCase* tc);
	
	double compDelay(int n){
		double fdcq = 0.826e-9;
		double lutd = 0.147e-9;
		double muxcy = 0.034e-9;
		double fdd = 0.017e-9;
		
		return (fdcq + lutd + (muxcy*n/2) + fdd);
	};
	double muxDelay(){
		double fdcq = 0.826e-9;
		double lutd = 0.147e-9;
		double fdd = 0.017e-9;
		
		return (fdcq + lutd + fdd);
	}


private:
	
	int          wIn_;                   /**< The number of bits of the input */
	int          wCount_;                /**< The number of bits of the count */
	int          wOut_;                  /**< The number of bits of the shifted output */
	bool         computeSticky_;         /**< If true, compute the sticky bit. If false, save this hardware */
	int          countType_;             /**< -1|0|1. If -1 is present then generic LZOC is instatiated */
	string       level_[42];             /**< The names of the signals, just to make code more readable */ 
	string       leveld_[42];            /**< Same but possibly delayed  */
	int          size_[42];              /**< Their size. Do we need to count more than 2^42 bits in FloPoCo? */      
	entityType_t entityType_;            /**< Entity type. Can be either generic or specific */
	bool         levelRegistered_ [42]; /**< if boolean true, the corresponding level signal is registered*/ 
	int          countDepth_[42];        /**< the depths for the levels of the architecture */	
	mpz_class    maxValue_;              /**< utilitary var */
};

#endif
