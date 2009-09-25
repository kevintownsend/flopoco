#ifndef SRT4STEP_HPP
#define SRT4STEP_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../Operator.hpp"


namespace flopoco{

	/** The SRT4Step class */
	class SRT4Step : public Operator
	{
	public:
		/**
		 * The SRT4Step constructor
		 * @param[in]		target		the target device
		 * @param[in]		wF			the the with of the fraction for the f-p number X
		 */
		SRT4Step(Target* target, int wF);

		/**
		 * SRT4Step destructor
		 */
		~SRT4Step();

		/**
		 * Method belonging to the Operator class overloaded by the SRT4Step class
		 * @param[in,out] o     the stream where the current architecture will be outputed to
		 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
		 **/
		void outputVHDL(std::ostream& o, std::string name);


	
	private:
		/** The width of the fraction for the input X */
		int wF; 

	};
}
#endif //SRT4STEP_HPP
