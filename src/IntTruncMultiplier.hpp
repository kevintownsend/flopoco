#ifndef IntTruncMultiplier_HPP
#define IntTruncMultiplier_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntAdder.hpp"
#include "IntNAdder.hpp"


namespace flopoco{

	class SoftDSP{
	public:
		SoftDSP(int topX, int topY, int bottomX, int bottomY):
			topX(topX), topY(topY), bottomX(bottomX), bottomY(bottomY){
			/*constructor*/
			
		}	
		
		~SoftDSP(){}
		
		void getTopRightCorner(int &xT, int &yT){
			xT = topX;
			yT = topY;
		}
		
		void getBottomLeftCorner(int &xB, int &yB){
			xB = bottomX;
			yB = bottomY;
		}
	
	int topX, topY, bottomX, bottomY;
	};


	class IntTruncMultiplier : public Operator
	{
	public:
		
		IntTruncMultiplier(Target* target, DSP** configuration, vector<SoftDSP*> softDSPs, int wX, int wY, int k);

		IntTruncMultiplier(Target* target,int wX, int wY,float ratio, int k);
		IntTruncMultiplier(Target* target,int wX,float ratio, int k);
	
		/** IntTruncMultiplier destructor */
		~IntTruncMultiplier();

	protected:

		int wX; /**< the width (in bits) of the input  X  */
		int wY; /**< the width (in bits) of the input  Y  */
		int wt; /**< the width (in bits) of the output R  */
		float ratio; /**<the ratio between slice and dsp */


		void printConfiguration(DSP** configuration, vector<SoftDSP*> softDSPs);
	private:

	};

}
#endif
