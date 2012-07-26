#include "Tiling.hpp"
using namespace std;


namespace flopoco{


		Tiling::Tiling(Target* target_, int wX_, int wY_, double ratio_, bool truncated_, int truncationSize_) :
			target(target_), wX(wX_), wY(wY_), ratio(ratio_), truncated(truncated_), truncationSize(truncationSize_)
		{
			if(truncated){
				xOrigin = wX-1;
				yOrigin = wY-1;
				
				xIncrement = -1;
				yIncrement = -1;
			}else{
				xOrigin = 0;
				yOrigin = 0;
				
				xIncrement = 1;
				yIncrement = 1;
			}
			
			xBorder = wX - xOrigin - 1;
			yBorder = wY - yOrigin - 1;
		}
		
		Tiling::~Tiling();

		vector<MultiplierBlock*> Tiling::createTiling()
		{
			vector<MultiplierBlock*> hTiling, vTiling, mixedTiling, bestConfiguration;
			
			hTiling = createHorizontalTiling();
			vTiling = createVerticalTiling();
			mixedTiling = createMixedTiling();
			
			bestConfiguration = hTiling;
			if(computeTilingCost(vTiling) > computeTilingCost(bestConfiguration))
				bestConfiguration = vTiling;
			if(computeTilingCost(mixedTiling) > computeTilingCost(bestConfiguration))
				bestConfiguration = mixedTiling;
				
			return bestConfiguration;
		}
		
		vector<MultiplierBlock*> Tiling::createHorizontalTiling()
		{
			vector<MultiplierBlock*> result;
			int xCurrent, yCurrent;
			bool tilingIncomplete = true;
			
			while(tilingIncomplete)
			{
				while(xCurrent < xBorder){
					// continue here
				}
			}
		}
		
		vector<MultiplierBlock*> Tiling::createVerticalTiling()
		{
			cerr << "Function createVerticalTiling not yet implemented. Exiting." << endl;
			exit(1);
		}
		
		vector<MultiplierBlock*> Tiling::createMixedTiling()
		{
			cerr << "Function createMixedTiling not yet implemented. Exiting." << endl;
			exit(1);
		}
		
		
		double Tiling::computeTilingCost(vector<MultiplierBlock*> configuration)
		{
			cerr << "Function computeTilingCost not yet implemented. Exiting." << endl;
			exit(1);
		}
		
		
		bool Tiling::validateTiling(vector<MultiplierBlock*> configuration)
		{
			cerr << "Function validateTiling not yet implemented. Exiting." << endl;
			exit(1);
		}
}
