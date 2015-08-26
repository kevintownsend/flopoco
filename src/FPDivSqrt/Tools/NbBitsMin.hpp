#include "../../UserInterface.hpp"

namespace flopoco
{
	//This file computes the minimum number of bits required as LUT's input (without prescaling)
	void plotPDDiagram(int delta, int t, int radix, int digitSet);
	bool checkDistrib(int delta, int t, int radix, int digitSet);
	float L(int k, float ro, float d);
	float U(int k, float ro, float d);
	float estimateCost(int nbBit, int radix, int digitSet);
	void computeNbBit(int radix, int digitSet);
	void NbBitsMinRegisterFactory();
}

