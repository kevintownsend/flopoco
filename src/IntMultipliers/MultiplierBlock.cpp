#include "MultiplierBlock.hpp"
using namespace std;


namespace flopoco{


		MultiplierBlock::MultiplierBlock(Operator* op, int wX, int wY, int topX, int topY, bool goToDSP,int cycle, MultiplierBlock* previous, MultiplierBlock* next) :
			op(op), wX(wX), wY(wY), topX(topX), topY(topY)
		{
			if(cycle==-1)
				cycle = op->getCurrentCycle();
			else
				cycle=cycle;
		}

		void MultiplierBlock::setSignalName(string name)
		{
			signalName=name;
		}

		void MultiplierBlock::setSignalLength(int length)
		{
			signalLength=length;
		}

		string MultiplierBlock::getSigName()
		{
			return signalName;
		}

		int MultiplierBlock::getSigLength()
		{
			return signalLength;
		}

		int MultiplierBlock::getwX()
		{
			return wX;
		}

		int MultiplierBlock::getwY()
		{
			return wY;
		}

		int MultiplierBlock::gettopX()
		{
			return topX;
		}

		int MultiplierBlock::gettopY()
		{
			return topY;
		}

		int MultiplierBlock::getCycle()
		{
			return cycle;
		}
	
		bool MultiplierBlock::getgoToDSP()
		{
			return goToDSP;
		}

		MultiplierBlock* MultiplierBlock::getNext()
		{
			return next;
		}

		MultiplierBlock* MultiplierBlock::getPrevious()
		{
			return previous;
		}
}
