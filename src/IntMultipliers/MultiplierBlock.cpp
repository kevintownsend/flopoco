#include "MultiplierBlock.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	
#include <cstdlib>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "../Operator.hpp"
#include "IntMultiplier.hpp"
#include "IntAdder.hpp"
#include "IntMultiAdder.hpp"
#include "IntAddition/NewCompressorTree.hpp"
#include "IntAddition/PopCount.hpp"
#include "utils.hpp"
#include<vector>
#include<list>
using namespace std;


namespace flopoco{


		MultiplierBlock::MultiplierBlock(int wX, int wY, int topX, int topY, string input1, string input2, int weightShift_, int cycle) :
			wX(wX), wY(wY), topX(topX), topY(topY), weightShift(weightShift_), inputName1(input1), inputName2(input2)
		{
			srcFileName=":MultiplierBlock"; 
			cycle=cycle;
			weight=topX+topY-weightShift;
			previous=NULL;
			next=NULL;
			
			
		}


		bool MultiplierBlock::operator <= (MultiplierBlock* b){
		if ((weight<=b->getWeight())) 
			return true;
		else
			return false;
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


		int MultiplierBlock::getbotX()
		{
			return topX+wX;
		}

		int MultiplierBlock::getbotY()
		{
			return topY+wY;
		}


		int MultiplierBlock::getCycle()
		{
			return cycle;
		}
	

		void MultiplierBlock::setNext(MultiplierBlock* b)
		{
			this->next=b;
		}

		MultiplierBlock* MultiplierBlock::getNext()
		{
			if(next!=NULL)
				return next;
			else
				return NULL;
		}

		

		MultiplierBlock* MultiplierBlock::getPrevious()
		{
			if(previous!=NULL)
				return previous;
			else
				return NULL;
		}

		void MultiplierBlock::setPrevious(MultiplierBlock* b)
		{
			this->previous=b;
		}


				//TODO improve the chaining
				// - pass the full target
				// - use the "double multiplier mode" on Altera
				// ...

		bool MultiplierBlock::canBeChained(MultiplierBlock* next, bool isXilinx)
		{
			//for now just the stupid chaining
			if (isXilinx)
			{
				if(neighbors(next))
					return true;
				else
					return false;
			}
			//Altera chaining
			else
			{
				if((next==NULL) && (previous==NULL) && (next->next==NULL) && (next->previous==NULL) 
						&& (wX==next->wX) && (wY==next->wY))
				{
					if((wX>18) || (wY>18))
						return false;
					else
					{
						if((getWeight() == next->getWeight()) && (this!=next))
							return true;
						else
							return false;
					}
				}
				else 
					return false;
			}

		}

		bool MultiplierBlock::neighbors(MultiplierBlock* next)
		{
			if((((this->topX==next->gettopX()) &&
					(this->topY==next->gettopY()+17)) ||
				((this->topY==next->gettopY()) &&
					(this->topX==next->gettopX()+17)) )  ||
				
				(((this->topX==next->gettopX()) &&
					(this->topY==next->gettopY()-17)) ||
				((this->topY==next->gettopY()) &&
					(this->topX==next->gettopX()-17)) )) 
				return true;
			else
				return false;

		}



	



     

}
