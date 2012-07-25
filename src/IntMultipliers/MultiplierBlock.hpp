#ifndef MultiplierBlock_HPP
#define MultiplierBlock_HPP
#include <stdio.h>
#include <stdlib.h>
#include "../Operator.hpp"

namespace flopoco {

	//FIXME: please comment your code!

	class MultiplierBlock 
	{
	public: 
	
		/**
		 * The default constructor
		 */
		MultiplierBlock(Operator* op, int wX, int wY, int topX, int topY, bool goToDSP=false,int cycle=-1, MultiplierBlock* previous=NULL, MultiplierBlock* next=NULL);
	
		
		/**
		 * Set the value of signalName
		 * @param the new value
		 */
		void setSignalName(string name);

		void setSignalLength(int length);

		string getSigName();

		int getSigLength();

		int getwX();

		int getwY();

		int gettopX();

		int gettopY();

		int getCycle();
	
		bool getgoToDSP();

		MultiplierBlock* getNext();

		MultiplierBlock* getPrevious();

		
	private:
	
		Operator* op;
		int wX; 							/**< x size */
		int wY; 							/**< y size */
		int topX; 							/**< x position (top right corner) */
		int topY; 							/**< y position (top right corner */
		bool goToDSP; 						/**< a bit saying if it should go into a DSP */
		int cycle;							/**< cycle */
		MultiplierBlock* previous;
		MultiplierBlock* next;
		string signalName;
		int signalLength;
	};

}
#endif
