/*
  A class to manage heaps of weighted bits in FloPoCo
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2012.
  All rights reserved.

*/
#ifndef __BITHEAP_HPP
#define __BITHEAP_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "Table.hpp"
#include "DualTable.hpp"
#include "IntAddition/BasicCompressor.hpp"
#include "IntMultipliers/IntTilingMult.hpp"



/* 
Each bit in the bit heap is flagged with the cycle at which it is produced.
Compression works as follows:

setCycle(0)
while needed
  manageCriticalPath() to advance the cycle  
  consider the subset of bits with cycles 0 to getCycle()
  build the longest possible heap of size 2 or less from the right,
     and feed them to an adder.
  compress the others bits greedily, 
    adding more bits, flagged with the current cycle, 
    setting the cycle of the consumed bits to -1 so they won't be considered in subsequent iterations
  
*/

namespace flopoco{





	class BitHeap 
	{

		class WeightedBit
		{
		public:
			/** constructor 
			    @param bh the parent bit heap */ 
			WeightedBit(BitHeap* bh, int weight, int cycle=-1,  double criticalPath=-1.0);

			/** destructor */ 
			~WeightedBit(){};
		
			/** mark this bit as already compressed */
			void done(){
				processed=true;
			};

			/** return the cycle at which this bit is defined */
			int getCycle(){
				return cycle;
			};

			/** return the critical path of this bit */
			double getCriticalPath(int atCycle);

			/** return the VHDL signal name of this bit */
			string getName(){
				return name;
			};

			/** returns true if this bit is still to compress, false if this bit was already compressed */
			bool todo(){
				return (!processed);
			};
		
			/** ordering by availability time */
			bool operator< (const WeightedBit& b); 
			/** ordering by availability time */
			bool operator<= (const WeightedBit& b); 


		private:
			int cycle;  /**< The cycle at which the bit was created */
			double criticalPath; /**< The cycle at which the bit was created */
			BitHeap* bh;
			bool processed;
			int weight;
			int uid;
			string name;
 			string srcFileName;
 			string uniqueName_;
		};

	public:

		/** The constructor
		 @param op         the operator in which this bit heap is beeing built
		 @param maxWeight  the maximum weight of the heap (it should be known statically, shouldn't it?) */
		BitHeap(Operator* op, int maxWeight);
		~BitHeap();
		
		/** add a bit to the bit heap. The bit will be added at the cycle op->currentCycle() with critical path op->getCriticalPath().
		 @param weight   the weight of the bit to be added
		 @param rhs      the right-hand VHDL side defining this bit.
		 @param comment  a VHDL comment for this bit*/
		void addBit(unsigned weight, string rhs, string comment="");

        void elemReduce(int i, BasicCompressor* bc);


        /** remove a bit from the bitheap.
         @param weight  the weight of the bit to be removed
         @param dir if dir==0 the bit will be removed from the begining of the list 
         			if dir==1 the bit will be removed from the end of the list
        */
        void removeBit(unsigned weight, int dir);

		/** get the parent operator */
		Operator* getOp() {return op;};

		/** generate the VHDL for the bit heap. This will involve adding more bits for intermediate results.*/
		void generateCompressorVHDL();

        /** generate the final adder for the bit heap (when the columns height is maximum 2*/
        void adderVHDL();


        /**is making the compression for the bitheap**/
		void compress();
		
		/** return the current height a column (bits not yet compressed) */
		unsigned currentHeight(unsigned w);

		int getUid(unsigned w);

		/** counts the bits not processed yet in wb */
		int count(list<WeightedBit*> wb, int cycle);

		/** marks the compressed bits as done*/
		void reduce(int c, int red);

        /** returns the maximum height list from the bitheap vector*/
		int getMaxHeight();

		void getMaxWeight();

 		string getResultVHDLName();
	private:
		Operator* op;
		vector<int> uid;   /**< unique id, per weight */
#if 0
		const static int consumed=-1;
		vector<vector<int> > cycle;   /**< external index is the weight (column). The int is the cycle of each bit. Consumed bits have their cycle set to -1 */
		vector<vector<double> > cpd;  /**< external index is the weight (column). The double is the critical path delay of each bit */
		string bit( int w, int h); /**< just provide the name of the bit of weight w and height h*/
#else
		vector<list<WeightedBit*> > bits; /**<  The list is ordered by arrival time of the bits, i.e. lexicographic order on (cycle, cp)*/
#endif
		unsigned maxWeight;     /**< The compressor tree will produce a result for weights < maxWeight (work modulo 2^maxWeight)*/
		unsigned minWeight; /**< bits smaller than this one are already compressed */         
		string srcFileName;	};


}
#endif
