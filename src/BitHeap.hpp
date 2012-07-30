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
#include "IntMultipliers/MultiplierBlock.hpp"



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

            /** returns the stage when this bit should be compressed */ 
            int computeStage(int stagesPerCycle, double elementaryTime);


			

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

		//adds a new MultiplierBlock in the list he already has
		void  addDSP(MultiplierBlock* m);

		
        void elemReduce(unsigned i, BasicCompressor* bc);
		void iterateDSP();

		//** computes the latest bit from a column, in order to compress just the bits which are smaller than that one*/
        BitHeap::WeightedBit* computeLatest(unsigned w, int c0, int c1);
        
		//** computes the latest bit from the bitheap, in order to manage the cycle before the final adding*/
        BitHeap::WeightedBit* getFinalLatestBit();

        BitHeap::WeightedBit* getFirstSoonestBit();

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

		

		/** search for the possible chainings and generates the VHDL code too*/
		 void doChaining();

        /**is making the compression for the bitheap**/
		void compress(int stage);
		
		/** return the current height a column (bits not yet compressed) */
		unsigned currentHeight(unsigned w);

		int getUid(unsigned w);

		/** counts the bits not processed yet in wb */
		int count(list<WeightedBit*> wb, int cycle);
		
		void printColumnInfo(int w);
		
		void generatePossibleCompressors();

		/** remove the compressed bits */
		void removeCompressedBits(int c, int red);

        /** returns the maximum height list from the bitheap vector*/
		unsigned getMaxHeight();

		void getMaxWeight();

		//generate the VHDL code for 1 dsp, 
		void generateVHDLforDSP(MultiplierBlock* m, int uid,int i);

 		void initializeDrawing();

        void closeDrawing(int offsetY);

        void drawConfiguration(int offsetY);

        void drawBit(int cnt, int w, int turnaroundX, int offsetY, int c);

	private:
		Operator* op;
		unsigned maxWeight;     /**< The compressor tree will produce a result for weights < maxWeight (work modulo 2^maxWeight)*/
		unsigned minWeight; /**< bits smaller than this one are already compressed */    
		bool usedCompressors[10]; /** the list of compressors which were used at least once*/
		unsigned chunkDoneIndex; 
		unsigned inConcatIndex; /** input index - to form the inputsignals of the compressor*/
		unsigned outConcatIndex; /** output index - to form the outputsignals of the compressor*/
		unsigned compressorIndex; /** the index of the instance of compressors*/
        unsigned cnt[100000]; /** number of bits which will be compressed in the current iteration*/
		vector<int> uid;   /**< unique id, per weight */
     	ofstream fileFig;
        ostringstream fig;
        bool drawCycleLine;
        int drawCycleNumber;    
        int stagesPerCycle;  
        double elementaryTime; 
        bool didCompress; 
#if 0
		const static int consumed=-1;
		vector<vector<int> > cycle;   /**< external index is the weight (column). The int is the cycle of each bit. Consumed bits have their cycle set to -1 */
		vector<vector<double> > cpd;  /**< external index is the weight (column). The double is the critical path delay of each bit */
		string bit( int w, int h); /**< just provide the name of the bit of weight w and height h*/
#else
		
		vector<list<WeightedBit*> > bits; /**<  The list is ordered by arrival time of the bits, i.e. lexicographic order on (cycle, cp)*/
		vector<MultiplierBlock*> mulBlocks; //the vector of multiplier blocks
#endif
     
		string srcFileName;	};


}
#endif
