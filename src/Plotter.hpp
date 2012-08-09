/*
  A class used for plotting various drawings in SVG format
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Kinga Illyes, Bogdan Popa

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2012.
  All rights reserved.

*/

#ifndef __PLOTTER_HPP
#define __PLOTTER_HPP

#include <vector>
#include <sstream>
#include "WeightedBit.hpp"
#include "IntMultipliers/MultiplierBlock.hpp"
//#include "BitHeap.hpp"


namespace flopoco
{

	class BitHeap;

	class Plotter
	{

		class Snapshot
		{
			public:

				Snapshot(vector<list<WeightedBit*> > bitheap, int maxWeight_, bool didCompress_, int stage_);

				~Snapshot(){};

				vector<list<WeightedBit*> > bits;
				int maxWeight;
				bool didCompress; 
				int stage;
		};

		public:

			/** constructor */
			Plotter();
			
			/** destructor */
			~Plotter();

			/** takes a snapshot of the bitheap's current state */
			void heapSnapshot(bool compress, int stage);

			/** plots all the bitheap's stages */
			void plotBitHeap();

			/** plots multiplier area and lozenge views */
			void plotMultiplierConfiguration(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g);

			void setBitHeap(BitHeap* bh_);

			void addSmallMult(int topX, int topY, int dx, int dy);

			stringstream ss;

		
		private:

			/** draws an area view of the DSP configuration */
			void drawAreaView(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g);

			/** draws a lozenge view of the DSP configuration */
			void drawLozengeView(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g);

			/** draws a line between the specified coordinates */
			void drawLine(int wX, int wY, int wRez, int offsetX, int offsetY, int scalingFactor, bool isRectangle);

			/** draws a DSP block */
			void drawDSP(int wXY,  int wY, int i, int xT, int yT, int xB, int yB, int offsetX, int offsetY, int scalingFactor,  bool isRectangle);

			/** draws the target rectangle or lozenge */
			void drawTargetFigure(int wX, int wY, int offsetX, int offsetY, int scalingFactor, bool isRectangle);

			/** draws a small multiplier table */
			void drawSmallMult(int wX, int wY, int xT, int yT, int xB, int yB, int offsetX, int offsetY, int scalingFactor,  bool isRectangle);
			
			void drawLittleClock(int x, int y, int cyclenumber, int currentcycle, int stageNumber, int currentStage);
			
			string romanNumber(int i);

			/** draws a single bit */
			void drawBit(int cnt, int w, int turnaroundX, int offsetY, int c);
			
			void addECMAFunction();

			ofstream fig;
			ofstream fig2;
//			vector<vector<list<WeightedBit*> > > snapshots;
			vector<Snapshot*> snapshots;

			int topX[10000];
			int topY[10000];
			int dx, dy;
			int smallMultIndex;

			string srcFileName;

			BitHeap* bh;


	};
}

#endif
