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

#include "IntMultipliers/MultiplierBlock.hpp"


namespace flopoco
{

	class Plotter
	{
		public:

			/** constructor */
			Plotter();
			
			/** destructor */
			~Plotter();

			/** takes a snapshot of the bitheap's current state */
			void heapSnapshot(int stage);

			/** plots all the bitheap's stages */
			void plotBitHeap();

			/** plots multiplier area and lozenge views */
			void plotMultiplierConfiguration(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g);

		
		private:

			/** draws an area view of the DSP configuration */
			void drawAreaView(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g);

			/** draws a lozenge view of the DSP configuration */
			void drawLozengeView(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g);

			/** draws a line between the specified coordinates */
			void drawLine(int wX, int wY, int wRez, int offsetX, int offsetY, int scalingFactor, bool isRectangle);

			/** draws a DSP block */
			void drawDSP(int wX,  int wY, int i, int xT, int yT, int xB, int yB, int offsetX, int offsetY, int scalingFactor,  bool isRectangle);

			/** draws the target rectangle or lozenge */
			void drawTargetFigure(int wX, int wY, int offsetX, int offsetY, int scalingFactor, bool isRectangle);
			
			void drawLittleClock(int x, int y, int cyclenumber, int currentcycle, int stageNumber, int currentStage);
			
			string romanNumber(int i);

			//int wX;
			//int wY;
			//int wOut;
			//int g;
			ofstream fig;
			ofstream fig2;

	};
}

#endif
