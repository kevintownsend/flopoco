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

#include "Plotter.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	
#include <vector>
#include <list>

using namespace std;

namespace flopoco
{
	Plotter::Plotter()
	{
		
	}



	Plotter::~Plotter()
	{

	}



	void Plotter::heapSnapshot(int stage)
	{

	}



	void Plotter::plotBitHeap()
	{
	
	}



	void Plotter::plotMultiplierConfiguration(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g)
	{
		drawAreaView(uid, mulBlocks, wX, wY, wOut, g);
		drawLozengeView(uid, mulBlocks, wX, wY, wOut, g);
	}



	void Plotter::drawAreaView(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g)
	{
		ostringstream figureFileName;
		figureFileName << "view_area_" << uid << ".svg";
		
		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fig.open (figureFileName.str().c_str(), ios::trunc);


		//scaling factor for the whole drawing
		int scalingFactor = 5;

		//offsets for the X and Y axes
		int offsetX = 200;
		int offsetY = 120;

		//file header
		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;

		//draw target rectangle	
		drawTargetFigure(wX, wY, offsetX, offsetY, scalingFactor, true);

		//draw DSPs
		int xT, yT, xB, yB;

		for(unsigned i=0; i<mulBlocks.size(); i++)
		{
			xT = mulBlocks[i]->gettopX();
			yT = mulBlocks[i]->gettopY();
			xB = mulBlocks[i]->getbotX();
			yB = mulBlocks[i]->getbotY();
			drawDSP(wX, wY, i, xT, yT, xB, yB, offsetX, offsetY, scalingFactor, true);
		}



		//draw truncation line
		if(wX+wY-wOut > 0)
		{
			drawLine(wX, wY, wOut, offsetX, offsetY, scalingFactor, true);    
		}

		//draw guard line
		if(g>0)
		{
			drawLine(wX, wY, wOut+g, offsetX, offsetY, scalingFactor, true);
		}

		fig << "</svg>" << endl;

		fig.close();

	}



	void Plotter::drawLozengeView(int uid, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g)
	{
		ostringstream figureFileName;
		figureFileName << "view_lozenge_" << uid << ".svg";
		
		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fig.open (figureFileName.str().c_str(), ios::trunc);


		//scaling factor for the whole drawing
		int scalingFactor = 5;

		//offsets for the X and Y axes
		int offsetX = 180 + wY*scalingFactor;
		int offsetY = 120;

		//file header
		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;

		//draw target lozenge	
		drawTargetFigure(wX, wY, offsetX, offsetY, scalingFactor, false);

		//draw DSPs
		int xT, yT, xB, yB;

		for(unsigned i=0; i<mulBlocks.size(); i++)
		{
			xT = mulBlocks[i]->gettopX();
			yT = mulBlocks[i]->gettopY();
			xB = mulBlocks[i]->getbotX();
			yB = mulBlocks[i]->getbotY();
			drawDSP(wX, wY, i, xT, yT, xB, yB, offsetX, offsetY, scalingFactor, false);
		}


		//draw truncation line
		if(wX+wY-wOut > 0)
		{
			drawLine(wX, wY, wOut, offsetX, offsetY, scalingFactor, false);    
		}

		//draw guard line
		if(g>0)
		{
			drawLine(wX, wY, wOut+g, offsetX, offsetY, scalingFactor, false);
		}


		fig << "</svg>" << endl;

		fig.close();

	}



	void Plotter::drawLine(int wX, int wY, int wRez, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
	{
		if(isRectangle)
			fig << "<line x1=\"" << offsetX + scalingFactor*(wRez - wX) << "\" y1=\"" << offsetY << "\" x2=\"" << offsetX + scalingFactor*wRez
			    << "\" y2=\"" << offsetY + scalingFactor*wY <<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>" << endl ;	
		else
			fig << "<line x1=\"" << offsetX + scalingFactor*(wRez - wY) << "\" y1=\"" << offsetY << "\" x2=\"" << offsetX + scalingFactor*(wRez - wY)
			    << "\" y2=\"" << offsetY + scalingFactor*wY <<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>" << endl ;	
	}



	void Plotter::drawDSP(int wX, int wY, int i, int xT, int yT, int xB, int yB, int offsetX, int offsetY, int scalingFactor,  bool isRectangle)
	{
		//because the X axis is opposing, all X coordinates have to be turned around
		int turnaroundX;

		int wxDSP = xB - xT;
		int wyDSP = yB - yT;

		if(isRectangle)
		{
			turnaroundX = offsetX + wX * scalingFactor;
			fig << "<rect x=\"" << turnaroundX - xB*scalingFactor << "\" y=\"" << yT*scalingFactor + offsetY 
			    << "\" height=\"" << (yB-yT)*scalingFactor << "\" width=\"" << (xB-xT)*scalingFactor
			    << "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;
			fig << "<text x=\"" << (2*turnaroundX - scalingFactor*(xT+xB))/2 -12 << "\" y=\"" << ((yT+yB)*scalingFactor)/2 + offsetY + 7
			    << "\" fill=\"blue\">D" <<  xT / wxDSP <<  yT / wyDSP  << "</text>" << endl;
		}   
		else 
		{
			turnaroundX = wX * scalingFactor;
			fig << "<polygon points=\"" << turnaroundX - 5*xB + offsetX - 5*yT << "," << 5*yT + offsetY << " "
			    << turnaroundX - 5*xT + offsetX - 5*yT << "," << 5*yT + offsetY << " " 
			    << turnaroundX - 5*xT + offsetX - 5*yB << "," << 5*yB + offsetY << " "
			    << turnaroundX - 5*xB + offsetX - 5*yB << "," << 5*yB + offsetY
			    << "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;

			fig << "<text x=\"" << (2*turnaroundX - xB*5 - xT*5 + 2*offsetX)/2 - 14 - (yT*5 + yB*5)/2 
			    << "\" y=\"" << ( yT*5 + offsetY + yB*5 + offsetY )/2 + 7 
			    << "\" fill=\"blue\">D" <<  xT / wxDSP <<  yT / wyDSP  << "</text>" << endl;	
	
		}	
	}



	void Plotter::drawTargetFigure(int wX, int wY, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
	{
		if(isRectangle)
			fig << "<rect x=\"" << offsetX << "\" y=\"" << offsetY << "\" height=\"" << wY * scalingFactor << "\" width=\"" << wX * scalingFactor 
			    <<"\" style=\"fill:rgb(255, 255, 255);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;
		else
			fig << "<polygon points=\"" << offsetX << "," << offsetY << " " 
			    << wX*scalingFactor + offsetX << "," << offsetY << " " 
			    << wX*scalingFactor + offsetX - scalingFactor*wY << "," << wY*scalingFactor + offsetY << " "
			    << offsetX - scalingFactor*wY << "," << wY*scalingFactor + offsetY 	
			    << "\" style=\"fill:rgb(255, 255, 255);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;

	}



}

