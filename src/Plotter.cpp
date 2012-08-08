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
#include "BitHeap.hpp"
#include "Plotter.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	
#include <vector>
#include <list>
#include <math.h>

using namespace std;

namespace flopoco
{

	Plotter::Plotter()
	{
		//srcFileName=bh->getOp()->getSrcFileName() + ":Plotter";

	}



	Plotter::~Plotter()
	{

	}


	void Plotter::heapSnapshot(bool compress, int stage)
//	void Plotter::heapSnapshot(vector<list<WeightedBit*> > heap, bool compress, int stage)
	{
		snapshots.push_back(bh);
		didCompress.push_back(compress);
		stages.push_back(stage);

		//srcFileName=bh->getOp()->getSrcFileName() + ":Plotter";
		

		//REPORT(DEBUG, "sdfsdfds");
	}



	void Plotter::plotBitHeap()
	{
	
	}



	void Plotter::setBitHeap(BitHeap* bh_)
	{
		bh = bh_;
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
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" onload=\"init(evt)\" >" << endl;

		addECMAFunction();
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
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" onload=\"init(evt)\" >" << endl;
		
		addECMAFunction();

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
      
		drawLittleClock(500,300,3,0,5,1);


		fig << "<text id=\"tooltip\" x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>" << endl;

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
			    << "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\"" 
				<< "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\" "
				<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\"X[" 
				<< xB << ":" << xT << "] * Y[" << yB << ":" << yT << "]\"/> " << endl;

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
			    << "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;stroke:rgb(0,0,0)\" "
				<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\"X[" 
				<< xB << ":" << xT << "] * Y[" << yB << ":" << yT << "]\"/> " << endl;

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
	
	string Plotter::romanNumber(int i)
	{
	if(i==0) return "0";
	if(i==1) return  "I";
	if(i==2) return "II";
	if(i==3) return "III";
	if(i==4) return"IV";
	if(i==5) return "V";
	if(i==6) return "VI";
	if(i==7) return"VII";
	if(i==8) return "VIII";
	if(i==9) return "IX";
	if(i==10) return"X";
	}
	
	void Plotter::drawLittleClock(int x, int y, int cyclenumber, int currentcycle, int stageNumber, int currentStage)
	{ //TODO correct fig2
	   ostringstream figureFileName;
		figureFileName << "littleclock"  << ".svg";
		
		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fig2.open (figureFileName.str().c_str(), ios::trunc);

		double scalefactor=3.0;
		

		//file header
		fig2 << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig2 << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig2 << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig2 << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
      	fig2<<"<g transform=\"scale("<<scalefactor<<") translate(1000,100)\">"<<endl;
		fig2<<" <circle cx=\""<<x<<"\" cy=\""<<y<<"\" r=\"100\" stroke=\"black\" stroke-width=\"2\" fill=\"darkcyan\" />"<<endl;
        fig2<<" <circle cx=\""<<x<<"\" cy=\""<<y<<"\" r=\"80\" stroke=\"black\" stroke-width=\"2\" fill=\"burlywood\" />"<<endl;
        
      //   fig << "<rect x=\"" << x-25<< "\" y=\"" << y-60
		//	    << "\" height=\"" << 10 <<"\" width=\"" << 60
		//	    << "\" style= \"fill:rgb(139, 69, 19);stroke:purple;stroke-width:1\" /> "<<endl;
		//fig<< "<polygon points= \" "<< x-22 <<","<< y-63 <<" "
		//						<< x-16 <<","<< y-68 <<" "
		///						<< x-10 <<","<< y-63 <<" "
		//						<< x-10 <<","<< y-45 <<" "
		//						<< x-22 <<","<< y-45 <<" \" "<< "style=\"fill:rgb(139, 134, 130);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;
        
		fig2<< "<polygon points= \" "<< x-35 <<","<< y-30 <<" "
								<< x+35 <<","<< y-30 <<" "
								<< x+50 <<","<< y-23 <<" "
								<< x+35 <<","<< y-15 <<" "
								<< x+25 <<","<< y-15 <<" "
								<< x+15 <<","<< y- 5 <<" "
								<< x+15 <<","<< y +5 <<" "
								<< x+25 <<","<< y+15 <<" "
								<< x+35 <<","<< y+15 <<" "
								<< x+35 <<","<< y+30 <<" "
								<< x-35 <<","<< y+30 <<" "
								<< x-35 <<","<< y+15 <<" "
								<< x-25 <<","<< y+15 <<" "
								<< x-15 <<","<< y+5 <<" "
								<< x-15 <<","<< y-5 <<" "
								<< x-25 <<","<< y-15 <<" "
								<< x-35 <<","<< y-15 <<" "
								<< x-35 <<","<< y-30 <<" \" "<< "style= \"fill:rgb(139, 134, 130);stroke:black;stroke-width:1\" /> "<<endl;
		fig2<<"<text x=\""<<x-25<<"\" y=\""<<y-16<<"\" stroke=\"white\"  fill=\"white\" >FloPoCo</text> "<<endl;
		
		
	
		double adjust=3.1415+3.1415/2;
		fig2<<" <circle cx=\""<<x<<"\" cy=\""<<y<<"\" r=\"4\" stroke=\"black\" stroke-width=\"2\" fill=\"black\" />"<<endl;
		
		for(int i=0;i<cyclenumber;i++)
		{
			double angle=360/cyclenumber*i *3.1415 / 180 + adjust;
			
			fig2<<" <circle cx=\""<<x+80*cos(angle)<<"\" cy=\""<<y+ 80*sin(angle)<<"\" r=\"4\" stroke=\"black\" stroke-width=\"2\" fill=\"red\" />"<<endl;
			fig2<<"<text x=\""<<x+65*cos(angle)<<"\" y=\""<<y+65*sin(angle)<<"\" stroke=\"red\" fill=\"red\" >"<<romanNumber(i)<<"</text> "<<endl;
			if(i==currentcycle)
			{
				fig2<<"<line x1=\""<<x<<"\" y1=\""<<y<<"\" x2=\""<<x+ 70*cos(angle)
					<<"\" y2=\""<<y+ 70*sin(angle)<<"\" style=\"stroke-width: 5px; stroke:rgb(139, 69, 19);\"  />"<<endl;
				
				fig2<< "<polygon points= \" "<< x+ 72*cos(angle) <<","<< y+ 72*sin(angle) <<" "
								<< x+60*cos(angle-0.2)<<","<< y+60*sin(angle-0.2) <<" "
								<< x+60*cos(angle+0.2) <<","<< y+60*sin(angle+0.2)<<" \" "
								<< "style= \"fill:rgb(139, 134, 130);stroke:black;stroke-width:1\" /> "<<endl;
			}
		}	
		
		
		for(int i=0;i<stageNumber;i++)
		{	double adjust=3.1415+3.1415/2;
			double angle=360/stageNumber*i *3.1415 / 180 + adjust;
			
			fig2<<" <circle cx=\""<<x+100*cos(angle)<<"\" cy=\""<<y+ 100*sin(angle)
				<<"\" r=\"4\" stroke=\"black\" stroke-width=\"2\" fill=\"yellowgreen\" />"<<endl;
			
			fig2<<"<text x=\""<<x+117*cos(angle)<<"\" y=\""<<y+117*sin(angle)
				<<"\" stroke=\"black\" fill=\"black\" >"<<i<<"</text> "<<endl;
		
			if(i==currentStage)
			{	
				fig2<<"<line x1=\""<<x<<"\" y1=\""<<y<<"\" x2=\""<<x+ 80*cos(angle)
					<<"\" y2=\""<<y+ 80*sin(angle)<<"\" style=\"stroke-width: 5px; stroke:rgb(139, 69, 19);\"  />"<<endl;
				fig2<< "<polygon points= \" "<< x+ 100*cos(angle) <<","<< y+ 100*sin(angle) <<" "
					<< x+82*cos(angle-0.2)<<","<< y+82*sin(angle-0.2) <<" "
					<< x+82*cos(angle+0.2) <<","<< y+82*sin(angle+0.2)<<" \" "
					<< "style= \"fill:rgb(139, 134, 130);stroke:black;stroke-width:1\" /> "<<endl;
			}
		}
		
		fig2<<"<text x=\""<<x-12<<"\" y=\""<<y+70<<"\" stroke=\"midnightblue\"  fill=\"midnightblue\" >cycle</text> "<<endl;
		fig2<<"<text x=\""<<x-12<<"\" y=\""<<y+95<<"\" stroke=\"ghostwhite\" fill=\"ghostwhite\" >stage</text> "<<endl;
		fig2<<"</g>"<<endl;
		fig2 << "</svg>" << endl;
		fig2.close();
		
	}

	void Plotter::drawBit(int cnt, int w, int turnaroundX, int offsetY, int c)
	{
		const std::string colors[] = { "#97bf04","#0f1af2", "#f5515c", "#3958ff","#f2eb8d", "indianred", "yellow", "lightgreen"};

		int index = c % 8;

		fig << "<circle cx=\"" << turnaroundX - w*10 - 5 << "\" cy=\"" << offsetY - cnt*10 - 5 << "\" r=\"3\" fill=\"" << colors[index] << "\"/>" << endl;
	}

	void Plotter::addECMAFunction()
	{
		fig << "<script> <![CDATA[  " << endl;
		fig << "function init(evt) {" << endl;
		fig << "if ( window.svgDocument == null ) {" << endl;
		fig << "svgDocument = evt.target.ownerDocument;}" << endl;
		fig << "tooltip = svgDocument.getElementById('tooltip');}" << endl;
		fig << "function ShowTooltip(evt) {" << endl;
		fig << "tooltip.setAttributeNS(null,\"x\",evt.clientX+10);" << endl;
		fig << "tooltip.setAttributeNS(null,\"y\",evt.clientY+30);" << endl;
		fig << "tooltip.firstChild.data = evt.target.getAttributeNS(null,\"mouseovertext\");" << endl;
		fig << "tooltip.setAttributeNS(null,\"visibility\",\"visible\");}" << endl;
		fig << "function HideTooltip(evt) {" << endl;
		fig << "tooltip.setAttributeNS(null,\"visibility\",\"hidden\");}]]></script>" << endl;

	}


#if 0
	void Plotter::initializeBitHeapPlotting()
	{
		ostringstream figureFileName;
		figureFileName << "bit_heap.svg";
	
		drawCycleNumber = 0;

		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fileFig.open (figureFileName.str().c_str(), ios::trunc);

		fileFig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fileFig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fileFig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
	}




	void BitHeap::drawConfiguration(int offsetY)
	{
		int turnaroundX = 1500;
		int color = 0;
		int tempCycle = 0;
		int cnt = 0;
		double tempCP = 0;
		if((drawCycleLine) && (drawCycleNumber==0))
			{
				fig << "<text x=\"" << turnaroundX + 85 << "\" y=\"" << 40
				    << "\" fill=\"midnightblue\">" << "Cycle"<< "</text>" << endl;
			}

		if(drawCycleLine)
			{
				drawCycleNumber++;
				fig << "<line x1=\"" << turnaroundX + 200 << "\" y1=\"" << offsetY +10 << "\" x2=\"" << turnaroundX - bits.size()*10 - 50
				    << "\" y2=\"" << offsetY +10 << "\" style=\"stroke:midnightblue;stroke-width:2\" />" << endl;

				fig << "<text x=\"" << turnaroundX + 100 << "\" y=\"" << offsetY + 3
				    << "\" fill=\"midnightblue\">" << drawCycleNumber - 1 << "</text>" << endl;

				fig << "<text x=\"" << turnaroundX + 100 << "\" y=\"" << offsetY + 27
				    << "\" fill=\"midnightblue\">" << drawCycleNumber  << "</text>" << endl;

				drawCycleLine = false;
			}
		else
			{
				fig << "<line x1=\"" << turnaroundX + 200 << "\" y1=\"" << offsetY +10 << "\" x2=\"" << turnaroundX - bits.size()*10 - 50
				    << "\" y2=\"" << offsetY +10 << "\" style=\"stroke:lightsteelblue;stroke-width:1\" />" << endl;
				drawCycleLine = false;
			}

		for(unsigned i=0; i<bits.size(); i++)
			{

				if(bits[i].size()>0)
					{
						color=0;
						tempCycle = 0;
						tempCP = 0;
						cnt = 0;
						for(list<WeightedBit*>::iterator it = bits[i].begin(); it!=bits[i].end(); ++it)
							{
								if(it==bits[i].begin())
									{
										tempCycle = (*it)->getCycle();
										tempCP = (*it)->getCriticalPath(tempCycle);
									}
								else
									{
										if((tempCycle!=(*it)->getCycle()) || 
										   ((tempCycle==(*it)->getCycle()) && (tempCP!=(*it)->getCriticalPath((*it)->getCycle()))))
											{
												tempCycle = (*it)->getCycle();
												tempCP = (*it)->getCriticalPath(tempCycle);
												color++;
											}
									}

								drawBit(cnt, i, turnaroundX, offsetY, (*it)->computeStage(stagesPerCycle, elementaryTime));
								cnt++;
							}
					}
			}

        
	}


	
	void Plotter::closeBitHeapPlotting(int offsetY)
	{
		int turnaroundX = 1500;
		fig << "<line x1=\"" << turnaroundX + 50 << "\" y1=\"" << 20 << "\" x2=\"" << turnaroundX + 50
		    << "\" y2=\"" << offsetY +30 << "\" style=\"stroke:midnightblue;stroke-width:1\" />" << endl;

		fig << "</g></svg>" << endl;

		fileFig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"100%\" height=\"100%\" viewBox=\"" <<turnaroundX - bits.size()*10 - 80 
		        << " " << 100 << " " << turnaroundX + 50 << " " << offsetY + 20 <<  "\">" << endl; 
		fileFig << "<g transform=\"rotate(45)\">" <<endl;
		//    fileFig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"0 0 500 500\">" << endl; 

		fileFig << fig.str();
       
               

		fileFig.close();
	}


#endif


}

