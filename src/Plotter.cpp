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

	Plotter::Snapshot::Snapshot(vector<list<WeightedBit*> > bitheap, int minWeight_, 
			int maxWeight_, unsigned maxHeight_, bool didCompress_,  int cycle_, double cp_):
		maxWeight(maxWeight_), minWeight(minWeight_), maxHeight(maxHeight_), didCompress(didCompress_) , 
		cycle(cycle_), cp(cp_)
	{
		for(int w=minWeight; w<maxWeight_; w++)
		{

			list<WeightedBit*> t;

			if(bitheap[w].size()>0)	
			{	
				for(list<WeightedBit*>::iterator it = bitheap[w].begin(); it!=bitheap[w].end(); ++it)	
				{
					WeightedBit* b = new WeightedBit(*it);
					t.push_back(b);
				}		

				bits.push_back(t);
			}
		}
	}



	bool Plotter::Snapshot::operator< (const Snapshot& b)
	{
		if ((this->cycle<b.cycle) || ((this->cycle==b.cycle) && (this->cp<b.cp))) 
			return true;
		else
			return false;
	} 	



	Plotter::Plotter(BitHeap* bh_):bh(bh_)
	{
		srcFileName=bh_->getOp()->getSrcFileName() + ":Plotter";
		smallMultIndex = 0;
	}



	Plotter::~Plotter()
	{

	}



	void Plotter::heapSnapshot(bool compress, int cycle, double cp)
	{
		if(!(cp==0.0 && cycle==0))
			if(compress)
			{
				unsigned size=snapshots.size();	 

				Snapshot* s = new Snapshot(bh->bits, bh->getMinWeight(), bh->getMaxWeight(), bh->getMaxHeight(),
						compress, cycle, cp);
				bool proceed=true;
				//snapshots.push_back(s);
				if (size==0)
				{
					snapshots.push_back(s);

				}
				else
				{
					vector<Snapshot*>::iterator it = snapshots.begin();

					it++;	

					while(proceed) 
					{
						if (it==snapshots.end() || (*s < **it))
						{ // test in this order to avoid segfault!

							snapshots.insert(it, s);
							proceed=false;

						}
						else 
						{
							it++;
						}
					}
				}



			}
	}



	void Plotter::plotBitHeap()
	{
		drawInitialHeap();
		drawCompressionHeap();
	}



	void Plotter::setBitHeap(BitHeap* bh_)
	{
		bh = bh_;
	}



	void Plotter::plotMultiplierConfiguration(string name, vector<MultiplierBlock*> mulBlocks, 
			int wX, int wY, int wOut, int g)
	{
		drawAreaView(name, mulBlocks, wX, wY, wOut, g);
		drawLozengeView(name, mulBlocks, wX, wY, wOut, g);

	}



	void Plotter::addSmallMult(int topX_, int topY_, int dx_, int dy_)
	{

		topX[smallMultIndex]=topX_;
		topY[smallMultIndex]=topY_;
		dx=dx_;
		dy=dy_;
		smallMultIndex++;

		//ss << topX.size() << "  " ;

	}



	void Plotter::drawInitialHeap()
	{
		initializeHeapPlotting(true);

		int offsetY = 0;
		int turnaroundX = snapshots[0]->bits.size() * 10 + 80;


		offsetY += 20 + snapshots[0]->maxHeight * 10;

		drawInitialConfiguration(snapshots[0]->bits, snapshots[0]->minWeight, offsetY, turnaroundX);

		fig << "<line x1=\"" << turnaroundX + 30 << "\" y1=\"" << 20 << "\" x2=\"" << turnaroundX + 30 
			<< "\" y2=\"" << offsetY +30 << "\" style=\"stroke:midnightblue;stroke-width:1\" />" << endl;

		fig << "<text id=\"tooltip\" x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>" << endl;

		fig << "</svg>" << endl;

		fig.close();

	}



	void Plotter::drawCompressionHeap()
	{
		initializeHeapPlotting(false);

		int offsetY = 0;
		int turnaroundX = snapshots[snapshots.size()-1]->maxWeight * 10 + 100;

		//lastStage=snapshots[0]->stage;

		bool timeCondition;


		for(unsigned i=0; i< snapshots.size(); i++)
		{

			if(snapshots[i]->didCompress)
			{

				timeCondition = true;
				if (i > snapshots.size()-3)
					timeCondition = false;

				offsetY += 15 + snapshots[i]->maxHeight * 10;
				drawConfiguration(snapshots[i]->bits, i,snapshots[i]->cycle, snapshots[i]->cp, 
						snapshots[i]->minWeight, offsetY, turnaroundX, timeCondition);

				if (i!=snapshots.size()-1)
				{
					int j=i+1;

					while(snapshots[j]->didCompress==false)
						j++;

					int c = snapshots[j]->cycle - snapshots[i]->cycle;

					fig << "<line x1=\"" << turnaroundX + 150 << "\" y1=\"" 
						<< offsetY + 10 << "\" x2=\"" << 50
						<< "\" y2=\"" << offsetY + 10 << "\" style=\"stroke:lightsteelblue;stroke-width:1\" />" << endl;

					while(c>0)
					{
						offsetY += 10;
						fig << "<line x1=\"" << turnaroundX + 150 << "\" y1=\"" 
							<< offsetY  << "\" x2=\"" << 50
							<< "\" y2=\"" << offsetY  << "\" style=\"stroke:midnightblue;stroke-width:2\" />" << endl;

						c--;

					}
				}


			}

		}

		fig << "<line x1=\"" << turnaroundX + 30 << "\" y1=\"" << 20 << "\" x2=\"" << turnaroundX + 30 
			<< "\" y2=\"" << offsetY +30 << "\" style=\"stroke:midnightblue;stroke-width:1\" />" << endl;

		fig << "<text id=\"tooltip\" x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>" << endl;

		fig << "</svg>" << endl;

		fig.close();

	}



	void Plotter::drawAreaView(string name, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g)
	{
		ostringstream figureFileName;
		figureFileName << "tiling_square_" << name << ".svg";

		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);

		fig.open (figureFileName.str().c_str(), ios::trunc);


		//scaling factor for the whole drawing
		int scalingFactor = 5;

		//offsets for the X and Y axes
		int offsetX = 150;
		int offsetY = 180;

		//file header
		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" onload=\"init(evt)\" >" << endl;

		addECMAFunction();
		//draw target rectangle	



		for(int i=0; i<smallMultIndex; i++)
		{
			drawSmallMult(wX, wY,  topX[i], topY[i], topX[i] + dx, topY[i] + dy, offsetX, offsetY, scalingFactor, true);
		}

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

		fig << "<text id=\"tooltip\" x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>" << endl;

		fig << "</svg>" << endl;

		fig.close();

	}



	void Plotter::drawLozengeView(string name, vector<MultiplierBlock*> mulBlocks, int wX, int wY, int wOut, int g)
	{
		ostringstream figureFileName;
		figureFileName << "tiling_sheared_" << name << ".svg";

		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);

		fig.open (figureFileName.str().c_str(), ios::trunc);


		//scaling factor for the whole drawing
		int scalingFactor = 5;

		//offsets for the X and Y axes
		int offsetX = 180 + wY*scalingFactor;
		int offsetY = 180;

		//file header
		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" onload=\"init(evt)\">" << endl;

		addECMAFunction();

		//draw target lozenge	





		for(int i=0; i<smallMultIndex; i++)
		{
			drawSmallMult(wX, wY,  topX[i], topY[i], topX[i] + dx, topY[i] + dy, offsetX, offsetY, scalingFactor, false);
		}

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



		//drawLittleClock(500,300,3,0,5,1);


		fig << "<text id=\"tooltip\" x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>" << endl;

		fig << "</svg>" << endl;

		fig.close();

	}



	void Plotter::drawLine(int wX, int wY, int wRez, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
	{
		if(isRectangle)
			fig << "<line x1=\"" << offsetX + scalingFactor * (wRez - wY)
				<< "\" y1=\"" << offsetY << "\" x2=\"" << offsetX + scalingFactor*wRez
				<< "\" y2=\"" << offsetY + scalingFactor*wY <<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>" << endl ;	
		else
			fig << "<line x1=\"" << offsetX + scalingFactor*(wRez - wY) 
				<< "\" y1=\"" << offsetY << "\" x2=\"" << offsetX + scalingFactor*(wRez - wY)
				<< "\" y2=\"" << offsetY + scalingFactor*wY <<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>" << endl ;	
	}



	void Plotter::drawDSP(int wX, int wY, int i, int xT, int yT, int xB, int yB, 
			int offsetX, int offsetY, int scalingFactor,  bool isRectangle)
	{
		//because the X axis is opposing, all X coordinates have to be turned around
		int turnaroundX;

		int wxDSP = xB - xT;
		int wyDSP = yB - yT;

		int xTT=xT;

		if(xTT<0)
			xTT=-1;
		else
			xTT=xT/wxDSP;

		int yTT=yT;

		if(yTT<0)
			yTT=-1;
		else
			yTT=yT/wyDSP;	

		if(isRectangle)
		{
			turnaroundX = offsetX + wX * scalingFactor;
			fig << "<rect x=\"" << turnaroundX - xB*scalingFactor << "\" y=\"" << yT*scalingFactor + offsetY 
				<< "\" height=\"" << (yB-yT)*scalingFactor << "\" width=\"" << (xB-xT)*scalingFactor
				<< "\" style=\"fill:rgb(200, 200, 200);fill-opacity:0.7;stroke-width:1;stroke:rgb(0,0,0)\"" 
				<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\"X[" 
				<< xB << ":" << xT << "] * Y[" << yB << ":" << yT << "]\"/> " << endl;

			fig << "<text x=\"" << (2*turnaroundX - scalingFactor*(xT+xB))/2 -12 
				<< "\" y=\"" << ((yT+yB)*scalingFactor)/2 + offsetY + 7
				<< "\" fill=\"blue\">D(" <<  xTT+1<<")("<<  yTT+1  << ")</text>" << endl;
		}
		else 
		{
			turnaroundX = wX * scalingFactor;
			fig << "<polygon points=\"" << turnaroundX - 5*xB + offsetX - 5*yT << "," << 5*yT + offsetY << " "
				<< turnaroundX - 5*xT + offsetX - 5*yT << "," << 5*yT + offsetY << " " 
				<< turnaroundX - 5*xT + offsetX - 5*yB << "," << 5*yB + offsetY << " "
				<< turnaroundX - 5*xB + offsetX - 5*yB << "," << 5*yB + offsetY
				<< "\" style=\"fill:rgb(200, 200, 200);stroke-width:1;fill-opacity:0.7;stroke:rgb(0,0,0)\" "
				<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\"X[" 
				<< xB << ":" << xT << "] * Y[" << yB << ":" << yT << "]\"/> " << endl;

			fig << "<text x=\"" << (2*turnaroundX - xB*5 - xT*5 + 2*offsetX)/2 - 14 - (yT*5 + yB*5)/2 
				<< "\" y=\"" << ( yT*5 + offsetY + yB*5 + offsetY )/2 + 7 
				<< "\" fill=\"blue\">D(" <<  xTT+1<<")("<<  yTT+1  << ")</text>" << endl;


		}	
	}




	void Plotter::drawTargetFigure(int wX, int wY, int offsetX, int offsetY, int scalingFactor, bool isRectangle)
	{
		if(isRectangle)
			fig << "<rect x=\"" << offsetX << "\" y=\"" << offsetY 
				<< "\" height=\"" << wY * scalingFactor << "\" width=\"" << wX * scalingFactor 
				<<"\" style=\"fill:rgb(255, 0, 0);stroke-width:1;fill-opacity:0.1;stroke:rgb(0,0,0)\"/>" << endl;
		else
			fig << "<polygon points=\"" << offsetX << "," << offsetY << " " 
				<< wX*scalingFactor + offsetX << "," << offsetY << " " 
				<< wX*scalingFactor + offsetX - scalingFactor*wY << "," << wY*scalingFactor + offsetY << " "
				<< offsetX - scalingFactor*wY << "," << wY*scalingFactor + offsetY 	
				<< "\" style=\"fill:rgb(255, 0, 0);stroke-width:1;fill-opacity:0.1;stroke:rgb(0,0,0)\"/>" << endl;

	}



	void Plotter::drawSmallMult(int wX, int wY, int xT, int yT, int xB, 
			int yB, int offsetX, int offsetY, int scalingFactor,  bool isRectangle)
	{

		int turnaroundX;


		if(isRectangle)
		{
			turnaroundX = offsetX + wX * scalingFactor;
			fig << "<rect x=\"" << turnaroundX - xB*scalingFactor << "\" y=\"" << yT*scalingFactor + offsetY 
				<< "\" height=\"" << (yB-yT)*scalingFactor << "\" width=\"" << (xB-xT)*scalingFactor
				<< "\" style=\"fill:rgb(255, 228, 196);stroke-width:1;stroke:rgb(165,42,42)\"" 
				<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\"X[" 
				<< xB << ":" << xT << "] * Y[" << yB << ":" << yT << "]\"/> " << endl;

		}   
		else 
		{
			turnaroundX = wX * scalingFactor;
			fig << "<polygon points=\"" << turnaroundX - 5*xB + offsetX - 5*yT << "," << 5*yT + offsetY << " "
				<< turnaroundX - 5*xT + offsetX - 5*yT << "," << 5*yT + offsetY << " " 
				<< turnaroundX - 5*xT + offsetX - 5*yB << "," << 5*yB + offsetY << " "
				<< turnaroundX - 5*xB + offsetX - 5*yB << "," << 5*yB + offsetY
				<< "\" style=\"fill:rgb(255, 228, 196);stroke-width:1;stroke:rgb(165,42,42)\" "
				<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\"X[" 
				<< xB << ":" << xT << "] * Y[" << yB << ":" << yT << "]\"/> " << endl;


		}
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
		return "0";
	}

	void Plotter::drawLittleClock(int x, int y, int cyclenumber, int currentcycle, int stageNumber, int currentStage)
	{ 
#if 0
		ostringstream figureFileName;
		figureFileName << "littleclock"  << ".svg";

		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);

		fig2.open (figureFileName.str().c_str(), ios::trunc);
#endif
		double scalefactor=0.2;


		//file header
#if 0
		fig2 << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig2 << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig2 << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig2 << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
#endif
		fig <<"<g transform=\"scale("<<scalefactor<<") translate(1000,100)\">"<<endl;
		fig <<" <circle cx=\""<<x<<"\" cy=\""<<y
			<<"\" r=\"100\" stroke=\"black\" stroke-width=\"2\" fill=\"darkcyan\" />"<<endl;
		fig <<" <circle cx=\""<<x<<"\" cy=\""<<y
			<<"\" r=\"80\" stroke=\"black\" stroke-width=\"2\" fill=\"burlywood\" />"<<endl;

		//   fig << "<rect x=\"" << x-25<< "\" y=\"" << y-60
		//	    << "\" height=\"" << 10 <<"\" width=\"" << 60
		//	    << "\" style= \"fill:rgb(139, 69, 19);stroke:purple;stroke-width:1\" /> "<<endl;
		//fig<< "<polygon points= \" "<< x-22 <<","<< y-63 <<" "
		//						<< x-16 <<","<< y-68 <<" "
		///						<< x-10 <<","<< y-63 <<" "
		//						<< x-10 <<","<< y-45 <<" "
		//						<< x-22 <<","<< y-45 <<" \" "
		//						<< "style=\"fill:rgb(139, 134, 130);stroke-width:1;stroke:rgb(0,0,0)\"/>" << endl;

		fig << "<polygon points= \" "<< x-35 <<","<< y-30 <<" "
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
			<< x-35 <<","<< y-30 <<" \" "
			<< "style= \"fill:rgb(139, 134, 130);stroke:black;stroke-width:1\" /> "<<endl;
		fig <<"<text x=\""<<x-25<<"\" y=\""<<y-16<<"\" stroke=\"white\"  fill=\"white\" >FloPoCo</text> "<<endl;



		double adjust=3.1415+3.1415/2;
		fig <<" <circle cx=\""<<x<<"\" cy=\""<<y<<"\" r=\"4\" stroke=\"black\" stroke-width=\"2\" fill=\"black\" />"<<endl;

		for(int i=0;i<cyclenumber;i++)
		{
			double angle=360/cyclenumber*i *3.1415 / 180 + adjust;

			fig <<" <circle cx=\""<<x+80*cos(angle)<<"\" cy=\""<<y+ 80*sin(angle)
				<<"\" r=\"4\" stroke=\"black\" stroke-width=\"2\" fill=\"red\" />"<<endl;
			fig <<"<text x=\""<<x+65*cos(angle)<<"\" y=\""<<y+65*sin(angle)
				<<"\" stroke=\"red\" fill=\"red\" >"<<romanNumber(i)<<"</text> "<<endl;
			if(i==currentcycle)
			{
				fig <<"<line x1=\""<<x<<"\" y1=\""<<y<<"\" x2=\""<<x+ 70*cos(angle)
					<<"\" y2=\""<<y+ 70*sin(angle)<<"\" style=\"stroke-width: 5px; stroke:rgb(139, 69, 19);\"  />"<<endl;

				fig << "<polygon points= \" "<< x+ 72*cos(angle) <<","<< y+ 72*sin(angle) <<" "
					<< x+60*cos(angle-0.2)<<","<< y+60*sin(angle-0.2) <<" "
					<< x+60*cos(angle+0.2) <<","<< y+60*sin(angle+0.2)<<" \" "
					<< "style= \"fill:rgb(139, 134, 130);stroke:black;stroke-width:1\" /> "<<endl;
			}
		}	


		for(int i=0;i<stageNumber;i++)
		{	double adjust=3.1415+3.1415/2;
			double angle=360/stageNumber*i *3.1415 / 180 + adjust;

			fig <<" <circle cx=\""<<x+100*cos(angle)<<"\" cy=\""<<y+ 100*sin(angle)
				<<"\" r=\"4\" stroke=\"black\" stroke-width=\"2\" fill=\"yellowgreen\" />"<<endl;

			fig <<"<text x=\""<<x+117*cos(angle)<<"\" y=\""<<y+117*sin(angle)
				<<"\" stroke=\"black\" fill=\"black\" >"<<i<<"</text> "<<endl;

			if(i==currentStage)
			{	
				fig <<"<line x1=\""<<x<<"\" y1=\""<<y<<"\" x2=\""<<x+ 80*cos(angle)
					<<"\" y2=\""<<y+ 80*sin(angle)<<"\" style=\"stroke-width: 5px; stroke:rgb(139, 69, 19);\"  />"<<endl;
				fig << "<polygon points= \" "<< x+ 100*cos(angle) <<","<< y+ 100*sin(angle) <<" "
					<< x+82*cos(angle-0.2)<<","<< y+82*sin(angle-0.2) <<" "
					<< x+82*cos(angle+0.2) <<","<< y+82*sin(angle+0.2)<<" \" "
					<< "style= \"fill:rgb(139, 134, 130);stroke:black;stroke-width:1\" /> "<<endl;
			}
		}

		fig <<"<text x=\""<<x-12<<"\" y=\""<<y+70<<"\" stroke=\"midnightblue\"  fill=\"midnightblue\" >cycle</text> "<<endl;
		fig <<"<text x=\""<<x-12<<"\" y=\""<<y+95<<"\" stroke=\"ghostwhite\" fill=\"ghostwhite\" >stage</text> "<<endl;
		fig <<"</g>"<<endl;
		//fig2 << "</svg>" << endl;
		//fig2.close();

		//fig << fig2.str();

	}

	void Plotter::drawBit(int cnt, int w, int turnaroundX, int offsetY, int color, int cycle, int cp, string name)
	{
		const std::string colors[] = { "#97bf04","#0f1af2", 
			"orange", "#f5515c",  "lightgreen", "fuchsia", "indianred"};

		int index = color % 7;

		//REPORT(INFO, "bit name " << name << "   color " << color );

		int ci,c1,c2,c3;//print cp as a number as a rational number, in nanoseconds

		c3 = cp % 10;
		cp = cp / 10;
		c2 = cp % 10;	
		cp = cp / 10;
		c1 = cp % 10;
		cp = cp / 10;
		ci = cp % 10;

		fig << "<circle cx=\"" << turnaroundX - w*10 - 5 << "\" cy=\"" 
			<< offsetY - cnt*10 - 5 << "\" r=\"3\" fill=\"" << colors[index] << "\" stroke=\"black\" stroke-width=\"0.5\"" 
			<< " onmousemove=\"ShowTooltip(evt)\" onmouseout=\"HideTooltip(evt)\" mouseovertext=\""
			<< name << ", " << cycle << " : " << ci << "." << c1 << c2 << c3 << " ns\"/> " << endl;
		//<< cycle << " : " << cp << "\"/> " << endl;

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



	void Plotter::initializeHeapPlotting(bool isInitial)
	{
		ostringstream figureFileName;
		if(isInitial)
			figureFileName << "BitHeap_initial_" << bh->getName()  << ".svg";
		else 
			figureFileName << "BitHeap_compression_" << bh->getName()  << ".svg";



		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);

		fig.open (figureFileName.str().c_str(), ios::trunc);

		fig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
		fig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" onload=\"init(evt)\" >" << endl;

		addECMAFunction(); 

	}



	void Plotter::drawInitialConfiguration(vector<list<WeightedBit*> > bits, int minWeight, int offsetY, int turnaroundX)
	{
		int color = 0;
		int cnt = 0;
		//int tempCycle = 0;
		//double tempCP = 0;

		fig << "<line x1=\"" << turnaroundX + 150 << "\" y1=\"" 
			<< offsetY +10 << "\" x2=\"" << turnaroundX - bits.size()*10 - 50
			<< "\" y2=\"" << offsetY +10 << "\" style=\"stroke:lightsteelblue;stroke-width:1\" />" << endl;

		//turnaroundX -= minWeight*10;
		//

		//REPORT(DEBUG, "printing final config");	
		vector<WeightedBit*> orderedBits;




		for(unsigned i=minWeight; i<bits.size(); i++)
		{
			if(bits[i].size()>0)
			{
				//REPORT(DEBUG, "in if " << i);
				for(list<WeightedBit*>::iterator bit = bits[i].begin(); bit!=bits[i].end(); ++bit)
				{
					if(orderedBits.size()==0)
					{
						orderedBits.push_back((*bit));
						//REPORT(DEBUG, "first element");
					}
					//REPORT(DEBUG, "in sec for");
					bool proceed=true;
					vector<WeightedBit*>::iterator iterBit = orderedBits.begin();


					while(proceed)
					{
						//REPORT(DEBUG, "in while");
						if (iterBit==orderedBits.end())
						{
							orderedBits.push_back((*bit));
							proceed=false;
						}
						else
						{
							if( (**bit) == (**iterBit))
							{
								proceed=false;
							}
							else
							{
								if( (**bit) < (**iterBit))
								{
									//REPORT(DEBUG, "here *");
									orderedBits.insert(iterBit, *bit);
									proceed=false;
								}
								else
								{
									iterBit++;
								}
							}

						}
					}

#if 0
					REPORT(DEBUG, "in while");
					if( (iterBit==orderedBits.end()) || (**bit) == (**iterBit))
					{
						proceed=false;
						REPORT(DEBUG, "proceed false");
					}
					else
					{
						REPORT(DEBUG, "first else");
						if ( (**bit) < (**iterBit))
						{
							REPORT(DEBUG, "here *");
							orderedBits.insert(iterBit, *bit);
							proceed=false;
						}
						else
						{
							if ( iterBit==orderedBits.end())
							{
								REPORT(DEBUG, "here **");
								orderedBits.insert(iterBit, *bit);
								proceed=false;
							}
							else
							{
								REPORT(DEBUG, "here *** "<<orderedBits.size());
								iterBit++;
							}
						}
					}
#endif

				}
			}

#if 0
			if (it==l.end() || (*bit <= **it))
			{ // test in this order to avoid segfault!

				l.insert(it, bit);
				proceed=false;
			}
			else
			{

				it++;
			}
#endif
		}

		// for(unsigned k=0; k<orderedBits.size(); k++)
		// 	REPORT(INFO, orderedBits[k]->getCycle() << "  " << orderedBits[k]->getCriticalPath(orderedBits[k]->getCycle()));


		for(unsigned i=0; i<bits.size(); i++)		{
			if(bits[i].size()>0)			{
				//				tempCycle = 0;
				//        tempCP = 0;
				cnt = 0;
				for(list<WeightedBit*>::iterator it = bits[i].begin(); it!=bits[i].end(); ++it)
				{
					color=0;

					for (unsigned j=0; j<orderedBits.size(); j++)
					{
						if ( (**it) == (*orderedBits[j]) )
							color=j;
					}


#if 0
					if(it==bits[i].begin())
					{
						tempCycle = (*it)->getCycle();
						tempCP = (*it)->getCriticalPath(tempCycle);
					}
					else
					{
						if((tempCycle!=(*it)->getCycle()) || 
								((tempCycle==(*it)->getCycle()) && 
								 (tempCP!=(*it)->getCriticalPath((*it)->getCycle()))))
						{
							tempCycle = (*it)->getCycle();
							tempCP = (*it)->getCriticalPath(tempCycle);
							color++;
						}
					}

#endif




					int cy = (*it)->getCycle();
					double cp = (*it)->getCriticalPath(cy)*100000000000;
					//	if(stage>=(*it)->computeStage(stagesPerCycle, elemTime))
					{
						drawBit(cnt, i, turnaroundX, offsetY, color, cy, cp, (*it)->getName());
						cnt++;
					}


				}
			}
		}


	}



	void Plotter::drawConfiguration(vector<list<WeightedBit*> > bits,unsigned nr, int cycle, double criticalPath,
			int minWeight, int offsetY, int turnaroundX, bool timeCondition)
	{


		int cnt = 0;


		int ci,c1,c2,c3;//print cp as a number as a rational number, in nanoseconds

		int cpint = criticalPath * 1000000000000;

		c3 = cpint % 10;
		cpint = cpint / 10;
		c2 = cpint % 10;	
		cpint = cpint / 10;
		c1 = cpint % 10;
		cpint = cpint / 10;
		ci = cpint % 10;

		//REPORT(INFO, snapshots.size()<< " nr= "<<nr);

		if(nr ==0 )
			fig << "<text x=\"" << turnaroundX + 50 << "\" y=\"" << offsetY + 3
				<< "\" fill=\"midnightblue\">" << "before first compression" << "</text>" << endl;
		else
			if (nr == snapshots.size()-1)
				fig << "<text x=\"" << turnaroundX + 50 << "\" y=\"" << offsetY + 3
					<< "\" fill=\"midnightblue\">" << "before final addition" << "</text>" << endl;
			else
				if (nr == snapshots.size()-2)
					fig << "<text x=\"" << turnaroundX + 50 << "\" y=\"" << offsetY + 3
						<< "\" fill=\"midnightblue\">" << "before 3-bit height additions" << "</text>" << endl;
				else
				{
					fig << "<text x=\"" << turnaroundX + 50 << "\" y=\"" << offsetY + 3
						<< "\" fill=\"midnightblue\">" << cycle << "</text>" << endl;

					fig << "<text x=\"" << turnaroundX + 80 << "\" y=\"" << offsetY + 3
						<< "\" fill=\"midnightblue\">" << ci << "." << c1 << c2 << c3 << " ns"  << "</text>" << endl;
				}



		turnaroundX -= minWeight*10;

		for(unsigned i=0; i<bits.size(); i++)
		{

			//REPORT(INFO, "wtf" << i);

			if(bits[i].size()>0)
			{

				cnt = 0;
				for(list<WeightedBit*>::iterator it = bits[i].begin(); it!=bits[i].end(); ++it)
				{



					int cy = (*it)->getCycle();
					double cp = (*it)->getCriticalPath(cy);
					if (timeCondition)
					{
						//if(stage>=(*it)->computeStage(stagesPerCycle, elemTime))
						//			REPORT(INFO, cy << "  " << cycle);
						//		REPORT(INFO, cp << "  " << criticalPath);
						if ((cy<cycle) || ((cy==cycle) && (cp<=criticalPath)))
						{
							cp = cp * 1000000000000;  //picoseconds
							drawBit(cnt, i, turnaroundX, offsetY, (*it)->getType(), cy, cp, (*it)->getName());
							cnt++;
						}
					}
					else
					{
						drawBit(cnt, i, turnaroundX, offsetY, (*it)->getType(), cy, cp, (*it)->getName());
						cnt++;
					}


				}
			}
		}






	}



}

