#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "NbBitsMin.hpp"
#include "../../UserInterface.hpp"

using namespace std;

namespace flopoco
{

	void NbBitsMinRegisterFactory(){
		UserInterface::add("NbBitsMin", // name
											 "A tool for FPDiv to compute where to truncate both partial remainder and divider.",
											 "Miscellaneous", // categories
											 "",
											 "radix(int): It has to be 2^n; \
digitSet(int): the range you allow for each digit [-digitSet, digitSet]",
											 "",
											 NULL
											 ) ;

	}

	void computeNbBit (int radix, int digitSet)
	{

		float ro = (float)digitSet/(radix-1);
		cout<<"Hence your rendundancy coefficiant is "<<ro<<endl<<endl;


		double exactDeltaMin = 1-log2((2*ro - 1)/(2*(digitSet-1)));

		int deltaMinMinus = exactDeltaMin;
		int deltaMinPlus = deltaMinMinus + 1;

		double nbBitMinus = -log2((2*ro-1)/2 - (digitSet-ro)*pow(2, 0-deltaMinMinus))+deltaMinMinus+log2(radix);
		int nbBitM = ceil(nbBitMinus);

		double nbBitPlus = -log2((2*ro-1)/2 - (digitSet-ro)*pow(2, 0-deltaMinPlus))+deltaMinPlus+log2(radix);
		int nbBitP = ceil(nbBitPlus);

		cout<<"There're 2 possibilities :"<<endl;
		cout<<"-Delta = "<<deltaMinPlus<<", nbBits = "<<nbBitP<<endl;
		cout<<"-Delta = "<<deltaMinMinus<<", nbBits = "<<nbBitM<<endl<<endl;

		int delta = (nbBitP-nbBitM < 0 ? deltaMinPlus : deltaMinMinus);
		int nbBit = (nbBitP-nbBitM < 0 ? nbBitP : nbBitM);

		cout<<"Therefore, the min is for Delta = "<<delta<<endl;
		cout<<"You'll need "<<nbBit<<" bits : ";
		cout<<delta-1<<" bits for D and ";
		cout<<nbBit-delta+1;
		cout<<" bits for the partial remainder"<<endl;

		cout<<"Checking for better configurations..."<<endl<<endl;

		//Optimized parameters for a digit set [-7,7] with a radix 8
		//delta = 3;
		//nbBit = 7;
		//radix = 8;
		//digitSet = 7;

		if(checkDistrib(delta, nbBit-delta+1, radix, digitSet))
		{
			int gain = 0;
			cout<<"Optimization found!"<<endl;
			while(checkDistrib(delta, nbBit-delta, radix, digitSet))
			{
				nbBit--;
				gain++;
			}
			while(checkDistrib(delta-1, nbBit-delta+1, radix, digitSet))
			{
				delta--;
				nbBit--;
				gain++;
			}

			cout<<"You lost "<<gain<<" extra bits"<<endl;
		}
		else
		{
			cout<<"No better configuration found"<<endl;
		}


		cout<<"Final values are : ";
		cout<<delta-1<<" bits for D and ";
		cout<<nbBit-delta+1;
		cout<<" bits for the partial remainder"<<endl;

		plotPDDiagram(delta-1, nbBit-delta+1, radix, digitSet);
		//plotPDDiagram(3, 5, 8, 7);

		cout<<"An implementation of this configuration is approximately "<<estimateCost(nbBit, radix, digitSet)<<" times larger than the actual implementation"<<endl;

	}

	//Produces the P-D Diagram corresponding to the previous analysis in svg format
	void plotPDDiagram(int delta, int t, int radix, int digitSet)
	{
		float ro = (float)digitSet/((float)radix-1);

		ofstream svg("PDDiagram.svg", ios::out|ios::trunc);
		const int width = 1024;
		const int height = 750;
		int maxH = (-digitSet-ro)*(height/2)/(digitSet+1); //The height doesn't always match with the extrema of Uk and Lk,
		int minH = (digitSet+ro)*(height/2)/(digitSet+1);  //so those value are used to draw a better graph

		//Doc info
		svg<<"<?xml version=\"1.0\"?>"<<endl;
		svg<<"<!DOCTYPE  svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"<<endl;
		svg<<"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" >"<<endl;
		svg<<"<text y=\"15\">Les cases incluent leurs bords bas et gauche</text>"<<endl;
		svg<<"<text y=\"40\">		radix"<<radix<<", [-"<<digitSet<<","<<digitSet<<"]"<<"</text>"<<endl;

		//Lk-Uk
		for(int i = -digitSet; i <= digitSet; i++)
		{
			//Lk
			svg<<"	<line x1=\""<<50<<"\" y1=\""<<50+height/2<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2-(i-ro)*(height/2)/(digitSet+1);
			svg<<"\" style=\"stroke:rgb(0,200,200);stroke-width:2\"/>"<<endl;
			//Uk
			svg<<"	<line x1=\""<<50<<"\" y1=\""<<50+height/2<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2-(i+ro)*(height/2)/(digitSet+1);
			svg<<"\" style=\"stroke:rgb(0,200,200);stroke-width:2\"/>"<<endl;
			//Overlapping
			svg<<"	<line x1=\""<<50+width+(i+digitSet+1)*10<<"\" y1=\""<<50+height/2-(i-ro)*(height/2)/(digitSet+1);
			svg<<"\" x2=\""<<50+width+(i+digitSet+1)*10<<"\" y2=\""<<50+height/2-(i+ro)*(height/2)/(digitSet+1);
			svg<<"\" style=\"stroke:rgb(0,0,200);stroke-width:2\"/>"<<endl;
		}


		//Horizontal axis subdiv
		int nbSubdiv = pow(2, delta);
		int pitch = width/nbSubdiv;

		for(int i = 0; i <= nbSubdiv; i++)
		{
			svg<<"	<line x1=\""<<50+width/2+i*pitch/2<<"\" y1=\""<<50+height/2+maxH;//50+3*height/4-((delta%2==0?i:i-1)*(-digitSet-ro)*(height/2)/(digitSet+1)/2/nbSubdiv);
			svg<<"\" x2=\""<<50+width/2+i*pitch/2<<"\" y2=\""<<50+height/2+minH;// /4-((delta%2==0?i:i-1)*(digitSet+ro)*(height/2)/(digitSet+1)/2/nbSubdiv);
			svg<<"\" style=\"stroke:rgb(100,100,100);stroke-width:2\"/>"<<endl;
		}


		//Vertical axis subdiv
		nbSubdiv = pow(2, t);

		for(int i = -nbSubdiv/2; i <= nbSubdiv/2; i++)
		{
			svg<<"	<line x1=\""<<50+width/2<<"\" y1=\""<<50+height/2+i*2*minH/nbSubdiv;
			svg<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2+i*2*minH/nbSubdiv;
			svg<<"\" style=\"stroke:rgb(100,100,100);stroke-width:1\"/>"<<endl;
		}


		//Axis
		svg<<"	<line x1=\""<<50<<"\" y1=\""<<50<<"\" x2=\""<<50<<"\" y2=\""<<50+height<<"\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>"<<endl;//vertical 0
		svg<<"	<line x1=\""<<50<<"\" y1=\""<<50+height/2<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2<<"\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>"<<endl;//horizontal

		svg<<"	<line x1=\""<<50+width<<"\" y1=\""<<50+height/2+minH;
		svg<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2+maxH<<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>"<<endl;//vertical 1

		svg<<"	<line x1=\""<<50+width/2<<"\" y1=\""<<50+height/2+minH/2;
		svg<<"\" x2=\""<<50+width/2<<"\" y2=\""<<50+height/2+maxH/2<<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>"<<endl;//vertical 0.5

		//Prescaling (cas particulier base 8, digitSet 7)
	//	svg<<"	<line x1=\""<<50+width/2<<"\" y1=\""<<50+height/2+minH/2;
	//	svg<<"\" x2=\""<<50+width/2<<"\" y2=\""<<50+height/2+maxH/2<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 0.5
	//
	//	svg<<"	<line x1=\""<<50+5*width/8-3<<"\" y1=\""<<50+height/2+5*minH/8;
	//	svg<<"\" x2=\""<<50+5*width/8-3<<"\" y2=\""<<50+height/2+5*maxH/8<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 0.625-
	//
	//	svg<<"	<line x1=\""<<50+5*width/8+3<<"\" y1=\""<<50+height/2+5*minH/8;
	//	svg<<"\" x2=\""<<50+5*width/8+3<<"\" y2=\""<<50+height/2+5*maxH/8<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 0.625+
	//
	//	svg<<"	<line x1=\""<<50+3*width/4+3<<"\" y1=\""<<50+height/2+3*minH/4;
	//	svg<<"\" x2=\""<<50+3*width/4+3<<"\" y2=\""<<50+height/2+3*maxH/4<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 0.75+
	//
	//	svg<<"	<line x1=\""<<50+3*width/4-3<<"\" y1=\""<<50+height/2+3*minH/4;
	//	svg<<"\" x2=\""<<50+3*width/4-3<<"\" y2=\""<<50+height/2+3*maxH/4<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 0.75-
	//
	//	svg<<"	<line x1=\""<<50+25*width/32<<"\" y1=\""<<50+height/2+25*minH/32;
	//	svg<<"\" x2=\""<<50+25*width/32<<"\" y2=\""<<50+height/2+25*maxH/32<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 25/32
	//
	//	svg<<"	<line x1=\""<<50+15*width/16-3<<"\" y1=\""<<50+height/2+15*minH/16;
	//	svg<<"\" x2=\""<<50+15*width/16-3<<"\" y2=\""<<50+height/2+15*maxH/16<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 15/16-
	//
	//	svg<<"	<line x1=\""<<50+15*width/16+3<<"\" y1=\""<<50+height/2+15*minH/16;
	//	svg<<"\" x2=\""<<50+15*width/16+3<<"\" y2=\""<<50+height/2+15*maxH/16<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 15/16+


		//Case pb, cas particulier base 8 digitSet 7
	//	svg<<"	<rect x=\""<<50+5*width/8<<"\" y=\""<<50+3*minH/2<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+5*width/8<<"\" y=\""<<50+minH/2-minH/16<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+minH/2<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+3*minH/2-minH/16<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+3*minH/2-minH/8<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+9*minH/16+1<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;

		//Fin
		svg<<"</svg>"<<endl;
	}


	bool checkDistrib(int delta, int t, int radix, int digitSet)
	{
		float ro = digitSet/(radix-1);
		float wMax = U(digitSet, ro, 1);

		for(int k = -digitSet+1; k <= digitSet; k++)
		{

			//the position of the left intersection with Lk
			float leftL = L(k, ro, 0.5);
			float leftCeiledL = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*leftL);
			int leftLineL = (wMax-leftCeiledL)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram
			int leftCol = 0;
			bool leftCornerL = (leftL-leftCeiledL == 0);
			vector<int> crossedBoxes;

			//the position of the left intersection with Uk-1
			float leftU = U(k, ro, 0.5);
			float leftCeiledU = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*leftU);
			int leftLineU = (wMax-leftCeiledU)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram
			bool leftCornerU = (leftU-leftCeiledU == 0);

			int index;

			for(float i = 0.5; i < 1 ; i += pow(2, -delta))
			{
				//the position of the right intersection with Lk
				float rightL = L(k, ro, i+pow(2, -delta));//value of Lk for d=i
				float rightCeiledL = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*rightL);
				int rightLineL = (wMax-rightCeiledL)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram
				int rightCol = (i-0.5+pow(2, -delta))*pow(2, delta-1)/0.5;

				//the position of the right intersection with Uk-1
				float rightU = U(k-1, ro, i+pow(2, -delta));//value of Uk-1 for d=i
				float rightCeiledU = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*rightU);
				int rightLineU = (wMax-rightCeiledU)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram

				bool rightCornerL = (rightL-rightCeiledL == 0);
				bool rightCornerU = (rightU-rightCeiledU == 0);

				if(leftCornerL && k-ro>0)
				{						//positive coeff =>the bottom-left corner, the box crossed is above the computed-one
										//       _______________
										//	    | /				|   <====== crossed box
					leftLineL--;		//		|/______________|
										//		/               |   <====== computed crossed box
										//	   /|_______________|
				}
				if(rightCornerL && k-ro<0)
				{
					rightLineL--;
				}


				if(leftCornerU && k-ro>0)
				{
					leftLineU--;
				}
				if(rightCornerU && k-ro<0)
				{
					rightLineU--;
				}




				for(int c = (k-ro<0 ? leftLineL : rightLineL); c <= (k-ro<0 ? rightLineL : leftLineL); c++)
				{
					index = leftCol + pow(2,delta-1)*c;
					crossedBoxes.push_back(index);
				}

				for(int c = (k-ro<0 ? leftLineU : rightLineU); c <= (k-ro<0 ? rightLineU : leftLineU); c++)
				{
					index = leftCol + pow(2,delta-1)*c;
					for(vector<int>::iterator it = crossedBoxes.begin(); it != crossedBoxes.end(); ++it)
					{
						if (index == *it)
						{
							return false;
						}
					}
				}



				rightLineL+=(rightCornerL && k-ro<0 ? 1 : 0);
				rightLineU+=(rightCornerU && k-ro<0 ? 1 : 0);

				leftL = rightL;
				leftCeiledL = rightCeiledL;
				leftLineL = rightLineL;
				leftU = rightU;
				leftCeiledU = rightCeiledU;
				leftLineU = rightLineU;
				leftCol = rightCol;

				crossedBoxes.clear();
			}

		}

		return true;
	}

	float L(int k, float ro, float d)
	{
		return (k-ro)*d;
	}

	float U(int k, float ro, float d)
	{
		return (k+ro)*d;
	}

	float estimateCost(int nbBit, int radix, int digitSet)
	{
		int nbIter = (23+6)>>1;
		float flopocoCost = nbIter*(2*27 + 3);
		cout<<"FlopocoCost = "<<flopocoCost<<endl;
		float res = (nbIter*2/log2(radix))*((ceil(log2(digitSet))+1)*ceil(pow(2, nbBit-6))+27*(1+ceil(pow(2, log2(digitSet)-7))));
		cout<<"Cost of this configuration = "<<res<<endl;
		return res/flopocoCost;
	}
}
