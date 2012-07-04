#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "BasicCompressor.hpp"


using namespace std;


// personalized parameter
//string BasicCompressor::operatorInfo = "UserDefinedInfo list parameter;



BasicCompressor::BasicCompressor(Target * target, vector<int> h)
:Operator(target)
{
	ostringstream name;
	stringstream nm;
	
	int w=0;
	int param=0;
	
	vector<int> height(h.size(), 0);
	
	for(unsigned i=0; i<h.size();i++)
		height[h.size()-i-1]=h[i];
	
	name << "Compressor_";
	
	for(unsigned i=0; i<height.size();i++)
	{
		w=w+height[i];
		param=param+intpow2(height.size()-i-1)*height[i];
		name<<height[i];
	}

	int wOut=intlog2(param);
	
	name << "_" << wOut;
	setName(name.str());
	setCopyrightString("Bogdan Popa & Kinga Illyes 2012");
	
	
	
	stringstream xs;

	
	for(unsigned i=0;i<height.size();i++)
	{
		addInput(join("X",i), h[i]);
		
		if(i!=0)
		{
			xs<<"& X"<<height.size()-i-1<<" ";
		}
		else
		{
			xs<<"X"<<height.size()-1<<" ";	
		}
	}	
	
	xs<<";\n";
	
	addOutput("R", wOut);
	
	vhdl << tab << declare("X", w) << " <=" << xs.str();
	
	vhdl << "with X select R <= \n";
	
	for (mpz_class i = 0; i < (1 << w); i++) 
	{
		
		mpz_class ppcnt=0;//popcnt(i);
		mpz_class ii=i;
		for(unsigned j=0;j<h.size();j++)
		{
			
		
			ppcnt+=popcnt(ii-((ii>>h[j])<<h[j]))*intpow2(j);
			ii=ii>>h[j];
		}
		
		vhdl << tab << "\"" << unsignedBinary(ppcnt,wOut) << "\" when \""
		<< unsignedBinary(i,w) << "\", \n";
		
		
		 }
		 
		 
		 
		 vhdl << tab << "\"" << std::string (wOut, '-') << "\" when others;\n" << endl;
		 
		 
	};
	
	
	void BasicCompressor::emulate(TestCase * tc)
	{
		mpz_class sx = tc->getInputValue("X");
		tc->addExpectedOutput("R", popcnt(sx));
	}
	
	
	