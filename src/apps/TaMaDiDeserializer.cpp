/*
  TaMaDiDeserializer evaluation of worst cases for rounding functions without 
  overheating the planet

  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author :   Bogdan Pasca

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
  All rights reserved.

 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <cstdlib>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"

#include "TaMaDiDeserializer.hpp"
#include "TaMaDiModule.hpp"


using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	TaMaDiDeserializer::TaMaDiDeserializer(Target* target, int wp, int d, int iterations, int wIntervalID, int n, int inFIFODepth, int peFIFODepth, int outFIFODepth):
	Operator(target), wp(wp), d(d), iterations(iterations), wIntervalID(wIntervalID), n(n) 
	{
		srcFileName="TaMaDiDeserializer";
		ostringstream name;

		name <<"TaMaDiDeserializer_wp"<<wp<<"_interations"<<iterations<<"_degree"<<d<<"_uid"<<getNewUId();
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca (2011)");		
		
		addInput   ("MainInput", 64);
		addInput   ("MainInputValid");
		addOutput  ("MainFIFOOutput", (d+1)*wp + wIntervalID);
		addOutput  ("ValidData");
		//FIXME DO proper Deserializer implementation
		
		vhdl << tab << declare("mainSignal",(d+1)*wp+wIntervalID) << " <= ";
		int lastC = ((d+1)*wp+wIntervalID) % 64;
		vhdl << "MainInput"<<range(lastC-1,0);
		for (int i=0; i< ((d+1)*wp+wIntervalID-lastC)/64; i++)
			vhdl << "& MainInput";
		vhdl << ";"<<endl;
		setCycle(1);
		vhdl << tab << "MainFIFOOutput <= mainSignal;" <<endl;
		vhdl << tab << "ValidData <= '1';"<<endl;
	}

	TaMaDiDeserializer::~TaMaDiDeserializer() {
	}

}



