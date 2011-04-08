/*
  TaMaDiShiftRegister evaluation of worst cases for rounding functions without 
  overheating the planet

  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author :   Bogdan Pasca

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

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

#include "TaMaDiShiftRegister.hpp"


using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	TaMaDiShiftRegister::TaMaDiShiftRegister(Target* target, int widthLocation, int n):
	Operator(target)
	{
		srcFileName="TaMaDiShiftRegister";
		ostringstream name;

		name <<"TaMaDiShiftRegister_width"<<widthLocation<<"_dimmension"<<n<<"_uid"<<getNewUId();
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca (2011)");		

		addInput   ("DataIn", widthLocation*n); //parallel load
		addInput   ("ParallelWrite");
		
		addOutput  ("DataOut",widthLocation);
		addOutput  ("Ready");
		
		int counterWidth = intlog2(n-1);

		declare("counter",counterWidth);
		declare("readys");
		for (int i=n-1;i>=0;i--)
			declare( join("reg",i), widthLocation );
		
		vhdl << tab << "process(clk, rst, ParallelWrite)"<<endl;
		vhdl << tab << "	begin"<<endl;
		vhdl << tab << " 		if rst='1' then "<<endl;
		vhdl << tab << "		   readys <= '1';"<<endl;
		vhdl << tab << "		   counter <= (others => '0');"<<endl;
		vhdl << tab << "		elsif clk'event and clk='1' then"<<endl;
		vhdl << tab << "			if (ParallelWrite ='1') then"<<endl;
		vhdl << tab << "				counter <= CONV_STD_LOGIC_VECTOR("<<n-1<<","<<counterWidth<<")"<<";"<<endl;
		for (int i=0;i<n;i++) //parallel Load
			vhdl << tab << join("reg",i) << "<= DataIn"<<range((i+1)*widthLocation-1,i*widthLocation)<<";"<<endl;
		vhdl << tab << "				readys <= '0';"<<endl;
		vhdl << tab << "			else"<<endl;
		for (int i=0;i<n-1;i++) //parallel Load
			vhdl << tab << join("reg",i) << "<= "<<join("reg",i+1)<<";"<<endl;

		vhdl << tab << join("reg",n-1) << " <= " << zg(widthLocation,0)<<";"<<endl;
		vhdl << tab << " if counter > 1 then"<<endl;
		vhdl << tab << "	counter <= counter - '1';"<<endl;
		vhdl << tab << " 	readys <= '0';"<<endl;
		vhdl << tab << " else "<<endl;
		vhdl << tab << " 	readys <= '1';"<<endl;
		vhdl << tab << " end if;"<<endl; 
		vhdl << tab << "			end if;"<<endl;
		vhdl << tab << " end if;"<<endl; 
		vhdl << tab << "end process;"<<endl;

		vhdl << tab << "DataOut <= reg0;"<<endl;
		vhdl << tab << "Ready <= readys;"<<endl;
	}

	TaMaDiShiftRegister::~TaMaDiShiftRegister() {
	}

}



