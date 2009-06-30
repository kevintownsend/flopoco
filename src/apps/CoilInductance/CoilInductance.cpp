/*
 * Floating Point Adder for FloPoCo
 *
 * Author :  Radu Tudoran, Bogdan Pasca
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>


#include <gmpxx.h>
#include "../../utils.hpp"
#include "../../Operator.hpp"

#include "CoilInductance.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <math.h>
#include <locale>

#include <stdio.h>
#include <mpfr.h>

using namespace std;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


CoilInductance::CoilInductance(Target* target, int LSBI, int MSBI, int LSBO, int MSBO, char *filepath) :
	Operator(target), MSBI(MSBI), LSBI(LSBI), LSBO(LSBO), MSBO(MSBO) ,filepath(filepath){
	
	if ((MSBI < LSBI)){
		cerr << 
			" CoilInductance: Input constraint LSBI <= MSBI not met."<<endl;
		exit (EXIT_FAILURE);
	}
	
	if ((MSBO < LSBO)){
		cerr << 
			" CoilInductance: Input constraint LSBO <= MSBO not met."<<endl;
		exit (EXIT_FAILURE);
	}
	
	wE=8;
	wF=23;
	
	ostringstream name; 
	name <<"CoilInductance_"<<abs(LSBI)<<"_"<<abs(MSBI)<<"_"<<abs(LSBO)<<"_"<<abs(MSBO);
	setName(name.str());
	
	setCopyrightString("Bogdan Pasca, Radu Tudoran (2009)");
	
	inputWidth= MSBI-LSBI;
	outputWidth= MSBO-LSBO;
	integratorWidth=3;
	
	addOutput("O",outputWidth);
	
	//Counters for addressing the memories and for frequency division
	
	vhdl<<endl;
	vhdl<<tab<<"process (clk)"<<endl<<tab<<"variable count:std_logic_vector"<<range(integratorWidth,0)<<":=(others=>'0');"<<endl<<tab<<"begin"<<endl<<tab<<tab<<"if clk'event and clk = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if count < "<<pow(2,integratorWidth)<<" then "<<endl<<tab<<tab<<tab<<tab<<"count:=count+'1';"<<endl<<tab<<tab<<tab<<tab<<declare("out_clk1",1)<<"<= '0';"<<endl<<tab<<tab<<tab<<tab<<declare("out_rst",1)<<"<= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<"elsif count = "<<pow(2,integratorWidth)<<" then "<<endl<<tab<<tab<<tab<<tab<<"count:=count+'1';"<<endl<<tab<<tab<<tab<<tab<<use("out_clk1")<<"<= '1'; "<<endl<<tab<<tab<<tab<<tab<<use("out_rst")<<"<= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<"else "<<endl<<tab<<tab<<tab<<tab<<"count:= CONV_STD_LOGIC_VECTOR(0,"<<integratorWidth<<");"<<endl<<tab<<tab<<tab<<tab<<use("out_clk1")<<"<= '0'; "<<endl<<tab<<tab<<tab<<tab<<use("out_rst")<<"<= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<declare("signal_t",integratorWidth)<<"<="<<" count"<<range(integratorWidth-1,0)<<";"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	//Memories instantiation
	
	//Computing the segments x1-x0 x3-x2 y1-y0 y3-y2 z1-z0 z3-z2
	
	//syncCycleFromSignal("????"); sincronization with memories
	
	
	
	vhdl<<tab<<declare("signal_x0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_x1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_x2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_x3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	
	vhdl<<tab<<declare("signal_y0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_y1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_y2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_y3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	vhdl<<tab<<declare("signal_z0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_z1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_z2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_z3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	
		
	
	//performing x1-x0
	
	setCycleFromSignal("signal_x0");
	
	vhdl<<tab<<declare("segmentX1mX0",inputWidth)<<"<="<<use("signal_x1")<<" - "<<use("signal_x0")<<";"<<endl;
	
	//performing x0-x1
	
	setCycleFromSignal("signal_x1");
	
	vhdl<<tab<<declare("segmentX0mX1",inputWidth)<<"<="<<use("signal_x0")<<" - "<<use("signal_x1")<<";"<<endl;
	
	//performing x0-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("segmentX0mX2",inputWidth)<<"<="<<use("signal_x0")<<" - "<<use("signal_x2")<<";"<<endl;
	
	
	//performing x3-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("segmentX3mX2",inputWidth)<<"<="<<use("signal_x3")<<" - "<<use("signal_x2")<<";"<<endl;
	
	
	//performing x3-x0
	
	setCycleFromSignal("signal_x0");
	
	vhdl<<tab<<declare("segmentX3mX0",inputWidth)<<"<="<<use("signal_x3")<<" - "<<use("signal_x0")<<";"<<endl;
	
	
	//performing y1-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("segmentY1mY0",inputWidth)<<"<="<<use("signal_y1")<<" - "<<use("signal_y0")<<";"<<endl;
	
	//performing y0-y1
	
	setCycleFromSignal("signal_y1");
	
	vhdl<<tab<<declare("segmentY0mY1",inputWidth)<<"<="<<use("signal_y0")<<" - "<<use("signal_y1")<<";"<<endl;
	
	//performing y0-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("segmentY0mY2",inputWidth)<<"<="<<use("signal_y0")<<" - "<<use("signal_y2")<<";"<<endl;
	
	
	
	//performing y3-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("segmentY3mY2",inputWidth)<<"<="<<use("signal_y3")<<" - "<<use("signal_y2")<<";"<<endl;
	
	
	//performing y3-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("segmentY3mY0",inputWidth)<<"<="<<use("signal_y3")<<" - "<<use("signal_y0")<<";"<<endl;
	
	
	//performing z1-z0
	
	setCycleFromSignal("signal_z0");

	vhdl<<tab<<declare("segmentZ1mZ0",inputWidth)<<"<="<<use("signal_z1")<<" - "<<use("signal_z0")<<";"<<endl;

		
	//performing z0-z1
	
	setCycleFromSignal("signal_z1");
	
	vhdl<<tab<<declare("segmentZ0mZ1",inputWidth)<<"<="<<use("signal_z0")<<" - "<<use("signal_z1")<<";"<<endl;
		
	
	//performing z0-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("segmentZ0mZ2",inputWidth)<<"<="<<use("signal_z0")<<" - "<<use("signal_z2")<<";"<<endl;
	
	
	//performing z3-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("segmentZ3mZ2",inputWidth)<<"<="<<use("signal_z3")<<" - "<<use("signal_z2")<<";"<<endl;
	
	
	//performing z3-z0
	
	setCycleFromSignal("signal_z0");
	
	vhdl<<tab<<declare("segmentZ3mZ0",inputWidth)<<"<="<<use("signal_z3")<<" - "<<use("signal_z0")<<";"<<endl;
	
	nextCycle();
	
	//referenceCycle1 = getCurrentCycle();
	
	//computing the var1 (x1-x0)*(x3-x2)+(y1-y0)*(y3-y2)+(z1-z0)*(z3-z2)
	
	
		//(x1-x0)*(x3-x2)
	//setCycleFromSignal("signX1mX0v1");
	//setCycle(referenceCycle1);
	
	vhdl<<tab<<declare("signX1mX0v1",1)<<"<= "<<use("segmentX1mX0")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorX1mX0v1",inputWidth-1)<<"<= ( others=> "<<use("signX1mX0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveX1mX0v1",inputWidth-1)<<"<="<<use("sign2vectorX1mX0v1")<<" xor "<<use("segmentX1mX0")<<range(inputWidth-2,0)<<";"<<endl;
	
	
	vhdl<<declare("converted2PositiveX1mX0v1",inputWidth-1)<<"<="<<use("partialConverted2PositiveX1mX0v1") << " + "<< use("signX1mX0v1")<<";"<<endl;
	
	
	
	vhdl<<tab<<declare("signX3mX2v1",1)<<"<= "<<use("segmentX3mX2")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorX3mX2v1",inputWidth-1)<<"<= ( others=> "<<use("signX3mX2v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveX3mX2v1",inputWidth-1)<<"<="<<use("sign2vectorX3mX2v1")<<" xor "<<use("segmentX3mX2")<<range(inputWidth-2,0)<<";"<<endl;
	
	vhdl<<tab<<declare("converted2PositiveX3mX2v1",inputWidth-1)<<"<="<<use("partialConverted2PositiveX3mX2v1") << " + "<< use("signX3mX2v1")<<";"<<endl;
	
	//referenceCycle2 = getCurrentCycle();
	
	vhdl<<tab<<declare("signSegmentXv1",1)<<" <= "<<use("signX1mX0v1")<<" xor "<<use("signX3mX2v1")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("converted2PositiveX1mX0v1temp", inputWidth-1)<<"<="<<use("converted2PositiveX1mX0v1")<<";"<<endl;
	vhdl<<tab<<declare("converted2PositiveX3mX2v1temp", inputWidth-1)<<"<="<<use("converted2PositiveX3mX2v1")<<";"<<endl;
	
	multiplierXSegv1 = new IntMultiplier(target, inputWidth-1, inputWidth-1);
	multiplierXSegv1->changeName(getName()+"multiplierXSegv1");
	oplist.push_back(multiplierXSegv1);
	inPortMap  (multiplierXSegv1, "X", use("converted2PositiveX1mX0v1temp"));
	inPortMap  (multiplierXSegv1, "Y", use("converted2PositiveX3mX2v1temp"));
	outPortMap (multiplierXSegv1, "R","partialComputedProductSX1v1");
	vhdl << instance(multiplierXSegv1, "multiplierXSegv1");
	
	syncCycleFromSignal("partialComputedProductSX1v1");
	
	inputWidthSegments=2*(inputWidth-1)+1;
	vhdl<<tab<<declare("partialComputedProductSX1v1Bit",inputWidthSegments)<<"<= '0' & "<<use("partialComputedProductSX1v1")<<";"<<endl;
	vhdl<<tab<<declare("sign2VectorSegmentXv1",inputWidthSegments)<<"<= (others => "<<use("signSegmentXv1")<<");"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSX1v1",inputWidthSegments)<<"<="<<use("partialComputedProductSX1v1Bit")<<" xor "<<use("sign2VectorSegmentXv1")<<";"<<endl;
	
	
	vhdl<<tab<<declare("convertedSegmentXv1",inputWidthSegments)<<"<="<<use("partialConvertedProductSX1v1")<<" + "<<use("signSegmentXv1") <<";"<<endl;
	
		
	
		
	//(y1-y0)*(y3-y2)
	setCycleFromSignal("signX1mX0v1");
	//setCycle(referenceCycle1);
	
	vhdl<<tab<<declare("signY1mY0v1",1)<<"<= "<<use("segmentY1mY0")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorY1mY0v1",inputWidth-1)<<"<= ( others=> "<<use("signY1mY0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveY1mY0v1",inputWidth-1)<<"<="<<use("sign2vectorY1mY0v1")<<" xor "<<use("segmentY1mY0")<<range(inputWidth-2,0)<<";"<<endl;
		
	vhdl<<declare("converted2PositiveY1mY0v1",inputWidth-1)<<"<="<<use("partialConverted2PositiveY1mY0v1") << " + "<< use("signY1mY0v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signY3mY2v1",1)<<"<= "<<use("segmentY3mY2")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorY3mY2v1",inputWidth-1)<<"<= ( others=> "<<use("signY3mY2v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveY3mY2v1",inputWidth-1)<<"<="<<use("sign2vectorY3mY2v1")<<" xor "<<use("segmentY3mY2")<<range(inputWidth-2,0)<<";"<<endl;
	
	vhdl<<tab<<declare("converted2PositiveY3mY2v1",inputWidth-1)<<"<="<<use("partialConverted2PositiveY3mY2v1") << " + "<< use("signY3mY2v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signSegmentYv1",1)<<" <= "<<use("signY1mY0v1")<<" xor "<<use("signY3mY2v1")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("converted2PositiveY1mY0v1temp", inputWidth-1)<<"<="<<use("converted2PositiveY1mY0v1")<<";"<<endl;
	vhdl<<tab<<declare("converted2PositiveY3mY2v1temp", inputWidth-1)<<"<="<<use("converted2PositiveY3mY2v1")<<";"<<endl;
	
	multiplierYSegv1 = new IntMultiplier(target, inputWidth-1, inputWidth-1);
	multiplierYSegv1->changeName(getName()+"multiplierYSegv1");
	oplist.push_back(multiplierYSegv1);
	inPortMap  (multiplierYSegv1, "X", use("converted2PositiveY1mY0v1temp"));
	inPortMap  (multiplierYSegv1, "Y", use("converted2PositiveY3mY2v1temp"));
	outPortMap (multiplierYSegv1, "R","partialComputedProductSY1v1");
	vhdl << instance(multiplierYSegv1, "multiplierYSegv1");
	
	syncCycleFromSignal("partialComputedProductSY1v1");
	
	inputWidthSegments=2*(inputWidth-1)+1;
	vhdl<<tab<<declare("partialComputedProductSY1v1Bit",inputWidthSegments)<<"<= '0' & "<<use("partialComputedProductSY1v1")<<";"<<endl;
	vhdl<<tab<<declare("sign2VectorSegmentYv1",inputWidthSegments)<<"<= (others => "<<use("signSegmentYv1")<<");"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSY1v1",inputWidthSegments)<<"<="<<use("partialComputedProductSY1v1Bit")<<" xor "<<use("sign2VectorSegmentYv1")<<";"<<endl;
	
	vhdl<<tab<<declare("convertedSegmentYv1",inputWidthSegments)<<"<="<<use("partialConvertedProductSY1v1")<<" + "<<use("signSegmentYv1") <<";"<<endl;
	
	
	
			
		//(z1-z0)*(z3-z2)
	setCycleFromSignal("signX1mX0v1");
	//setCycle(referenceCycle1);
	
	vhdl<<tab<<declare("signZ1mZ0v1",1)<<"<= "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorZ1mZ0v1",inputWidth-1)<<"<= ( others=> "<<use("signZ1mZ0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveZ1mZ0v1",inputWidth-1)<<"<="<<use("sign2vectorZ1mZ0v1")<<" xor "<<use("segmentZ1mZ0")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("zeroInput4ConvertionZv1",inputWidth-1)<<"<= (others => '0' );"<<endl;
	
	vhdl<<declare("converted2PositiveZ1mZ0v1",inputWidth-1)<<"<="<<use("partialConverted2PositiveZ1mZ0v1") << " + "<< use("signZ1mZ0v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signZ3mZ2v1",1)<<"<= "<<use("segmentZ3mZ2")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorZ3mZ2v1",inputWidth-1)<<"<= ( others=> "<<use("signZ3mZ2v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveZ3mZ2v1",inputWidth-1)<<"<="<<use("sign2vectorZ3mZ2v1")<<" xor "<<use("segmentZ3mZ2")<<range(inputWidth-2,0)<<";"<<endl;
	
	vhdl<<tab<<declare("converted2PositiveZ3mZ2v1",inputWidth-1)<<"<="<<use("partialConverted2PositiveZ3mZ2v1") << " + "<< use("signZ3mZ2v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signSegmentZv1",1)<<" <= "<<use("signZ1mZ0v1")<<" xor "<<use("signZ3mZ2v1")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("converted2PositiveZ1mZ0v1temp", inputWidth-1)<<"<="<<use("converted2PositiveZ1mZ0v1")<<";"<<endl;
	vhdl<<tab<<declare("converted2PositiveZ3mZ2v1temp", inputWidth-1)<<"<="<<use("converted2PositiveZ3mZ2v1")<<";"<<endl;
	
	multiplierZSegv1 = new IntMultiplier(target, inputWidth-1, inputWidth-1);
	multiplierZSegv1->changeName(getName()+"multiplierZSegv1");
	oplist.push_back(multiplierZSegv1);
	inPortMap  (multiplierZSegv1, "X", use("converted2PositiveZ1mZ0v1temp"));
	inPortMap  (multiplierZSegv1, "Y", use("converted2PositiveZ3mZ2v1temp"));
	outPortMap (multiplierZSegv1, "R","partialComputedProductSZ1v1");
	vhdl << instance(multiplierZSegv1, "multiplierZSegv1");
	
	syncCycleFromSignal("partialComputedProductSZ1v1");
	
	inputWidthSegments=2*(inputWidth-1)+1;
	vhdl<<tab<<declare("partialComputedProductSZ1v1Bit",inputWidthSegments)<<"<= '0' & "<<use("partialComputedProductSZ1v1")<<";"<<endl;
	vhdl<<tab<<declare("sign2VectorSegmentZv1",inputWidthSegments)<<"<= (others => "<<use("signSegmentZv1")<<");"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZ1v1",inputWidthSegments)<<"<="<<use("partialComputedProductSZ1v1Bit")<<" xor "<<use("sign2VectorSegmentZv1")<<";"<<endl;
	
	vhdl<<tab<<declare("convertedSegmentZv1",inputWidthSegments)<<"<="<<use("partialConvertedProductSZ1v1")<<" + "<<use("signSegmentZv1") <<";"<<endl;
	
	
	nextCycle();
	
		//ading the 3 results for segments
	
	vhdl<<tab<<declare("convertedSegmentYv1temp",inputWidthSegments)<<"<="<<use("convertedSegmentYv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentZv1temp",inputWidthSegments)<<"<="<<use("convertedSegmentZv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentXv1temp",inputWidthSegments)<<"<="<<use("convertedSegmentXv1")<<";"<<endl;
	
	adder4var1 = new IntNAdder(target,inputWidthSegments,3);
	adder4var1->changeName(getName()+"adder4var1");
	oplist.push_back(adder4var1);
	inPortMap  (adder4var1, "X0", use("convertedSegmentYv1temp"));
	inPortMap  (adder4var1, "X1", use("convertedSegmentZv1temp"));
	inPortMap  (adder4var1, "X2", use("convertedSegmentXv1temp") );
	inPortMapCst(adder4var1,"Cin","'0'");
	outPortMap (adder4var1, "R","result4Var1");
	vhdl << instance(adder4var1, "adder4var1");
	
	syncCycleFromSignal("result4Var1");
	
	//int signofLSBI=LSBI>=0?(+1):(-1);
	int signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2<<endl;
	
	convert2FPv1 = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2,1,wE,wF);
	convert2FPv1->changeName(getName()+"convert2FPv1");
	oplist.push_back(convert2FPv1);
	inPortMap  (convert2FPv1, "I", use("result4Var1"));
	outPortMap (convert2FPv1, "O","Var1");
	vhdl << instance(convert2FPv1, "convert2FPv1");
	
	syncCycleFromSignal("Var1");
	
		
	//computing the var2 sqrt( (x3-x2)^2 +  (y3-y2)^2 +(z3-z2)^2 )
	
		//computing (x3-x2)^2
	
	setCycleFromSignal("converted2PositiveX3mX2v1temp");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("segmentX3mX2ns",inputWidth-1)<<"<="<<use("converted2PositiveX3mX2v1")<<";"<<endl;
	
	
	sqrXv2 = new IntSquarer(target, inputWidth-1);
	sqrXv2->changeName(getName()+"sqrXv2");
	oplist.push_back(sqrXv2);
	inPortMap  (sqrXv2, "X", use("segmentX3mX2ns"));
	outPortMap (sqrXv2, "R","sqrX3mX2ns");
	vhdl << instance(sqrXv2, "sqrXv2");
	
	syncCycleFromSignal("sqrX3mX2ns");
	
	
		//computing (y3-y2)^2
	
	setCycleFromSignal("converted2PositiveY3mY2v1temp");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("segmentY3mY2ns",inputWidth-1)<<"<="<<use("converted2PositiveY3mY2v1")<<";"<<endl;
	
	
	sqrYv2 = new IntSquarer(target, inputWidth-1);
	sqrYv2->changeName(getName()+"sqrYv2");
	oplist.push_back(sqrYv2);
	inPortMap  (sqrYv2, "X", use("segmentY3mY2ns"));
	outPortMap (sqrYv2, "R","sqrY3mY2ns");
	vhdl << instance(sqrYv2, "sqrYv2");
	
	syncCycleFromSignal("sqrY3mY2ns");
	
		//computing (z3-z2)^2
	
	setCycleFromSignal("converted2PositiveZ3mZ2v1temp");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("segmentZ3mZ2ns",inputWidth-1)<<"<="<<use("converted2PositiveZ3mZ2v1")<<";"<<endl;
	
	
	sqrZv2 = new IntSquarer(target, inputWidth-1);
	sqrZv2->changeName(getName()+"sqrZv2");
	oplist.push_back(sqrZv2);
	inPortMap  (sqrZv2, "X", use("segmentZ3mZ2ns"));
	outPortMap (sqrZv2, "R","sqrZ3mZ2ns");
	vhdl << instance(sqrZv2, "sqrZv2");
	
	syncCycleFromSignal("sqrZ3mZ2ns");
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	
	adder4SQRTv2 = new IntNAdder(target,(inputWidth-1)*2,3);
	adder4SQRTv2->changeName(getName()+"adder4SQRTv2");
	oplist.push_back(adder4SQRTv2);
	inPortMap  (adder4SQRTv2, "X0", use("sqrY3mY2ns"));
	inPortMap  (adder4SQRTv2, "X1", use("sqrZ3mZ2ns"));
	inPortMap  (adder4SQRTv2, "X2", use("sqrX3mX2ns") );
	inPortMapCst(adder4SQRTv2,"Cin","'0'");
	outPortMap (adder4SQRTv2, "R","result4SQRTv2");
	vhdl << instance(adder4SQRTv2, "adder4SQRTv2");
	
	syncCycleFromSignal("result4SQRTv2");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	convert2FP4sqrtv2 = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2-1,0,wE,wF);
	convert2FP4sqrtv2->changeName(getName()+"convert2FP4sqrtv2");
	oplist.push_back(convert2FP4sqrtv2);
	inPortMap  (convert2FP4sqrtv2, "I", use("result4SQRTv2"));
	outPortMap (convert2FP4sqrtv2, "O","fpSQRVar2");
	vhdl << instance(convert2FP4sqrtv2, "convert2FP4sqrtv2");
	
	syncCycleFromSignal("fpSQRVar2");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar2temp",wE+wF+3)<<"<="<<use("fpSQRVar2")<<";"<<endl;
	
	sqrt4var2 = new  FPSqrt(target, wE, wF, 1, 1);
	sqrt4var2->changeName(getName()+"sqrt4var2");
	oplist.push_back(sqrt4var2);
	inPortMap  (sqrt4var2, "X", use("fpSQRVar2temp"));
	outPortMap (sqrt4var2, "R","Var2");
	vhdl << instance(sqrt4var2, "sqrt4var2");
	
	syncCycleFromSignal("Var2");
	
	
	//computing the var3 sart(((x3-x0)+(x0-x1)*t)^2+((y3-y0)+(y0-y1)*t)^2+((z3-z0)+(z0-z1)*t)^2)
	
	
	
	
	
	
		//(x3-x0)+(x0-x1)*t
	setCycleFromSignal("signX1mX0v1");
	
	vhdl<<tab<<declare("signX0mX1",1)<<"<="<<use("segmentX0mX1")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorX0mX1",inputWidth-1)<<" <= (others => "<<use("signX0mX1")<<");"<<endl;
	
	vhdl<<tab<<declare("segmentX0mX1PartialPositive",inputWidth-1)<<"<="<<use("sign2vectorX0mX1")<<" xor "<<use("segmentX0mX1")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("segmentX0mX12Positive",inputWidth-1)<<"<="<<use("segmentX0mX1PartialPositive")<<" + "<<use("signX0mX1")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("segmentX0mX12Positivetemp",inputWidth-1)<<"<="<<use("segmentX0mX12Positive")<<";"<<endl;
	vhdl<<tab<<declare("signal_t_temp",integratorWidth)<<"<="<<use("signal_t")<<";"<<endl;
	
	multiplierXv3 = new IntMultiplier(target, inputWidth-1,integratorWidth);
	multiplierXv3->changeName(getName()+"multiplierXv3");
	oplist.push_back(multiplierXv3);
	inPortMap  (multiplierXv3, "X", use("segmentX0mX12Positivetemp"));
	inPortMap  (multiplierXv3, "Y", use("signal_t_temp"));
	outPortMap (multiplierXv3, "R","partialComputedProductSXv3");
	vhdl << instance(multiplierXv3, "multiplierXv3");
	
	syncCycleFromSignal("partialComputedProductSXv3");
	
	vhdl<<tab<<declare("sign2vectorX0mX1_2",inputWidth)<<"<= (others => "<<use("signX0mX1")<<");"<<endl;
	vhdl<<tab<<declare("partialComputedProductSXv3Bit",inputWidth)<<"<= '0' & "<<use("partialComputedProductSXv3")<<range(inputWidth-2+integratorWidth,integratorWidth)<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSXv3",inputWidth)<<"<="<<use("partialComputedProductSXv3Bit")<<" xor "<< use("sign2vectorX0mX1_2")<<";"<<endl;

	vhdl<<tab<<declare("sumXPartv3",inputWidth)<<"<="<<use("partialConvertedProductSXv3")<<" + "<<use("segmentX3mX0")<<" + "<<use("signX0mX1")<<";"<<endl;
	
	vhdl<<tab<<declare("signsumXPartv3",1)<<"<="<<use("sumXPartv3")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vector4sumXPartv3",inputWidth-1)<<"<= (others => "<<use("signsumXPartv3")<<");"<<endl;
	vhdl<<tab<<declare("sumX2PartialPositivePartv3",inputWidth-1)<<"<="<<use("sumXPartv3")<<range(inputWidth-2,0)<<" xor "<<use("sign2vector4sumXPartv3")<<";"<<endl;
	vhdl<<tab<<declare("sumX2PositivePartv3",inputWidth-1)<<"<="<<use("sumX2PartialPositivePartv3")<<" + "<<use("signsumXPartv3")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sumX2PositivePartv3temp",inputWidth-1)<<"<="<<use("sumX2PositivePartv3")<<";"<<endl;
	
	sqrXv3 = new IntSquarer(target, inputWidth-1);
	sqrXv3->changeName(getName()+"sqrXv3");
	oplist.push_back(sqrXv3);
	inPortMap  (sqrXv3, "X", use("sumX2PositivePartv3temp"));
	outPortMap (sqrXv3, "R","sqrXnsv3");
	vhdl << instance(sqrXv3, "sqrXv3");
	
	syncCycleFromSignal("sqrXnsv3");
	
	
	
		//(y3-y0)+(y0-y1)*t
	setCycleFromSignal("signY1mY0v1");
	
	vhdl<<tab<<declare("signY0mY1",1)<<"<="<<use("segmentY0mY1")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorY0mY1",inputWidth-1)<<" <= (others => "<<use("signY0mY1")<<");"<<endl;
	
	vhdl<<tab<<declare("segmentY0mY1PartialPositive",inputWidth-1)<<"<="<<use("sign2vectorY0mY1")<<" xor "<<use("segmentY0mY1")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("segmentY0mY12Positive",inputWidth-1)<<"<="<<use("segmentY0mY1PartialPositive")<<" + "<<use("signY0mY1")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("segmentY0mY12Positivetemp",inputWidth-1)<<"<="<<use("segmentY0mY12Positive")<<";"<<endl;
	//vhdl<<tab<<declare("signal_t_temp",integratorWidth)<<"<="<<use("signal_t")<<";"<<endl;
	
	multiplierYv3 = new IntMultiplier(target, inputWidth-1,integratorWidth);
	multiplierYv3->changeName(getName()+"multiplierYv3");
	oplist.push_back(multiplierYv3);
	inPortMap  (multiplierYv3, "X", use("segmentY0mY12Positivetemp"));
	inPortMap  (multiplierYv3, "Y", use("signal_t_temp"));
	outPortMap (multiplierYv3, "R","partialComputedProductSYv3");
	vhdl << instance(multiplierYv3, "multiplierYv3");
	
	syncCycleFromSignal("partialComputedProductSYv3");
	
	vhdl<<tab<<declare("sign2vectorY0mY1_2",inputWidth)<<"<= (others => "<<use("signY0mY1")<<");"<<endl;
	vhdl<<tab<<declare("partialComputedProductSYv3Bit",inputWidth)<<"<= '0' & "<<use("partialComputedProductSYv3")<<range(inputWidth-2+integratorWidth,integratorWidth)<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSYv3",inputWidth)<<"<="<<use("partialComputedProductSYv3Bit")<<" xor "<< use("sign2vectorY0mY1_2")<<";"<<endl;
		
	vhdl<<tab<<declare("sumYPartv3",inputWidth)<<"<="<<use("partialConvertedProductSYv3")<<" + "<<use("segmentY3mY0")<<" + "<<use("signY0mY1")<<";"<<endl;
	
	vhdl<<tab<<declare("signsumYPartv3",1)<<"<="<<use("sumYPartv3")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vector4sumYPartv3",inputWidth-1)<<"<= (others => "<<use("signsumYPartv3")<<");"<<endl;
	vhdl<<tab<<declare("sumY2PartialPositivePartv3",inputWidth-1)<<"<="<<use("sumYPartv3")<<range(inputWidth-2,0)<<" xor "<<use("sign2vector4sumYPartv3")<<";"<<endl;
	vhdl<<tab<<declare("sumY2PositivePartv3",inputWidth-1)<<"<="<<use("sumY2PartialPositivePartv3")<<" + "<<use("signsumYPartv3")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sumY2PositivePartv3temp",inputWidth-1)<<"<="<<use("sumY2PositivePartv3")<<";"<<endl;
	
	sqrYv3 = new IntSquarer(target, inputWidth-1);
	sqrYv3->changeName(getName()+"sqrYv3");
	oplist.push_back(sqrYv3);
	inPortMap  (sqrYv3, "X", use("sumY2PositivePartv3temp"));
	outPortMap (sqrYv3, "R","sqrYnsv3");
	vhdl << instance(sqrYv3, "sqrYv3");
	
	syncCycleFromSignal("sqrYnsv3");
	
	
	
	//(z3-z0)+(z0-z1)*t
	setCycleFromSignal("signZ1mZ0v1");
	
	vhdl<<tab<<declare("signZ0mZ1",1)<<"<="<<use("segmentZ0mZ1")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorZ0mZ1",inputWidth-1)<<" <= (others => "<<use("signZ0mZ1")<<");"<<endl;
	
	vhdl<<tab<<declare("segmentZ0mZ1PartialPositive",inputWidth-1)<<"<="<<use("sign2vectorZ0mZ1")<<" xor "<<use("segmentZ0mZ1")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("segmentZ0mZ12Positive",inputWidth-1)<<"<="<<use("segmentZ0mZ1PartialPositive")<<" + "<<use("signZ0mZ1")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("segmentZ0mZ12Positivetemp",inputWidth-1)<<"<="<<use("segmentZ0mZ12Positive")<<";"<<endl;
	//vhdl<<tab<<declare("signal_t_temp",integratorWidth)<<"<="<<use("signal_t")<<";"<<endl;
	
	multiplierZv3 = new IntMultiplier(target, inputWidth-1,integratorWidth);
	multiplierZv3->changeName(getName()+"multiplierZv3");
	oplist.push_back(multiplierZv3);
	inPortMap  (multiplierZv3, "X", use("segmentZ0mZ12Positivetemp"));
	inPortMap  (multiplierZv3, "Y", use("signal_t_temp"));
	outPortMap (multiplierZv3, "R","partialComputedProductSZv3");
	vhdl << instance(multiplierZv3, "multiplierZv3");
	
	syncCycleFromSignal("partialComputedProductSZv3");
	
	vhdl<<tab<<declare("sign2vectorZ0mZ1_2",inputWidth)<<"<= (others => "<<use("signZ0mZ1")<<");"<<endl;
	vhdl<<tab<<declare("partialComputedProductSZv3Bit",inputWidth)<<"<= '0' & "<<use("partialComputedProductSZv3")<<range(inputWidth-2+integratorWidth,integratorWidth)<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZv3",inputWidth)<<"<="<<use("partialComputedProductSZv3Bit")<<" xor "<< use("sign2vectorZ0mZ1_2")<<";"<<endl;
		
	vhdl<<tab<<declare("sumZPartv3",inputWidth)<<"<="<<use("partialConvertedProductSZv3")<<" + "<<use("segmentZ3mZ0")<<" + "<<use("signZ0mZ1")<<";"<<endl;
	
	vhdl<<tab<<declare("signsumZPartv3",1)<<"<="<<use("sumZPartv3")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vector4sumZPartv3",inputWidth-1)<<"<= (others => "<<use("signsumZPartv3")<<");"<<endl;
	vhdl<<tab<<declare("sumZ2PartialPositivePartv3",inputWidth-1)<<"<="<<use("sumZPartv3")<<range(inputWidth-2,0)<<" xor "<<use("sign2vector4sumZPartv3")<<";"<<endl;
	vhdl<<tab<<declare("sumZ2PositivePartv3",inputWidth-1)<<"<="<<use("sumZ2PartialPositivePartv3")<<" + "<<use("signsumZPartv3")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sumZ2PositivePartv3temp",inputWidth-1)<<"<="<<use("sumZ2PositivePartv3")<<";"<<endl;
	
	sqrZv3 = new IntSquarer(target, inputWidth-1);
	sqrZv3->changeName(getName()+"sqrZv3");
	oplist.push_back(sqrZv3);
	inPortMap  (sqrZv3, "X", use("sumZ2PositivePartv3temp"));
	outPortMap (sqrZv3, "R","sqrZnsv3");
	vhdl << instance(sqrZv3, "sqrZv3");
	
	syncCycleFromSignal("sqrZnsv3");
	
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	
	adder4SQRTv3 = new IntNAdder(target,(inputWidth-1)*2,3);
	adder4SQRTv3->changeName(getName()+"adder4SQRTv3");
	oplist.push_back(adder4SQRTv3);
	inPortMap  (adder4SQRTv3, "X0", use("sqrYnsv3"));
	inPortMap  (adder4SQRTv3, "X1", use("sqrZnsv3"));
	inPortMap  (adder4SQRTv3, "X2", use("sqrXnsv3") );
	inPortMapCst(adder4SQRTv3,"Cin","'0'");
	outPortMap (adder4SQRTv3, "R","result4SQRTv3");
	vhdl << instance(adder4SQRTv3, "adder4SQRTv3");
	
	syncCycleFromSignal("result4SQRTv3");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	convert2FP4sqrtv3 = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2-1,0,wE,wF);
	convert2FP4sqrtv3->changeName(getName()+"convert2FP4sqrtv3");
	oplist.push_back(convert2FP4sqrtv3);
	inPortMap  (convert2FP4sqrtv3, "I", use("result4SQRTv3"));
	outPortMap (convert2FP4sqrtv3, "O","fpSQRVar3");
	vhdl << instance(convert2FP4sqrtv3, "convert2FP4sqrtv3");
	
	syncCycleFromSignal("fpSQRVar3");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar3temp",wE+wF+3)<<"<="<<use("fpSQRVar3")<<";"<<endl;
	
	sqrt4var3 = new  FPSqrt(target, wE, wF, 1, 1);
	sqrt4var3->changeName(getName()+"sqrt4var3");
	oplist.push_back(sqrt4var3);
	inPortMap  (sqrt4var3, "X", use("fpSQRVar3temp"));
	outPortMap (sqrt4var3, "R","Var3");
	vhdl << instance(sqrt4var3, "sqrt4var3");
	
	syncCycleFromSignal("Var3");
	
	
	//computing the var4 sart(((x0-x2)+(x1-x0)*t)^2+((y0-y2)+(y1-y0)*t)^2+((z0-z2)+(z1-z0)*t)^2)
	
	
	
	
	
		//((x0-x2)+(x1-x0)*t)^2
	setCycleFromSignal("converted2PositiveX1mX0v1temp");
	
	
		
	vhdl<<tab<<declare("segmentX1mX02Positivetemp",inputWidth-1)<<"<="<<use("converted2PositiveX1mX0v1")<<";"<<endl;
	//vhdl<<tab<<declare("signal_t_temp",integratorWidth)<<"<="<<use("signal_t")<<";"<<endl;
	
	multiplierXv4 = new IntMultiplier(target, inputWidth-1,integratorWidth);
	multiplierXv4->changeName(getName()+"multiplierXv4");
	oplist.push_back(multiplierXv4);
	inPortMap  (multiplierXv4, "X", use("segmentX1mX02Positivetemp"));
	inPortMap  (multiplierXv4, "Y", use("signal_t_temp"));
	outPortMap (multiplierXv4, "R","partialComputedProductSXv4");
	vhdl << instance(multiplierXv4, "multiplierXv4");
	
	syncCycleFromSignal("partialComputedProductSXv4");
	
	vhdl<<tab<<declare("sign2vectorX1mX0_2",inputWidth)<<"<= (others => "<<use("signX1mX0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialComputedProductSXv4Bit",inputWidth)<<"<= '0' & "<<use("partialComputedProductSXv4")<<range(inputWidth-2+integratorWidth,integratorWidth)<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSXv4",inputWidth)<<"<="<<use("partialComputedProductSXv4Bit")<<" xor "<< use("sign2vectorX1mX0_2")<<";"<<endl;

	vhdl<<tab<<declare("sumXPartv4",inputWidth)<<"<="<<use("partialConvertedProductSXv4")<<" + "<<use("segmentX0mX2")<<" + "<<use("signX1mX0v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signsumXPartv4",1)<<"<="<<use("sumXPartv4")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vector4sumXPartv4",inputWidth-1)<<"<= (others => "<<use("signsumXPartv4")<<");"<<endl;
	vhdl<<tab<<declare("sumX2PartialPositivePartv4",inputWidth-1)<<"<="<<use("sumXPartv4")<<range(inputWidth-2,0)<<" xor "<<use("sign2vector4sumXPartv4")<<";"<<endl;
	vhdl<<tab<<declare("sumX2PositivePartv4",inputWidth-1)<<"<="<<use("sumX2PartialPositivePartv4")<<" + "<<use("signsumXPartv4")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sumX2PositivePartv4temp",inputWidth-1)<<"<="<<use("sumX2PositivePartv4")<<";"<<endl;
	
	sqrXv4 = new IntSquarer(target, inputWidth-1);
	sqrXv4->changeName(getName()+"sqrXv4");
	oplist.push_back(sqrXv4);
	inPortMap  (sqrXv4, "X", use("sumX2PositivePartv4temp"));
	outPortMap (sqrXv4, "R","sqrXnsv4");
	vhdl << instance(sqrXv4, "sqrXv4");
	
	syncCycleFromSignal("sqrXnsv4");
	
	
	
		//((y0-y2)+(y1-y0)*t)^2
	setCycleFromSignal("converted2PositiveY1mY0v1temp");
	
			
	vhdl<<tab<<declare("segmentY1mY02Positivetemp",inputWidth-1)<<"<="<<use("converted2PositiveY1mY0v1")<<";"<<endl;
	//vhdl<<tab<<declare("signal_t_temp",integratorWidth)<<"<="<<use("signal_t")<<";"<<endl;
	
	multiplierYv4 = new IntMultiplier(target, inputWidth-1,integratorWidth);
	multiplierYv4->changeName(getName()+"multiplierYv4");
	oplist.push_back(multiplierYv4);
	inPortMap  (multiplierYv4, "X", use("segmentY1mY02Positivetemp"));
	inPortMap  (multiplierYv4, "Y", use("signal_t_temp"));
	outPortMap (multiplierYv4, "R","partialComputedProductSYv4");
	vhdl << instance(multiplierYv4, "multiplierYv4");
	
	syncCycleFromSignal("partialComputedProductSYv4");
	
	vhdl<<tab<<declare("sign2vectorY1mY0_2",inputWidth)<<"<= (others => "<<use("signY1mY0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialComputedProductSYv4Bit",inputWidth)<<"<= '0' & "<<use("partialComputedProductSYv4")<<range(inputWidth-2+integratorWidth,integratorWidth)<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSYv4",inputWidth)<<"<="<<use("partialComputedProductSYv4Bit")<<" xor "<< use("sign2vectorY1mY0_2")<<";"<<endl;

	vhdl<<tab<<declare("sumYPartv4",inputWidth)<<"<="<<use("partialConvertedProductSYv4")<<" + "<<use("segmentY0mY2")<<" + "<<use("signY1mY0v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signsumYPartv4",1)<<"<="<<use("sumYPartv4")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vector4sumYPartv4",inputWidth-1)<<"<= (others => "<<use("signsumYPartv4")<<");"<<endl;
	vhdl<<tab<<declare("sumY2PartialPositivePartv4",inputWidth-1)<<"<="<<use("sumYPartv4")<<range(inputWidth-2,0)<<" xor "<<use("sign2vector4sumYPartv4")<<";"<<endl;
	vhdl<<tab<<declare("sumY2PositivePartv4",inputWidth-1)<<"<="<<use("sumY2PartialPositivePartv4")<<" + "<<use("signsumYPartv4")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sumY2PositivePartv4temp",inputWidth-1)<<"<="<<use("sumY2PositivePartv4")<<";"<<endl;
	
	sqrYv4 = new IntSquarer(target, inputWidth-1);
	sqrYv4->changeName(getName()+"sqrYv4");
	oplist.push_back(sqrYv4);
	inPortMap  (sqrYv4, "X", use("sumY2PositivePartv4temp"));
	outPortMap (sqrYv4, "R","sqrYnsv4");
	vhdl << instance(sqrYv4, "sqrYv4");
	
	syncCycleFromSignal("sqrYnsv4");
	
	
			//((z0-z2)+(z1-z0)*t)^2
	setCycleFromSignal("converted2PositiveZ1mZ0v1temp");
	
	
		
	vhdl<<tab<<declare("segmentZ1mZ02Positivetemp",inputWidth-1)<<"<="<<use("converted2PositiveZ1mZ0v1")<<";"<<endl;
	//vhdl<<tab<<declare("signal_t_temp",integratorWidth)<<"<="<<use("signal_t")<<";"<<endl;
	
	multiplierZv4 = new IntMultiplier(target, inputWidth-1,integratorWidth);
	multiplierZv4->changeName(getName()+"multiplierZv4");
	oplist.push_back(multiplierZv4);
	inPortMap  (multiplierZv4, "X", use("segmentZ1mZ02Positivetemp"));
	inPortMap  (multiplierZv4, "Y", use("signal_t_temp"));
	outPortMap (multiplierZv4, "R","partialComputedProductSZv4");
	vhdl << instance(multiplierZv4, "multiplierZv4");
	
	syncCycleFromSignal("partialComputedProductSZv4");
	
	vhdl<<tab<<declare("sign2vectorZ1mZ0_2",inputWidth)<<"<= (others => "<<use("signZ1mZ0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialComputedProductSZv4Bit",inputWidth)<<"<= '0' & "<<use("partialComputedProductSZv4")<<range(inputWidth-2+integratorWidth,integratorWidth)<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZv4",inputWidth)<<"<="<<use("partialComputedProductSZv4Bit")<<" xor "<< use("sign2vectorZ1mZ0_2")<<";"<<endl;

	vhdl<<tab<<declare("sumZPartv4",inputWidth)<<"<="<<use("partialConvertedProductSZv4")<<" + "<<use("segmentZ0mZ2")<<" + "<<use("signZ1mZ0v1")<<";"<<endl;
	
	vhdl<<tab<<declare("signsumZPartv4",1)<<"<="<<use("sumZPartv4")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vector4sumZPartv4",inputWidth-1)<<"<= (others => "<<use("signsumZPartv4")<<");"<<endl;
	vhdl<<tab<<declare("sumZ2PartialPositivePartv4",inputWidth-1)<<"<="<<use("sumZPartv4")<<range(inputWidth-2,0)<<" xor "<<use("sign2vector4sumZPartv4")<<";"<<endl;
	vhdl<<tab<<declare("sumZ2PositivePartv4",inputWidth-1)<<"<="<<use("sumZ2PartialPositivePartv4")<<" + "<<use("signsumZPartv4")<<";"<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sumZ2PositivePartv4temp",inputWidth-1)<<"<="<<use("sumZ2PositivePartv4")<<";"<<endl;
	
	sqrZv4 = new IntSquarer(target, inputWidth-1);
	sqrZv4->changeName(getName()+"sqrZv4");
	oplist.push_back(sqrZv4);
	inPortMap  (sqrZv4, "X", use("sumZ2PositivePartv4temp"));
	outPortMap (sqrZv4, "R","sqrZnsv4");
	vhdl << instance(sqrZv4, "sqrZv4");
	
	syncCycleFromSignal("sqrZnsv4");
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	
	adder4SQRTv4 = new IntNAdder(target,(inputWidth-1)*2,3);
	adder4SQRTv4->changeName(getName()+"adder4SQRTv4");
	oplist.push_back(adder4SQRTv4);
	inPortMap  (adder4SQRTv4, "X0", use("sqrYnsv4"));
	inPortMap  (adder4SQRTv4, "X1", use("sqrZnsv4"));
	inPortMap  (adder4SQRTv4, "X2", use("sqrXnsv4") );
	inPortMapCst(adder4SQRTv4,"Cin","'0'");
	outPortMap (adder4SQRTv4, "R","result4SQRTv4");
	vhdl << instance(adder4SQRTv4, "adder4SQRTv4");
	
	syncCycleFromSignal("result4SQRTv4");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	convert2FP4sqrtv4 = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2-1,0,wE,wF);
	convert2FP4sqrtv4->changeName(getName()+"convert2FP4sqrtv4");
	oplist.push_back(convert2FP4sqrtv4);
	inPortMap  (convert2FP4sqrtv4, "I", use("result4SQRTv4"));
	outPortMap (convert2FP4sqrtv4, "O","fpSQRVar4");
	vhdl << instance(convert2FP4sqrtv4, "convert2FP4sqrtv4");
	
	syncCycleFromSignal("fpSQRVar4");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar4temp",wE+wF+3)<<"<="<<use("fpSQRVar4")<<";"<<endl;
	
	sqrt4var4 = new  FPSqrt(target, wE, wF, 1, 1);
	sqrt4var4->changeName(getName()+"sqrt4var4");
	oplist.push_back(sqrt4var4);
	inPortMap  (sqrt4var4, "X", use("fpSQRVar4temp"));
	outPortMap (sqrt4var4, "R","Var4");
	vhdl << instance(sqrt4var4, "sqrt4var4");
	
	syncCycleFromSignal("Var4");
	
	
	
	
	
	
	
	
	vhdl<<tab<<"O<="<<use("Var1")<<range(outputWidth-1,0)<< "or "<<use("Var2")<<range(outputWidth-1,0)<<" or "<<use("Var3")<<range(outputWidth-1,0)<<" or "<<use("Var4")<<range(outputWidth-1,0)<<";"<<endl;

	}

CoilInductance::~CoilInductance() {
}

void CoilInductance::emulate(TestCase * tc)
{
}

void CoilInductance::buildStandardTestCases(TestCaseList* tcl){
	
}

