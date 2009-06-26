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
	
	addOutput("O",outputWidth);
	
	//Counters for addressing the memories and for frequency division
	
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
	
	vhdl<<tab<<declare("not_x0",inputWidth)<<"<= not( "<<use("signal_x0")<<");"<<endl;
	
	segment1X = new IntAdder(target,inputWidth);
	segment1X->changeName(getName()+"segment1X");
	oplist.push_back(segment1X);
	inPortMap  (segment1X, "X", use("signal_x1"));
	inPortMap  (segment1X, "Y", use("not_x0"));
	inPortMapCst(segment1X, "Cin", "'1'");
	outPortMap (segment1X, "R","segmentX1mX0");
	vhdl << instance(segment1X, "segment1X");
	
	syncCycleFromSignal("segmentX1mX0");
	
	
	//performing x3-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("not_x2",inputWidth)<<"<= not( "<<use("signal_x2")<<");"<<endl;
	
	segment2X = new IntAdder(target,inputWidth);
	segment2X->changeName(getName()+"segment2X");
	oplist.push_back(segment2X);
	inPortMap  (segment2X, "X", use("signal_x3"));
	inPortMap  (segment2X, "Y", use("not_x2"));
	inPortMapCst(segment2X, "Cin", "'1'");
	outPortMap (segment2X, "R","segmentX3mX2");
	vhdl << instance(segment2X, "segment2X");
	
	syncCycleFromSignal("segmentX3mX2");
	
	
	
	//performing y1-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("not_y0",inputWidth)<<"<= not( "<<use("signal_y0")<<");"<<endl;
	
	segment1Y = new IntAdder(target,inputWidth);
	segment1Y->changeName(getName()+"segment1Y");
	oplist.push_back(segment1Y);
	inPortMap  (segment1Y, "X", use("signal_y1"));
	inPortMap  (segment1Y, "Y", use("not_y0"));
	inPortMapCst(segment1Y, "Cin", "'1'");
	outPortMap (segment1Y, "R","segmentY1mY0");
	vhdl << instance(segment1Y, "segment1Y");
	
	syncCycleFromSignal("segmentY1mY0");
	
	
	//performing y3-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("not_y2",inputWidth)<<"<= not( "<<use("signal_y2")<<");"<<endl;
	
	segment2Y = new IntAdder(target,inputWidth);
	segment2Y->changeName(getName()+"segment2Y");
	oplist.push_back(segment2Y);
	inPortMap  (segment2Y, "X", use("signal_y3"));
	inPortMap  (segment2Y, "Y", use("not_y2"));
	inPortMapCst(segment2Y, "Cin", "'1'");
	outPortMap (segment2Y, "R","segmentY3mY2");
	vhdl << instance(segment2Y, "segment2Y");
	
	syncCycleFromSignal("segmentY3mY2");
	
	
	//performing z1-z0
	
	setCycleFromSignal("signal_z0");
	
	vhdl<<tab<<declare("not_z0",inputWidth)<<"<= not( "<<use("signal_z0")<<");"<<endl;
	
	segment1Z = new IntAdder(target,inputWidth);
	segment1Z->changeName(getName()+"segment1Z");
	oplist.push_back(segment1Z);
	inPortMap  (segment1Z, "X", use("signal_z1"));
	inPortMap  (segment1Z, "Y", use("not_z0"));
	inPortMapCst(segment1Z, "Cin", "'1'");
	outPortMap (segment1Z, "R","segmentZ1mZ0");
	vhdl << instance(segment1Z, "segment1Z");
	
	syncCycleFromSignal("segmentZ1mZ0");
	
	
	//performing z3-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("not_z2",inputWidth)<<"<= not( "<<use("signal_z2")<<");"<<endl;
	
	segment2Z = new IntAdder(target,inputWidth);
	segment2Z->changeName(getName()+"segment2Z");
	oplist.push_back(segment2Z);
	inPortMap  (segment2Z, "X", use("signal_z3"));
	inPortMap  (segment2Z, "Y", use("not_z2"));
	inPortMapCst(segment2Z, "Cin", "'1'");
	outPortMap (segment2Z, "R","segmentZ3mZ2");
	vhdl << instance(segment2Z, "segment2Z");
	
	syncCycleFromSignal("segmentZ3mZ2");
	
	
	//computing the var1 (x1-x0)*(x3-x2)+(y1-y0)*(y3-y2)+(z1-z0)*(z3-z2)
	
	
		//(x1-x0)*(x3-x2)
	setCycleFromSignal("segmentX1mX0");
	
	vhdl<<tab<<declare("signX1mX0v1",1)<<"<= "<<use("segmentX1mX0")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorX1mX0v1",inputWidth-1)<<"<= ( others=> "<<use("signX1mX0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveX1mX0v1",inputWidth-1)<<"<="<<use("sign2vectorX1mX0v1")<<" xor "<<use("segmentX1mX0")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("zeroInput4ConvertionXv1",inputWidth-1)<<"<= (others => '0' );"<<endl;
	
	vhdl<<tab<<declare("signX3mX2v1",1)<<"<= "<<use("segmentX3mX2")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorX3mX2v1",inputWidth-1)<<"<= ( others=> "<<use("signX3mX2v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveX3mX2v1",inputWidth-1)<<"<="<<use("sign2vectorX3mX2v1")<<" xor "<<use("segmentX3mX2")<<range(inputWidth-2,0)<<";"<<endl;
	
	
	vhdl<<tab<<declare("signSegmentXv1",1)<<" <= "<<use("signX1mX0v1")<<" xor "<<use("signX3mX2v1")<<";"<<endl;
	
	
	covertInitialXS1v1 = new IntAdder(target,inputWidth-1);
	covertInitialXS1v1->changeName(getName()+"covertInitialXS1v1");
	oplist.push_back(covertInitialXS1v1);
	inPortMap  (covertInitialXS1v1, "X", use("partialConverted2PositiveX1mX0v1"));
	inPortMap  (covertInitialXS1v1, "Y", use("zeroInput4ConvertionXv1"));
	inPortMap  (covertInitialXS1v1, "Cin", use("signX1mX0v1") );
	outPortMap (covertInitialXS1v1, "R","converted2PositiveX1mX0v1");
	vhdl << instance(covertInitialXS1v1, "covertInitialXS1v1");
	
	syncCycleFromSignal("converted2PositiveX1mX0v1");
	
	setCycleFromSignal("zeroInput4ConvertionXv1");
	
	covertInitialXS2v1 = new IntAdder(target,inputWidth-1);
	covertInitialXS2v1->changeName(getName()+"covertInitialXS2v1");
	oplist.push_back(covertInitialXS2v1);
	inPortMap  (covertInitialXS2v1, "X", use("partialConverted2PositiveX3mX2v1"));
	inPortMap  (covertInitialXS2v1, "Y", use("zeroInput4ConvertionXv1"));
	inPortMap  (covertInitialXS2v1, "Cin", use("signX3mX2v1") );
	outPortMap (covertInitialXS2v1, "R","converted2PositiveX3mX2v1");
	vhdl << instance(covertInitialXS2v1, "covertInitialXS2v1");
	
	syncCycleFromSignal("converted2PositiveX3mX2v1");
	
	multiplierXSegv1 = new IntMultiplier(target, inputWidth-1, inputWidth-1);
	multiplierXSegv1->changeName(getName()+"multiplierXSegv1");
	oplist.push_back(multiplierXSegv1);
	inPortMap  (multiplierXSegv1, "X", use("converted2PositiveX1mX0v1"));
	inPortMap  (multiplierXSegv1, "Y", use("converted2PositiveX3mX2v1"));
	outPortMap (multiplierXSegv1, "R","partialComputedProductSX1v1");
	vhdl << instance(multiplierXSegv1, "multiplierXSegv1");
	
	syncCycleFromSignal("partialComputedProductSX1v1");
	
	inputWidthSegments=2*(inputWidth-1)+1;
	vhdl<<tab<<declare("partialComputedProductSX1v1Bit",inputWidthSegments)<<"<= '0' & "<<use("partialComputedProductSX1v1")<<";"<<endl;
	vhdl<<tab<<declare("sign2VectorSegmentXv1",inputWidthSegments)<<"<= (others => "<<use("signSegmentXv1")<<");"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSX1v1",inputWidthSegments)<<"<="<<use("partialComputedProductSX1v1Bit")<<" xor "<<use("sign2VectorSegmentXv1")<<";"<<endl;
	
	//~ if(target->frequencyMHz()>=250)
		//~ nextCycle();
		
	vhdl<<tab<<declare("zeroInput4FinalConvertionXv1",inputWidthSegments)<<"<= (others => '0' );"<<endl;
	vhdl<<tab<<declare("tempsignSegmentXv1",1)<<"<="<<use("signSegmentXv1")<<";"<<endl;
	
	covertFinalXv1 = new IntAdder(target,inputWidthSegments);
	covertFinalXv1->changeName(getName()+"covertFinalXv1");
	oplist.push_back(covertFinalXv1);
	inPortMap  (covertFinalXv1, "X", use("partialConvertedProductSX1v1"));
	inPortMap  (covertFinalXv1, "Y", use("zeroInput4FinalConvertionXv1"));
	inPortMap  (covertFinalXv1, "Cin", use("tempsignSegmentXv1") );
	outPortMap (covertFinalXv1, "R","convertedSegmentXv1");
	vhdl << instance(covertFinalXv1, "covertFinalXv1");
	
	syncCycleFromSignal("convertedSegmentXv1");
	
	
		
	//(y1-y0)*(y3-y2)
	setCycleFromSignal("segmentY1mY0");
	
	vhdl<<tab<<declare("signY1mY0v1",1)<<"<= "<<use("segmentY1mY0")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorY1mY0v1",inputWidth-1)<<"<= ( others=> "<<use("signY1mY0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveY1mY0v1",inputWidth-1)<<"<="<<use("sign2vectorY1mY0v1")<<" xor "<<use("segmentY1mY0")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("zeroInput4ConvertionYv1",inputWidth-1)<<"<= (others => '0' );"<<endl;
	
	vhdl<<tab<<declare("signY3mY2v1",1)<<"<= "<<use("segmentY3mY2")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorY3mY2v1",inputWidth-1)<<"<= ( others=> "<<use("signY3mY2v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveY3mY2v1",inputWidth-1)<<"<="<<use("sign2vectorY3mY2v1")<<" xor "<<use("segmentY3mY2")<<range(inputWidth-2,0)<<";"<<endl;
	
	
	vhdl<<tab<<declare("signSegmentYv1",1)<<" <= "<<use("signY1mY0v1")<<" xor "<<use("signY3mY2v1")<<";"<<endl;
	
	
	covertInitialYS1v1 = new IntAdder(target,inputWidth-1);
	covertInitialYS1v1->changeName(getName()+"covertInitialYS1v1");
	oplist.push_back(covertInitialYS1v1);
	inPortMap  (covertInitialYS1v1, "X", use("partialConverted2PositiveY1mY0v1"));
	inPortMap  (covertInitialYS1v1, "Y", use("zeroInput4ConvertionYv1"));
	inPortMap  (covertInitialYS1v1, "Cin", use("signY1mY0v1") );
	outPortMap (covertInitialYS1v1, "R","converted2PositiveY1mY0v1");
	vhdl << instance(covertInitialYS1v1, "covertInitialYS1v1");
	
	syncCycleFromSignal("converted2PositiveY1mY0v1");
	
	setCycleFromSignal("zeroInput4ConvertionYv1");
	
	covertInitialYS2v1 = new IntAdder(target,inputWidth-1);
	covertInitialYS2v1->changeName(getName()+"covertInitialYS2v1");
	oplist.push_back(covertInitialYS2v1);
	inPortMap  (covertInitialYS2v1, "X", use("partialConverted2PositiveY3mY2v1"));
	inPortMap  (covertInitialYS2v1, "Y", use("zeroInput4ConvertionYv1"));
	inPortMap  (covertInitialYS2v1, "Cin", use("signY3mY2v1") );
	outPortMap (covertInitialYS2v1, "R","converted2PositiveY3mY2v1");
	vhdl << instance(covertInitialYS2v1, "covertInitialYS2v1");
	
	syncCycleFromSignal("converted2PositiveY3mY2v1");
	
	multiplierYSegv1 = new IntMultiplier(target, inputWidth-1, inputWidth-1);
	multiplierYSegv1->changeName(getName()+"multiplierYSegv1");
	oplist.push_back(multiplierYSegv1);
	inPortMap  (multiplierYSegv1, "X", use("converted2PositiveY1mY0v1"));
	inPortMap  (multiplierYSegv1, "Y", use("converted2PositiveY3mY2v1"));
	outPortMap (multiplierYSegv1, "R","partialComputedProductSY1v1");
	vhdl << instance(multiplierYSegv1, "multiplierYSegv1");
	
	syncCycleFromSignal("partialComputedProductSY1v1");
	
	inputWidthSegments=2*(inputWidth-1)+1;
	vhdl<<tab<<declare("partialComputedProductSY1v1Bit",inputWidthSegments)<<"<= '0' & "<<use("partialComputedProductSY1v1")<<";"<<endl;
	vhdl<<tab<<declare("sign2VectorSegmentYv1",inputWidthSegments)<<"<= (others => "<<use("signSegmentYv1")<<");"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSY1v1",inputWidthSegments)<<"<="<<use("partialComputedProductSY1v1Bit")<<" xor "<<use("sign2VectorSegmentYv1")<<";"<<endl;
	
	//~ if(target->frequencyMHz()>=250)
		//~ nextCycle();
		
	vhdl<<tab<<declare("zeroInput4FinalConvertionYv1",inputWidthSegments)<<"<= (others => '0' );"<<endl;
	vhdl<<tab<<declare("tempsignSegmentYv1",1)<<"<="<<use("signSegmentYv1")<<";"<<endl;
	
	covertFinalYv1 = new IntAdder(target,inputWidthSegments);
	covertFinalYv1->changeName(getName()+"covertFinalYv1");
	oplist.push_back(covertFinalYv1);
	inPortMap  (covertFinalYv1, "X", use("partialConvertedProductSY1v1"));
	inPortMap  (covertFinalYv1, "Y", use("zeroInput4FinalConvertionYv1"));
	inPortMap  (covertFinalYv1, "Cin", use("tempsignSegmentYv1") );
	outPortMap (covertFinalYv1, "R","convertedSegmentYv1");
	vhdl << instance(covertFinalYv1, "covertFinalYv1");
	
	syncCycleFromSignal("convertedSegmentYv1");
	
			
		//(z1-z0)*(z3-z2)
	setCycleFromSignal("segmentZ1mZ0");
	
	vhdl<<tab<<declare("signZ1mZ0v1",1)<<"<= "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorZ1mZ0v1",inputWidth-1)<<"<= ( others=> "<<use("signZ1mZ0v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveZ1mZ0v1",inputWidth-1)<<"<="<<use("sign2vectorZ1mZ0v1")<<" xor "<<use("segmentZ1mZ0")<<range(inputWidth-2,0)<<";"<<endl;
	vhdl<<tab<<declare("zeroInput4ConvertionZv1",inputWidth-1)<<"<= (others => '0' );"<<endl;
	
	vhdl<<tab<<declare("signZ3mZ2v1",1)<<"<= "<<use("segmentZ3mZ2")<<of(inputWidth-1)<<";"<<endl;
	vhdl<<tab<<declare("sign2vectorZ3mZ2v1",inputWidth-1)<<"<= ( others=> "<<use("signZ3mZ2v1")<<");"<<endl;
	vhdl<<tab<<declare("partialConverted2PositiveZ3mZ2v1",inputWidth-1)<<"<="<<use("sign2vectorZ3mZ2v1")<<" xor "<<use("segmentZ3mZ2")<<range(inputWidth-2,0)<<";"<<endl;
	
	
	vhdl<<tab<<declare("signSegmentZv1",1)<<" <= "<<use("signZ1mZ0v1")<<" xor "<<use("signZ3mZ2v1")<<";"<<endl;
	
	
	covertInitialZS1v1 = new IntAdder(target,inputWidth-1);
	covertInitialZS1v1->changeName(getName()+"covertInitialZS1v1");
	oplist.push_back(covertInitialZS1v1);
	inPortMap  (covertInitialZS1v1, "X", use("partialConverted2PositiveZ1mZ0v1"));
	inPortMap  (covertInitialZS1v1, "Y", use("zeroInput4ConvertionZv1"));
	inPortMap  (covertInitialZS1v1, "Cin", use("signZ1mZ0v1") );
	outPortMap (covertInitialZS1v1, "R","converted2PositiveZ1mZ0v1");
	vhdl << instance(covertInitialZS1v1, "covertInitialZS1v1");
	
	syncCycleFromSignal("converted2PositiveZ1mZ0v1");
	
	setCycleFromSignal("zeroInput4ConvertionZv1");
	
	covertInitialZS2v1 = new IntAdder(target,inputWidth-1);
	covertInitialZS2v1->changeName(getName()+"covertInitialZS2v1");
	oplist.push_back(covertInitialZS2v1);
	inPortMap  (covertInitialZS2v1, "X", use("partialConverted2PositiveZ3mZ2v1"));
	inPortMap  (covertInitialZS2v1, "Y", use("zeroInput4ConvertionZv1"));
	inPortMap  (covertInitialZS2v1, "Cin", use("signZ3mZ2v1") );
	outPortMap (covertInitialZS2v1, "R","converted2PositiveZ3mZ2v1");
	vhdl << instance(covertInitialZS2v1, "covertInitialZS2v1");
	
	syncCycleFromSignal("converted2PositiveZ3mZ2v1");
	
	multiplierZSegv1 = new IntMultiplier(target, inputWidth-1, inputWidth-1);
	multiplierZSegv1->changeName(getName()+"multiplierZSegv1");
	oplist.push_back(multiplierZSegv1);
	inPortMap  (multiplierZSegv1, "X", use("converted2PositiveZ1mZ0v1"));
	inPortMap  (multiplierZSegv1, "Y", use("converted2PositiveZ3mZ2v1"));
	outPortMap (multiplierZSegv1, "R","partialComputedProductSZ1v1");
	vhdl << instance(multiplierZSegv1, "multiplierZSegv1");
	
	syncCycleFromSignal("partialComputedProductSZ1v1");
	
	inputWidthSegments=2*(inputWidth-1)+1;
	vhdl<<tab<<declare("partialComputedProductSZ1v1Bit",inputWidthSegments)<<"<= '0' & "<<use("partialComputedProductSZ1v1")<<";"<<endl;
	vhdl<<tab<<declare("sign2VectorSegmentZv1",inputWidthSegments)<<"<= (others => "<<use("signSegmentZv1")<<");"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZ1v1",inputWidthSegments)<<"<="<<use("partialComputedProductSZ1v1Bit")<<" xor "<<use("sign2VectorSegmentZv1")<<";"<<endl;
	
	//~ if(target->frequencyMHz()>=250)
		//~ nextCycle();
		
	vhdl<<tab<<declare("zeroInput4FinalConvertionZv1",inputWidthSegments)<<"<= (others => '0' );"<<endl;
	vhdl<<tab<<declare("tempsignSegmentZv1",1)<<"<="<<use("signSegmentZv1")<<";"<<endl;
	
	covertFinalZv1 = new IntAdder(target,inputWidthSegments);
	covertFinalZv1->changeName(getName()+"covertFinalZv1");
	oplist.push_back(covertFinalZv1);
	inPortMap  (covertFinalZv1, "X", use("partialConvertedProductSZ1v1"));
	inPortMap  (covertFinalZv1, "Y", use("zeroInput4FinalConvertionZv1"));
	inPortMap  (covertFinalZv1, "Cin", use("tempsignSegmentZv1") );
	outPortMap (covertFinalZv1, "R","convertedSegmentZv1");
	vhdl << instance(covertFinalZv1, "covertFinalZv1");
	
	syncCycleFromSignal("convertedSegmentZv1");
	
		//ading the 3 results for segments
	
	vhdl<<tab<<declare("convertedSegmentYv1temp",inputWidthSegments)<<"<="<<use("convertedSegmentYv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentZv1temp",inputWidthSegments)<<"<="<<use("convertedSegmentZv1")<<";"<<endl;
	
	adderPartialResult4Var1 = new IntAdder(target,inputWidthSegments);
	adderPartialResult4Var1->changeName(getName()+"adderPartialResult4Var1");
	oplist.push_back(adderPartialResult4Var1);
	inPortMap  (adderPartialResult4Var1, "X", use("convertedSegmentYv1temp"));
	inPortMap  (adderPartialResult4Var1, "Y", use("convertedSegmentZv1temp"));
	inPortMapCst(adderPartialResult4Var1, "Cin", "'0'" );
	outPortMap (adderPartialResult4Var1, "R","partialResult4Var1");
	vhdl << instance(adderPartialResult4Var1, "adderPartialResult4Var1");
	
	syncCycleFromSignal("partialResult4Var1");
	
	vhdl<<tab<<declare("convertedSegmentXv1temp",inputWidthSegments)<<"<="<<use("convertedSegmentXv1")<<";"<<endl;
	
	adderResult4Var1 = new IntAdder(target,inputWidthSegments);
	adderResult4Var1->changeName(getName()+"adderResult4Var1");
	oplist.push_back(adderResult4Var1);
	inPortMap  (adderResult4Var1, "X", use("partialResult4Var1"));
	inPortMap  (adderResult4Var1, "Y", use("convertedSegmentXv1temp"));
	inPortMapCst(adderResult4Var1, "Cin", "'0'" );
	outPortMap (adderResult4Var1, "R","result4Var1");
	vhdl << instance(adderResult4Var1, "adderResult4Var1");
	
	syncCycleFromSignal("result4Var1");
	
	//int signofLSBI=LSBI>=0?(+1):(-1);
	int signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2+1<<endl;
	
	convert2FPv1 = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2,1,wE,wF);
	convert2FPv1->changeName(getName()+"convert2FPv1");
	oplist.push_back(convert2FPv1);
	inPortMap  (convert2FPv1, "I", use("result4Var1"));
	outPortMap (convert2FPv1, "O","Var1");
	vhdl << instance(convert2FPv1, "convert2FPv1");
	
	syncCycleFromSignal("Var1");
	
	
	//computing the var2 sqrt( (x3-x2)^2 +  (y3-y2)^2 +(z3-z2)^2 )
	
		//computing (x3-x2)^2
	
	setCycleFromSignal("converted2PositiveX3mX2v1");
	
	vhdl<<tab<<declare("segmentX3mX2ns",inputWidth-1)<<"<="<<use("converted2PositiveX3mX2v1")<<";"<<endl;
	
	
	sqrXv2 = new IntSquarer(target, inputWidth-1);
	sqrXv2->changeName(getName()+"sqrXv2");
	oplist.push_back(sqrXv2);
	inPortMap  (sqrXv2, "X", use("segmentX3mX2ns"));
	outPortMap (sqrXv2, "R","sqrX3mX2ns");
	vhdl << instance(sqrXv2, "sqrXv2");
	
	syncCycleFromSignal("sqrX3mX2ns");
	
	
		//computing (y3-y2)^2
	
	setCycleFromSignal("converted2PositiveY3mY2v1");
	
	vhdl<<tab<<declare("segmentY3mY2ns",inputWidth-1)<<"<="<<use("converted2PositiveY3mY2v1")<<";"<<endl;
	
	
	sqrYv2 = new IntSquarer(target, inputWidth-1);
	sqrYv2->changeName(getName()+"sqrYv2");
	oplist.push_back(sqrYv2);
	inPortMap  (sqrYv2, "X", use("segmentY3mY2ns"));
	outPortMap (sqrYv2, "R","sqrY3mY2ns");
	vhdl << instance(sqrYv2, "sqrYv2");
	
	syncCycleFromSignal("sqrY3mY2ns");
	
		//computing (z3-z2)^2
	
	setCycleFromSignal("converted2PositiveZ3mZ2v1");
	
	vhdl<<tab<<declare("segmentZ3mZ2ns",inputWidth-1)<<"<="<<use("converted2PositiveZ3mZ2v1")<<";"<<endl;
	
	
	sqrZv2 = new IntSquarer(target, inputWidth-1);
	sqrZv2->changeName(getName()+"sqrZv2");
	oplist.push_back(sqrZv2);
	inPortMap  (sqrZv2, "X", use("segmentZ3mZ2ns"));
	outPortMap (sqrZv2, "R","sqrZ3mZ2ns");
	vhdl << instance(sqrZv2, "sqrZv2");
	
	syncCycleFromSignal("sqrZ3mZ2ns");
	
	
	//ading the 3 results for segments to be feed to sqrt
	
		
	adderPartialResult4SQRTv2 = new IntAdder(target,(inputWidth-1)*2);
	adderPartialResult4SQRTv2->changeName(getName()+"adderPartialResult4SQRTv2");
	oplist.push_back(adderPartialResult4SQRTv2);
	inPortMap  (adderPartialResult4SQRTv2, "X", use("sqrY3mY2ns"));
	inPortMap  (adderPartialResult4SQRTv2, "Y", use("sqrZ3mZ2ns"));
	inPortMapCst(adderPartialResult4SQRTv2, "Cin", "'0'" );
	outPortMap (adderPartialResult4SQRTv2, "R","partialResult4SQRTv2");
	vhdl << instance(adderPartialResult4SQRTv2, "adderPartialResult4SQRTv2");
	
	syncCycleFromSignal("partialResult4SQRTv2");
	
	vhdl<<tab<<declare("sqrX3mX2nstemp",(inputWidth-1)*2)<<"<="<<use("sqrX3mX2ns")<<";"<<endl;
	
	adderResult4SQRTv2 = new IntAdder(target,(inputWidth-1)*2);
	adderResult4SQRTv2->changeName(getName()+"adderResult4SQRTv2");
	oplist.push_back(adderResult4SQRTv2);
	inPortMap  (adderResult4SQRTv2, "X", use("sqrX3mX2nstemp"));
	inPortMap  (adderResult4SQRTv2, "Y", use("partialResult4SQRTv2"));
	inPortMapCst(adderResult4SQRTv2, "Cin", "'0'" );
	outPortMap (adderResult4SQRTv2, "R","result4SQRTv2");
	vhdl << instance(adderResult4SQRTv2, "adderResult4SQRTv2");
	
	syncCycleFromSignal("result4SQRTv2");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2+1<<endl;
	
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
	
	
	
	vhdl<<tab<<"O<="<<use("Var1")<<range(outputWidth-1,0)<< "or "<<use("Var2")<<range(outputWidth-1,0)<<";"<<endl;

	}

CoilInductance::~CoilInductance() {
}

void CoilInductance::emulate(TestCase * tc)
{
}

void CoilInductance::buildStandardTestCases(TestCaseList* tcl){
	
}

