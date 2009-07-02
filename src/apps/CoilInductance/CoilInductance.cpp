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


CoilInductance::CoilInductance(Target* target, int LSBI, int MSBI, int MaxMSBO,int LSBO, int MSBO, char *filepath) :
	Operator(target), MSBI(MSBI), LSBI(LSBI), MaxMSBO(MaxMSBO),LSBO(LSBO), MSBO(MSBO) ,filepath(filepath){
	
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
	inputWidthSegments=2*(inputWidth-1)+1;
	
	addOutput("O",outputWidth);
	
	//Counters for addressing the memories and for frequency division
	
	vhdl<<endl;
	vhdl<<tab<<"process (clk)"<<endl<<tab<<"variable count:std_logic_vector"<<range(integratorWidth,0)<<":=(others=>'0');"<<endl<<tab<<"begin"<<endl<<tab<<tab<<"if clk'event and clk = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if count < "<<pow(2,integratorWidth)<<" then "<<endl<<tab<<tab<<tab<<tab<<"count:=count+'1';"<<endl<<tab<<tab<<tab<<tab<<declare("out_clk1",1)<<"<= '0';"<<endl<<tab<<tab<<tab<<tab<<declare("out_rst",1)<<"<= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<"elsif count = "<<pow(2,integratorWidth)<<" then "<<endl<<tab<<tab<<tab<<tab<<"count:=count+'1';"<<endl<<tab<<tab<<tab<<tab<<use("out_clk1")<<"<= '1'; "<<endl<<tab<<tab<<tab<<tab<<use("out_rst")<<"<= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<"else "<<endl<<tab<<tab<<tab<<tab<<"count:= CONV_STD_LOGIC_VECTOR(0,"<<integratorWidth<<");"<<endl<<tab<<tab<<tab<<tab<<use("out_clk1")<<"<= '0'; "<<endl<<tab<<tab<<tab<<tab<<use("out_rst")<<"<= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	
	vhdl<<tab<<declare("out_clk2",1)<<"<="<<"clk;"<<endl;
	//vhdl<<tab<<declare("signal_t",integratorWidth)<<" <= "<<" count"<<range(integratorWidth-1,0)<<";"<<endl;
	vhdl<<tab<<declare("signal_tp",integratorWidth+1)<<" <= '0' & "<<" count"<<range(integratorWidth-1,0)<<";"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	//Memories instantiation
	
	//Computing the segments x1-x0 x3-x2 y1-y0 y3-y2 z1-z0 z3-z2
	
	//syncCycleFromSignal("????"); sincronization with memories
	
	addInput("signal_x0",inputWidth);
	addInput("signal_x1",inputWidth);
	addInput("signal_x2",inputWidth);
	addInput("signal_x3",inputWidth);
	
	addInput("signal_y0",inputWidth);
	addInput("signal_y1",inputWidth);
	addInput("signal_y2",inputWidth);
	addInput("signal_y3",inputWidth);
	
	addInput("signal_z0",inputWidth);
	addInput("signal_z1",inputWidth);
	addInput("signal_z2",inputWidth);
	addInput("signal_z3",inputWidth);
	
	/*
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
	*/
	
		
	
	//performing x1-x0
	
	setCycleFromSignal("signal_x0");
	
	vhdl<<tab<<declare("segmentX1mX0",inputWidth)<<" <= "<<use("signal_x1")<<" - "<<use("signal_x0")<<";"<<endl;
	
	//performing x0-x1
	
	setCycleFromSignal("signal_x1");
	
	vhdl<<tab<<declare("segmentX0mX1",inputWidth)<<" <= "<<use("signal_x0")<<" - "<<use("signal_x1")<<";"<<endl;
	
	//performing x0-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("segmentX0mX2",inputWidth)<<" <= "<<use("signal_x0")<<" - "<<use("signal_x2")<<";"<<endl;
	
	
	//performing x3-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("segmentX3mX2",inputWidth)<<" <= "<<use("signal_x3")<<" - "<<use("signal_x2")<<";"<<endl;
	
	
	//performing x3-x0
	
	setCycleFromSignal("signal_x0");
	
	vhdl<<tab<<declare("segmentX3mX0",inputWidth)<<" <= "<<use("signal_x3")<<" - "<<use("signal_x0")<<";"<<endl;
	
	
	//performing y1-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("segmentY1mY0",inputWidth)<<" <= "<<use("signal_y1")<<" - "<<use("signal_y0")<<";"<<endl;
	
	//performing y0-y1
	
	setCycleFromSignal("signal_y1");
	
	vhdl<<tab<<declare("segmentY0mY1",inputWidth)<<" <= "<<use("signal_y0")<<" - "<<use("signal_y1")<<";"<<endl;
	
	//performing y0-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("segmentY0mY2",inputWidth)<<" <= "<<use("signal_y0")<<" - "<<use("signal_y2")<<";"<<endl;
	
	
	
	//performing y3-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("segmentY3mY2",inputWidth)<<" <= "<<use("signal_y3")<<" - "<<use("signal_y2")<<";"<<endl;
	
	
	//performing y3-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("segmentY3mY0",inputWidth)<<" <= "<<use("signal_y3")<<" - "<<use("signal_y0")<<";"<<endl;
	
	
	//performing z1-z0
	
	setCycleFromSignal("signal_z0");

	vhdl<<tab<<declare("segmentZ1mZ0",inputWidth)<<" <= "<<use("signal_z1")<<" - "<<use("signal_z0")<<";"<<endl;

		
	//performing z0-z1
	
	setCycleFromSignal("signal_z1");
	
	vhdl<<tab<<declare("segmentZ0mZ1",inputWidth)<<" <= "<<use("signal_z0")<<" - "<<use("signal_z1")<<";"<<endl;
		
	
	//performing z0-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("segmentZ0mZ2",inputWidth)<<" <= "<<use("signal_z0")<<" - "<<use("signal_z2")<<";"<<endl;
	
	
	//performing z3-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("segmentZ3mZ2",inputWidth)<<" <= "<<use("signal_z3")<<" - "<<use("signal_z2")<<";"<<endl;
	
	
	//performing z3-z0
	
	setCycleFromSignal("signal_z0");
	
	vhdl<<tab<<declare("segmentZ3mZ0",inputWidth)<<" <= "<<use("signal_z3")<<" - "<<use("signal_z0")<<";"<<endl;
	
	nextCycle();
	
	//referenceCycle1 = getCurrentCycle();
	
	//computing the var1 (x1-x0)*(x3-x2)+(y1-y0)*(y3-y2)+(z1-z0)*(z3-z2)
	
	
		//(x1-x0)*(x3-x2)
	//setCycleFromSignal("signX1mX0v1");
	//setCycle(referenceCycle1);
	
	
	vhdl<<tab<<declare("convertedSegmentXv1",inputWidthSegments)<<" <= "<<use("segmentX1mX0")<<" * "<<use("segmentX3mX2")<<";"<<endl;
	
			
	//(y1-y0)*(y3-y2)
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle1);
	
	
	vhdl<<tab<<declare("convertedSegmentYv1",inputWidthSegments)<<" <= "<<use("segmentY1mY0")<<" * "<<use("segmentY3mY2")<<";"<<endl;
	
	
			
		//(z1-z0)*(z3-z2)
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle1);
	
	
	vhdl<<tab<<declare("convertedSegmentZv1",inputWidthSegments)<<" <= "<<use("segmentZ1mZ0")<<" * "<<use("segmentZ3mZ2")<<";"<<endl;
	
	
		
	nextCycle();
	
	
		//ading the 3 results for segments
	
	vhdl<<tab<<declare("convertedSegmentYv1temp",inputWidthSegments)<<" <= "<<use("convertedSegmentYv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentZv1temp",inputWidthSegments)<<" <= "<<use("convertedSegmentZv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentXv1temp",inputWidthSegments)<<" <= "<<use("convertedSegmentXv1")<<";"<<endl;
	
	adder4var = new IntNAdder(target,inputWidthSegments,3);
	adder4var->changeName(getName()+"adder4var");	//aici
	oplist.push_back(adder4var);
	inPortMap  (adder4var, "X0", use("convertedSegmentYv1temp"));
	inPortMap  (adder4var, "X1", use("convertedSegmentZv1temp"));
	inPortMap  (adder4var, "X2", use("convertedSegmentXv1temp") );
	inPortMapCst(adder4var,"Cin","'0'");
	outPortMap (adder4var, "R","result4Var1");
	vhdl << instance(adder4var, "adder4var1");
	
	syncCycleFromSignal("result4Var1");
	
	//int signofLSBI=LSBI>=0?(+1):(-1);
	int signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2<<endl;
	
	convert2FP = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2,1,wE,wF);
	convert2FP->changeName(getName()+"convert2FPv");	//aici
	oplist.push_back(convert2FP);
	inPortMap  (convert2FP, "I", use("result4Var1"));
	outPortMap (convert2FP, "O","Var1");
	vhdl << instance(convert2FP, "convert2FPv1");
	
	syncCycleFromSignal("Var1");
	
		
	//computing the var2 sqrt( (x3-x2)^2 +  (y3-y2)^2 +(z3-z2)^2 )
	
		//computing (x3-x2)^2
	
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("sqrX3mX2s",inputWidthSegments)<<" <= "<<use("segmentX3mX2")<<" * "<<use("segmentX3mX2")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrX3mX2ns",(inputWidth-1)*2)<<" <= "<<use("sqrX3mX2s")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
		
	
		//computing (y3-y2)^2
	
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("sqrY3mY2s",inputWidthSegments)<<" <= "<<use("segmentY3mY2")<<" * "<<use("segmentY3mY2")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrY3mY2ns",(inputWidth-1)*2)<<" <= "<<use("sqrY3mY2s")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	
		//computing (z3-z2)^2
	
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("sqrZ3mZ2s",inputWidthSegments)<<" <= "<<use("segmentZ3mZ2")<<" * "<<use("segmentZ3mZ2")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrZ3mZ2ns",(inputWidth-1)*2)<<" <= "<<use("sqrZ3mZ2s")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrX3mX2nstemp",(inputWidth-1)*2)<<" <= "<<use("sqrX3mX2ns")<<";"<<endl;
	vhdl<<tab<<declare("sqrY3mY2nstemp",(inputWidth-1)*2)<<" <= "<<use("sqrY3mY2ns")<<";"<<endl;
	vhdl<<tab<<declare("sqrZ3mZ2nstemp",(inputWidth-1)*2)<<" <= "<<use("sqrZ3mZ2ns")<<";"<<endl;
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	
	adder4SQRTv = new IntNAdder(target,(inputWidth-1)*2,3);
	adder4SQRTv->changeName(getName()+"adder4SQRTv");	//aici
	oplist.push_back(adder4SQRTv);
	inPortMap  (adder4SQRTv, "X0", use("sqrY3mY2nstemp"));
	inPortMap  (adder4SQRTv, "X1", use("sqrZ3mZ2nstemp"));
	inPortMap  (adder4SQRTv, "X2", use("sqrX3mX2nstemp") );
	inPortMapCst(adder4SQRTv,"Cin","'0'");
	outPortMap (adder4SQRTv, "R","result4SQRTv2");
	vhdl << instance(adder4SQRTv, "adder4SQRTv2");
	
	syncCycleFromSignal("result4SQRTv2");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	convert2FP4sqrtv = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2-1,0,wE,wF);
	convert2FP4sqrtv->changeName(getName()+"convert2FP4sqrtv");	//aici
	oplist.push_back(convert2FP4sqrtv);
	inPortMap  (convert2FP4sqrtv, "I", use("result4SQRTv2"));
	outPortMap (convert2FP4sqrtv, "O","fpSQRVar2");
	vhdl << instance(convert2FP4sqrtv, "convert2FP4sqrtv2");
	
	syncCycleFromSignal("fpSQRVar2");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar2temp",wE+wF+3)<<" <= "<<use("fpSQRVar2")<<";"<<endl;
	
	sqrt4var = new  FPSqrt(target, wE, wF, 1, 0);
	sqrt4var->changeName(getName()+"sqrt4var");	//aici
	oplist.push_back(sqrt4var);
	inPortMap  (sqrt4var, "X", use("fpSQRVar2temp"));
	outPortMap (sqrt4var, "R","Var2");
	vhdl << instance(sqrt4var, "sqrt4var2");
	
	syncCycleFromSignal("Var2");
	
	
	//computing the var3 sqrt(((x3-x0)+(x0-x1)*t)^2+((y3-y0)+(y0-y1)*t)^2+((z3-z0)+(z0-z1)*t)^2)
	
	
	
	
	
	
		//(x3-x0)+(x0-x1)*t
	setCycleFromSignal("convertedSegmentXv1");
	
	
	
	vhdl<<tab<<declare("partialComputedProductSXMv3", inputWidth+integratorWidth)<<" <= "<<use("segmentX0mX1")<<" * "<<use("signal_tp")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSXv3",inputWidth)<<" <= "<<use("partialComputedProductSXMv3")<<range(inputWidth-1+integratorWidth,integratorWidth)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sumXPartv3",inputWidth)<<" <= "<<use("partialConvertedProductSXv3")<<" + "<<use("segmentX3mX0")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrXsv3",inputWidthSegments)<<" <= "<<use("sumXPartv3")<<" * "<<use("sumXPartv3")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrXnsv3",(inputWidth-1)*2)<<" <= "<<use("sqrXsv3")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	
	
		//(y3-y0)+(y0-y1)*t
	setCycleFromSignal("convertedSegmentXv1");
	
	
	
	vhdl<<tab<<declare("partialComputedProductSYMv3", inputWidth+integratorWidth)<<" <= "<<use("segmentY0mY1")<<" * "<<use("signal_tp")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSYv3",inputWidth)<<" <= "<<use("partialComputedProductSYMv3")<<range(inputWidth-1+integratorWidth,integratorWidth)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sumYPartv3",inputWidth)<<" <= "<<use("partialConvertedProductSYv3")<<" + "<<use("segmentY3mY0")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrYsv3",inputWidthSegments)<<" <= "<<use("sumYPartv3")<<" * "<<use("sumYPartv3")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrYnsv3",(inputWidth-1)*2)<<" <= "<<use("sqrYsv3")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	
	//(z3-z0)+(z0-z1)*t
	setCycleFromSignal("convertedSegmentXv1");
	
	
	
	vhdl<<tab<<declare("partialComputedProductSZMv3", inputWidth+integratorWidth)<<" <= "<<use("segmentZ0mZ1")<<" * "<<use("signal_tp")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZv3",inputWidth)<<" <= "<<use("partialComputedProductSZMv3")<<range(inputWidth-1+integratorWidth,integratorWidth)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sumZPartv3",inputWidth)<<" <= "<<use("partialConvertedProductSZv3")<<" + "<<use("segmentZ3mZ0")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrZsv3",inputWidthSegments)<<" <= "<<use("sumZPartv3")<<" * "<<use("sumZPartv3")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrZnsv3",(inputWidth-1)*2)<<" <= "<<use("sqrZsv3")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrXnsv3temp",(inputWidth-1)*2)<<" <= "<<use("sqrXnsv3")<<";"<<endl;
	vhdl<<tab<<declare("sqrYnsv3temp",(inputWidth-1)*2)<<" <= "<<use("sqrYnsv3")<<";"<<endl;
	vhdl<<tab<<declare("sqrZnsv3temp",(inputWidth-1)*2)<<" <= "<<use("sqrZnsv3")<<";"<<endl;
	
	
	
	inPortMap  (adder4SQRTv, "X0", use("sqrYnsv3temp"));
	inPortMap  (adder4SQRTv, "X1", use("sqrZnsv3temp"));
	inPortMap  (adder4SQRTv, "X2", use("sqrXnsv3temp") );
	inPortMapCst(adder4SQRTv,"Cin","'0'");
	outPortMap (adder4SQRTv, "R","result4SQRTv3");
	vhdl << instance(adder4SQRTv, "adder4SQRTv3");
	
	syncCycleFromSignal("result4SQRTv3");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	inPortMap  (convert2FP4sqrtv, "I", use("result4SQRTv3"));
	outPortMap (convert2FP4sqrtv, "O","fpSQRVar3");
	vhdl << instance(convert2FP4sqrtv, "convert2FP4sqrtv3");
	
	syncCycleFromSignal("fpSQRVar3");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar3temp",wE+wF+3)<<" <= "<<use("fpSQRVar3")<<";"<<endl;
	
	inPortMap  (sqrt4var, "X", use("fpSQRVar3temp"));
	outPortMap (sqrt4var, "R","Var3");
	vhdl << instance(sqrt4var, "sqrt4var3");
	
	syncCycleFromSignal("Var3");
	
	
	target->setNotPipelined();
	acc4var = new FPAdder(target, wE, wF, wE, wF, wE, wF);
	target->setPipelined();
	acc4var->changeName(getName()+"accumulator4var");	
	oplist.push_back(acc4var);
	inPortMap  (acc4var, "X", use("Var3"));
	inPortMapCst (acc4var, "Y","accVar3");
	outPortMap (acc4var, "R","tempVar3");
	vhdl << instance(acc4var, "accumulator4var3");
	
	syncCycleFromSignal("tempVar3");
	
	vhdl<<tab<<"process(clk)"<<endl<<tab<<"variable temp: std_logic_vector( "<<wE+wF+3-1<<" downto 0):=(others=>'0');"<<endl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk2'event and out_clk2 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if out_rst = '0' then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"temp:="<<use("tempVar3")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"else"<<endl<<tab<<tab<<tab<<tab<<"temp:=(others=>'0');"<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("accVar3",wF+wE+3)<<"<= temp;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
	
	
	//computing the var4 sqrt(((x0-x2)+(x1-x0)*t)^2+((y0-y2)+(y1-y0)*t)^2+((z0-z2)+(z1-z0)*t)^2)
	
	
	
	
	
		//((x0-x2)+(x1-x0)*t)^2
	setCycleFromSignal("convertedSegmentXv1");
	
	vhdl<<tab<<declare("partialComputedProductSXMv4", inputWidth+integratorWidth)<<" <= "<<use("segmentX1mX0")<<" * "<<use("signal_tp")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSXv4",inputWidth)<<" <= "<<use("partialComputedProductSXMv4")<<range(inputWidth-1+integratorWidth,integratorWidth)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sumXv4",inputWidth)<<" <= "<<use("partialConvertedProductSXv4")<<" + "<<use("segmentX0mX2")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrXsv4",inputWidthSegments)<<" <= "<<use("sumXv4")<<" * "<<use("sumXv4")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrXnsv4",(inputWidth-1)*2)<<" <= "<<use("sqrXsv4")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
		
		
	
		//((y0-y2)+(y1-y0)*t)^2
	setCycleFromSignal("convertedSegmentXv1");
	
	vhdl<<tab<<declare("partialComputedProductSYMv4", inputWidth+integratorWidth)<<" <= "<<use("segmentY1mY0")<<" * "<<use("signal_tp")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSYv4",inputWidth)<<" <= "<<use("partialComputedProductSYMv4")<<range(inputWidth-1+integratorWidth,integratorWidth)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sumYv4",inputWidth)<<" <= "<<use("partialConvertedProductSYv4")<<" + "<<use("segmentY0mY2")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrYsv4",inputWidthSegments)<<" <= "<<use("sumYv4")<<" * "<<use("sumYv4")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrYnsv4",(inputWidth-1)*2)<<" <= "<<use("sqrYsv4")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
			//((z0-z2)+(z1-z0)*t)^2
	setCycleFromSignal("convertedSegmentXv1");
	
	vhdl<<tab<<declare("partialComputedProductSZMv4", inputWidth+integratorWidth)<<" <= "<<use("segmentZ1mZ0")<<" * "<<use("signal_tp")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZv4",inputWidth)<<" <= "<<use("partialComputedProductSZMv4")<<range(inputWidth-1+integratorWidth,integratorWidth)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sumZv4",inputWidth)<<" <= "<<use("partialConvertedProductSZv4")<<" + "<<use("segmentZ0mZ2")<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrZsv4",inputWidthSegments)<<" <= "<<use("sumZv4")<<" * "<<use("sumZv4")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrZnsv4",(inputWidth-1)*2)<<" <= "<<use("sqrZsv4")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrXnsv4temp",(inputWidth-1)*2)<<" <= "<<use("sqrXnsv4")<<";"<<endl;
	vhdl<<tab<<declare("sqrYnsv4temp",(inputWidth-1)*2)<<" <= "<<use("sqrYnsv4")<<";"<<endl;
	vhdl<<tab<<declare("sqrZnsv4temp",(inputWidth-1)*2)<<" <= "<<use("sqrZnsv4")<<";"<<endl;
	
	
	inPortMap  (adder4SQRTv, "X0", use("sqrYnsv4temp"));
	inPortMap  (adder4SQRTv, "X1", use("sqrZnsv4temp"));
	inPortMap  (adder4SQRTv, "X2", use("sqrXnsv4temp") );
	inPortMapCst(adder4SQRTv,"Cin","'0'");
	outPortMap (adder4SQRTv, "R","result4SQRTv4");
	vhdl << instance(adder4SQRTv, "adder4SQRTv4");
	
	syncCycleFromSignal("result4SQRTv4");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	
	inPortMap  (convert2FP4sqrtv, "I", use("result4SQRTv4"));
	outPortMap (convert2FP4sqrtv, "O","fpSQRVar4");
	vhdl << instance(convert2FP4sqrtv, "convert2FP4sqrtv4");
	
	syncCycleFromSignal("fpSQRVar4");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar4temp",wE+wF+3)<<" <= "<<use("fpSQRVar4")<<";"<<endl;
	
	inPortMap  (sqrt4var, "X", use("fpSQRVar4temp"));
	outPortMap (sqrt4var, "R","Var4");
	vhdl << instance(sqrt4var, "sqrt4var4");
	
	syncCycleFromSignal("Var4");
	
	
	inPortMap  (acc4var, "X", use("Var4"));
	inPortMapCst (acc4var, "Y","accVar4");
	outPortMap (acc4var, "R","tempVar4");
	vhdl << instance(acc4var, "accumulator4var4");
	
	syncCycleFromSignal("tempVar4");
	
	nextCycle(); //in order for the whole pipeline to be syncronized with the accumulator -> used as a reference signal for syncronization the accVar4
	
	vhdl<<tab<<"process(clk)"<<endl;
	vhdl<<tab<<"variable temp: std_logic_vector( "<<wE+wF+3-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk2'event and out_clk2 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if out_rst = '0' then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"temp:="<<use("tempVar4")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"else"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"temp:=(others=>'0');"<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("accVar4",wF+wE+3)<<"<= temp;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
		
	
	//computing the var5 ( ((x0-x2)+(x1-x0)*t)*(x3-x2) +((y0-y2)+(y1-y0)*t)*(y3-y2)+((z0-z2)+(z1-z0)*t)*(z3-z2))
	
	
	
	
	
		//((x0-x2)+(x1-x0)*t)*(x3-x2)
	setCycleFromSignal("sqrXsv4");
	
	vhdl<<tab<<declare("convertedProductSXv5",inputWidthSegments)<<" <= "<<use("segmentX3mX2")<<" * "<<use("sumXv4")<<";"<<endl;
	
	
	
	
		//((y0-y2)+(y1-y0)*t)*(y3-y2)
	setCycleFromSignal("sqrXsv4");
	
	vhdl<<tab<<declare("convertedProductSYv5",inputWidthSegments)<<" <= "<<use("segmentY3mY2")<<" * "<<use("sumYv4")<<";"<<endl;
	
	
	
		//((z0-z2)+(z1-z0)*t)*(z3-z2)
	setCycleFromSignal("sqrXsv4");
	
	vhdl<<tab<<declare("convertedProductSZv5",inputWidthSegments)<<" <= "<<use("segmentZ3mZ2")<<" * "<<use("sumZv4")<<";"<<endl;
	
	
	
	
	
	//Adding the 3 results and the caries that are neded to finalize the transformation of the numbers in 2's complement
	
	nextCycle();
	
	vhdl<<tab<<declare("partialConvertedProductSYv5temp",inputWidthSegments)<<" <= "<<use("convertedProductSYv5")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZv5temp",inputWidthSegments)<<" <= "<<use("convertedProductSZv5")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSXv5temp",inputWidthSegments)<<" <= "<<use("convertedProductSXv5")<<";"<<endl;
	
	
	
	
	inPortMap  (adder4var, "X0", use("partialConvertedProductSYv5temp"));
	inPortMap  (adder4var, "X1", use("partialConvertedProductSZv5temp"));
	inPortMap  (adder4var, "X2", use("partialConvertedProductSXv5temp") );
	inPortMapCst(adder4var,"Cin","'0'");
	outPortMap (adder4var, "R","result4Var5");
	vhdl << instance(adder4var, "adder4var5");
	
	syncCycleFromSignal("result4Var5");
	
	//int signofLSBI=LSBI>=0?(+1):(-1);
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2<<endl;
	
	
	nextCycle();
	
	vhdl<<tab<<declare("tempVar5",inputWidthSegments)<<" <= "<<use("result4Var5")<<" + accVar5fix"<<";"<<endl; //the value is hardcoded because otherwise some problems would occur because of the later declaration
	
	vhdl<<tab<<"process(clk)"<<endl;
	vhdl<<tab<<"variable temp: std_logic_vector( "<<inputWidthSegments-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk2'event and out_clk2 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if out_rst = '0' then"<<endl<<tab<<tab<<tab<<tab<<"temp:="<<use("tempVar5")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"else"<<endl<<tab<<tab<<tab<<tab<<"temp:=(others=>'0');"<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("accVar5fix",inputWidthSegments)<<"<= temp;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
		
	inPortMap  (convert2FP, "I", use("accVar5fix"));
	outPortMap (convert2FP, "O","accVar5");
	vhdl << instance(convert2FP, "convert2FPv5");
	
	syncCycleFromSignal("accVar5");
	
	//de aici in jos toate componentele ar trebui sa functioneze pe clockul out_clk1 (este o diviziune cu 9+1 a clk)
	
	setCycleFromSignal("accVar4");
	
	
	vhdl<<tab<<declare("Var1temp",wE+wF+3)<<" <= "<<use("Var1")<<";"<<endl;
	vhdl<<tab<<declare("Var2temp1",wE+wF+3)<<" <= "<<use("Var2")<<";"<<endl;
	
	div4Log =new FPDiv(target, wE, wF);
	div4Log->changeName(getName()+"div4Acc");
	oplist.push_back(div4Log);
	inPortMap  (div4Log, "X", use("Var1temp"));
	inPortMap (div4Log, "Y",use("Var2temp1"));
	outPortMap (div4Log, "R","var1divvar2");
	vhdl << instance(div4Log, "var1divvar24acc");
	
	syncCycleFromSignal("var1divvar2");
	
	
	setCycleFromSignal("accVar4");
		
	
	vhdl<<tab<<declare("Var2temp",wE+wF+3)<<" <= "<<use("Var2")<<";"<<endl;
	vhdl<<tab<<declare("accVar3temp",wE+wF+3)<<" <= "<<use("accVar3")<<";"<<endl;
	
	inPortMap  (acc4var, "X", use("Var2temp"));
	inPortMap (acc4var, "Y",use("accVar3temp"));
	outPortMap (acc4var, "R","var3pvar2");
	vhdl << instance(acc4var, "var3plusvar2");
	
	syncCycleFromSignal("var3pvar2");
	
	nextCycle();
	
	setCycleFromSignal("accVar4");
	
	vhdl<<tab<<declare("accVar5temp",wE+wF+3)<<" <= "<<use("accVar5")<<";"<<endl;
	vhdl<<tab<<declare("Var2temp2",wE+wF+3)<<" <= "<<use("Var2")<<";"<<endl;
	
	inPortMap  (div4Log, "X", use("accVar5temp"));
	inPortMap (div4Log, "Y",use("Var2temp2"));
	outPortMap (div4Log, "R","var5divvar2");
	vhdl << instance(div4Log, "var5divvar24log");
	
	syncCycleFromSignal("var5divvar2");
	
	vhdl<<tab<<declare("var3pvar2temp",wE+wF+3)<<" <= "<<use("var3pvar2")<<";"<<endl;
	vhdl<<tab<<declare("minusvar5divvar2",wE+wF+3)<<" <= "<<use("var5divvar2")<<range(wE+wF+2,wE+wF+1)<<" & "<<"( not ("<<use("var5divvar2")<<of(wE+wF)<<")) & "<<use("var5divvar2")<<range(wE+wF-1,0)<<";"<<endl;	
	
	inPortMap  (acc4var, "X", use("var3pvar2temp"));
	inPortMap (acc4var, "Y",use("minusvar5divvar2"));
	outPortMap (acc4var, "R","numerator4Log");
	vhdl << instance(acc4var, "numerator4LogAdder");
	
	syncCycleFromSignal("numerator4Log");
	
	
	setCycleFromSignal("var5divvar2");
	
	vhdl<<tab<<declare("accVar4temp",wE+wF+3)<<" <= "<<use("accVar4")<<";"<<endl;
	
	
	inPortMap  (acc4var, "X", use("accVar4temp"));
	inPortMap (acc4var, "Y",use("minusvar5divvar2"));
	outPortMap (acc4var, "R","denominator4Log");
	vhdl << instance(acc4var, "denominator4LogAdder");
	
	syncCycleFromSignal("denominator4Log");
	
	
	inPortMap  (div4Log, "X", use("numerator4Log"));
	inPortMap (div4Log, "Y",use("denominator4Log"));
	outPortMap (div4Log, "R","result4Log");
	vhdl << instance(div4Log, "div4log");
	
	syncCycleFromSignal("result4Log");
	
	

	log4Acc = new FPLog(target, wE, wF);
	log4Acc->changeName(getName()+"log4Acc");
	oplist.push_back(log4Acc);
	inPortMap  (log4Acc, "X", use("result4Log"));
	outPortMap (log4Acc, "R","resultLog");
	vhdl << instance(log4Acc, "log4acc");
	
	
	syncCycleFromSignal("resultLog");
	
	vhdl<<tab<<declare("var1divvar2temp",wE+wF+3)<<" <= "<<use("var1divvar2")<<";"<<endl;
	
	target->setNotPipelined();
	mult4Acc = new FPMultiplier(target, wE, wF, wE, wF, wE, wF, 1);
	target->setPipelined();
	mult4Acc->changeName(getName()+"mult4Acc");
	oplist.push_back(mult4Acc);
	inPortMap  (mult4Acc, "X", use("var1divvar2temp"));
	inPortMap (mult4Acc, "Y",use("resultLog"));
	outPortMap (mult4Acc, "R","value4LongAcc");
	vhdl << instance(mult4Acc, "mult4acc");
	
	
	syncCycleFromSignal("value4LongAcc");
	
	nextCycle();
	
	vhdl<<tab<<declare("value4LongAcctemp",wE+wF+3)<<" <= "<<use("value4LongAcc")<<";"<<endl;
	
	finalAcc = new LongAcc(target, wE, wF, MaxMSBO, LSBO, MSBO);
	finalAcc->changeName(getName()+"finalAcc");
	oplist.push_back(finalAcc);
	inPortMap  (finalAcc, "X", use("value4LongAcctemp"));
	outPortMap (finalAcc, "A","finalResult");
	outPortMap (finalAcc, "XOverflow","XOverflow");
	outPortMap (finalAcc, "XUnderflow","XUnderflow");
	outPortMap (finalAcc, "AccOverflow","AccOverflow");
	vhdl << instance(finalAcc, "finalAcc");
		
	syncCycleFromSignal("finalResult");
		
	
	vhdl<<tab<<"O <= "<<use("finalResult")<<range(outputWidth-1,0)<<";"<<endl;
	
	//vhdl<<tab<<"O<="<<use("Var1")<<range(outputWidth-1,0)<<" or "<<use("accVar4")<<range(outputWidth-1,0)<<" or "<<use("numerator4Log")<<range(outputWidth-1,0)<<";"<<endl;

	}

void CoilInductance::outputVHDL(std::ostream& o, std::string name) {
  
	licence(o);
	o << "library ieee; " << endl;
	o << "use ieee.std_logic_1164.all;" << endl;
	o << "use ieee.std_logic_arith.all;" << endl;
	o << "use ieee.std_logic_signed.all;" << endl;
	o << "library work;" << endl;
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLComponentDeclarations();	
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);		
	o<<buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);

}
	
	
CoilInductance::~CoilInductance() {
}

void CoilInductance::emulate(TestCase * tc)
{
}

void CoilInductance::buildStandardTestCases(TestCaseList* tcl){
	
}

