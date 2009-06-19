/*
 * Floating Point Adder for FloPoCo
 *
 * Author : Bogdan Pasca, Florent de Dinechin, Radu Tudoran
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

// TODO move close path prenormalization up to the Swap Difference box
//   if it becomes a part of the critical path
// TODO remove pipeline stage after finalRoundAdd if slack allows

// TODO clean up X propagation to remove warnings

// TODO Single path adder

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>


#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "Fix2FP.hpp"

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


Fix2FP::Fix2FP(Target* target, int MSBI, int LSBI, int wER, int wFR) :
	Operator(target), MSBI(MSBI), LSBI(LSBI), wER(wER), wFR(wFR) {

	ostringstream name;

	long absMSB=MSBI>=0?MSBI:-MSBI;
	long absLSB=LSBI>=0?LSBI:-LSBI;
	name<<"Fix2FP_"<< absMSB <<"_"<<absLSB<<"_"<<wER<<"_"<<wFR; 
	setName(name.str()); 

	setCopyrightString("Bogdan Pasca, Florent de Dinechin, Radu Tudoran (2009)");		

	//parameter set up
	
	wF = wFR;
	wE = wER;
		if(MSBI>0)
			MSB=MSBI+1;
		else
			MSB=MSBI;
	LSB=LSBI;
		
			
		
	/* Set up the IO signals */
		
	addInput ("I", MSB-LSB);
	addFPOutput("O", wE,wF);
	
	/*	VHDL code description	*/
	
	inputWidth=MSB-LSB;
	
	vhdl<<declare("input",inputWidth) << " <= I;"<<endl;
	
	// code for the LZOCShifter part
	
	
	//ostringstream signIndex;
	//ostringstream magnitudeRange;
	//signIndex<<"input("<<MSB-1<<")";
	//magnitudeRange<<"input("<<MSB-1<<" downto "<<LSB<< ")";
	
	
	
	//lzocs =  new LZOCShifter(target,MSB-LSB-1,MSB-LSB-1);
	
	vhdl<<declare("signSignal",1)<<"<="<<use("input")<<of(MSB-1-LSB)<<";"<<endl;
	vhdl<<declare("passedInput",inputWidth)<<"<="<<use("input")<<range(MSB-1 -LSB,0)<<";"<<endl;
	vhdl<<declare("input2LZOC",inputWidth-1)<<"<="<<use("passedInput")<<range(MSB-2 -LSB,0)<<";"<<endl;
	
	//int maximalOutputValue = (MSB-LSB)>(wF+4)?MSB-LSB:(wF+4);
	if((MSB-LSB)>(wF))
	{
		
		//code for the case when a rounding is required
	int maximalOutputValue = (MSB-LSB)>wF+4? MSB-LSB:wF+4;
	
	lzocs		= new LZOCShifterSticky(target,inputWidth-1 , maximalOutputValue, intlog2(inputWidth-1), 0, -1);
	lzocs->changeName(getName()+"_LZCS");
	oplist.push_back(lzocs);
	//inPortMap  (lzocs, "I", magnitudeRange.str());
	inPortMap  (lzocs, "I", use("input2LZOC"));
	//inPortMap	(lzocs,"OZB",signIndex.str());
	inPortMap	(lzocs,"OZb",use("signSignal"));
	outPortMap (lzocs, "Count","temporalExponent");
	outPortMap (lzocs, "O","temporalFraction");
	vhdl << instance(lzocs, "LZOC_component");
	
	sizeExponentValue=lzocs->getCountWidth();
	sizeFractionPlusOne=wF+1;
	
	syncCycleFromSignal("temporalExponent");
	
	
	//Code for creating the exponent
	
	
	vhdl<<declare("MSB2Signal",wE)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<MSB-2<<","<<wE<<");"<<endl;
	vhdl<<declare("zeroPadding4Exponent",wE- intlog2(inputWidth-1))<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE- intlog2(inputWidth-1)<<");"<<endl;
	vhdl<<declare("valueExponent",wE)<<"<= not("<<use("zeroPadding4Exponent")<<" & "<<use("temporalExponent")<<");"<<endl;
	//vhdl<<declare("oneBit",1)<<"<='1';"<<endl;
	
	exponentConvertion = new IntAdder(target,wE);
	exponentConvertion->changeName(getName()+"exponentConvertion");
	oplist.push_back(exponentConvertion);
	inPortMap  (exponentConvertion, "X", use("MSB2Signal"));
	inPortMap  (exponentConvertion, "Y", use("valueExponent"));
	inPortMapCst(exponentConvertion, "Cin", "'1'");
	//inPortMap  (exponentConvertion, "Cin", "oneBit");
	outPortMap (exponentConvertion, "R","partialConvertedExponent");
	vhdl << instance(exponentConvertion, "exponentConvertion");
	
	syncCycleFromSignal("partialConvertedExponent");
	
	vhdl<<declare("biassOfOnes",wE-1)<<"<=CONV_STD_LOGIC_VECTOR("<<pow(2,wE)-1<<","<<wE-1<<");"<<endl;
	vhdl<<declare("biassSignal",wE)<<"<="<<"'0' &"<<use("biassOfOnes")<<";"<<endl;
	vhdl<<declare("biassSignalBit",wE+1)<<"<="<<"'0' &"<<use("biassSignal")<<";"<<endl;
	//vhdl<<declare("zeroBitExponent",1)<<"<='0';"<<endl;
	vhdl<<declare("partialConvertedExponentBit",wE+1)<<"<= '0' & "<<use("partialConvertedExponent")<<";"<<endl;
	vhdl<<declare("sign4OU",1)<<"<="<<use("partialConvertedExponent")<<of(wE-1)<<";"<<endl;
	
	exponentFinal = new IntAdder(target,wE+1);
	exponentFinal->changeName(getName()+"exponentFinal");
	oplist.push_back(exponentFinal);
	inPortMap  (exponentFinal, "X", use("partialConvertedExponentBit"));
	inPortMap  (exponentFinal, "Y", use("biassSignalBit"));
	inPortMapCst(exponentFinal, "Cin", "'0'");
	//inPortMap  (exponentFinal, "Cin", "zeroBitExponent");
	outPortMap (exponentFinal, "R","convertedExponentBit");
	vhdl << instance(exponentFinal, "exponentFinal");
	
	syncCycleFromSignal("convertedExponentBit");
	
	vhdl<<declare("convertedExponent",wE)<<"<= "<<use("convertedExponentBit")<<range(wE-1,0)<<";"<<endl;
	
	
	vhdl<<declare("underflowSignal",1)<<"<= '1' when ("<<use("sign4OU")<<" = '1' and "<<use("convertedExponentBit")<<range(wE,wE-1)<<" = \"01\" ) else '0' ;"<<endl;
	vhdl<<declare("overflowSignal",1)<<"<= '1' when ("<<use("sign4OU")<<" = '0' and "<<use("convertedExponentBit")<<range(wE,wE-1)<<" = \"10\" ) else '0' ;"<<endl;
	
	//code for verifing if the number is zero
	
	setCycleFromSignal("passedInput");
	
	vhdl<<declare("minusOne4ZD",MSB -LSB)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<MSB-LSB<<");"<<endl;
	//vhdl<<declare("carry4ZD",1)<<" <= '0'; "<<endl;
	
	zeroD = new IntAdder(target, MSB-LSB );
	zeroD->changeName(getName()+"zeroD");
	oplist.push_back(zeroD);
	inPortMap  (zeroD, "X", use("passedInput"));
	inPortMap  (zeroD, "Y", use("minusOne4ZD"));
	inPortMapCst(zeroD, "Cin", "'0'");
	//inPortMap  (zeroD, "Cin", "carry4ZD");
	outPortMap (zeroD, "R","zeroDS");
	vhdl << instance(zeroD, "zeroD");
	
	syncCycleFromSignal("zeroDS");
	//selection of mux 2
	vhdl<<declare("zeroInput",1)<<"<= "<<use("zeroDS")<<of(MSB-LSB-1)<<" and not("<<use("signSignal")<<");"<<endl;
	
	
	
	//code for the Convertion of the fraction
	
	//syncCycleFromSignal("tempFractionResult");
	setCycleFromSignal("temporalFraction");
	
	//vhdl<<declare("tfr",sizeFractionPlusOne) <<"<="<<use("temporalFraction")<<range(maximalOutputValue-1,maximalOutputValue-sizeFractionPlusOne)<<";"<<endl;
		
	//vhdl<<declare("semn1",23)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<5<<","<<bbb<<");"<<endl;
	
	vhdl<<declare("sign2vector",maximalOutputValue)<<"<="<<"(others=>"<<use("signSignal")<<");"<<endl;
	vhdl<<declare("tempConvert",maximalOutputValue)<<"<="<<use("sign2vector")<<" xor "<<use("temporalFraction")<<";"<<endl;
	vhdl<<declare("tempConvert0",maximalOutputValue+1)<<"<= '0' & "<<use("tempConvert")<<";"<<endl;
	//vhdl<<declare("tempAddSign",sizeFractionPlusOne)<<"<="<<"CONV_STD_LOGIC_VECTOR("<< use("input")<<of(MSB-1)<<","<<sizeFractionPlusOne<<");"<<endl;
	vhdl<<declare("tempPaddingAddSign",maximalOutputValue)<<"<=(others=>'0');"<<endl;
	vhdl<<declare("tempAddSign",maximalOutputValue+1)<<"<="<<use("tempPaddingAddSign")<<" & "<<use("signSignal")<<";"<<endl;
	
	
		//Integer adder for obtaining the fraction value
		
	//vhdl<<declare("zeroBit",1)<<"<="<<" '0'; "<<endl;
	
	fractionConvert = new IntAdder(target,maximalOutputValue+1);
	fractionConvert->changeName(getName()+"_fractionConvert");
	oplist.push_back(fractionConvert);
	inPortMap  (fractionConvert, "X", use("tempConvert0"));
	inPortMap  (fractionConvert, "Y", use("tempAddSign"));
	inPortMapCst(fractionConvert, "Cin", "'0'");
	//inPortMap  (fractionConvert, "Cin", "zeroBit");
	outPortMap (fractionConvert, "R","tempFractionResult");
	vhdl << instance(fractionConvert, "fractionConverter");
	
	syncCycleFromSignal("tempFractionResult");
	vhdl<<declare("correctingExponent",1)<<"<="<<use("tempFractionResult")<<of(maximalOutputValue)<<";"<<endl;
	vhdl<<declare("fractionConverted",wF)<<"<="<<use("tempFractionResult")<<range(maximalOutputValue-2,maximalOutputValue -wF-1)<<";"<<endl;
	
	//Zero Compare of the bits from the remainder of magnitude
	
	int sizeOfRemainder=maximalOutputValue-sizeFractionPlusOne-1;
	
	setCycleFromSignal("tempFractionResult");
	vhdl<<declare("minusOne",sizeOfRemainder)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<sizeOfRemainder<<");"<<endl;
	vhdl<<declare("fractionRemainder",sizeOfRemainder)<<"<="<<use("tempFractionResult")<<range(sizeOfRemainder-1,0)<<";"<<endl;
	//vhdl<<declare("zeroBit2",1)<<"<="<<" '0'; "<<endl;
	
	//cout<<"aici "<<sizeOfRemainder<<" tot aici";
	
	oneSubstracter = new IntAdder(target,sizeOfRemainder);
	oneSubstracter->changeName(getName()+"_oneSubstracter");
	oplist.push_back(oneSubstracter);
	inPortMap  (oneSubstracter, "X", use("fractionRemainder"));
	inPortMap  (oneSubstracter, "Y", use("minusOne"));
	inPortMapCst(oneSubstracter, "Cin", "'0'");
	//inPortMap  (oneSubstracter, "Cin", "zeroBit2");
	outPortMap (oneSubstracter, "R","zeroFractionResult");
	vhdl << instance(oneSubstracter, "oneSubstracter");
	
	syncCycleFromSignal("zeroFractionResult");
	//selection of mux 2
	vhdl<<declare("zeroRemainder",1)<<"<= not( "<<"not ("<<use("tempFractionResult")<<of(sizeOfRemainder-1)<<") and "<<use("zeroFractionResult")<<of(sizeOfRemainder-1)<<");"<<endl;
	
	
	// signals for Muxes
	//selection of mux 1
	//de luat semnalele dupa conversie
	vhdl<<declare("firstBitofRest",1)<<"<="<<use("tempFractionResult")<<of(maximalOutputValue-wF-2)<<";"<<endl;
	
	//selection of mux 3
	vhdl<<declare("lastBitOfFraction")<<"<="<<use("tempFractionResult")<<of(maximalOutputValue-wF-1)<<";"<<endl;
	
	//vhdl<<"with "<<use("lastBitOfFraction")<<"select "<<endl<<declare("outputOfMux3",1)<<" <= '0' when '0', "<<endl<<"'1' when others;";
	vhdl<<declare("outputOfMux3",1)<<"<="<<use("lastBitOfFraction")<<";"<<endl;
	vhdl<<"with "<<use("zeroRemainder")<<" select "<<endl<<declare("outputOfMux2",1)<<" <= "<< use("outputOfMux3")<<" when '0', "<<endl<<"\t'1' when others;"<<endl;
	vhdl<<"with "<<use("firstBitofRest")<<" select "<<endl<<declare("outputOfMux1",1)<<" <= "<< use("outputOfMux2")<<" when '1', "<<endl<<"\t'0' when others;"<<endl;
	//vhdl<<declare("zeroInput4Rounding",wF+wE+1)<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wF+wE+1<<");"<<endl;
	vhdl<<declare("possibleCorrector4Rounding",wF+wE+1)<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE<<") & "<<use("correctingExponent")<<" & "<<"CONV_STD_LOGIC_VECTOR(0,"<<wF<<");"<<endl;
	vhdl<<declare("concatenationForRounding",wE+wF+1)<<"<= '0' &"<<use("convertedExponent")<<" & "<<use("fractionConverted")<<";"<<endl;
	
	
	roundingApproximator = new IntAdder(target,wF+wE+1);
	roundingApproximator->changeName(getName()+"roundingApproximator");
	oplist.push_back(roundingApproximator);
	inPortMap  (roundingApproximator, "X", use("concatenationForRounding"));
	inPortMap  (roundingApproximator, "Y", use("possibleCorrector4Rounding"));
	inPortMap  (roundingApproximator, "Cin", use("outputOfMux1"));
	outPortMap (roundingApproximator, "R","roundedResult");
	vhdl << instance(roundingApproximator, "roundingApproximator");
	
	syncCycleFromSignal("roundedResult");
	
	
	vhdl<<declare("convertedExponentAfterRounding",wE)<<"<="<<use("roundedResult")<<range(wE+wF-1,wF)<<";"<<endl;
	vhdl<<declare("convertedFractionAfterRounding",wF)<<"<="<<use("roundedResult")<<range(wF-1,0)<<";"<<endl;
	
	vhdl<<declare("MSBSelection",1)<<"<= "<<use("overflowSignal")<<" or "<<use("roundedResult")<<of(wF+wE)<<";"<<endl;
	vhdl<<declare("LSBSelection",1)<<"<= not("<<use("underflowSignal")<<" and not ( "<<use("zeroInput")<<" )) ;"<<endl;
	vhdl<<declare("Selection",2)<<"<="<<use("MSBSelection")<<" & "<<use("LSBSelection")<<";"<<endl;
	
	//vhdl<< "with "<<use("Selection")<<"select"<<endl << declare("specialBits",2)<<" <= "<< " \"00\" when "..
	vhdl<<declare("specialBits",2)<<" <= "<<use("Selection")<<";"<<endl;
	
	//assembling the result

	vhdl<<"O<="<<use("specialBits")<<"&"<<use("signSignal")<<"&"<<use("convertedExponentAfterRounding")<<"&"<<use("convertedFractionAfterRounding")<<";"<<endl;
	
	
	}
	else
	{
		//code for the case when we do not need a rounding
	
	int maximalOutputValue = wF+1;
		
		
	//code for zero detector of the input

	
	vhdl<<declare("minusOne4ZD",MSB -LSB)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<MSB-LSB<<");"<<endl;
	//vhdl<<declare("carry4ZD",1)<<" <= '0'; "<<endl;
	
	zeroD = new IntAdder(target, MSB-LSB );
	zeroD->changeName(getName()+"zeroD");
	oplist.push_back(zeroD);
	inPortMap  (zeroD, "X", use("passedInput"));
	inPortMap  (zeroD, "Y", "minusOne4ZD");
	inPortMapCst(zeroD, "Cin", "'0'");
	//inPortMap  (zeroD, "Cin", "carry4ZD");
	outPortMap (zeroD, "R","zeroDS");
	vhdl << instance(zeroD, "zeroD");
	
	setCycleFromSignal("zeroDS");
	
	vhdl<<declare("zeroInput",1)<<"<= "<<use("zeroDS")<<of(MSB-LSB-1)<<" and not ("<<use("signSignal")<<");"<<endl;
		
	//code for the leading zeros/ones
	
	setCycleFromSignal("input2LZOC");
	
	lzocs		= new LZOCShifterSticky(target,inputWidth-1 , maximalOutputValue, intlog2(inputWidth-1), 0, -1);
	lzocs->changeName(getName()+"_LZCS");
	oplist.push_back(lzocs);
	inPortMap  (lzocs, "I", use("input2LZOC"));
	inPortMap	(lzocs,"OZb",use("signSignal"));
	outPortMap (lzocs, "Count","temporalExponent");
	outPortMap (lzocs, "O","temporalFraction");
	vhdl << instance(lzocs, "LZOC_component");
	
	sizeExponentValue=lzocs->getCountWidth();
	sizeFractionPlusOne=wF+1;
	
	syncCycleFromSignal("temporalExponent");
		
	//code for the fraction
		
	//vhdl<<declare("convertedFraction",wF)<<"<= "<<use("temporalFraction")<<range(maximalOutputValue-1,maximalOutputValue-wF)<<";"<<endl;
	
	
	setCycleFromSignal("temporalFraction");
	
	int sizeFractionPlusOne=wF+1;
	
	vhdl<<declare("tfr",sizeFractionPlusOne) <<"<="<<use("temporalFraction")<<range(maximalOutputValue-1,maximalOutputValue-sizeFractionPlusOne)<<";"<<endl;
		
	vhdl<<declare("sign2vector",sizeFractionPlusOne)<<"<="<<"(others=>"<<use("signSignal")<<");"<<endl;
	vhdl<<declare("tempConvert",sizeFractionPlusOne)<<"<="<<use("sign2vector")<<" xor "<<use("tfr")<<";"<<endl;
	//vhdl<<declare("tempAddSign",sizeFractionPlusOne)<<"<="<<"CONV_STD_LOGIC_VECTOR("<< use("signSignal")<<","<<sizeFractionPlusOne<<");"<<endl;
	vhdl<<declare("tempPaddingAddSign",sizeFractionPlusOne)<<"<=(others=>'0');"<<endl;
	vhdl<<declare("tempAddSign",sizeFractionPlusOne+1)<<"<="<<use("tempPaddingAddSign")<<" & "<<use("signSignal") <<";"<<endl;
	vhdl<<declare("tempConvert0",sizeFractionPlusOne+1)<<"<= '0' & "<<use("tempConvert")<<";"<<endl;
	
	
		//Integer adder for obtaining the fraction value
	
	fractionConvert = new IntAdder(target,sizeFractionPlusOne+1);
	fractionConvert->changeName(getName()+"_fractionConvert");
	oplist.push_back(fractionConvert);
	inPortMap  (fractionConvert, "X", use("tempConvert0"));
	inPortMap  (fractionConvert, "Y", use("tempAddSign"));
	//inPortMap  (fractionConvert, "X", "tempConvert");
	//inPortMap  (fractionConvert, "Y", "tempAddSign");
	inPortMapCst(fractionConvert, "Cin", "'0'"); 
	outPortMap (fractionConvert, "R","tempFractionResult");
	vhdl << instance(fractionConvert, "fractionConverter");
	
	syncCycleFromSignal("tempFractionResult");
	
	vhdl<<declare("correctingExponent",1)<<"<="<<use("tempFractionResult")<<of(sizeFractionPlusOne)<<";"<<endl;
	vhdl<<declare("convertedFraction",wF)<<"<="<<use("tempFractionResult")<<range(wF-1,0)<<";"<<endl;
	
	
	
	//code for creating the exponent

	setCycleFromSignal("temporalExponent");

	vhdl<<declare("MSB2Signal",wE)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<MSB-2<<","<<wE<<");"<<endl;
	vhdl<<declare("zeroPadding4Exponent",wE- intlog2(inputWidth-1))<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE- intlog2(inputWidth-1)<<");"<<endl;
	vhdl<<declare("valueExponent",wE)<<"<= not("<<use("zeroPadding4Exponent")<<" & "<<use("temporalExponent")<<");"<<endl;
	//vhdl<<declare("oneBit",1)<<"<='1';"<<endl;
	
	exponentConvertion = new IntAdder(target,wE);
	exponentConvertion->changeName(getName()+"exponentConvertion");
	oplist.push_back(exponentConvertion);
	inPortMap  (exponentConvertion, "X", use("MSB2Signal"));
	inPortMap  (exponentConvertion, "Y", use("valueExponent"));
	inPortMapCst(exponentConvertion, "Cin", "'1'");
	//inPortMap  (exponentConvertion, "Cin", "oneBit");
	outPortMap (exponentConvertion, "R","partialConvertedExponent");
	vhdl << instance(exponentConvertion, "exponentConvertion");
	
	syncCycleFromSignal("partialConvertedExponent");
	
	vhdl<<declare("biassOfOnes",wE-1)<<"<=CONV_STD_LOGIC_VECTOR("<<pow(2,wE)-1<<","<<wE-1<<");"<<endl;
	vhdl<<declare("biassSignal",wE)<<"<="<<"'0' &"<<use("biassOfOnes")<<";"<<endl;
	vhdl<<declare("biassSignalBit",wE+1)<<"<="<<"'0' &"<<use("biassSignal")<<";"<<endl;
	vhdl<<declare("zeroBitExponent",1)<<"<='0';"<<endl;
	vhdl<<declare("partialConvertedExponentBit",wE+1)<<"<= '0' & "<<use("partialConvertedExponent")<<";"<<endl;
	vhdl<<declare("sign4OU",1)<<"<="<<use("partialConvertedExponent")<<of(wE-1)<<";"<<endl;
	
	exponentFinal = new IntAdder(target,wE+1);
	exponentFinal->changeName(getName()+"exponentFinal");
	oplist.push_back(exponentFinal);
	inPortMap  (exponentFinal, "X", use("partialConvertedExponentBit"));
	inPortMap  (exponentFinal, "Y", use("biassSignalBit"));
	inPortMapCst(exponentFinal, "Cin", "'0'");
	//inPortMap  (exponentFinal, "Cin", "zeroBitExponent");
	outPortMap (exponentFinal, "R","convertedExponentBit");
	vhdl << instance(exponentFinal, "exponentFinal");
	
	syncCycleFromSignal("convertedExponentBit");
	
	vhdl<<declare("OUflowSignal1",2)<<"<="<<use("convertedExponentBit")<<range(wE,wE-1)<<";"<<endl;
	
	vhdl<<declare("underflowSignal",1)<<"<= '1' when ("<<use("sign4OU")<<" = '1' and "<<use("OUflowSignal1")<<" = \"01\" ) else '0' ;"<<endl;
	vhdl<<declare("overflowSignal1",1)<<"<= '1' when ("<<use("sign4OU")<<" = '0' and "<<use("OUflowSignal1")<<" = \"10\" ) else '0' ;"<<endl;
	
	
	vhdl<<declare("zeroInput4Exponent",wE+1)<<"<=(others=>'0');"<<endl;
	vhdl<<declare("possibleConvertedExponent2",wE)<<"<= "<<use("convertedExponentBit")<<range(wE-1,0)<<";"<<endl;
	vhdl<<declare("possibleConvertedExponent20",wE+1)<<"<= '0' & "<<use("possibleConvertedExponent2")<<";"<<endl;
	vhdl<<declare("sign4OU2",1)<<"<="<<use("possibleConvertedExponent2")<<of(wE-1)<<";"<<endl;
	
	expCorrect = new IntAdder(target,wE+1);
	expCorrect->changeName(getName()+"expCorrect");
	oplist.push_back(expCorrect);
	inPortMap  (expCorrect, "X", use("possibleConvertedExponent20"));
	inPortMap  (expCorrect, "Y", use("zeroInput4Exponent"));
	inPortMap  (expCorrect, "Cin", use("correctingExponent"));
	outPortMap (expCorrect, "R","finalConvertedExponent");
	vhdl << instance(expCorrect, "expCorrect");
	
	syncCycleFromSignal("finalConvertedExponent");
	
	
	vhdl<<declare("convertedExponent",wE)<<"<= "<<use("finalConvertedExponent")<<range(wE-1,0)<<";"<<endl;
	
	vhdl<<declare("overflowSignal2",1)<<"<= '1' when ("<<use("sign4OU2")<<" = '0' and "<<use("finalConvertedExponent")<<range(wE,wE-1)<<" = \"10\" ) else '0' ;"<<endl;
	
	vhdl<<declare("overflowSignal",1)<<"<="<<use("overflowSignal2")<<" or "<<use("overflowSignal1")<<";"<<endl;
	
	
		
	
	
	//code for the special bits
	
	vhdl<<declare("MSBSelection",1)<<"<= "<<use("overflowSignal")<<";"<<endl;
	vhdl<<declare("LSBSelection",1)<<"<= not("<<use("underflowSignal")<<" and not ( "<<use("zeroInput")<<" ) ) ;"<<endl;
	vhdl<<declare("Selection",2)<<"<="<<use("MSBSelection")<<" & "<<use("LSBSelection")<<";"<<endl;
	
	//vhdl<< "with "<<use("Selection")<<"select"<<endl << declare("specialBits",2)<<" <= "<< " \"00\" when "..
	vhdl<<declare("specialBits",2)<<" <= "<<use("Selection")<<";"<<endl;
	
	
	//assembling the result

	vhdl<<"O<="<<use("specialBits")<<"&"<<use("signSignal")<<"&"<<use("convertedExponent")<<"&"<<use("convertedFraction")<<";"<<endl;
	
		
	}


}

Fix2FP::~Fix2FP() {
}


void Fix2FP::emulate(TestCase * tc)
{
	
	
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("I");
	//cout << " " << unsignedBinary(svX, (MSBI-LSBI+1))<< endl;
	
	mpz_class tmpSUB = (1 << (MSBI-LSBI+1));
	mpz_class tmpCMP = (1 << (MSBI-LSBI))-1;
	//cout << " " << unsignedBinary(tmpCMP, (MSBI-LSBI+1))<< endl;
	//cout << unsignedBinary(tmpSUB, (MSBI-LSBI+1)+1)<< endl;


	if (svX > tmpCMP){ //negative number 
		svX = svX - tmpSUB;
		//cout << "-" << unsignedBinary(-svX, (MSBI-LSBI+1))<< endl;
	}
	
	mpfr_t x;
	mpfr_init2(x, 10000); //init to infinite prec
	mpfr_set_z(x, svX.get_mpz_t(), GMP_RNDN);
	
	mpfr_t cst, tmp2;
	mpfr_init2(cst, 10000); //init to infinite prec
	mpfr_init2(tmp2, 10000); //init to infinite prec
	
	
	mpfr_set_ui(cst, 2 , GMP_RNDN);
	mpfr_set_si(tmp2, LSBI , GMP_RNDN);
	mpfr_pow(cst, cst, tmp2, GMP_RNDN);
	
	mpfr_mul(x, x, cst, GMP_RNDN);
	
	//cout << "x: ";
        //mpfr_out_str (stdout, 2,  40, x, GMP_RNDN);
	//cout << endl;
   
	
	mpfr_t myFP;
	mpfr_init2(myFP, wFR+1);
	mpfr_set(myFP, x, GMP_RNDN);
	
	//cout << "my: x";
        //mpfr_out_str (stdout, 2,  40, myFP, GMP_RNDN);
	//cout << endl;
	
	
	FPNumber  fpr(wER, wFR, myFP);
	mpz_class svR = fpr.getSignalValue();
	tc->addExpectedOutput("O", svR);

	// clean-up
	mpfr_clears(x, myFP, NULL);

}





void Fix2FP::buildStandardTestCases(TestCaseList* tcl){
	
}

