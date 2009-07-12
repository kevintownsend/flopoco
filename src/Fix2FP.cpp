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


Fix2FP::Fix2FP(Target* target, int LSBI, int MSBI, int Signed,int wER, int wFR) :
	Operator(target), MSBI(MSBI), LSBI(LSBI), Signed(Signed),wER(wER), wFR(wFR) {

	ostringstream name;
	
	if ((MSBI < LSBI)){
		cerr << " Fix2FP: Input constraint LSB <= MSB not met."<<endl;
		exit (EXIT_FAILURE);
	}
	
	mpz_class maxExpWE = mpz_class(1)<<(wER-1);
	mpz_class minExpWE = 1 - maxExpWE;
	
	if (( maxExpWE < MSBI ) || ( minExpWE > LSBI)){
		cerr << " The exponent is too small for full coverage. Try increasing the exponent !"<<endl;
		exit (EXIT_FAILURE);
	}
	

	int absMSB = MSBI>=0?MSBI:-MSBI;
	int absLSB = LSBI>=0?LSBI:-LSBI;
	name<<"Fix2FP_"<< (LSBI<0?"M":"")<<absLSB<<"_"<<(MSBI<0?"M":"")<<absMSB <<"_"<< (Signed==1?"S":"US") << "_" <<wER<<"_"<<wFR; 
	setName(name.str()); 

	setCopyrightString("Radu Tudoran, Bogdan Pasca (2009)");		

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
	
	vhdl << tab << declare("input",inputWidth) << " <= I;"<<endl;
	
	// code for the LZOCShifter part
	
	if(Signed!=0)
		vhdl << tab << declare("signSignal",1)<<"<="<<use("input")<<of(MSB-1-LSB)<<";"<<endl;
	else
		vhdl << tab << declare("signSignal",1)<<"<= '0';"<<endl;
	
	vhdl << tab << declare("passedInput",inputWidth)<<"<="<<use("input")<<range(MSB-1 -LSB,0)<<";"<<endl;
	vhdl << tab << declare("input2LZOC",inputWidth-1)<<"<="<<use("passedInput")<<range(MSB-2 -LSB,0)<<";"<<endl;
	
	if ( MSB-LSB > wF ){ //ROUNDING NEEDED
		
		//code for the case when a rounding is required
		int maximalOutputValue = (MSB-LSB)>wF+4? MSB-LSB:wF+4;
	
		if(Signed!=0){	
			lzocs = new LZOCShifterSticky(target,inputWidth-1 , maximalOutputValue, intlog2(inputWidth-1), 0, -1);
			lzocs->changeName(getName()+"_LZOCS");
			oplist.push_back(lzocs);

			inPortMap  (lzocs, "I", use("input2LZOC"));
			inPortMap  (lzocs,"OZb",use("signSignal"));
			outPortMap (lzocs, "Count","temporalExponent");
			outPortMap (lzocs, "O","temporalFraction");
			vhdl << instance(lzocs, "LZOC_component");
	
			sizeExponentValue=lzocs->getCountWidth();
		}else{
		
			lzcs = new LZOCShifterSticky(target, inputWidth , maximalOutputValue, intlog2(inputWidth), 0, 0);
			lzcs->changeName(getName()+"_LZCS");
			oplist.push_back(lzcs);
			
			inPortMap  (lzcs, "I", use("passedInput"));
			outPortMap (lzcs, "Count","temporalExponent");
			outPortMap (lzcs, "O","temporalFraction");
			vhdl << instance(lzcs, "LZC_component");
	
			sizeExponentValue=lzcs->getCountWidth();
		}
		sizeFractionPlusOne=wF+1;
		
		syncCycleFromSignal("temporalExponent");
	
		//Code for creating the exponent
		if(Signed!=0)
			vhdl << tab << declare("MSB2Signal",wE)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<MSB-2<<","<<wE<<");"<<endl;
		else
			vhdl << tab << declare("MSB2Signal",wE)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<MSB-1<<","<<wE<<");"<<endl;
	
		if(Signed!=0)
			vhdl << tab << declare("zeroPadding4Exponent",wE- intlog2(inputWidth-1))<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE- intlog2(inputWidth-1)<<");"<<endl;
		else
			vhdl << tab << declare("zeroPadding4Exponent",wE- intlog2(inputWidth))<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE- intlog2(inputWidth)<<");"<<endl;
	
		vhdl << tab << declare("valueExponent",wE)<<"<= not("<<use("zeroPadding4Exponent")<<" & "<<use("temporalExponent")<<");"<<endl;
	
		exponentConversion = new IntAdder(target,wE);
		exponentConversion->changeName(getName()+"exponentConversion");
		oplist.push_back(exponentConversion);
		inPortMap  (exponentConversion, "X", use("MSB2Signal"));
		inPortMap  (exponentConversion, "Y", use("valueExponent"));
		inPortMapCst(exponentConversion, "Cin", "'1'");
		outPortMap (exponentConversion, "R","partialConvertedExponent");
		vhdl << instance(exponentConversion, "exponentConversion");
	
		syncCycleFromSignal("partialConvertedExponent");
	
		vhdl << tab << declare("biassOfOnes",wE-1)<<"<=CONV_STD_LOGIC_VECTOR("<<pow(2,wE)-1<<","<<wE-1<<");"<<endl;
		vhdl << tab << declare("biassSignal",wE)<<"<="<<"'0' &"<<use("biassOfOnes")<<";"<<endl;
		vhdl << tab << declare("biassSignalBit",wE+1)<<"<="<<"'0' &"<<use("biassSignal")<<";"<<endl;
		vhdl << tab << declare("partialConvertedExponentBit",wE+1)<<"<= '0' & "<<use("partialConvertedExponent")<<";"<<endl;
		vhdl << tab << declare("sign4OU",1)<<"<="<<use("partialConvertedExponent")<<of(wE-1)<<";"<<endl;
	
		exponentFinal = new IntAdder(target,wE+1);
		exponentFinal->changeName(getName()+"exponentFinal");
		oplist.push_back(exponentFinal);
		inPortMap  (exponentFinal, "X", use("partialConvertedExponentBit"));
		inPortMap  (exponentFinal, "Y", use("biassSignalBit"));
		inPortMapCst(exponentFinal, "Cin", "'0'");
		outPortMap (exponentFinal, "R","convertedExponentBit");
		vhdl << instance(exponentFinal, "exponentFinal");
	
		syncCycleFromSignal("convertedExponentBit");
	
		vhdl << tab << declare("convertedExponent",wE)<<"<= "<<use("convertedExponentBit")<<range(wE-1,0)<<";"<<endl;
	
		vhdl << tab << declare("underflowSignal",1)<<"<= '1' when ("<<use("sign4OU")<<" = '1' and "<<use("convertedExponentBit")<<range(wE,wE-1)<<" = \"01\" ) else '0' ;"<<endl;
		vhdl << tab << declare("overflowSignal",1)<<"<= '1' when ("<<use("sign4OU")<<" = '0' and "<<use("convertedExponentBit")<<range(wE,wE-1)<<" = \"10\" ) else '0' ;"<<endl;
	
		//code for verifing if the number is zero
	
		setCycleFromSignal("passedInput");
	
		if(Signed!=0){
			vhdl << tab << declare("minusOne4ZD",MSB -LSB)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<MSB-LSB<<");"<<endl;
	
			zeroD = new IntAdder(target, MSB-LSB );
			zeroD->changeName(getName()+"zeroD");
			oplist.push_back(zeroD);
			inPortMap  (zeroD, "X", use("passedInput"));
			inPortMap  (zeroD, "Y", use("minusOne4ZD"));
			inPortMapCst(zeroD, "Cin", "'0'");
			outPortMap (zeroD, "R","zeroDS");
			vhdl << instance(zeroD, "zeroD");
	
			syncCycleFromSignal("zeroDS");
	
			vhdl << tab << declare("zeroInput",1)<<"<= "<<use("zeroDS")<<of(MSB-LSB-1)<<" and not("<<use("signSignal")<<");"<<endl;
		}else{

			vhdl << tab << declare("minusOne4ZD",MSB -LSB+1)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<MSB-LSB+1<<");"<<endl;
			vhdl << tab << declare("passedInputBit",MSB-LSB+1)<<"<= '0' & "<<use("passedInput")<<";"<<endl;
			zeroD = new IntAdder(target, MSB-LSB +1);
			zeroD->changeName(getName()+"zeroD");
			oplist.push_back(zeroD);
			inPortMap  (zeroD, "X", use("passedInputBit"));
			inPortMap  (zeroD, "Y", "minusOne4ZD");
			inPortMapCst(zeroD, "Cin", "'0'");
			outPortMap (zeroD, "R","zeroDS");
			vhdl << instance(zeroD, "zeroD");
	
			setCycleFromSignal("zeroDS");
	
			vhdl << tab << declare("zeroInput",1)<<"<= "<<use("zeroDS")<<of(MSB-LSB)<<" and not ("<<use("signSignal")<<");"<<endl;
		}

		//code for the Convertion of the fraction
		setCycleFromSignal("temporalFraction");
	
		if(Signed!=0){
	
			vhdl << tab << declare("sign2vector",maximalOutputValue)<<"<="<<"(others=>"<<use("signSignal")<<");"<<endl;
			vhdl << tab << declare("tempConvert",maximalOutputValue)<<"<="<<use("sign2vector")<<" xor "<<use("temporalFraction")<<";"<<endl;
			vhdl << tab << declare("tempConvert0",maximalOutputValue+1)<<"<= '0' & "<<use("tempConvert")<<";"<<endl;
			vhdl << tab << declare("tempPaddingAddSign",maximalOutputValue)<<"<=(others=>'0');"<<endl;
			vhdl << tab << declare("tempAddSign",maximalOutputValue+1)<<"<="<<use("tempPaddingAddSign")<<" & "<<use("signSignal")<<";"<<endl;

			//Integer adder for obtaining the fraction value
			fractionConvert = new IntAdder(target,maximalOutputValue+1);
			fractionConvert->changeName(getName()+"_fractionConvert");
			oplist.push_back(fractionConvert);
			inPortMap  (fractionConvert, "X", use("tempConvert0"));
			inPortMap  (fractionConvert, "Y", use("tempAddSign"));
			inPortMapCst(fractionConvert, "Cin", "'0'");
			outPortMap (fractionConvert, "R","tempFractionResult");
			vhdl << instance(fractionConvert, "fractionConverter");
	
			syncCycleFromSignal("tempFractionResult");
		}else{
			vhdl << tab << declare("tempFractionResult",maximalOutputValue+1)<<"<= '0' &"<<use("temporalFraction")<<";"<<endl;
		}
	
		vhdl << tab << declare("correctingExponent",1)<<"<="<<use("tempFractionResult")<<of(maximalOutputValue)<<";"<<endl;
		vhdl << tab << declare("fractionConverted",wF)<<"<="<<use("tempFractionResult")<<range(maximalOutputValue-2,maximalOutputValue -wF-1)<<";"<<endl;
	
		vhdl << tab << declare("firstBitofRest",1)<<"<="<<use("tempFractionResult")<<of(maximalOutputValue-wF-2)<<";"<<endl;
	
		//selection of mux 3
		vhdl << tab << declare("lastBitOfFraction")<<"<="<<use("tempFractionResult")<<of(maximalOutputValue-wF-1)<<";"<<endl;
	
		//Zero Compare of the bits from the remainder of magnitude
	
		int sizeOfRemainder=maximalOutputValue-sizeFractionPlusOne-1;
	
		setCycleFromSignal("tempFractionResult");
		vhdl << tab << declare("minusOne",sizeOfRemainder)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<sizeOfRemainder<<");"<<endl;
		vhdl << tab << declare("fractionRemainder",sizeOfRemainder)<<"<="<<use("tempFractionResult")<<range(sizeOfRemainder-1,0)<<";"<<endl;
		oneSubstracter = new IntAdder(target,sizeOfRemainder);
		oneSubstracter->changeName(getName()+"_oneSubstracter");
		oplist.push_back(oneSubstracter);
		inPortMap  (oneSubstracter, "X", use("fractionRemainder"));
		inPortMap  (oneSubstracter, "Y", use("minusOne"));
		inPortMapCst(oneSubstracter, "Cin", "'0'");
		outPortMap (oneSubstracter, "R","zeroFractionResult");
		vhdl << instance(oneSubstracter, "oneSubstracter");
	
		syncCycleFromSignal("zeroFractionResult");
	
		vhdl << tab << declare("zeroRemainder",1)<<"<= not( "<<"not ("<<use("tempFractionResult")<<of(sizeOfRemainder-1)<<") and "<<use("zeroFractionResult")<<of(sizeOfRemainder-1)<<");"<<endl;
	
		// signals for Muxes
	
		if(inputWidth>32&&target->frequencyMHz()>=250)
			nextCycle();
	
		//selection of mux 3
		vhdl << tab << declare("outputOfMux3",1)<<"<="<<use("lastBitOfFraction")<<";"<<endl;
		vhdl << tab << "with "<<use("zeroRemainder")<<" select "<<endl<<tab<<declare("outputOfMux2",1)<<" <= "<< use("outputOfMux3")<<" when '0', "<<endl<<tab<<"\t'1' when others;"<<endl;
		vhdl << tab << "with "<<use("firstBitofRest")<<" select "<<endl<<tab<<declare("outputOfMux1",1)<<" <= "<< use("outputOfMux2")<<" when '1', "<<endl<<tab<<"\t'0' when others;"<<endl;
		vhdl << tab << declare("possibleCorrector4Rounding",wF+wE+1)<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE<<") & "<<use("correctingExponent")<<" & "<<"CONV_STD_LOGIC_VECTOR(0,"<<wF<<");"<<endl;
		vhdl << tab << declare("concatenationForRounding",wE+wF+1)<<"<= '0' &"<<use("convertedExponent")<<" & "<<use("fractionConverted")<<";"<<endl;
	
	  	vhdl << tab << declare("testC",wE+wF+1)<<"<="<<use("concatenationForRounding")<<";"<<endl;
	  	vhdl << tab << declare("testR",wE+wF+1)<<"<="<<use("possibleCorrector4Rounding")<<";"<<endl;
	 	vhdl << tab << declare("testM",1)<<"<="<<use("outputOfMux1")<<";"<<endl;
	 
	
	
		roundingAdder = new IntAdder(target,wF+wE+1);
		roundingAdder->changeName(getName()+"roundingAdder");
		oplist.push_back(roundingAdder);
		inPortMap  (roundingAdder, "X", use("concatenationForRounding"));
		inPortMap  (roundingAdder, "Y", use("possibleCorrector4Rounding"));
		inPortMap  (roundingAdder, "Cin", use("outputOfMux1"));
		outPortMap (roundingAdder, "R","roundedResult");
		vhdl << instance(roundingAdder, "roundingAdder");
	
		syncCycleFromSignal("roundedResult");
	
	
		vhdl << tab << declare("convertedExponentAfterRounding",wE)<<"<="<<use("roundedResult")<<range(wE+wF-1,wF)<<";"<<endl;
		vhdl << tab << declare("convertedFractionAfterRounding",wF)<<"<="<<use("roundedResult")<<range(wF-1,0)<<";"<<endl;
	
		vhdl << tab << declare("MSBSelection",1)<<"<= "<<use("overflowSignal")<<" or "<<use("roundedResult")<<of(wF+wE)<<";"<<endl;
		vhdl << tab << declare("LSBSelection",1)<<"<= not("<<use("underflowSignal")<<" and not ( "<<use("zeroInput")<<" )) ;"<<endl;
		vhdl << tab << declare("Selection",2)<<"<="<<use("MSBSelection")<<" & "<<use("LSBSelection") <<" when "<<use("zeroInput")<<"='0' else \"00\"" <<";"<<endl;
		vhdl << tab << declare("specialBits",2)<<" <= "<<use("Selection")<<";"<<endl;
	
		//assembling the result

		vhdl << tab << "O<="<<use("specialBits")<<"&"<<use("signSignal")<<"&"<<use("convertedExponentAfterRounding")<<"&"<<use("convertedFractionAfterRounding")<<";"<<endl;
	} 
	//==========================================================================
	//==========================================================================
	else //ROUNDING NOT NEEDED
	{
	
		int maximalOutputValue = wF+1;
		
		//code for zero detector of the input

		if(Signed!=0){
			vhdl << tab << declare("minusOne4ZD",MSB -LSB)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<MSB-LSB<<");"<<endl;
	
			zeroD = new IntAdder(target, MSB-LSB );
			zeroD->changeName(getName()+"zeroD");
			oplist.push_back(zeroD);
			inPortMap  (zeroD, "X", use("passedInput"));
			inPortMap  (zeroD, "Y", "minusOne4ZD");
			inPortMapCst(zeroD, "Cin", "'0'");
			outPortMap (zeroD, "R","zeroDS");
			vhdl << instance(zeroD, "zeroD");
	
			setCycleFromSignal("zeroDS");
	
			vhdl << tab << declare("zeroInput",1)<<"<= "<<use("zeroDS")<<of(MSB-LSB-1)<<" and not ("<<use("signSignal")<<");"<<endl;
		}else{
			vhdl << tab << declare("minusOne4ZD",MSB -LSB+1)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<-1<<","<<MSB-LSB+1<<");"<<endl;
			vhdl << tab << declare("passedInputBit",MSB-LSB+1)<<"<= '0' & "<<use("passedInput")<<";"<<endl;
			zeroD = new IntAdder(target, MSB-LSB +1);
			zeroD->changeName(getName()+"zeroD");
			oplist.push_back(zeroD);
			inPortMap  (zeroD, "X", use("passedInputBit"));
			inPortMap  (zeroD, "Y", "minusOne4ZD");
			inPortMapCst(zeroD, "Cin", "'0'");
			outPortMap (zeroD, "R","zeroDS");
			vhdl << instance(zeroD, "zeroD");
	
			setCycleFromSignal("zeroDS");
	
			vhdl << tab << declare("zeroInput",1)<<"<= "<<use("zeroDS")<<of(MSB-LSB)<<" and not ("<<use("signSignal")<<");"<<endl;
		}

		//code for the leading zeros/ones
		setCycleFromSignal("input2LZOC");
	
		if(Signed!=0){
			lzocs		= new LZOCShifterSticky(target,inputWidth-1 , maximalOutputValue, intlog2(inputWidth-1), 0, -1);
			lzocs->changeName(getName()+"_LZCS");
			oplist.push_back(lzocs);
			inPortMap  (lzocs, "I", use("input2LZOC"));
			inPortMap	(lzocs,"OZb",use("signSignal"));
			outPortMap (lzocs, "Count","temporalExponent");
			outPortMap (lzocs, "O","temporalFraction");
			vhdl << instance(lzocs, "LZOC_component");
	
			sizeExponentValue=lzocs->getCountWidth();
		}else{
			lzcs = new LZOCShifterSticky(target, inputWidth , maximalOutputValue, intlog2(inputWidth), 0, 0);
			lzcs->changeName(getName()+"_LZCS");
			oplist.push_back(lzcs);
			inPortMap  (lzcs, "I", use("passedInput"));
			outPortMap (lzcs, "Count","temporalExponent");
			outPortMap (lzcs, "O","temporalFraction");
			vhdl << instance(lzcs, "LZC_component");
	
			sizeExponentValue=lzcs->getCountWidth();
		}	
		sizeFractionPlusOne=wF+1;
	
		syncCycleFromSignal("temporalExponent");

		//code for the fraction
		setCycleFromSignal("temporalFraction");
	
		int sizeFractionPlusOne=wF+1;

		if(Signed!=0){
			vhdl << tab << declare("tfr",sizeFractionPlusOne) <<"<="<<use("temporalFraction")<<range(maximalOutputValue-1,maximalOutputValue-sizeFractionPlusOne)<<";"<<endl;
		
			vhdl << tab << declare("sign2vector",sizeFractionPlusOne)<<"<="<<"(others=>"<<use("signSignal")<<");"<<endl;
			vhdl << tab << declare("tempConvert",sizeFractionPlusOne)<<"<="<<use("sign2vector")<<" xor "<<use("tfr")<<";"<<endl;
			vhdl << tab << declare("tempPaddingAddSign",sizeFractionPlusOne)<<"<=(others=>'0');"<<endl;
			vhdl << tab << declare("tempAddSign",sizeFractionPlusOne+1)<<"<="<<use("tempPaddingAddSign")<<" & "<<use("signSignal") <<";"<<endl;
			vhdl << tab << declare("tempConvert0",sizeFractionPlusOne+1)<<"<= '0' & "<<use("tempConvert")<<";"<<endl;
	
			//Integer adder for obtaining the fraction value
	
			fractionConvert = new IntAdder(target,sizeFractionPlusOne+1);
			fractionConvert->changeName(getName()+"_fractionConvert");
			oplist.push_back(fractionConvert);
			inPortMap  (fractionConvert, "X", use("tempConvert0"));
			inPortMap  (fractionConvert, "Y", use("tempAddSign"));
			inPortMapCst(fractionConvert, "Cin", "'0'"); 
			outPortMap (fractionConvert, "R","tempFractionResult");
			vhdl << instance(fractionConvert, "fractionConverter");
	
			syncCycleFromSignal("tempFractionResult");
		}else{
			vhdl << tab << declare("tempFractionResult",maximalOutputValue+1)<<"<= '0' &"<<use("temporalFraction")<<";"<<endl;
		}
	
		vhdl << tab << declare("correctingExponent",1)<<"<="<<use("tempFractionResult")<<of(sizeFractionPlusOne)<<";"<<endl;
		vhdl << tab << declare("convertedFraction",wF)<<"<="<<use("tempFractionResult")<<range(wF-1,0)<<";"<<endl;
	
		//code for creating the exponent
		setCycleFromSignal("temporalExponent");
		if(Signed!=0)
			vhdl << tab << declare("MSB2Signal",wE,true)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<MSB-2<<","<<wE<<");"<<endl;
		else
			vhdl << tab << declare("MSB2Signal",wE,true)<<"<="<<"CONV_STD_LOGIC_VECTOR("<<MSB-1<<","<<wE<<");"<<endl;
	
		if(Signed!=0)
			vhdl << tab << declare("zeroPadding4Exponent",wE- intlog2(inputWidth-1),true)<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE- intlog2(inputWidth-1)<<");"<<endl;
		else
			vhdl << tab << declare("zeroPadding4Exponent",wE- intlog2(inputWidth),true)<<"<="<<"CONV_STD_LOGIC_VECTOR(0,"<<wE- intlog2(inputWidth)<<");"<<endl;
	
		vhdl << tab << declare("valueExponent",wE,true)<<"<= not("<<use("zeroPadding4Exponent")<<" & "<<use("temporalExponent")<<");"<<endl;
	
		exponentConversion = new IntAdder(target,wE);
		exponentConversion->changeName(getName()+"exponentConversion");
		oplist.push_back(exponentConversion);
		inPortMap  (exponentConversion, "X", use("MSB2Signal"));
		inPortMap  (exponentConversion, "Y", use("valueExponent"));
		inPortMapCst(exponentConversion, "Cin", "'1'");
		outPortMap (exponentConversion, "R","partialConvertedExponent");
		vhdl << instance(exponentConversion, "exponentConversion");
	
		syncCycleFromSignal("partialConvertedExponent");
	
		vhdl << tab << declare("biassOfOnes",wE-1)<<"<=CONV_STD_LOGIC_VECTOR("<<pow(2,wE)-1<<","<<wE-1<<");"<<endl;
		vhdl << tab << declare("biassSignal",wE)<<"<="<<"'0' &"<<use("biassOfOnes")<<";"<<endl;
		vhdl << tab << declare("biassSignalBit",wE+1)<<"<="<<"'0' &"<<use("biassSignal")<<";"<<endl;
		vhdl << tab << declare("zeroBitExponent",1)<<"<='0';"<<endl;
		vhdl << tab << declare("partialConvertedExponentBit",wE+1)<<"<= '0' & "<<use("partialConvertedExponent")<<";"<<endl;
		vhdl << tab << declare("sign4OU",1)<<"<="<<use("partialConvertedExponent")<<of(wE-1)<<";"<<endl;
	
		exponentFinal = new IntAdder(target,wE+1);
		exponentFinal->changeName(getName()+"exponentFinal");
		oplist.push_back(exponentFinal);
		inPortMap  (exponentFinal, "X", use("partialConvertedExponentBit"));
		inPortMap  (exponentFinal, "Y", use("biassSignalBit"));
		inPortMapCst(exponentFinal, "Cin", "'0'");
		outPortMap (exponentFinal, "R","convertedExponentBit");
		vhdl << instance(exponentFinal, "exponentFinal");
	
		syncCycleFromSignal("convertedExponentBit");
	
		vhdl << tab << declare("OUflowSignal1",2)<<"<="<<use("convertedExponentBit")<<range(wE,wE-1)<<";"<<endl;
	
		vhdl << tab << declare("underflowSignal",1)<<"<= '1' when ("<<use("sign4OU")<<" = '1' and "<<use("OUflowSignal1")<<" = \"01\" ) else '0' ;"<<endl;
	
		vhdl << tab << declare("overflowSignal1",1)<<"<= '1' when ("<<use("sign4OU")<<" = '0' and "<<use("OUflowSignal1")<<" = \"10\" ) else '0' ;"<<endl;
	
		syncCycleFromSignal("correctingExponent");
	
		vhdl << tab << declare("zeroInput4Exponent",wE+1)<<"<=(others=>'0');"<<endl;
		vhdl << tab << declare("possibleConvertedExponent2",wE)<<"<= "<<use("convertedExponentBit")<<range(wE-1,0)<<";"<<endl;
		vhdl << tab << declare("possibleConvertedExponent20",wE+1)<<"<= '0' & "<<use("possibleConvertedExponent2")<<";"<<endl;
		vhdl << tab << declare("sign4OU2",1)<<"<="<<use("possibleConvertedExponent2")<<of(wE-1)<<";"<<endl;
	
		expCorrect = new IntAdder(target,wE+1);
		expCorrect->changeName(getName()+"expCorrect");
		oplist.push_back(expCorrect);
		inPortMap  (expCorrect, "X", use("possibleConvertedExponent20"));
		inPortMap  (expCorrect, "Y", use("zeroInput4Exponent"));
		inPortMap  (expCorrect, "Cin", use("correctingExponent"));
		outPortMap (expCorrect, "R","finalConvertedExponent");
		vhdl << instance(expCorrect, "expCorrect");
	
		syncCycleFromSignal("finalConvertedExponent");
	
		vhdl << tab << declare("convertedExponent",wE)<<"<= "<<use("finalConvertedExponent")<<range(wE-1,0)<<";"<<endl;
		vhdl << tab << declare("overflowSignal2",1)<<"<= '1' when ("<<use("sign4OU2")<<" = '0' and "<<use("finalConvertedExponent")<<range(wE,wE-1)<<" = \"10\" ) else '0' ;"<<endl;
		vhdl << tab << declare("overflowSignal",1)<<"<="<<use("overflowSignal2")<<" or "<<use("overflowSignal1")<<";"<<endl;

		//code for the special bits
	
		vhdl << tab << declare("MSBSelection",1)<<"<= "<<use("overflowSignal")<<";"<<endl;
		vhdl << tab << declare("LSBSelection",1)<<"<= not("<<use("underflowSignal")<<" or  "<<use("zeroInput")<<"  ) ;"<<endl;
		vhdl << tab << declare("Selection",2)<<"<="<<use("MSBSelection")<<" & "<<use("LSBSelection") <<" when "<<use("zeroInput")<<"='0' else \"00\"" <<";"<<endl;
		vhdl << tab << declare("specialBits",2)<<" <= "<<use("Selection")<<";"<<endl;
		
		//assembling the result
		vhdl << tab << "O<="<<use("specialBits")<<"&"<<use("signSignal")<<"&"<<use("convertedExponent")<<"&"<<use("convertedFraction")<<";"<<endl;
	}
}


Fix2FP::~Fix2FP() {
}


void Fix2FP::emulate(TestCase * tc)
{
	
	/* Get I/O values */
	mpz_class svX = tc->getInputValue("I");
 	
	mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI+1));
	mpz_class tmpCMP = (mpz_class(1)  << (MSBI-LSBI))-1;

	if (Signed != 0)
		if (svX > tmpCMP){ //negative number 
			svX = svX - tmpSUB;
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

	mpfr_t myFP;
	mpfr_init2(myFP, wFR+1);
	mpfr_set(myFP, x, GMP_RNDN);

	FPNumber  fpr(wER, wFR, myFP);
	mpz_class svR = fpr.getSignalValue();
	tc->addExpectedOutput("O", svR);

	// clean-up
	mpfr_clears(x, myFP, NULL);

}


void Fix2FP::buildStandardTestCases(TestCaseList* tcl){
	
}

