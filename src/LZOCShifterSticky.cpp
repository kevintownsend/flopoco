/*
 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
 *
 * Authors : Florent de Dinechin, Bogdan Pasca
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
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "LZOCShifterSticky.hpp"

using namespace std;

LZOCShifterSticky::LZOCShifterSticky(Target* target, int wIn, int wOut, bool computeSticky, const int countType) :
	Operator(target), wIn_(wIn), wOut_(wOut), computeSticky_(computeSticky), countType_(countType) {
	
	int i, p2i;

	if (countType_ < 0)
		setEntityType(generic);
	else
		setEntityType(specific);

	setOperatorName();
	setOperatorType();

	/* Set up the internal architecture signals */

	/* Terminology: 
	   - level i is the signal on which we will test wether the 2^(i-1) leading bits are zeroes.
	   - i goes from wCount_ to 0 */
	
	/*  from the end of the pipeline to the beginning, 
	    determine the sizes of all register levels */
	i=0;
	size_[0] = wOut_; /* size of the result */
	while(size_[i]-wOut_ < wIn_){
		i++;
		p2i = 1<<(i-1);
		size_[i] = size_[i-1] + p2i;
		/* Invariant: size_[i] = wOut_ + 2^i -1 */
	}
	/* the attribute that gives the number of bits of the LZO count */
	wCount_ = i;
	/* the first stage doesn't need to register zeroes */
	size_[wCount_] = 1<<wCount_;
	
	/* should be identical to : wCount_ = intlog2(wIn_+1); 
	   +1 for the case all zeroes */
	
	/* Set up the IO signals */
	addInput ("I", wIn_);
	/* The width of the lzo count output will be floor(log2(wIn_+1)). */		
	addOutput("Count", wCount_);
	addOutput("O", wOut_);
	/* if we generate a generic LZOC */
	if (entityType_==generic) 
		addInput ("OZb"); 
	/* if we require a sticky bit computation */
	if(computeSticky_) 
		addOutput("Sticky"); 
		
	if(verbose){
		cout <<endl<<"  wCount_=" << wCount_ << "    Level sizes: ";
		for(i=0; i<=wCount_; i++) 
			cout << size_[i]<<" ";
		cout <<endl;
	}

	double criticalPath = 0.0;
	for (int i=wCount_; i>=0; i--){
		ostringstream levelName, leveldName, stickyName;
		levelName << "level"  << i;
		stickyName << "sticky" << i;
		level_[i] = levelName.str();

	double stageDelay;
	if (entityType_==generic) 
		stageDelay =  2* target->localWireDelay() * (1<<i);
	else		
		stageDelay =  1.2 * target->localWireDelay() * (1<<i);		//TODO


		if (i==wCount_)
		{
			criticalPath=stageDelay;
			levelRegistered_[i] = false;
		}	
		else
			if (criticalPath + stageDelay > 1/target->frequency()) {
				criticalPath=stageDelay;
				levelRegistered_[i] = true;
				incrementPipelineDepth();
			}
	}
	
	if (isSequential())
	{
		for (int i=wCount_; i>=0; i--){
			for (int j=0; j<=i; j++)
				if (levelRegistered_[j])
					countDepth_[i]++;	
			if (levelRegistered_[i]){
				leveld_[i] = level_[i] + "_d" ;
				ostringstream levelName;
				levelName << "level"<<i;		
				addRegisteredSignalWithoutReset(levelName.str(), size_[i]);	
			}
			else{
				leveld_[i] = level_[i];
				ostringstream levelName;
				levelName << "level"<<i;		
				addSignal(levelName.str(), size_[i]);	
			}
		}
	
		for (int i=wCount_-1; i>=0; i--){
			ostringstream countName;
			countName << "count"<<i;		
			addDelaySignal(countName.str(), 1, countDepth_[i]);
		}
			
		if(computeSticky_){
			for (int j=wCount_; j>=0; j--){
				ostringstream stickyName;
				ostringstream eqVer;
				ostringstream newSticky;
				eqVer<<"eqVer"<<j;
				newSticky<<"newSticky"<<j;
				stickyName<<"sticky"<<j;
				
				if (levelRegistered_[j])			
					addRegisteredSignalWithoutReset(stickyName.str());
				else
					addSignal(stickyName.str());
			
				addSignal(eqVer.str());
				addSignal(newSticky.str());	
			}
		}						
			
		if (entityType_==generic)
			addDelaySignalNoReset("sozb",1, getPipelineDepth());
	}else /* combinatorial version */
	{
		setPipelineDepth(0);
		for (int i=wCount_; i>=0; i--){
				leveld_[i] = level_[i];
				ostringstream levelName;
				levelName << "level"<<i;		
				addSignal(levelName.str(), size_[i]);
				if(computeSticky_){
					ostringstream stickyName;
					ostringstream eqVer;
					ostringstream newSticky;
					stickyName<<"sticky"<<i;
					eqVer<<"eqVer"<<i;
					newSticky<<"newSticky"<<i;
					addSignal(stickyName.str());
					addSignal(eqVer.str());
					addSignal(newSticky.str());
				}					
		}
		
		for (int i=wCount_-1; i>=0; i--){
			ostringstream countName;
			countName << "count"<<i;		
			addSignal(countName.str(), 1);
		}
		addSignal("sozb",1);
	}
	addSignal("preCount", wCount_);	
}

LZOCShifterSticky::~LZOCShifterSticky() {
}

void LZOCShifterSticky::setEntityType(entityType_t eType){
	entityType_ = eType; 
}

int LZOCShifterSticky::getCountWidth() const{
	return wCount_;
}

void LZOCShifterSticky::setOperatorName(){
	ostringstream name; 
	ostringstream computationInterface;
	/* The computationInterface refers to the VHDL entity ports. 
	   If computation is generic, then an extra input port is 
	   available on the entity for specifying the count type */ 
	if (countType_ < 0) 
		computationInterface<<"gen";
	else 
		computationInterface<<"spec";
	
	name<<"LZOCShifterSticky_"<<wIn_<<"_"<<wOut_<<"_"<<computeSticky_<<"_" //semantic name
	    <<computationInterface.str();
	
	uniqueName_=name.str();
}

void LZOCShifterSticky::outputVHDL(std::ostream& o, std::string name) {
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	outputVHDLSignalDeclarations(o);	
	beginArchitecture(o);
	
	/* connect input to level wCount_, possibly with 0 padding */
	o << tab << "level" << wCount_<< " <= " ;
	if(wIn_ == size_[wCount_]) // no padding needed
		o << "I;" <<endl;
	else 
		o << "I & (" << size_[wCount_]-wIn_-1 << " downto 0 => '0');" << endl ;

	if (entityType_==generic)
		o<<tab<<"sozb <= OZb;"<<endl;
	
	if(computeSticky_)
		o<<tab<<"sticky" << wCount_ << " <= '0';"<<endl;

	for  (int i = wCount_; i>=1; i--){
		int p2i = 1 << i;
		int p2io2 = 1 << (i-1);

		//======================================================================
		o << tab << "count" << i-1 << " <= '1' when " << leveld_[i] << "(" << size_[i]-1 << " downto " << size_[i]-p2io2 << ") = (" << p2io2-1 << " downto 0 => "; 
		if (entityType_==generic)
			o<<getDelaySignalName("sozb",getPipelineDepth()-countDepth_[i-1]);
		else
			o<<"'"<<countType_<<"'";
		o<<" )   else '0';" << endl; 
		//======================================================================
		// REM in the followIn_g,  size_[i]-size_[i-1]-1 is almost always equal to p2io2, except in the first stage
		if (i==wCount_){
			o << tab << level_[i-1] << " <= " ;
			o << " " << "("<<leveld_[i] << "(" << size_[i]/2 -1 << " downto " << 0 << ") & "<<zeroGenerator(size_[i-1]-size_[i]/2, 0)<<" )  when count" << i-1 << "='1'" << endl;
			
			if (size_[i-1]-size_[i]>0)
				o << tab << tab << tab <<  "else (" << leveld_[i] << "(" << size_[i]-1 << " downto 0) &  "<<zeroGenerator(size_[i-1]-size_[i],0)<<") ;" << endl; 
			else
				o << tab << tab << tab <<  "else " << leveld_[i] << "(" << size_[i]-1 << " downto "<<size_[i]-size_[i-1]<<");" << endl; 
		}
		else
		{
			o << tab << level_[i-1] << " <= " ;
			o << " " << leveld_[i] << "(" << size_[i]-1 << " downto " << size_[i]-size_[i-1] << ")   when count" << i-1 << "='0'" << endl;
			o << tab << tab << tab <<  "else " << leveld_[i] << "(" << size_[i-1]-1 << " downto 0);" << endl; 
		}
		
		if(computeSticky_){
			ostringstream eqVer;
			ostringstream newSticky;
			eqVer << "eqVer"<<i;
			newSticky <<"newSticky"<<i;
			
			o << tab << eqVer.str() <<" <= '1' when " << leveld_[i] << "(" << size_[i]-size_[i-1]-1 << " downto 0) = (" << size_[i]-size_[i-1]-1 << " downto 0 => '0') else '0';"
			  << endl;
			o << tab << newSticky.str() <<" <= sticky" << i << " or not("<<eqVer.str()<<");"<<endl;  
		
			if (size_[i]-size_[i-1]-1>=0){
				o << tab << "sticky" << i-1 << " <= sticky"<<i;
				if(isSequential() && i<wCount_ && levelRegistered_[i]) 
					o << "_d";
				o <<" when (count" << i-1 << "='1') else "<<newSticky.str();
			}
			else{
				o << tab << "sticky" << i-1 << " <= sticky"<<i<<" when count" << i-1 << "='1'  else sticky" << i;	
				if(isSequential() && i<wCount_ && levelRegistered_[i]) 
					o << "_d";
			}
			o << ";" << endl;
		}
			
		o << endl;
	}

	if(computeSticky_)
		if (levelRegistered_[0])
			o << tab << "sticky <= sticky0_d;" << endl;
		else
			o << tab << "sticky <= sticky0;" << endl;
			
	o << tab << "O      <= "<<leveld_[0]<<";" << endl;
	o << tab << "preCount  <= ";
	for (int i=wCount_-1; i >=0; i--){
		ostringstream name;
		name << "count" << i;
		if (isSequential())
			o << getDelaySignalName(name.str(), countDepth_[i]);
		else
			o << name.str();
		if (i>0) o << " & ";
	}
	
	o << ";" << endl;

	ostringstream binInputWidth;
	
	printBinNum(binInputWidth, wIn_, wCount_);	

	o << tab << "Count <= \""<<binInputWidth.str()<<"\" when preCount = ("<<wCount_-1<<" downto  0 => '1')";
	o << tab << "         else preCount;"<<endl;
	 
	if (isSequential()){		
		outputVHDLRegisters(o); o<<endl;
	}
	endArchitecture(o);
}



TestIOMap LZOCShifterSticky::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("I"));
	tim.add(*getSignalByName("Count"));
	tim.add(*getSignalByName("O"));
	if (computeSticky_)
		tim.add(*getSignalByName("Sticky"));
	if (entityType_ == generic )
		tim.add(*getSignalByName("OZb"));	
	
	return tim;
}


void LZOCShifterSticky::fillTestCase(mpz_class a[])
{
	if (entityType_==specific)
	{
		mpz_class& si     = a[0];
		mpz_class& scount = a[1];
		mpz_class& so     = a[2];
		mpz_class& ssticky = a[3];
		
		int sticky=0;
		
		/* Count the leading zero/one s */
		int j = wIn_-1;                      /* the index of the MSB of the input */
		int bit = (countType_ == 0) ? 0 : 1; /* what are we counting in the specific case */
		/* from the MSB towards the LSB, check if current bit of input = to the bit we test against */
		for (j = wIn_-1; j >= 0; j--)
				if (mpz_tstbit(si.get_mpz_t(), j) != bit)
					break;
		
		/* the number of bits is then equal to the index of the MSB - the index where we stoped previously */
		int icount = wIn_ - 1 - j;
		
		/* assign this value to the scount output */
		scount = icount; 
	
		/* compute the max value on wOut_ bits */
		maxValue_=2;
		for (int i=2;i<=wOut_;i++)
			maxValue_=maxValue_*2;
		maxValue_--;
	
		mpz_class inputValue;
		inputValue=si;
		sticky=0;

		if (countType_==0) /* if we are counting zeros */
			if (inputValue!=0) 
				while (!((inputValue<=maxValue_)&&(2*inputValue>maxValue_)))
					if (inputValue>maxValue_){
						if(mpz_tstbit(inputValue.get_mpz_t(), 0)==1)
							sticky=1;
						inputValue=inputValue/2;
					}
					else
						inputValue=inputValue*2;
			else {}
		else /* if we are counting ones */
		{
			int restOfBits = wIn_ - icount;
			if (icount>0)
			{
				mpz_class ones = 1;
				for (int i=1;i<=icount;i++)
					ones = ones*2;
			
				ones=ones-1;
						
				for (int i=1;i<=restOfBits;i++)
					ones=ones*2;
				inputValue=inputValue-ones; // the input without the leading ones
			}

			if ((wIn_<=wOut_) || ((wIn_>wOut_) && (restOfBits<wOut_) ))	//shift result in place	
				for (int i=1;i<=(wOut_-restOfBits);i++)
					inputValue=inputValue*2;
			else
				for (int i=1;i<=restOfBits-wOut_;i++){
					if(mpz_tstbit(inputValue.get_mpz_t(), 0)==1)
						sticky=1;
					inputValue=inputValue/2;
				}
		}
				
		/* assign the shifted result to the testCase output */		
		so=inputValue;
				
		if (computeSticky_)
			/* assign the sticky bit computation to the output if required*/	
			ssticky = sticky;
				
		if (verbose)
		{
			cout<<"TestCase report:"<<endl;
			cout<<tab<<"Input value:"<<si<<endl;
			cout<<tab<<"Count value:"<<icount<<endl;
			cout<<tab<<"Output value:"<<so<<endl;
			if (computeSticky_)
			cout<<tab<<"Sitcky value:"<<ssticky<<endl;
		}		
			
	}
	else
	{ /* if entity type is generic */
		if (! computeSticky_)
		{
			mpz_class& si     = a[0];
			mpz_class& scount = a[1];
			mpz_class& so     = a[2];
			mpz_class& sozb   = a[3];
		
			int j=wIn_-1;
			/* Count the leading zero/one s */
			for (j = (wIn_-1); j >= 0; j--)
				if (mpz_tstbit(si.get_mpz_t(), j) != sozb)
					break;
		
			int icount =(wIn_-1)-j;
			scount = icount; 	
			
			/* compute max value on wOut_ bits */
			maxValue_=2;
			for (int i=2;i<=wOut_;i++)
				maxValue_=maxValue_*2;
			maxValue_--;
	
			mpz_class inputValue;
			inputValue=si;
			
			if (sozb==0) /* if we are counting zeros */
				if (inputValue!=0) 
					while (!((inputValue<=maxValue_)&&(2*inputValue>maxValue_)))
						if (inputValue>maxValue_)
							inputValue=inputValue/2;
						else
							inputValue=inputValue*2;
				else {}
			else /* if we are counting ones */
			{
				int restOfBits = wIn_ - icount;
				if (icount>0)
				{
					mpz_class ones=1;
					for (int i=1;i<=icount;i++)
						ones = ones*2;
			
					ones=ones-1;
						
					for (int i=1;i<=restOfBits;i++)
						ones=ones*2;
					inputValue=inputValue-ones; /* the input without the leading ones */
				}

				if ((wIn_<=wOut_) || ((wIn_>wOut_) && (restOfBits<wOut_) ))	/* shift result in place */	
					for (int i=1;i<=(wOut_-restOfBits);i++)
						inputValue=inputValue*2;
				else
					for (int i=1;i<=restOfBits-wOut_;i++)
						inputValue=inputValue/2;
			}
			
			so=inputValue;	
		}
		else{
		/* if we also require a sticky bit computation */
			mpz_class& si     = a[0];
			mpz_class& scount = a[1];
			mpz_class& so     = a[2];
			mpz_class& ssticky= a[3];
			mpz_class& sozb   = a[4];
			
			
			int sticky = 0;
			int j=wIn_-1;
			/* Count the leading zero/one s */
			for (j = (wIn_-1); j >= 0; j--)
				if (mpz_tstbit(si.get_mpz_t(), j) != sozb)
					break;
		
			int icount =(wIn_-1)-j;
			scount = icount; //the count result
			
			//compute max value on wOut_ bits
			maxValue_=2;
			for (int i=2;i<=wOut_;i++)
				maxValue_=maxValue_*2;
			maxValue_--;
	
			mpz_class inputValue;
			inputValue=si;
			sticky=0;

			if (sozb==0) /* if we are counting zeros */
				if (inputValue!=0) 
					while (!((inputValue<=maxValue_)&&(2*inputValue>maxValue_)))
						if (inputValue>maxValue_){
							if(mpz_tstbit(inputValue.get_mpz_t(), 0)==1)
								sticky=1;
							inputValue=inputValue/2;
						}
						else
							inputValue=inputValue*2;
				else {}
			else /* if we are counting ones */
			{
				int restOfBits = wIn_ - icount;
				if (icount>0)
				{
					mpz_class ones=1;
					for (int i=1;i<=icount;i++)
						ones = ones*2;
			
					ones=ones-1;
						
					for (int i=1;i<=restOfBits;i++)
						ones=ones*2;
					inputValue=inputValue-ones; // the input without the leading ones
				}

				if ((wIn_<=wOut_) || ((wIn_>wOut_) && (restOfBits<wOut_) ))	//shift result in place	
					for (int i=1;i<=(wOut_-restOfBits);i++)
						inputValue=inputValue*2;
				else
					for (int i=1;i<=restOfBits-wOut_;i++){
						if(mpz_tstbit(inputValue.get_mpz_t(), 0)==1)
							sticky=1;
						inputValue=inputValue/2;
					}
			}
			so=inputValue;
			ssticky = sticky;
		}
	}
}



