/*
 * An integer adder for FloPoCo
 *
 * It may be pipelined to arbitrary frequency.
 * Also useful to derive the carry-propagate delays for the subclasses of Target
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
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntAdder.hpp"

using namespace std;


IntAdder::IntAdder(Target* target, int wIn, const int p) :
	Operator(target), wIn_(wIn), forcePipeline_(p) {

	setOperatorName();
	setOperatorType();
	// Set up the IO signals
	addInput ("X"  , wIn_);
	addInput ("Y"  , wIn_);
	addInput ("Cin", 1  );
	addOutput("R"  , wIn_);


	if (isSequential()){
		//the maximum chunk size so that (optimisically) the target freqency can be reached. 
	bool status = target->suggestSubaddSize(chunkSize_ ,wIn_);
	
	if (!status)
		cout << "Warning: the desired frequency is not possible; optimizing to maximum frequency"<<endl;

		//the pipeLevels_ gives th number of chunks that the addition must be split in so that the frequency will ve reached
		
		//set the pipeline depth 
		if (wIn_ <= chunkSize_)
			pipeLevels_ = 0;
		else
			pipeLevels_ = wIn_ / chunkSize_;
			
		setPipelineDepth(pipeLevels_+p); 
		if(verbose) {
			cout << tab <<"Estimated delay will be " << target->adderDelay(wIn) <<endl; 
			cout << tab << "chunk="<<chunkSize_ << " freq=" << 1e-6/target->adderDelay(chunkSize_) <<" levels="<<pipeLevels_ <<endl;
		}
		if(chunkSize_ > (wIn/(pipeLevels_+1))+2)
			chunkSize_ = (wIn/(pipeLevels_+1))+2;
		if(verbose)
			cout << tab << "after chunk balancing, chunk=="<<chunkSize_ << " freq=" << 1e-6/target->adderDelay(chunkSize_) <<"  levels="<<pipeLevels_ <<endl;
		lastChunkSize_ = wIn - (pipeLevels_)*chunkSize_;
		if(verbose)
			cout << tab << "last chunk=="<<lastChunkSize_ <<endl;

		//if no pipelining is required - then these signals are not present
		if (pipeLevels_>0)
			for(int i=0; i<=pipeLevels_; i++){
				ostringstream snamex, snamey, snamer, cOutR, iSum;
				int size;
				
				snamex <<"ix_"<<i;
				snamey <<"iy_"<<i;
				iSum<<"iSum"<<i;	//isum will contain the 1 bit extended sum of ix and iy
				snamer <<"ir_"<<i;
				cOutR<<"cOutR"<<i;
				
				if(i<pipeLevels_)
					size = chunkSize_;
				else
					size = lastChunkSize_;
							
				addDelaySignalBus(snamex.str(), size, i);
				addDelaySignalBus(snamey.str(), size, i);
				
				if (i<pipeLevels_)	
					addSignal(iSum.str(),size+1); //is larger because inputs are padded with one 0 each in MSB position
				else
					addSignal(iSum.str(),size);
				
				addDelaySignalBus(snamer.str(), size, pipeLevels_-i);
				
				if (i<pipeLevels_)
					addRegisteredSignalWithSyncReset(cOutR.str(), 1);
			}	
		
	}
	if (p==1)
		addRegisteredSignalWithSyncReset("r1reg",wIn_);
	
}

IntAdder::~IntAdder() {
}

void IntAdder::setOperatorName(){
	ostringstream name;
	name << "IntAdder_" << wIn_;
	uniqueName_ = name.str();
}

void IntAdder::outputVHDL(std::ostream& o, std::string name) {
	ostringstream signame;
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	if(isSequential()){
		if(pipeLevels_==0){
			//addition simplification
			if (forcePipeline_==0)
			o << tab << "R <= X + Y + Cin;" <<endl;      
			else
			{
			o<< tab <<" r1reg <= X + Y + Cin;"<<endl;
			o<< tab <<" R<=r1reg_d;"<<endl;
			}
			
		}
		else{
			// Initialize the chunks
			for(int i=0; i<=pipeLevels_; i++){
				if (!((i==pipeLevels_)&&(lastChunkSize_==0)))
				{
					int maxIndex;
					if(i==pipeLevels_)
						maxIndex = wIn_-1;
					else 
						maxIndex = i*chunkSize_+chunkSize_-1;
						
					ostringstream snamex, snamey, snamer;
					snamex <<"ix_"<<i;
					snamey <<"iy_"<<i;
					snamer <<"ir_"<<i;
				
					o << tab << snamex.str() << " <= ";
					if(i < pipeLevels_)
						o << "";
						o << "X(" << maxIndex<< " downto " << i*chunkSize_ << ");" << endl;
				
					o << tab << snamey.str() << " <= ";
					if(i < pipeLevels_)
						o << "";
						o << "Y(" << maxIndex<< " downto " << i*chunkSize_ << ");" << endl;
				}
			}
			// then the pipe_level adders of size chunkSize_
			for(int i=0; i<=pipeLevels_; i++){
				if (!((i==pipeLevels_)&&(lastChunkSize_==0)))
				{
					int size;
					if(i==pipeLevels_) 
						size= lastChunkSize_ -1 ;
					else  
						size = chunkSize_;
					
					ostringstream snamex, snamey, snamer, iSum, cOutR, cOutRM;
					snamex <<"ix_"<<i;
					snamey <<"iy_"<<i;
					snamer <<"ir_"<<i;
					iSum << "iSum"<<i;
					cOutR << "cOutR"<<i;
					cOutRM<< "cOutR"<<i-1;
					
					o<<endl<<endl;
					
					o << tab << iSum.str() << " <= ";
					if(i!=pipeLevels_)
						o<< "(\"0\" & ";
						o<< getDelaySignalName(snamex.str(), i) ; 
					if(i!=pipeLevels_)	
						o<< ")";
						o<< " + ";
					if(i!=pipeLevels_)	
						o<< "(\"0\" & ";
						o<< getDelaySignalName(snamey.str(), i);
					if(i!=pipeLevels_)	
						o<< ")";
					
					// add the carry in
					if(i>0) {
						ostringstream carry;
						carry <<"ir_"<<i-1;
						o  << " + " << getDelaySignalName(cOutRM.str() , 1);
					}else
						o  << " + Cin";
					
					o << ";" << endl;
					
					if(i<pipeLevels_) 
						o << tab<<snamer.str() <<" <= "<< iSum.str()<<"("<<size-1<<" downto 0);"<<endl;
					else 
						o << tab<<snamer.str() <<" <= "<< iSum.str()<<"("<<size<<" downto 0);"<<endl;
						
					if(i<pipeLevels_) 
						o << tab<< cOutR.str() <<" <= "<< iSum.str()<<"("<<size<<");"<<endl;
				}
			}
			// Then the output to R 
			for(int i=0; i<=pipeLevels_; i++){
				if (!((i==pipeLevels_)&&(lastChunkSize_==0)))
				{
					int maxIndex, size;
					if(i==pipeLevels_) {
						maxIndex = wIn_-1;
						size=lastChunkSize_;
					}
					else  {
						maxIndex = i*chunkSize_+chunkSize_-1;
						size = chunkSize_;
					}
					if (forcePipeline_==0)
					{
						ostringstream snamer;
						snamer <<"ir_"<<i;
						o << tab << "R(" << maxIndex << " downto " << i*chunkSize_ << ")  <=  "
							<<  getDelaySignalName(snamer.str(), pipeLevels_ -i) << "(" << size-1 << " downto 0);" << endl; 
					}
					else
					{
					ostringstream snamer;
						snamer <<"ir_"<<i;
						o << tab << "r1reg_d(" << maxIndex << " downto " << i*chunkSize_ << ")  <=  "
							<<  getDelaySignalName(snamer.str(), pipeLevels_ -i) << "(" << size-1 << " downto 0);" << endl; 
					o<<tab<<"R<=r1reg_d;"<<endl;
					}
					
				}
			}
		}
		outputVHDLRegisters(o);
	}
	else{
		o << tab << "R <= X + Y + Cin;" <<endl;
	}
	o << "end architecture;" << endl << endl;
}

TestIOMap IntAdder::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("X"));
	tim.add(*getSignalByName("Y"));
	tim.add(*getSignalByName("Cin"));
	tim.add(*getSignalByName("R"));
	return tim;
}

void IntAdder::fillTestCase(mpz_class a[])
{
	mpz_class& svX = a[0];
	mpz_class& svY = a[1];
	mpz_class& svC = a[2];
	mpz_class& svR = a[3];

	svR = svX + svY + svC;
	// Don't allow overflow
	mpz_clrbit(svR.get_mpz_t(),wIn_); 
}


