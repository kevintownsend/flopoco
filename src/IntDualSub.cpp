/*
  An operator which performes x-y and y-x in parallel for FloPoCo

  Author : Bogdan Pasca

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

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
#include "IntDualSub.hpp"

using namespace std;


namespace flopoco{

	IntDualSub::IntDualSub(Target* target, int wIn, int opType, map<string, double> inputDelays):
		Operator(target), wIn_(wIn), opType_(opType), inputDelays_(inputDelays)
	{
		ostringstream name;
		name << "IntDualSub_" << wIn_;
		setName(name.str());

		if (opType==0) 
			son_ = "yMx";
		else
			son_ = "xPy";
	
		// Set up the IO signals
		addInput ("X"  , wIn_);
		addInput ("Y"  , wIn_);
		addOutput("RxMy", wIn_);
		addOutput("R"+son_, wIn_);
	
		if (verbose){
			cout <<"delay for X is   "<< inputDelays["X"]<<endl;	
			cout <<"delay for Y is   "<< inputDelays["Y"]<<endl;
		}

		if (isSequential()){
			//compute the maximum input delay
			maxInputDelay = 0;
			map<string, double>::iterator iter;
			for (iter = inputDelays.begin(); iter!=inputDelays.end();++iter)
				if (iter->second > maxInputDelay)
					maxInputDelay = iter->second;
	
			if (verbose)
				cout << "The maximum input delay is "<<	maxInputDelay<<endl;
	
			double	objectivePeriod;
			objectivePeriod = 1/ target->frequency();
		
			if (verbose)
				cout << "Objective period is "<< objectivePeriod<<" at an objective frequency of "<<target->frequency() << endl;

			if (objectivePeriod<maxInputDelay){
				//It is the responsability of the previous components to not have a delay larger than the period
				cout << "Warning, the combinatorial delay at the input of "<<this->getName()<<"is above limit"<<endl;
				maxInputDelay = objectivePeriod;
			}
		
			if (((objectivePeriod - maxInputDelay) - target->lutDelay())<0)	{
				bufferedInputs = 1;
				maxInputDelay=0;
				target->suggestSubaddSize(chunkSize_ ,wIn_);
				nbOfChunks = ceil(double(wIn_)/double(chunkSize_));
				cSize = new int[nbOfChunks+1];
				cSize[nbOfChunks-1]=( ((wIn_%chunkSize_)==0)?chunkSize_:wIn_-(nbOfChunks-1)*chunkSize_);
				for (int i=0;i<=nbOfChunks-2;i++)
					cSize[i]=chunkSize_;				
			}
			else{
				int cS0; 
				bufferedInputs=0;
				int maxInAdd;
				target->suggestSlackSubaddSize(maxInAdd, wIn_, maxInputDelay);
				//int maxInAdd = ceil(((objectivePeriod - maxInputDelay) - target->lutDelay())/target->carryPropagateDelay()); 			
				cS0 = (maxInAdd<=wIn_?maxInAdd:wIn_);
				if ((wIn_-cS0)>0)
					{
						int newWIn = wIn_-cS0;
						target->suggestSubaddSize(chunkSize_,newWIn);
						cout << "chunckSize" << chunkSize_ << endl;
						cout << "newWIn" << newWIn << endl;
						nbOfChunks = ceil( double(newWIn)/double(chunkSize_));
				
						cSize = new int[nbOfChunks+1];
						cSize[0] = cS0;
						cSize[nbOfChunks]=( (( (wIn_-cSize[0])%chunkSize_)==0)?chunkSize_:(wIn_-cSize[0])-(nbOfChunks-1)*chunkSize_);
						for (int i=1;i<=nbOfChunks-1;i++)
							cSize[i]=chunkSize_;				
						nbOfChunks++;			
					}
				else{
					nbOfChunks=1;
					cSize = new int[1];
					cSize[0] = cS0;
				}
			}
		
			if (verbose){
				cout<<endl<<endl<<"Buffered Inputs "<<(bufferedInputs?"yes":"no")<<endl;
				cout<<endl;
				for (int i=nbOfChunks-1;i>=0;i--)
					cout<<cSize[i]<<" ";
				cout<<endl; 
			}	
		
			for (int i=0;i<nbOfChunks;i++){
				ostringstream t; t<<"X"<<i;	
				addDelaySignalBus(t.str(),cSize[i],bufferedInputs); 
				t.str(""); 
				t<<"Y"<<i;
				addDelaySignalBus(t.str(),cSize[i],bufferedInputs); 
				t.str(""); 
			}
		
			for (int i=0; i<nbOfChunks-1;i++){
				ostringstream t;
				t<<"xMycin"<<i+1<<"r"<<i;
				addDelaySignal(t.str(),cSize[i]+1,1);
				t.str("");
				t<<son_<<"cin"<<i+1<<"r"<<i;
				addDelaySignal(t.str(),cSize[i]+1,1);
			}
		
			for (int i=0; i<nbOfChunks;i++){
				ostringstream t;
				t<<"xMyr"<<i;
				addDelaySignalBus(t.str(),cSize[i],nbOfChunks-2-i);
				t.str("");
				t<<son_<<"r"<<i;
				addDelaySignalBus(t.str(),cSize[i],nbOfChunks-2-i);
			}	
		
			for (int i=0; i<nbOfChunks;i++){
				ostringstream t; t<<"sX"<<i;
				addDelaySignalBus(t.str(),cSize[i],i); t.str(""); 
				t<<"sY"<<i;
				addDelaySignalBus(t.str(),cSize[i],i); t.str("");
				if (i==0)
					addSignal("cin0",1);
			}	

			setPipelineDepth(nbOfChunks-1+bufferedInputs);

			outDelayMap["RxMy"] = target->adderDelay(cSize[nbOfChunks-1]);
			outDelayMap["R"+son_] = target->adderDelay(cSize[nbOfChunks-1]);  
			if (verbose)
				cout<< "Last addition size is "<<cSize[nbOfChunks-1]<< " having a delay of "<<target->adderDelay(cSize[nbOfChunks-1])<<endl;

		}
	
	}

	IntDualSub::~IntDualSub() {
	}


	void IntDualSub::outputVHDL(std::ostream& o, std::string name) {
		ostringstream signame;
		licence(o,"Bogdan Pasca (2008)");
		Operator::stdLibs(o);
		outputVHDLEntity(o);
		newArchitecture(o,name);
		outputVHDLSignalDeclarations(o);
		beginArchitecture(o);
	
		if(isSequential()){
			//first assignments. This level might be transformed from signals to registers
			// if the (T-inputDelay)<lutDelay
			for (int i=0;i<nbOfChunks;i++){
				//if (i==0)
				//	o << tab << "carry <= Cin; "<<endl;
				int sum=0;
				for (int j=0;j<=i;j++) sum+=cSize[j];
				o << tab << "X"<<i<<" <= X("<<sum-1<<" downto "<<sum-cSize[i]<<");"<<endl;
				o << tab << "Y"<<i<<" <= Y("<<sum-1<<" downto "<<sum-cSize[i]<<");"<<endl;
			}
			//connect first assignments to second level of signals
			for (int i=0;i<nbOfChunks;i++){
				ostringstream xi,yi;
				xi << "X"<<i;
				yi << "Y"<<i;
				o << tab << "sX"<<i<<" <= " << delaySignal(xi.str(), bufferedInputs) << ";" << endl;
				o << tab << "sY"<<i<<" <= " << delaySignal(yi.str(), bufferedInputs)<<";"<<endl;
				if (i==0)
					o << tab << "cin0 <= '"<< !(opType_) <<"';"<<endl;
			}

			//additions	for x - y	
			for (int i=0;i<nbOfChunks;i++){
				ostringstream sxi,syi;
				sxi << "sX"<<i;
				syi << "sY"<<i;
				if (i==0 && nbOfChunks>1)
					o << tab << "xMycin"<<i+1<<"r"<<i<<" <= (\"0\" & sX"<<i<<") + (\"0\" & not(sY"<<i<<")) + cin0;"<<endl;
				else 
					if (i<nbOfChunks-1)
						o << tab << "xMycin"<<i+1<<"r"<<i<<" <= ( \"0\" & " << delaySignal(sxi.str(), i)<< ")"
						  << " + ( \"0\" & not(" << delaySignal(syi.str(), i)<< "))"
						  << " + xMycin"<<i<<"r"<<i-1<<"_d("<<cSize[i-1]<<");"<<endl;
			}

			//additions	for y - x	or x + y
			for (int i=0;i<nbOfChunks;i++){
				ostringstream sxi,syi;
				sxi << "sX"<<i;
				syi << "sY"<<i;

				if (i==0 && nbOfChunks>1){
					o << tab << son_<<"cin"<<i+1<<"r"<<i<<" <= ";
					if (opType_==0)
						o<<"(\"0\" & not(sX"<<i<<")) + (\"0\" & sY"<<i<<") + cin0;"<<endl;
					else
						o<<"(\"0\" & sX"<<i<<") + (\"0\" & sY"<<i<<") + cin0;"<<endl;

				}else 
					if (i<nbOfChunks-1){
						o << tab << "yMxcin"<<i+1<<"r"<<i<<" <=";
						if (opType_==0)	
							o<<" ( \"0\" & not(" << delaySignal(sxi.str(), i)<< "))";
						else
							o<<" ( \"0\" & " << delaySignal(sxi.str(), i) << ")";
						o << " + ( \"0\" & " << delaySignal(syi.str(), i)<< ")"
						  << " + yMxcin"<<i<<"r"<<i-1<<"_d("<<cSize[i-1]<<");"<<endl;
					}
			}

			//assign the partial additions which will propagate to the result for x-y
			for (int i=0;i<nbOfChunks;i++){
				ostringstream sxi,syi;
				sxi << "sX"<<i;
				syi << "sY"<<i;
				if (i<nbOfChunks-1)
					o << tab << "xMyr"<<i<<" <= xMycin"<<i+1<<"r"<<i<<"_d("<<cSize[i]-1<<" downto 0);"<<endl;
				else{
					o << tab << "xMyr"<<i<<" <= " << delaySignal(sxi.str(), i)
					  << " + not(" << delaySignal(syi.str(), i)<<")";
					if (nbOfChunks>1)				
						o << " + xMycin"<<i<<"r"<<i-1<<"_d("<<cSize[i-1]<<");"<<endl;
					else
						o << " + cin0;"<<endl;	
				}
			}

			//assign the partial additions which will propagate to the result for y-x || x +y
			for (int i=0;i<nbOfChunks;i++){
				ostringstream sxi,syi;
				sxi << "sX"<<i;
				syi << "sY"<<i;
				if (i<nbOfChunks-1)
					o << tab << son_<<"r"<<i<<" <= "<<son_<<"cin"<<i+1<<"r"<<i<<"_d("<<cSize[i]-1<<" downto 0);"<<endl;
				else{
					o << tab << son_<<"r"<<i<<" <= not(" << delaySignal(sxi.str(), i)<<")"<<
						" + " << delaySignal(syi.str(), i);
					if (nbOfChunks>1)				
						o << " + "<<son_<<"cin"<<i<<"r"<<i-1<<"_d("<<cSize[i-1]<<");"<<endl;
					else
						o << " + cin0;"<<endl;	
				}
			}

			//assign output by composing the result for x - y
			o << tab << "RxMy <= ";
			for (int i=nbOfChunks-1;i>=0;i--){
				ostringstream xmyri;
				xmyri << "xMyr"<<i;
				if (i==0)
					o << delaySignal(xmyri.str(), nbOfChunks-2-i)<<";"<<endl;
				else
					o << delaySignal(xmyri.str(), nbOfChunks-2-i)<<" & ";			
			} 
			o<<endl;

			//assign output by composing the result for y - x || x + y
			o << tab << "R" << son_ << " <= ";
			for (int i=nbOfChunks-1;i>=0;i--){
				ostringstream ri;
				ri << son_ << "r"<<i;
				if (i==0)
					o <<  delaySignal(ri.str(), nbOfChunks-2-i)<<";"<<endl;
				else
					o << delaySignal(ri.str(), nbOfChunks-2-i)<<" & ";			
			} o<<endl;

			outputVHDLRegisters(o);
		}
		else{
			o << tab << "RxMy <= X + not(Y) + '1';" <<endl;
			o << tab << "R"<<son_<<" <= "<< (opType_==0?"not(X)":"X")<<" + Y + '1';" <<endl;
		}
		o << "end architecture;" << endl << endl;
	}





	void IntDualSub::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		mpz_class svRxMy = svX - svY;
		tc->addExpectedOutput("RxMy", svRxMy);

		mpz_class svR2;
		if (opType_==0)
			svR2=svY-svX;  
		else {
			svR2=svX+svY;
			// Don't allow overflow
			mpz_clrbit(svR2.get_mpz_t(),wIn_); 
		}
		tc->addExpectedOutput("R"+son_, svR2);
	}


	void IntDualSub::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		tc = new TestCase(this); 
		tc->addInput("X", mpz_class(0) );
		tc->addInput("Y", mpz_class(1));
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this); 
		tc->addInput("X", mpz_class(0) );
		tc->addInput("Y", mpz_class(-1));
		emulate(tc);
		tcl->add(tc);

	
	}






#if 0 // to keep the FIXME below

	// FIXME doesn't work for:    flopoco  -frequency=500 IntDualSub 26 0 TestBench 10000
	void IntDualSub::fillTestCase(mpz_class a[])
	{
		mpz_class& svX = a[0];
		mpz_class& svY = a[1];
		mpz_class& svR1 = a[2];
		mpz_class& svR2 = a[3];

		svX = getLargeRandom(wIn_-1);	
		svY = getLargeRandom(wIn_-1);

		svR1 = svX - svY;
		svR2 = svY - svX;
#ifndef _WIN32
		if(verbose)
			cout<<endl<< "x is "<< svX <<" y is "<<svY << " x-y is "<<	svR1 << " y-x is "<<svR2<<endl;
#endif
		// Don't allow overflow
		mpz_clrbit(svR1.get_mpz_t(),wIn_);
		mpz_clrbit(svR2.get_mpz_t(),wIn_);  
	}
#endif




}
