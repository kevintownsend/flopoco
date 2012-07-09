/*
  A class to manage heaps of weighted bits in FloPoCo
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2012.
  All rights reserved.

*/
#include "BitHeap.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	

#include "utils.hpp"
#include<vector>
#include<list>

using namespace std;

/*

*/

namespace flopoco{

	BitHeap::WeightedBit::WeightedBit(BitHeap* bh, int weight_, int cycle_, double criticalPath_) : 
		bh(bh), weight(weight_)
	{
		srcFileName=bh->getOp()->getSrcFileName() + ":BitHeap:WeightedBit";
		if(cycle_==-1)
			cycle=bh->getOp()->getCurrentCycle();
		else
			cycle=cycle_;
		if(criticalPath_==-1)
			criticalPath=bh->getOp()->getCriticalPath();
		else
			criticalPath=criticalPath_;
		std::ostringstream p;
		p  << "heap_w" << weight << "_" << bh->getUid(weight);
		name=p.str();

	}

	
	/* which bit was defined earlier */
	bool BitHeap::WeightedBit::operator< (const WeightedBit& b){
		if ((cycle<b.cycle) || (cycle==b.cycle && criticalPath<b.criticalPath)) 
			return true;
		else
			return false;
	} 

	bool BitHeap::WeightedBit::operator<= (const WeightedBit& b){
		if ((cycle<b.cycle) || (cycle==b.cycle && criticalPath<=b.criticalPath)) 
			return true;
		else
			return false;
	} 
	
	double BitHeap::WeightedBit::getCriticalPath(int atCycle){
		if (cycle>atCycle){
			THROWERROR("For bit " << name << ", getCriticalPath() called for cycle "<<atCycle<< " but this bit is defined only at cycle "<< cycle);
		}
		if (cycle==atCycle)
			return criticalPath;
		if (cycle<atCycle)
			return 0.0;
		
		return 0.0;  //because it returned no value on this branch
	}


	BitHeap::BitHeap(Operator* op, int maxWeight) :
		op(op), maxWeight(maxWeight)
	{
		// Set up the vector of lists of weighted bits, and the vector of uids
		srcFileName=op->getSrcFileName() + ":BitHeap"; // for REPORT to work
		REPORT(DEBUG, "Creating BitHeap of size " << maxWeight);
		for (int i=0; i< maxWeight; i++) {
			uid.push_back(0);
			list<WeightedBit*> t;
			bits.push_back(t);
		}
	}	


	BitHeap::~BitHeap()
	{
		for(unsigned i=0; i<bits.size(); i++){
			list<WeightedBit*>& l = bits[i];
			list<WeightedBit*>::iterator it;
			for(it=l.begin(); it!=l.end(); it++){
				delete (*it);
			}
		}
	}	


	int BitHeap::getUid(unsigned w){
		return uid[w]++;
	}	


	double BitHeap::computeMaxCP(int w, int c0, int c1)
	{
		double max=0.0;

		REPORT(DEBUG, "Crash here "<<w<<" and "<<bits.size());

		if(w>=bits.size())
		{
			return 0.0;
		}
		else
		{
			list<WeightedBit*>::iterator it1, it2;
			int j=0,k=0;
			//it=bits[w].begin();

			for(it1=bits[w].begin(); it1!=bits[w].end(); it1++)
			{
				if(j<c0)
					if ((*it1)->getCriticalPath(op->getCurrentCycle())>max)
						max=(*it1)->getCriticalPath(op->getCurrentCycle());
				j++;
			}

			for(it2=bits[w+1].begin(); it2!=bits[w+1].end(); it2++)
			{
				if(k<c1)
					if ((*it2)->getCriticalPath(op->getCurrentCycle())>max)
						max=(*it2)->getCriticalPath(op->getCurrentCycle());
				k++;
			}


			/*

			//A little more efficient

			while((j<c0)&&(it!=bits[w].end()))
			{
				if ((*it)->getCriticalPath(op->getCurrentCycle())>max)
					max=(*it)->getCriticalPath(op->getCurrentCycle());
				j++;
				it++;
			}
		
			if(c1!=0)
			{
				it=bits[w+1].begin();
				while((k<c1)&&(it!=bits[w+1].end()))
				{
					if ((*it)->getCriticalPath(op->getCurrentCycle())>max)
						max=(*it)->getCriticalPath(op->getCurrentCycle());
					k++;	
					it++;
				}
			}
			*/

			return max;

		}
	}







	void  BitHeap::addBit(unsigned w, string rhs, string comment)
	{
		list<WeightedBit*> t;
			
		if(bits.size()==w)
		{
			bits.push_back(t);
		}

		WeightedBit* bit= new WeightedBit(this, w) ;
		list<WeightedBit*>& l=bits[w];

		//insert so that the list is sorted by bit cycle/delay
		list<WeightedBit*>::iterator it=l.begin(); 
		bool proceed=true;
		while(proceed) {
			if (it==l.end() || (**it <= *bit)){ // test in this order to avoid segfault!
				l.insert(it, bit);
				proceed=false;
			}
			else 
				it++;
		}
		// now generate VHDL
		op->vhdl << tab << op->declare(bit->getName()) << " <= " << rhs << ";";
		if(comment.size())
			op->vhdl <<  " -- " << comment;
		op->vhdl <<  endl;		

		REPORT(DEBUG, "added bit on column " << w);	
	};





	void BitHeap::addBit(unsigned weight, int cycle,  double criticalPath, string rhs, string comment)
	{
		list<WeightedBit*> t;
			
		if(bits.size()==weight)
		{
			bits.push_back(t);
		}

		//????????????????????????????????????????????


		//REPORT(DEBUG,"CYCLE before if in addBit="<<cycle);
		//REPORT(DEBUG,"criticalPath="<<criticalPath <<" (1/op->getTarget()->frequency()) " <<1/op->getTarget()->frequency());

		if (criticalPath > (1/op->getTarget()->frequency()))
		{
			criticalPath=0.0;
			cycle++;
		}

		//REPORT(DEBUG,"CYCLE in addBit="<<cycle);
		//////////////////////////////////////////////

		WeightedBit* bit= new WeightedBit(this, weight, cycle, criticalPath);
		list<WeightedBit*>& l=bits[weight];

		//insert so that the list is sorted by bit cycle/delay
		list<WeightedBit*>::iterator it=l.begin(); 
		bool proceed=true;
		while(proceed) {
			if (it==l.end() || (**it <= *bit)){ // test in this order to avoid segfault!
				l.insert(it, bit);
				proceed=false;
			}
			else 
				it++;
		}
		// now generate VHDL
		op->vhdl << tab << op->declare(bit->getName()) << " <= " << rhs << ";";
		if(comment.size())
			op->vhdl <<  " -- " << comment;
		op->vhdl <<  endl;		

		REPORT(DEBUG, "added bit on column " << weight);	
	};






	void BitHeap::removeBit(unsigned weight, int dir){
		
		list<WeightedBit*>& l=bits[weight];
    

    //if dir=0 the bit will be removed from the begining of the list, else from the end of the list of weighted bits
		if(dir==0)
			l.pop_front();
		else if(dir==1)
			l.pop_back();

		


		REPORT(DEBUG,"remove bit from column " << weight);


	}

  
    void BitHeap::elemReduce(int i, BasicCompressor* bc)
    {

    	reduce(i,bc->getColumn(0));
    	if(bc->getColumn(1)!=0)
    		reduce(i+1,bc->getColumn(1));

    	double maxCP = computeMaxCP(i,bc->getColumn(0),bc->getColumn(1));

    	addBit(i, op->getCurrentCycle(), op->getTarget()->lutDelay() + maxCP, "rhs","comm");
    	addBit(i+1, op->getCurrentCycle(), op->getTarget()->lutDelay() + maxCP, "rhs","comm");
    	addBit(i+2, op->getCurrentCycle(), op->getTarget()->lutDelay() + maxCP, "rhs","comm");
    }





	string adderBitName(int c, int bit, int h) {
		ostringstream p;
		p  << "adderBit" << c << "_" << bit << "_" << h;



		return p.str();
	};



	unsigned BitHeap::currentHeight(unsigned w) {
		int h=0;
		list<WeightedBit*>& l = bits[w];
		h=l.size();
		return h;
	}
	


	int BitHeap::count(list<WeightedBit*> wb, int cycle)
	{
		int bitsToCompress=0;
		for (list<WeightedBit*>::iterator iter = wb.begin(), end = wb.end(); iter != end; ++iter)
		{
			if (((*iter)->todo()==true)&&((*iter)->getCycle()<=cycle)){
    				bitsToCompress++;
    			}
		}
		return bitsToCompress;
	}

	void BitHeap::reduce(int c, int red)
	{

       while(red>0)
       {
       	removeBit(c,0);
       	red--;
       }

	}


	int BitHeap::getMaxHeight()

 	{
		int max=0; 
			for(int i=0; i<bits.size();i++)
			{
              if(bits[i].size()>max)
              	max=bits[i].size();
			}
     return max;

	}
	// There are two distinct notions: compression stages, and compression cycles. Several stages can fit in one cycle.


	
     vector<BasicCompressor *> possibleCompressors (0, (BasicCompressor*) 0);
     

     void BitHeap::generatePossibleCompressors()
    {
    	//Build a vector of basic compressors that are fit for this target
		//Size 10 should cover all possible compressors at this time
		//Not very elegant, but should work
		vector<BasicCompressor *>& pC=possibleCompressors;
		int maxCompressibleBits, col0, col1;

		maxCompressibleBits = op->getTarget()->lutInputs();
		
		
		//Generate all "workable" compressors for 2 columns, on this target, descending on fitness function
		for(col0=maxCompressibleBits; col0>=3; col0--)
			for(col1=col0; col1>=0; col1--)
				if(col0 + 2*col1<=maxCompressibleBits)
				{
					vector<int> newVect;
					newVect.push_back(col1);
					newVect.push_back(col0);
					pC.push_back(new BasicCompressor(op->getTarget(), newVect));
				}
        
    }


	void BitHeap::generateCompressorVHDL()
	{

		// Should be while (maxheight>2) compress();
		// then a second loop which does the final addition, pipelined

		REPORT(DEBUG, "in generateCompressorVHDL");

		generatePossibleCompressors();

		int heapCount=0;

         while (getMaxHeight()>2)
         	{
         		compress(heapCount);
         		++heapCount;
         	}

         adderVHDL();
	

	}

   
   void BitHeap::adderVHDL()
   {
     
     REPORT(DEBUG,"FinalAdder");

   }



	void BitHeap::compress(int heapCount)
	{
		Target* target=op->getTarget();

		op->vhdl << tab << "-- Beginning of code generated by BitHeap::generateCompressorVHDL" << endl;
		
		unsigned n = target->lutInputs();
		unsigned w;
		
		vector<BasicCompressor*> compressors (n+1, (BasicCompressor*) 0);
		
		unsigned chunkDoneIndex=0;
		//unsigned stage=0;
		//bool isPipelined=target->isPipelined();

		REPORT(DEBUG, "maxWeight="<< maxWeight);
		REPORT(DEBUG, "Column height before compression");
		for (w=0; w<bits.size(); w++) {
			REPORT(DEBUG, "   w=" << w << ":\t height=" << bits[w].size()); 
		}		

		
		REPORT(DEBUG, "Checking whether rightmost columns need compressing");
		bool alreadyCompressed=true;
		w=0;
		while(w<bits.size() && alreadyCompressed) 
		{
			//REPORT(DEBUG, "Current height for level " << w << " is " << currentHeight(w));
			
			if(currentHeight(w)>1)
			{
				alreadyCompressed=false;
			}
			else
			{
				REPORT(DEBUG, "Level " << w << " is already compressed; will go directly to the final result");
				w++;
			}
		}
		
		
		
		//FIXME Two poorly done concatenations here, must be checked after the compression is fully done
		if(w!=0)
		{
			//op->vhdl << tab << op->declare(join("tempR", chunkDoneIndex), w, true) << " <= " ;
			for(int i=w-1; i>=0; i--) {
				//REPORT(DEBUG, "Current height " << currentHeight(i)<<" of " << i);
				if(currentHeight(i)==1) 
				{
					op->vhdl << (bits[i].front())->getName()   << " & ";
					bits[i].pop_front();
				}
				else
				{
					op->vhdl << "'0' & ";
				}	
			}
			op->vhdl << " \"\"; -- already compressed" << endl; 
			chunkDoneIndex++;			
		}
		
		
		
		minWeight=w;

		REPORT(DEBUG, "minWeight="<< minWeight);
		

		
		
		REPORT(DEBUG,"start compressing");
		//Remaining structure must be compressed
		unsigned j;

		// i is the column index
		for(unsigned i=minWeight; i<maxWeight; i++)
		{

			j=0; // index of the compressors among the possible ones
			

			// map the best compressor as many times as possible 
			while(count(bits[i], op->getCurrentCycle()) >= possibleCompressors[j]->getColumn(0))
			{
				
				
				compressors.push_back(possibleCompressors[j]) ;
				REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
				elemReduce(i,possibleCompressors[j]);
				
			}

			++j;

			// The remaining bits will fit one non-best compressor, find it and apply it.
			while((j<possibleCompressors.size())&&(count(bits[i], op->getCurrentCycle())>2))
			{

				if(count(bits[i], op->getCurrentCycle())==possibleCompressors[j]->getColumn(0))
				{
					if(possibleCompressors[j]->getColumn(1)!=0)
					{

						if((count(bits[i+1], op->getCurrentCycle())>=possibleCompressors[j]->getColumn(1))&&(i<bits.size()-1))
						{
							elemReduce(i,possibleCompressors[j]);
						}
						else 
							++j;
					}
					else
					{
						elemReduce(i,possibleCompressors[j]);
						++j;
					}

				}
				else
					++j;
			}



			/*while((j<possibleCompressors.size())&&(count(bits[i], op->getCurrentCycle())>2))
			{
								
				if(count(bits[i], op->getCurrentCycle())==possibleCompressors[j]->getColumn(0))
				{
					if(possibleCompressors[j]->getColumn(1)!=0)
					{
						//REPORT(DEBUG, "equality " << count(bits[i]) << possibleCompressors[j]->getColumn(0));

						if((count(bits[i+1], op->getCurrentCycle())>=possibleCompressors[j]->getColumn(1))&&(i<bits.size()-1))
						{
					
						REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
						compressors.push_back(possibleCompressors[j]);
						elemReduce(i,possibleCompressors[j]);
						
						++j;
					

				}
				else
					++j;
			}*/
		}
	};

}

