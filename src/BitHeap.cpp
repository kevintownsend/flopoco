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

namespace flopoco
{

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
		chunkDoneIndex=0;
		inConcatIndex=0;
		outConcatIndex=0;
		compressorIndex=0;
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
			maxWeight++;
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
			maxWeight++;
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

        list<WeightedBit*>::iterator it = bits[i].begin();
        stringstream signal[2];

        for(unsigned j=0; j<bc->getColumnSize(0); j++)
        {
        	signal[0] << (*it)->getName();
        	if(j!=bc->getColumnSize(0)-1)
        		signal[0] << " & ";
        	it++;
        }



    	reduce(i,bc->getColumnSize(0));
    	if(bc->getColumnSize(1)!=0)
    	{
    		it = bits[i+1].begin();
    		for(unsigned j=0; j<bc->getColumnSize(1); j++)
      		{
        		signal[1] << (*it)->getName();
        		if(j!=bc->getColumnSize(1)-1)
        			signal[1] << " & ";
        		it++;
        	}
    		reduce(i+1,bc->getColumnSize(1));
    	}
    	double maxCP = 0.0;//computeMaxCP(i,bc->getColumnSize(0),bc->getColumnSize(1));

    	op->vhdl << endl;

    	for(unsigned j=0; j<bc->height.size(); j++)
    	{
    		if (bc->getColumnSize(j)==1)
    			op->vhdl << tab << op->declare(join("in_concat_", compressorIndex,"_", inConcatIndex), bc->getColumnSize(j)) 
    							<< "(0)<=" << signal[j].str() << ";" << endl;
    		else
    			op->vhdl << tab << op->declare(join("in_concat_", compressorIndex,"_", inConcatIndex), bc->getColumnSize(j)) 
    							<< "<=" << signal[j].str() << ";" << endl;
    		op->inPortMap(bc, join("X",j), join("in_concat_", compressorIndex,"_", inConcatIndex));
    		++inConcatIndex;
    	}
    	REPORT(DEBUG, "outIndex=" << outConcatIndex);
    	//op->vhdl << tab << op->declare(join("out_concat_", compressorIndex,"_", outConcatIndex), bc->getOutputSize());
    	op->vhdl << endl;
    	op->outPortMap(bc, "R", join("out_concat_", compressorIndex,"_", outConcatIndex));
    
    	op->vhdl << tab << op->instance(bc, join("Compressor_",compressorIndex))<<endl;
    	
    	addBit(i, op->getCurrentCycle(), op->getTarget()->lutDelay() + maxCP, 
    					join("out_concat_", compressorIndex,"_", outConcatIndex, "(0)"),"");
    	addBit(i+1, op->getCurrentCycle(), op->getTarget()->lutDelay() + maxCP, 
    					join("out_concat_", compressorIndex,"_", outConcatIndex, "(1)"),"");
    	if(!((bc->getColumnSize(0)==3) && (bc->getColumnSize(1)==0)))
    		addBit(i+2, op->getCurrentCycle(), op->getTarget()->lutDelay() + maxCP, 
    						join("out_concat_", compressorIndex,"_", outConcatIndex, "(2)"),"");

    	++compressorIndex;
    	++outConcatIndex;


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
		if(&wb == NULL )
			return 0;
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
				if((col0 + col1<=maxCompressibleBits) && (intlog2(col0 + 2*col1)<=3))
				{
					REPORT(DEBUG, "col0 " << col0 <<", col1 " << col1);
					vector<int> newVect;
					newVect.push_back(col0);
					newVect.push_back(col1);
					pC.push_back(new BasicCompressor(op->getTarget(), newVect));
					REPORT(DEBUG, "size col0 "<< pC[pC.size()-1]->getColumnSize(0));
				}
        
    }


	void BitHeap::generateCompressorVHDL()
	{

		// Should be while (maxheight>2) compress();
		// then a second loop which does the final addition, pipelined

		REPORT(DEBUG, "in generateCompressorVHDL");

		generatePossibleCompressors();

		op->vhdl << tab << "-- Beginning of code generated by BitHeap::generateCompressorVHDL" << endl;

        while (getMaxHeight()>2)
        	{
         		compress();
        	}



        REPORT(DEBUG, "Column height after all compressions");
		for (int w=0; w<bits.size(); w++) 
		{
			REPORT(DEBUG, "   w=" << w << ":\t height=" << bits[w].size()); 
		}	

        for(int i=0;i<10;i++)
        {
         	if (usedCompressors[i]==true)
         	{
         		op->getOpListR().push_back(possibleCompressors[i]);
         	}
        }
        
        adderVHDL();
	

	}

	void BitHeap::adderVHDL()
	{
		stringstream inAdder0, inAdder1, outAdder;
		REPORT(DEBUG, maxWeight << "  " << minWeight);

		for(unsigned i=maxWeight-1; i>=minWeight; i--)
		{
			list<WeightedBit*>::iterator it = bits[i].begin();
			REPORT(DEBUG, bits[i].size());
			if(bits[i].size()==2)
			{
				inAdder0 << (*it)->getName();
				it++;
				inAdder1 << (*it)->getName();
			}
			else
			{
				inAdder0 << (*it)->getName();
				inAdder1 << "\'0\'";
			}

			if (i!=minWeight)
			{
				inAdder0<<" & "; 
				inAdder1<<" & ";
			}
		}

		inAdder0 << ";";
		inAdder1 << ";";

		op->vhdl << tab << op->declare("inAdder0", maxWeight-minWeight) << " <= " << inAdder0.str() << endl;
		op->vhdl << tab << op->declare("inAdder1", maxWeight-minWeight) << " <= " << inAdder1.str() << endl;

		op->vhdl << tab << op->declare("outAdder", maxWeight-minWeight+1) << " <= ('0' & inAdder0) + ('0' & inAdder1);" << endl;





   		op->vhdl << tab << "-- concatenate all the compressed chunks" << endl;
		op->vhdl << tab << op->declare("CompressionResult", maxWeight+1) << " <= outAdder " ;
		for(int i=chunkDoneIndex-1; i>=0; i--)
			op->vhdl <<  " & " << join("tempR", i);
		op->vhdl << ";" << endl;
		op->vhdl << tab << "-- End of code generated by BitHeap::generateCompressorVHDL" << endl;

		REPORT(DEBUG,"FinalAdder");

	}



	void BitHeap::compress()
	{
		Target* target=op->getTarget();

		
		
		unsigned n = target->lutInputs();
		unsigned w;
		
		vector<BasicCompressor*> compressors (n+1, (BasicCompressor*) 0);
		
		
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
			op->vhdl << tab << op->declare(join("tempR", chunkDoneIndex), w, true) << " <= " ;
			for(int i=w-1; i>=0; i--) {
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
			while(count(bits[i], op->getCurrentCycle()) >= possibleCompressors[j]->getColumnSize(0))
			{
				
				
				//compressors.push_back(possibleCompressors[j]) ;
				REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
				
				elemReduce(i,possibleCompressors[j]);
				usedCompressors[j]=true;


				
			}

			++j;

			// The remaining bits will fit one non-best compressor, find it and apply it.
			/*while(   (j < possibleCompressors.size())   &&   ( count(bits[i],op->getCurrentCycle()) > 2 )  )
			{

				if (  ( count(bits[i],op->getCurrentCycle())  ==  possibleCompressors[j]->getColumnSize(0) )  && (  count(bits[i+1],op->getCurrentCycle())  >=  possibleCompressors[j]->getColumnSize(1)  ))
				{
					elemReduce(i,possibleCompressors[j]);
					usedCompressors[j]=true;

					if (possibleCompressors[j]->getColumnSize(1)!=0)
					{
						REPORT(DEBUG,"Using Compressor " << j <<", reduce columns "<<i<<" and "<<i+1);
						++j;
					}
					else
					{
						REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
						REPORT(DEBUG, "j "<<j);
						++j;
					}
				}
				else
				{
					REPORT(DEBUG, "j "<<j);
					++j;
					
				}
			}*/





			while(   (j < possibleCompressors.size())   &&   ( count(bits[i],op->getCurrentCycle()) > 2 )   )
			{

				if(   count(bits[i],op->getCurrentCycle())  >=  possibleCompressors[j]->getColumnSize(0)   )
				{
					if(possibleCompressors[j]->getColumnSize(1)!=0)
					{

						if((count(bits[i+1], op->getCurrentCycle())>=possibleCompressors[j]->getColumnSize(1))&&(i<bits.size()-1))
						{
							
							elemReduce(i,possibleCompressors[j]);
							REPORT(DEBUG,"Using Compressor " << j <<", reduce columns "<<i<<" and "<<i+1);
							usedCompressors[j]=true;
						}
						else 
							++j;
					}
					else
					{
						
						elemReduce(i,possibleCompressors[j]);
						REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
						usedCompressors[j]=true;
						++j;
					}

				}
				else
					++j;
			}



		}
	}

}

