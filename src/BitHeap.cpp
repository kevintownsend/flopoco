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
			cycle = bh->getOp()->getCurrentCycle();
		else
			cycle=cycle_;
		if(criticalPath_==-1)
			criticalPath = bh->getOp()->getCriticalPath();
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
	
	double BitHeap::WeightedBit::getCriticalPath(int atCycle)
	{

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



	BitHeap::WeightedBit* BitHeap::computeLatest(int w, int c0, int c1)
	{

		if(w>=bits.size())
			{	
				return NULL;
			}

		if(c1==0)
			{
				int k=1;
				WeightedBit **bb;
				for(list<WeightedBit*>::iterator it = bits[w].begin(); it!=bits[w].end(); ++it)
					{
						if (k==c0)
							{
								bb = &*it;
							}
						k++;
					}
				return *bb;
			}
		else
			{
				int i=1, j=1;
				WeightedBit **b0;
				for(list<WeightedBit*>::iterator it = bits[w].begin(); it!=bits[w].end(); ++it)
					{
						if (i==c0)
							{
								b0 = &*it;
							}
						i++;
					}
				WeightedBit **b1;
				for(list<WeightedBit*>::iterator it = bits[w+1].begin(); it!=bits[w+1].end(); ++it)
					{
						if (j==c1)
							{
								b1 = &*it;
							}
						j++;
					}


				if((**b0) <= (**b1))
					return *b1;
				else 
					return *b0;
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
			if (it==l.end() || (*bit <= **it)){ // test in this order to avoid segfault!
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
		op->vhdl <<  " -- " << "cycle= "<< bit->getCycle() <<" cp= "<<bit->getCriticalPath(bit->getCycle());		
		op->vhdl <<  endl;		
		
		REPORT(DEBUG, "added bit on column " << w <<" cycle= "<< bit->getCycle() <<" cp= "<<bit->getCriticalPath(bit->getCycle()));
		
		printColumnInfo(w);	
	};



	void BitHeap::removeBit(unsigned weight, int dir)
	{
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

		op->vhdl << endl;

		// set time to the instant of creation of the latest bit, then advance of lutDelay()+localWireDelay()
		WeightedBit* b =  computeLatest(i, bc->getColumnSize(0), bc->getColumnSize(1)) ; 
		// REPORT(DEBUG, "name " << b->getName() << " cycle " <<b->getCycle() << " cp " << b->getCriticalPath(b->getCycle()));
		// op->vhdl << "-- name " << b->getName() << " cycle " <<b->getCycle() << endl;
		op->setCycle(  b ->getCycle()  );
		op->setCriticalPath(  b ->getCriticalPath(op->getCurrentCycle()));
		op->manageCriticalPath( op->getTarget()->lutDelay() + op->getTarget()->localWireDelay() );

		// build the first input of the compressor
		for(unsigned j=0; j<bc->getColumnSize(0); j++)
			{
				signal[0] << (*it)->getName();
				if(j!=bc->getColumnSize(0)-1)
					signal[0] << " & ";
				it++;
			}

		// build the second input of the compressor, if any
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
			}
        

		for(unsigned j=0; j<bc->height.size(); j++)
			{
				if (bc->getColumnSize(j)==1)
					op->vhdl << tab << op->declare(join("in_concat_", compressorIndex,"_", inConcatIndex), bc->getColumnSize(j)) 
					         << "(0) <= " << signal[j].str() << ";" << endl;
				else
					op->vhdl << tab << op->declare(join("in_concat_", compressorIndex,"_", inConcatIndex), bc->getColumnSize(j)) 
					         << " <= " << signal[j].str() << ";" << endl;
				op->inPortMap(bc, join("X",j), join("in_concat_", compressorIndex,"_", inConcatIndex));
				++inConcatIndex;
			}

		op->outPortMap(bc, "R", join("out_concat_", compressorIndex,"_", outConcatIndex));
    
		op->vhdl << tab << op->instance(bc, join("Compressor_",compressorIndex));

		removeCompressedBits(i,bc->getColumnSize(0));
		if(bc->getColumnSize(1)!=0)
			removeCompressedBits(i+1,bc->getColumnSize(1));
		
		// add the bits, at the current (global) instant.
		addBit(i, join("out_concat_", compressorIndex,"_", outConcatIndex, "(0)"),"");
		addBit(i+1,	join("out_concat_", compressorIndex,"_", outConcatIndex, "(1)"),"");
		if(!((bc->getColumnSize(0)==3) && (bc->getColumnSize(1)==0)))
			addBit(i+2, join("out_concat_", compressorIndex,"_", outConcatIndex, "(2)"),"");


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
		list<WeightedBit*>::iterator it = wb.begin();
		
		if(op->getTarget()->isPipelined())
			{
				for (list<WeightedBit*>::iterator iter = wb.begin(), end = wb.end(); iter != end; ++iter)
					{
						if ((*iter)->getCycle()<=cycle)
							{
								bitsToCompress++;
							}
					}
			}
		else
			{
				for (list<WeightedBit*>::iterator iter = wb.begin(), end = wb.end(); iter != end; ++iter)
					{
						if (  ((*iter)->getCycle()<cycle) ||  
						      (((*iter)->getCycle()==cycle)  &&  ((*iter)->getCriticalPath(cycle) <= (*it)->getCriticalPath(cycle))) )
							{
								bitsToCompress++;
							}
					}
			
				if(wb.size()>2)
					{
						double delta = op->getTarget()->lutDelay();
						while(bitsToCompress<=2)
							{
								bitsToCompress=0;
								it = wb.begin();
								for (list<WeightedBit*>::iterator iter = wb.begin(), end = wb.end(); iter != end; ++iter)
									{
										if (  ((*iter)->getCycle()<cycle) ||  
										      (((*iter)->getCycle()==cycle)  &&  ((*iter)->getCriticalPath(cycle) <= (*it)->getCriticalPath(cycle) + delta)) )
											{
												bitsToCompress++;
											}
									}
								delta += op->getTarget()->lutDelay();
							}
					}
			}
		
		return bitsToCompress;
	}

	
	void BitHeap::removeCompressedBits(int c, int red)
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
					}
        
	}

	
	void BitHeap::printColumnInfo(int w)
	{
		int i=0;
		
		for(list<WeightedBit*>::iterator it = bits[w].begin(); it!=bits[w].end(); ++it)
			{
				REPORT(DEBUG, "element "<<i<<" cycle = "<<(*it)->getCycle() << " and cp = "<<(*it)->getCriticalPath((*it)->getCycle()));
				i++; 
			}
	}



	void BitHeap::generateCompressorVHDL()
	{
	
		generatePossibleCompressors();

		op->vhdl << tab << endl << "-- Beginning of code generated by BitHeap::generateCompressorVHDL" << endl;
		minWeight=0;

		while (getMaxHeight()>2)
			{
				compress();
			}
         
		REPORT(DEBUG, "Column height after all compressions");
		for (int w=0; w<bits.size(); w++) 
			{
				REPORT(DEBUG, "   w=" << w << ":\t height=" << bits[w].size());
				printColumnInfo(w);  
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
	



	BitHeap::WeightedBit* BitHeap::getFinalLatestBit()
	{
	
		double maxCycle=0;
		double maxCP = 0;
		WeightedBit **b;
		for(int w=minWeight; w<maxWeight; w++)
			{
			
			
				if (bits[w].size() >= 1)
					{
						for(list<WeightedBit*>::iterator it = bits[w].begin(); it!=bits[w].end(); ++it)
							{
								if (maxCycle < (*it)->getCycle())
									{
										maxCycle = (*it)->getCycle();
										maxCP = (*it)->getCriticalPath(maxCycle);
										b = &*it;
									}
								else
									if ((maxCycle == (*it)->getCycle())  &&  (maxCP < (*it)->getCriticalPath(maxCycle)))
										{
											maxCP = (*it)->getCriticalPath(maxCycle);
											b = &*it;
										}
							}
					}
			}	
		
		return *b;
	}

	void BitHeap::adderVHDL()
	{


		stringstream inAdder0, inAdder1, outAdder;

		unsigned i=maxWeight-1;
	
		while((i>=minWeight)&&(i<maxWeight))
			{

				if(i>=0)
					{  
						list<WeightedBit*>::iterator it = bits[i].begin();
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

				--i;
			}

		

		inAdder0 << ";";
		inAdder1 << ";";
		
		WeightedBit* b = getFinalLatestBit();
		
		op->setCycle(  b ->getCycle()  );
		op->setCriticalPath(  b->getCriticalPath(op->getCurrentCycle()));
		op->manageCriticalPath(b->getCriticalPath(op->getCurrentCycle()) + op->getTarget()->adderDelay(maxWeight-minWeight));

		op->vhdl << tab << op->declare("inAdder0", maxWeight-minWeight) << "<= " << inAdder0.str() << endl;
		op->vhdl << tab << op->declare("inAdder1", maxWeight-minWeight) << " <= " << inAdder1.str() << endl;

		op->vhdl << tab << op->declare("outAdder", maxWeight-minWeight+1) << " <= ('0' & inAdder0) + ('0' & inAdder1);" << endl;





		op->vhdl << tab << "-- concatenate all the compressed chunks" << endl;


		op->vhdl << tab << op->declare("CompressionResult", (maxWeight+1)) << " <= outAdder ";



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
		


		REPORT(DEBUG, "maxWeight="<< maxWeight);
		REPORT(DEBUG, "Column height before compression");
		for (w=0; w<bits.size(); w++) {
			REPORT(DEBUG, "   w=" << w << ":\t height=" << bits[w].size());
			printColumnInfo(w); 
		}		

		
		REPORT(DEBUG, "Checking whether rightmost columns need compressing");
		bool alreadyCompressed=true;
		w=minWeight;
		while(w<bits.size() && alreadyCompressed) 
			{
			
			
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
		
		


		if(w!=minWeight)
			{
				if (w-minWeight>1)
					op->vhdl << tab << op->declare(join("tempR", chunkDoneIndex), w-minWeight, true) << " <= " ;
				else
					op->vhdl << tab << op->declare(join("tempR", chunkDoneIndex), 1, false) << " <= " ;
				unsigned i=w-1;
				while((i>=minWeight)&&(i<w)) 
					{
						REPORT(DEBUG,"crash "<<i);
						if(currentHeight(i)==1) 
							{
								op->vhdl << (bits[i].front())->getName();
								bits[i].pop_front();
							}
						else
							{
								op->vhdl << "'0'";
							}
						if (i!=minWeight)
							op->vhdl << " & ";	
						i--;
					}
				op->vhdl << "; -- already compressed" << endl; 
				chunkDoneIndex++;			
			}
		
		
		
		minWeight=w;

		REPORT(DEBUG, "minWeight="<< minWeight);
		
		int cnt[maxWeight-minWeight];
		
		for(int i=minWeight; i<maxWeight; i++)
			{
				cnt[i] = bits[i].size();
			}
		
		
		REPORT(DEBUG,"start compressing");
		//Remaining structure must be compressed
		unsigned j;

		// i is the column index
		for(unsigned i=minWeight; i<maxWeight; i++)
			{

				j=0; // index of the compressors among the possible ones
			

				// map the best compressor as many times as possible 
				while(cnt[i] >= possibleCompressors[j]->getColumnSize(0))
					{
				
				
						//compressors.push_back(possibleCompressors[j]) ;
						REPORT(DEBUG,endl);
						REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
						cnt[i]-=possibleCompressors[j]->getColumnSize(0);
						elemReduce(i,possibleCompressors[j]);
						usedCompressors[j]=true;


				
					}

				++j;

				while(   (j < possibleCompressors.size())   &&   ( cnt[i] > 2 )   )
					{

						if(   cnt[i]  >=  possibleCompressors[j]->getColumnSize(0)   )
							{
								if(possibleCompressors[j]->getColumnSize(1)!=0)
									{

										if((cnt[i+1]>=possibleCompressors[j]->getColumnSize(1))&&(i<bits.size()-1))
											{
												REPORT(DEBUG,endl);
												REPORT(DEBUG,"Using Compressor " << j <<", reduce columns "<<i<<" and "<<i+1);
												cnt[i]-=possibleCompressors[j]->getColumnSize(0);
												cnt[i+1]-=possibleCompressors[j]->getColumnSize(1);
												elemReduce(i,possibleCompressors[j]);
							
												usedCompressors[j]=true;
											}
										else 
											++j;
									}
								else
									{
										REPORT(DEBUG,endl);
										REPORT(DEBUG,"Using Compressor " << j <<", reduce column "<<i);
										elemReduce(i,possibleCompressors[j]);
										cnt[i]-=possibleCompressors[j]->getColumnSize(0);
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

