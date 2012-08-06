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
#include <vector>
#include <list>

using namespace std;


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
		p  << "heap" << bh->getGUid() << "_w" << weight << "_" << bh->newUid(weight);
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

	int BitHeap::WeightedBit::computeStage(int stagesPerCycle, double elementaryTime)
	{
		return (getCycle()*stagesPerCycle + getCriticalPath(getCycle())/elementaryTime);        
	}


	BitHeap::BitHeap(Operator* op, int maxWeight) :
		op(op), maxWeight(maxWeight)
	{
		// Set up the vector of lists of weighted bits, and the vector of uids
		srcFileName=op->getSrcFileName() + ":BitHeap"; // for REPORT to work
		uniqueName_=op->getName() + ":BitHeap"; // for REPORT to work
		guid = Operator::getNewUId();
		REPORT(DEBUG, "Creating BitHeap of size " << maxWeight);
		chunkDoneIndex=0;
		inConcatIndex=0;
		outConcatIndex=0;
		compressorIndex=0;
		adderIndex=0;
		for(int i=0; i<10;++i)
			usedCompressors[i]=false;
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


	int BitHeap::newUid(unsigned w){
		return uid[w]++;
	}	

	int BitHeap::getGUid(){
		return guid;
	}	



	BitHeap::WeightedBit* BitHeap::computeLatest(unsigned w, int c0, int c1)
	{

		if(w>=maxWeight)
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



	void BitHeap::addUnsignedBitVector(unsigned weight, string x, unsigned size){
#if 0 // Should be a warning
		if(weight+size>maxWeight)
			THROWERROR("in addUnsignedBitVector: Size of signal " << x << " is " << size << 
			           ", adding it at weight " << weight << " overflows the bit heap (maxWeight=" << maxWeight << ")");
#endif
		
		op->setCycleFromSignal(x);
		for (unsigned i=0; i<size; i++) {
			ostringstream s;
			s << x << of(i);
			addBit(weight+i, s.str());
		}
	}

	void BitHeap::subtractUnsignedBitVector(unsigned weight, string x, unsigned size){
#if 0 // Should be a warning
		if(weight+size>maxWeight)
			THROWERROR("in subtractUnsignedBitVector: Size of signal " << x << " is " << size << 
			           ", adding it at weight " << weight << " overflows the bit heap (maxWeight=" << maxWeight << ")");
#endif

		for (unsigned i=0; i<size; i++) {
			ostringstream s;
			s << "not " << x << of(i);
			addBit(weight+i, s.str());
		}
		addConstantOneBit(weight);
		if (size+weight < maxWeight) {
			// complement all the way to maxWeight
			for (unsigned i=size+weight; i<maxWeight; i++) {
				addConstantOneBit(i);
			}
		}
	}


	void BitHeap::addSignedBitVector(unsigned weight, string x, unsigned size){
		//TODO
	}


	void BitHeap::subtractSignedBitVector(unsigned weight, string x, unsigned size){
		// TODO
	}


	void  BitHeap::addDSP(MultiplierBlock* m)
	{
		
		//now, I can insert in any position, because the supertile chaining will be done in an ascending way(hopefully)
		mulBlocks.push_back(m);

		
		//this is the old code, inserting increasingly by weight - maybe not needed anymore
		/*
		//insert to proper position, the vector of multiplierBlocks should be ascending by weigth
		unsigned i=0;
		unsigned size=mulBlocks.size();
		

		while((i<size)&&(mulBlocks[i]<=m))
		{	
			REPORT(DETAILED,"size "<<size<<" i "<<i);
			++i;
		}
		
		mulBlocks.insert(mulBlocks.begin()+i,m);
		*/

	}



	void BitHeap::buildSupertiles()
	{
		for(unsigned i=0;i<mulBlocks.size();i++)
			for(unsigned j=0;j<mulBlocks.size();j++)
			{ 
				//if 2 blocks can be chained, the the chaining is done ascending by weight.
				//TODO improve the chaining		
			

				if(mulBlocks[i]->canBeChained(mulBlocks[j]))
				{
					
					if(mulBlocks[j]->getWeight()<=mulBlocks[i]->getWeight())
						if((mulBlocks[j]->getNext()==NULL)&&(mulBlocks[i]->getPrevious()==NULL))
						{	REPORT(DETAILED,"mulblocks[j]= "<<mulBlocks[j]->getWeight()<<" mulBlock[i]= "<<mulBlocks[i]->getWeight());
							REPORT(DETAILED,"j<=i");
							mulBlocks[j]->setNext(mulBlocks[i]);
							mulBlocks[i]->setPrevious(mulBlocks[j]);
						}	
											
					else
						if((mulBlocks[i]->getNext()==NULL)&&(mulBlocks[j]->getPrevious()==NULL))	
						{	REPORT(DETAILED,"mulblocks[j]= "<<mulBlocks[j]->getWeight()<<" mulBlock[i]= "<<mulBlocks[i]->getWeight());
							REPORT(DETAILED,"j>i");
							mulBlocks[i]->setNext(mulBlocks[j]);
							mulBlocks[j]->setPrevious(mulBlocks[i]);
						}			
				}

			}


		//just for debugging
		for(unsigned i=0;i<mulBlocks.size();i++)
		{
			if(mulBlocks[i]->getNext()!=NULL)
			REPORT(DETAILED, mulBlocks[i]->getWeight()<<" chained with "<<mulBlocks[i]->getNext()->getWeight());
		}
	}



	void BitHeap::generateSupertileVHDL()
	{
		//making all the possible supertiles
		buildSupertiles();


		//generate the VHDL code for each supertile
		for(unsigned i=0;i<mulBlocks.size();i++)
		{	
			//take just the blocks which are roots
			if(mulBlocks[i]->getPrevious()==NULL)
			{
				int uid=0;
				MultiplierBlock* next;
				MultiplierBlock* current=mulBlocks[i];
				int newLength=0;
			
				//the first DSP from the supertile(it has the smallest weight in the supertile)
				generateVHDLforDSP(current,uid,i);

				//iterate on the of the supertile
				while(current->getNext()!=NULL)
				{	
					uid++;
					next=current->getNext();
					newLength=current->getSigLength();					
					generateVHDLforDSP(next,uid,i);
					//TODO ! replace 17 with multiplierblock->getshift~ something like that
					
					//******pipeline*******//	
					op->setCycleFromSignal(next->getSigName());
					op->syncCycleFromSignal(current->getSigName());
					op->manageCriticalPath(  op->getTarget()->DSPAdderDelay() ) ; 

					//addition, the 17lsb-s from the first block will go directly to bitheap
					op->vhdl << tab <<	op->declare(join("DSPch",i,"_",uid),newLength)<< "<= " <<next->getSigName() 
						<< " +  ("<< zg(17)<<" & "<<current->getSigName()<<range(newLength-1,17)<<" );"<<endl ; 


					//sending the 17 lsb to the bitheap
           			for(int k=16;k>=0;k--)
					{
						int weight=current->getWeight()+k;	
						if(weight>=0)
						{	
							stringstream s;
							s<<current->getSigName()<<"("<<k<<")";
							addBit(weight,s.str());
						}
					}


					//setting the name and length of the current block, to be used properly in the next iteration
					stringstream q;
					q<<join("DSPch",i,"_",uid);
					next->setSignalName(q.str());		
					next->setSignalLength(newLength);
					//next
					current=next;
				}
		
				// adding the result to the bitheap (in the last block from the supertile)
					string name=current->getSigName();
				int length=current->getSigLength();
				for(int k=length-1;k>=0;k--)
				{
					int weight=current->getWeight()+k;	
					REPORT(DETAILED,"k= "<<k <<" weight= "<<weight);				
					if(weight>=0)
					{	
						stringstream s;
						s<<name<<"("<<k<<")";
						addBit(weight,s.str());
					}
				}
			}
		}
	}


	void  BitHeap::addBit(unsigned w, string rhs, string comment)
	{
			
		// ignore bits beyond the declared maxWeight
		if(w >= maxWeight)
			return; 

		WeightedBit* bit= new WeightedBit(this, w) ;
		// created at (op->getCycle(), opt-getCriticalPath())

		list<WeightedBit*>& l=bits[w];

		//insert so that the list is sorted by bit cycle/delay
		list<WeightedBit*>::iterator it=l.begin(); 
		bool proceed=true;
		unsigned count=0;
		while(proceed) 
			{
				if (it==l.end() || (*bit <= **it))
					{ // test in this order to avoid segfault!

						l.insert(it, bit);
						proceed=false;
					}
				else 
					{
						count++;
						it++;
					}
			}

        //not sure if we need this anymore; also, if used, will kill new version of bitheap
        //
		//the cnt[w] should increase, because the bit introduced is smaller then the biggest bit which will be compressed in current iteration
		//if(count<cnt[w])
		//	cnt[w]++;

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

  
	void BitHeap::elemReduce(unsigned i, BasicCompressor* bc)
	{

		list<WeightedBit*>::iterator it = bits[i].begin();
		stringstream signal[2];

		op->vhdl << endl;

		
		WeightedBit* b =  computeLatest(i, bc->getColumnSize(0), bc->getColumnSize(1)) ; 
		if(b)
			{
				op->setCycle(  b ->getCycle()  );
				op->setCriticalPath(  b ->getCriticalPath(op->getCurrentCycle()));
				op->manageCriticalPath( op->getTarget()->lutDelay() + op->getTarget()->localWireDelay() );
				if(b->getCycle() < op->getCurrentCycle())
					{
						drawCycleLine = true;
					}
			}

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



	unsigned BitHeap::currentHeight(unsigned w) {
		int h=0;
		list<WeightedBit*>& l = bits[w];
		h=l.size();
		return h;
	}
	

	//c - the current column, red- number of bits to be reduced
	void BitHeap::removeCompressedBits(int c, int red)
	{
		while(red>0)
			{
				removeBit(c,0);
				red--;
			}
	}


	unsigned BitHeap::getMaxHeight()

	{
		unsigned max=0; 
		for(unsigned i=0; i<maxWeight; i++)
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
				//REPORT(DEBUG, "element "<<i<<" cycle = "<<(*it)->getCycle() << " and cp = "<<(*it)->getCriticalPath((*it)->getCycle()));
				i++; 
			}
	}



	void BitHeap::generateCompressorVHDL()
	{
		op->vhdl << tab << endl << "-- Beginning of code generated by BitHeap::generateCompressorVHDL" << endl;

		// add the constant bits to the actual bit heap 
		op->setCycle(0); 
		op->setCriticalPath(0.0);
		op->vhdl << tab << endl << "-- Adding the constant bits" << endl;
		for (unsigned w=0; w<maxWeight; w++)
			if (1 == ((constantBits>>w) & 1) )
				addBit(w, "'1'");

		
		generateSupertileVHDL();
	

		
		initializeDrawing();
		generatePossibleCompressors();
        
		WeightedBit* firstBit = getFirstSoonestBit();
        
		int minCycle = firstBit->getCycle();
		double minCP = firstBit->getCriticalPath(minCycle);

		elementaryTime = op->getTarget()->lutDelay() + op->getTarget()->localWireDelay();

		stagesPerCycle = (1/op->getTarget()->frequency()) / elementaryTime;

		int stage = minCycle*stagesPerCycle + minCP/elementaryTime ;

		REPORT(DETAILED, "minimum stage " << stage);

		minWeight=0;

		int offsetY=0;

		didCompress = true;

		//compressing until the maximum height of the columns is 3
		//FIXME change 2 to 3 in order to use extra additions 
		while (getMaxHeight()>3)
			{
                
				REPORT(DETAILED, "didCompress " << didCompress);
				if(didCompress)
					{
						offsetY += 20 + getMaxHeight()*10;
						drawConfiguration(offsetY);
					}
				compress(stage);
				stage++;
			}
		offsetY += 20 + getMaxHeight()*10;
		drawConfiguration(offsetY);

        REPORT(DEBUG, endl);
        REPORT(DEBUG, "only three levels left");
        
		//do additional compressions or additions
		if (getMaxHeight()>2)
			{
				unsigned i = 0;
				// remember the initial heights
				for(unsigned i=0; i<maxWeight; i++)	{
						cnt[i]=bits[i].size();
					}

				//find the first column with 3 bits (starting from LSB); the rest go to the final adder
				while(cnt[i]<3)
						i++;

    	    	REPORT(DEBUG, "Column height after all compressions");
	        	for (unsigned w=0; w<bits.size(); w++) 
		    	{
			    	REPORT(DEBUG, "   w=" << w << ":\t height=" << bits[w].size());
			    	printColumnInfo(w);  
		    	}

				int first2, first3;
				first3=i;
				bool doLoop=true;

				WeightedBit* latestBit;
				while(doLoop)
					{
						// Now we are sure cnt[i] is 3;
						//look for the first two
                        REPORT(DEBUG, i << "  " << cnt[i]);
						while (i<maxWeight && cnt[i]==3) 
							i++;
                        
                        REPORT(DEBUG, i << "  " << cnt[i]);
						
                        if (i==maxWeight) {
							applyCompressor3_2(maxWeight-1);
							doLoop=false;
						}
						else {
							// Now we are sure there is at least one 2
							first2=i;
						
							if(first3<first2-1) {
								for(int j=first3; j<first2-1; j++)
									applyCompressor3_2(j);
							}
							
							//look for the first three again
							while (i<maxWeight &&  cnt[i]<3) 
								i++;
							
							if (i==maxWeight) 
							{	
								latestBit = getLatestBit(first2-1, maxWeight-1); 
								if(latestBit)
								{
									op->setCycle(  latestBit ->getCycle()  );
									op->setCriticalPath(   latestBit ->getCriticalPath(op->getCurrentCycle()));
									op->manageCriticalPath( op->getTarget()->localWireDelay() + op->getTarget()->adderDelay(maxWeight-first2+1) );
									if( latestBit ->getCycle() < op->getCurrentCycle())
									{
										drawCycleLine = true;
									}
								}


								applyAdder(first2-1,maxWeight-1);							
								doLoop=false;
							}
							else 
							{
								first3=i;
								latestBit = getLatestBit(first2-1, first3-1); 
								if(latestBit)
								{
									op->setCycle(   latestBit  ->getCycle()  );
									op->setCriticalPath( latestBit->getCriticalPath(op->getCurrentCycle()));
									op->manageCriticalPath(op->getTarget()->localWireDelay() + op->getTarget()->adderDelay(first3 - first2 + 1));
									if(latestBit->getCycle() < op->getCurrentCycle())
									{
										drawCycleLine = true;
									}
								}


								applyAdder(first2-1, first3-1);	
                                //i++;                        
							}
						}
						

#if 0
						//REPORT(DEBUG,"print i "<<i);
						if ((cnt[i]==3) && (cnt[i+1]==3))
							{
								applyCompressor3_2(i);
								i++;
								//REPORT(DEBUG, "crash here");
							}
						else
							{
								unsigned j = i+1; 
								while((j<maxWeight) && (cnt[j]<3))
									{
										j++;

										REPORT(DEBUG, j);
									}

								j--;
								applyAdder(i,j);
								j++;
								i=j;
							}
#endif

						}
					}
         
		REPORT(DEBUG, "Column height after all compressions");
		for (unsigned w=0; w<bits.size(); w++) 
			{
				REPORT(DEBUG, "   w=" << w << ":\t height=" << bits[w].size());
				printColumnInfo(w);  
			}	

       
		//checking the used compressors 
		for(int i=0;i<10;i++)
			{
				if (usedCompressors[i]==true)
					{	
						op->getOpListR().push_back(possibleCompressors[i]);
					}
			}
	    
		
		offsetY += 20 + getMaxHeight() * 10;
		drawConfiguration(offsetY);
		closeDrawing(offsetY);
		
		//final addition
		generateFinalAddVHDL();
		
	}


	//returns a pointer to the first bit which has the smallest cycle & CP combination
	BitHeap::WeightedBit* BitHeap::getFirstSoonestBit()
	{
		int minCycle= 9999;
		double minCP = 100;
		WeightedBit *b=0;
		for(unsigned w=minWeight; w<maxWeight; w++)
			{
				if(bits[w].size()==0)
					{

					}
				else{			
					if (bits[w].size() >= 1)
						{
							for(list<WeightedBit*>::iterator it = bits[w].begin(); it!=bits[w].end(); ++it)
								{
									if (minCycle > (*it)->getCycle())
										{
											minCycle = (*it)->getCycle();
											minCP = (*it)->getCriticalPath(minCycle);
											b = *it;
										}
									else
										if ((minCycle == (*it)->getCycle())  &&  (minCP > (*it)->getCriticalPath(minCycle)))
											{
												minCP = (*it)->getCriticalPath(minCycle);
												b = *it;
											}
								}
						}
				}
			}	
		
		return b;
      
	}

	void BitHeap::applyCompressor3_2(int col)
	{
		for(unsigned i=0; i<possibleCompressors.size(); i++)
			{
				if((possibleCompressors[i]->getColumnSize(0)==3) && (possibleCompressors[i]->getColumnSize(1)==0))
					{
						REPORT(DEBUG, endl);
						REPORT(DEBUG, "Using Compressor3_2, reduce column " << col);
						elemReduce(col, possibleCompressors[i]);
						usedCompressors[i]=true;
					}
			}
	}




	// addition stretching from weights col0(LSB) to msb included
	// assumes cnt has been set up
	void BitHeap::applyAdder(int lsb, int msb)
	{
		stringstream inAdder0, inAdder1, cin, outAdder;

	    REPORT(DEBUG, "Applying an adder from columns " << lsb << " to " << msb);					
		for(int i = msb; i>=lsb+1; i--)
			{
				list<WeightedBit*>::iterator it = bits[i].begin();
				
				if(cnt[i]>=2)
					{
						inAdder0 << (*it)->getName();
						it++;
						inAdder1 << (*it)->getName();
					}
				else
					{
						//if(bits[i].size()==1)
						if(cnt[i]==1)
							{
								inAdder0 << (*it)->getName();
								inAdder1 << "\'0\'";
							}
						else
							{
								inAdder0 << "\'0\'";                       
								inAdder1 << "\'0\'";
							}
					}

				inAdder0 << " & ";
				inAdder1 << " & ";

				removeCompressedBits(i, cnt[i]);
			}

		// We know the LSB col is of size 3
		list<WeightedBit*>::iterator it = bits[lsb].begin();
		inAdder0 << (*it)->getName();
		it++;
		inAdder1 << (*it)->getName();
		it++;
		cin << (*it)->getName() << ";";
		removeCompressedBits(lsb, cnt[lsb]);


		inAdder0 << ";";
		inAdder1 << ";";

		op->vhdl << tab << op->declare(join("inAdder0_",adderIndex), msb-lsb+2) << " <= \'0\' & " << inAdder0.str() << endl;
		op->vhdl << tab << op->declare(join("inAdder1_",adderIndex), msb-lsb+2) << " <= \'0\' & " << inAdder1.str() << endl;
		op->vhdl << tab << op->declare(join("cin_",adderIndex)) << " <= " << cin.str() << endl;

#if 0

		IntAdder* adder = new IntAdder(op->getTarget(), msb-lsb+2);
		op->getOpListR().push_back(adder);
		
		op->inPortMap(adder, "X", join("inAdder0_",adderIndex));
		op->inPortMap(adder, "Y", join("inAdder1_",adderIndex));
		op->inPortMap(adder, "Cin", join("cin_",adderIndex));
		
		op->outPortMap(adder, "R", join("outAdder_",adderIndex));

		op->vhdl << tab << op->instance(adder, join("Adder_",adderIndex));
				
#else
		op->vhdl << tab << op->declare(join("outAdder_",adderIndex),msb-lsb+2) << " <= " 
                        <<  join("inAdder0_",adderIndex) << " + " << join("inAdder1_",adderIndex) << " + " << join("cin_",adderIndex) << ";" <<endl ;

#endif

		for(int i=lsb; i<msb + 2 ; i++)	{
			addBit(i, join("outAdder_", adderIndex,"(",i-lsb,")"),"");
		}
 

		adderIndex++;

	}


	//returns a pointer to the bit which has the biggest criticalPath+cycle combination, considering the given column limits (both included)
	BitHeap::WeightedBit* BitHeap::getLatestBit(int lsbColumn, int msbColumn)
	{
	
		double maxCycle=0;
		double maxCP = 0;
		WeightedBit *b=0;
		for(int w=lsbColumn; w<=msbColumn; w++)
			{
				if(bits[w].size()==0)
					{

					}
				else{			
					if (bits[w].size() >= 1)
						{
							for(list<WeightedBit*>::iterator it = bits[w].begin(); it!=bits[w].end(); ++it)
								{
									if (maxCycle < (*it)->getCycle())
										{
											maxCycle = (*it)->getCycle();
											maxCP = (*it)->getCriticalPath(maxCycle);
											b = *it;
										}
									else
										if ((maxCycle == (*it)->getCycle())  &&  (maxCP <= (*it)->getCriticalPath(maxCycle)))
											{
												maxCP = (*it)->getCriticalPath(maxCycle);
												b = *it;
											}
								}
						}
				}
			}	
		
		return b;
	}
	


	//the final addition
	void BitHeap::generateFinalAddVHDL()
	{
		stringstream inAdder0, inAdder1, outAdder;

		unsigned i=maxWeight-1;
	
		//forming the input signals for the first and second line
		while((i>=minWeight)&&(i<maxWeight))
			{
				REPORT(DEBUG,"i=   "<<i);
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
								if (bits[i].size()==1)
									{
										inAdder0 << (*it)->getName();
										inAdder1 << "\'0\'";
									}
								else
									{
										inAdder0 << "\'0\'";
										inAdder1 << "\'0\'";
									}
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
		
		WeightedBit* b = getLatestBit(minWeight, maxWeight-1);
		//managing the pipeline
		if(b)
			{
				op->setCycle(  b ->getCycle()  );
				op->setCriticalPath(  b->getCriticalPath(op->getCurrentCycle()));
				op->manageCriticalPath(op->getTarget()->localWireDelay() + op->getTarget()->adderDelay(maxWeight-minWeight));
				if(b->getCycle() < op->getCurrentCycle())
					{
						drawCycleLine = true;
					}
			}

		op->vhdl << tab << op->declare(join("inAdder0_", guid), maxWeight-minWeight) << "<= " << inAdder0.str() << endl;
		op->vhdl << tab << op->declare(join("inAdder1_", guid), maxWeight-minWeight) << " <= " << inAdder1.str() << endl;
		op->vhdl << tab << op->declare(join("outAdder_", guid), maxWeight-minWeight+1) 
		         << " <= ('0' & "<< join("inAdder0_", guid) << ") + ('0' & " << join("inAdder1_", guid) << ");" << endl;
		op->vhdl << tab << "-- concatenate all the compressed chunks" << endl;
		//result		
		op->vhdl << tab << op->declare(join("CompressionResult", guid), (maxWeight+1)) << " <= " << join("outAdder", guid);


		//adding the rightmost bits
		for(int i=chunkDoneIndex-1; i>=0; i--)
			op->vhdl <<  " & " << join("tempR", i, "_", guid);

		op->vhdl << ";" << endl;
		op->vhdl << tab << "-- End of code generated by BitHeap::generateCompressorVHDL" << endl;

	}


	string BitHeap::getSumName() {
		return join("CompressionResult", guid);
	}




	//the compression
	void BitHeap::compress(int stage)
	{
		Target* target=op->getTarget();

		
		
		unsigned n = target->lutInputs();
		unsigned w;
		didCompress = false;
		
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
					op->vhdl << tab << op->declare(join("tempR", chunkDoneIndex, "_", guid), w-minWeight, true) << " <= " ;
				else
					op->vhdl << tab << op->declare(join("tempR", chunkDoneIndex, "_", guid), 1, false) << " <= " ;
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
		
		
		
		for(unsigned i=minWeight; i<maxWeight; i++)
			{
				cnt[i]=0;
				for(list<WeightedBit*>::iterator it = bits[i].begin(); it!=bits[i].end(); it++)
					{
						if((*it)->computeStage(stagesPerCycle, elementaryTime)<=stage)
							{
								cnt[i]++;
							}
					}
				//cnt[i]=bits[i].size();
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
						didCompress = true;


				
					}

				++j;
				//search for the next best compressor which fits the remaining bits
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
												didCompress = true;
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
										didCompress = true;
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


	void BitHeap::initializeDrawing()
	{
		ostringstream figureFileName;
		figureFileName << "bit_heap.svg";
	
		drawCycleNumber = 0;

		FILE* pfile;
		pfile  = fopen(figureFileName.str().c_str(), "w");
		fclose(pfile);
		
		fileFig.open (figureFileName.str().c_str(), ios::trunc);

		fileFig << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
		fileFig << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
		fileFig << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
	}

	void BitHeap::drawConfiguration(int offsetY)
	{
		int turnaroundX = 1500;
		int color = 0;
		int tempCycle = 0;
		int cnt = 0;
		double tempCP = 0;
		if((drawCycleLine) && (drawCycleNumber==0))
			{
				fig << "<text x=\"" << turnaroundX + 85 << "\" y=\"" << 40
				    << "\" fill=\"midnightblue\">" << "Cycle"<< "</text>" << endl;
			}

		if(drawCycleLine)
			{
				drawCycleNumber++;
				fig << "<line x1=\"" << turnaroundX + 200 << "\" y1=\"" << offsetY +10 << "\" x2=\"" << turnaroundX - bits.size()*10 - 50
				    << "\" y2=\"" << offsetY +10 << "\" style=\"stroke:midnightblue;stroke-width:2\" />" << endl;

				fig << "<text x=\"" << turnaroundX + 100 << "\" y=\"" << offsetY + 3
				    << "\" fill=\"midnightblue\">" << drawCycleNumber - 1 << "</text>" << endl;

				fig << "<text x=\"" << turnaroundX + 100 << "\" y=\"" << offsetY + 27
				    << "\" fill=\"midnightblue\">" << drawCycleNumber  << "</text>" << endl;

				drawCycleLine = false;
			}
		else
			{
				fig << "<line x1=\"" << turnaroundX + 200 << "\" y1=\"" << offsetY +10 << "\" x2=\"" << turnaroundX - bits.size()*10 - 50
				    << "\" y2=\"" << offsetY +10 << "\" style=\"stroke:lightsteelblue;stroke-width:1\" />" << endl;
				drawCycleLine = false;
			}

		for(unsigned i=0; i<bits.size(); i++)
			{

				if(bits[i].size()>0)
					{
						color=0;
						tempCycle = 0;
						tempCP = 0;
						cnt = 0;
						for(list<WeightedBit*>::iterator it = bits[i].begin(); it!=bits[i].end(); ++it)
							{
								if(it==bits[i].begin())
									{
										tempCycle = (*it)->getCycle();
										tempCP = (*it)->getCriticalPath(tempCycle);
									}
								else
									{
										if((tempCycle!=(*it)->getCycle()) || 
										   ((tempCycle==(*it)->getCycle()) && (tempCP!=(*it)->getCriticalPath((*it)->getCycle()))))
											{
												tempCycle = (*it)->getCycle();
												tempCP = (*it)->getCriticalPath(tempCycle);
												color++;
											}
									}

								drawBit(cnt, i, turnaroundX, offsetY, (*it)->computeStage(stagesPerCycle, elementaryTime));
								cnt++;
							}
					}
			}

        
	}

	void BitHeap::drawBit(int cnt, int w, int turnaroundX, int offsetY, int c)
	{
		const std::string colors[] = { "#97bf04","#0f1af2", "#f5515c", "#3958ff","#f2eb8d", "indianred", "yellow", "lightgreen"};

		int index = c % 8;

		fig << "<circle cx=\"" << turnaroundX - w*10 - 5 << "\" cy=\"" << offsetY - cnt*10 - 5 << "\" r=\"3\" fill=\"" << colors[index] << "\"/>" << endl;
	}

	void BitHeap::closeDrawing(int offsetY)
	{
		int turnaroundX = 1500;
		fig << "<line x1=\"" << turnaroundX + 50 << "\" y1=\"" << 20 << "\" x2=\"" << turnaroundX + 50
		    << "\" y2=\"" << offsetY +30 << "\" style=\"stroke:midnightblue;stroke-width:1\" />" << endl;

		fig << "</svg>" << endl;

		fileFig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"100%\" height=\"100%\" viewBox=\"" <<turnaroundX - bits.size()*10 - 80 
		        << " " << 100 << " " << turnaroundX + 50 << " " << offsetY + 20 <<  "\">" << endl; 
		//    fileFig << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"0 0 500 500\">" << endl; 

		fileFig << fig.str();
       
               

		fileFig.close();
	}

	void BitHeap::generateVHDLforDSP(MultiplierBlock* m, int uid,int i)
	{
		REPORT(DETAILED,"dsp");
			
		stringstream s;
		int topX=m->gettopX();
		int topY=m->gettopY();
		int botX=topX+m->getwX()-1;
		int botY=topY+m->getwY()-1;
		string input1=m->getInputName1();
		string input2=m->getInputName2();
		int zerosX;
		int zerosY;
		if(m->getwX()>m->getwY())
		{
			zerosX=25-m->getwX();
			zerosY=18-m->getwY();
		}		
		else 
		{		
			zerosX=25-m->getwY();
			zerosY=18-m->getwX();
		}	


		//int zerosY=18-m->getwY();
		
		//PIPELINE!!!!
		op->setCycleFromSignal(m->getInputName1());
		op->syncCycleFromSignal(m->getInputName2());
		op->manageCriticalPath(  op->getTarget()->DSPMultiplierDelay() ) ;
	//	op->setCycle(uid);
		if(uid==0)	
			op->vhdl << tab << op->declare(join("DSPch",i,"_",uid), m->getwX()+m->getwY()+zerosX+zerosY) 
				<< " <= (" <<zg(zerosX)<<" & " << input1<<range(botX,topX)<<") * (" <<zg(zerosY) <<" & "
			    << input2 <<range(botY,topY)<<");"<<endl;
		else
			op->vhdl << tab << op->declare(join("DSP",i,"_",uid), m->getwX()+m->getwY()+zerosX+zerosY) 
				<< " <= (" <<zg(zerosX)<<" & " << input1<<range(botX,topX)<<") * (" <<zg(zerosY) <<" & "
			    << input2 <<range(botY,topY)<<");"<<endl;

	
		if(uid==0)
			s<<join("DSPch",i,"_",uid);
		else
			
			s<<join("DSP",i,"_",uid);

		m->setSignalName(s.str());
		m->setSignalLength(m->getwX()+m->getwY()+zerosX+zerosY);
		REPORT(DETAILED,"dspout");
	}


}

