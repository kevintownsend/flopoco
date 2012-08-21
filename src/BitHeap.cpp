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
#include "Plotter.hpp"
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

	void BitHeap::setPlotter(Plotter* plotter_)
	{
		plotter = plotter_;
	}

	Plotter* BitHeap::getPlotter()
	{
		return plotter;
	}




	int BitHeap::newUid(unsigned w){
		return uid[w]++;
	}	

	int BitHeap::getGUid(){
		return guid;
	}	



	WeightedBit* BitHeap::computeLatest(unsigned w, int c0, int c1)
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

			addBit(weight+i, s.str(),"",1);
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

			addBit(weight+i, s.str(),"", 1);
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
					{
						if((mulBlocks[j]->getNext()==NULL)&&(mulBlocks[i]->getPrevious()==NULL))
						{	REPORT(DETAILED,"mulblocks[j]= "<<mulBlocks[j]->getWeight()<<" mulBlock[i]= "<<mulBlocks[i]->getWeight());
							REPORT(DETAILED,"j<=i");
							mulBlocks[j]->setNext(mulBlocks[i]);
							mulBlocks[i]->setPrevious(mulBlocks[j]);
						}	
					}						
					else
					{
						if((mulBlocks[i]->getNext()==NULL)&&(mulBlocks[j]->getPrevious()==NULL))	
						{	REPORT(DETAILED,"mulblocks[j]= "<<mulBlocks[j]->getWeight()<<" mulBlock[i]= "<<mulBlocks[i]->getWeight());
							REPORT(DETAILED,"j>i");
							mulBlocks[i]->setNext(mulBlocks[j]);
							mulBlocks[j]->setPrevious(mulBlocks[i]);
						}	
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

							addBit(weight,s.str(),"",1);
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

						addBit(weight,s.str(),"",1);
					}
				}
			}
		}
	}

	int BitHeap::computeStage()
	{
		return (op->getCurrentCycle()*stagesPerCycle + op->getCriticalPath()/elementaryTime); 
	
	}


	void  BitHeap::addBit(unsigned w, string rhs, string comment, int type)
	{
			
		// ignore bits beyond the declared maxWeight
		if(w >= maxWeight)
			return; 

		WeightedBit* bit= new WeightedBit(getGUid(), newUid(w), w, type, op->getCurrentCycle(), op->getCriticalPath()) ;
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
		REPORT(DEBUG, "name " << bit->getName());
		
		printColumnInfo(w);	
	};



	void BitHeap::removeBit(unsigned weight, int dir)
	{
		list<WeightedBit*>& l=bits[weight];

		WeightedBit* bit;
    
		//if dir=0 the bit will be removed from the begining of the list, else from the end of the list of weighted bits
		if(dir==0)
			l.pop_front();
		else if(dir==1)
			l.pop_back();

		

		REPORT(DEBUG,"remove bit from column " << weight);
	}

  
	void BitHeap::elemReduce(unsigned i, BasicCompressor* bc, int type)
	{
		REPORT(DEBUG, "Entering elemReduce() ");
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
        
		string in_concat=join("in_concat_bh", getGUid(), "_");
		string out_concat=join("out_concat_bh", getGUid(), "_");
		string compressor=join("Compressor_bh", getGUid(), "_");


		for(unsigned j=0; j<bc->height.size(); j++)
			{
				if (bc->getColumnSize(j)==1)
					op->vhdl << tab << op->declare(join(in_concat, compressorIndex,"_", inConcatIndex), bc->getColumnSize(j)) 
					         << "(0) <= " << signal[j].str() << ";" << endl;
				else
					op->vhdl << tab << op->declare(join(in_concat, compressorIndex,"_", inConcatIndex), bc->getColumnSize(j)) 
					         << " <= " << signal[j].str() << ";" << endl;
				op->inPortMap(bc, join("X",j), join(in_concat, compressorIndex,"_", inConcatIndex));
				++inConcatIndex;
			}

		op->outPortMap(bc, "R", join(out_concat, compressorIndex,"_", outConcatIndex));

		op->vhdl << tab << op->instance(bc, join(compressor, compressorIndex));

		removeCompressedBits(i,bc->getColumnSize(0));
		if(bc->getColumnSize(1)!=0)
			removeCompressedBits(i+1,bc->getColumnSize(1));
		
		// add the bits, at the current (global) instant.
		addBit(i, join(out_concat, compressorIndex,"_", outConcatIndex, "(0)"),"",type);
		addBit(i+1,	join(out_concat, compressorIndex,"_", outConcatIndex, "(1)"),"",type);
		if(!((bc->getColumnSize(0)==3) && (bc->getColumnSize(1)==0)))
			addBit(i+2, join(out_concat, compressorIndex,"_", outConcatIndex, "(2)"),"",type);

		REPORT(DEBUG, "Exiting elemReduce() ");

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
				REPORT(DEBUG, "element "<<i<<" cycle = "<<(*it)->getCycle() 
					<< " and cp = "<<(*it)->getCriticalPath((*it)->getCycle()));
				i++; 
			}
	}



	void BitHeap::generateCompressorVHDL()
	{
		op->vhdl << tab << endl << "-- Beginning of code generated by BitHeap::generateCompressorVHDL" << endl;

		// add the constant bits to the actual bit heap 
		//
		WeightedBit* firstBit = getFirstSoonestBit();
        
		int minCycle = firstBit->getCycle();
		double minCP = firstBit->getCriticalPath(minCycle);
		op->setCycle(minCycle); 
		op->setCriticalPath(minCP);
		op->vhdl << endl << tab << "-- Adding the constant bits" << endl;
		for (unsigned w=0; w<maxWeight; w++)
			if (1 == ((constantBits>>w) & 1) )

				addBit(w, "'1'","",2);

		
		generateSupertileVHDL();
		
			

		
		//initializeDrawing();
		generatePossibleCompressors();
        


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
                
				REPORT(DEBUG, "didCompress " << didCompress);

				plotter->heapSnapshot(didCompress, stage);
				compress(stage);
				stage++;
			}

		plotter->heapSnapshot(true, stage);
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


				while(i<maxWeight)
				{
					//REPORT(INFO, "i= "<< i << " cnt= " << cnt[i]);
					// Now we are sure cnt[i] is 3
					if (i==maxWeight-1)
					{
						applyCompressor3_2(i);
					}
					else
					{
						if((cnt[i+1]==3) || (cnt[i+1]==1))
						{
							applyCompressor3_2(i);
							do
							{
								i++;
							}
							while(cnt[i]!=3);

						}
						else
						{
							if(cnt[i+1]==2)
							{
								int j=i;
								do
								{
									i++;
								}
								while(cnt[i]==2);
								REPORT(INFO, "j= "<< j << " i-1= " << i-1);
								latestBit = getLatestBit(j, i-1); 
								if(latestBit)
								{
									op->setCycle(  latestBit ->getCycle()  );
									op->setCriticalPath(   latestBit ->getCriticalPath(op->getCurrentCycle()));
									op->manageCriticalPath( op->getTarget()->localWireDelay() + op->getTarget()->adderDelay(i-j) );

									stage = computeStage();
								}
								applyAdder(j, i-1);
								
								while((i<=maxWeight) && (cnt[i]!=3))
								{
									i++;
								}
							}
						}
					}
				}



#if 0
				while(doLoop)
					{
						// Now we are sure cnt[i] is 3;
						//look for the first two
						REPORT(DEBUG, " looking for the first two: before, i=" << i << " cnt[i]= " << cnt[i]);
						while (i<maxWeight && cnt[i]==3) 
							i++;
						REPORT(DEBUG, "                             after, i=" << i << " cnt[i]= " << cnt[i]);
                        
						)
						if (i==maxWeight) {
							applyCompressor3_2(maxWeight-1);
							doLoop=false;
						}
						else {
							// Now we are sure there is at least one 2
							first2=i;
						
							REPORT(DEBUG, " first3=" << first3 <<" first2=" << first2);
							if(first3<first2-1) {
								for(int j=first3; j<first2-1; j++){
									applyCompressor3_2(j);
								}
							}
							REPORT(DEBUG, "look for the first three again... ");
							
							//look for the first three again
							while (i<maxWeight &&  cnt[i]<3) 
								i++;
							REPORT(DEBUG, " found i= " << i);
							
							if (i==maxWeight) 
							{	
								latestBit = getLatestBit(first2-1, maxWeight-1); 
								if(latestBit)
								{
									op->setCycle(  latestBit ->getCycle()  );
									op->setCriticalPath(   latestBit ->getCriticalPath(op->getCurrentCycle()));
									op->manageCriticalPath( op->getTarget()->localWireDelay() + op->getTarget()->adderDelay(maxWeight-first2+1) );

									stage = computeStage();
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

									stage = computeStage();
								}


								applyAdder(first2-1, first3-1);	                       
							}
						}
						

					}
#endif

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
	  
	  	stage = computeStage();

		plotter->heapSnapshot(true, stage);	
#if 0
		offsetY += 20 + getMaxHeight() * 10;
		drawConfiguration(offsetY);
		closeDrawing(offsetY);
#endif
		
		//final addition
		generateFinalAddVHDL();

		plotter->plotBitHeap();
		
	}


	//returns a pointer to the first bit which has the smallest cycle & CP combination
	WeightedBit* BitHeap::getFirstSoonestBit()
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
		REPORT(DEBUG, "Entering applyCompressor3_2(" << col << ") ");
		unsigned i;
		for(i=0; i<possibleCompressors.size(); i++)
			{
				if((possibleCompressors[i]->getColumnSize(0)==3) && (possibleCompressors[i]->getColumnSize(1)==0))
					break; // exit the loop
			}
		REPORT(DEBUG, "Using Compressor3_2, reduce column " << col);
		elemReduce(col, possibleCompressors[i], 4);
		usedCompressors[i]=true;
		REPORT(DEBUG, "Exiting applyCompressor3_2(" << col << ") ");
	}




	// addition stretching from weights col0(LSB) to msb included
	// assumes cnt has been set up
	void BitHeap::applyAdder(int lsb, int msb)
	{

		stringstream inAdder0, inAdder1, cin;

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
		string inAdder0Name = join("inAdder0_bh", getGUid(), "_");
		string inAdder1Name = join("inAdder1_bh", getGUid(), "_");
		string cinName = join("cin_bh", getGUid(), "_");
		string outAdder = join("outAdder_bh", getGUid(), "_");

		op->vhdl << tab << op->declare(join(inAdder0Name, adderIndex), msb-lsb+2) << " <= \'0\' & " << inAdder0.str() << endl;
		op->vhdl << tab << op->declare(join(inAdder1Name, adderIndex), msb-lsb+2) << " <= \'0\' & " << inAdder1.str() << endl;
		op->vhdl << tab << op->declare(join(cinName, adderIndex)) << " <= " << cin.str() << endl;

#if 0

		IntAdder* adder = new IntAdder(op->getTarget(), msb-lsb+2);
		op->getOpListR().push_back(adder);
		
		op->inPortMap(adder, "X", join(inAdder0Name,adderIndex));
		op->inPortMap(adder, "Y", join(inAdder1Name,adderIndex));
		op->inPortMap(adder, "Cin", join(cinName,adderIndex));
		
		op->outPortMap(adder, "R", join(outAdder,adderIndex));

		op->vhdl << tab << op->instance(adder, join("Adder_bh", getGUid(), "_", adderIndex));
				
#else
		op->vhdl << tab << op->declare(join(outAdder, adderIndex), msb-lsb+2) << " <= " 
                        <<  join(inAdder0Name, adderIndex) << " + " << join(inAdder1Name,adderIndex) << " + " << join(cinName, adderIndex) << ";" <<endl ;

#endif


		for(int i=lsb; i<msb + 2 ; i++)	
		{
			addBit(i, join(outAdder, adderIndex,"(",i-lsb,")"),"",3); //adder working as a compressor = type 3 for added bit
		}
 

		adderIndex++;

	}


	//returns a pointer to the bit which has the biggest criticalPath+cycle combination, considering the given column limits (both included)
	WeightedBit* BitHeap::getLatestBit(int lsbColumn, int msbColumn)
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

		string inAdder0Name = join("finaldderIn0_bh", getGUid());
		string inAdder1Name = join("finalAdderIn1_bh", getGUid());
		string outAdderName = join("finalAdderOut_bh", getGUid());

		op->vhdl << tab << op->declare(inAdder0Name, maxWeight-minWeight) << "<= " << inAdder0.str() << endl;
		op->vhdl << tab << op->declare(inAdder1Name, maxWeight-minWeight) << " <= " << inAdder1.str() << endl;
		op->vhdl << tab << op->declare(outAdderName, maxWeight-minWeight+1) 
		         << " <= ('0' & "<< inAdder0Name << ") + ('0' & " << inAdder1Name << ");" << endl;
		op->vhdl << tab << "-- concatenate all the compressed chunks" << endl;
		//result		
		op->vhdl << tab << op->declare(join("CompressionResult", guid), (maxWeight+1)) << " <= " << outAdderName;


		//adding the rightmost bits
		for(int i=chunkDoneIndex-1; i>=0; i--)
			op->vhdl <<  " & " << join("tempR_bh", guid, "_", i);

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
					op->vhdl << tab << op->declare(join("tempR_bh", guid, "_", chunkDoneIndex), w-minWeight, true) << " <= " ;
				else
					op->vhdl << tab << op->declare(join("tempR_bh", guid, "_", chunkDoneIndex), 1, false) << " <= " ;
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
		
		int addx=0;
		int addy=0;
		
		if(topX<0)
		addx=0-topX;
		
		
		if(topY<0)
		addy=0-topY;	


		//int zerosY=18-m->getwY();
		
		//PIPELINE!!!!
		op->setCycleFromSignal(m->getInputName1());
		op->syncCycleFromSignal(m->getInputName2());
		op->manageCriticalPath(  op->getTarget()->DSPMultiplierDelay() ) ;
	//	op->setCycle(uid);
		if(uid==0)	
			op->vhdl << tab << op->declare(join("DSPch",i,"_",uid), m->getwX()+m->getwY()+zerosX+zerosY) 
				<< " <= (" <<zg(zerosX)<<" & "<<zg(addx)<<" & " << input1<<range(botX,topX-addx)<<") * (" <<zg(zerosY) <<" & "<<zg(addy)<<" & "
			    << input2 <<range(botY,topY-addy)<<");"<<endl;
		else
			op->vhdl << tab << op->declare(join("DSP",i,"_",uid), m->getwX()+m->getwY()+zerosX+zerosY) 
					<< " <= (" <<zg(zerosX)<<" & "<<zg(addx)<<" & " << input1<<range(botX,topX-addx)<<") * (" <<zg(zerosY) <<" & "<<zg(addy)<<" & "
			    << input2 <<range(botY,topY-addy)<<");"<<endl;

	
		if(uid==0)
			s<<join("DSPch",i,"_",uid);
		else
			
			s<<join("DSP",i,"_",uid);

		m->setSignalName(s.str());
		m->setSignalLength(m->getwX()+m->getwY()+zerosX+zerosY);
		REPORT(DETAILED,"dspout");
	}


}

