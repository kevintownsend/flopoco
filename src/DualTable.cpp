/*
  A generic class for DualTables of values

  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author :Radu Tudoran

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

*/


/* 
 The class contains a primitive Dual Port Memory block, that is nested inside the DualTable class. This is done in order to prevent the miss used of the inner object.
 The Memories are splitied in fundamental units(blocks) that can be recognized by the simulator as dual port memories. The limit of a memory block until it is recognized as
 a dual port one is 2^11*9 bits. However this value could be moved in the target in the eventuality that some other targets would have others such limits. However, on all models
 of Virtex 4 on which the design was tested this value holds. This class is an abstract one, and in order to be used it must be extended by a concret class which will implement the
 "function" method. It is very similar to the Table class.
*/

#include <iostream>
#include "utils.hpp"
#include "DualTable.hpp"


using namespace std;



namespace flopoco{

	extern vector<Operator*> oplist;



	int DualTable::double2input(double x){
		cerr << "Error, double2input is being used and has not been overriden";
		return 1;
	}

	double DualTable::input2double(int x) {
		cerr << "Error, input2double is being used and has not been overriden";
		return 1;
	}

	mpz_class DualTable::double2output(double x){
		cerr << "Error, double2output is being used and has not been overriden";
		return 0;
	}

	double DualTable::output2double(mpz_class x) {
		cerr << "Error, output2double is being used and has not been overriden";
		return 1;
	}



	DualTable::DualTable(Target* target, int _wIn, int _wOut, int _minIn, int _maxIn) : 
		Operator(target),
		wIn(_wIn), wOut(_wOut), minIn(_minIn), maxIn(_maxIn),target(target)
	{
	
		setCopyrightString("Radu Tudoran (2009)");

		//limitSingleMemory = intpow2(11)*9;
		limitSingleMemory = target->sizeOfMemoryBlock();
		
		// Set up the IO signals
	
		addInput ("X1"  , wIn);
		addOutput ("Y1"  , wOut);
		addInput ("X2"  , wIn);
		addOutput ("Y2"  , wOut);			
		
		if(intpow2(wIn)*wOut<=limitSingleMemory)
			{
		
				if(maxIn==-1) maxIn=(1<<wIn)-1;
				if(minIn<0) {
					cerr<<"ERROR in DualTable::DualTable, minIn<0\n";
					exit(EXIT_FAILURE);
				}
				if(maxIn>=(1<<wIn)) {
					cerr<<"ERROR in DualTable::DualTable, maxIn too large\n";
					exit(EXIT_FAILURE);
				}
				if((minIn==0) && (maxIn==(1<<wIn)-1)) 
					full=true;
				else
					full=false;
	
				nrOfMemBlocks=1;
	
	
	
			}
		else
			{
		
				maxIn= intlog2(limitSingleMemory / wOut);
				if(intpow2(maxIn)*wOut>limitSingleMemory)
					maxIn--;
	
	
		
				nrOfMemBlocks = intpow2(wIn-maxIn);
	
				if(nrOfMemBlocks * maxIn<wIn)
					nrOfMemBlocks++;
		
			
				if(minIn<0) {
					cerr<<"ERROR in DualTable::DualTable, minIn<0\n";
					exit(EXIT_FAILURE);
				}
	
				full=true;
	
				cerr << "WARNING : FloPoCo is building a DualTable with " << wIn <<" X "<< wOut<< " , it will be large." << endl;
		
			}
	
	
	}
	
	DualTable::DualTable(Target* target) : 
		Operator(target)
	{
		setCopyrightString("Radu Tudoran (2010)");
	}


	// We have to define this method because the constructor of DualTable cannot use the (pure virtual) function()
	void DualTable::outputVHDL(std::ostream& o, std::string name) {
		int i;
		mpz_class y;
		ostringstream data;

		if( intpow2(wIn) * wOut<=limitSingleMemory)
			{
	
		
				for(int v=minIn;v<=maxIn;v++)
					{data<<function(v)<<" ";
					}
				
		
				primitiveBlock  = new primitiveDualMemory(target,wIn,wOut,minIn,maxIn);
				primitiveBlock->setInputData(data);
				primitiveBlock  ->changeName(getName()+"primitiveBlock");	
				oplist.push_back(primitiveBlock  );
				inPortMapCst  (primitiveBlock  , "X1","X1");
				inPortMapCst  (primitiveBlock  , "X2","X2");
				outPortMap (primitiveBlock  , "Y1","data1");
				outPortMap (primitiveBlock  , "Y2","data2");
				vhdl << instance(primitiveBlock  ,"primitiveBlock");
		
				vhdl<<tab<<"Y1 <=data1;"<<endl;
				vhdl<<tab<<"Y2 <=data2;"<<endl;
		
			}
		else
			{
	
				ostringstream name1,name2,name3;
	
				int rangeHighAddress = wIn - maxIn;
		
				vhdl<<tab<<declare("address1L",maxIn)<<" <= X1"<<range(maxIn-1,0)<<";"<<endl;
				vhdl<<tab<<declare("address1H",rangeHighAddress,true)<<" <= X1"<<range(wIn-1,maxIn)<<";"<<endl;
				vhdl<<tab<<declare("address2L",maxIn)<<" <= X2"<<range(maxIn-1,0)<<";"<<endl;
				vhdl<<tab<<declare("address2H",rangeHighAddress,true)<<" <= X2"<<range(wIn-1,maxIn)<<";"<<endl;
	
	
				primitiveBlocks = (primitiveDualMemory**) calloc (nrOfMemBlocks,sizeof(primitiveDualMemory*));
		
	
				for(int c1=0;c1<nrOfMemBlocks;c1++)
					{
						data.str("");
						name1.str("");
						name1<<"primitiveBlock_"<<c1;
						name2.str("");
						name2<<"data1FromM_"<<c1;
						name3.str("");
						name3<<"data2FromM_"<<c1;
		
						int maxInInner= intpow2(maxIn);
		
						for(int v=minIn;v<=maxInInner-1-minIn;v++)
							{
								data<<function(v +maxInInner*c1)<<" ";
							}
		
						primitiveBlocks[c1]  = new primitiveDualMemory(target,maxIn,wOut,minIn,maxInInner-1-minIn);
						primitiveBlocks[c1]->setInputData(data);
						primitiveBlocks[c1]  ->changeName(getName()+name1.str());	
						oplist.push_back(primitiveBlocks[c1] );
						inPortMapCst  (primitiveBlocks[c1]  , "X1","address1L");
						inPortMapCst  (primitiveBlocks[c1]  , "X2","address2L");
						outPortMap (primitiveBlocks[c1]  , "Y1",name2.str());
						outPortMap (primitiveBlocks[c1]  , "Y2",name3.str());
						vhdl << instance(primitiveBlocks[c1]  ,name1.str());
		
			
					}
	
				vhdl<<endl;
				vhdl <<tab<< "  with address1H select  Y1 <= " << endl;
				for (int c1 = 0; c1 < nrOfMemBlocks; c1++) {
					name2.str("");
					name2<<"data1FromM_"<<c1;
					vhdl  << tab << name2.str()<<" when \"" << unsignedBinary(c1,rangeHighAddress) << "\"," << endl;
				}
				vhdl << tab << "\"";
				for (i = 0; i < wOut; i++) 
					vhdl << "-";
				vhdl <<  "\" when others;" << endl;
	
				vhdl<<endl;
				vhdl <<tab<< "  with address2H select  Y2 <= " << endl;
				for (int c1 = 0; c1 < nrOfMemBlocks; c1++) {
					name3.str("");
					name3<<"data2FromM_"<<c1;
					vhdl 	<<tab <<tab << name3.str()<<" when \"" << unsignedBinary(c1,rangeHighAddress) << "\"," << endl;
				}
				vhdl << tab<<tab << "\"";
				for (i = 0; i < wOut; i++) 
					vhdl << "-";
				vhdl <<  "\" when others;" << endl;
	
		
			}
	
		Operator::outputVHDL(o,  name);
	}






	int DualTable::size_in_LUTs() {
		return wOut*int(intpow2(wIn-target_->lutInputs()));
	}


}

	
	
