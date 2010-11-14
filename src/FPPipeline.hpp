#ifndef FPPipeline_HPP
#define FPPipeline_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>


#include "FPAdder.hpp"
#include "FPAdderSinglePath.hpp"

#include "FPMultiplier.hpp"
#include "FPSquarer.hpp"
#include "FPSqrt.hpp"
#include "FPExp.hpp"
#include "FPLog.hpp"

#include "Operator.hpp"
// #include "HOTBM/sollya.h"	// Do NOT use libsollya from user's environment
// #include "UtilSollya.hh"

namespace flopoco{
	
	extern vector<Operator*> oplist;
	
	/** The FPPipeline class.  */
	class FPPipeline : public Operator {
		
		
		
		class FPNode {
			public: 
				typedef enum {
					adder, 
					subtracter,
					multiplier,
					sqr,
					sqrt
				} FPOperators;
				
				FPNode(FPOperators fop, vector<FPNode*> argInNodes, string argOName = "", bool outNode = false){
					outputNode = outNode;
					fun = fop;
					for (unsigned i=0; i< argInNodes.size(); i++){
						inNodes.push_back(argInNodes[i]);
						inNodesNames.push_back(argInNodes[i]->oName);
					}	
					
					if (argOName=="")
						oName = setNodeName();
					else
						oName = argOName;

					if (outNode){
						ostringstream u;
						u << Operator::getNewUId();
						oNameReal = oName;
						oName = "tmp_"+u.str()+oName;
					}
					
				};
				
				FPNode(FPOperators fop, vector<string> argInNodesNames, string argOName = "", bool outNode = false){
					outputNode = outNode;
					fun = fop;
					for (unsigned i=0; i< argInNodesNames.size(); i++)
						inNodesNames.push_back(argInNodesNames[i]);
					// 					inNodesNames(argInNodesNames);
					if (argOName=="")
						oName = setNodeName();
					else
						oName = argOName;
					if (outNode){
						ostringstream u;
						u << Operator::getNewUId();
						oNameReal = oName;
						oName = "tmp_"+u.str()+oName;
					}
				};
				
				string setNodeName(){
					ostringstream varName;
					varName << "var";
					if (inputNode()){
						for (unsigned i=0; i< inNodesNames.size(); i++)
							if (inNodesNames[i] != ""){
								varName << "_"<<inNodesNames[i];
							}else{}
					}else{
						for (unsigned i=0; i< inNodes.size(); i++)
							varName << "_"<<inNodes[i]->oName;
					}
					return varName.str();
				};
				
				bool inputNode(){
					return ( inNodes.size()==0?true:false);
				}
				
				
				FPOperators fun;
				string inName[2];
				string oName;
				string oNameReal;
				int inputs;
				vector<FPNode*> inNodes;
				vector<string> inNodesNames;
				bool outputNode;
		};
		
		public:
			
			/* TODO: Doxygen parameters*/ 
			FPPipeline(Target* target, string func, int wE, int wF);
			
			/**
			* FPPipeline destructor
			*/
			~FPPipeline();
			
			void printTree(FPNode *node){
				if (node->inputNode()){
					cout << node->oName << "=operation"<<node->fun<<"(";
					for (unsigned i=0; i<node->inNodesNames.size(); i++){
						cout << node->inNodesNames[i];
						if (i<node->inNodesNames.size()-1)
							cout<<",";
					}
					cout << ");"<<endl;
				}else{
					for (unsigned i=0; i< node->inNodes.size(); i++){
						printTree(node->inNodes[i]);
					}
					
					cout << node->oName << "=operation"<< node->fun<<"(";
					for (unsigned i=0; i< node->inNodes.size(); i++){
						cout << node->inNodes[i]->oName;
						if (i< node->inNodes.size()-1)
							cout << ",";
					}
					cout << ");"<<endl;
				}
			};
			
			
			void generateVHDL(FPNode *node, bool top){
				cout << "Generating VHDL ... " << endl;
				
				if (node->inputNode()){
					cout << "Input node detected " << endl;
					//we start at cycle 0, for now
					setCycle(0);
					//check if inputs are already declared. if not declare the inputs
					for (unsigned i=0; i< node->inNodesNames.size();i++){
						if (!isSignalDeclared(node->inNodesNames[i])){
							cout << "signal " << node->inNodesNames[i] << "   declared" << endl;
							addFPInput(node->inNodesNames[i], wE, wF);
						}
					}
					
				}else{
					for (unsigned i=0; i< node->inNodes.size(); i++){
						generateVHDL(node->inNodes[i], false);
					}
					//sync with all inputs. 
					for (unsigned i=0; i< node->inNodes.size(); i++)
						syncCycleFromSignal(node->inNodes[i]->oName);
				}
				
				if (node->outputNode)
					addFPOutput(node->oNameReal, wE, wF);
				
// 				cout << "    code " << endl;
				Operator* op1;
				//let's instantiate the proper operator

				switch (node->fun) {
					
					case FPNode::adder:{ 
						cout << " instance adder. Oplistsize =" <<oplist.size() << endl;
						
						op1 = new FPAdderSinglePath(target_, wE, wF, wE, wF, wE, wF);
						oplist.push_back(op1);

						inPortMap( op1, "X", node->inNodesNames[0]);
						inPortMap( op1, "Y", node->inNodesNames[1]);
						outPortMap( op1, "R", node->oName);
						
						ostringstream tmp;
						tmp << "adder" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}
					case FPNode::sqr:{
						cout << " instance squarer Oplistsize =" <<oplist.size()<<endl;
						op1 = new FPSquarer(target_, wE, wF, wF);
						oplist.push_back(op1);
						
						inPortMap( op1, "X", node->inNodesNames[0]);
						outPortMap( op1, "R", node->oName);
						
						ostringstream tmp;
						tmp << "squarer" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						cout << "    generated squarer instance" << endl;
						break;
					}
					default:{
						cout << "nothing else implemented yet" << endl;
						exit(-1);
					}
				}
				
				if (node->outputNode){
					syncCycleFromSignal(node->oName);
 					vhdl << tab << node->oNameReal << " <= " << node->oName << ";" << endl;
				}
				
			};
			
		protected:
			int wE;   /**< Exponent size*/ 
			int wF;  /**< Significand fraction size */
	};
	
}
#endif
