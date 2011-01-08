#ifndef FPPipeline_HPP
#define FPPipeline_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "FPAdderSinglePath.hpp"

#include "FPMultiplier.hpp"
#include "FPSquarer.hpp"
#include "FPSqrt.hpp"
#include "FPExp.hpp"
#include "FPLog.hpp"
#include "FPSqrtPoly.hpp"
#include "FPSqrt.hpp"

#include "Operator.hpp"
#include "FPExpressions/ExpressionParserData.h"
// #include "HOTBM/sollya.h"	// Do NOT use libsollya from user's environment
// #include "UtilSollya.hh"

namespace flopoco{
	
	extern vector<Operator*> oplist;
	
	/** The FPPipeline class.  */
	class FPPipeline : public Operator {
		
		
		
		class FPNode {
			public: 
				typedef enum {
					input,
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
					
					if (fop==input){
						vinputNode = true;
						oName = argInNodes[0]->oName;
					}else{
						vinputNode = false;
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

					if (fop==input){
						vinputNode = true;
						oName = argInNodesNames[0];
					}else {
						vinputNode = false;
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
					return ( inNodes.size()==0?true:false) or (vinputNode);
				}
				
				
				FPOperators fun;
				bool vinputNode;
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
					cout << node->oName << "=InputOperation"<<node->fun<<"(";
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
// 					cout << "Input node detected " << endl;
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
					
					bool nodesAreInputs = true;
					for (unsigned i=0; i< node->inNodes.size(); i++)
						if (!node->inNodes[i]->inputNode()) 
							nodesAreInputs = false;

					if (!nodesAreInputs)
						nextCycle();
				}
				
				if (node->outputNode)
					addFPOutput(node->oNameReal, wE, wF);
				
// 				cout << "    code " << endl;
				Operator* op1;
				//let's instantiate the proper operator

				switch (node->fun) {

					case FPNode::input:{ 
						break;
					}
					
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

					case FPNode::multiplier:{ 
						cout << " instance adder. Oplistsize =" <<oplist.size() << endl;
						
						op1 = new FPMultiplier(target_, wE, wF, wE, wF, wE, wF);
						oplist.push_back(op1);

						inPortMap( op1, "X", node->inNodesNames[0]);
						inPortMap( op1, "Y", node->inNodesNames[1]);
						outPortMap( op1, "R", node->oName);
						
						ostringstream tmp;
						tmp << "multiplier" << getNewUId();
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
					case FPNode::sqrt:{
						cout << " instance sqrt Oplistsize =" <<oplist.size()<<endl;
#ifdef ha
						int degree = int ( floor ( double(wF) / 10.0) );
						op1 = new FPSqrtPoly(target_, wE, wF, 0, degree);
#else
						op1 = new FPSqrt(target_, wE, wF);//, 1, degree);
#endif
						oplist.push_back(op1);
						
						inPortMap( op1, "X", node->inNodesNames[0]);
						outPortMap( op1, "R", node->oName);
						
						ostringstream tmp;
						tmp << "squarer" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						cout << "    generated square root instance" << endl;
						break;
					}

					default:{
						cout << "nothing else implemented yet" << endl;
						exit(-1);
					}
				}
				
				if (node->outputNode){
					syncCycleFromSignal(node->oName);
					nextCycle();
 					vhdl << tab << node->oNameReal << " <= " << node->oName << ";" << endl;
				}
				
			};

			/* --------------------------------------------------------------- */
			void generateVHDL_c(node* n, bool top){
				cout << "Generating VHDL ... " << endl;
				
				if (n->type == 0){
					//we start at cycle 0, for now
					setCycle(0);
					//check if inputs are already declared. if not declare the inputs
					if (!isSignalDeclared(n->name)){
						cout << "signal " << n->name << "   declared" << endl;
						addFPInput(n->name, wE, wF);
					}
				}else{
					//iterate on all inputs
					nodeList* lh = n->nodeArray;
					while (lh!=NULL){
						generateVHDL_c(lh->n, false);
						lh=lh->next;
					}
					lh = n->nodeArray;
					while (lh!=NULL){
						syncCycleFromSignal(lh->n->name);
						lh=lh->next;
					}
					cout << "finished with node" << endl;
				}
				
				if (n->name==NULL){
					//assign a unique name;
					ostringstream t;
					t << "tmp_var_"<<getNewUId();
					cout << " new temporary variable created "<<t.str()<<" size="<<t.str().length() <<endl; 
					string w = t.str();
					char *c  = new char[t.str().length()+1];
//					c=(char*)malloc(t.str().length()*sizeof(char));
					c = strncpy(c, t.str().c_str(), t.str().length() );
					c[t.str().length()]=NULL;
					cout << " new temporary variable created "<< c <<" size="<<t.str().length() <<endl; 
					n->name = c;
				}

				ostringstream t;				
				if (n->isOutput){
					t << "out_" << n->name;
					addFPOutput(t.str(), wE, wF);
				}

				
				Operator* op1;
				//let's instantiate the proper operator

				switch (n->type) {

					case 0:{  //input
						break;
					}
					
					case 1:{ //adder 
						cout << " instance adder" << endl;
						
						op1 = new FPAdderSinglePath(target_, wE, wF, wE, wF, wE, wF);
						oplist.push_back(op1);

						inPortMap( op1, "X", n->nodeArray->n->name);
						inPortMap( op1, "Y", n->nodeArray->next->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "adder" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}
					case 2:{ //subtracter 
						cout << " instance subtracter"<< endl;
						
						op1 = new FPAdderSinglePath(target_, wE, wF, wE, wF, wE, wF);
						oplist.push_back(op1);

						ostringstream temp;
						temp << "minus_"<<n->nodeArray->next->n->name;
						vhdl<<tab<<declare(temp.str(),wE+wF+3) << " <= " << n->nodeArray->next->n->name <<range(wE+wF+2,wE+wF+1)
						                                       << " & not(" << n->nodeArray->next->n->name <<of(wE+wF)<<")"
						                                       << " & " << n->nodeArray->next->n->name << range(wE+wF-1,0)<<";"<<endl;  

						inPortMap( op1, "X", n->nodeArray->n->name);
						inPortMap( op1, "Y", temp.str());
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "adder" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}

					case 3:{ //multiplier
						cout << " instance adder"<< endl;
						
						op1 = new FPMultiplier(target_, wE, wF, wE, wF, wE, wF);
						oplist.push_back(op1);

						inPortMap( op1, "X", n->nodeArray->n->name);
						inPortMap( op1, "Y", n->nodeArray->next->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "multiplier" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}
					case 5:{ //squarer
						cout << " instance squarer" << endl;
						
						op1 = new FPSquarer(target_, wE, wF, wF);
						oplist.push_back(op1);

						inPortMap( op1, "X", n->nodeArray->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "squarer" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}
					case 6:{ //sqrt
						cout << " instance sqrt"<<endl;
#ifdef ha
						int degree = int ( floor ( double(wF) / 10.0) );
						op1 = new FPSqrtPoly(target_, wE, wF, 0, degree);
#else
						op1 = new FPSqrt(target_, wE, wF);//, 1, degree);
#endif
						oplist.push_back(op1);
						
						inPortMap( op1, "X", n->nodeArray->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "squarer" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						cout << "    generated square root instance" << endl;
						break;
					}
					case 7:{ //exponential
						cout << " instance exp" << endl;
						
						op1 = new FPExp(target_, wE, wF, 0, 0);
						oplist.push_back(op1);

						inPortMap( op1, "X", n->nodeArray->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "exponential" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}
					case 8:{ //logarithm
						cout << " instance log" << endl;
						
						op1 = new FPLog(target_, wE, wF, 9);
						oplist.push_back(op1);

						inPortMap( op1, "X", n->nodeArray->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "logarithm" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}


					default:{
						cout << "nothing else implemented yet" << endl;
						exit(-1);
					}
				}
				
				if (n->isOutput){
					syncCycleFromSignal(n->name);
					nextCycle();
 					vhdl << tab << "out_"<<n->name << " <= " << n->name << ";" << endl;
				}
				
			};


			
		protected:
			int wE;   /**< Exponent size*/ 
			int wF;  /**< Significand fraction size */
	};
	
}
#endif
