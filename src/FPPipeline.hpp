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
#include "ConstMult/FPRealKCM.hpp"

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
				REPORT(INFO, "Generating VHDL ... ");
				
				if (node->inputNode()){
					//we start at cycle 0, for now
					setCycle(0);
					//check if inputs are already declared. if not declare the inputs
					for (unsigned i=0; i< node->inNodesNames.size();i++){
						if (!isSignalDeclared(node->inNodesNames[i])){
							REPORT(DETAILED, "signal " << node->inNodesNames[i] << "   declared");
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
				
				Operator* op1;
				//let's instantiate the proper operator

				switch (node->fun) {

					case FPNode::input:{ 
						break;
					}
					
					case FPNode::adder:{ 
						REPORT(DETAILED, " instance adder. Oplistsize =" <<oplist.size());
						
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
						REPORT(DETAILED, " instance adder. Oplistsize =" <<oplist.size());
						
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
						REPORT(DETAILED, " instance squarer Oplistsize =" <<oplist.size());
						op1 = new FPSquarer(target_, wE, wF, wF);
						oplist.push_back(op1);
						
						inPortMap( op1, "X", node->inNodesNames[0]);
						outPortMap( op1, "R", node->oName);
						
						ostringstream tmp;
						tmp << "squarer" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						REPORT(INFO, "generated squarer instance");
						break;
					}
					case FPNode::sqrt:{
						REPORT(DETAILED, " instance sqrt Oplistsize =" <<oplist.size());
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
						REPORT(DETAILED, "    generated square root instance");
						break;
					}

					default:{
						cerr << "nothing else implemented yet" << endl;
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
				REPORT(DETAILED, "Generating VHDL ... ");
				
				if (n->type == 0){
					//we start at cycle 0, for now
					setCycle(0);
					//check if inputs are already declared. if not declare the inputs
					if (n->name!=NULL){
						if (!isSignalDeclared(n->name)){
							REPORT(DETAILED, "signal " << n->name << "   declared");
							addFPInput(n->name, wE, wF);
						}
					}else{
						//this is a constant, so it has no name, and is not declared
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
						if (lh->n->name!=NULL)
							syncCycleFromSignal(lh->n->name);
						lh=lh->next;
					}
					REPORT(DETAILED, "finished with node");
				}
				
				bool hadNoName = (n->name==NULL);
				
				if (n->name==NULL){
					//assign a unique name;
					ostringstream t;
					t << "tmp_var_"<<getNewUId();
					string w = t.str();
					char *c  = new char[t.str().length()+1];
					c = strncpy(c, t.str().c_str(), t.str().length() );
					c[t.str().length()]=NULL;
					REPORT(DETAILED, " new temporary variable created "<< c <<" size="<<t.str().length());
					n->name = c;
					REPORT(DETAILED, " the value was created for the constant " << n->value);
				}
				
				if ((hadNoName)&&(n->type == 0)){
					//it is a constant_expr
					mpfr_t mpx;
					mpfr_init2 (mpx, wF+1);
					mpfr_set_str (mpx, n->s_value, 10, GMP_RNDN);
					vhdl << tab << declare(n->name, wE+wF+3) << " <= \""<<fp2bin(mpx, wE, wF)<< "\";"<<endl;
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
						REPORT(DETAILED, " instance adder");
						
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
						REPORT(DETAILED, " instance subtracter");
						
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
						REPORT(DETAILED, " instance multiplier");
						if (((n->nodeArray->n->type==0)&&(n->nodeArray->n->s_value!=NULL)) || 
						    ((n->nodeArray->next->n->type==0)&&(n->nodeArray->next->n->s_value!=NULL))){
							REPORT(INFO, "constant node detected");
							ostringstream constant_expr, operand_name;
							if ((n->nodeArray->n->type==0)&&(n->nodeArray->n->s_value!=NULL)){
								//the first one is the constant
								constant_expr << n->nodeArray->n->s_value;
								operand_name << n->nodeArray->next->n->name;		
							}else{
								constant_expr << n->nodeArray->next->n->s_value;
								operand_name << n->nodeArray->n->name;		
							}
							
							REPORT(INFO, "Constant is "<< constant_expr.str());

							op1 = new FPRealKCM(target_,wE, wF, constant_expr.str());
							oplist.push_back(op1);
							
							inPortMap( op1, "X", operand_name.str());
							outPortMap( op1, "R", n->name);

							ostringstream tmp;
							tmp << "constant_multiplier" << getNewUId();
							vhdl << instance(op1, tmp.str())<<endl;
						}else{
							//we just plug-in a regular multiplier
							op1 = new FPMultiplier(target_, wE, wF, wE, wF, wE, wF);
							oplist.push_back(op1);

							inPortMap( op1, "X", n->nodeArray->n->name);
							inPortMap( op1, "Y", n->nodeArray->next->n->name);
							outPortMap( op1, "R", n->name);
						
							ostringstream tmp;
							tmp << "multiplier" << getNewUId();
							vhdl << instance(op1, tmp.str())<<endl;
						}
						break;
					}
					case 5:{ //squarer
						REPORT(DETAILED, " instance squarer");
						
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
						REPORT(DETAILED, " instance sqrt");
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
						break;
					}
					case 7:{ //exponential
						REPORT(DETAILED, " instance exp");
						
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
						REPORT(DETAILED, " instance log");
						
						op1 = new FPLog(target_, wE, wF, 9);
						oplist.push_back(op1);

						inPortMap( op1, "X", n->nodeArray->n->name);
						outPortMap( op1, "R", n->name);
						
						ostringstream tmp;
						tmp << "logarithm" << getNewUId();
						vhdl << instance(op1, tmp.str())<<endl;
						break;
					}
					case 17:{ //assignement
						vhdl << tab << declare( n->name, wE+wF+3) << "<= " << n->nodeArray->n->name <<";"<<endl;
						break;
					}

					default:{
						
						cerr << "nothing else implemented yet for operation code: "<<n->type << endl;
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
