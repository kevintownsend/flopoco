/*
Floating-point pipeline generator for FloPoCo

Author : Florent de Dinechin

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
CeCILL license, 2008-2010.

All rights reserved
*/

#ifdef HAVE_SOLLYA

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <cstdlib>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "FPPipeline.hpp"

using namespace std;

namespace flopoco{
	
	extern vector<Operator*> oplist;
	
	FPPipeline::FPPipeline(Target* target, string func, int wE_, int wF_): 
	Operator(target), wE(wE_), wF(wF_) {
		// Name HAS to be unique!
		// will cause weird bugs otherwise
		ostringstream complete_name;
		complete_name << "Pipeline" << getNewUId(); 
		setName(complete_name.str());
		// r = x^2 + y^2 + z^2 example
		
		vector<FPNode*> fpNodeList;
		
		vector<string> xSqrInNames, ySqrInNames, zSqrInNames;
		xSqrInNames.push_back("x");
		ySqrInNames.push_back("y");
		zSqrInNames.push_back("z");
		
		FPNode *xsqr = new FPNode(FPNode::sqr, xSqrInNames);
		// 		fpNodeList.push_back(xsqr);
		// 		cleanupNodeList(fpNodeList,xsqr);
		FPNode *ysqr = new FPNode(FPNode::sqr, ySqrInNames);
		// 		fpNodeList.push_back(ysqr);
		FPNode *zsqr = new FPNode(FPNode::sqr, zSqrInNames);
		// 		fpNodeList.push_back(zsqr);
		
		vector<FPNode*> xsqrysqr;
		xsqrysqr.push_back(xsqr);
		xsqrysqr.push_back(ysqr);
		
		FPNode *addxsqrysqr = new FPNode(FPNode::adder, xsqrysqr);
		// 		fpNodeList.push_back(addxsqrysqr);
		
		vector<FPNode*> xyz;
		xyz.push_back(addxsqrysqr);
		xyz.push_back(zsqr);
		
		FPNode *addaddxsqrysqrzsqr = new FPNode(FPNode::adder, xyz, "r", true);
		fpNodeList.push_back(addaddxsqrysqrzsqr);
		
		for (unsigned i=0; i<fpNodeList.size(); i++){
			cout << "------------------" << endl;
			printTree( fpNodeList[i]);
		}
		
		for (unsigned i=0; i<fpNodeList.size(); i++){
			cout << "------------------" << endl;
			generateVHDL( fpNodeList[i], true);
		}
		
		
		
	}
	
	FPPipeline::~FPPipeline() {
	}
	
	
	
}
#endif //HAVE_SOLLYA


