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

//#define sumeofsquares
//#define polynomial
//#define sqrtx2y2
#define sqrtx2y2z2		

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
		
		ifstream f;
		f.open(func.c_str());//
		//print file contents
		
		string line;
		if (f.is_open())
		{
			while ( f.good() )
			{
				getline (f,line);
				cout << line << endl;
			}
			f.close();
		}

#ifdef sumeofsquares		
		vector<FPNode*> fpNodeList;
		
		vector<string> inx, iny, inz;
		inx.push_back("x");
		iny.push_back("y");
		inz.push_back("z");

		FPNode *nx = new FPNode(FPNode::input, inx);
		FPNode *ny = new FPNode(FPNode::input, iny);
		FPNode *nz = new FPNode(FPNode::input, inz);
		
		vector<FPNode*> vsqrx, vsqry,vsqrz;
		vsqrx.push_back(nx);
		vsqry.push_back(ny);
		vsqrz.push_back(nz);

		FPNode *nxsqr = new FPNode(FPNode::sqr, vsqrx);
		FPNode *nysqr = new FPNode(FPNode::sqr, vsqry);
		FPNode *nzsqr = new FPNode(FPNode::sqr, vsqrz);
		
		vector<FPNode*> vadd1;
		vadd1.push_back(nxsqr);
		vadd1.push_back(nysqr);
		
		FPNode *nadd1 = new FPNode(FPNode::adder, vadd1);
		
		vector<FPNode*> vadd2;
		vadd2.push_back(nadd1);
		vadd2.push_back(nzsqr);
		
		FPNode *nadd2 = new FPNode(FPNode::adder, vadd2, "r", true);
		fpNodeList.push_back(nadd2);
#endif

#ifdef polynomial		
		vector<FPNode*> fpNodeList;
		
		vector<string> inx;
		inx.push_back("x");
		FPNode *nx = new FPNode(FPNode::input, inx);

		vector<FPNode*> vsqrx;
		vsqrx.push_back(nx);
		FPNode *sqrx = new FPNode(FPNode::sqr, vsqrx);
		
		vector<string> ina2;
		ina2.push_back("a2");
		FPNode *na2 = new FPNode(FPNode::input, ina2);

		vector<FPNode*> vprod1;
		vprod1.push_back(sqrx);
		vprod1.push_back(na2);
		
		FPNode *nprod1 = new FPNode(FPNode::multiplier, vprod1);




		vector<string> ina1;
		ina1.push_back("a1");
		FPNode *na1 = new FPNode(FPNode::input, ina1);

		vector<FPNode*> vprod2;
		vprod2.push_back(nx);
		vprod2.push_back(na1);
		
		FPNode *nprod2 = new FPNode(FPNode::multiplier, vprod2);

		vector<string> ina0;
		ina0.push_back("a0");
		FPNode *na0 = new FPNode(FPNode::input, ina0);

		vector<FPNode*> vadd1;
		vadd1.push_back(na0);
		vadd1.push_back(nprod2);
		
		FPNode *nadd1 = new FPNode(FPNode::adder, vadd1);


		vector<FPNode*> vadd2;
		vadd2.push_back(nadd1);
		vadd2.push_back(nprod1);
		
		FPNode *nadd2 = new FPNode(FPNode::adder, vadd2, "r", true);
		fpNodeList.push_back(nadd2);

#endif


#ifdef sqrtx2y2		
		vector<FPNode*> fpNodeList;
		
		vector<string> inx;
		inx.push_back("x");
		FPNode *nx = new FPNode(FPNode::input, inx);

		vector<FPNode*> vsqrx;
		vsqrx.push_back(nx);
		FPNode *sqrx = new FPNode(FPNode::sqr, vsqrx);

		vector<string> iny;
		iny.push_back("y");
		FPNode *ny = new FPNode(FPNode::input, iny);

		vector<FPNode*> vsqry;
		vsqry.push_back(ny);
		FPNode *sqry = new FPNode(FPNode::sqr, vsqry);

		vector<FPNode*> vadd1;
		vadd1.push_back(sqrx);
		vadd1.push_back(sqry);
		
		FPNode *nadd1 = new FPNode(FPNode::adder, vadd1);

		vector<FPNode*> vsqrt1;
		vsqrt1.push_back(nadd1);

		FPNode *nsqrt1 = new FPNode(FPNode::sqrt, vsqrt1, "r", true);
		fpNodeList.push_back(nsqrt1);

#endif

#ifdef sqrtx2y2z2		
		vector<FPNode*> fpNodeList;
		
		vector<string> inx;
		inx.push_back("x");
		FPNode *nx = new FPNode(FPNode::input, inx);

		vector<FPNode*> vsqrx;
		vsqrx.push_back(nx);
		FPNode *nsqrx = new FPNode(FPNode::sqr, vsqrx);

		vector<string> iny;
		iny.push_back("y");
		FPNode *ny = new FPNode(FPNode::input, iny);

		vector<FPNode*> vsqry;
		vsqry.push_back(ny);
		FPNode *nsqry = new FPNode(FPNode::sqr, vsqry);

		vector<string> inz;
		inz.push_back("z");
		FPNode *nz = new FPNode(FPNode::input, inz);

		vector<FPNode*> vsqrz;
		vsqrz.push_back(nz);
		FPNode *nsqrz = new FPNode(FPNode::sqr, vsqrz);

		vector<FPNode*> vadd1;
		vadd1.push_back(nsqrx);
		vadd1.push_back(nsqry);
		FPNode *nadd1 = new FPNode(FPNode::adder, vadd1);

		vector<FPNode*> vadd2;
		vadd2.push_back(nadd1);
		vadd2.push_back(nsqrz);
		FPNode *nadd2 = new FPNode(FPNode::adder, vadd2);

		vector<FPNode*> vsqrt1;
		vsqrt1.push_back(nadd2);

		FPNode *nsqrt1 = new FPNode(FPNode::sqrt, vsqrt1, "r", true);
		fpNodeList.push_back(nsqrt1);

#endif


		
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


