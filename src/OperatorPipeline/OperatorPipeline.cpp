

#include "OperatorPipeline.hpp"

//*
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
//*/
//#ifdef HAVE_SOLLYA

//*
//#include "ExpressionParser.h"//this was not found idk why
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

program* p;

///TODO I used an implementation of this here so it can compile. Not sure if it realy have to be done here
int FlopocoExpressionlex(void){

}

int FlopocoExpressionparse(void){

}

int FlopocoExpressionerror(char *){

}

/** This functions instanciates and returns a pointer to an the proper
	operator for a given node.
  */
//Operator* instanciateOp(node n){
//	Operator* op;
//	if ( !strcmp(n->type, "") )
//		op=new Operator();
//	else if ( !strcmp(n->type, "") )
//		op=new Operator();
//	else if ( !strcmp(n->type, "") )
//		op=new Operator();
//	else
//		return null;
//
//	return op;
//}

//*/


//#define sumeofsquaresNaive
//#define sumeofsquares
//#define polynomial
//#define sqrtx2y2
//#define sqrtx2y2z2

using namespace std;

namespace flopoco{

OperatorPipeline::OperatorPipeline(Target* target, string filename, int wE_, int wF_):
    Operator(target), wE(wE_), wF(wF_) {
    // Name HAS to be unique!
    // will cause weird bugs otherwise
	
    ostringstream complete_name;
    complete_name << "OperatorPipeline" << getNewUId();
    setName(complete_name.str());
    // r = x^2 + y^2 + z^2 example
    srcFileName = "OperatorPipeline";

    // redirect stdin to the file pointer
    int my_stdin = dup(0);
    close(0);
    int fp = open(filename.c_str(), O_RDONLY, "r");

    dup2(fp, 0);
    FlopocoExpressionparse();
    close(fp);
    dup2(0, my_stdin);

    REPORT(DEBUG, "-----------------------------------");
    nodeList* head = p->assignList;
    while (head!=NULL){
        printExpression(head->n);
        REPORT(DEBUG,endl);
        head = head->next;
    }
    REPORT(DEBUG, "-----------------------------------");
    varList* headv = p->outVariableList;
    while (headv != NULL){
        REPORT(DEBUG, "out: variable " << headv->name	<< ";");
        headv = headv->next;
    }
    REPORT(DEBUG, "-----------------------------------");
    head = p->assignList;
    /* creates the computational tree our of the assignment list, by linking
           all variables already declared to their use */
    makeComputationalTree(NULL, head, head);

    REPORT(DEBUG, "NEW NODES: ------------------------");
    head = p->assignList;
    while (head!=NULL){
        printExpression(head->n);
        REPORT(DEBUG,endl);
        head = head->next;
    }
    REPORT(DEBUG, "-----------------------------------");

    /* create a node output list, where each node is a computational datapath.
           if in the user-provided outputList, some of the variables are part of
           intermediary computations for one, say larger, node these will not be added
           to the new list */

    nodeList* outList = createOuputList(p->assignList, p->outVariableList);

    REPORT(DEBUG, "PROPER OUT LIST: ------------------");
    nodeList* outListHead = outList;
    while (outListHead != NULL){
        printExpression( outListHead->n);
        outListHead = outListHead->next;
    }
    REPORT(DEBUG,endl);

    nodeList* oh = outList;



    optimise_tree();


    while (oh!=NULL){
        generateVHDL_c( oh->n, true);
        oh = oh->next;
    }
}

OperatorPipeline::~OperatorPipeline() {
}

void OperatorPipeline::optimise_tree(){

}

void OperatorPipeline::generateVHDL_c(node* n, bool top){
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
        c[t.str().length()]=0;
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
    //here we add the transformations from our tree to vhdl
    //I left the fpipeline case 1 as an exemple

    switch (n->type)
    {
        case 0:{  //input
            break;
        }
        case 1:{ //adder
            /*
        REPORT(DETAILED, " instance adder");

        op1 = new FPAddSinglePath(target_, wE, wF, wE, wF, wE, wF);
        oplist.push_back(op1);

        inPortMap( op1, "X", n->nodeArray->n->name);
        inPortMap( op1, "Y", n->nodeArray->next->n->name);
        outPortMap( op1, "R", n->name);

        ostringstream tmp;
        tmp << "adder" << getNewUId();
        vhdl << instance(op1, tmp.str())<<endl;
        //*/
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


}
//#endif //HAVE_SOLLYA


