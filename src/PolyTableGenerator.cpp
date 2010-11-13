/*
 * Table Generator unit for FloPoCo
 *
 * Author : Mioara Joldes
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon

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
#include "PolyTableGenerator.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	PolyTableGenerator::PolyTableGenerator(Target* target, PiecewiseFunction* pf, int wInX, int wOutX, int n): 
		Table(target),	wInX_(wInX), wOutX_(wOutX), pwf(pf){
		
		setCopyrightString("Mioara Joldes (2010)");
		srcFileName = "PolyTableGenerator";
		/* Start initialization */
		setToolPrecision(165);
		int nrMaxIntervals=1024*1024;

		int	nrFunctions=(pwf->getPiecewiseFunctionArray()).size();
		int	iter;
		int guardBits =1;
		
		/*this is the size of the multipliers to be prefered*/
		/*in the current version this is not used yet */
		int msize =1;
		
		mpfr_t a;
		mpfr_t b;
		mpfr_t eps;
		mpfr_init2(a,getToolPrecision());
		mpfr_init2(b,getToolPrecision());
		mpfr_init2(eps,getToolPrecision());
	
		/* Convert parameters to their required type */
		mpfr_set_d(a,0.,GMP_RNDN);
		mpfr_set_d(b,1.,GMP_RNDN);
		mpfr_set_ui(eps, 1, GMP_RNDN);
		mpfr_mul_2si(eps, eps, -wOutX-1-guardBits, GMP_RNDN); /* eps< 2^{-woutX-1} */
	
	
		vector<sollya_node_t> polys;
		vector<mpfr_t*> errPolys;
	
		sollya_node_t tempNode,tempNode2, tempNode3, nDiff, sX,sY,aiNode;
		sollya_chain_t tempChain, tempChain2;
	
		mpfr_t ai;
		mpfr_t bi;
		mpfr_t zero;
		mpfr_t* mpErr;
		int nrIntCompleted=-1;
	
		mpfr_init2(ai,getToolPrecision());
		mpfr_init2(bi,getToolPrecision());
		mpfr_init2(zero,getToolPrecision());
 
		int k, errBoundBool;
		tempNode2=parseString("0"); /*ct part*/
		tempChain = makeIntPtrChainFromTo(0,n); /*monomials*/
 
		int sizeList[n+1];
		for (k=0;k<=n;k++)	
			sizeList[k]=0;

		for (iter=0; iter<nrFunctions;iter++){
			Function *fi;
			fi=(pwf->getPiecewiseFunctionArray(iter));
			
			REPORT(INFO, "Now approximating "<<fi->getName());
			/*start with one interval; subdivide until the error is satisfied for all functions involved*/
			int nrIntervals = 1;
			int precShift=0;
			errBoundBool =0;
		
			/*Convert the input string into a sollya evaluation tree*/ 
			tempNode = copyTree(fi->getSollyaNode());
		
			while((errBoundBool==0)&& (nrIntervals <=nrMaxIntervals)){
				errBoundBool=1; 
				REPORT(DETAILED, "Now trying with "<< nrIntervals <<" intervals");
				for (k=0; k<nrIntervals; k++){
					mpfr_set_ui(ai,k,GMP_RNDN);
					mpfr_set_ui(bi,1,GMP_RNDN);
					mpfr_div_ui(ai, ai, nrIntervals, GMP_RNDN);
					mpfr_div_ui(bi, bi, nrIntervals, GMP_RNDN);		
					mpfr_set_ui(zero,0,GMP_RNDN);
			
					aiNode = makeConstant(ai);
					sX = makeAdd(makeVariable(),aiNode);
					sY = substitute(tempNode, sX);
					if (sY == 0)
						REPORT(DEBUG, "Sollya error when performing range mapping.");
					if(verbose>=DEBUG){
						cout<<"\n-------------"<<endl;
						printTree(sY);
						cout<<"\nover: "<<sPrintBinary(zero)<<" "<< sPrintBinary(bi)<<"withprecshift:"<<precShift<<endl;
					}
					tempChain2 = makeIntPtrChainCustomized(wOutX+1+guardBits,n+1, precShift,msize ); //precision
				
					tempNode3 = FPminimax(sY, tempChain ,tempChain2, NULL, zero, bi, FIXED, ABSOLUTESYM, tempNode2,NULL);
					polys.push_back(tempNode3);
				
					if (verbose>=DEBUG){
						printTree(tempNode3);
						printf("\n");
					}
					/*Compute the error*/
					nDiff = makeSub(sY, tempNode3);
					mpErr= (mpfr_t *) safeMalloc(sizeof(mpfr_t));
					mpfr_init2(*mpErr,getToolPrecision());	
					uncertifiedInfnorm(*mpErr, nDiff, zero, bi, 501, getToolPrecision()); 
					if (verbose>=DEBUG){
						cout<< "infinite norm:"<<sPrintBinary(*mpErr)<<endl;
						cout<< "eps:"<<sPrintBinary(eps)<<endl;
					}
					errPolys.push_back(mpErr);
					if (mpfr_cmp(*mpErr, eps)>0) {
						errBoundBool=0; 
						REPORT(DETAILED, tab << "we have found an interval where the error is not small enough");
						REPORT(DETAILED, tab << "failed at interval "<< k+1 << "/" << nrIntervals);
						REPORT(DETAILED, tab << "proceed to splitting"); 
						/*erase the polys and the errors put so far for this function; keep the previous good ones intact*/
						polys.resize(nrIntCompleted+1);
						errPolys.resize(nrIntCompleted+1);
						nrIntervals=2 * nrIntervals;
						precShift=precShift+1;
						break;
					}
				}
			}
			if (errBoundBool==1){ 
				/*we have the good polynomials for one function*/
					REPORT(DETAILED, "the number of intervals is:"<< nrIntervals); 
					REPORT(DETAILED, "We proceed to the next function");
				/*we have to update the total size of the coefficients*/
				int ii=0;
				while (tempChain2!=NULL){
					int sizeNew=*((int *)first(tempChain2));
					tempChain2=tail(tempChain2);
					if (sizeNew>sizeList[ii] )sizeList[ii]= sizeNew;
					ii++;
				}
			nrIntCompleted=nrIntCompleted+nrIntervals;	
			nrIntArray.push_back(intlog2(nrIntervals-1));
			} 
		}
		/*here we have the good polynomials for all the piecewise functions*/
		/*clear some vars used*/

		mpfr_clear(ai);
		mpfr_clear(bi);
		mpfr_clear(zero);
	
		if (errBoundBool==1){
			if(verbose>=DEBUG){
				cout<< "the total number of intervals is:"<< nrIntCompleted<<endl; 
				cout<< "We proceed to the extraction of the coefficients:"<<endl; 
			}

			/*Get the maximum error*/
			mpfr_t *mpErrMax;
			mpErrMax=(mpfr_t*) safeMalloc(sizeof(mpfr_t));
			mpfr_init2(*mpErrMax, getToolPrecision());
			mpfr_set(*mpErrMax,*errPolys[0], GMP_RNDN);
		
			for (k=1;(unsigned)k<errPolys.size();k++){
				if (mpfr_cmp(*mpErrMax, *(errPolys[k]))<0)
					mpfr_set(*mpErrMax,*(errPolys[k]), GMP_RNDN);
			}
		
			maxError=(mpfr_t*) safeMalloc(sizeof(mpfr_t));
			mpfr_init2(*maxError,getToolPrecision());
			mpfr_set(*maxError,*mpErrMax,GMP_RNDN);
		
			mpfr_clear(*mpErrMax);
			free(mpErrMax);
			if (verbose>=DEBUG){
				cout<< "maximum error="<<sPrintBinary(*maxError)<<endl;
			}
			/*Extract coefficients*/
			vector<FixedPointCoefficient*> fpCoeffVector;
		
			k=0;
			for (k=0;(unsigned)k<polys.size();k++){
				if (verbose>=DEBUG){
					cout<<"\n----"<< k<<"th polynomial:----"<<endl;
					printTree(polys[k]);
				}
			
				fpCoeffVector = getPolynomialCoefficients(polys[k], sizeList);
				polyCoeffVector.push_back(fpCoeffVector);
			}
			/*setting of Table parameters*/
			wIn=intlog2(polys.size()-1);
			minIn=0;
			maxIn=polys.size()-1;
			wOut=0;
			for(k=0; (unsigned)k<coeffParamVector.size();k++){
				wOut=wOut+(*coeffParamVector[k]).getSize()+(*coeffParamVector[k]).getWeight()+1; /*a +1 is necessary for the sign*/
			}
		
			ostringstream name;
			/*Set up the name of the entity */
			name <<"PolyTableGenerator_"<<wIn<<"_"<<wOut;
			setName(name.str());
		
			/*Set up the IO signals*/
			addInput ("X"	, wIn, true);
			addOutput ("Y"	, wOut, true,1);
				
			/*This operator is combinatorial (in fact is just a ROM.*/
			//setCombinatorial();

			//generateDebugPwf();
			if (verbose>=INFO){	
				printPolynomialCoefficientsVector();
				REPORT(DETAILED, "Parameters for polynomial evaluator:");
				printCoeffParamVector();
			}
		}else{
			/*if we didn't manage to have the good polynomials for all the piecewise functions*/
			cout<< "PolyTableGenerator error: Could not approximate the given function(s)"<<endl; 
		}
	
		mpfr_clear(a);
		mpfr_clear(b);
		mpfr_clear(eps);
		free_memory(tempNode);
		free_memory(tempNode2);
		free_memory(nDiff);
		free_memory(sX);
		freeChain(tempChain,freeIntPtr);
		freeChain(tempChain2,freeIntPtr);
	}


	/*This constructor receives the function to be approximated as a string
	Look above to find one that receives the function as a Piecewise function already parsed*/
	PolyTableGenerator::PolyTableGenerator(Target* target, string func, int wInX, int wOutX, int n):
	Table(target) {
		
		setCopyrightString("Mioara Joldes (2010)");		

		/*parse the string, create a list of functions, create an array of f's, compute an approximation on each interval*/
		PiecewiseFunction *pf=new PiecewiseFunction(func);
		PolyTableGenerator(target, pf, wInX, wOutX, n);
	}


	PolyTableGenerator::PolyTableGenerator(Target* target, string func, int wInX, int wOutX, int n, double xmin, double xmax, double scale ): 
	Table(target),	wInX_(wInX), wOutX_(wOutX), f(new Function(func, xmin, xmax, scale)){
	
		/* Start initialization */
		setToolPrecision(165);
	
		/* End of initialization */
		//int verbose=1;
		int nrMaxIntervals=1024*1024;	
		/* Convert the input string into a sollya evaluation tree */
		sollya_node_t tempNode = f->getSollyaNode(); //function
		int guardBits =1;
		mpfr_t a;
		mpfr_t b;
		mpfr_t eps;
	
		mpfr_init2(a,getToolPrecision());
		mpfr_init2(b,getToolPrecision());
		mpfr_init2(eps,getToolPrecision());
	
		/* Convert parameters to their required type */
		mpfr_set_d(a,0.,GMP_RNDN);
		mpfr_set_d(b,1.,GMP_RNDN);
		mpfr_set_ui(eps, 1, GMP_RNDN);
		mpfr_mul_2si(eps, eps, -wOutX-1-guardBits, GMP_RNDN); // eps< 2^{-woutX-1}
	
		//sollya_node_t w=parseString("1"); //weight
		// sollya_range_t r;
		//r=guessDegree(tempNode, w, a, b, eps);
	
		sollya_node_t tempNode2=parseString("0"); //ct part
		sollya_chain_t tempChain = makeIntPtrChainFromTo(0,n); //monomials

		//start with one interval; subdivide until the error is satisfied
	
		int nrIntervals = 1;
		int precShift=0;
		vector<sollya_node_t> polys;
	
		vector<mpfr_t*> errPolys;
	
		sollya_node_t tempNode3, nDiff;
		sollya_chain_t tempChain2;
	
		mpfr_t ai;
		mpfr_t bi;
		mpfr_t zero;
		mpfr_t* mpErr;
	

		mpfr_init2(ai,getToolPrecision());
		mpfr_init2(bi,getToolPrecision());
	
		mpfr_init2(zero,getToolPrecision());
	
		int k;
		int errBoundBool =0;
		sollya_node_t sX,sY,aiNode;
	
		while((errBoundBool==0)&& (nrIntervals <=nrMaxIntervals)){
			errBoundBool=1; //suppose the nr of intervals is good
			//polys.reserve(nrIntervals);
			//errPolys.reserve(nrIntervals);
			if (verbose>=DEBUG){
			cout<<"trying with "<< nrIntervals <<" intervals"<<endl;
			}
			for (k=0; k<nrIntervals; k++){
				mpfr_set_ui(ai,k,GMP_RNDN);
				mpfr_set_ui(bi,1,GMP_RNDN);
				mpfr_div_ui(ai, ai, nrIntervals, GMP_RNDN);
				mpfr_div_ui(bi, bi, nrIntervals, GMP_RNDN);		
		
				mpfr_set_ui(zero,0,GMP_RNDN);
			
				aiNode = makeConstant(ai);
				sX = makeAdd(makeVariable(),aiNode);
				//sY = simplifyTreeErrorfree(substitute(tempNode, sX));
				sY = substitute(tempNode, sX);
				if (sY == 0)
				cout<<"Sollya error when performing range mapping."<<endl;
			
				if (verbose>=DEBUG){
				cout<<"\n-------------"<<endl;	
				printTree(sY);
				cout<<"\nover: "<<sPrintBinary(zero)<<" "<< sPrintBinary(bi)<<"withprecshift:"<<precShift<<endl;	
				}
			
				tempChain2 = makeIntPtrChainToFromBy(wOutX+1+guardBits,n+1, precShift); //precision
			
				//tempNode3 = FPminimax(firstArg, tempChain, tempChain2, tempChain3, a, b, resB, resC, tempNode, tempNode2);
				tempNode3 = FPminimax(sY, tempChain ,tempChain2, NULL,			zero, bi, FIXED, ABSOLUTESYM, tempNode2,NULL);
			
				polys.push_back(tempNode3);
				if (verbose>=DEBUG){
					printTree(tempNode3);
					printf("\n");
				}
			
				//Compute the error 
				nDiff = makeSub(sY, tempNode3);
				mpErr= (mpfr_t *) safeMalloc(sizeof(mpfr_t));
				mpfr_init2(*mpErr,getToolPrecision());	
				uncertifiedInfnorm(*mpErr, nDiff, zero, bi, 501/*default in sollya*/, getToolPrecision()); 
				if (verbose>=DEBUG){
				cout<< "infinite norm:"<<sPrintBinary(*mpErr)<<endl;
				cout<< "eps:"<<sPrintBinary(eps)<<endl;
				}
				errPolys.push_back(mpErr);
				if (mpfr_cmp(*mpErr, eps)>0) {
					errBoundBool=0; //we have found an interval where the error is not good
					if(verbose>=DEBUG){
						cout<< "we have found an interval where the error is not good, proceed to splitting"<<endl;
					}
					//k=nrIntervals;
					polys.clear();
					errPolys.clear();
					nrIntervals=2 * nrIntervals;
					precShift=precShift+1;
					break;
				}
			}
		} 
		if (errBoundBool==1){
			if(verbose==1){
				cout<< "the number of intervals is:"<< nrIntervals<<endl; 
				cout<< "We proceed to the extraction of the coefficients:"<<endl; 
			}
			//Get the maximum error
	
			mpfr_t *mpErrMax;
			mpErrMax=(mpfr_t*) safeMalloc(sizeof(mpfr_t));
			mpfr_init2(*mpErrMax, getToolPrecision());
			mpfr_set(*mpErrMax,*errPolys[0], GMP_RNDN);
		
			for (k=1;(unsigned)k<errPolys.size();k++){
				if (mpfr_cmp(*mpErrMax, *(errPolys[k]))<0)
					mpfr_set(*mpErrMax,*(errPolys[k]), GMP_RNDN);
			}
		
			maxError=(mpfr_t*) safeMalloc(sizeof(mpfr_t));
			mpfr_init2(*maxError,getToolPrecision());
			mpfr_set(*maxError,*mpErrMax,GMP_RNDN);
		
			mpfr_clear(*mpErrMax);
			free(mpErrMax);
			if (verbose>=DEBUG){
				cout<< "maximum error="<<sPrintBinary(*maxError)<<endl;
			}
			//Extract coefficients
				vector<FixedPointCoefficient*> fpCoeffVector;
		
				k=0;
			for (k=0;k<nrIntervals;k++){
					if (verbose>=DEBUG){
						cout<<"\n----"<< k<<"th polynomial:----"<<endl;
						printTree(polys[k]);
					}
			
					fpCoeffVector = getPolynomialCoefficients(polys[k], tempChain2);
					polyCoeffVector.push_back(fpCoeffVector);
			}

			/*****setting of Table parameters**/
			//int wInZ, minInZ, maxInZ, wOutZ;
			wIn=intlog2(nrIntervals-1);
			minIn=0;
			maxIn=nrIntervals-1;
			wOut=0;
			for(k=0; (unsigned)k<coeffParamVector.size();k++){
				wOut=wOut+(*coeffParamVector[k]).getSize()+(*coeffParamVector[k]).getWeight()+1; //a +1 is necessary for the sign
			}
	
			ostringstream name;
			/* Set up the name of the entity */
			name <<"PolyTableGenerator_"<<wIn<<"_"<<wOut;
			setName(name.str());
	
			// Set up the IO signals
			addInput ("X"	, wIn);
			addOutput ("Y"	, wOut);
			
			/* This operator is combinatorial (in fact is just a ROM.*/
	//		setCombinatorial();

			generateDebugPwf();
			/**********************************/
		
			if (verbose>=DEBUG){	
				printPolynomialCoefficientsVector();
				cout<<"Parameters for polynomial evaluator:"<<endl;
				printCoeffParamVector();
			}
		}
		else{
			cout<< "something went wrong"<<endl; 
		}
		
		mpfr_clear(a);
		mpfr_clear(b);

		//free_memory(tempNode);
		free_memory(tempNode2);
		//free_memory(tempNode3);

		freeChain(tempChain,freeIntPtr);
		freeChain(tempChain2,freeIntPtr);
		//finishTool();
	}

	PolyTableGenerator::~PolyTableGenerator() {
	}


	MPPolynomial* PolyTableGenerator::getMPPolynomial(sollya_node_t t){
		int degree=1,i;
		sollya_node_t *nCoef;
		mpfr_t *coef;
		
		//printTree(t);
		getCoefficients(&degree, &nCoef, t);
		//cout<<degree<<endl;
		coef = (mpfr_t *) safeCalloc(degree+1,sizeof(mpfr_t));
		
			
		for (i = 0; i <= degree; i++){
			mpfr_init2(coef[i], getToolPrecision());
			//cout<< i<<"th coeff:"<<endl;
			//printTree(getIthCoefficient(t, i));
			evaluateConstantExpression(coef[i], getIthCoefficient(t, i), getToolPrecision());
			if (verbose>=DEBUG){
				cout<< i<<"th coeff:"<<sPrintBinary(coef[i])<<endl;
			}
		}
			MPPolynomial* mpPx = new MPPolynomial(degree, coef);
			//Cleanup 
			for (i = 0; i <= degree; i++)
				mpfr_clear(coef[i]);
			free(coef);
		
		return mpPx;
	}

	vector<FixedPointCoefficient*> PolyTableGenerator::getPolynomialCoefficients(sollya_node_t t, sollya_chain_t c){
		int degree=1,i,size, weight;
		sollya_node_t *nCoef;
		mpfr_t *coef;
		sollya_chain_t cc;
		vector<FixedPointCoefficient*> coeffVector;
		FixedPointCoefficient* zz;
		
		//printTree(t);
		getCoefficients(&degree, &nCoef, t);
		//cout<<degree<<endl;
		coef = (mpfr_t *) safeCalloc(degree+1,sizeof(mpfr_t));
		cc=c;
			
		for (i = 0; i <= degree; i++)
		{
			mpfr_init2(coef[i], getToolPrecision());
			//cout<< i<<"th coeff:"<<endl;
			//printTree(getIthCoefficient(t, i));
			evaluateConstantExpression(coef[i], getIthCoefficient(t, i), getToolPrecision());
			if (verbose>=DEBUG){
				cout<< i<<"th coeff:"<<sPrintBinary(coef[i])<<endl;
			}
			size=*((int *)first(cc));
			cc=tail(cc);
			if (mpfr_sgn(coef[i])==0){
				weight=0;
				size=1;
			} 
			else 
				weight=mpfr_get_exp(coef[i]);
			
			zz= new FixedPointCoefficient(size, weight, coef[i]);
			coeffVector.push_back(zz);
			updateMinWeightParam(i,zz);
		}
			
			//Cleanup 
			for (i = 0; i <= degree; i++)
				mpfr_clear(coef[i]);
			free(coef);
			
		return coeffVector;
	}


	vector<FixedPointCoefficient*> PolyTableGenerator::getPolynomialCoefficients(sollya_node_t t, int* sizeList){
		int degree=1,i,size, weight;
		sollya_node_t *nCoef;
		mpfr_t *coef;
		
		vector<FixedPointCoefficient*> coeffVector;
		FixedPointCoefficient* zz;
		
		//printTree(t);
		getCoefficients(&degree, &nCoef, t);
		//cout<<degree<<endl;
		coef = (mpfr_t *) safeCalloc(degree+1,sizeof(mpfr_t));
		
		for (i = 0; i <= degree; i++)
		{
			mpfr_init2(coef[i], getToolPrecision());
			//cout<< i<<"th coeff:"<<endl;
			//printTree(getIthCoefficient(t, i));
			evaluateConstantExpression(coef[i], getIthCoefficient(t, i), getToolPrecision());
			if (verbose>=DEBUG){
				cout<< i<<"th coeff:"<<sPrintBinary(coef[i])<<endl;
			}
			size=sizeList[i];
			
			if (mpfr_sgn(coef[i])==0){
				weight=0;
				size=1;
			} 
			else 
			weight=mpfr_get_exp(coef[i]);
			
			zz= new FixedPointCoefficient(size, weight, coef[i]);
			coeffVector.push_back(zz);
			updateMinWeightParam(i,zz);
		}
			
		//Cleanup 
		for (i = 0; i <= degree; i++)
			mpfr_clear(coef[i]);
		free(coef);
		
		return coeffVector;
	}


	void PolyTableGenerator::updateMinWeightParam(int i, FixedPointCoefficient* zz)
	{
		if (coeffParamVector.size()<=(unsigned)i) {
			coeffParamVector.push_back(zz);
		}
		else if (mpfr_sgn((*(*coeffParamVector[i]).getValueMpfr()))==0) coeffParamVector[i]=zz;
		else if ( ((*coeffParamVector[i]).getWeight() <(*zz).getWeight()) && (mpfr_sgn(*((*zz).getValueMpfr()))!=0) )
		coeffParamVector[i]=zz;
	}


	vector<vector<FixedPointCoefficient*> > PolyTableGenerator::getPolynomialCoefficientsVector(){
		return polyCoeffVector;
	}


	void PolyTableGenerator::printPolynomialCoefficientsVector(){
		int i,j,nrIntervals, degree;
		vector<FixedPointCoefficient*> pcoeffs;
		nrIntervals=polyCoeffVector.size();

		for (i=0; i<nrIntervals; i++){	
			pcoeffs=polyCoeffVector[i];
			degree= pcoeffs.size();
			REPORT(DEBUG, "polynomial "<<i<<": ");
			for (j=0; j<degree; j++){
				REPORT(DEBUG, " "<<(*pcoeffs[j]).getSize()<< " "<<(*pcoeffs[j]).getWeight());
			}
		}
	}


	vector<FixedPointCoefficient*> PolyTableGenerator::getCoeffParamVector(){
		return coeffParamVector;
	}


	void PolyTableGenerator::printCoeffParamVector(){
		int j, degree;
		degree= coeffParamVector.size();
		for (j=0; j<degree; j++){		
			REPORT(DETAILED, " "<<(*coeffParamVector[j]).getSize()<< " "<<(*coeffParamVector[j]).getWeight()); 
		}
	}


	mpfr_t * PolyTableGenerator::getMaxApproxError(){
		return maxError;
	}


	void PolyTableGenerator::generateDebug(){
		int j;
		cout<<"f=";
		printTree(simplifyTreeErrorfree(f->getSollyaNode()));
		cout<<" wIn="<<wInX_<<" wOut="<<(-1)*wOutX_<<endl;
		cout<<"k="<<polyCoeffVector.size()<<" d="<<coeffParamVector.size()<<endl;
		cout<<"The size of the coefficients is:"<<endl;
		for (j=0; (unsigned)j<coeffParamVector.size(); j++){
			cout<<"c"<<j<<":"<<(*coeffParamVector[j]).getSize()+(*coeffParamVector[j]).getWeight()+1<<endl; 
		}
	}


	vector<int> PolyTableGenerator::getNrIntArray(){
		return nrIntArray;
	}


	void PolyTableGenerator::generateDebugPwf(){
		int j;
		cout<<"pwf=";
		cout<<pwf->getName()<<endl;
		cout<<" wIn="<<wInX_<<" wOut="<<(-1)*wOutX_<<endl;
		cout<<"k="<<polyCoeffVector.size()<<" d="<<coeffParamVector.size()<<endl;
		cout<<"The size of the branch is:"<<endl;
		for (j=0; (unsigned)j<getNrIntArray().size(); j++){
			cout<<j<<":"<<(getNrIntArray())[j]<<endl;
		}

		cout<<"The size of the coefficients is:"<<endl;
		for (j=0; (unsigned)j<coeffParamVector.size(); j++){
			cout<<"c"<<j<<":"<<(*coeffParamVector[j]).getSize()+(*coeffParamVector[j]).getWeight()+1<<endl; 
		}
	}


	sollya_chain_t PolyTableGenerator::makeIntPtrChainCustomized(int m, int n, int precshift, int msize) {
		int i,j, temp;
		int *elem;
		sollya_chain_t c;
		int tempTable[n+1];
		tempTable[0]=m;
		for (i=1; i<n; i++){
			temp=(tempTable[i-1]-precshift)/msize;
			if (temp!=0)
				tempTable[i]=temp*msize ;
			else	
				tempTable[i]=(tempTable[i-1]-precshift);
		}

//		tempTable[0]+=2;		
//		tempTable[1]+=2;
//		tempTable[2]+=2;
//		tempTable[4]+=2;
		tempTable[3]+=1;
		
	
		c = NULL;
	 
		for(j=n-1;j>=0;j--) {
			elem = (int *) malloc(sizeof(int));
			*elem = tempTable[j];
			c = addElement(c,elem);
		}
		return c;
	}


/****************************************************************************************/
/************Implementation of virtual methods of Class Table***************************/
	int PolyTableGenerator::double2input(double x){
		int result;
		cerr << "???	PolyTableGenerator::double2input not yet implemented ";
		exit(1);
		return result;
	}


	double	PolyTableGenerator::input2double(int x) {
		double y;
		cerr << "??? PolyTableGenerator::double2input not yet implemented ";
		exit(1);
		return(y);
	}


	mpz_class	PolyTableGenerator::double2output(double x){
		cerr << "???	PolyTableGenerator::double2input not yet implemented ";
		exit(1);
		return 0;
	}


	double	PolyTableGenerator::output2double(mpz_class x) {
		double y;
		cerr << "???	PolyTableGenerator::double2input not yet implemented ";
		exit(1);
	
		return(y);
	}


/***************************************************************************************/
/********************This is the implementation of the actual mapping*******************/
	mpz_class	PolyTableGenerator::function(int x)
	{
		mpz_class r=0; mpfr_t *cf; mpz_t c;char * z;
		int amount,j,nrIntervals, degree,trailingZeros,numberSize;
		vector<FixedPointCoefficient*> pcoeffs;
		nrIntervals=polyCoeffVector.size();
		if ((x<0) ||(x>=nrIntervals)) {}//x is not in the good range
		else{
			pcoeffs=polyCoeffVector[(unsigned)x];
			degree= pcoeffs.size();
			amount=0;
			//cout<<"polynomial "<<x<<": "<<endl;
			//r=mpz_class( 133955332 )+(mpz_class( 65664 )<< 27 )+(mpz_class( 64 )<< 44 )
			for (j=0; j<degree; j++){		
				//cout<<" "<<(*pcoeffs[j]).getSize()<< " "<<(*pcoeffs[j]).getWeight()<<endl; 
				r=r+(mpz_class(0)<<amount);
				
				
				cf=(*pcoeffs[j]).getValueMpfr();
				
				//cout<< j<<"th coeff:"<<sPrintBinary(*cf)<<endl;
				z=sPrintBinaryZ(*cf);
				//cout<< j<<"th coeff:"<<z<<" "<<strlen(z)<<endl;
				mpz_init(c);
				if (mpfr_sgn(*cf)!=0) {
					
				trailingZeros=(*coeffParamVector[j]).getSize()+(*pcoeffs[j]).getWeight()-strlen(z);
				//cout<<"Trailing zeros="<< trailingZeros<<endl;
				numberSize= (*coeffParamVector[j]).getSize()+(*coeffParamVector[j]).getWeight()+1 ;
					//mpz_set_str (mpz_t rop, char *str, int base) 
					mpz_set_str (c,z,2);
					if (mpfr_sgn(*cf)<0) {
						if (j==0)
								r=r+(((mpz_class(1)<<numberSize) -	((mpz_class(c)<<trailingZeros) - (mpz_class(1)<<2) ) )<<amount);
						else
							r=r+(((mpz_class(1)<<numberSize) -	(mpz_class(c)<<trailingZeros)	)<<amount);
					}
					else{
						if (j==0)
							r=r+(((mpz_class(c)<<trailingZeros) + (mpz_class(1)<<2) )<<amount);
						else
							r=r+((mpz_class(c)<<trailingZeros)<<amount);							
					}
					
				}
				else{
					if (j==0)
 						r=r+(((mpz_class(0)) + (mpz_class(1)<<2))<<amount);
 				else
 					r=r+(mpz_class(0)<<amount);
				} 
				
				amount=amount+(*coeffParamVector[j]).getSize()+(*coeffParamVector[j]).getWeight()+1;
			}
		}
		return r;
	}
/***************************************************************************************/

}
#endif //HAVE_SOLLYA


