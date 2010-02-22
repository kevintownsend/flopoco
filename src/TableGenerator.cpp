/*
 * Table Generator unit for FloPoCo
 *
 * Author : Mioara Joldes
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
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
#include "TableGenerator.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	TableGenerator::TableGenerator(Target* target, string func, int wInX, int wOutX, int n, double xmin, double xmax, double scale ): // TODO extend list
	Table(target),	 wInX_(wInX), wOutX_(wOutX), f(*new Function(func, xmin, xmax, scale))   {
	
	/* Start initialization */
  
  
  setToolPrecision(165);
  
  /* End of initialization */
  //int verbose=0;
  int nrMaxIntervals=1024*1024;	
	/* Convert the input string into a sollya evaluation tree */
  sollya_node_t tempNode = f.getSollyaNode(); //function
  
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
  mpfr_mul_2si(eps, eps, -wOutX-2, GMP_RNDN); // eps< 2^{-woutX-1}
  
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
    if(verbose){
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
      
      if(verbose){
      cout<<"\n-------------"<<endl;	
      printTree(sY);
      cout<<"\nover: "<<sPrintBinary(zero)<<" "<< sPrintBinary(bi)<<"withprecshift:"<<precShift<<endl;	
      }
      
      tempChain2 = makeIntPtrChainToFromBy(wOutX+1,n+1, precShift); //precision
      
      //tempNode3 = FPminimax(firstArg, tempChain, tempChain2, tempChain3, a, b, resB, resC, tempNode, tempNode2);
      tempNode3 = FPminimax(sY, tempChain ,tempChain2, NULL,       zero, bi, FIXED, ABSOLUTESYM, tempNode2,NULL);
      
      polys.push_back(tempNode3);
      if (verbose){
      
       printTree(tempNode3);
       printf("\n");
      }
      
      //Compute the error 
		  nDiff = makeSub(sY, tempNode3);
		  mpErr= (mpfr_t *) safeMalloc(sizeof(mpfr_t));
		  mpfr_init2(*mpErr,getToolPrecision());  
			uncertifiedInfnorm(*mpErr, nDiff, zero, bi, 501/*default in sollya*/, getToolPrecision()); 
      if (verbose){
      cout<< "infinite norm:"<<sPrintBinary(*mpErr)<<endl;
      }
		  errPolys.push_back(mpErr);
		  if (mpfr_cmp(*mpErr, eps)>0) {
		    errBoundBool=0; //we have found an interval where the error is not good
		    if(verbose){
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
  
  if(verbose){
    cout<< "the number of intervals is:"<< nrIntervals<<endl; 
    cout<< "We proceed to the extraction of the coefficients:"<<endl; 
  }
  //Get the maximum error
  
   mpfr_t *mpErrMax;
   mpErrMax=(mpfr_t*) safeMalloc(sizeof(mpfr_t));
   mpfr_init2(*mpErrMax, getToolPrecision());
   mpfr_set(*mpErrMax,*errPolys[0], GMP_RNDN);
   
   for (k=1;k<errPolys.size();k++){
    if (mpfr_cmp(*mpErrMax, *(errPolys[k]))<0)
      mpfr_set(*mpErrMax,*(errPolys[k]), GMP_RNDN);
   }
   
   maxError=(mpfr_t*) safeMalloc(sizeof(mpfr_t));
   mpfr_init2(*maxError,getToolPrecision());
   mpfr_set(*maxError,*mpErrMax,GMP_RNDN);
   
   mpfr_clear(*mpErrMax);
   free(mpErrMax);
   if (verbose){
    cout<< "maximum error="<<sPrintBinary(*maxError)<<endl;
   }
  //Extract coefficients
		vector<FixedPointCoefficient*> fpCoeffVector;
		
    k=0;
   for (k=0;k<nrIntervals;k++){
      if (verbose){
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
	  for(k=0; k<coeffParamVector.size();k++){
	    wOut=wOut+(*coeffParamVector[k]).getSize()+(*coeffParamVector[k]).getWeight()+1; //a +1 is necessary for the sign
	  }
	  
		ostringstream name;
		/* Set up the name of the entity */
		name <<"TableGenerator_"<<wIn<<"_"<<wOut;
		setName(name.str());
		
		// Set up the IO signals
		addInput ("X"  , wIn);
		addOutput ("Y"  , wOut);
				
		/* This operator is combinatorial (in fact is just a ROM.*/
		setCombinatorial();
	
	
	  /**********************************/
	  
	 if (verbose){  
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

	TableGenerator::~TableGenerator() {
	}


MPPolynomial* TableGenerator::getMPPolynomial(sollya_node_t t){
	  int degree=1,i;
		sollya_node_t *nCoef;
		mpfr_t *coef;
		
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
				if (verbose){
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

vector<FixedPointCoefficient*> TableGenerator::getPolynomialCoefficients(sollya_node_t t, sollya_chain_t c){
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
				if (verbose){
		      cout<< i<<"th coeff:"<<sPrintBinary(coef[i])<<endl;
		    }
		    size=*((int *)first(cc));
		    cc=tail(cc);
		    //if (mpfr_sgn(coef[i])==0) weight=0;
		    //else 
		    
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

void TableGenerator::updateMinWeightParam(int i, FixedPointCoefficient* zz)
{
  if (coeffParamVector.size()<=(unsigned)i) {
    coeffParamVector.push_back(zz);
  }
  else if ((*coeffParamVector[i]).getWeight() <(*zz).getWeight()) 
  coeffParamVector[i]=zz;

}

vector<vector<FixedPointCoefficient*> > TableGenerator::getPolynomialCoefficientsVector(){
return polyCoeffVector;
}
void TableGenerator::printPolynomialCoefficientsVector(){
  int i,j,nrIntervals, degree;
  vector<FixedPointCoefficient*> pcoeffs;
  nrIntervals=polyCoeffVector.size();

  for (i=0; i<nrIntervals; i++){  
    pcoeffs=polyCoeffVector[i];
    degree= pcoeffs.size();
    cout<<"polynomial "<<i<<": "<<endl;
    for (j=0; j<degree; j++){     
      cout<<" "<<(*pcoeffs[j]).getSize()<< " "<<(*pcoeffs[j]).getWeight()<<endl; 
    }
  }
}

vector<FixedPointCoefficient*> TableGenerator::getCoeffParamVector(){
return coeffParamVector;
}
void TableGenerator::printCoeffParamVector(){
  int j, degree;
  
    degree= coeffParamVector.size();
    
    for (j=0; j<degree; j++){     
      cout<<" "<<(*coeffParamVector[j]).getSize()<< " "<<(*coeffParamVector[j]).getWeight()<<endl; 
    }
  
}
mpfr_t * TableGenerator::getMaxApproxError(){
return maxError;
}

void TableGenerator::generateDebug(){
 cout<<"f="<<endl;
 printTree(f.getSollyaNode());
 cout<<"wIn="<<wInX_<<"wOut="<<wOutX_<<endl;
 
}


/****************************************************************************************/
/************Implementation of virtual methods of Class Table***************************/

int TableGenerator::double2input(double x){
		int result;
		cerr << "???  TableGenerator::double2input not yet implemented ";
		exit(1);
		return result;
	}


	double  TableGenerator::input2double(int x) {
		double y;
		cerr << "??? TableGenerator::double2input not yet implemented ";
		exit(1);
		return(y);
	}

	mpz_class  TableGenerator::double2output(double x){
		cerr << "???  TableGenerator::double2input not yet implemented ";
		exit(1);
		return 0;
	}

	double  TableGenerator::output2double(mpz_class x) {
		double y;
		cerr << "???  TableGenerator::double2input not yet implemented ";
		exit(1);
  
		return(y);
	}

/***************************************************************************************/
/********************This is the implementation of the actual mapping*******************/
	mpz_class  TableGenerator::function(int x)
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
        //r=r+(mpz_class(1)<<amount);
        
        
        cf=(*pcoeffs[j]).getValueMpfr();
        
        //cout<< j<<"th coeff:"<<sPrintBinary(*cf)<<endl;
        z=sPrintBinaryZ(*cf);
        //cout<< j<<"th coeff:"<<z<<" "<<strlen(z)<<endl;
        mpz_init(c);
        if (mpfr_sgn(*cf)!=0) {
          
         trailingZeros=(*coeffParamVector[j]).getSize()+(*pcoeffs[j]).getWeight()-strlen(z);
         numberSize= (*coeffParamVector[j]).getSize()+(*coeffParamVector[j]).getWeight()+1 ;
          //mpz_set_str (mpz_t rop, char *str, int base) 
          mpz_set_str (c,z,2);
          if (mpfr_sgn(*cf)<0) {
            r=r+(((mpz_class(1)<<numberSize) -   (mpz_class(c)<<trailingZeros)   )<<amount);
          }
          else{
             r=r+((mpz_class(c)<<trailingZeros)<<amount);
          }
          
        }
        else{
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


