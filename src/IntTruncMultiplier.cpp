/*
  Truncated Integer Multiplier for FloPoCo
 
  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
 * Authors : Sebastian Banescu, Bogdan Pasca, Radu Tudoran

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.
*/

#include <typeinfo>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <limits.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntMultiplier.hpp"
#include "IntTruncMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;


	IntTruncMultiplier::IntTruncMultiplier(Target* target, int wX, float ratio, int k,int uL):
		Operator(target), wX(wX), wY(wX),ratio(ratio),targetPrecision(k),useLimits(uL){
		isSquarer = true;
		ostringstream name;
		name <<"IntTruncSquarer_"<<wX;
		setName(name.str());
	
		setCopyrightString("Sebastian Banescu, Bogdan Pasca, Radu Tudoran 2010");
	
		addInput ("X", wX);
		addOutput("R", wX + wX- k); 
		
		wInX = wX;
		wInY = wX;
		vn=wInX + 2* getExtraWidth();
		vm=wInY + 2* getExtraHeight();
		vnme = vn-getExtraWidth();		
		vmme = vm - getExtraHeight() ;
		
		nrDSPs = estimateDSPs();
		nrSoftDSPs = 0;
		cout << "Nr of estimated DSP blocks " << nrDSPs << endl;
		truncationOffset = estimateNrOfDiscardedCols(k);
		cout << "Nr of discarded cols " << truncationOffset << endl;
		nrOfShifts4Virtex=4;
		float const scale=100.0;
		costDSP = ( (1.0+scale) - scale * ratio );
		//~ cout<<"Cost DSP is "<<costDSP<<endl;
		costLUT = ( (1.0+scale) - scale * (1-ratio) ) /  ((float) target_->getEquivalenceSliceDSP() );
	
		runAlgorithm();
		//test for insertSoftDSPs 		
			
			//~ int n=vn;
			//~ int m=vm;
			//~ int exw=getExtraWidth();
			//~ int exh=getExtraHeight();
	
			//~ mat = new int*[m];
			//~ for(int i=0;i<m;i++)
			//~ {
				//~ mat[i] = new int [n];
				//~ for(int j=0;j<n;j++)
					//~ mat[i][j]=0;
			//~ }
			
			
			//~ int w, h;
			//~ target->getDSPWidths(w, h);
			
			//~ cout<<w<<" "<<h<<endl;
			//~ int shift = 17;
			
			
			//~ cout<<"ExtraWidth="<<exw<<" ExtraHeight="<<exh<<endl;
			
			//~ nrDSPs=2;
			//~ globalConfig = new DSP*[nrDSPs];
			//~ globalConfig[0] = new DSP(shift, w, h*2);	
			//~ globalConfig[0]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[0]->setTopRightCorner(exw, exh);
			//~ globalConfig[0]->setBottomLeftCorner(exw+w-1, exh+h*2-1);
			
			/*
			globalConfig[1] = new DSP(shift, w, h*2);	
			globalConfig[1]->setNrOfPrimitiveDSPs(2);
			globalConfig[1]->rotate();
			globalConfig[1]->setTopRightCorner(exw , exh +2*h  );
			globalConfig[1]->setBottomLeftCorner(exw + h*2-1, exh +2*h+ w-1);
			*/
			//~ globalConfig[1] = new DSP(shift, w, h*2);	
			//~ globalConfig[1]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[1]->rotate();
			//~ globalConfig[1]->setTopRightCorner(exw +w, exh   );
			//~ globalConfig[1]->setBottomLeftCorner(exw +w+ h*2-1, exh + w-1);
			
			
			//~ target_->getDSPWidths(dspWidth,dspHeight);
						
			
			//~ int test;
			//~ insertSoftDSPswithLimits(globalConfig);
			//~ cout<<"Nr of part "<<test<<endl;
			
			
		}
		
	IntTruncMultiplier::IntTruncMultiplier(Target* target, int wX, int wY, float ratio, int k,int uL):
		Operator(target), wX(wX), wY(wY), ratio(ratio),targetPrecision(k),useLimits(uL){
		
		isSquarer = false;	
		ostringstream name;
		name <<"IntTruncMultiplier_"<<wX<<"_"<<wY;
		setName(name.str());
	
		setCopyrightString("Sebastian Banescu, Bogdan Pasca, Radu Tudoran 2010");
	
		addInput ("X", wX);
		addInput ("Y", wY);
		addOutput("R", wX + wY- k); 
	
		wInX = wX;
		wInY = wY;
		vn=wInX + 2* getExtraWidth();
		vm=wInY + 2* getExtraHeight();
		vnme = vn - getExtraWidth();		
		vmme = vm - getExtraHeight() ;
		nrDSPs = estimateDSPs();
		nrSoftDSPs = 0;
		cout << "Nr of estimated DSP blocks " << nrDSPs << endl;
		truncationOffset = estimateNrOfDiscardedCols(k);
		cout << "Nr of discarded cols " << truncationOffset << endl;
		nrOfShifts4Virtex=4;
		float const scale=100.0;
		costDSP = ( (1.0+scale) - scale * ratio );
		//~ cout<<"Cost DSP is "<<costDSP<<endl;
		costLUT = ( (1.0+scale) - scale * (1-ratio) ) /  ((float) target_->getEquivalenceSliceDSP() );
		
		runAlgorithm();
		
		/*
		DSP** configuration;
		//configuration = (DSP**) malloc(sizeof(DSP*));
		//configuration[0] = (DSP*)malloc(sizeof(DSP));
		configuration = new DSP*[2];
		
		configuration[0] = target->createDSP();
		configuration[0]->setTopRightCorner(9,6);
		configuration[0]->setBottomLeftCorner(32,22);
		
		configuration[1] = target->createDSP();
		configuration[1]->setTopRightCorner(9,23);
		configuration[1]->setBottomLeftCorner(32,39);
		configuration[1]->setShiftIn(configuration[0]);
		configuration[0]->setShiftOut(configuration[1]);
		
		nrDSPs = 2;
		multiplicationInDSPs(configuration);
		vector<SoftDSP*> softDSPs;
		SoftDSP* sd = new SoftDSP(39, 15, 41, 20);
		softDSPs.push_back(sd);
		multiplicationInSlices(softDSPs);
		displayAll(configuration, softDSPs);
		
		// FOLLOWING CODE VALIDATES THE FINAL TILING CONFIGURATION
		if (isTilingValid(configuration, softDSPs, k)){
			cerr << "This is a good tiling" << endl;
		}else{
			cerr << "This is NOT a good tiling. Redo tiling " << endl;
		}
		
		
		//from here it will start the algorithm
		nrOfShifts4Virtex=4;
		int n=vn;
		int m=vm;
		mat = new int*[m];
			for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
			
			
		 

		//SoftDSP** softDSPs = new SoftDSP*[100];
		//SoftDSP* d1 = new SoftDSP(9, 9, 11, 29);
		//softDSPs[0] = d1;
		//nrSoftDSPs++;
		
		
		//runAlgorithm();
		*/
		/* FOLLOWING CODE VALIDATES THE FINAL TILING CONFIGURATION
		if (isTilingValid(configuration, softDSPs, k)){
			cerr << "This is a good tiling" << endl;
		}else{
			cerr << "This is NOT a good tiling. Redo tiling " << endl;
		}
		*/
		/*
		//test radu
			
			wInX=wX;
			wInY=wY;
			vn=wX + 2* getExtraWidth();
			vm=wY + 2* getExtraHeight();
			vnme = vn-getExtraWidth();		
			vmme = vm - getExtraHeight();
			nrOfShifts4Virtex=4;
			
			int n=vn;
			int m=vm;
			int exw=getExtraWidth();
			int exh=getExtraHeight();
	
			mat = new int*[m];
			for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
			
			nrDSPs=4;
			globalConfig = new DSP*[nrDSPs];
			int w, h;
			target->getDSPWidths(w, h);
			
			cout<<w<<" "<<h<<endl;
			int shift = 17;
			*/
			
				//test for -- configurations
			//~ globalConfig[0] = new DSP(shift, w, h*2);	
			//~ globalConfig[0]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[0]->setTopRightCorner(exw, exh);
			//~ globalConfig[0]->setBottomLeftCorner(exw+w-1, exh+h*2-1);
			
			//~ globalConfig[1] = new DSP(shift, w, h*2);	
			//~ globalConfig[1]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[1]->setTopRightCorner(exw , exh +2*h );
			//~ globalConfig[1]->setBottomLeftCorner(exw +w-1, exh + 4*h-1);
			
			
			
			//~ globalConfig[2] = new DSP(shift, w, h*2);	
			//~ globalConfig[2]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[2]->rotate();
			//~ globalConfig[2]->setTopRightCorner(vnme - h*4, exh  );
			//~ globalConfig[2]->setBottomLeftCorner(vnme -2*h-1 , exh +w-1);
			
			
			//~ globalConfig[3] = new DSP(shift, w, h*2);	
			//~ globalConfig[3]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[3]->rotate();
			//~ globalConfig[3]->setTopRightCorner(vnme - h*2, exh );
			//~ globalConfig[3]->setBottomLeftCorner(vnme-1, exh +w-1);
			
			
					//test for Z configurations
			//~ globalConfig[0] = new DSP(shift, w, h*2);	
			//~ globalConfig[0]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[0]->setTopRightCorner(exw, exh);
			//~ globalConfig[0]->setBottomLeftCorner(exw+w-1, exh+h*2-1);
			
			//~ globalConfig[1] = new DSP(shift, w, h*2);	
			//~ globalConfig[1]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[1]->setTopRightCorner(exw +w, exh +h );
			//~ globalConfig[1]->setBottomLeftCorner(exw +w+ w-1, exh + h+2*h-1);
			
			
			
			//~ globalConfig[2] = new DSP(shift, w, h*2);	
			//~ globalConfig[2]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[2]->rotate();
			//~ globalConfig[2]->setTopRightCorner(vnme - h*3, exh  );
			//~ globalConfig[2]->setBottomLeftCorner(vnme -h-1 , exh +w-1);
			
			
			//~ globalConfig[3] = new DSP(shift, w, h*2);	
			//~ globalConfig[3]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[3]->rotate();
			//~ globalConfig[3]->setTopRightCorner(vnme - h*2, exh +w);
			//~ globalConfig[3]->setBottomLeftCorner(vnme-1, exh +w+w-1);
						
						
						//test for L configurations
			//~ globalConfig[0] = new DSP(shift, w, h*2);	
			//~ globalConfig[0]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[0]->setTopRightCorner(exw, exh);
			//~ globalConfig[0]->setBottomLeftCorner(exw+w-1, exh+h*2-1);
			
			//~ globalConfig[1] = new DSP(shift, w, h*2);	
			//~ globalConfig[1]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[1]->rotate();
			//~ globalConfig[1]->setTopRightCorner(exw, exh +h*2 );
			//~ globalConfig[1]->setBottomLeftCorner(exw + h*2-1, exh + h*2+w-1);
			
			
			
			//~ globalConfig[2] = new DSP(shift, w, h*2);	
			//~ globalConfig[2]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[2]->rotate();
			//~ globalConfig[2]->setTopRightCorner(vnme - h*2, exh  );
			//~ globalConfig[2]->setBottomLeftCorner(vnme -1 , exh +w-1);
			
			
			//~ globalConfig[3] = new DSP(shift, w, h*2);	
			//~ globalConfig[3]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[3]->setTopRightCorner(vnme - w, exh +w);
			//~ globalConfig[3]->setBottomLeftCorner(vnme-1, exh +w+h*2-1);
			
						
			//~ //display(globalConfig);
			
			//~ bindDSPs(globalConfig);
			
			
			//test for insertSoftDSPs 			
			
			//~ wInX=wX;
			//~ wInY=wY;
			//~ vn=wX + 2* getExtraWidth();
			//~ vm=wY + 2* getExtraHeight();
			//~ vnme = vn-getExtraWidth();		
			//~ vmme = vm - getExtraHeight();
			//~ nrOfShifts4Virtex=4;
			
			//~ int n=vn;
			//~ int m=vm;
			//~ int exw=getExtraWidth();
			//~ int exh=getExtraHeight();
	
			//~ mat = new int*[m];
			//~ for(int i=0;i<m;i++)
			//~ {
				//~ mat[i] = new int [n];
				//~ for(int j=0;j<n;j++)
					//~ mat[i][j]=0;
			//~ }
			
			
			//~ int w, h;
			//~ target->getDSPWidths(w, h);
			
			//~ cout<<w<<" "<<h<<endl;
			//~ int shift = 17;
			
			
			//~ cout<<"ExtraWidth="<<exw<<" ExtraHeight="<<exh<<endl;
			
			//~ nrDSPs=2;
			//~ globalConfig = new DSP*[nrDSPs];
			//~ globalConfig[0] = new DSP(shift, w, h*2);	
			//~ globalConfig[0]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[0]->setTopRightCorner(exw, exh);
			//~ globalConfig[0]->setBottomLeftCorner(exw+w-1, exh+h*2-1);
			
			/*
			globalConfig[1] = new DSP(shift, w, h*2);	
			globalConfig[1]->setNrOfPrimitiveDSPs(2);
			globalConfig[1]->rotate();
			globalConfig[1]->setTopRightCorner(exw , exh +2*h  );
			globalConfig[1]->setBottomLeftCorner(exw + h*2-1, exh +2*h+ w-1);
			*/
			//~ globalConfig[1] = new DSP(shift, w, h*2);	
			//~ globalConfig[1]->setNrOfPrimitiveDSPs(2);
			//~ globalConfig[1]->rotate();
			//~ globalConfig[1]->setTopRightCorner(exw +w, exh   );
			//~ globalConfig[1]->setBottomLeftCorner(exw +w+ h*2-1, exh + w-1);
			
			
			//~ target_->getDSPWidths(dspWidth,dspHeight);
						
			
			//~ int test;
			//~ insertSoftDSPswithLimits(globalConfig);
			//~ cout<<"Nr of part "<<test<<endl;
	}

	IntTruncMultiplier::~IntTruncMultiplier() {
	}
	
	bool IntTruncMultiplier::isTilingValid(DSP** configuration, vector<SoftDSP*> softDSPs, int k){
		//mpfr_t targetError;
		mpfr_init2(targetError,1000);
		mpfr_set_ui(targetError, 2, GMP_RNDN);
		mpfr_pow_si(targetError, targetError, k, GMP_RNDN);

		mpfr_t roundingError;
		mpfr_init2(roundingError,1000);
		mpfr_set_ui(roundingError, 2, GMP_RNDN);
		mpfr_pow_si(roundingError, roundingError, k-1, GMP_RNDN);
				
		
		//cout << "The TargetError is : " << mpfr_get_d(targetError, GMP_RNDN) << endl;


		mpfr_t *approxError;
		approxError = evalTruncTilingErrorInverted(configuration, softDSPs);
		double apxError = mpfr_get_d(*approxError, GMP_RNDN);
		//cout << "The trunc error is : " << apxError << endl;
		
		if (apxError < 0) // then next loop will never finish (REASON: some multipliers overlap)
		{
			cout << "Unexpected value for trunc error. Execution stopped." << endl;
			return false;
		}
		/* error after compensation */
		mpfr_t t;
		mpfr_init2(t, 1000);
		mpfr_t t2;
		mpfr_init2(t2, 1000);
		
		bool found = false;
		
		mpfr_set_ui( t, 0, GMP_RNDN);
		if (mpfr_cmp(t, *approxError) == 0)
		{
			mpfr_set(*approxError, t, GMP_RNDN);
			found = true;
		}
		
		int i=0;
		do{
			mpfr_set_ui( t, 2, GMP_RNDN);
			mpfr_pow_si(  t, t, i, GMP_RNDN); 
			mpfr_mul_ui( t2, t, 2, GMP_RNDN);
		
			if  ((mpfr_cmp(t, *approxError) <= 0) && (mpfr_cmp(*approxError, t2) <= 0)){
				found = true;
				mpfr_set(*approxError, t, GMP_RNDN);		
			}else{
				i++;
			} 
		
		}while(! found);

		//cout << "The trunc after compensation : " << mpfr_get_d(*approxError, GMP_RNDN) << endl;
		//cout << "Rounding error is : " << mpfr_get_d(roundingError, GMP_RNDN) << endl;


		//mpfr_t errorSum;
		mpfr_init2(errorSum, 1000);
		mpfr_add( errorSum, roundingError, *approxError, GMP_RNDN);
		
		
		//cout << "The total error is : " << mpfr_get_d(errorSum, GMP_RNDN) << endl;
		
		if ( mpfr_cmp( errorSum, targetError) <= 0){
			return true;
		}else{
			return false;
		}
	}
	
	void IntTruncMultiplier::printConfiguration(DSP** configuration, vector<SoftDSP*> softDSPs){
		if (configuration!=NULL){
			int i=0;
			int xB,xT,yB,yT;
			for(i=0; i<nrDSPs; i++){
				configuration[i]->getTopRightCorner(xT,yT);
				configuration[i]->getBottomLeftCorner(xB,yB);
				cout << "HARD DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
			}
		}

		int xB,xT,yB,yT;
		for (unsigned k=0; k < softDSPs.size(); k++){
			softDSPs[k]->trim(vnme, vmme);
			softDSPs[k]->getTopRightCorner(xT,yT);
			softDSPs[k]->getBottomLeftCorner(xB,yB);
			cout << "SOFT DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
		}
	}
	
	mpfr_t *IntTruncMultiplier::evalTruncTilingError(DSP** configuration, SoftDSP** softDSPs){
		/* fist we get the maximal sum */
		mpfr_t *fullSum;
		cout << "wX wY are " << wX << "    " << wY << endl;
		fullSum = evalMaxValue(wX, wY);
		cout << "Full Sum is :" << mpfr_get_d(*fullSum, GMP_RNDN) << endl;
		int extW = getExtraWidth();
		int extH = getExtraHeight(); 
		
		if (configuration!=NULL){
			int i=0;
			int xB,xT,yB,yT;
			
			for (i=0; i<nrDSPs; i++){
				configuration[i]->getTopRightCorner(xT,yT);
				configuration[i]->getBottomLeftCorner(xB,yB);
				xT -= extW;
				xB -= extW;
				yT -= extH;
				yB -= extH;
				int power = xT + yT;
				mpfr_t* currentSum;
				currentSum = evalMaxValue(min(xB,wX) - xT + 1, min(yB,wY) - yT +1);
				mpfr_t s;
				mpfr_init2(s, 1000);
				mpfr_set_ui( s, 2, GMP_RNDN);
				mpfr_pow_si( s, s, power, GMP_RNDN);
				mpfr_mul( *currentSum, *currentSum, s, GMP_RNDN);
				mpfr_sub( *fullSum, *fullSum, *currentSum, GMP_RNDN);
				cout << "HARD DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
				cout << "This one is :" << mpfr_get_d(*currentSum, GMP_RNDN) << endl;
				i++;
			}
		}
		
		int xB,xT,yB,yT;
		for (int k=0; k < nrSoftDSPs; k++){
			softDSPs[k]->getTopRightCorner(xT,yT);
			softDSPs[k]->getBottomLeftCorner(xB,yB);
			xT -= extW;
			xB -= extW;
			yT -= extH;
			yB -= extH;
			softDSPs[k]->trim(vnme, vmme);
			int power = xT + yT;
			mpfr_t* currentSum;
			if ((xB>xT) && (yB > yT))
				currentSum = evalMaxValue(xB - xT + 1, yB - yT + 1);
			else{
				currentSum = (mpfr_t *) malloc( sizeof( mpfr_t));
				mpfr_init2( *currentSum, 1000);
				mpfr_set_si( *currentSum, 0, GMP_RNDN); 
			}
			
			
			mpfr_t s;
			mpfr_init2(s, 1000);
			mpfr_set_ui( s, 2, GMP_RNDN);
			mpfr_pow_si( s, s, power, GMP_RNDN);
			mpfr_mul( *currentSum, *currentSum, s, GMP_RNDN);
			mpfr_sub( *fullSum, *fullSum, *currentSum, GMP_RNDN);
			cout << "SOFT DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
			cout << "This one is :" << mpfr_get_d(*currentSum, GMP_RNDN) << endl;
		}
		return fullSum;
	}
	
	mpfr_t *IntTruncMultiplier::evalTruncTilingErrorInverted(DSP** configuration, vector<SoftDSP*> softDSPs){
		/* fist we get the maximal sum */
		mpfr_t *fullSum;
		cout << endl << endl << "wX wY are " << wX << "    " << wY << endl;
		fullSum = evalMaxValue(wX, wY);
		if (isSquarer) // only half of tiling grid
			mpfr_div_si(*fullSum, *fullSum, 2, GMP_RNDU);
		int extW = getExtraWidth();
		int extH = getExtraHeight(); 
		
		if (configuration!=NULL){
			int i=0;
			int xB,xT,yB,yT;
			for(i=0; i<nrDSPs; i++){
				configuration[i]->getTopRightCorner(xT,yT);
				configuration[i]->getBottomLeftCorner(xB,yB);
				xT -= extW;
				xB -= extW;
				yT -= extH;
				yB -= extH;
				xB = min(xB,wX-1); 
				xT = max(xT,0);
				yB = min(yB,wY-1);
				yT = max(yT,0);
				//int power = xT + yT;
				int power = (wX - xB-1) + (wY - yB-1);
				cout << "Power " << power << endl; 
				mpfr_t* currentSum;
				currentSum = evalMaxValue(min(xB,wX-1) - max(xT,0) + 1, min(yB,wY-1) - max(yT,0) +1);
				double area=0.0;
				
				if (isSquarer) // subtract parts that are out of triangle
				{
					int z, q;
					for (z=xT; z<=xB; z++)
						for (q=yT; q<=yB; q++)
							if ((wX-z) > q)
								area += pow(2.0, (double)((wX - z-1) + (wY - q-1)));
				}
				
				mpfr_t s;
				mpfr_init2(s, 1000);
				mpfr_set_ui( s, 2, GMP_RNDN);
				mpfr_pow_si( s, s, power, GMP_RNDN);
				mpfr_mul( *currentSum, *currentSum, s, GMP_RNDN);
				mpfr_sub_d( *currentSum, *currentSum, area, GMP_RNDN);
				mpfr_sub( *fullSum, *fullSum, *currentSum, GMP_RNDN);
				cout << "HARD DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
				cout << "This one is :" << mpfr_get_d(*currentSum, GMP_RNDN) << endl;
			}
		}
		
		int xB,xT,yB,yT;
		for (unsigned k=0; k < softDSPs.size(); k++){
			softDSPs[k]->trim(vnme, vmme);
			softDSPs[k]->getTopRightCorner(xT,yT);
			softDSPs[k]->getBottomLeftCorner(xB,yB);
			xT -= extW;
			xB -= extW;
			yT -= extH;
			yB -= extH;
			xB = min(xB,wX-1); 
			xT = max(xT,0);
			yB = min(yB,wY-1);
			yT = max(yT,0);
			//int power = xT + yT;
			int power = (wX - xB-1) + (wY - yB-1);
			//cout << "Power " << power << endl; 
			mpfr_t* currentSum;
			if ((xB>=xT) && (yB >=yT))
				currentSum = evalMaxValue(xB - xT + 1, yB - yT + 1);
			else{
				currentSum = (mpfr_t *) malloc( sizeof( mpfr_t));
				mpfr_init2( *currentSum, 1000);
				mpfr_set_si( *currentSum, 0, GMP_RNDN); 
			}
			double area=0.0;
				if (isSquarer) // subtract parts that are out of triangle
				{
					int z, q;
					for (z=xT; z<=xB; z++)
						for (q=yT; q<=yB; q++)
							if ((wX-z) > q)
								area += pow(2.0, (double)((wX - z-1) + (wY - q-1)));
				}
			
			mpfr_t s;
			mpfr_init2(s, 1000);
			mpfr_set_ui( s, 2, GMP_RNDN);
			mpfr_pow_si( s, s, power, GMP_RNDN);
			mpfr_mul( *currentSum, *currentSum, s, GMP_RNDN);
			mpfr_sub_d( *currentSum, *currentSum, area, GMP_RNDN);
			mpfr_sub( *fullSum, *fullSum, *currentSum, GMP_RNDN);
			cout << "SOFT DSP Top right = " << xT << ", " << yT << " and bottom left = " << xB << ", " <<yB << endl;
			cout << "This one is :" << mpfr_get_d(*currentSum, GMP_RNDN) << endl;
		}
		return fullSum;
	}
	
	mpfr_t* IntTruncMultiplier::evalMaxValue(int w, int h){
		mpfr_t *l, *r;
		l = (mpfr_t *) malloc( sizeof( mpfr_t) );
		r = (mpfr_t *) malloc( sizeof( mpfr_t) );
		mpfr_init2( *l, 1000);
		mpfr_set_ui( *l, 2, GMP_RNDN);
		mpfr_pow_ui( *l, *l, w, GMP_RNDN);
		mpfr_add_si( *l, *l, -1, GMP_RNDN);
		mpfr_init2( *r, 1000);
		mpfr_set_ui( *r, 2, GMP_RNDN);
		mpfr_pow_ui( *r, *r, h, GMP_RNDN);
		mpfr_add_si( *r, *r, -1, GMP_RNDN);
		mpfr_mul( *l, *l, *r, GMP_RNDN);
		return l; 
	}		

//	void IntTruncMultiplier::emulate(TestCase* tc)
//	{
//		mpz_class svX = tc->getInputValue("X");
//		mpz_class svY = tc->getInputValue("Y");

//		mpz_class svR = svX * svY;

//		tc->addExpectedOutput("R", svR);
//	}

	int IntTruncMultiplier::estimateNrOfDiscardedCols(int k){
		int d = k+1;
		double limit = k;
		double apxError;
		double aux;
		do{
			d--;
			aux = pow(2.0, (double)d);
			aux *= (d-1);
			aux += 2.0;
			aux = log2(aux);
			apxError = ceil(aux);
		}while (limit < apxError);
		return d;
	}

	/** The movement of the DSP blocks with values belonging to their widths and heights still needs to be done. Now it only runs with one type of move on one of the directions, which is not ok for the cases when the DSPs are not squares.
	 */
	void IntTruncMultiplier::tilingAlgorithm(int i, int n,bool repl,int lastMovedDSP)
	{

		if(i==n)
			{
				
				if(repl==true) // if previous DSPs were moved this one needs to recompute all positions 
					{
						//cout<<" Pas 4_1 "<<i<<endl;
						if(replace(globalConfig,i)) // repostioned the DSP
							{
						//		cout<<"Pas 4_0_1 "<<i<<endl;
								compareCost();
								//display(globalConfig);
								rot[i]=false;
								tilingAlgorithm(i,n,false,lastMovedDSP);	
							}
						else // could not reposition the DSP in the bounds of the tiling board
							{
						//		cout<<"Pas 4_5_1 "<<i<<endl;
								rot[i]=false;
								if( lastMovedDSP>=0) // go one level up the backtracking stack
									tilingAlgorithm(lastMovedDSP,n,false, lastMovedDSP);
							}
					}
				else // the last DSP is being moved on the tiling board
					{
						//	cout<<"Pas __1 "<<i<<endl;
						if(move(globalConfig,i)) // successfuly moved the last block
							{
								//cout<<" Pas 1_1 "<<i<<endl;
								compareCost();
								tilingAlgorithm(i,n,repl,i);		//repl should be false
							}
						//~ else
						//~ if(move(globalConfig,i,DSPh,DSPw))
						//~ {
						//~ cout<<" Pas 1_1 "<<i<<endl;
						//~ compareCost();
						//~ tilingAlgorithm(i,n,repl,i);		//repl should be false
						//~ }
 						else // could not find a position for the last block
							{
								if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() ))
									{ // if the DSP was not rotated and is not sqare then roteate it
										//display(globalConfig);
										//~ cout<<" Pas 2_1 "<<i<<endl;
										globalConfig[i]->rotate();
										//display(globalConfig);
										rot[i]=true;
										if(replace(globalConfig,i)) // the DSP could be repositioned
											{
												//display(globalConfig);
												compareCost();
												tilingAlgorithm(i,n,repl,i);		//repl should be false
											}
										else // go to the previous block 
											{
												if(i-1>=0)
													tilingAlgorithm(i-1,n,false,i);
											}
									}
								else // the DSP was either rotated already or is square
									{
										//~ cout<<" Pas 3_1 "<<i<<endl;
										if(i-1>=0)
											tilingAlgorithm(i-1,n,repl,i);		//repl should be false
									}
							}
					}
			}
		else // we are not at the last DSP
			{
				if(repl==true) // the previuos DSPs were successfuly repositioned
					{
					//	cout<<" Pas 4_2 "<<i<<endl;
						if(replace(globalConfig,i)) // the current DSP was successfuly repositioned
							{
								rot[i]=false;
								tilingAlgorithm(i+1,n,repl, lastMovedDSP);
							}
						else // the current DSP could not be repositioned
							{// go to the DSP block that was moved (not repostioned) the last time
								rot[i]=false;
								if( lastMovedDSP>=0) 
									tilingAlgorithm( lastMovedDSP,n,false, lastMovedDSP);
							}
			
		
					}
				else // the folling DSP could not be moved or repositioned 
					{	
						//~ if(i==0)
						//~ display(globalConfig);
						if(move(globalConfig,i)) // the current DSP was successfuly moved
							{
					//			cout<<"Pas 1_2 "<<i<<endl;
								if(i==0){
									cout<<"DSP #1 has made another step!"<<endl;
									//~ display(globalConfig);
									//~ cout<<endl<<endl<<endl;
				
								}
								tilingAlgorithm(i+1,n,true,i);
							}
						//~ if(counterfirst%100==0)
						//~ cout<<counterfirst<<"DSP #1 has made 100 steps!"<<endl;
						//~ display(globalConfig);
						//~ cout<<endl<<endl<<endl;
				
						//~ }
						//~ tilingAlgorithm(i+1,n,true,i);
						//~ }
						else // the current DSP was not moved successfuly
							{
								if(rot[i]==false && (globalConfig[i]->getMaxMultiplierWidth() != globalConfig[i]->getMaxMultiplierHeight() ))
									{// if the DSP was not rotated and is not sqare then roteate it
					//					cout<<" Pas 2_2 "<<i<<endl;
										globalConfig[i]->rotate();
										if(replace(globalConfig,i)) // the current DSP was successfuly repositioned
											{
												rot[i]=true;
												tilingAlgorithm(i+1,n,true,i);
											}
										else // the current DSP was not successfuly repositioned
											{
												if(i-1>=0)
													tilingAlgorithm(i-1,n,repl,i);
											}
									}
								else // the DSP is either square or has been already rotated
									{
					//					cout<<" Pas 3_2 "<<i<<endl;
										if(i-1>=0)
											tilingAlgorithm(i-1,n,repl,i);		//repl should be false
									}
							}
					}
			}

	
	}

	void IntTruncMultiplier::runAlgorithm()
	{	
		int n=vn;
		int m=vm;

		mat = new int*[m];
		for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		
		//if the estimated number of DSPs is grather then 0 then we should run the algorithm
		if(nrDSPs>0){
			tempc= new DSP*[nrDSPs];
			for(int ii=0;ii<nrDSPs;ii++)
				tempc[ii]= new DSP();
			
			target_->getDSPWidths(dspWidth,dspHeight);
			
			/*we will try the algorithm with 2 values of nrDSPs	
			One will be the estimated value(nrDSPs) and the second one will be nrDSPs-1	
			*/
			rot = new bool[nrDSPs];
			for(int i =0;i<nrDSPs;i++)
				rot[i]=false;
		
			//The second
			numberDSP4Overlap=nrDSPs;
			initTiling(globalConfig,nrDSPs);
			cout << "NRDSPs = " << nrDSPs << endl;
			//this will initialize the bestConfig with the first configuration
			bestCost = FLT_MAX ;
			cout<<"Max score is"<<bestCost<<endl;
			//bestConfig = (DSP**)malloc(nrDSPs * sizeof(DSP*));
			bestConfig = new DSP*[nrDSPs];
			for(int i=0;i<nrDSPs;i++)
				bestConfig[i]= new DSP();
			compareCost();
			cout<<"New best score is"<<bestCost<<endl;
			//getchar();
			//display(bestConfig);
			//the best configuration should be consider initially the first one. So the bestConfig parameter will be initialized with global config and hence the bestCost will be initialized with the first cost
		
			//the one
			numberDSP4Overlap=nrDSPs;
			tilingAlgorithm(nrDSPs-1,nrDSPs-1,false,nrDSPs-1);
			bindDSPs(bestConfig);
			vector<SoftDSP*> configSoft;
			if(useLimits==0)
				configSoft = insertSoftDSPs(bestConfig);
			else
				configSoft = insertSoftDSPswithLimits(bestConfig);
			
			
			displayAll(bestConfig, configSoft);
			printConfiguration(bestConfig, configSoft);
			cout<<"Best cost is "<<bestCost<<endl;
			generateVHDLCode(bestConfig, configSoft);
			/*
			globalConfig[2]->setTopRightCorner(2,26);
			globalConfig[2]->setBottomLeftCorner(25,59);
			globalConfig[1]->setTopRightCorner(36,2);
			globalConfig[1]->setBottomLeftCorner(59,35);
			globalConfig[3]->setTopRightCorner(26,36);
			globalConfig[3]->setBottomLeftCorner(59,59);
			bestCost = computeCost(globalConfig);
			display(globalConfig);
			*/
			// After all configurations with the nrDSPs number of DSPs were evaluated then a new search is carryed with one DSP less
			// After the initialization of the new configuration with nrDSPs-1, the cost must be evaluated and confrunted with the best score obtained so far.
		}
		else
		{
			n=vn;
			m=vm;
	
			mat = new int*[m];
			for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
			tempc= new DSP*[0];
			bestConfig = new DSP*[1];
			globalConfig = new DSP*[1];
			tempc[0]=bestConfig[0]=globalConfig[0]=NULL;
			
			bestCost = FLT_MAX ;
			cout<<"Max score is"<<bestCost<<endl;
			compareCost();
			cout<<"New best score is"<<bestCost<<endl;
			display(bestConfig);
		}
	
	}
	
	bool IntTruncMultiplier::compareOccupation(DSP** config)
	{
		int totalSize = wInX * wInY;
		int sizeBest = totalSize;
		int sizeConfig = totalSize;
		//~ cout<<"Total size is "<<totalSize<<endl;
		int c1X,c2X,c1Y,c2Y;
		int nj,ni,njj,nii;
		int ii,jj,i,j;
		int n,m;
		//~ n=wInX + 2* getExtraWidth();
		//~ m=wInY + 2* getExtraHeight();
		n=vm;
		m=vm;
	
		int nmew = n-getExtraWidth();
		int ew = getExtraWidth();
		int mmeh = m - getExtraHeight() ;
		int eh = getExtraHeight();
	
		for(int ti=0;ti<nrDSPs;ti++)
			if(config[ti]!=NULL)
				{
			
			
					config[ti]->getTopRightCorner(c1X,c1Y);
					config[ti]->getBottomLeftCorner(c2X,c2Y);
					c1X=n-c1X-1;
					c2X=n-c2X-1;
					//fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
					j = c2X;
					i = c1Y;
					jj = c1X;
					ii = c2Y;
			
					//cout<<" New coordinates are ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
			
			
						
					if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
						{
							nj=ni=njj=nii=0;
						}
					else
						{
							if( j < getExtraWidth() )
								nj = getExtraWidth() ;
							else
								nj = j;
							if( jj >= n - getExtraWidth() )
								njj = n -getExtraWidth() -1;
							else
								njj = jj;
							
							if( i < getExtraHeight() )
								ni = getExtraHeight() ;
							else
								ni = i;
							if( ii >= m - getExtraHeight() )
								nii = m -getExtraHeight() -1;
							else
								nii = ii;
							//costSlice +=target->getIntMultiplierCost(njj-nj+1,nii-ni+1);
							sizeConfig-=(njj-nj+1)*(nii-ni+1);
							
						}
			
				}
		
		//cout<<"Size config is"<<sizeConfig<<endl;
		//display(config);
		
		
		
		for(int ti=0;ti<nrDSPs;ti++)
			if(bestConfig[ti]!=NULL)
				{
			
			
					bestConfig[ti]->getTopRightCorner(c1X,c1Y);
					bestConfig[ti]->getBottomLeftCorner(c2X,c2Y);
					c1X=n-c1X-1;
					c2X=n-c2X-1;
					//fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
					j = c2X;
					i = c1Y;
					jj = c1X;
					ii = c2Y;
			
					//cout<<" New coordinates are ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
			
			
						
					if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
						{
							nj=ni=njj=nii=0;
						}
					else
						{
							if( j < getExtraWidth() )
								nj = getExtraWidth() ;
							else
								nj = j;
							if( jj >= n - getExtraWidth() )
								njj = n -getExtraWidth() -1;
							else
								njj = jj;
							
							if( i < getExtraHeight() )
								ni = getExtraHeight() ;
							else
								ni = i;
							if( ii >= m - getExtraHeight() )
								nii = m -getExtraHeight() -1;
							else
								nii = ii;
							//costSlice +=target->getIntMultiplierCost(njj-nj+1,nii-ni+1);
							sizeBest-=(njj-nj+1)*(nii-ni+1);
							
						}
			
				}
		
		
		//cout<<"Size best is "<<sizeBest<<endl;
		//display(bestConfig);
		
		
		if(sizeBest >= sizeConfig)		
			return true;
		else
			return false;
	
	}

	void IntTruncMultiplier::fillMatrix(int **&matrix,int lw,int lh,int topleftX,int topleftY,int botomrightX,int botomrightY,int value)
	{
		for(int j=topleftX;j<=botomrightX;j++)
			{
				for(int i=topleftY;i<=botomrightY;i++)
					{
						if(j>-1&&i>-1&&i<lh&&j<lw)
							matrix[i][j]=value;
					}
			}
	}

	void IntTruncMultiplier::display(DSP** config)
	{
	
	
		int **mat;
		int n,m;
		int count=1;
		n=wInX + 2* getExtraWidth();
		m=wInY + 2* getExtraHeight();
		cout<<"real width"<<wInX<<"real height"<<wInY<<endl;
		cout<<"width "<<n<<"height "<<m<<endl;
		mat = new int*[m];
	
		int nmew = n-getExtraWidth();
		int ew = getExtraWidth();
		int mmeh = m - getExtraHeight() ;
		int eh = getExtraHeight();
		int nj,ni,njj,nii;
		for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
	
		count++;
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					{
						if(mat[i][j]==0)
							{
								int ver =0;
								int ii=i,jj=j;
								while(ver<6&&(ii<m-1||jj<n-1))
									{
							
							
										if(ver<3)
											{
								
												if(ver==0||ver==1)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=2;							
													}
							
												if(ver==0||ver==2)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=1;
													}
								
								
							
												for(int k=ii,q=jj;k>i-1&&(ver==0||ver==2);k--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=1;
															else
																ver=3;
															jj--;
														}
									
												for(int k=ii,q=jj;q>j-1&&(ver==0||ver==1);q--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=2;
															else
																ver=3;
															ii--;
														}
								
							
											}
										else
											{
												if(ver==3||ver==5)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=4;
													}
								
								
												if(ver==3||ver==4)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=5;							
													}
							
								
															
												for(int k=ii,q=jj;q>j-1&&(ver==3||ver==4);q--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=5;
															else
																ver=6;
															ii--;
														}
								
												for(int k=ii,q=jj;k>i-1&&(ver==3||ver==5);k--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=4;
															else
																ver=6;
															jj--;
														}
												if(ver==5&&jj==n-1)
													ver=6;
												if(ver==4&&ii==m-1)
													ver=6;
										
															
								
											}
									}
						
						
					
						
						
								if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
									{
										cout<<"Partition number "<<count<<" is totally out of the real multiplication bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
									}
								else
									{
										if( j < getExtraWidth() )
											nj = getExtraWidth() ;
										else
											nj = j;
										if( jj >= n - getExtraWidth() )
											njj = n -getExtraWidth() -1;
										else
											njj = jj;
							
										if( i < getExtraHeight() )
											ni = getExtraHeight() ;
										else
											ni = i;
										if( ii >= m - getExtraHeight() )
											nii = m -getExtraHeight() -1;
										else
											nii = ii;
										cout<<"Partition number "<<count<<" with bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<" has now bounds ("<<nj<<" , "<<ni<<" , "<<njj<<" , "<<nii<<")"<<endl;
									}
						
								cout<<j<<" "<<i<<" "<<jj<<" "<<ii<<endl;
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
						
							}
					}
			
			}
		
		
		
		char af;
		int afi;
		for(int i=0;i<m;i++)
			{
				if(i==getExtraHeight())
					cout<<endl;
			
				for(int j=0;j<n;j++)
					{
						if(j==getExtraWidth())
							cout<<" ";
						if(j==n-getExtraWidth())
							cout<<" ";
				
						if(mat[i][j]<10)
							afi=mat[i][j];
						else
							afi=mat[i][j]+7;
						af=(int)afi+48;
						cout<<af;
					}
				cout<<endl;
				if(i==m-getExtraHeight()-1)
					cout<<endl;
			}
		
		for(int ii=0;ii<m;ii++)
		   delete [] (mat[ii]);
	
		delete[] (mat);
		
		
	
	}

	void IntTruncMultiplier::displayAll(DSP** config, vector<SoftDSP*> softConfig)
	{
	
	
		int **mat;
		int n,m;
		int count=1;
		n=wInX + 2* getExtraWidth();
		m=wInY + 2* getExtraHeight();
		cout<<"real width"<<wInX<<"real height"<<wInY<<endl;
		cout<<"width "<<n<<"height "<<m<<endl;
		mat = new int*[m];
	
		int nmew = n-getExtraWidth();
		int ew = getExtraWidth();
		int mmeh = m - getExtraHeight() ;
		int eh = getExtraHeight();
		int nj,ni,njj,nii;
		for(int i=0;i<m;i++)
			{
				mat[i] = new int [n];
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
			
		for(unsigned i=0;i<softConfig.size();i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				softConfig[i]->getTopRightCorner(c1X,c1Y);
				softConfig[i]->getBottomLeftCorner(c2X,c2Y);
				cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
			
		count++;
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					{
						if(mat[i][j]==0)
							{
								int ver =0;
								int ii=i,jj=j;
								while(ver<6&&(ii<m-1||jj<n-1))
									{
							
							
										if(ver<3)
											{
								
												if(ver==0||ver==1)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=2;							
													}
							
												if(ver==0||ver==2)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=1;
													}
								
								
							
												for(int k=ii,q=jj;k>i-1&&(ver==0||ver==2);k--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=1;
															else
																ver=3;
															jj--;
														}
									
												for(int k=ii,q=jj;q>j-1&&(ver==0||ver==1);q--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=2;
															else
																ver=3;
															ii--;
														}
								
							
											}
										else
											{
												if(ver==3||ver==5)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=4;
													}
								
								
												if(ver==3||ver==4)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=5;							
													}
							
								
															
												for(int k=ii,q=jj;q>j-1&&(ver==3||ver==4);q--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=5;
															else
																ver=6;
															ii--;
														}
								
												for(int k=ii,q=jj;k>i-1&&(ver==3||ver==5);k--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=4;
															else
																ver=6;
															jj--;
														}
												if(ver==5&&jj==n-1)
													ver=6;
												if(ver==4&&ii==m-1)
													ver=6;
										
															
								
											}
									}
						
						
					
						
						
								if( j>= nmew || jj< ew || i >= mmeh || ii < eh)
									{
										cout<<"Partition number "<<count<<" is totally out of the real multiplication bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<endl;
									}
								else
									{
										if( j < getExtraWidth() )
											nj = getExtraWidth() ;
										else
											nj = j;
										if( jj >= n - getExtraWidth() )
											njj = n -getExtraWidth() -1;
										else
											njj = jj;
							
										if( i < getExtraHeight() )
											ni = getExtraHeight() ;
										else
											ni = i;
										if( ii >= m - getExtraHeight() )
											nii = m -getExtraHeight() -1;
										else
											nii = ii;
										cout<<"Partition number "<<count<<" with bounds. ("<<j<<" , "<<i<<" , "<<jj<<" , "<<ii<<")"<<" has now bounds ("<<nj<<" , "<<ni<<" , "<<njj<<" , "<<nii<<")"<<endl;
									}
						
								cout<<j<<" "<<i<<" "<<jj<<" "<<ii<<endl;
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
						
							}
					}
			
			}
		
		
		
		char af;
		int afi;
		for(int i=0;i<m;i++)
			{
				if(i==getExtraHeight())
					cout<<endl;
			
				for(int j=0;j<n;j++)
					{
						if(j==getExtraWidth())
							cout<<" ";
						if(j==n-getExtraWidth())
							cout<<" ";
				
						if(mat[i][j]<10)
							afi=mat[i][j];
						else
							afi=mat[i][j]+7;
						af=(int)afi+48;
						cout<<af;
					}
				cout<<endl;
				if(i==m-getExtraHeight()-1)
					cout<<endl;
			}
		
		for(int ii=0;ii<m;ii++)
		   delete [] (mat[ii]);
	
		delete[] (mat);
		
		
	
	}

	int IntTruncMultiplier::partitionOfGridSlices(DSP** config,int &partitions)
	{
		//~ cout<<"Incepe"<<endl;
		int costSlice=0;
	
		int n,m;
		int count=1;
		n=vn;
		m=vm;
	
		int nmew = vnme;
		int ew = getExtraWidth();
		int mmeh = vmme;
		int eh = getExtraHeight();
		int nj,ni,njj,nii;
		
		
		//~ cout<<"width "<<n<<"height "<<m<<endl;
	
		for(int i=0;i<m;i++)
			{
			
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				//~ cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
		//partitions = count;
		partitions = 0;
	
		//~ cout<<"Partea 2"<<endl;
		
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					{
						if(mat[i][j]==0)
							{
								int ver =0;
								int ii=i,jj=j;
								while(ver<6&&(ii<m-1||jj<n-1))
									{
							
							
										if(ver<3)
											{
								
												if(ver==0||ver==1)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=2;							
													}
							
												if(ver==0||ver==2)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=1;
													}
												//~ if(count==3)
												//~ cout<<"P0  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
								
								
							
												for(int k=ii,q=jj;k>i-1&&(ver==0||ver==2);k--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=1;
															else
																ver=3;
															jj--;
														}
												//~ if(count==3)
												//~ {
												//~ cout<<"P1   ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl; 
												//~ }
									
												for(int k=ii,q=jj;q>j-1&&(ver==0||ver==1);q--)
													if(mat[k][q]!=0)
														{
															if(ver==0)
																ver=2;
															else
																ver=3;
															ii--;
														}
												//~ if(count==3)
												//~ {cout<<"P2  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
												//~ cout<<endl;
												//~ }
							
											}
										else
											{
												if(ver==3||ver==5)
													jj++;
							
												if(jj>n-1)
													{
														jj=n-1;
														ver=4;
													}
								
												//~ if(count==3)
												//~ cout<<"P3  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
								
												if(ver==3||ver==4)
													ii++;
												if(ii>m-1)
													{
														ii=m-1;
														ver=5;							
													}
							
								
												//~ if(count==3)
												//~ cout<<"P3  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;

								
												for(int k=ii,q=jj;q>j-1&&(ver==3||ver==4);q--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=5;
															else
																ver=6;
															ii--;
														}
												//~ if(count==3)
												//~ {
												//~ cout<<"P4   ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl; 
												//~ }
								
												for(int k=ii,q=jj;k>i-1&&(ver==3||ver==5);k--)
													if(mat[k][q]!=0)
														{
															if(ver==3)
																ver=4;
															else
																ver=6;
															jj--;
														}
												if(ver==5&&jj==n-1)
													ver=6;
												if(ver==4&&ii==m-1)
													ver=6;
										
												//~ if(count==3)
												//~ {cout<<"P5  ii:="<<ii<<" jj:="<<jj<<" ver:="<<ver<<endl;
												//~ cout<<endl;
												//~ }
							
								
											}
									}
						
								//~ cout<<count<<endl;
						
					
												
						
								if( j>=nmew || jj< ew || i >= mmeh || ii < eh)
									{
							
									}
								else
									{
										if( j < ew )
											nj = ew ;
										else
											nj = j;
										if( jj >=nmew )
											njj = nmew -1;
										else
											njj = jj;
							
										if( i < eh )
											ni = eh ;
										else
											ni = i;
										if( ii >=mmeh)
											nii = mmeh -1;
										else
											nii = ii;
							
										partitions++;
										//cout << "IntMultiplierCost ("<<nj<<", "<<njj<<") ("<<ni<<", "<<nii<<") cost="<< target_->getIntMultiplierCost(njj-nj+1,nii-ni+1) << endl;
										costSlice += target_->getIntMultiplierCost(njj-nj+1,nii-ni+1);//(njj-nj+1)*(nii-ni+1);
							
							
														
									}
						
						
						
								fillMatrix(mat,n,m,j,i,jj,ii,count);
								count++;
						
							}
					}
			
			}
		
		//de verificat
		
		//cout<<"Count "<<count<<" Partitions "<<partitions<<endl;
		
		//partitions =count -partitions;
		 
		
		//~ char af;
		//~ int afi;
		//~ for(int i=0;i<m;i++)
		//~ {
		//~ for(int j=0;j<n;j++)
		//~ {
		//~ if(mat[i][j]<10)
		//~ afi=mat[i][j];
		//~ else
		//~ afi=mat[i][j]+7;
		//~ af=(int)afi+48;
		//~ cout<<af;
		//~ }
		//~ cout<<endl;
		//~ }
	
		//~ cout<<"gata"<<endl;
	
		//~ for(int ii=0;ii<m;ii++)
		//~ delete[](mat[ii]);
	
		//~ delete[] (mat);
	
		return costSlice;
	}	
	
	int IntTruncMultiplier::bindDSPs4Virtex(DSP** &config)
	{
		int nrOfUsedDSPs=0;
		int dx,dy;
		
		for(int i=0;i<nrDSPs;i++){
			
			config[i]->getTopRightCorner(dx,dy);
			if(config[i]!=NULL )  // && (dx<=vnme && dy<vmme)  // with this condition in this if the algorithm will be encorege to try to compute with less dsps that the user has given
				{
					nrOfUsedDSPs++;
				}
			//countsShift[i]=0;	
			}
			
			
		DSP* ref;
	
		sortDSPs(config);
		
		//display(config);
			
		int itx,ity,jtx,jty,ibx,iby,jbx,jby;
		//int prev=0;
			
		//cout<<endl<<endl;	
		int count;
		
		for(int i=0;i<nrDSPs;i++)
			{
				
				if(config[i]!=NULL)
					{
						
						ref=config[i];
						count=ref->getNrOfPrimitiveDSPs();
						bool ver=true;
						int rw,rh;
						int sa;
						//while(ver==true&&ref->getShiftOut()==NULL && countsShift[prev] <nrOfShifts4Virtex-1)
						while(ver==true&&ref->getShiftOut()==NULL && count <nrOfShifts4Virtex)
							{
								ver=false;
								ref->getTopRightCorner(itx,ity);
								ref->getBottomLeftCorner(ibx,iby);
								rw=ref->getMaxMultiplierWidth();
								rh=ref->getMaxMultiplierHeight();
					
								for(int j=0;j<nrDSPs&&ver==false;j++)
									{
										
										if(config[j]!=NULL &&j!=i && (count+ config[j]->getNrOfPrimitiveDSPs()<=nrOfShifts4Virtex))
											{
												config[j]->getTopRightCorner(jtx,jty);
												if((jtx<=vnme && jty<vmme))
												{
												
												sa = config[j]->getShiftAmount();
												//~ cout<<"Now considering taking(in left) dsp nr. "<<i<<" with tx:="<<itx<<" ty:="<<ity<<" bx:="<<ibx<<"by:="<<iby<<" with dsp nr. "<<j<<" with tx:="<<jtx<<" ty:="<<jty<<endl;
												//config[j]->getBottomLeftCorner(jbx,jby);
												
												if(rw!=34 && rh!=34)
												{
												if(jtx==ibx+1&&jty==ity&&rw==sa&&config[j]->getShiftIn()==NULL)
													{
														//~ cout<<"DSP #"<<i<<" bind with DSP# "<<j<<endl;
														ver=true;
														ref->setShiftOut(config[j]);
														config[j]->setShiftIn(ref);
														nrOfUsedDSPs--;
														ref=config[j];
														count+=ref->getNrOfPrimitiveDSPs();
														//~ countsShift[prev]++;
														//~ countsShift[j] = countsShift[prev];
														//~ prev = j;								
													}
												}
												else
												{
													config[j]->getBottomLeftCorner(jbx,jby);
													//~ cout<<config[i]->getMaxMultiplierHeight()<<" "<<rh<<" ";
													//~ cout<<(config[i]->getMaxMultiplierHeight()% rh)<<endl;

													//if( jtx==ibx+1 && sa!=0 &&   config[j]->getMaxMultiplierWidth() %sa==0 && ( (config[j]->getMaxMultiplierWidth() == 34 && jty==ity)   || ( config[j]->getMaxMultiplierWidth()==17 && jty==ity+sa   )  ))
					
													
													
													if( (jtx==ibx+1) && sa!=0 && ( (config[j]->getMaxMultiplierWidth() == 34 && jby ==iby ) || (config[j]->getMaxMultiplierWidth() == 17  && (iby + sa == jby)) ) )
													{
														//~ cout<<" case 1_2 DSP #"<<i<<" bind with DSP# "<<j<<endl;
														ver=true;
														ref->setShiftOut(config[j]);
														config[j]->setShiftIn(ref);
														nrOfUsedDSPs--;
														ref=config[j];
														count+=ref->getNrOfPrimitiveDSPs();
													}
													
												}
											}
											}
									}
					
								for(int j=0;j<nrDSPs&&ver==false;j++)
									{
										if(config[j]!=NULL &&j!=i&& count+ config[j]->getNrOfPrimitiveDSPs()<=nrOfShifts4Virtex)
											{
												config[j]->getTopRightCorner(jtx,jty);
												if((jtx<=vnme && jty<vmme))
												{
												sa = config[j]->getShiftAmount();
												//~ cout<<"Now considering taking(down) dsp nr. "<<i<<" with tx:="<<itx<<" ty:="<<ity<<" bx:="<<ibx<<"by:="<<iby<<" with dsp nr. "<<j<<" with tx:="<<jtx<<" ty:="<<jty<<endl;
												//config[j]->getBottomLeftCorner(jbx,jby);
												if(rw!=34 && rh!=34)
												{
												if(iby+1==jty&&itx==jtx&&rh==sa&&config[j]->getShiftIn()==NULL)
													{
														//~ cout<<"DSP #"<<i<<" bind with DSP# "<<j<<endl;
														ver=true;
														ref->setShiftOut(config[j]);
														config[j]->setShiftIn(ref);
														nrOfUsedDSPs--;
														ref=config[j];								
														count+=ref->getNrOfPrimitiveDSPs();
														//~ countsShift[prev]++;
														//~ countsShift[j] = countsShift[prev];
														//~ prev = j;
													}
												}
												else
												{
													//~ cout<<config[i]->getMaxMultiplierWidth()<<" "<<rw<<endl;
													//~ cout<<(config[i]->getMaxMultiplierWidth()% rw)<<endl;

													//&&  (config[j]->getMaxMultiplierWidth()% rw==0)
													config[j]->getTopRightCorner(jtx,jty);
													config[j]->getBottomLeftCorner(jbx,jby);
													//if( iby+1==jty &&sa!=0&& config[j]->getMaxMultiplierHeight()% sa==0 && ( (config[j]->getMaxMultiplierHeight() == 34 && jtx==itx      )   || ( config[j]->getMaxMultiplierHeight()==17 && jtx==itx+sa)  ))
													
													if(   iby+1==jty &&sa!=0&& (  ((config[j]->getMaxMultiplierHeight() == 34  &&  jbx ==ibx) ||  (config[j]->getMaxMultiplierHeight() == 17  &&  (jbx ==ibx+sa  )))     ) )
													
													{
														//~ cout<<"case 2_2 DSP #"<<i<<" bind with DSP# "<<j<<endl;
														ver=true;
														ref->setShiftOut(config[j]);
														config[j]->setShiftIn(ref);
														nrOfUsedDSPs--;
														ref=config[j];								
														count+=ref->getNrOfPrimitiveDSPs();
													}
												}
											}
											}						
									}					
							}
				
					}
			
			}
	
		//~ cout<<" nr of dsps after bind "<<nrOfUsedDSPs<<endl;
		return nrOfUsedDSPs;
	
	}
	
	void IntTruncMultiplier::sortDSPs(DSP** &config)
	{
		int ix,iy,jx,jy;
		DSP* temp;
		for(int i=0;i<nrDSPs-1;i++)
			{		
				for(int j=i+1;j<nrDSPs;j++)
					{
						config[i]->getTopRightCorner(ix,iy);
						config[j]->getTopRightCorner(jx,jy);
						if(iy>jy)
							{
								temp=config[i];
								config[i]=config[j];
								config[j]=temp;				
							}
						else
							if(iy==jy)
								{
									if(ix>jx)
										{
											temp=config[i];
											config[i]=config[j];
											config[j]=temp;						
										}
								}
					}
			}
	
	}

	int IntTruncMultiplier::bindDSPs(DSP** &config)
	{
		if  ( (target_->getID() == "Virtex4") ||
		      (target_->getID() == "Virtex5") ||
		      (target_->getID() == "Spartan3"))  // then the target is a Xilinx FPGA
			{
				return bindDSPs4Virtex(config);
			}
		else // the target is Stratix
			{
				return bindDSPs4Stratix(config);
			}
	}

	void IntTruncMultiplier::compareCost()
	{
		//~ cout<<"Inta la cost"<<endl;
	
		//~ DSP** tempc;
		
		//~ tempc= new DSP*[nrDSPs];
		//~ for(int ii=0;ii<nrDSPs;ii++)
		//~ tempc[ii]= new DSP();
		
		//display(globalConfig);
		//getchar();
		//memcpy(tempc,globalConfig,sizeof(DSP*) *nrDSPs );
		for(int ii=0;ii<nrDSPs;ii++)
			memcpy(tempc[ii],globalConfig[ii],sizeof(DSP) );
	
		//display(tempc);
		//getchar();
		//~ cout<<"intra la display cost"<<endl;
	
		float temp = computeCost(tempc);
	
		//cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
	
		if(temp < bestCost)
			{
				//cout<<"Costul e mai bun la cel curent!Schimba"<<endl;
				
				//cout<<"Interchange! Score for current is"<<temp<<" and current best is"<<bestCost<<endl;
				//getchar();
				//~ int c1X,c2X,c1Y,c2Y;
				//~ int n=wInX + 2* getExtraWidth();
				//~ tempc[0]->getTopRightCorner(c1X,c1Y);
				//~ tempc[0]->getBottomLeftCorner(c2X,c2Y);
				//~ cout<<"DSP #"<<1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				//~ c1X=n-c1X-1;
				//~ c2X=n-c2X-1;
				//~ //cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
				//~ cout<<"matrix DSP #"<<1<<"has topleft ("<<c2X<<","<<c1Y<<") and botomright ("<<c1X<<","<<c2Y<<")"<<endl;
		
		
		
				bestCost=temp;
				//memcpy(bestConfig,tempc,sizeof(DSP*) *nrDSPs );	
				for(int ii=0;ii<nrDSPs;ii++)
					memcpy(bestConfig[ii],tempc[ii],sizeof(DSP) );
				//display(bestConfig);
			}
		else
			if(temp == bestCost )
				{
					//cout<<"Cost egal!!!"<<endl;
					//cout<<"Rezult compare is"<<compareOccupation(tempc)<<endl;
					//~ cout<<"score temp is"<<temp<<" and current best is"<<bestCost<<endl;
					//display(bestConfig);
					if(compareOccupation(tempc)==true)
						{
							//cout<<"Interchange for equal cost. Now best has cost "<<temp<<endl;
							
							bestCost=temp;
			
							for(int ii=0;ii<nrDSPs;ii++)
								memcpy(bestConfig[ii],tempc[ii],sizeof(DSP) );
							//	display(bestConfig);
						}
				}
	
	
		//~ for(int ii=0;ii<nrDSPs;ii++)
		//~ free(tempc[ii]);
	
		//~ delete[] (tempc);
	
	}

	float IntTruncMultiplier::computeCost(DSP** &config)
	{
	
		
		float acc=0.0;
		
		//costLUT = ( (1.0+scale) - scale * (1-ratio) ) /  ((float)100);
	
		//cout<<"Cost of a DSP is "<<costDSP<<endl<<"Cost of a Slice is "<<costLUT<<endl;
	
		int nrOfUsedDSPs=0;
		
		for(int i=0;i<nrDSPs;i++)
			if(config[i]!=NULL)
				{
					int nrPrimitiveDSPs = config[i]->getNrOfPrimitiveDSPs();
					//acc += (costDSP*float(nrPrimitiveDSPs));
					nrOfUsedDSPs += nrPrimitiveDSPs;
				}
		
		//cout<<"Number of used DSP blocks is "<<nrOfUsedDSPs<<endl;
		vector<SoftDSP*> configSoft;
		if(useLimits==0)
			configSoft = insertSoftDSPs(config);
		else
			configSoft = insertSoftDSPswithLimits(config);
		
		int partitions = (int)configSoft.size();
		float LUTs4Multiplication =  0;
		int stx, sty, sbx, sby;
		
		for (int i=0; i<partitions; i++)
		{
			configSoft[i]->getTopRightCorner(stx,sty);
			configSoft[i]->getBottomLeftCorner(sbx,sby); 
			LUTs4Multiplication += target_->getIntMultiplierCost(sbx-stx+1, sby-sty+1);
		}
		
		//cout<<"Number of slices 4 multiplication of the rest is "<<LUTs4Multiplication<<endl;
		
		acc =((float)nrOfUsedDSPs)*costDSP + costLUT * LUTs4Multiplication;
		
		//cout << "Accum = " << acc << endl;
		//~ cout<<"Number of partitions for LUTs is "<<partitions<<endl;
		nrOfUsedDSPs = bindDSPs(config);
		//~ cout<<"Number of operands coming from DSPs is "<<nrOfUsedDSPs<<endl;
		
	
		float LUTs4NAdder=((float)target_->getIntNAdderCost(wInX + wInY,nrOfUsedDSPs+partitions) );
		//float LUTs4NAdder=((float)200);
				
				
	
		//~ cout<<"LUTs used for last "<<nrOfUsedDSPs+partitions<<" adder are"<<LUTs4NAdder<<endl;
		
		acc +=  LUTs4NAdder* costLUT;	
		//cout << "Accum = " << acc << endl;
		
		//~ Substracting the cost of different additions that can be done in DSPs(Virtex) or the cost of a DSP if they can be combined(Altera)
		
		
		return acc;
			
	}

	int IntTruncMultiplier::estimateDSPs()
	{
		if (ratio>1)
			ratio =1;
		float t1,t2, t3, t4;
		int Xd, Yd; //the dimension of the multiplier on X and on Y
		//int multX1, multY1, mult1, multX2, multY2, mult2;
		bool fitMultiplicaication = false;
		//int tempDSP;
		target_->getDSPWidths(Xd,Yd);
		int wInX = vnme-getExtraWidth();
		int wInY = vmme-getExtraHeight();
		int maxDSP, mDSP = target_->getNumberOfDSPs();
		int wInXt = wInX;
		int wInYt = wInY;
		if (wInY > wInX)
		{
			wInYt = wInX;
			wInXt = wInY;
		}
		//~ tempDSP =    int(ceil( t1) * ceil(t2) );
		
		//t1 = wInX*wInY;
		//t2 = Xd*Yd; 
		//tempDSP =    int(ceil(   ((float) t1) / ((float) t2) ));
		t1 = ((float) wInXt) / ((float) Xd);
		t2 = ((float) wInYt) / ((float) Yd);
		t3 = ((float) wInXt) / ((float) Yd);
		t4 = ((float) wInYt) / ((float) Xd);
		if (t1 < 1) // DSP width > multiplication width
		{
			if (t3 < 1) // DSP height > multiplication width
				maxDSP = int (ceil(t4));
			else
				maxDSP = int (ceil(t2));
		}
		else
		{
			if (t2 < 1) // DSP height > multiplication height
				maxDSP = int (ceil(t1));
			else if (t4 < 1) // DSP width > multiplication height
				maxDSP = int (ceil(t3));
			else // none of the above
			{
				int rw = wInXt % Xd;
				int rh = wInXt % Yd;
				
				if ((rw == 0) || (rh == 0))// divide by width or height
					maxDSP = int (ceil(t1)*ceil(t2));
				else if ((rh-Xd) <= 0) // can pad with width
				{
					if ((rw-Yd) <= 0) // can pad with height
					{	
						if ((rw-Yd) >= (rh-Xd)) // pad with height
							maxDSP = int (floor(t1)*ceil(t2)+ceil(t4));
						else // pad with width
							maxDSP = int (floor(t3)*ceil(t4)+ceil(t2));
					}
					else // can't pad with height
						maxDSP = int (floor(t3)*ceil(t4)+ceil(t2));
				}
				else if ((rw-Yd) <= 0) // can pad with height
					maxDSP = int (floor(t1)*ceil(t2)+ceil(t4));
				else
					maxDSP = int (ceil(t1)*ceil(t2));
			}
		}
		int limitDSP = maxDSP;
		/* After estimating the maximal number of DSPs for the entire 
		 * board subtract the number of DSPs not needed to cover the 
		 * part of the tiling board that is of no interest */
		int areaTrunc = 0;
		int areaCommon = 0;
		int areaSq = 0;
		
		// Compute area not needed after truncation
		for(int i=0,j=(wInY-truncationOffset);i<wInX&&j<wInY&&i<truncationOffset;i++,j++)	
			{
				for(int k=j;k<wInY;k++)
				{
					if(k>=0)
						areaTrunc++;
					if(k>(wInX-i))
						areaCommon++;
				}
			}
			
		if(isSquarer)	
		{
			// Area above second diagonal minus max surface of standing out ridges
			areaSq = wInX*wInY/2 - (wInXt/Xd+1)*(Xd*Yd/2);
			maxDSP -= floor((double)areaSq/(Xd*Yd));
			areaTrunc -= areaCommon;
		}
		// (truncationOffset/Xd)*(Xd*Yd/2)
		areaTrunc -= (truncationOffset*Yd/2);
		if (areaTrunc < 0)
			areaTrunc = 0;
			
		maxDSP -= floor((double)areaTrunc/(Xd*Yd));
		
		if(maxDSP > mDSP)
		{
			fitMultiplicaication = true;
			//maxDSP = tempDSP;
			maxDSP = mDSP; //set the maximum number of DSPs to the multiplication size
		}
			
		if (ratio == 1){
			if (! fitMultiplicaication){
				REPORT(INFO, "Warning!!! The number of existing DSPs on this FPGA is not enough to cover the whole multiplication!");
			}else{
				REPORT(INFO, "Warning!!! The minimum number of DSP that is neccessary to cover the whole addition will be used!");
			}
			
			return maxDSP;
		}else{	
			//cout<<" Eq Slice DSP "<<target_->getEquivalenceSliceDSP()<<endl;
			//cout<<" IntM cost "<<target_->getIntMultiplierCost(wInX,wInY)<<endl;
			//float temp = ( float(target_->getIntMultiplierCost(wInX,wInY)) * ratio)  /   ((1.-ratio)*float(target_->getEquivalenceSliceDSP())) ;
			float scaleDSPs=1.0;
			float temp;
			if(isSquarer || (truncationOffset>0))
			{
				int trCornerMultX = truncationOffset/2;
				int trCornerMultY = truncationOffset - trCornerMultX;
				float intMultCost = float(target_->getIntMultiplierCost(wInX,wInY)-target_->getIntMultiplierCost(trCornerMultX, trCornerMultY));
				temp = (intMultCost * ratio * maxDSP)  /   (float(target_->getEquivalenceSliceDSP()*limitDSP)) ;
				cout<<"temp= "<<temp<<endl;
			}
			else
			{
				if(maxDSP>4)
					scaleDSPs= ((float)maxDSP) / ((float)maxDSP-2);
				temp = ( float(target_->getIntMultiplierCost(wInX,wInY)) * ratio* scaleDSPs)  /   (float(target_->getEquivalenceSliceDSP())) ;
				//float temp = ( float(target_->getIntMultiplierCost(wInX,wInY)) )  /   ((1.-ratio)*float(target_->getEquivalenceSliceDSP())) ;
				cout<<"val calc "<<temp<<endl;
			}
			int i_tmp = int(ceil(temp));
			cout<<" rounded "<<i_tmp<<endl;
	
			if(i_tmp > maxDSP){
				if (fitMultiplicaication){
					REPORT(INFO, "Warning!!! The number of estimated DSPs with respect with this ratio of preference is grather then the needed number of DSPs to perform this multiplication!");
				}else{
					REPORT(INFO, "Warning!!! The number of estimated DSPs with respect with this ratio of preference is grather then the total number of DSPs that exist on this board!");
				}
				
				i_tmp = maxDSP;
				//cout<<" final val is(the max dsps) "<<i_tmp<<endl;
			}
	
			return i_tmp ;
		}
	}

	int  IntTruncMultiplier::getExtraHeight()
	{	
		int x,y;	
		target_->getDSPWidths(x,  y);
		float temp = ratio * 0.75 * ((float) y);
		return ((int)temp);
		//return 4;

	}
	
	int  IntTruncMultiplier::getExtraWidth()
	{
		int x,y;	
		target_->getDSPWidths(x,y);
		float temp = ratio * 0.75 * ((float) x);
		return ((int)temp);
		//return 4;

	}
	
	bool IntTruncMultiplier::isValidPosition(int x, int y)
	{
		if (x>=vn || y>=vm || x<0 || y<0) // then not in tiling grid bounds
		{
			cout << "then not in tiling grid bounds" << endl;
			return false;
		}
			
		if (isSquarer && (x > y)) // then the position is above the secondary diagonal
		{
			cout << "then the position: ("<<x<<", "<<y<<") is above the secondary diagonal" << endl;
			return false;
		}	
		if ((truncationOffset>0) && ((x-truncationOffset) > y)) // then the position is above the truncation boundary
		{
			cout << "then the position is above the truncation boundary" << endl;
			return false;
		}
			
		return true;
	}

	bool IntTruncMultiplier::checkOverlap(DSP** config, int index)
	{			
		int x,y,w,h;
		h = config[index]->getMaxMultiplierHeight();
		w = config[index]->getMaxMultiplierWidth();
		int poscur = config[index]->getCurrentPosition();
		int posav = config[index]->getAvailablePositions();
		int minY,minX;
		config[index]->getTopRightCorner(minX,minY);
		//~ cout<<index<<" "<< minX<<" "<<minY<<" "<<poscur<<" "<<posav<<" ";
		minX+=w; 
		for(;poscur<posav;poscur++)
		{
			if(minY>config[index]->Ypositions[poscur])
			{
				minY=config[index]->Ypositions[poscur];
				minX=config[index]->Xpositions[poscur];
			}
		}
		//~ cout<<index<<" "<< minX<<" "<<minY<<" ";
		
		config[index]->getBottomLeftCorner(x,y);
		
		
		
		long area = (vn - minX+6) * (vm -minY+6) + w * (vm- y);

		//~ cout<<" area "<<area<<" dspArea "<<dsplimit<<" ";
		int dsplimittemp = (int)ceil( ((double)area) / dsplimit );
		//~ cout<<" limit "<<dsplimittemp<<" nrrest "<<(numberDSP4Overlap -index-1)<<endl;
		//I have added +1 to the number of dsps that could be place because when we work with combined dsps there exists the case when an are of a dsp is not added
		if( dsplimittemp+1 < (numberDSP4Overlap -index-1))

			return true;
		
		//return false;
		
		//not used now
		
		int xtr1, ytr1, xbl1, ybl1, xtr2, ytr2, xbl2, ybl2;
		//int nrdsp = sizeof(config)/sizeof(DSP*);
		config[index]->getTopRightCorner(xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
	

		
		bool a1 ;
		bool a2 ;
		bool a3 ;
		bool a4 ;
		bool a5 ;
		bool a6 ;
				
		bool b1 ;
		bool b2 ;
		bool b3 ;
		bool b4 ;
		bool b5 ;
		bool b6 ;
		
	
		if(verbose)
			cout << tab << tab << "checkOverlap: ref is block #" << index << ". Top-right is at (" << xtr1 << ", " << ytr1 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
	
		for (int i=0; i<index; i++)
			if (config[i] != NULL)
				{
					config[i]->getTopRightCorner(xtr2, ytr2);
					config[i]->getBottomLeftCorner(xbl2, ybl2);
					//cout<<index<<" "<<i<<" "<<xbl1<<" "<<xtr2<<endl;
					if (((xbl1 < xbl2) && (ytr2 > ybl1)) || 	// config[index] is above and to the right of config[i]
						  ((xbl1 < xtr2)))// && (ybl1 < ybl2))) 							// config[index] is to the right of config[i]
						return true;
			
				
					if(verbose)
						cout << tab << tab << "checkOverlap: comparing with block #" << i << ". Top-right is at (" << xtr2 << ", " << ytr2 << ") and Bottom-right is at (" << xbl1 << ", " << ybl1 << ")" << endl;
			
					//~ if (((xtr2 <= xbl1) && (ytr2 <= ybl1) && (xtr2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the upper and/or right part(s) of config[i]
					//~ ((xbl2 <= xbl1) && (ybl2 <= ybl1) && (xbl2 >= xtr1) && (ybl2 >= ytr1)) || // config[index] overlaps the bottom and/or left part(s) of config[i]
					//~ ((xtr2 <= xbl1) && (ybl2 <= ybl1) && (xtr2 >= xtr1) && (ybl2 >= ytr1)) || // config[index] overlaps the upper and/or left part(s) of config[i]
					//~ ((xbl2 <= xbl1) && (ytr2 <= ybl1) && (xbl2 >= xtr1) && (ytr2 >= ytr1)) || // config[index] overlaps the bottom and/or right part(s) of config[i]
					//~ ((xbl2 >= xbl1) && (ybl2 <= ybl1) && (ytr2 >= ytr1) && (xtr2 <= xtr1)) || // config[index] overlaps the center part of config[i]
				
					//~ ((xtr1 <= xbl2) && (ytr1 <= ybl2) && (xtr1 >= xtr2) && (ytr1 >= ytr2)) || // config[i] overlaps the upper and/or right part(s) of config[index]
					//~ ((xbl1 <= xbl2) && (ybl1 <= ybl2) && (xbl1 >= xtr2) && (ybl1 >= ytr2)) || // config[i] overlaps the bottom and/or left part(s) of config[index]
					//~ ((xtr1 <= xbl2) && (ybl1 <= ybl2) && (xtr1 >= xtr2) && (ybl1 >= ytr2)) || // config[i] overlaps the upper and/or left part(s) of config[index]
					//~ ((xbl1 <= xbl2) && (ytr1 <= ybl2) && (xbl1 >= xtr2) && (ytr1 >= ytr2)) || // config[i] overlaps the bottom and/or right part(s) of config[index]
					//~ ((xbl1 >= xbl2) && (ybl1 <= ybl2) && (ytr1 >= ytr2) && (xtr1 <= xtr2))    // config[i] overlaps the center part of config[index]
					//~ )
					//~ return true;
					
					
					
					// the optimisation of the above if
					a1 = (xtr2 <= xbl1);
					a2 = (xtr2 >= xtr1);
					a3 = (xbl2 <= xbl1);
					a4 = (xbl2 >= xtr1);
					a5 = (xbl2 >= xbl1);
					a6 = (xtr1 >= xtr2);
				
					b1 = (ytr2 <= ybl1);
					b2 = (ytr2 >= ytr1);
					b3 = (ybl2 <= ybl1);
					b4 = (ybl2 >= ytr1);
					b5 = (ytr1 >= ytr2);
					b6 = (ybl1 <= ybl2);
				
					if ((((a1 && a2)||(a3 && a4)) && ((b1 && b2)||(b3 && b4))) || 
					    (((a4 && a6)||(a5 && a1)) && ((b6 && b1)||(b4 && b5))) || 
					    (((a5 && b3) && ( b2 && a6)) || ((a3 && b6) && (b5 && a2))))
						return true;
					
				
			
			
				}	
		if(verbose)
			cout << tab << tab << "checkOverlap: return false" << endl;	
		return false;
	}

	/**
		There is one case that is not resolved w=yet. For DSP with widths different then height the algorithm should move the dsp with both values
	*/

	bool IntTruncMultiplier::move(DSP** config, int index)
	{
		int xtr1, ytr1, xbl1, ybl1;
		int w, h;
		//target->getDSPWidths(w,h);
		w= config[index]->getMaxMultiplierWidth();
		h= config[index]->getMaxMultiplierHeight();
	
		config[index]->getTopRightCorner(xtr1, ytr1);
		config[index]->getBottomLeftCorner(xbl1, ybl1);
	
		if(verbose)
			cout << tab << "replace : DSP #" << index << " width is " << w << ", height is " << h << endl; 
		/*
		 * Initialized in the constructor:
		 * vn = wInX+2*getExtraWidth(); width of the tiling grid including extensions
		 * vm = wInY+2*getExtraHeight(); height of the tiling grid including extensions
		 * */
		//int exh = getExtraHeight();
		//int exw = getExtraWidth();
		int pos; // index for list of positions of a DSP
		
		
		if(index==0) // the first DSP block can move freely on the tiling grid
		{
			return false;
			//if ((xtr1 > 0) && (ytr1 > 0) && (xbl1 < vn-1) && (ybl1 < vm-1))
					{// then the DSP block is placed outside the bounds 		
	
						//do{
							// move down one unit
							/*
							ytr1++;
							ybl1++;
			
							if (ytr1 > exh) // the DSP block has reached the bottom limit of the tiling grid
								{
									// move to top of grid and one unit to the left 
									xtr1++;
									xbl1++;
									
									ytr1 = exh -1; //0
									ybl1 = ytr1 + h-1;
										
									if (xtr1 > exw) // the DSP block has reached the left limit of the tiling grid
										return false;
								}						
							config[index]->setTopRightCorner(xtr1, ytr1);
							config[index]->setBottomLeftCorner(xbl1, ybl1);
							*/ 
						//}while (checkOverlap(config, index));
					}
		}
		else // all DSP blocks except the first one can move only in fixed positions
		{
						do{
							// move to next position
							pos = config[index]->pop();
							if(pos >= 0)
							{
									ytr1 = config[index]->Ypositions[pos];
									ybl1 = ytr1 + h-1;
									xtr1 = config[index]->Xpositions[pos];
									xbl1 = xtr1 + w-1;
							}
							else
								return false;
								
							//if ((ytr1 > gh) || (xtr1 > gw) || (ybl1 < exh) || (xbl1 < exw)) // the DSP block is out of the tiling grid
							//	continue;
														
							config[index]->setTopRightCorner(xtr1, ytr1);
							config[index]->setBottomLeftCorner(xbl1, ybl1);
						}while (checkOverlap(config, index));
		}
		return true;
	}

	bool IntTruncMultiplier::replace(DSP** config, int index)
	{
		//~ cout<<"DSP nr "<<index<<endl;
		int xtr1, ytr1, xbl1, ybl1;
		int w, h;
		string targetID = target_->getID();
		
		//target->getDSPWidths(w,h);
		w= config[index]->getMaxMultiplierWidth();
		h= config[index]->getMaxMultiplierHeight();
		config[index]->setPosition(0);
		
		if (index > 1)
		{// take all positions from the previous DSP
			int curpos =config[index-1]->getCurrentPosition();
			int avpos = config[index-1]->getAvailablePositions();
				memcpy(config[index]->Xpositions, config[index-1]->Xpositions + curpos, sizeof(int)*(avpos - curpos));	
				memcpy(config[index]->Ypositions, config[index-1]->Ypositions + curpos, sizeof(int)*(avpos - curpos));	
				config[index]->setPosition((avpos - curpos));
		}
		
		if (index > 0)
		{
			w= config[index-1]->getMaxMultiplierWidth();
			h= config[index-1]->getMaxMultiplierHeight();
			//int dif = abs(h-w);
			//int maxd = h;
			int mind = w;
			
			if (w > h)
			{
				mind = h;
			}
			int mindX = mind;
			int mindY = mind;
			

			if (targetID == "Virtex5")
			{ // align bottom-left corner of current with X-possition of previous to catch ideal case
				mindX = abs(w-h);
				mindY = h;
			}
			else if ((targetID == "StratixIV") || (targetID == "StratixIV"))
			{ // align on diagonal
				mindX = -w;
				mindY = h;	
			}
			//int positionDisplacementX[] = {0, dif, mind, maxd, w, w, w, w};
			int positionDisplacementX[] = {0, w, mindX};
			//int positionDisplacementY[] = {h, h, h, h, 0, dif, mind, maxd};
			int positionDisplacementY[] = {h, 0, mindY};

			int x,y,x1,y1, pos;
			config[index-1]->getTopRightCorner(x, y);
			
			//* Loop unrolling
			bool extraPosition = ((w!=h) || ((targetID == "StratixIV") || (targetID == "StratixIV")));
			for (int i=0; i<3; i++)
			{
				x1 =x+positionDisplacementX[i];
				y1 =y+positionDisplacementY[i];
				//if(  (x1<vn&&y1<vm)||(x1>vn&&y1>vm)  ) //allows dsp to get out of the board in the diagonal of the bottom left corner
				if (isValidPosition(x1, y1)) 
					if((i!=2) || extraPosition)
					config[index]->push(x1, y1);
				
			}
			/*
			cout<<endl<<"index "<<index<<" ";
			config[index]->resetPosition();
			do
			{
				pos = config[index]->pop();
				if(pos>=0)
				cout<<" ("<<config[index]->Xpositions[pos]<<" , "<<config[index]->Ypositions[pos]<<")";	
			}while(pos>=0);
			cout<<endl;
			getchar();
			*/
			w= config[index]->getMaxMultiplierWidth();
			h= config[index]->getMaxMultiplierHeight();

			config[index]->resetPosition();
			
			do{// go to next position in list
				pos = config[index]->pop();
				if(pos >= 0)
				{
					ytr1 = config[index]->Ypositions[pos];
					ybl1 = ytr1 + h-1;
					xtr1 = config[index]->Xpositions[pos];
					xbl1 = xtr1 + w-1;
				}
				else
				{
					config[index]->setTopRightCorner(vn, vm);
					config[index]->setBottomLeftCorner(vn+w, vm+h);
					return false;
				}
				
				//~ cout<<index<<" (X,Y)="<<xtr1<<" "<<ytr1<<endl;
					
				//if ((ytr1 > gh) && (xtr1 > gw) && (ybl1 < exh) && (xbl1 < exw)) // the DSP block is out of the tiling grid
				//	continue;
											
				config[index]->setTopRightCorner(xtr1, ytr1);
				config[index]->setBottomLeftCorner(xbl1, ybl1);
			}while (checkOverlap(config, index));
				
			return true;
		}
		else
		{	
		// TODO if first DSP block	
		//xtr1 = ytr1 = 0;
		//xbl1 = w-1;
		//ybl1 = h-1;
		
		//cout<<"index"<<endl;
			
		int exh = getExtraHeight();
		int exw = getExtraWidth();
		
		xtr1 = exw;
		ytr1 = exh;	
		ybl1 = ytr1 + h-1;
		xbl1 = xtr1 + w-1;
	
		config[index]->setTopRightCorner(xtr1, ytr1);
		config[index]->setBottomLeftCorner(xbl1, ybl1);
		
		if(verbose)
			cout << tab << "replace : DSP width is " << w << ", height is " << h << endl; 
		// try yo place the DSP block inside the extended region of the tiling grid
		}
		return true;
	}

	void IntTruncMultiplier::initTiling(DSP** &config, int dspCount)
	{
		int w,h; 
		target_->getDSPWidths(w, h);
		dsplimit = w*h;
		config = new DSP*[nrDSPs];
		//nrOfShifts4Virtex=4;
		//countsShift = new int[nrDSPs];
		for (int i=0; i<nrDSPs; i++)
		{
			config[i] = NULL;
		}
		for (int i=0; i<dspCount; i++)
		{
			if(verbose)
				cout << "initTiling : iteration #" << i << endl; 
			config[i] = target_->createDSP();						
			config[i]->setNrOfPrimitiveDSPs(1);
			
			config[i]->allocatePositions(3*i); // each DSP offers 8 positions
			if(!replace(config, i))
			{
				w=config[i]->getMaxMultiplierWidth();
				h= config[i]->getMaxMultiplierHeight();
				config[i]->setTopRightCorner(vn, vm);
				config[i]->setBottomLeftCorner(vn+w, vm+h);
			}
		}
		
	}
	
	void IntTruncMultiplier::initTiling2(DSP** &config, int dspCount)
	{
		//nrOfShifts4Virtex =2; 
		int w, h;
		target_->getDSPWidths(w, h);
		int min = h;
		int mw,mh;
		mw = vn - 2*getExtraWidth();
		mh = vm - 2*getExtraHeight();
		if(w < h)
		{
			min = w;
			w = h;
			h = min;
		}	
		
		// allocate the maximum number of DSP objects
		config = new DSP*[nrDSPs];
		for (int i=0; i<nrDSPs; i++)
		{	
			config[i] = NULL;
		}
		
		/* NOTE: In case the multiplication is narrower than twice the minimum dimension of a DSP 
		 * then we do not form paris of multipliers
		 */
		if ((vm < 2*min) || (vn < 2*min))
		{
			for (int i=0; i<dspCount; i++)
			{
				config[i] = target_->createDSP();
				config[i]->setNrOfPrimitiveDSPs(1);
				config[i]->allocatePositions(3*i); // each DSP offers 3 positions
				if(!replace(config, i))
				{
					w=config[i]->getMaxMultiplierWidth();
					h= config[i]->getMaxMultiplierHeight();
					config[i]->setTopRightCorner(vn, vm);
					config[i]->setBottomLeftCorner(vn+w, vm+h);
				}
			}
			dsplimit = w*h;
			return; // force exit
		}
		dsplimit = w*h*2;
		// else form paris
		// verify if we have an even or odd nr of DSPs
		bool singleDSP = false;
		if (nrDSPs % 2 == 1)
			singleDSP = true;
		
		int nrDSPsprim = nrDSPs;
			
		// compute new number of DSPs

		numberDSP4Overlap = dspCount = nrDSPs = (int) ceil((double) dspCount/2);

		
		// set shift amount according to target
		int shift = 0;
		if ((target_->getID() == "Virtex4") || (target_->getID() == "Virtex5"))
			shift = 17;
		
		int start = 0; // starting position
		
		
		//~ cout<< "mw "<<mw<<" mh "<<mh<<" w "<<w<<" h "<<h<<"estimate dsps "<<nrDSPsprim<<endl;
		//~ cout<<"h*0.9 + h*4="<<h*0.9 + h*4<<endl;
		//~ cout<<"cond "<<(mw == h*4 || mh ==h*4)<<" tout "<<((mw == h*4 || mh ==h*4) ||(  (h*0.9 + h*4 <=mw) &&  mh>=2*w  ) || (   (h*0.9 + h*4 <=mh) &&  mw>=2*w  ))<<endl;
		//~ cout<<" area= "<<mw*mh<<" area covered by DSPs= "<<h*w*nrDSPsprim<<" condition "<<( mw*mh*0.8>= h*w*nrDSPsprim )<<endl;
		
		
		
		
		if (dspCount >= 16 &&
			((mw == h*4 || mh ==h*4) ||
			(  (h*0.9 + h*4 <=mw) &&  mh>=2*w  ) ||
			(   (h*0.9 + h*4 <=mh) &&  mw>=2*w  ) ||
			( mw*mh>= h*w*nrDSPsprim )
			)) 
		{ // we have more than 16 paris of multipliers we can group 4 such pairs into a super-block
			cout<<"A super DSP was created"<<endl;
			config[start] = new DSP(0, w*2, h*4);
			config[start]->setNrOfPrimitiveDSPs(4);
			config[start]->allocatePositions(3*start); // each DSP offers 3 positions
			if(!replace(config, start))
			{
				int w=config[start]->getMaxMultiplierWidth();
				int h= config[start]->getMaxMultiplierHeight();
				config[start]->setTopRightCorner(vn, vm);
				config[start]->setBottomLeftCorner(vn+w, vm+h);
			}
			start++;
			dspCount -= 3;

			numberDSP4Overlap = nrDSPs = dspCount;

			
			int i;
			for (i=start; i<3; i++)
			{
				config[i] = new DSP(shift, w, h*2);	
				config[i]->setNrOfPrimitiveDSPs(2);
				config[i]->allocatePositions(3*i); // each DSP offers 3 positions
				if(!replace(config, i))
				{
					int w=config[i]->getMaxMultiplierWidth();
					int h= config[i]->getMaxMultiplierHeight();
					config[i]->setTopRightCorner(vn, vm);
					config[i]->setBottomLeftCorner(vn+w, vm+h*2);
				}
			}
			start = i;	
		}

		/*NOTE: if the program entered the previous if clause it will also enter the following 
		 * 	if clause. This is the effect we want to obtain */
		if (dspCount >= 10 		&&
			((mw == h*4 || mh ==h*4) ||
			(  (h*0.9 + h*4 <=mw) &&  mh>=2*w  ) ||
			(   (h*0.9 + h*4 <=mh) &&  mw>=2*w  ) ||
			( mw*mh>= h*w*nrDSPsprim )
			)	
		)
		{ // we have more than 10 paris of multipliers we can group 4 such pairs into a super-block
			cout<<"A super DSP was created"<<endl;
			config[start] = new DSP(0, w*2, h*4);
			config[start]->setNrOfPrimitiveDSPs(4);
			config[start]->allocatePositions(3*start); // each DSP offers 3 positions
			dspCount -= 3;
			numberDSP4Overlap = nrDSPs = dspCount;
			if(!replace(config, start))
			{
				int w=config[start]->getMaxMultiplierWidth();
				int h= config[start]->getMaxMultiplierHeight();
				config[start]->setTopRightCorner(vn, vm);
				config[start]->setBottomLeftCorner(vn+w, vm+h);
			}
			start++;

		}
		
		// initialize all DSP paris except first one
		dspCount--;
		for (int i=start; i<dspCount; i++)
		{
			config[i] = new DSP(shift, w, h*2);	
			config[i]->setNrOfPrimitiveDSPs(2);
			config[i]->allocatePositions(3*i); // each DSP offers 3 positions
			if(!replace(config, i))
			{
				int w=config[i]->getMaxMultiplierWidth();
				int h= config[i]->getMaxMultiplierHeight();
				config[i]->setTopRightCorner(vn, vm);
				config[i]->setBottomLeftCorner(vn+w, vm+h*2);
			}
		}
		
		//initializa last DSP or DSP pair
		if (singleDSP) // then the last position is a single DSP 
		{ // allocate single DSP
			config[dspCount] = target_->createDSP();
			config[dspCount]->setNrOfPrimitiveDSPs(1);
			config[dspCount]->allocatePositions(3*dspCount); // each DSP offers 3 positions
			if(!replace(config, dspCount))
			{
				int w=config[dspCount]->getMaxMultiplierWidth();
				int h= config[dspCount]->getMaxMultiplierHeight();
				config[dspCount]->setTopRightCorner(vn, vm);
				config[dspCount]->setBottomLeftCorner(vn+w, vm+h);
			}
			
		}
		else // then the last position is a DSP pair
		{ // allocate DSP pair
			config[dspCount] = new DSP(shift, w, h*2);
			config[dspCount]->setNrOfPrimitiveDSPs(2);
			config[dspCount]->allocatePositions(3*dspCount); // each DSP offers 3 positions
			if(!replace(config, dspCount))
			{
				int w=config[dspCount]->getMaxMultiplierWidth();
				int h= config[dspCount]->getMaxMultiplierHeight();
				config[dspCount]->setTopRightCorner(vn, vm);
				config[dspCount]->setBottomLeftCorner(vn+w, vm+h*2);
			}
			
		}	
		
	}
	
	int IntTruncMultiplier::bindDSPs4Stratix(DSP** config)
	{
		int nrDSPs_ = nrDSPs;
		int DSPcount = nrDSPs_;
	
		for (int i=0; i<nrDSPs_; i++)
			if (config[i] == NULL)
				DSPcount--;
		
			
		for (int i=0; i<DSPcount; i++)
			if (config[i]->getNumberOfAdders() < 3) // serach for an aligned DSP block that can be added with this one
				{
					int xtri, ytri, wx, wy, nrOp;
					bool bound;
					config[i]->getTopRightCorner(xtri, ytri);
					target_->getDSPWidths(wx,wy);
				
					if (verbose)
						cout << "bindDSP4Stratix: DSP #" << i << " has less than 3 operands. Top-right is at ( " << xtri << ", " << ytri << ") width is " << wx << endl;
					DSP** operands;
					DSP** operandsj;
					DSP** opsk;
					for (int j=0; j<DSPcount; j++)
						if (i != j)
							{
								nrOp = config[i]->getNumberOfAdders();
								if (nrOp == 3)
									break;
			
								operands = config[i]->getAdditionOperands();
								bound = false;
								// check if the DSP blocks are bound
								for (int k=0; k<nrOp; k++)
									if (operands[k] == config[j]) // the DSPs are allready bound
										{
											if (verbose)
												cout << "bindDSP4Stratix: DSP #" << j << " has #" << i << " as one of its operands allready." << endl;
											bound = true;
											break;
										}
				
								if (bound)
									continue;
								// if they are not bound
								int xtrj, ytrj, nrOpj;
								config[j]->getTopRightCorner(xtrj, ytrj);
								nrOpj = config[j]->getNumberOfAdders();
				
								if (verbose)
									cout << "bindDSP4Stratix:" << tab << " Checking against DSP #" << j << " Top-right is at ( " << xtrj << ", " << ytrj << ") width is " << wx << endl;
			
								if ((((xtrj-xtri == -wx) && (ytrj-ytri == wy)) || ((xtrj-xtri == wx) && (ytrj-ytri == -wy))) && // if they have the same alignment
									 (nrOpj + nrOp < 3)) // if the DSPs we want to bind have less than 3 other DSPs to add together
									{ // copy the operands from one another and bind them 
					
										if (verbose)
											cout << "bindDSP4Stratix : DSP #" << j << " together with #" << i << " have fewer than 3 operands. We can bind them." << endl;
										operandsj = config[j]->getAdditionOperands();
					
										for (int k=0; k<nrOp; k++)
											{
												operandsj[nrOpj+k] = operands[k];
												// each operand of congif[i] also gets bounded 
												int opcntk = operands[k]->getNumberOfAdders();
												operands[k]->setNumberOfAdders(nrOp+nrOpj+1);
												opsk = operands[k]->getAdditionOperands();
												for (int l=0; l < nrOpj; l++)
													opsk[l+opcntk] = operandsj[l];
												opsk[nrOpj+opcntk] = config[j];
												operands[k]->setAdditionOperands(opsk);
											}
										operandsj[nrOp+nrOpj] = config[i];
										config[j]->setAdditionOperands(operandsj);
										config[j]->setNumberOfAdders(nrOp+nrOpj+1);
					
										for (int k=0; k<nrOpj; k++)
											{
												operands[nrOp+k] = operandsj[k];
												// each operand of congif[j] also gets bounded 
												int opcntk = operandsj[k]->getNumberOfAdders();
												operandsj[k]->setNumberOfAdders(nrOp+nrOpj+1);
												opsk = operandsj[k]->getAdditionOperands();
												for (int l=0; l < nrOp; l++)
													opsk[l+opcntk] = operands[l];
												opsk[nrOp+opcntk] = config[i];
												operandsj[k]->setAdditionOperands(opsk);
											}
										operands[nrOp+nrOpj] = config[j];
										config[i]->setAdditionOperands(operands);
										config[i]->setNumberOfAdders(nrOp+nrOpj+1);
									}
							}
				}
		
		/* We now have:
		 * 	- pairs of DSP objects which have the number of operands set to 1
		 * 	- triplets of DSP objects which have the number of operands set to 2
		 * 	- quadruplets of DSP objects which have the number of operands set to 3 
		 * We keep a counter for each possible number of operands. */
		int pair[3] = {0, 0, 0};
	
	
		for (int i=0; i<DSPcount; i++)
			{
				if (verbose)
					cout << "bindDSP4Stratix : DSP #" << i << " has " << config[i]->getNumberOfAdders() << " adders" << endl;
				pair[config[i]->getNumberOfAdders()-1]++;
			}
		if (verbose)
			{
				cout << "bindDSP4Stratix : one " << pair[0] << endl;
				cout << "bindDSP4Stratix : two " << pair[1] << endl;
				cout << "bindDSP4Stratix : three " << pair[2] << endl;
			}
		return (DSPcount - pair[0]/2 - pair[1]*2/3 - pair[2]*3/4);
	}

	DSP** IntTruncMultiplier::splitLargeBlocks(DSP** config, int &numberOfDSPs)
	{
		int h, w, dspH, dspW, tmp, nrDSPonHeight, nrDSPonWidth, shiftAmount, newNrDSPs=0;
		getTarget()->getDSPWidths(dspW, dspH);
		
		// count total number of DSPs
		for (int i=0; i<numberOfDSPs; i++)
		{
			h = config[i]->getMaxMultiplierHeight();
			w = config[i]->getMaxMultiplierWidth();
			
			if ((h % dspH) != 0) // match width and height
			{
				tmp = dspH;
				dspH = dspW;
				dspW = tmp;
			}
			
			nrDSPonHeight = h/dspH;
			nrDSPonWidth = w/dspW;
			
			newNrDSPs += (nrDSPonHeight*nrDSPonWidth);
		}
		
		DSP** returnConfig = new DSP*[newNrDSPs];
		int index = 0;
		int xtr, xbl, ytr, ybl;
		for (int i=0; i<numberOfDSPs; i++)
		{
			h = config[i]->getMaxMultiplierHeight();
			w = config[i]->getMaxMultiplierWidth();
			shiftAmount = config[i]->getShiftAmount();
			config[i]->getTopRightCorner(xtr, ytr);
			config[i]->getBottomLeftCorner(xbl, ybl);
			
			if ((h % dspH) != 0) // match width and height
			{
				tmp = dspH;
				dspH = dspW;
				dspW = tmp;
			}
			
			nrDSPonHeight = h/dspH;
			nrDSPonWidth = w/dspW;
			// create DSP blocks
			for (int j=0; j<nrDSPonHeight; j++)
				for (int k=0; k<nrDSPonWidth; k++)
				{
					returnConfig[index] = getTarget()->createDSP();
					returnConfig[index]->setTopRightCorner(xtr+k*dspW, ytr+j*dspH);
					returnConfig[index]->setBottomLeftCorner(xtr+(k+1)*dspW-1, ytr+(j+1)*dspH-1);
					returnConfig[index]->setMultiplierHeight(dspH);
					returnConfig[index]->setMultiplierWidth(dspW);
					// take care of shiftings between DSPs
					if (shiftAmount == dspH)
					{
						if (j > 0) 
						{
							returnConfig[index]->setShiftIn(returnConfig[index-(j*nrDSPonWidth)]);
							returnConfig[index-(j*nrDSPonWidth)]->setShiftOut(returnConfig[index]);	
						}
					}
					else if (shiftAmount == dspW)
					{
						if (k > 0)
						{
							returnConfig[index]->setShiftIn(returnConfig[index-1]);
							returnConfig[index-1]->setShiftOut(returnConfig[index]);
						}
					}
					index++;
				}
		}
		numberOfDSPs = newNrDSPs;
		return returnConfig;
	}
	
	void IntTruncMultiplier::convertCoordinates(int &tx, int &ty, int &bx, int &by)
	{
		int ax = bx;
		int ay = by;
		bx = vnme - tx-1;
		by = vmme - ty-1;
		tx = vnme - ax-1;
		ty = vmme - ay-1;	
	}

	int IntTruncMultiplier::multiplicationInDSPs(DSP** config)
	{
		int nrOp = 0;		 			// number of resulting adder operands
		int trx1, try1, blx1, bly1; 	// coordinates of the two corners of a DSP block 
		int fpadX, fpadY, bpadX, bpadY;	// zero padding on both axis
		int multW, multH; 				// width and height of the multiplier the DSP block is using
		ostringstream xname, yname, mname, cname, sname;
		DSP** tempc = new DSP*[nrDSPs];	// buffer that with hold a copy of the global configuration
			
		//memcpy(tempc, config, sizeof(DSP*) * nrDSPs );
		for (int i=0; i<nrDSPs; i++)
		{
			tempc[i] = config[nrDSPs-i-1];
			DSP* si = tempc[i]->getShiftIn();
			DSP* so = tempc[i]->getShiftOut();
			tempc[i]->setShiftIn(so);
			tempc[i]->setShiftOut(si);
		}
		
		if ( ( target_->getID() == "Virtex4") ||
			 ( target_->getID() == "Virtex5") ||
			 ( target_->getID() == "Spartan3"))  // then the target is A Xilinx FPGA 
			{
				for (int i=0; i<nrDSPs; i++)
					if (tempc[i] != NULL)
						{
							//~ cout << "At DSP#"<< i+1 << " tempc["<<i<<"]" << endl; 
							DSP* d = tempc[i];
							int j=0;
							int connected = 0;
			
							while (d != NULL)
								{
									connected++;
									d = d->getShiftOut();
								}
			
							d = tempc[i];
			
							while (d != NULL)
								{
									d->getTopRightCorner(trx1, try1);
									d->getBottomLeftCorner(blx1, bly1);
									convertCoordinates(trx1, try1, blx1, bly1);
									int sh = 0;
									
									if (blx1 < try1)
										sh = 1;
										
									fpadX = wInX-blx1-1;
									//~ cout << "fpadX = " << fpadX << endl;
									fpadX = (fpadX<0)?0:fpadX;
									
									fpadY = wInY-bly1-1;
									//~ cout << "fpadY = " << fpadY << endl;
									fpadY = (fpadY<0)?0:fpadY;
									
									bpadX = trx1;
									bpadX = (bpadX<0)?0:bpadX;
									
									bpadY = try1;
									bpadY = (bpadY<0)?0:bpadY;
									
									multW = blx1-trx1+1;
									multH = bly1-try1+1;
			
									setCycle(0);
									
									xname.str("");
									xname << "x" << i << "_" << j;
									vhdl << tab << declare(xname.str(), multW) << " <= X" << range(blx1, trx1) << ";" << endl;
									
									yname.str("");
									yname << "y" << i << "_" << j;
									vhdl << tab << declare(yname.str(), multH) << " <= Y" << range(bly1, try1) << ";" << endl;
				
									if ((d->getShiftIn() != NULL) && (j>0)) // multiply accumulate
										{
											mname.str("");
											mname << "pxy" << i;
											cname.str("");
											cname << "txy" << i << j;
											setCycle(j);
											vhdl << tab << declare(cname.str(), multW+multH) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
											vhdl << tab << declare(join(mname.str(),j), multW+multH+1) << " <= (\"0\" & " << use(cname.str()) << ") + " << use(join(mname.str(), j-1)) << range(multW+multH-1, d->getShiftAmount()) << ";" << endl;	
											if (d->getShiftOut() == NULL) // concatenate the entire partial product
												{
													setCycle(connected);
													sname.seekp(ios_base::beg);
													//sname << zg(wInX+wInY+extW+extH-blx1-bly1-3, 0) << " & " << use(join(mname.str(),j)) << range(multW-fpadX + multH-fpadY-1, 0) << " & " << sname.str();
													sname << zg(fpadX+fpadY-1-sh, 0) << " & " << use(join(mname.str(),j)) << range(multW+multH, 0) << " & " << sname.str();
												}
											else // concatenate only the lower portion of the partial product
												{
													setCycle(connected);
													sname.seekp(ios_base::beg);
													sname << use(join(mname.str(),j)) << range(d->getShiftAmount()-1, 0) << " & " << sname.str();
												}
										}
									else // only multiplication
										{
											mname.str("");
											mname << "pxy" << i << j;
											vhdl << tab << declare(mname.str(), multW+multH) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
											sname.str("");
											if (d->getShiftOut() == NULL) // concatenate the entire partial product
												{
													setCycle(connected);
													//sname << zg(wInX+wInY+extW+extH-blx1-bly1-2, 0) << " & " << use(mname.str()) << range(multH+multW-1, bpadX+bpadY)<< " & " << zg(trx1-extW,0) << " & " << zg(try1-extH,0) <<  ";" << endl;
													sname << zg(fpadX+fpadY, 0) << " & " << use(mname.str()) << " & " << zg(trx1,0) << " & " << zg(try1+sh,0) <<  ";" << endl;
												}
											else // concatenate only the lower portion of the partial product
												{
													setCycle(connected);
													sname << use(mname.str()) << range(d->getShiftAmount()-1, 0) << " & " << zg(trx1,0) << " & " << zg(try1+sh,0) << ";" << endl;
												}
										}
				
									// erase d from the tempc buffer to avoid handleing it twice
									for (int k=i+1; k<nrDSPs; k++)
										{
											if ((tempc[k] != NULL) && (tempc[k] == d))
												{
													//~ cout << "tempc[" << k << "] deleted" << endl;
													tempc[k] = NULL;
													break;
												}
										}
				
				
									d = d->getShiftOut();
									j++;
								}	
							sname.seekp(ios_base::beg);
							sname << tab << declare(join("addOpDSP", nrOp),wInX+wInY) << " <= " << sname.str();
							vhdl << sname.str();
							nrOp++;		
						}
		
				return nrOp;
			}
		else // the target is Stratix
			{
				int boundDSPs;  				// number of bound DSPs in a group
				DSP** addOps;					// addition operands bound to a certain DSP
				int mPadX, mPadY, minPad;		// minimum padding of addition operands
	
				for (int i=0; i<nrDSPs; i++)
					if (tempc[i] != NULL)
						{
							mPadX = INT_MIN;
							mPadY = INT_MIN;
							cout << "At DSP#"<< i+1 << " tempc["<<i<<"]" << endl; 
							tempc[i]->getTopRightCorner(trx1, try1);
							tempc[i]->getBottomLeftCorner(blx1, bly1);
							convertCoordinates(trx1, try1, blx1, bly1);
							fpadX = blx1-wInX-getExtraWidth()+1;
							fpadX = (fpadX<0)?0:fpadX;
							fpadY = bly1-wInY-getExtraHeight()+1;
							fpadY = (fpadY<0)?0:fpadY;
							bpadX = getExtraWidth()-trx1;
							bpadX = (bpadX<0)?0:bpadX;
							mPadX = (bpadX>mPadX)?bpadX:mPadX;
							bpadY = getExtraHeight()-try1;
							bpadY = (bpadY<0)?0:bpadY;
							mPadY = (bpadY>mPadY)?bpadY:mPadY;
							minPad = bpadY+bpadX;
							multW = tempc[i]->getMaxMultiplierWidth();
							multH = tempc[i]->getMaxMultiplierHeight();
			
							setCycle(0);
							xname.str("");
							xname << "x" << i << "_0";
							vhdl << tab << declare(xname.str(), multW, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadX,0) << " & " << "X" << range(blx1-fpadX, trx1+bpadX) << " & " << zg(bpadX,0) << ";" << endl;
							yname.str("");
							yname << "y" << i << "_0";
							vhdl << tab << declare(yname.str(), multH, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadY,0) << " & " << "Y" << range(bly1-fpadY, try1+bpadY) << " & " << zg(bpadY,0) << ";" << endl;
			
							boundDSPs = tempc[i]->getNumberOfAdders();
							int ext = 0; // the number of carry bits of the addtion
							if (boundDSPs > 0) // need to traverse the addition operands list and perform addtion
								{
									ext = (boundDSPs>1)?2:1;
									cout << "boundDSPs = " << boundDSPs << endl;
									nextCycle();
									mname.str("");
									mname << "mult_" << i << "_0";
									vhdl << tab << declare(mname.str(), multW+multH, true, Signal::registeredWithAsyncReset) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
				
									addOps = tempc[i]->getAdditionOperands();
				
									//
									for (int j=0; j<3; j++)
										if (addOps[j] == NULL)
											cout << "addOps["<< j << "]=NULL" << endl;
										else
											cout << "addOps["<< j << "]=not null" << endl;
				
				
					
									//
				
				
									for (int j=0; j<boundDSPs; j++)
										{
											cout << "j = " << j << endl;
											// erase addOps[j] from the tempc buffer to avoid handleing it twice
											for (int k=i+1; k<nrDSPs; k++)
												{
													if ((tempc[k] != NULL) && (tempc[k] == addOps[j]))
														{
															cout << "tempc[" << k << "] deleted" << endl;
															tempc[k] = NULL;
															break;
														}
												}
					
											addOps[j]->getTopRightCorner(trx1, try1);
											addOps[j]->getBottomLeftCorner(blx1, bly1);
											fpadX = blx1-wInX-getExtraWidth()+1;
											fpadX = (fpadX<0)?0:fpadX;
											fpadY = bly1-wInY-getExtraHeight()+1;
											fpadY = (fpadY<0)?0:fpadY;
											bpadX = getExtraWidth()-trx1;
											bpadX = (bpadX<0)?0:bpadX;
											mPadX = (bpadX>mPadX)?bpadX:mPadX;
											bpadY = getExtraHeight()-try1;
											bpadY = (bpadY<0)?0:bpadY;
											mPadY = (bpadY>mPadY)?bpadY:mPadY;
											minPad = (minPad<(bpadY+bpadX))?minPad:(bpadY+bpadX);
											multW = addOps[j]->getMaxMultiplierWidth();
											multH = addOps[j]->getMaxMultiplierHeight();
					
											setCycle(0);
											xname.str("");
											xname << "x" << i << "_" << j+1;
											vhdl << tab << declare(xname.str(), multW, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadX,0) << " & " << "X" << range(blx1-fpadX, trx1+bpadX) << " & " << zg(bpadX,0) << ";" << endl;
											yname.str("");
											yname << "y" << i << "_" << j+1;
											vhdl << tab << declare(yname.str(), multH, true, Signal::registeredWithAsyncReset) << " <= " << zg(fpadY,0) << " & " << "Y" << range(bly1-fpadY, try1+bpadY) << " & " << zg(bpadY,0) << ";" << endl;
					
											nextCycle();
											mname.str("");
											mname << "mult_" << i << "_" << j+1;
											vhdl << tab << declare(mname.str(), multW+multH, true, Signal::registeredWithAsyncReset) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
										}
				
									nextCycle();
									vhdl << tab << declare(join("addDSP", nrOp), multW+multH+ext, true, Signal::registeredWithAsyncReset) << " <= ";
				
									int j=0;
				
									for (j=0; j<boundDSPs; j++)
										{
											mname.str("");
											mname << "mult_" << i << "_" << j;
											vhdl << "(" << zg(ext,0) << " & " << use(mname.str()) << ") + "; 
										}
				
									mname.str("");
									mname << "mult_" << i << "_" << j;
									vhdl << "(" << zg(ext,0) << " & " << use(mname.str()) << ");" << endl; 
				
								} 
							else // multiply the two terms and you're done
								{
									nextCycle();
									vhdl << tab << declare(join("addDSP", nrOp), multW+multH, true, Signal::registeredWithAsyncReset) << " <= " << use(xname.str()) << " * " << use(yname.str()) << ";" << endl;
								}
							vhdl << tab << declare(join("addOpDSP", nrOp), wInX+wInY) << " <= " << zg(wInX+minPad-blx1-ext,0) << " & " << zg(wInY-bly1,0) << " & " << use(join("addDSP", nrOp)) << range(multW+multH+ext-1, minPad) << " & " << zg(trx1+try1+minPad-mPadX-mPadY,0) << ";" << endl;
							nrOp++;
						}
				return nrOp;
			}

	}

	int IntTruncMultiplier::multiplicationInSlices(vector<SoftDSP*> config)
	{
		unsigned partitions;
		for (partitions=0; partitions<config.size(); partitions++)
		{	
			int njj, nj, ni, nii;
			config[partitions]->getCoordinates(nj, ni, njj, nii);
			convertCoordinates(nj, ni, njj, nii);
			int sh = 0;
			if (ni < njj)
				sh = 1;
			// CODE GENERATING VHDL 
			setCycle(0);	
			target_->setUseHardMultipliers(false);
			IntMultiplier* mult =  new IntMultiplier(target_, njj-nj+1, nii-ni+1);
			ostringstream cname;
			cname << mult->getName() << "_" << partitions;
			mult->changeName(cname.str());
			oplist.push_back(mult);
			// TODO: compute width of x and y + corretc range for X and Y
			vhdl << tab << declare(join("x_",partitions), njj-nj+1) << " <= " << "X" << range(njj, nj) << ";" << endl;
			inPortMap(mult, "X", join("x_",partitions));
			vhdl << tab << declare(join("y_",partitions), nii-ni+1) << " <= " << "Y" << range(nii, ni) << ";" << endl;
			inPortMap(mult, "Y", join("y_",partitions));
		
			outPortMap(mult, "R", join("result", partitions));
		
			vhdl << instance(mult, join("Mult", partitions));
		
			syncCycleFromSignal(join("result", partitions));
			
			vhdl << tab << declare(join("addOpSlice", partitions), wInX+wInY) << " <= " << zg(wInX-njj-1+wInY-nii-1-sh, 0) << " & " << join("result", partitions) << " & " << zg(nj+ni+sh, 0) << ";" << endl;
			cout<<"partitions " << partitions << " @ cycle " << getCurrentCycle() << endl;
		}
				
		return (int)partitions;
	}	

	void IntTruncMultiplier::generateVHDLCode(DSP** config, vector<SoftDSP*> softConfig){
		
		DSP** bestConfig = config;
		bindDSPs(bestConfig);      
		bestConfig = splitLargeBlocks(bestConfig, nrDSPs);
		int nrDSPOperands = multiplicationInDSPs(bestConfig);
		int nrSliceOperands = multiplicationInSlices(softConfig);
		map<string, double> inMap;
		
		for (int j=0; j<nrDSPOperands; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOpDSP" << j;
				syncCycleFromSignal(concatPartialProd.str());
			}	
			
		for (int j=0; j<nrSliceOperands; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOpSlice" << j;
				syncCycleFromSignal(concatPartialProd.str());
			}	
		
		IntNAdder* add =  new IntNAdder(getTarget(), wInX+wInY, nrDSPOperands+nrSliceOperands, inMap);
		//IntCompressorTree* add =  new IntCompressorTree(target, adderWidth, opCount);
		oplist.push_back(add);

		for (int j=0; j<nrDSPOperands; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOpDSP" << j;

				inPortMap (add, join("X",j) , concatPartialProd.str());
			}	
			
		for (int j=0; j<nrSliceOperands; j++)
			{
				ostringstream concatPartialProd;
				concatPartialProd  << "addOpSlice" << j;
				cout << "@ In Port Map Current Cycle is " << getCurrentCycle() << endl;
				inPortMap (add, join("X",j+nrDSPOperands) , concatPartialProd.str());
			}	
	
		inPortMapCst(add, "Cin", "'0'");
		outPortMap(add, "R", "addRes");
		vhdl << instance(add, "adder");

		syncCycleFromSignal("addRes");

		vhdl << tab << "R <= " << "addRes" << range(wX+wY-1, targetPrecision) << ";" << endl;
	}
	
	vector<SoftDSP*> IntTruncMultiplier::insertSoftDSPs(DSP** config)
	{
		
		//int costSlice=0;
	
		int n,m;
		int count=1;
		n=vn;
		m=vm;
	
		//int nmew = vnme;
		int ew = getExtraWidth();
		//int mmeh = vmme;
		int eh = getExtraHeight();
		//int nj,ni,njj,nii;
		vector<SoftDSP*> configSoft;
		configSoft.clear();
		SoftDSP *tempSoft;
		//configSoft.push_back(temp);
		//configSoft[1];
		//configSoft.size();
		
		
		
		//~ cout<<"width "<<n<<"height "<<m<<endl;
	
		for(int i=0;i<m;i++)
			{
			
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			}
		
		if(isSquarer)
		{
			for (int j=0; j<m; j++)
				{
					for (int i=0; i<n-j-1; i++)
						mat[j][i] = 1369;
				}
				
			for (int i=0; i<ew+truncationOffset; i++)
				for (int j=m-truncationOffset+i-eh; j<m; j++)
					mat[j][i] = 1369;
		}			
			
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					printf("%5d",mat[i][j]);
				printf("\n");
			}	
			
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				//~ cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
		
		
		int deepX=vn,deepY=vm,rdeepX=vm;
			int ti;
		int sbx,sby;
		//~ for(int i=40;i<vmme;i++)
			//~ mat[i][vnme-1]=nrDSPs+1;
		//mat[vmme-1][vnme-2]=4;	
		
		//~ for(int i=35;i<vnme - 24;i++)
			//~ mat[eh][i]=nrDSPs+1;
		//~ mat[eh+1][35]=4;
		
		//mat[40][vnme-1]=nrDSPs+1;
		int lastModifiedIndex = -1;
		int lastDeepX, lastDeepY;	
			
		while(!isTilingValid(config,configSoft,targetPrecision) )
		{
		deepX=vn,deepY=vm,rdeepX=vm;
		bool found = false;
		//search the deapest position
		for(int i=ew+1;i<=vnme;i++)	
		{
			ti = vn-i;
			for(int j=eh;j<vmme;j++)
			{
				if(mat[j][ti]==0&&((deepX+deepY)>(i+j))  )
				{
					found = true;
					deepX=i;
					rdeepX=ti;
					deepY=j;					
				}				
			}			
		}	
		
		
		if (!found)
			break;
			
		//cout<<"Deepest position is X="<<deepX<<" Y="<<deepY<<" realX"<<rdeepX<<endl;
		mat[deepY][rdeepX]=count++;
		
		//small display
		//~ for(int i=0;i<vm;i++)	
		//~ {
			//~ for(int j=0;j<vn;j++)
			//~ {
				//~ cout<<mat[i][j];
			//~ }
			//~ cout<<endl;
		//~ }
		
				
		//puts the current 1x1 multiplication on the board and will try to combine and extend the neighbours
		//bool isextentionok=false;
		int ref;
		//try in right
		if(  mat[deepY][rdeepX+1]!=1369 && mat[deepY][rdeepX+1]>nrDSPs && mat[deepY-1][rdeepX+1]!=mat[deepY][rdeepX+1])
		{
			//cout<<"Has neighbor in left of it"<<endl;
			ref = mat[deepY][rdeepX+1];
			int i;
			for(i=deepY+1; mat[i][rdeepX]==0 && mat[i][rdeepX+1]==ref;i++ )
			{
				mat[i][rdeepX]=ref;
			}
			if(  mat[i][rdeepX+1]==ref  )
			{
				//new multiplier should be created
				for(i=deepY+1;mat[i][rdeepX]==ref;i++)
					mat[i][rdeepX]=0;
				//cout<<"New SoftDSP was created"<<endl;
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 0;
				lastDeepY = 0;
				tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
				configSoft.push_back(tempSoft);
			}
			else
			{
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 1;
				lastDeepY = 0;
				//the previous multiplier was extended
				mat[deepY][rdeepX]=ref;
				count--;		
				configSoft[(ref-nrDSPs-1)]->getBottomLeftCorner(sbx,sby);
				sbx++;
				configSoft[(ref-nrDSPs-1)]->setBottomLeftCorner(sbx,sby);
				
						
			}			
		}
		else
		{
			//try top
			if(  mat[deepY-1][rdeepX]!=1369 && mat[deepY-1][rdeepX]>nrDSPs &&  mat[deepY-1][rdeepX+1]!=mat[deepY-1][rdeepX])
			{
				//cout<<"Has neighbor in top of it"<<endl;
				ref = mat[deepY-1] [rdeepX];
				int i;
				for(i=rdeepX-1;mat[deepY][i]==0 && mat[deepY-1][i]==ref;i--)
				{
					mat[deepY][i]=ref;
				}
				if(mat[deepY-1][i]==ref)
				{
					//new multiplier should be created
					for(i=rdeepX-1;mat[deepY][i]==ref;i--)
						mat[deepY][i]=0;
					//cout<<"New SoftDSP was created"<<endl;
					lastModifiedIndex = ref-nrDSPs-1;
					lastDeepX = 0;
					lastDeepY = 0;
					tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
					configSoft.push_back(tempSoft);
				}
				else
				{
					lastModifiedIndex = ref-nrDSPs-1;
					lastDeepX = 0;
					lastDeepY = 1;
					//the previous multiplier was extended
					mat[deepY][rdeepX]=ref;
					count--;
					configSoft[(ref-nrDSPs-1)]->getBottomLeftCorner(sbx,sby);
					sby++;
					configSoft[(ref-nrDSPs-1)]->setBottomLeftCorner(sbx,sby);
				}				
				
			}
			else
			{
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 0;
				lastDeepY = 0;
				//no SoftDSP exists in the neighbourhood
				//cout<<"New SoftDSP was created"<<endl;
				tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
				configSoft.push_back(tempSoft);
				
			}
		}
		
		
		//displayAll(config,configSoft);
		
		//small display
		//~ for(int i=0;i<vm;i++)	
		//~ {
			//~ for(int j=0;j<vn;j++)
			//~ {
				//~ cout<<mat[i][j];
			//~ }
			//~ cout<<endl;
		//~ }
		
		}
		
		// cut last expanded multiplier and fit to reduce resource use
		if  ((mpfr_cmp(errorSum, targetError) < 0) && (lastDeepX != lastDeepY))
		{
			configSoft[lastModifiedIndex]->getBottomLeftCorner(sbx, sby);
			configSoft[lastModifiedIndex]->setBottomLeftCorner(sbx-lastDeepX, sby-lastDeepY);
			
			int stx, sty;
			configSoft[lastModifiedIndex]->getTopRightCorner(stx, sty);
			sbx = stx = sbx*lastDeepX + stx*((lastDeepX==1)?0:1);
			sby = sty = sby*lastDeepY + sty*((lastDeepY==1)?0:1);
			
			lastDeepX = ((lastDeepX==1)?0:1);
			lastDeepY = ((lastDeepY==1)?0:1);
			tempSoft =  new SoftDSP(stx, sty, sbx, sby);
			configSoft.push_back(tempSoft);
			int idx = configSoft.size()-1;
			
			while (!isTilingValid(config,configSoft,targetPrecision))
			{
				sbx += lastDeepX;
				sby += lastDeepY;
				configSoft[idx]->setBottomLeftCorner(sbx, sby);
			}	
		}
		
		vector<SoftDSP*> newSoft;
		newSoft.clear();
		
		if (isSquarer)
		{
			for (int i=m-1; i>=0; i--)
			{
				if ((mat[i][n-i-1] > nrDSPs) && (mat[i][n-i-1] != 1369) && (mat[i-1][n-i-1] == 1369) && (mat[i][n-i-2] == 1369))
				{
					//cout << "in " << i << " " << mat[i-1][n-i-1] << " " << mat[i-1][n-i-2] << mat[i][n-i-1] <<endl;
					mat[i][n-i-1] = -1;//int ref = mat[i][n-i-1];
				
					int j = n-i-1;
					// expand ref
					//cout << "ref=" << ref << endl;
					while ((mat[i][j+1] > nrDSPs) && (mat[i][j+1] < count) && (mat[n-2-j][j+1] > nrDSPs))
					{
						//cout << "under " << mat[i][j+1] << " over " << mat[n-1-j][j] << endl;
						j++;
						for (int k=n-i-1; k<=j; k++)
						{
							mat[n-1-j][k] = -1;//ref;
							mat[n-1-k][j] = -1;//ref;
						}
					}
					j -= n-i-1;
					tempSoft = new SoftDSP(i-j, i-j, i, i);
					newSoft.push_back(tempSoft);
				}
			}	
			bool modified;
			do{
			modified = false;
			for (int j=0; j<n; j++)
			{
				for (int i=m-j-1; i<(m-1); i++)
				{
					if ((mat[i][j] > nrDSPs) && (mat[i][j] != 1369) && (mat[i][j-1] == 1369) && (mat[i+1][j] == 1369) && (j!=n-i-1)) //&& (mat[i-1][j-1] == 1369))
					{
						mat[i][j] = -1;//int ref = mat[i][j];
						// expand ref
						//cout << "ref=" << ref << endl;	
						int oldi = i;
						int oldj = j;
						// go down
						while ((i<(m-1)) && (mat[i+1][j] > nrDSPs) && (mat[i+1][j] != 1369))
						{
							i++;
							mat[i][j] = -1;
						}
						// go right & down
						while ((i<(m-1)) && (j<(n-1)) && (mat[oldi][j+1] > nrDSPs) && (mat[oldi][j+1] < count) && (mat[i+1][j+1] == mat[oldi][j+1]))// && (mat[n-2-j][j+1] > nrDSPs)
						{
							j++;
							i++;
							for (int k=oldi; k<=i; k++)
								mat[k][j] = -1;//ref;
							for (int k=oldj; k<=j; k++)	
								mat[i][k] = -1;//ref;
							
						}
						// go right
						while ((j<(n-1)) && (mat[oldi][j+1] > nrDSPs) && (mat[i][j+1] == mat[oldi][j+1]))
						{
							j++;
							for (int k=oldi; k<=i; k++)
								mat[k][j] = -1;//ref;
						}
						
						tempSoft = new SoftDSP(m-j-1, oldi, n-oldj-1, i);
						newSoft.push_back(tempSoft);
						modified = true;
					}
					
					cout << "i="<<i<<" j="<<j<< endl;
					if ((mat[i][j] > nrDSPs) && (mat[i][j] != 1369))
					{
						int ref = mat[i][j];
						mat[i][j] = -1;
						// expand ref
						cout << "ref=" << ref <<" i="<<i<<" j="<<j<< endl;	
						int oldi = i;
						int oldj = j;
						// go down
						while ((i<m) && (mat[i+1][j] > nrDSPs))
						{
							i++;	
							mat[i][j] = -1;
						}
						
						// go right
						while ((j<(n-1)) && (mat[oldi][j+1] > nrDSPs) && (mat[i][j+1] == mat[oldi][j+1]))
						{
							j++;
							for (int k=oldi; k<=i; k++)
								mat[k][j] = -1;
						}
						
						tempSoft = new SoftDSP(m-j-1, oldi, n-oldj-1, i);
						newSoft.push_back(tempSoft);
						modified = true;
					}
					
				}
			}
			}while(modified); 
			configSoft = newSoft;
		}
		
		for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
					printf("%5d",mat[i][j]);
				printf("\n");
			}	
			
		/* 
		
		 
		if (isSquarer)
		{
			int tx1, ty1, bx1, by1, tx2, ty2, bx2, by2;
			for (int i=(configSoft.size()-1); i>=0; i--)
			{
				configSoft[i]->getTopRightCorner(tx2, ty2);
				configSoft[i]->getBottomLeftCorner(bx2, by2);
				
				if ((tx2==bx2) && (ty2==by2) && (tx2==ty2)) // single 1x1 on the secondary diagonal
				{
					tx1 = tx2;
					ty1 = ty2;
					bx1 = bx2;
					by1 = by2;
					tempSoft = new SoftDSP(tx2, ty2, bx2, by2);
					newSoft.push_back(tempSoft);
				}
				else if ((by1 == by2) && ((tx1-1)==tx2)) // expand previously found multiplier
				{
					tx1 = tx2;
					tempSoft->setTopRightCorner(tx2, ty2);	
				}
				else
				{
					tempSoft = new SoftDSP(tx2, ty2, bx2, by2);
					newSoft.push_back(tempSoft);
				}
			}
			
			configSoft = newSoft;
		}
		*/

		//~ int stx,sty,
		//~ cout<<"Number of SoftDSPs created is "<<configSoft.size()<<endl;
		//~ for(int i=0;i<configSoft.size();i++)
		//~ {
			//~ configSoft[i]->getTopRightCorner(stx,sty);
			//~ configSoft[i]->getBottomLeftCorner(sbx,sby);
			//~ cout<<"DSP# "<<i+nrDSPs+1<<"has top-right corner X="<<stx<<" and Y= "<<sty<<" and bottom-left corner X"<<sbx<<" and Y="<<sby<<endl;
		//~ }
		//~ cout<<targetPrecision<<endl;
	
		return configSoft;
	}
	
	vector<SoftDSP*>IntTruncMultiplier::insertSoftDSPswithLimits(DSP** config)
	{
		
		
		//int costSlice=0;
	
		int n,m;
		int count=1;
		n=vn;
		m=vm;
	
		//int nmew = vnme;
		int ew = getExtraWidth();
		//int mmeh = vmme;
		int eh = getExtraHeight();
		//int nj,ni,njj,nii;
		vector<SoftDSP*> configSoft;
		configSoft.clear();
		SoftDSP *tempSoft;
		//configSoft.push_back(temp);
		//configSoft[1];
		//configSoft.size();
		
		
		
		//~ cout<<"width "<<n<<"height "<<m<<endl;
		
		
		for(int i=0;i<m;i++)
			{
			
				for(int j=0;j<n;j++)
					mat[i][j]=0;
			
			}
		if(isSquarer)
		{
			for(int j=eh;j<vmme;j++)
				{
					for(int i=ew;i<vnme-j-1;i++)
						mat[j][i] = 1369;
				}
		}
		
			
			
		for(int i=0;i<nrDSPs;i++)
			{
				int c1X,c2X,c1Y,c2Y;
			
				config[i]->getTopRightCorner(c1X,c1Y);
				config[i]->getBottomLeftCorner(c2X,c2Y);
				//~ cout<<"DSP #"<<i+1<<"has toprigh ("<<c1X<<","<<c1Y<<") and botomleft ("<<c2X<<","<<c2Y<<")"<<endl;
				c1X=n-c1X-1;
				c2X=n-c2X-1;
				//~ cout<<"new x1 "<<c1X<<" new x2 "<<c2X<<endl;
			
				fillMatrix(mat,n,m,c2X,c1Y,c1X,c2Y,count);
				count++;			
			}
		
		
		int deepX=vn,deepY=vm,rdeepX=vm;
			int ti;
		int sbx,sby;
		int stx,sty;
		//~ for(int i=40;i<vmme;i++)
			//~ mat[i][vnme-1]=nrDSPs+1;
		//mat[vmme-1][vnme-2]=4;	
		
		//~ for(int i=35;i<vnme - 24;i++)
			//~ mat[eh][i]=nrDSPs+1;
		//~ mat[eh+1][35]=4;
		
		//mat[40][vnme-1]=nrDSPs+1;
		int lastModifiedIndex = -1;
		int lastDeepX, lastDeepY;	
			
			
			
			
		while(!isTilingValid(config,configSoft,targetPrecision) )
		{
		deepX=vn,deepY=vm,rdeepX=vm;
		bool found = false;
		//search the deapest position
		for(int i=ew+1;i<=vnme;i++)	
		{
			ti = vn-i;
			for(int j=eh;j<vmme;j++)
			{
				if(mat[j][ti]==0&&((deepX+deepY)>(i+j))  )
				{
					found = true;
					deepX=i;
					rdeepX=ti;
					deepY=j;					
				}				
			}			
		}	
		
		
		if (!found)
			return configSoft;
			
		//~ cout<<"Deepest position is X="<<deepX<<" Y="<<deepY<<" realX"<<rdeepX<<endl;
		mat[deepY][rdeepX]=count++;
		
		//small display
		//~ for(int i=0;i<vm;i++)	
		//~ {
			//~ for(int j=0;j<vn;j++)
			//~ {
				//~ cout<<mat[i][j];
			//~ }
			//~ cout<<endl;
		//~ }
		
				
		//puts the current 1x1 multiplication on the board and will try to combine and extend the neighbours
		//bool isextentionok=false;
		int ref;
		int sdspw,sdsph;
		//try in right
		if(   mat[deepY][rdeepX+1]!=1369 &&   mat[deepY][rdeepX+1]>nrDSPs && mat[deepY-1][rdeepX+1]!=mat[deepY] [rdeepX+1])
		{
			//~ cout<<"Has neighbor in left of it"<<endl;
			ref = mat[deepY] [rdeepX+1];
			int i;
			
			
			configSoft[(ref-nrDSPs-1)]->getCoordinates(stx,sty,sbx,sby);
			sdspw = sbx-stx+1;
			sdsph = sby-sty+1;
			
			
			//~ for(i=deepX+1;i<vn&&mat[deepY][i]==ref;i++)
				//~ countlimit++;
				
			if( ( ((float)dspWidth*dspHeight)*0.5)< (sdspw*sdsph) &&  sdspw>=dspWidth ){
				
				//~ cout<<"New SoftDSP was created because the previous one exeeds the limit size"<<endl;
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 0;
				lastDeepY = 0;
				tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
				configSoft.push_back(tempSoft);
				
			}
			else
			{
			
			for(i=deepY+1; mat[i][rdeepX]==0 && mat[i][rdeepX+1]==ref;i++ )
			{
				mat[i][rdeepX]=ref;
			}
			if(  mat[i][rdeepX+1]==ref  )
			{
				//new multiplier should be created
				for(i=deepY+1;mat[i][rdeepX]==ref;i++)
					mat[i][rdeepX]=0;
				//~ cout<<"New SoftDSP was created"<<endl;
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 0;
				lastDeepY = 0;
				tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
				configSoft.push_back(tempSoft);
			}
			else
			{
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 1;
				lastDeepY = 0;
				//the previous multiplier was extended
				mat[deepY][rdeepX]=ref;
				count--;		
				configSoft[(ref-nrDSPs-1)]->getBottomLeftCorner(sbx,sby);
				sbx++;
				configSoft[(ref-nrDSPs-1)]->setBottomLeftCorner(sbx,sby);
				
						
			}
			}			
		}
		else
		{
			//try top
			if(   mat[deepY-1][rdeepX]!=1369 && mat[deepY-1][rdeepX]>nrDSPs &&  mat[deepY-1][rdeepX+1]!=mat[deepY-1][rdeepX])
			{
				//~ cout<<"Has neighbor in top of it"<<endl;
				ref = mat[deepY-1] [rdeepX];
				int i;
				//~ int countlimit=0;
				//~ for(i=deepY-1;i>0&&mat[i][rdeepX]==ref;i--)
					//~ countlimit++;
				
				
			configSoft[(ref-nrDSPs-1)]->getCoordinates(stx,sty,sbx,sby);
			sdspw = sbx-stx+1;
			sdsph = sby-sty+1;
			
				
			if( ( ((float)dspWidth*dspHeight)*0.5)< (sdspw*sdsph) &&  sdsph>=dspHeight )
				{
					//~ cout<<"New SoftDSP was created because the previous one exeeds the limit size"<<endl;
					lastModifiedIndex = ref-nrDSPs-1;
					lastDeepX = 0;
					lastDeepY = 0;
					tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
					configSoft.push_back(tempSoft);
				}
				else
				{
				
				for(i=rdeepX-1;  mat[deepY][i]==0 && mat[deepY-1][i]==ref;i--)
				{
					mat[deepY][i]=ref;
				}
				if( mat[deepY-1][i]==ref)
				{
					//new multiplier should be created
					for(i=rdeepX-1;mat[deepY][i]==ref;i--)
						mat[deepY][i]=0;
					//~ cout<<"New SoftDSP was created"<<endl;
					lastModifiedIndex = ref-nrDSPs-1;
					lastDeepX = 0;
					lastDeepY = 0;
					tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
					configSoft.push_back(tempSoft);
				}
				else
				{
					lastModifiedIndex = ref-nrDSPs-1;
					lastDeepX = 0;
					lastDeepY = 1;
					//the previous multiplier was extended
					mat[deepY][rdeepX]=ref;
					count--;
					configSoft[(ref-nrDSPs-1)]->getBottomLeftCorner(sbx,sby);
					sby++;
					configSoft[(ref-nrDSPs-1)]->setBottomLeftCorner(sbx,sby);
				}	
				}				
				
			}
			else
			{
				lastModifiedIndex = ref-nrDSPs-1;
				lastDeepX = 0;
				lastDeepY = 0;
				//no SoftDSP exists in the neighbourhood
				//cout<<"New SoftDSP was created"<<endl;
				tempSoft =  new SoftDSP(deepX-1,deepY,deepX-1,deepY);
				configSoft.push_back(tempSoft);
				
			}
		}
		
		
		//displayAll(config,configSoft);
		
		//small display
		//~ for(int i=0;i<vm;i++)	
		//~ {
			//~ for(int j=0;j<vn;j++)
			//~ {
				//~ cout<<mat[i][j];
			//~ }
			//~ cout<<endl;
		//~ }
		
		}
		
		
		
		if  ((mpfr_cmp(errorSum, targetError) < 0) && (lastDeepX != lastDeepY))
		{
			configSoft[lastModifiedIndex]->getBottomLeftCorner(sbx, sby);
			configSoft[lastModifiedIndex]->setBottomLeftCorner(sbx-lastDeepX, sby-lastDeepY);
			
			int stx, sty;
			configSoft[lastModifiedIndex]->getTopRightCorner(stx, sty);
			sbx = stx = sbx*lastDeepX + stx*((lastDeepX==1)?0:1);
			sby = sty = sby*lastDeepY + sty*((lastDeepY==1)?0:1);
			
			lastDeepX = ((lastDeepX==1)?0:1);
			lastDeepY = ((lastDeepY==1)?0:1);
			tempSoft =  new SoftDSP(stx, sty, sbx, sby);
			configSoft.push_back(tempSoft);
			int idx = configSoft.size()-1;
			
			while (!isTilingValid(config,configSoft,targetPrecision))
			{
				sbx += lastDeepX;
				sby += lastDeepY;
				configSoft[idx]->setBottomLeftCorner(sbx, sby);
			}	
		}
		 
		//~ int stx,sty,
		//~ cout<<"Number of SoftDSPs created is "<<configSoft.size()<<endl;
		//~ for(int i=0;i<configSoft.size();i++)
		//~ {
			//~ configSoft[i]->getTopRightCorner(stx,sty);
			//~ configSoft[i]->getBottomLeftCorner(sbx,sby);
			//~ cout<<"DSP# "<<i+nrDSPs+1<<"has top-right corner X="<<stx<<" and Y= "<<sty<<" and bottom-left corner X"<<sbx<<" and Y="<<sby<<endl;
		//~ }
		//~ cout<<targetPrecision<<endl;
	
		
		
		//small display
		//~ for(int i=0;i<vm;i++)	
		//~ {
			//~ for(int j=0;j<vn;j++)
			//~ {
				//~ cout<<mat[i][j];
			//~ }
			//~ cout<<endl;
		//~ }
		
		
		 //~ cout<<"Number of SoftDSPs created is "<<configSoft.size()<<endl;
		//~  cout<<"DSP width="<<dspWidth<<" height="<<dspHeight<<endl;
		return configSoft;
		
	}
	
	
}
