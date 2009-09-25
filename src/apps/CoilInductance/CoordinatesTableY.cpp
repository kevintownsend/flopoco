#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "../../utils.hpp"
#include "CoordinatesTableY.hpp"
#include <stdio.h>
#include <mpfr.h>



using namespace std;
using std::ifstream;



namespace flopoco{

 
	CoordinatesTableY::CoordinatesTableY(Target* target, int wIn, int LSBI,int MSBI,char *filename) : 
		DualTable(target, wIn, (MSBI-LSBI)), MSBI(MSBI), LSBI(LSBI), adrWidth(wIn) , wOutm(MSBI-LSBI),filepath(filename)
	{
		if ((MSBI < LSBI)){
			cerr << 
				" CoordinatesTableY: Input constraint LSBI <= MSBI not met."<<endl;
			exit (EXIT_FAILURE);
		}
	   
		ostringstream name; 
		name <<"CoordinatesX_"<<wIn<<"_"<<abs(LSBI)<<"_"<<abs(MSBI);
		setName(name.str());
	
		readParams();
	
	}

	CoordinatesTableY::~CoordinatesTableY() {}
	

	int    CoordinatesTableY::double2input(double x){
		int result;
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
		return result;
	}


	double CoordinatesTableY::input2double(int x) {
		double y;
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
		return(y);
	}

	mpz_class CoordinatesTableY::double2output(double x){
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
		return 0;
	}

	double CoordinatesTableY::output2double(mpz_class x) {
		double y;
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
  
		return(y);
	}


	void CoordinatesTableY::readParams()
	{
		ifstream indata; // indata is like cin
		indata.open(filepath);
	
	
		indata>>L>>NrSVert;
	   
		config = new int*[L];
		for(int j=0;j<L;j++)
			config[j]= new int[NrSVert];
	 

		nrTurns=0;
		for(int j=0;j<L;j++)
			for(int i=0;i<NrSVert;i++)
				{indata>>config[j][i];
					nrTurns++;
				}
	
		indata>>out_radius>>turn_radius>>insulation>>N;

		indata.close();
		initParameters();	
	}
	void CoordinatesTableY::initParameters()
	{
		//declaratii si initializari
		mpfr_set_default_prec(1000);
	
	
		mpfr_inits(out_rad_m,turn_rad_m,insu_m,nr_m,(mpfr_ptr) 0);
		mpfr_set_d(out_rad_m,out_radius,GMP_RNDN);
		mpfr_set_d(turn_rad_m,turn_radius,GMP_RNDN);
		mpfr_set_d(insu_m,insulation,GMP_RNDN);
	
		mpfr_t temp1,temp2;
		mpfr_inits( temp1,temp2,(mpfr_ptr) 0);
	
		//float h= ((float)(NrSVert-1))*(2.0*turn_radius+insulation);
	
	
		mpfr_init2(h,1000);
	
		mpfr_mul_si(temp1,turn_rad_m,2,GMP_RNDN);
		mpfr_add(temp1,temp1,insu_m,GMP_RNDN);
		mpfr_mul_si(h,temp1,(NrSVert-1),GMP_RNDN);
	
		//float fi=2*M_PI/N;
	
	
	
		mpfr_init2(fi,1000);
		mpfr_const_pi(temp2,GMP_RNDN);
		mpfr_div_si(temp1,temp2,N,GMP_RNDN);
		mpfr_mul_ui(fi,temp1,2,GMP_RNDN);
		
	
	}

	mpz_class CoordinatesTableY::function(int x)
	{
	
		
		mpz_class r;
	
		mpfr_t temp1,temp2;
		mpfr_inits2(1000, temp1,temp2,(mpfr_ptr) 0);
	
		mpfr_t ddj,deltak,rk;
		mpfr_init2(ddj,1000);
		mpfr_init2(deltak,1000);
		mpfr_init2(rk,1000);
		mpfr_t temp,tauqj;
		mpfr_init2(temp,1000);
		mpfr_init2(tauqj,1000);


		mpfr_t yv;
		mpfr_init2(yv,1000);
	
		int q;
	
		//set the coefficients that traverse the coil system
		int myx=x;
		q=myx%N;
		myx/=N;
		int k=0;
		int j=0;
		int i=0;
		int ver=1;
	

		for(i=0;(i<L)&&(ver==1);i++)
			for(j=0;(j<NrSVert)&&(ver==1);j++)
				for(k=0;(k<config[i][j])&&(ver==1);k++)
					{
						if(myx>0)
							myx--;
						else
							ver=0;
				
					}
				
		if(myx==0){
					
			//~ //computing the mpfr value of this coordinate
				
	
			//ddj=((float)j)*(2.0*turn_radius+insulation);
			mpfr_mul_ui(temp1,turn_rad_m,2,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_m,GMP_RNDN);
			mpfr_mul_ui(ddj,temp1,j,GMP_RNDN);
			
				
			//float deltak=((float)k)*(2.0*turn_radius+insulation);
			mpfr_mul_ui(temp1,turn_rad_m,2,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_m,GMP_RNDN);
			mpfr_mul_ui(deltak,temp1,k,GMP_RNDN);
			//float rk=out_radius-deltak;
			mpfr_sub(rk,out_rad_m,deltak,GMP_RNDN);
			
				
										
			//float yv=rk*sin(q*fi);
					
			mpfr_mul_ui(temp1,fi,q,GMP_RNDN);
			mpfr_sin(temp1,temp1,GMP_RNDN);
			mpfr_mul(yv,rk,temp1,GMP_RNDN);
					
		
		
			ver =1;
		
			if(mpfr_sgn(yv)<0)
				{ver=0;
					mpfr_mul_si(yv,yv,(-1),GMP_RNDN);
				}
		
			mpfr_set_ui(temp1, 2 , GMP_RNDN);
			mpfr_set_si(temp2, abs(LSBI) , GMP_RNDN);
			mpfr_pow(temp1, temp1, temp2, GMP_RNDN);
				
				
			mpfr_mul(yv, yv, temp1, GMP_RNDN);	
		
			mpz_t xvz;
			mpz_init2(xvz,1000);
			mpfr_get_z(xvz,yv,GMP_RNDN);
			r=mpz_class(xvz);
		
			if(ver==0)//if the number was negative..covert to binary
				{	
					mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI));
					r = tmpSUB -r;
					//r= r*(mpz_class(-1));
				}	
		
		}
	
	
		return r;
	}

}
