#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "../../utils.hpp"
#include "CoordinatesTableZ.hpp"
#include <stdio.h>
#include <mpfr.h>



using namespace std;
using std::ifstream;


namespace flopoco{



	CoordinatesTableZ::CoordinatesTableZ(Target* target, int wIn, int LSBI,int MSBI,char *filename) : 
		DualTable(target, wIn, (MSBI-LSBI)), MSBI(MSBI), LSBI(LSBI), adrWidth(wIn) , wOutm(MSBI-LSBI),filepath(filename)
   {
		if ((MSBI < LSBI)){
			cerr << 
				" CoordinatesTableZ: Input constraint LSBI <= MSBI not met."<<endl;
			exit (EXIT_FAILURE);
		}
	   
		ostringstream name; 
		name <<"CoordinatesX_"<<wIn<<"_"<<abs(LSBI)<<"_"<<abs(MSBI);
		setName(name.str());
	   
		readParams();
	
	}

	CoordinatesTableZ::~CoordinatesTableZ() {}
	

	int    CoordinatesTableZ::double2input(double x){
		int result;
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
		return result;
	}


	double CoordinatesTableZ::input2double(int x) {
		double y;
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
		return(y);
	}

	mpz_class CoordinatesTableZ::double2output(double x){
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
		return 0;
	}

	double CoordinatesTableZ::output2double(mpz_class x) {
		double y;
		cerr << "??? PolynomialTable::double2input not yet implemented ";
		exit(1);
  
		return(y);
	}


	void CoordinatesTableZ::readParams()
	{
		ifstream indata; // indata is like cin
		indata.open(filepath);
	
	
		indata>>L>>NrSVert;
		//cout<<L<<NrSVert;
	   
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
	void CoordinatesTableZ::initParameters()
	{
		//declaratii si initializari
		mpfr_set_default_prec(1000);
	
	
		mpfr_inits(out_rad_m,turn_rad_m,insu_m,nr_m,(mpfr_ptr) 0);
		mpfr_set_d(out_rad_m,out_radius,GMP_RNDN);
		mpfr_set_d(turn_rad_m,turn_radius,GMP_RNDN);
		mpfr_set_d(insu_m,insulation,GMP_RNDN);
	
		mpfr_t temp1,temp2,temp3,temp4;
		mpfr_inits( temp1,temp2,temp3,temp4,(mpfr_ptr) 0);
	
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
		
		for(int j=0;j<5;j++)
			{
				//mpfr_inits(displacemenX[j],displacementZ[j],alfa[j],(mpfr_ptr) 0);
				mpfr_init_set_si(displacementZ[j],(+0),GMP_RNDN);
				mpfr_init_set_si(alfa[j],(+0),GMP_RNDN);
			}
		switch(L){
		case 1:	
			//all were set to zero
			break;
		case 2:	
			
			//deplasamentZ[1]=(NrSVert-1)*(2.0*turn_radius+insulation);
			mpfr_mul_si(temp1,turn_rad_m,2,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_m,GMP_RNDN);
			mpfr_mul_si(displacementZ[1],temp1,(NrSVert-1),GMP_RNDN);
			
			//alfa[1]=M_PI;
			mpfr_const_pi(alfa[1],GMP_RNDN);
		
			break;
		case 3:	
			//deplasamentZ[1]=h+insulation;deplasamentZ[2]=h;     
			mpfr_add(displacementZ[1],h,insu_m,GMP_RNDN);
			mpfr_set(displacementZ[2],h,GMP_RNDN);
				
			//alfa[1]=M_PI/2.0;		alfa[2]=M_PI;
			mpfr_const_pi(temp1,GMP_RNDN);
			mpfr_div_ui(alfa[1],temp1,2,GMP_RNDN);		
			mpfr_const_pi(alfa[2],GMP_RNDN);
		
			break;
		case 4:	
			//deplasamentZ[1]=h+insulation;deplasamentZ[2]=h*3.0/2.0+insulation;         deplasamentZ[3]=h;   
			mpfr_add(displacementZ[1],h,insu_m,GMP_RNDN);
			mpfr_mul_ui(temp1,h,3,GMP_RNDN);
			mpfr_div_ui(temp1,temp1,2,GMP_RNDN);
			mpfr_add(displacementZ[2],temp1,insu_m,GMP_RNDN);
			mpfr_set(displacementZ[3],h,GMP_RNDN);
			
			//alfa[1]=M_PI/3.0;		alfa[2]=2.0*M_PI/3.0;		alfa[3]=M_PI;
			mpfr_const_pi(temp1,GMP_RNDN);
			mpfr_div_ui(temp2,temp1,3,GMP_RNDN);		
			mpfr_set(alfa[1],temp2,GMP_RNDN);
			mpfr_mul_ui(alfa[2],temp2,2,GMP_RNDN);
			mpfr_const_pi(alfa[3],GMP_RNDN);
		
			break;
		case 5:	
			//deplasamentZ[1]=h+insulation; deplasamentZ[2]=h*(1.0+sin(M_PI/4.0))+insulation;
			//deplasamentZ[3]=h*(1.0+sin(M_PI/4.0))+insulation;      deplasamentZ[4]=h;
			mpfr_add(displacementZ[1],h,insu_m,GMP_RNDN);
			mpfr_const_pi(temp1,GMP_RNDN);
			mpfr_div_ui(temp1,temp1,4,GMP_RNDN);	
			mpfr_sin(temp1,temp1,GMP_RNDN);
			mpfr_add_ui(temp1,temp1,1,GMP_RNDN);
			mpfr_mul(temp1,temp1,h,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_m,GMP_RNDN);
			mpfr_set(displacementZ[2],temp1,GMP_RNDN);
			mpfr_set(displacementZ[3],temp1,GMP_RNDN);
			mpfr_set(displacementZ[4],h,GMP_RNDN);
		
					
			//alfa[1]=M_PI/4.0;		alfa[2]=M_PI/2.0;		alfa[3]=3.0*M_PI/4.0;		alfa[4]=M_PI;
				
			mpfr_const_pi(temp1,GMP_RNDN);
			mpfr_div_ui(temp2,temp1,4,GMP_RNDN);		
			mpfr_set(alfa[1],temp2,GMP_RNDN);
				
			mpfr_mul_ui(alfa[2],temp2,2,GMP_RNDN);
			mpfr_mul_ui(alfa[3],temp2,3,GMP_RNDN);
			mpfr_const_pi(alfa[4],GMP_RNDN);
				
			break;
		}
		
		mpfr_clear(temp1);
		mpfr_clear(temp2);
		mpfr_clear(temp3);
		mpfr_clear(temp4);
	}

	mpz_class CoordinatesTableZ::function(int x)
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


		mpfr_t zv;
		mpfr_init2(zv,1000);
	
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
			
				
			//float temp=rk *(1.0+cos((q)*fi))+deltak ;
			mpfr_mul_ui(temp1,fi,q,GMP_RNDN);
			mpfr_cos(temp1,temp1,GMP_RNDN);
			mpfr_add_ui(temp1,temp1,1,GMP_RNDN);
			mpfr_mul(temp1,temp1,rk,GMP_RNDN);
			mpfr_add(temp,temp1,deltak,GMP_RNDN);
			/*float tauqj;
			  if(temp==0)
			  tauqj=M_PI/2;
			  else
			  tauqj=atan(ddj/temp);*/
			if(mpfr_sgn(temp)==0)
				{
					mpfr_const_pi(temp1,GMP_RNDN);
					mpfr_div_ui(tauqj,temp1,2,GMP_RNDN);	
				}
			else
				{
					mpfr_div(temp1,ddj,temp,GMP_RNDN);
					mpfr_atan(tauqj,temp1,GMP_RNDN);
				}
					
			//float zv=deplasamentZ[i]+sqrt(pow(ddj,2)+pow(temp,2))*sin(alfa[i]+tauqj);
			mpfr_hypot(temp1,ddj,temp,GMP_RNDN);
			mpfr_add(temp2,alfa[i],tauqj,GMP_RNDN);
			mpfr_sin(temp2,temp2,GMP_RNDN);
			mpfr_mul(temp1,temp1,temp2,GMP_RNDN);
			mpfr_add(zv,displacementZ[i],temp1,GMP_RNDN);
					
					
			ver =1;
		
			if(mpfr_sgn(zv)<0)
				{ver=0;
					mpfr_mul_si(zv,zv,(-1),GMP_RNDN);
				}
		
			mpfr_set_ui(temp1, 2 , GMP_RNDN);
			mpfr_set_si(temp2, abs(LSBI) , GMP_RNDN);
			mpfr_pow(temp1, temp1, temp2, GMP_RNDN);
				
				
			mpfr_mul(zv, zv, temp1, GMP_RNDN);	
		
			mpz_t xvz;
			mpz_init2(xvz,1000);
			mpfr_get_z(xvz,zv,GMP_RNDN);
			r=mpz_class(xvz);
		
			if(ver==0)//if the number was negative..covert to binary
				{	
					mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI));
					r = tmpSUB - r;
					//~ r= r*(mpz_class(-1));
				}	
		
		}
	
	
		return r;
	}

}
