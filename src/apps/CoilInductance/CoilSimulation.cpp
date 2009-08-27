#include <gmpxx.h>
#include <iostream>   
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include<math.h>

#include <stdio.h>
#include <mpfr.h>



using namespace std;
using std::ofstream;


void initParameters();
void initParameterst();
void initParametersz();
void initParametersy();

void readParams();
//mpz_class functionX(int x,mpfr_t &xv);
//void functionZ(int x,mpfr_t &zv);
//void functionY(int x,mpfr_t &yv);
mpz_class  functiont(int x,mpfr_t &xv);
mpz_class  functionzt(int x,mpfr_t &zv);
mpz_class  functionyt(int x,mpfr_t &yv);

string unsignedBinary(mpz_class x, int size);


mpfr_t out_rad_m,turn_rad_m,insu_m,nr_m;
mpfr_t fi;
mpfr_t h;
mpfr_t displacementX[5],alfa[5],displacementZ[5];


 mpfr_t out_rad_mt,turn_rad_mt,insu_mt,nr_mt;

 mpfr_t displacementXt[5],alfat[5];

mpfr_t out_rad_mz,turn_rad_mz,insu_mz,nr_mz;

	mpfr_t displacementZz[5],alfaz[5];

mpfr_t out_rad_my,turn_rad_my,insu_my,nr_my;
	mpfr_t fiy;
	mpfr_t hy;


	/** The MSB for the input */
	 int MSBI; 
	/** The LSB for the input */
	 int LSBI; 


	/** configuration file path from which to initialize the parameters of the coil*/
	 char *filepath;
	/** number of coils	*/
	 int L;
	/** number of vertical stages	*/
	 int NrSVert;
	/** configuration */
	//int  config [L][NrSVert];
	 int **config;
	/** number of turns */
	 int nrTurns;
	/** outer radius */
	 float out_radius;
	/* radius of a turn */
	 float turn_radius;
	/**	insulation	*/
	 float insulation;
	/** number of points to divide a turn */
	 int N;
	 
	 int prec;
	 int limit1;
	 int limit2;
	 int limit4t; 


int main(int argc, char* argv[] )
{
int i=1;
	
	LSBI=atoi(argv[i++]);
	MSBI=atoi(argv[i++]);
	
	filepath=argv[i++]; //filepath
	limit1=atoi(argv[i++]);
	limit2=atoi(argv[i++]);
	limit4t=atoi(argv[i++]);
	
	prec= MSBI-LSBI;
	//prec=atoi(argv[i++]); //precision
	
	
	mpfr_set_default_prec(prec);
	
	readParams();

	initParameterst();
	initParametersz();
	initParametersy();
	//initParametersX();
	//initParametersZ();
	
	
	
	mpfr_t cst2,cst1;
	mpfr_init(cst1);
	mpfr_init2(cst2,4);
	
	mpfr_set_si(cst2,2,GMP_RNDN);
	mpfr_set_si(cst1,2,GMP_RNDN);
	
	mpfr_pow_si(cst2,cst2,-3,GMP_RNDN);
	mpfr_pow_si(cst1,cst1,-9,GMP_RNDN);
	
	
	mpfr_t segX1mX0,segX3mX2,segX3mX0,segX0mX1,segX0mX2;
	mpfr_inits(segX1mX0,segX3mX2,segX3mX0,segX0mX1,segX0mX2,(mpfr_ptr) 0 );
	mpfr_t sigX0,sigX1,sigX2,sigX3;
	mpfr_inits(sigX0,sigX1,sigX2,sigX3,(mpfr_ptr) 0 );
	
	mpfr_t segY1mY0,segY3mY2,segY3mY0,segY0mY1,segY0mY2;
	mpfr_inits(segY1mY0,segY3mY2,segY3mY0,segY0mY1,segY0mY2,(mpfr_ptr) 0 );
	mpfr_t sigY0,sigY1,sigY2,sigY3;
	mpfr_inits(sigY0,sigY1,sigY2,sigY3,(mpfr_ptr) 0 );
	
	mpfr_t segZ1mZ0,segZ3mZ2,segZ3mZ0,segZ0mZ1,segZ0mZ2;
	mpfr_inits(segZ1mZ0,segZ3mZ2,segZ3mZ0,segZ0mZ1,segZ0mZ2,(mpfr_ptr) 0 );
	mpfr_t sigZ0,sigZ1,sigZ2,sigZ3;
	mpfr_inits(sigZ0,sigZ1,sigZ2,sigZ3,(mpfr_ptr) 0 );
	
	
	mpz_class temmpz;
	
	mpfr_t temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	mpfr_inits2(prec*2,temp1,temp2,temp3,temp4,temp5,temp6,temp7,(mpfr_ptr) 0);
	
	mpfr_t var1,var2,var3,var4,var5;
	mpfr_inits2(24,var1,var2,var3,var4,var5,(mpfr_ptr) 0);
	
	mpfr_t temp1_23;
	mpfr_init2(temp1_23,24);
	mpfr_t temp2_23;
	mpfr_init2(temp2_23,24);
	mpfr_t temp3_23;
	mpfr_init2(temp3_23,24);
	mpfr_t sigt;
	mpfr_init2(sigt,4);
	
	
	
	mpfr_t acc;
	mpfr_inits2(94,acc,(mpfr_ptr)0);
	mpfr_set_si(acc,0,GMP_RNDN);



	

	/*
	//temmpz=functionX(1,sigX0);
        cout<<"The value in mpz class is "<<unsignedBinary(temmpz,prec);
	mpfr_printf("   %Rf  \n",sigX0);
	temmpz=functiont(1,sigX0);
        cout<<"The value in mpz class is "<<unsignedBinary(temmpz,prec);
	mpfr_printf("   %Rf  \n",sigX0);

	
	//functionZ(1,sigZ0);
        
	mpfr_printf("\n\n   %Rf  \n",sigZ0);
	temmpz=functionzt(1,sigZ0);
        cout<<"The value in mpz class is "<<unsignedBinary(temmpz,prec);
	mpfr_printf("   %Rf  \n",sigZ0);


	//functionY(1,sigY0);
        
	mpfr_printf("\n\n   %Rf  \n",sigY0);
	temmpz=functionyt(1,sigY0);
        cout<<"The value in mpz class is "<<unsignedBinary(temmpz,prec);
	mpfr_printf("   %Rf  \n",sigY0);

	*/

	for(int t=0;t<limit4t;t++)
	{
		mpfr_set_ui(sigt,t,GMP_RNDN);
		mpfr_mul(sigt,sigt,cst2,GMP_RNDN);
		
		mpfr_printf("t:= %Rf   %.4Rb \n",sigt,sigt);
	
		for(int a1=0;a1<limit1;a1++)
		{	
			

			 temmpz=functiont(a1,sigX2);
			 //cout<<"The value in mpz class is "<<unsignedBinary(temmpz,prec);
			 //mpfr_printf("   %Rf  \n",sigX0);

			
			
			if(a1%64!=63)
			{
				temmpz=functiont(a1+1,sigX3);
				//cout<<"The value in mpz class is "<<unsignedBinary(temmpz,prec);
				//mpfr_printf("   %Rf  \n",sigX1);

			}
			else
			{
				 functiont(a1-63,sigX3);
				
			}
			
			
			 functionyt(a1,sigY2);
			
			if(a1%64!=63)
			{
				 functionyt(a1+1,sigY3);
				
			}
			else
			{
				 functionyt(a1-63,sigY3);
				
			}
			
			 functionzt(a1,sigZ2);
			
			if(a1%64!=63)
			{
				 functionzt(a1+1,sigZ3);
				
			}
			else
			{
				 functionzt(a1-63,sigZ3);
				
			}
			
			//cout<<"Intra la loop 2"<<endl;
			
			for(int a2=0;a2<limit2;a2++)
				if(abs(a1-a2)>1)
				{
					 functiont(a2,sigX0);
					if(a2%64!=63)
					{
						functiont(a2+1,sigX1 );
						
					}
					else
					{
						 functiont(a2-63,sigX1 );
						
					}
					
					//cout<<"gata cu x2,3"<<endl;
					
					//~ mpfr_mul(sigX0,sigX0,cst1,GMP_RNDN);
					//~ mpfr_mul(sigX1,sigX1,cst1,GMP_RNDN);
					//~ mpfr_mul(sigX2,sigX2,cst1,GMP_RNDN);
					//~ mpfr_mul(sigX3,sigX3,cst1,GMP_RNDN);
			
			
					mpfr_sub(segX1mX0,sigX1,sigX0,GMP_RNDN);
					mpfr_sub(segX0mX1,sigX0,sigX1,GMP_RNDN);
					mpfr_sub(segX3mX2,sigX3,sigX2,GMP_RNDN);
					mpfr_sub(segX3mX0,sigX3,sigX0,GMP_RNDN);
					mpfr_sub(segX0mX2,sigX0,sigX2,GMP_RNDN);
					
					
					 functionyt(a2,sigY0);
					if(a2%64!=63)
					{
						 functionyt(a2+1,sigY1);
						
					}
					else
					{
						 functionyt(a2-63,sigY1);
						
					}
										
					//~ mpfr_mul(sigY0,sigY0,cst1,GMP_RNDN);
					//~ mpfr_mul(sigY1,sigY1,cst1,GMP_RNDN);
					//~ mpfr_mul(sigY2,sigY2,cst1,GMP_RNDN);
					//~ mpfr_mul(sigY3,sigY3,cst1,GMP_RNDN);
			
					mpfr_sub(segY1mY0,sigY1,sigY0,GMP_RNDN);
					mpfr_sub(segY0mY1,sigY0,sigY1,GMP_RNDN);
					mpfr_sub(segY3mY2,sigY3,sigY2,GMP_RNDN);
					mpfr_sub(segY3mY0,sigY3,sigY0,GMP_RNDN);
					mpfr_sub(segY0mY2,sigY0,sigY2,GMP_RNDN);		
					
					
					 functionzt(a2,sigZ0);
					if(a2%64!=63)
					{
						functionzt(a2+1,sigZ1 );
						
					}
					else
					{
						 functionzt(a2-63,sigZ1);
						
					}
					
					//~ mpfr_mul(sigZ0,sigZ0,cst1,GMP_RNDN);
					//~ mpfr_mul(sigZ1,sigZ1,cst1,GMP_RNDN);
					//~ mpfr_mul(sigZ2,sigZ2,cst1,GMP_RNDN);
					//~ mpfr_mul(sigZ3,sigZ3,cst1,GMP_RNDN);
			
					mpfr_sub(segZ1mZ0,sigZ1,sigZ0,GMP_RNDN);
					mpfr_sub(segZ0mZ1,sigZ0,sigZ1,GMP_RNDN);
					mpfr_sub(segZ3mZ2,sigZ3,sigZ2,GMP_RNDN);
					mpfr_sub(segZ3mZ0,sigZ3,sigZ0,GMP_RNDN);
					mpfr_sub(segZ0mZ2,sigZ0,sigZ2,GMP_RNDN);
							
					
					
					//var1
			
					//~ mpfr_printf("%Rf %Rf  ",segX1mX0,segX3mX2);
					//~ mpfr_printf("%Rf %Rf  ",segY1mY0,segY3mY2);
					//~ mpfr_printf("%Rf %Rf \n",segZ1mZ0,segZ3mZ2);
					
					mpfr_mul(temp1,segX1mX0,segX3mX2,GMP_RNDN);
					//mpfr_printf("%Rf ",temp1);
					mpfr_mul(temp2,segY1mY0,segY3mY2,GMP_RNDN);
					//mpfr_printf("%Rf ",temp2);
					mpfr_mul(temp3,segZ1mZ0,segZ3mZ2,GMP_RNDN);
					//mpfr_printf("%Rf \n",temp3);
					mpfr_add(temp4,temp1,temp2,GMP_RNDN);
					mpfr_add(var1,temp4,temp3,GMP_RNDN);
			
					
					
					
					//var2
			
					mpfr_sqr(temp1,segX3mX2,GMP_RNDN);
					mpfr_sqr(temp2,segY3mY2,GMP_RNDN);
					mpfr_sqr(temp3,segZ3mZ2,GMP_RNDN);
			
					mpfr_add(temp4,temp1,temp2,GMP_RNDN);
					mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
					mpfr_sqrt(var2,temp1_23,GMP_RNDN);
			
			
					//var3
				mpfr_mul(temp1,segX0mX1,sigt,GMP_RNDN);
				mpfr_add(temp1,temp1,segX3mX0,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY0mY1,sigt,GMP_RNDN);
				mpfr_add(temp2,temp2,segY3mY0,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ0mZ1,sigt,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ3mZ0,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_set(var3,temp1_23,GMP_RNDN);
				
				
				
					//var4
				mpfr_mul(temp1,segX1mX0,sigt,GMP_RNDN);
				mpfr_add(temp1,temp1,segX0mX2,GMP_RNDN);
				mpfr_set(temp5,temp1,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY1mY0,sigt,GMP_RNDN);
				mpfr_add(temp2,temp2,segY0mY2,GMP_RNDN);
				mpfr_set(temp6,temp2,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ1mZ0,sigt,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ0mZ2,GMP_RNDN);
				mpfr_set(temp7,temp3,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_set(var4,temp1_23,GMP_RNDN);
				
				
					//var5
				
				mpfr_mul(temp1,temp5,segX3mX2,GMP_RNDN);
				mpfr_mul(temp2,temp6,segY3mY2,GMP_RNDN);
				mpfr_mul(temp3,temp7,segZ3mZ2,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_set(var5,temp1_23,GMP_RNDN);
			
				
			//mpfr_printf("%Rf  %Rf  %Rf  %Rf  %Rf\n",var1,var2,var3,var4,var5);
				
				//formula finala
			
			mpfr_div(temp3_23,var5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,var3,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp2_23,var4,temp3_23,GMP_RNDN);
			
			//mpfr_printf("%Rf  ",temp2_23);
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf ",temp1_23);
			mpfr_mul(temp1_23,temp1_23,cst2,GMP_RNDN);
			//mpfr_printf("%Rf\n\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
				mpfr_add(acc,acc,temp1_23,GMP_RNDN);
					
					
				}
		}
	
	}
	
	mpfr_printf("Final value with t on 3 bits is:=%Rf\n",acc);
	
	
}



void initParametersX()
{
	mpfr_t temp1,temp2,temp3,temp4;
	mpfr_inits( temp1,temp2,temp3,temp4,(mpfr_ptr) 0);
	
	for(int j=0;j<5;j++)
		{//mpfr_inits(displacemenX[j],displacementZ[j],alfa[j],(mpfr_ptr) 0);
		mpfr_init_set_si(displacementX[j],(+0),GMP_RNDN);
		mpfr_init_set_si(alfa[j],(+0),GMP_RNDN);
		}
	switch(L){
		case 1:	
				//all were set to zero
				break;
		case 2:	
			
				//deplasamentX[1]=-2*insulation;
				mpfr_mul_si(displacementX[1],insu_m,(-2),GMP_RNDN);
			
				//alfa[1]=M_PI;
				mpfr_const_pi(alfa[1],GMP_RNDN);
		
				break;
		case 3:	
				//deplasamentX[1]=-insulation;   deplasamentX[2]=-h-2*insulation;
				mpfr_set(displacementX[1],insu_m,GMP_RNDN);
				mpfr_mul_si(temp1,insu_m,2,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_sub(displacementX[2],temp2,temp1,GMP_RNDN);
				
				//alfa[1]=M_PI/2.0;		alfa[2]=M_PI;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(alfa[1],temp1,2,GMP_RNDN);		
				mpfr_const_pi(alfa[2],GMP_RNDN);
		
				break;
		case 4:	
				//deplasamentX[2]=-h*cos(M_PI/6.0)-insulation; deplasamentX[3]=-2.0*h*cos(M_PI/6.0)-insulation;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,6,GMP_RNDN);	
				mpfr_cos(temp1,temp1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp3,temp1,temp2,GMP_RNDN);
				mpfr_sub(displacementX[2],temp3,insu_m,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-2),GMP_RNDN);
				mpfr_mul(displacementX[3],temp2,temp1,GMP_RNDN);
			
				//alfa[1]=M_PI/3.0;		alfa[2]=2.0*M_PI/3.0;		alfa[3]=M_PI;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp2,temp1,3,GMP_RNDN);		
				mpfr_set(alfa[1],temp2,GMP_RNDN);
				mpfr_mul_ui(alfa[2],temp2,2,GMP_RNDN);
				mpfr_const_pi(alfa[3],GMP_RNDN);
		
				break;
		case 5:	
				//deplasamentX[1]=-insulation;    deplasamentX[2]=-h*cos(M_PI/4.0)-2*insulation; 
				//deplasamentX[3]=-h*(1.0 + cos(M_PI/4.0))-3*insulation;deplasamentX[4]=-h*(1.0 + 2.0*cos(M_PI/4.0))-4*insulation;
				mpfr_mul_si(temp3,insu_m,(-1),GMP_RNDN);
				mpfr_set(displacementX[1],temp3,GMP_RNDN);
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,4,GMP_RNDN);	
				mpfr_cos(temp1,temp1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp2,temp2,temp1,GMP_RNDN);
				mpfr_sub(temp3,temp3,insu_m,GMP_RNDN);
				mpfr_add(displacementX[2],temp2,temp3,GMP_RNDN);
				mpfr_add_ui(temp4,temp1,1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp2,temp2,temp4,GMP_RNDN);
				mpfr_sub(temp3,temp3,insu_m,GMP_RNDN);
				mpfr_add(displacementX[3],temp2,temp3,GMP_RNDN);
				mpfr_mul_ui(temp1,temp1,2,GMP_RNDN);
				mpfr_add_ui(temp4,temp1,1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp2,temp2,temp4,GMP_RNDN);
				mpfr_sub(temp3,temp3,insu_m,GMP_RNDN);
				mpfr_add(displacementX[4],temp2,temp3,GMP_RNDN);
		
					
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

void initParametersZ()
{
		mpfr_t temp1,temp2,temp3,temp4;
	mpfr_inits( temp1,temp2,temp3,temp4,(mpfr_ptr) 0);
	
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
			
						
				break;
		case 3:	
				//deplasamentZ[1]=h+insulation;deplasamentZ[2]=h;     
				mpfr_add(displacementZ[1],h,insu_m,GMP_RNDN);
				mpfr_set(displacementZ[2],h,GMP_RNDN);
			
		
				break;
		case 4:	
				//deplasamentZ[1]=h+insulation;deplasamentZ[2]=h*3.0/2.0+insulation;         deplasamentZ[3]=h;   
				mpfr_add(displacementZ[1],h,insu_m,GMP_RNDN);
				mpfr_mul_ui(temp1,h,3,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,2,GMP_RNDN);
				mpfr_add(displacementZ[2],temp1,insu_m,GMP_RNDN);
				mpfr_set(displacementZ[3],h,GMP_RNDN);
			
		
		
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
		
					
								
				break;
		}
		
		mpfr_clear(temp1);
		mpfr_clear(temp2);
		mpfr_clear(temp3);
		mpfr_clear(temp4);	
	
}


void initParameters()
{
	//declaratii si initializari
	mpfr_set_default_prec(prec);
	
	
	mpfr_inits(out_rad_m,turn_rad_m,insu_m,nr_m,(mpfr_ptr) 0);
	mpfr_set_d(out_rad_m,out_radius,GMP_RNDN);
	mpfr_set_d(turn_rad_m,turn_radius,GMP_RNDN);
	mpfr_set_d(insu_m,insulation,GMP_RNDN);
	
	mpfr_t temp1,temp2;
	mpfr_inits( temp1,temp2,(mpfr_ptr) 0);
	
	//float h= ((float)(NrSVert-1))*(2.0*turn_radius+insulation);
	
	
	mpfr_init2(h,prec);
	
	mpfr_mul_si(temp1,turn_rad_m,2,GMP_RNDN);
	mpfr_add(temp1,temp1,insu_m,GMP_RNDN);
	mpfr_mul_si(h,temp1,(NrSVert-1),GMP_RNDN);
	
	//float fi=2*M_PI/N;
	
	
	
	mpfr_init2(fi,prec);
	mpfr_const_pi(temp2,GMP_RNDN);
	mpfr_div_si(temp1,temp2,N,GMP_RNDN);
	mpfr_mul_ui(fi,temp1,2,GMP_RNDN);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	
	
}


void readParams()
{
ifstream indata; 
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


void functionY(int x,mpfr_t &yv)
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
					
		
		
			
		
	}
	
	mpfr_clear(ddj);
	mpfr_clear(deltak);
	mpfr_clear(rk);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	//mpfr_clear(yv);
	mpfr_clear(temp);
	mpfr_clear(tauqj);
		mpfr_free_cache();
	
}



void functionZ(int x,mpfr_t &zv)
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
		
		
		
	}
	

	mpfr_clear(ddj);
	mpfr_clear(deltak);
	mpfr_clear(rk);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	//mpfr_clear(zv);
	mpfr_clear(temp);
	mpfr_clear(tauqj);
		mpfr_free_cache();
	
}


	
mpz_class functionX(int x,mpfr_t &xv)
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
					
					//float xv=deplasamentX[i]+sqrt(pow(ddj,2)+pow(temp,2))*cos(alfa[i]+tauqj);
					mpfr_hypot(temp1,ddj,temp,GMP_RNDN);
					
					
					mpfr_add(temp2,alfa[i],tauqj,GMP_RNDN);
					mpfr_cos(temp2,temp2,GMP_RNDN);
					mpfr_mul(temp1,temp1,temp2,GMP_RNDN);
					mpfr_add(xv,displacementX[i],temp1,GMP_RNDN);
					
					
		ver =1;
		
		if(mpfr_sgn(xv)<0)
		{ver=0;
		mpfr_mul_si(xv,xv,(-1),GMP_RNDN);
		}
		
		mpfr_set_ui(temp1, 2 , GMP_RNDN);
		mpfr_set_si(temp2, abs(LSBI) , GMP_RNDN);
		mpfr_pow(temp1, temp1, temp2, GMP_RNDN);
		
		mpz_t xvz;
		mpz_init2(xvz,1000);		
				
		mpfr_mul(temp1, xv, temp1, GMP_RNDN);	
		
		
		mpfr_get_z(xvz,temp1,GMP_RNDN);
		r=mpz_class(xvz);
		mpz_clear(xvz);
		
		if(ver==0)//if the number was negative..covert to binary
		{	
		mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI+1));
		r =r - tmpSUB;
		r= r*(mpz_class(-1));
		}
		
		
	}
	
	mpfr_clear(ddj);
	mpfr_clear(deltak);
	mpfr_clear(rk);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	//mpfr_clear(xv);
	mpfr_free_cache();
	mpfr_clear(temp);
	mpfr_clear(tauqj);
	
	mpfr_free_cache();
	return r;
}

string unsignedBinary(mpz_class x, int size){
	string s;
	mpz_class po2, number;
	char bit;

	if(x<0) {
		cerr<<"unsigned_binary: Positive number expected, got x="<<x.get_d()<<endl;
		exit(EXIT_FAILURE);
	}
	po2 = ((mpz_class) 1)<<size;
	number=x;
		
	for (int i = 0; i < size ; i++) {
		po2 = po2>>1;
		if (number >= po2) {
			bit = '1';
			number -= po2;
		}
		else {
			bit = '0';
		}
		s +=  bit;
	}
	return s;
}

void initParameterst()
{
	
	//cout<<"initializeaza";
	
	//declaratii si initializari
	mpfr_set_default_prec(1000);
	
	
	mpfr_inits(out_rad_mt,turn_rad_mt,insu_mt,nr_mt,(mpfr_ptr) 0);
	mpfr_set_d(out_rad_mt,out_radius,GMP_RNDN);
	mpfr_set_d(turn_rad_mt,turn_radius,GMP_RNDN);
	mpfr_set_d(insu_mt,insulation,GMP_RNDN);
	
	mpfr_t temp1,temp2,temp3,temp4;
	mpfr_inits( temp1,temp2,temp3,temp4,(mpfr_ptr) 0);
	
	//float h= ((float)(NrSVert-1))*(2.0*turn_radius+insulation);
	
	
	mpfr_init2(h,1000);
	
	mpfr_mul_si(temp1,turn_rad_mt,2,GMP_RNDN);
	mpfr_add(temp1,temp1,insu_mt,GMP_RNDN);
	mpfr_mul_si(h,temp1,(NrSVert-1),GMP_RNDN);
	
	//float fi=2*M_PI/N;
	
	
	
	mpfr_init2(fi,1000);
	mpfr_const_pi(temp2,GMP_RNDN);
	mpfr_div_si(temp1,temp2,N,GMP_RNDN);
	mpfr_mul_ui(fi,temp1,2,GMP_RNDN);
		
	for(int j=0;j<5;j++)
		{//mpfr_inits(displacemenX[j],displacementZ[j],alfat[j],(mpfr_ptr) 0);
		mpfr_init_set_si(displacementXt[j],(+0),GMP_RNDN);
		mpfr_init_set_si(alfat[j],(+0),GMP_RNDN);
		}
	switch(L){
		case 1:	
				//all were set to zero
				break;
		case 2:	
			
				//deplasamentX[1]=-2*insulation;
				mpfr_mul_si(displacementXt[1],insu_mt,(-2),GMP_RNDN);
			
				//alfat[1]=M_PI;
				mpfr_const_pi(alfat[1],GMP_RNDN);
		
				break;
		case 3:	
				//deplasamentX[1]=-insulation;   deplasamentX[2]=-h-2*insulation;
				mpfr_set(displacementXt[1],insu_mt,GMP_RNDN);
				mpfr_mul_si(temp1,insu_mt,2,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_sub(displacementXt[2],temp2,temp1,GMP_RNDN);
				
				//alfat[1]=M_PI/2.0;		alfat[2]=M_PI;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(alfat[1],temp1,2,GMP_RNDN);		
				mpfr_const_pi(alfat[2],GMP_RNDN);
		
				break;
		case 4:	
				//deplasamentX[2]=-h*cos(M_PI/6.0)-insulation; deplasamentX[3]=-2.0*h*cos(M_PI/6.0)-insulation;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,6,GMP_RNDN);	
				mpfr_cos(temp1,temp1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp3,temp1,temp2,GMP_RNDN);
				mpfr_sub(displacementXt[2],temp3,insu_mt,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-2),GMP_RNDN);
				mpfr_mul(displacementXt[3],temp2,temp1,GMP_RNDN);
			
				//alfat[1]=M_PI/3.0;		alfat[2]=2.0*M_PI/3.0;		alfat[3]=M_PI;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp2,temp1,3,GMP_RNDN);		
				mpfr_set(alfat[1],temp2,GMP_RNDN);
				mpfr_mul_ui(alfat[2],temp2,2,GMP_RNDN);
				mpfr_const_pi(alfat[3],GMP_RNDN);
		
				break;
		case 5:	
				//deplasamentX[1]=-insulation;    deplasamentX[2]=-h*cos(M_PI/4.0)-2*insulation; 
				//deplasamentX[3]=-h*(1.0 + cos(M_PI/4.0))-3*insulation;deplasamentX[4]=-h*(1.0 + 2.0*cos(M_PI/4.0))-4*insulation;
				mpfr_mul_si(temp3,insu_mt,(-1),GMP_RNDN);
				mpfr_set(displacementXt[1],temp3,GMP_RNDN);
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,4,GMP_RNDN);	
				mpfr_cos(temp1,temp1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp2,temp2,temp1,GMP_RNDN);
				mpfr_sub(temp3,temp3,insu_mt,GMP_RNDN);
				mpfr_add(displacementXt[2],temp2,temp3,GMP_RNDN);
				mpfr_add_ui(temp4,temp1,1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp2,temp2,temp4,GMP_RNDN);
				mpfr_sub(temp3,temp3,insu_mt,GMP_RNDN);
				mpfr_add(displacementXt[3],temp2,temp3,GMP_RNDN);
				mpfr_mul_ui(temp1,temp1,2,GMP_RNDN);
				mpfr_add_ui(temp4,temp1,1,GMP_RNDN);
				mpfr_mul_si(temp2,h,(-1),GMP_RNDN);
				mpfr_mul(temp2,temp2,temp4,GMP_RNDN);
				mpfr_sub(temp3,temp3,insu_mt,GMP_RNDN);
				mpfr_add(displacementXt[4],temp2,temp3,GMP_RNDN);
		
					
				//alfat[1]=M_PI/4.0;		alfat[2]=M_PI/2.0;		alfat[3]=3.0*M_PI/4.0;		alfat[4]=M_PI;
				
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp2,temp1,4,GMP_RNDN);		
				mpfr_set(alfat[1],temp2,GMP_RNDN);
				
				mpfr_mul_ui(alfat[2],temp2,2,GMP_RNDN);
				mpfr_mul_ui(alfat[3],temp2,3,GMP_RNDN);
				mpfr_const_pi(alfat[4],GMP_RNDN);
				
				break;
		}
		
		mpfr_clear(temp1);
		mpfr_clear(temp2);
		mpfr_clear(temp3);
		mpfr_clear(temp4);
		mpfr_free_cache();
}



mpz_class  functiont(int x,mpfr_t &xv)
{
	
	//cout<<"aiai";	
	
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


	
	//mpfr_t xv;
	//mpfr_init2(xv,1000);
	
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
				
			//mpfr_printf(" turnrad %Rf insu %Rf\n",turn_rad_mt,insu_mt);
			//ddj=((float)j)*(2.0*turn_radius+insulation);
			mpfr_mul_ui(temp1,turn_rad_mt,2,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_mt,GMP_RNDN);
			mpfr_mul_ui(ddj,temp1,j,GMP_RNDN);
							
				
				//float deltak=((float)k)*(2.0*turn_radius+insulation);
				mpfr_mul_ui(temp1,turn_rad_mt,2,GMP_RNDN);
				mpfr_add(temp1,temp1,insu_mt,GMP_RNDN);
				mpfr_mul_ui(deltak,temp1,k,GMP_RNDN);
				//float rk=out_radius-deltak;
				mpfr_sub(rk,out_rad_mt,deltak,GMP_RNDN);
			
					
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
					
					//float xv=deplasamentX[i]+sqrt(pow(ddj,2)+pow(temp,2))*cos(alfa[i]+tauqj);
					mpfr_hypot(temp1,ddj,temp,GMP_RNDN);
							
					
					//pt test
					//mpfr_sqr(xv,temp,GMP_RNDN);
					//mpfr_set(xv,temp1,GMP_RNDN);
					
					//mpfr_fprintf(stdout,"A:= %f  \n",mpfr_get_d(alfa[i],GMP_RNDN));
					//mpfr_fprintf(stdout,"T:= %f  \n",mpfr_get_d(tauqj,GMP_RNDN));
					
					mpfr_add(temp2,alfat[i],tauqj,GMP_RNDN);
					mpfr_cos(temp2,temp2,GMP_RNDN);
					mpfr_mul(temp1,temp1,temp2,GMP_RNDN);
					mpfr_add(xv,displacementXt[i],temp1,GMP_RNDN);
					
					
		
		ver =1;
		

		if(mpfr_sgn(xv)<0)
		{ver=0;
		mpfr_mul_si(xv,xv,(-1),GMP_RNDN);
		}
		
		mpfr_set_ui(temp1, 2 , GMP_RNDN);
		mpfr_set_si(temp2, abs(LSBI) , GMP_RNDN);
		mpfr_pow(temp1, temp1, temp2, GMP_RNDN);
				
				
		mpfr_mul(temp1, xv, temp1, GMP_RNDN);	
		
		mpz_t xvz;
		mpz_init2(xvz,1000);
		mpfr_get_z(xvz,temp1,GMP_RNDN);
		r=mpz_class(xvz);
		
		if(ver==0)//if the number was negative..covert to binary
		{	
		mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI+1));
		r =r - tmpSUB;
		r= r*(mpz_class(-1));
		}	
		
	}
	
	mpfr_clear(ddj);
	mpfr_clear(deltak);
	mpfr_clear(rk);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	//mpfr_clear(xv);
	//mpfr_free_cache();
	mpfr_clear(temp);
	mpfr_clear(tauqj);
	
	mpfr_free_cache();
	return r;
}

void initParametersz()
{
	//declaratii si initializari
	mpfr_set_default_prec(1000);
	
	
	mpfr_inits(out_rad_mz,turn_rad_mz,insu_mz,nr_mz,(mpfr_ptr) 0);
	mpfr_set_d(out_rad_mz,out_radius,GMP_RNDN);
	mpfr_set_d(turn_rad_mz,turn_radius,GMP_RNDN);
	mpfr_set_d(insu_mz,insulation,GMP_RNDN);
	
	mpfr_t temp1,temp2,temp3,temp4;
	mpfr_inits( temp1,temp2,temp3,temp4,(mpfr_ptr) 0);
	
	//float h= ((float)(NrSVert-1))*(2.0*turn_radius+insulation);
	
	
	mpfr_init2(h,1000);
	
	mpfr_mul_si(temp1,turn_rad_mz,2,GMP_RNDN);
	mpfr_add(temp1,temp1,insu_mz,GMP_RNDN);
	mpfr_mul_si(h,temp1,(NrSVert-1),GMP_RNDN);
	
	//float fi=2*M_PI/N;
	
	
	
	mpfr_init2(fi,1000);
	mpfr_const_pi(temp2,GMP_RNDN);
	mpfr_div_si(temp1,temp2,N,GMP_RNDN);
	mpfr_mul_ui(fi,temp1,2,GMP_RNDN);
		
	for(int j=0;j<5;j++)
		{
		//mpfr_inits(displacemenX[j],displacementZz[j],alfaz[j],(mpfr_ptr) 0);
		mpfr_init_set_si(displacementZz[j],(+0),GMP_RNDN);
		mpfr_init_set_si(alfaz[j],(+0),GMP_RNDN);
		}
	switch(L){
		case 1:	
				//all were set to zero
				break;
		case 2:	
			
				//deplasamentZ[1]=(NrSVert-1)*(2.0*turn_radius+insulation);
				mpfr_mul_si(temp1,turn_rad_mz,2,GMP_RNDN);
				mpfr_add(temp1,temp1,insu_mz,GMP_RNDN);
				mpfr_mul_si(displacementZz[1],temp1,(NrSVert-1),GMP_RNDN);
			
				//alfaz[1]=M_PI;
				mpfr_const_pi(alfaz[1],GMP_RNDN);
		
				break;
		case 3:	
				//deplasamentZ[1]=h+insulation;deplasamentZ[2]=h;     
				mpfr_add(displacementZz[1],h,insu_mz,GMP_RNDN);
				mpfr_set(displacementZz[2],h,GMP_RNDN);
				
				//alfaz[1]=M_PI/2.0;		alfaz[2]=M_PI;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(alfaz[1],temp1,2,GMP_RNDN);		
				mpfr_const_pi(alfaz[2],GMP_RNDN);
		
				break;
		case 4:	
				//deplasamentZ[1]=h+insulation;deplasamentZ[2]=h*3.0/2.0+insulation;         deplasamentZ[3]=h;   
				mpfr_add(displacementZz[1],h,insu_mz,GMP_RNDN);
				mpfr_mul_ui(temp1,h,3,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,2,GMP_RNDN);
				mpfr_add(displacementZz[2],temp1,insu_mz,GMP_RNDN);
				mpfr_set(displacementZz[3],h,GMP_RNDN);
			
				//alfaz[1]=M_PI/3.0;		alfaz[2]=2.0*M_PI/3.0;		alfaz[3]=M_PI;
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp2,temp1,3,GMP_RNDN);		
				mpfr_set(alfaz[1],temp2,GMP_RNDN);
				mpfr_mul_ui(alfaz[2],temp2,2,GMP_RNDN);
				mpfr_const_pi(alfaz[3],GMP_RNDN);
		
				break;
		case 5:	
				//deplasamentZ[1]=h+insulation; deplasamentZ[2]=h*(1.0+sin(M_PI/4.0))+insulation;
				//deplasamentZ[3]=h*(1.0+sin(M_PI/4.0))+insulation;      deplasamentZ[4]=h;
				mpfr_add(displacementZz[1],h,insu_mz,GMP_RNDN);
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp1,temp1,4,GMP_RNDN);	
				mpfr_sin(temp1,temp1,GMP_RNDN);
				mpfr_add_ui(temp1,temp1,1,GMP_RNDN);
				mpfr_mul(temp1,temp1,h,GMP_RNDN);
				mpfr_add(temp1,temp1,insu_mz,GMP_RNDN);
				mpfr_set(displacementZz[2],temp1,GMP_RNDN);
				mpfr_set(displacementZz[3],temp1,GMP_RNDN);
				mpfr_set(displacementZz[4],h,GMP_RNDN);
		
					
				//alfaz[1]=M_PI/4.0;		alfaz[2]=M_PI/2.0;		alfaz[3]=3.0*M_PI/4.0;		alfaz[4]=M_PI;
				
				mpfr_const_pi(temp1,GMP_RNDN);
				mpfr_div_ui(temp2,temp1,4,GMP_RNDN);		
				mpfr_set(alfaz[1],temp2,GMP_RNDN);
				
				mpfr_mul_ui(alfaz[2],temp2,2,GMP_RNDN);
				mpfr_mul_ui(alfaz[3],temp2,3,GMP_RNDN);
				mpfr_const_pi(alfaz[4],GMP_RNDN);
				
				break;
		}
		
		mpfr_clear(temp1);
		mpfr_clear(temp2);
		mpfr_clear(temp3);
		mpfr_clear(temp4);

	mpfr_free_cache();
}


mpz_class functionzt(int x,mpfr_t &zv)
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
			mpfr_mul_ui(temp1,turn_rad_mz,2,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_mz,GMP_RNDN);
			mpfr_mul_ui(ddj,temp1,j,GMP_RNDN);
			
				
				//float deltak=((float)k)*(2.0*turn_radius+insulation);
				mpfr_mul_ui(temp1,turn_rad_mz,2,GMP_RNDN);
				mpfr_add(temp1,temp1,insu_mz,GMP_RNDN);
				mpfr_mul_ui(deltak,temp1,k,GMP_RNDN);
				//float rk=out_radius-deltak;
				mpfr_sub(rk,out_rad_mz,deltak,GMP_RNDN);
			
				
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
					
					//float zv=deplasamentZ[i]+sqrt(pow(ddj,2)+pow(temp,2))*sin(alfaz[i]+tauqj);
					mpfr_hypot(temp1,ddj,temp,GMP_RNDN);
					mpfr_add(temp2,alfaz[i],tauqj,GMP_RNDN);
					mpfr_sin(temp2,temp2,GMP_RNDN);
					mpfr_mul(temp1,temp1,temp2,GMP_RNDN);
					mpfr_add(zv,displacementZz[i],temp1,GMP_RNDN);
					
					
		ver =1;
		
		if(mpfr_sgn(zv)<0)
		{ver=0;
		mpfr_mul_si(zv,zv,(-1),GMP_RNDN);
		}
		
		mpfr_set_ui(temp1, 2 , GMP_RNDN);
		mpfr_set_si(temp2, abs(LSBI) , GMP_RNDN);
		mpfr_pow(temp1, temp1, temp2, GMP_RNDN);
				
				
		mpfr_mul(temp1, zv, temp1, GMP_RNDN);	
		
		mpz_t xvz;
		mpz_init2(xvz,1000);
		mpfr_get_z(xvz,temp1,GMP_RNDN);
		r=mpz_class(xvz);
		
		if(ver==0)//if the number was negative..covert to binary
		{	
		mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI+1));
		r =r - tmpSUB;
		r= r*(mpz_class(-1));
		}	
		
	}
	
mpfr_clear(ddj);
	mpfr_clear(deltak);
	mpfr_clear(rk);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	//mpfr_clear(zv);
	mpfr_clear(temp);
	mpfr_clear(tauqj);
		mpfr_free_cache();
	
	return r;
}


void initParametersy()
{
	//declaratii si initializari
	mpfr_set_default_prec(1000);
	
	
	mpfr_inits(out_rad_my,turn_rad_my,insu_my,nr_my,(mpfr_ptr) 0);
	mpfr_set_d(out_rad_my,out_radius,GMP_RNDN);
	mpfr_set_d(turn_rad_my,turn_radius,GMP_RNDN);
	mpfr_set_d(insu_my,insulation,GMP_RNDN);
	
	mpfr_t temp1,temp2;
	mpfr_inits( temp1,temp2,(mpfr_ptr) 0);
	
	//float h= ((float)(NrSVert-1))*(2.0*turn_radius+insulation);
	
	
	mpfr_init2(h,1000);
	
	mpfr_mul_si(temp1,turn_rad_my,2,GMP_RNDN);
	mpfr_add(temp1,temp1,insu_my,GMP_RNDN);
	mpfr_mul_si(h,temp1,(NrSVert-1),GMP_RNDN);
	
	//float fi=2*M_PI/N;
	
	
	
	mpfr_init2(fi,1000);
	mpfr_const_pi(temp2,GMP_RNDN);
	mpfr_div_si(temp1,temp2,N,GMP_RNDN);
	mpfr_mul_ui(fi,temp1,2,GMP_RNDN);

	mpfr_clear(temp1);
	mpfr_clear(temp2);
		mpfr_free_cache();
		
	
}

mpz_class functionyt(int x,mpfr_t &yv)
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
			mpfr_mul_ui(temp1,turn_rad_my,2,GMP_RNDN);
			mpfr_add(temp1,temp1,insu_my,GMP_RNDN);
			mpfr_mul_ui(ddj,temp1,j,GMP_RNDN);
			
				
				//float deltak=((float)k)*(2.0*turn_radius+insulation);
				mpfr_mul_ui(temp1,turn_rad_my,2,GMP_RNDN);
				mpfr_add(temp1,temp1,insu_my,GMP_RNDN);
				mpfr_mul_ui(deltak,temp1,k,GMP_RNDN);
				//float rk=out_radius-deltak;
				mpfr_sub(rk,out_rad_my,deltak,GMP_RNDN);
			
				
										
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
				
				
		mpfr_mul(temp1, yv, temp1, GMP_RNDN);	
		
		mpz_t xvz;
		mpz_init2(xvz,1000);
		mpfr_get_z(xvz,temp1,GMP_RNDN);
		r=mpz_class(xvz);
		
		if(ver==0)//if the number was negative..covert to binary
		{	
		mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI+1));
		r =r - tmpSUB;
		r= r*(mpz_class(-1));
		}	
		
	}
	mpfr_clear(ddj);
	mpfr_clear(deltak);
	mpfr_clear(rk);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	//mpfr_clear(yv);
	mpfr_clear(temp);
	mpfr_clear(tauqj);
		mpfr_free_cache();
	
	return r;
}


