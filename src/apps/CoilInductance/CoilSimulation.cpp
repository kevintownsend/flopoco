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
void initParametersX();
void initParametersZ();
void readParams();
mpz_class functionX(int x);
mpz_class functionY(int x);
mpz_class functionZ(int x);
string unsignedBinary(mpz_class x, int size);

	
mpfr_t out_rad_m,turn_rad_m,insu_m,nr_m;
mpfr_t fi;
mpfr_t h;
mpfr_t displacementX[5],alfa[5],displacementZ[5];


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
	 int limit;



int main(int argc, char* argv[] )
{
int i=1;
	
	LSBI=atoi(argv[i++]);
	MSBI=atoi(argv[i++]);
	
	filepath=argv[i++]; //filepath
	limit=atoi(argv[i++]);
	
	prec= MSBI-LSBI;
	//prec=atoi(argv[i++]); //precision
	
	
	mpfr_set_default_prec(prec);
	
	readParams();

	//initParameters();
	initParametersX();
	initParametersZ();
	
	mpz_class tmpSUB = (mpz_class(1) << (MSBI-LSBI+1));
	mpz_class tmpCMP = (mpz_class(1)  << (MSBI-LSBI))-1;
	
	
	//~ mpfr_t te1;
	//~ mpfr_init2(te1,24);
	//~ mpfr_set_ld(te1,1.1,GMP_RNDN);
	//~ cout<<"ceva";
	//~ mpz_class temmpz1;
	//~ for(int oo1=0;oo1<limit;oo1++)
		//~ for(int oo2=0;oo2<limit;oo2++)
		//~ {
			//~ mpfr_mul(te1,te1,te1,GMP_RNDN);
			//~ mpfr_pow_si(te1,te1,3,GMP_RNDN);
			//~ mpfr_sqrt(te1,te1,GMP_RNDN);
			
			//~ temmpz1 = functionX(oo1);
			//~ if (temmpz1  > tmpCMP){ //negative number 
				//~ temmpz1  = temmpz1  - tmpSUB;
			//~ }
			//~ mpfr_set_z(te1,temmpz1.get_mpz_t(),GMP_RNDN);
			
			//~ }
	//~ cout<<"gata";
	//~ return 0;
			
			
	mpfr_t cst1,cst2,cst3,cst4,cst5,cst6;
	mpfr_init(cst1);
	mpfr_init2(cst2,3);
	mpfr_init2(cst3,3);
	mpfr_init2(cst4,3);
	mpfr_init2(cst5,4);
	
	mpfr_set_si(cst1,2,GMP_RNDN);
	mpfr_set_si(cst2,2,GMP_RNDN);
	mpfr_set_si(cst3,2,GMP_RNDN);
	mpfr_set_si(cst4,2,GMP_RNDN);
	mpfr_set_si(cst5,2,GMP_RNDN);
	
	mpfr_pow_si(cst1,cst1,-9,GMP_RNDN);
	mpfr_pow_si(cst2,cst2,-3,GMP_RNDN);
	mpfr_pow_si(cst3,cst3,-2,GMP_RNDN);
	mpfr_pow_si(cst4,cst4,-1,GMP_RNDN);
	mpfr_pow_si(cst5,cst5,-4,GMP_RNDN);
	
	


	
	
	//~ srand(prec);
	//~ for(i=0;i<21;i++)
	//~ {
		//~ int tt=rand()%8192;
		//~ cout<<tt<<"X:= ";
		//~ mpz_class v;
		//~ v=functionX(tt);
		//~ //cout<<unsignedBinary(functionX(tt),prec);
		//~ if (v > tmpCMP){ //negative number 
			//~ v = v - tmpSUB;
		//~ }
		
		//~ mpfr_t x;
		//~ mpfr_init(x);
		//~ mpfr_set_z(x,v.get_mpz_t(),GMP_RNDN);
		//~ mpfr_mul(x,x,cst1,GMP_RNDN);
		//~ mpfr_printf("%Rf",x);
		
		
		//~ cout<<" ; Y:= ";
		//~ v=functionY(tt);
		//~ if (v > tmpCMP){ //negative number 
			//~ v = v - tmpSUB;
		//~ }
		//~ mpfr_set_z(x,v.get_mpz_t(),GMP_RNDN);
		//~ mpfr_mul(x,x,cst1,GMP_RNDN);
		//~ mpfr_printf("%Rf",x);
		//~ cout<<" ; Z:= ";
		//~ v=functionZ(tt);
		//~ if (v > tmpCMP){ //negative number 
			//~ v = v - tmpSUB;
		//~ }
		//~ mpfr_set_z(x,v.get_mpz_t(),GMP_RNDN);
		//~ mpfr_mul(x,x,cst1,GMP_RNDN);
		//~ mpfr_printf("%Rf",x);
		
		//~ cout<<endl;
		//~ }
		
		
	/// Segments computation
	
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
	mpfr_t var1m,var2m,var3m,var4m,var5m;
	mpfr_inits2(24,var1m,var2m,var3m,var4m,var5m,(mpfr_ptr) 0);
	mp_exp_t minV1 = 255,maxV1 = -255,minV2= 255,maxV2= -255,minV3= 255,maxV3= -255,minV4= 255,maxV4= -255,minV5= 255,maxV5= -255,texp;
	mpfr_t miV1,maV1,miV2,maV2,miV3,maV3,miV4,maV4,miV5,maV5;
	mpfr_inits2(24,miV1,maV1,miV2,maV2,miV3,maV3,miV4,maV4,miV5,maV5,(mpfr_ptr)0);
	mpfr_set_si(var1m,-255,GMP_RNDN);
	mpfr_set_si(var2m,-255,GMP_RNDN);
	mpfr_set_si(var3m,-255,GMP_RNDN);
	mpfr_set_si(var4m,-255,GMP_RNDN);
	mpfr_set_si(var5m,-255,GMP_RNDN);
	
	mpfr_t temp1_23;
	mpfr_init2(temp1_23,24);
	mpfr_t temp2_23;
	mpfr_init2(temp2_23,24);
	mpfr_t temp3_23;
	mpfr_init2(temp3_23,24);
	mpfr_t sigt;
	mpfr_init2(sigt,3);
	mpfr_t sigt2;
	mpfr_init2(sigt2,4);
	mpfr_set_si(sigt2,7,GMP_RNDN);
	mpfr_mul(sigt2,sigt2,cst5,GMP_RNDN);
	//mpfr_printf("%Rf\n",sigt2);
	
	mpfr_init2(cst6,16);
	mpfr_set_si(cst6,2,GMP_RNDN);
	mpfr_pow_si(cst6,cst6,-15,GMP_RNDN);
	
	mpfr_t sigt3;
	mpfr_init2(sigt3,16);
	mpfr_set_si(sigt3,14349,GMP_RNDN);
	mpfr_mul(sigt3,sigt3,cst6,GMP_RNDN);
	
	mpfr_t sigt4;
	mpfr_init2(sigt4,5);
	
	//mpfr_printf("%Rf   %Rb\n",sigt3,sigt3);
	
	
	
	mpfr_t acc,acc1,acc2,acc31,acc21,acc4,accn1,accn2,accn3,accn4;
	mpfr_inits2(64,acc,acc1,acc2,acc31,acc21,acc4,accn1,accn2,accn3,accn4,(mpfr_ptr)0);
	mpfr_set_si(acc,0,GMP_RNDN);
	mpfr_set_si(acc1,0,GMP_RNDN);
	mpfr_set_si(acc2,0,GMP_RNDN);
	mpfr_set_si(acc4,0,GMP_RNDN);
	mpfr_set_si(acc31,0,GMP_RNDN);
	mpfr_set_si(acc21,0,GMP_RNDN);
	mpfr_set_si(accn1,0,GMP_RNDN);
	mpfr_set_si(accn2,0,GMP_RNDN);
	mpfr_set_si(accn3,0,GMP_RNDN);
	mpfr_set_si(accn4,0,GMP_RNDN);
	
	
	
	//~ minV1=mpfr_get_exp(cst1);
	//~ if(minV1<2)
	//~ printf("\n%d\n",((int)minV1));
	//mpfr_printf("aaaaaa");
	mpfr_t accvar3,accvar4,accvar5,accvar3_2,accvar3_1,maxD1,maxD2,accvar4_2,accvar4_1,accvar5_2,accvar5_1,accvar3_d2,accvar4_d2,var3p1,var4p1,accvar4p4,accvar3p4,accvar4p5,accvar3p5;
	mpfr_inits2(24,accvar3,accvar4,accvar5,accvar3_2,accvar3_1,maxD1,maxD2,accvar4_2,accvar4_1,accvar5_2,accvar5_1,accvar3_d2,accvar4_d2,var3p1,accvar4p4,accvar3p4,var4p1,accvar4p5,accvar3p5,(mpfr_ptr) 0);
	mpfr_set_si(maxD1,0,GMP_RNDN);
	mpfr_set_si(maxD2,0,GMP_RNDN);
	
	mpfr_t cevar3,cevar4,cevar5;
	mpfr_inits2(24,cevar3,cevar4,cevar5,(mpfr_ptr)0);
	
	mpfr_t avar3_d2 , avar4_d2;
	mpfr_inits2(24,avar3_d2, avar4_d2,(mpfr_ptr)0);
	
	mpfr_t mem,mem4;
	mpfr_init2(mem,24);
	mpfr_init2(mem4,24);
	
	ostringstream sv3,sv4;
	int verv3,verv4;
	//int mmcc=0;
	
	int cc1=2;
	for(int a1=0;a1<limit;a1++)
	{		

		if(a1%cc1==0)
		{
		cout<<(a1+1)*100/limit<<"%"<<endl;
		if(cc1<32)
		cc1*=2;
		}
	

		for(int a2=0;a2<limit;a2++)
		if(abs(a1-a2)>1)
		{
			
			
			temmpz = functionX(a1);
			if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
			}
			mpfr_set_z(sigX0,temmpz.get_mpz_t(),GMP_RNDN);
			if(a1%64!=63)
			{
				temmpz = functionX(a1+1);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigX1,temmpz.get_mpz_t(),GMP_RNDN);
			}
			else
			{
				temmpz = functionX(a1-63);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigX1,temmpz.get_mpz_t(),GMP_RNDN);
			}
			temmpz = functionX(a2);
			if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
			}
			mpfr_set_z(sigX2,temmpz.get_mpz_t(),GMP_RNDN);
			if(a2%64!=63)
			{
				temmpz = functionX(a2+1);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigX3,temmpz.get_mpz_t(),GMP_RNDN);
			}
			else
			{
				temmpz = functionX(a2-63);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigX3,temmpz.get_mpz_t(),GMP_RNDN);
			}
			
			mpfr_mul(sigX0,sigX0,cst1,GMP_RNDN);
			mpfr_mul(sigX1,sigX1,cst1,GMP_RNDN);
			mpfr_mul(sigX2,sigX2,cst1,GMP_RNDN);
			mpfr_mul(sigX3,sigX3,cst1,GMP_RNDN);
			
			
			mpfr_sub(segX1mX0,sigX1,sigX0,GMP_RNDN);
			mpfr_sub(segX0mX1,sigX0,sigX1,GMP_RNDN);
			mpfr_sub(segX3mX2,sigX3,sigX2,GMP_RNDN);
			mpfr_sub(segX3mX0,sigX3,sigX0,GMP_RNDN);
			mpfr_sub(segX0mX2,sigX0,sigX2,GMP_RNDN);
			
			
			
			temmpz = functionY(a1);
			if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
			}
			mpfr_set_z(sigY0,temmpz.get_mpz_t(),GMP_RNDN);
			if(a1%64!=63)
			{
				temmpz = functionY(a1+1);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigY1,temmpz.get_mpz_t(),GMP_RNDN);
			}
			else
			{
				temmpz = functionY(a1-63);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigY1,temmpz.get_mpz_t(),GMP_RNDN);
			}
			temmpz = functionY(a2);
			if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
			mpfr_set_z(sigY2,temmpz.get_mpz_t(),GMP_RNDN);
			if(a2%64!=63)
			{
				temmpz = functionY(a2+1);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigY3,temmpz.get_mpz_t(),GMP_RNDN);
			}
			else
			{
				temmpz = functionY(a2-63);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigY3,temmpz.get_mpz_t(),GMP_RNDN);
			}
			
			mpfr_mul(sigY0,sigY0,cst1,GMP_RNDN);
			mpfr_mul(sigY1,sigY1,cst1,GMP_RNDN);
			mpfr_mul(sigY2,sigY2,cst1,GMP_RNDN);
			mpfr_mul(sigY3,sigY3,cst1,GMP_RNDN);
			
			mpfr_sub(segY1mY0,sigY1,sigY0,GMP_RNDN);
			mpfr_sub(segY0mY1,sigY0,sigY1,GMP_RNDN);
			mpfr_sub(segY3mY2,sigY3,sigY2,GMP_RNDN);
			mpfr_sub(segY3mY0,sigY3,sigY0,GMP_RNDN);
			mpfr_sub(segY0mY2,sigY0,sigY2,GMP_RNDN);
			
			
			temmpz = functionZ(a1);
			if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
			}
			mpfr_set_z(sigZ0,temmpz.get_mpz_t(),GMP_RNDN);
			if(a1%64!=63)
			{
				temmpz = functionZ(a1+1);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigZ1,temmpz.get_mpz_t(),GMP_RNDN);
			}
			else
			{
				temmpz = functionZ(a1-63);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigZ1,temmpz.get_mpz_t(),GMP_RNDN);
			}
			temmpz = functionZ(a2);
			if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
			}
			mpfr_set_z(sigZ2,temmpz.get_mpz_t(),GMP_RNDN);
			if(a2%64!=63)
			{
				temmpz = functionZ(a2+1);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigZ3,temmpz.get_mpz_t(),GMP_RNDN);
			}
			else
			{
				temmpz = functionZ(a2-63);
				if (temmpz  > tmpCMP){ //negative number 
				temmpz  = temmpz  - tmpSUB;
				}
				mpfr_set_z(sigZ3,temmpz.get_mpz_t(),GMP_RNDN);
			}
			
			mpfr_mul(sigZ0,sigZ0,cst1,GMP_RNDN);
			mpfr_mul(sigZ1,sigZ1,cst1,GMP_RNDN);
			mpfr_mul(sigZ2,sigZ2,cst1,GMP_RNDN);
			mpfr_mul(sigZ3,sigZ3,cst1,GMP_RNDN);
			
			mpfr_sub(segZ1mZ0,sigZ1,sigZ0,GMP_RNDN);
			mpfr_sub(segZ0mZ1,sigZ0,sigZ1,GMP_RNDN);
			mpfr_sub(segZ3mZ2,sigZ3,sigZ2,GMP_RNDN);
			mpfr_sub(segZ3mZ0,sigZ3,sigZ0,GMP_RNDN);
			mpfr_sub(segZ0mZ2,sigZ0,sigZ2,GMP_RNDN);
			
			//var1
			
			mpfr_mul(temp1,segX1mX0,segX3mX2,GMP_RNDN);
			mpfr_mul(temp2,segY1mY0,segY3mY2,GMP_RNDN);
			mpfr_mul(temp3,segZ1mZ0,segZ3mZ2,GMP_RNDN);
			
			mpfr_add(temp4,temp1,temp2,GMP_RNDN);
			mpfr_add(var1,temp4,temp3,GMP_RNDN);
			
			mpfr_abs(temp1_23,var1,GMP_RNDN);
			if(mpfr_cmp(temp1_23,var1m)>0)
				mpfr_set(var1m,temp1_23,GMP_RNDN);
			texp=mpfr_get_exp(var1);
			if(texp<minV1&&(mpfr_sgn(var1)!=0)&&(mpfr_number_p(var1)!=0))
			{minV1=texp;
				mpfr_set(miV1,var1,GMP_RNDN);
			}
			if(texp>maxV1&&(mpfr_sgn(var1)!=0)&&(mpfr_number_p(var1)!=0))
			{maxV1=texp;
			mpfr_set(maV1,var1,GMP_RNDN);
				}
			
			
			//var2
			
			mpfr_sqr(temp1,segX3mX2,GMP_RNDN);
			mpfr_sqr(temp2,segY3mY2,GMP_RNDN);
			mpfr_sqr(temp3,segZ3mZ2,GMP_RNDN);
			
			mpfr_add(temp4,temp1,temp2,GMP_RNDN);
			mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
			mpfr_sqrt(var2,temp1_23,GMP_RNDN);
			
			mpfr_abs(temp1_23,var2,GMP_RNDN);
			if(mpfr_cmp(temp1_23,var2m)>0)
				mpfr_set(var2m,temp1_23,GMP_RNDN);
			
			texp=mpfr_get_exp(var2);
			if(texp<minV2&&(mpfr_sgn(var2)!=0)&&(mpfr_number_p(var2)!=0))
			{minV2=texp;
				mpfr_set(miV2,var2,GMP_RNDN);
			}
			if(texp>maxV2&&(mpfr_sgn(var2)!=0)&&(mpfr_number_p(var2)!=0))
			{maxV2=texp;
				mpfr_set(maV2,var2,GMP_RNDN);
			}
			
			mpfr_set_si(accvar3,0,GMP_RNDN);
			mpfr_set_si(accvar3_2,0,GMP_RNDN);
			mpfr_set_si(accvar3_1,0,GMP_RNDN);
			mpfr_set_si(accvar4,0,GMP_RNDN);
			mpfr_set_si(accvar4p4,0,GMP_RNDN);
			mpfr_set_si(accvar3p4,0,GMP_RNDN);
			mpfr_set_si(accvar3p5,0,GMP_RNDN);
			mpfr_set_si(accvar4p5,0,GMP_RNDN);
			mpfr_set_si(accvar4_2,0,GMP_RNDN);
			mpfr_set_si(accvar4_1,0,GMP_RNDN);
			mpfr_set_si(accvar5,0,GMP_RNDN);
			mpfr_set_si(accvar5_1,0,GMP_RNDN);
			mpfr_set_si(accvar5_2,0,GMP_RNDN);
			
			mpfr_set_si(accvar3_d2,0,GMP_RNDN);
			mpfr_set_si(accvar4_d2,0,GMP_RNDN);
			mpfr_set_si(var3p1,0,GMP_RNDN);
			mpfr_set_si(var4p1,0,GMP_RNDN);
			mpfr_set_si(avar3_d2,0,GMP_RNDN);
			mpfr_set_si(avar4_d2,0,GMP_RNDN);
			
					
			sv3.str("");
			sv4.str("");
			verv3=0;
			verv4=0;
			
			for(int t=0;t<8;t++)
			{
				mpfr_set_ui(sigt,t,GMP_RNDN);
				mpfr_mul(sigt,sigt,cst2,GMP_RNDN);
				//mpfr_printf("%Rf\n",sigt);
				
				
				
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
				mpfr_add(accvar3,accvar3,temp1_23,GMP_RNDN);
				
				if(t==0||t==4)
				{
					mpfr_add(accvar3_2,accvar3_2,temp1_23,GMP_RNDN);
					mpfr_add(accvar3_1,accvar3_1,temp1_23,GMP_RNDN);
					
				}
				if(t==2||t==6)
				{
					mpfr_add(accvar3_2,accvar3_2,temp1_23,GMP_RNDN);
					mpfr_add(accvar3_d2,accvar3_d2,temp1_23,GMP_RNDN);
				}
				
				if(t==4||t==3)
				{
					mpfr_add(avar3_d2,avar3_d2,temp1_23,GMP_RNDN);
				}
				
				if(t==4)
				{
				mpfr_set(var3p1,temp1_23,GMP_RNDN);	
				}
				
				
				
				
				//mpfr_printf("V3:%Rf ",temp1_23);
				
				if(t>0)
				{
					mpfr_sub(temp2_23,temp1_23,mem,GMP_RNDN);
					//mpfr_printf("D3:%Rf ",temp2_23);
					if(mpfr_sgn(temp2_23)<0)
					{
						sv3<<"-";						
					}
					else
					{
						sv3<<"+";
						verv3++;
					}
				}
				mpfr_set(mem,temp1_23,GMP_RNDN);
				
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
				mpfr_add(accvar4,accvar4,temp1_23,GMP_RNDN);
				
				//mpfr_printf("V4:%Rf ",temp1_23);
				
				if(t==0||t==4)
				{
					mpfr_add(accvar4_2,accvar4_2,temp1_23,GMP_RNDN);
					mpfr_add(accvar4_1,accvar4_1,temp1_23,GMP_RNDN);
					
				}
				if(t==2||t==6)
				{
					mpfr_add(accvar4_2,accvar4_2,temp1_23,GMP_RNDN);
					mpfr_add(accvar4_d2,accvar4_d2,temp1_23,GMP_RNDN);
				}
				
				if(t==4)
				{
				mpfr_set(var4p1,temp1_23,GMP_RNDN);	
				}
				
				if(t==4||t==3)
				{
					mpfr_add(avar4_d2,avar4_d2,temp1_23,GMP_RNDN);
				}
				
				if(t>0)
				{
					mpfr_sub(temp2_23,temp1_23,mem4,GMP_RNDN);
					//mpfr_printf("D4:%Rf ",temp2_23);
					if(mpfr_sgn(temp2_23)<0)
					{
						sv4<<"-";
					}
					else
					{
						sv4<<"+";
						verv4++;
						
					}
				}
				mpfr_set(mem4,temp1_23,GMP_RNDN);
				
				//var5
				
				mpfr_mul(temp1,temp5,segX3mX2,GMP_RNDN);
				mpfr_mul(temp2,temp6,segY3mY2,GMP_RNDN);
				mpfr_mul(temp3,temp7,segZ3mZ2,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_add(accvar5,accvar5,temp1_23,GMP_RNDN);
				
				if(t==0||t==4)
				{
					mpfr_add(accvar5_2,accvar5_2,temp1_23,GMP_RNDN);
					mpfr_add(accvar5_1,accvar5_1,temp1_23,GMP_RNDN);
					
				}
				if(t==2||t==6)
				{
					mpfr_add(accvar5_2,accvar5_2,temp1_23,GMP_RNDN);
				}
				
				
			}
			
			mpfr_set_si(cevar3,0,GMP_RNDN);
			mpfr_set_si(cevar4,0,GMP_RNDN);
			mpfr_set_si(cevar5,0,GMP_RNDN);
			
			for(int t=1;t<16;t=t+2)
			{
				
				mpfr_set_ui(sigt4,t,GMP_RNDN);
				mpfr_mul(sigt4,sigt4,cst5,GMP_RNDN);
				//mpfr_printf("%Rb\n",sigt4);
				
				
				
				
				//var3
				mpfr_mul(temp1,segX0mX1,sigt4,GMP_RNDN);
				mpfr_add(temp1,temp1,segX3mX0,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY0mY1,sigt4,GMP_RNDN);
				mpfr_add(temp2,temp2,segY3mY0,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ0mZ1,sigt4,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ3mZ0,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_add(cevar3,cevar3,temp1_23,GMP_RNDN);
				
				
				
				
				
				
				//var4
				mpfr_mul(temp1,segX1mX0,sigt4,GMP_RNDN);
				mpfr_add(temp1,temp1,segX0mX2,GMP_RNDN);
				mpfr_set(temp5,temp1,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY1mY0,sigt4,GMP_RNDN);
				mpfr_add(temp2,temp2,segY0mY2,GMP_RNDN);
				mpfr_set(temp6,temp2,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ1mZ0,sigt4,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ0mZ2,GMP_RNDN);
				mpfr_set(temp7,temp3,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_add(cevar4,cevar4,temp1_23,GMP_RNDN);
				
				
				//var5
				
				mpfr_mul(temp1,temp5,segX3mX2,GMP_RNDN);
				mpfr_mul(temp2,temp6,segY3mY2,GMP_RNDN);
				mpfr_mul(temp3,temp7,segZ3mZ2,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_add(cevar5,cevar5,temp1_23,GMP_RNDN);
				
			}
			mpfr_mul(cevar3,cevar3,cst2,GMP_RNDN);
			mpfr_mul(cevar4,cevar4,cst2,GMP_RNDN);
			mpfr_mul(cevar5,cevar5,cst2,GMP_RNDN);
			
			
			//var3 reload with t on 4 bits
				mpfr_mul(temp1,segX0mX1,sigt2,GMP_RNDN);
				mpfr_add(temp1,temp1,segX3mX0,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY0mY1,sigt2,GMP_RNDN);
				mpfr_add(temp2,temp2,segY3mY0,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ0mZ1,sigt2,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ3mZ0,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_set(accvar3p4,temp1_23,GMP_RNDN);
			
			//var4 reloaded with t on 4 bits
				mpfr_mul(temp1,segX1mX0,sigt2,GMP_RNDN);
				mpfr_add(temp1,temp1,segX0mX2,GMP_RNDN);
				mpfr_set(temp5,temp1,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY1mY0,sigt2,GMP_RNDN);
				mpfr_add(temp2,temp2,segY0mY2,GMP_RNDN);
				mpfr_set(temp6,temp2,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ1mZ0,sigt2,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ0mZ2,GMP_RNDN);
				mpfr_set(temp7,temp3,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_set(accvar4p4,temp1_23,GMP_RNDN);
				
				
				
			//var3 reload with t on 5 bits
				mpfr_mul(temp1,segX0mX1,sigt3,GMP_RNDN);
				mpfr_add(temp1,temp1,segX3mX0,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY0mY1,sigt3,GMP_RNDN);
				mpfr_add(temp2,temp2,segY3mY0,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ0mZ1,sigt3,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ3mZ0,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_set(accvar3p5,temp1_23,GMP_RNDN);
			
			//var4 reloaded with t on 5 bits
				mpfr_mul(temp1,segX1mX0,sigt3,GMP_RNDN);
				mpfr_add(temp1,temp1,segX0mX2,GMP_RNDN);
				mpfr_set(temp5,temp1,GMP_RNDN);
				mpfr_sqr(temp1,temp1,GMP_RNDN);
				mpfr_mul(temp2,segY1mY0,sigt3,GMP_RNDN);
				mpfr_add(temp2,temp2,segY0mY2,GMP_RNDN);
				mpfr_set(temp6,temp2,GMP_RNDN);
				mpfr_sqr(temp2,temp2,GMP_RNDN);
				mpfr_mul(temp3,segZ1mZ0,sigt3,GMP_RNDN);
				mpfr_add(temp3,temp3,segZ0mZ2,GMP_RNDN);
				mpfr_set(temp7,temp3,GMP_RNDN);
				mpfr_sqr(temp3,temp3,GMP_RNDN);
				
				mpfr_add(temp4,temp1,temp2,GMP_RNDN);
				mpfr_add(temp1_23,temp4,temp3,GMP_RNDN);
				mpfr_sqrt(temp1_23,temp1_23,GMP_RNDN);
				mpfr_set(accvar4p5,temp1_23,GMP_RNDN);
			
			
			
			//~ cout<<sv3.str();
			//~ if(verv3==7||verv3==0)
				//~ cout<<" OK         ";
			//~ else
				//~ cout<<" NOT        ";
			//~ cout<<sv4.str();
			//~ if(verv4==7||verv4==0)
				//~ cout<<" OK         ";
			//~ else
				//~ cout<<" NOT         ";
			 //~ cout<<endl;
			
			mpfr_set(var3,accvar3,GMP_RNDN);
			mpfr_mul(var3,var3,cst2,GMP_RNDN);
			mpfr_mul(accvar3_2,accvar3_2,cst3,GMP_RNDN);
			mpfr_mul(accvar3_d2,accvar3_d2,cst4,GMP_RNDN);
			mpfr_mul(accvar3_1,accvar3_1,cst4,GMP_RNDN);
			mpfr_mul(avar3_d2,avar3_d2,cst4,GMP_RNDN);
			
			mpfr_abs(temp3_23,var3,GMP_RNDN);
			mpfr_abs(temp2_23,accvar3_2,GMP_RNDN);
			mpfr_sub(temp1_23,temp3_23,temp2_23,GMP_RNDN);
			mpfr_abs(temp1_23,temp1_23,GMP_RNDN);
			if(mpfr_cmp(temp1_23,maxD2)>0)
				mpfr_set(maxD2,temp1_23,GMP_RNDN);
			mpfr_abs(temp2_23,accvar3_1,GMP_RNDN);
			mpfr_sub(temp1_23,temp3_23,temp2_23,GMP_RNDN);
			mpfr_abs(temp1_23,temp1_23,GMP_RNDN);
			if(mpfr_cmp(temp1_23,maxD1)>0)
				mpfr_set(maxD1,temp1_23,GMP_RNDN);
			
			
			
			texp=mpfr_get_exp(var3);
			if(texp<minV3&&(mpfr_sgn(var3)!=0)&&(mpfr_number_p(var3)!=0))
			{minV3=texp;
				mpfr_set(miV3,var3,GMP_RNDN);
			}
			if(texp>maxV3&&(mpfr_sgn(var3)!=0)&&(mpfr_number_p(var3)!=0))
			{maxV3=texp;
				mpfr_set(maV3,var3,GMP_RNDN);
			}
			
			
			mpfr_set(var4,accvar4,GMP_RNDN);
			mpfr_mul(var4,var4,cst2,GMP_RNDN);
			mpfr_mul(accvar4_2,accvar4_2,cst3,GMP_RNDN);
			mpfr_mul(accvar4_d2,accvar4_d2,cst4,GMP_RNDN);
			mpfr_mul(accvar4_1,accvar4_1,cst4,GMP_RNDN);
			mpfr_mul(avar4_d2,avar4_d2,cst4,GMP_RNDN);
			
			texp=mpfr_get_exp(var4);
			if(texp<minV4&&(mpfr_sgn(var4)!=0)&&(mpfr_number_p(var4)!=0))
			{minV4=texp;
				//cout<<"pt min"<<endl;
				//mpfr_printf("%Rf %d %d ",var4,a1,a2);
				mpfr_set(miV4,var4,GMP_RNDN);
				//mpfr_printf(" %Rf\n",miV4);
			}
			if(texp>maxV4&&(mpfr_sgn(var4)!=0)&&(mpfr_number_p(var4)!=0))
			{maxV4=texp;
				//cout<<"pt max"<<endl;
				//mpfr_printf("%Rf %d %d ",var4,a1,a2);
				mpfr_set(maV4,var4,GMP_RNDN);
				//mpfr_printf(" %Rf\n",maV4);
			}
			
			
			mpfr_set(var5,accvar5,GMP_RNDN);
			mpfr_mul(var5,var5,cst2,GMP_RNDN);
			mpfr_mul(accvar5_2,accvar5_2,cst3,GMP_RNDN);
			mpfr_mul(accvar5_1,accvar5_1,cst4,GMP_RNDN);
			
			texp=mpfr_get_exp(var5);
			if(texp<minV5&&(mpfr_sgn(var5)!=0)&&(mpfr_number_p(var5)!=0))
			{minV5=texp;
				mpfr_set(miV5,var5,GMP_RNDN);
			}
			if(texp>maxV5&&(mpfr_sgn(var5)!=0)&&(mpfr_number_p(var5)!=0))
			{maxV5=texp;
				mpfr_set(maV5,var5,GMP_RNDN);
			}
			
			if(mpfr_cmp(var3,var3m)>0)
				mpfr_set(var3m,var3,GMP_RNDN);
			if(mpfr_cmp(var4,var4m)>0)
				mpfr_set(var4m,var4,GMP_RNDN);
			if(mpfr_cmp(var5,var5m)>0)
				mpfr_set(var5m,var5,GMP_RNDN);
			
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
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(acc,acc,temp1_23,GMP_RNDN);
			
			
			
			mpfr_div(temp3_23,accvar5_2,var2,GMP_RNDN);
			mpfr_add(temp1_23,accvar3_2,var2,GMP_RNDN);
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			mpfr_sub(temp2_23,accvar4_2,temp3_23,GMP_RNDN);
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(acc2,acc2,temp1_23,GMP_RNDN);
			
			
			mpfr_div(temp3_23,accvar5_1,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,accvar3_1,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp2_23,accvar4_1,temp3_23,GMP_RNDN);
			
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(acc1,acc1,temp1_23,GMP_RNDN);
			
			//de aici se modifica
			
			mpfr_div(temp3_23,cevar5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,accvar3p5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
                        //mpfr_printf("%Rf  ",temp1_23);		
			
			
			mpfr_sub(temp2_23,accvar4p5,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(acc31,acc31,temp1_23,GMP_RNDN);
			
			
			
			
			
			
			
			
			
			mpfr_div(temp3_23,cevar5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,accvar3_d2,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
                        //mpfr_printf("%Rf  ",temp1_23);		
			
			
			mpfr_sub(temp2_23,accvar4_d2,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(acc4,acc4,temp1_23,GMP_RNDN);
			
			
			
						
			mpfr_div(temp3_23,cevar5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,cevar3,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
                        //mpfr_printf("%Rf  ",temp1_23);		
			
			
			mpfr_sub(temp2_23,cevar4,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(acc21,acc21,temp1_23,GMP_RNDN);
			
			
						
			
			mpfr_div(temp3_23,cevar5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,var3p1,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
                        //mpfr_printf("%Rf  ",temp1_23);		
			
			
			mpfr_sub(temp2_23,var4p1,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(accn1,accn1,temp1_23,GMP_RNDN);
			
			
			
			
			mpfr_div(temp3_23,var5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,accvar3p4,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
                        //mpfr_printf("%Rf  ",temp1_23);		
			
			
			mpfr_sub(temp2_23,accvar4p4,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(accn3,accn3,temp1_23,GMP_RNDN);
			
						
			
			
			mpfr_div(temp3_23,var5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_add(temp1_23,accvar3p5,var2,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			
			mpfr_sub(temp1_23,temp1_23,temp3_23,GMP_RNDN);
                        //mpfr_printf("%Rf  ",temp1_23);		
			
			
			mpfr_sub(temp2_23,accvar4p5,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp2_23);
			
			mpfr_div(temp1_23,temp1_23,temp2_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp1_23);
			
			mpfr_log(temp3_23,temp1_23,GMP_RNDN);
			//mpfr_printf("%Rf  ",temp3_23);
			mpfr_div(temp1_23,var1,var2,GMP_RNDN);
			mpfr_mul(temp1_23,temp1_23,temp3_23,GMP_RNDN);
			//mpfr_printf("%Rf\n",temp1_23);
			if(mpfr_number_p(temp1_23)!=0)
			mpfr_add(accn4,accn4,temp1_23,GMP_RNDN);
			
			
			//~ mpfr_printf("var3:%Rf %Rf %Rf  %Rf %Rf %Rf\nvar4: %Rf %Rf %Rf %Rf %Rf %Rf\n",var3,cevar3,var3p1,avar3_d2,accvar3p4,accvar3p5,var4,cevar4,var4p1,avar4_d2,accvar4p4,accvar4p5);
			//~ mpfr_printf("var5:%Rf %Rf\n",var5,cevar5);
			//~ mpfr_printf("var1:%Rf var2:%Rf\n",var1,var2);
			
			
			
		}
	}
	
	mpfr_printf("maxABS var1:=+/-%Rf\n",var1m);
	mpfr_printf("maxABS var2:=+/-%Rf\n",var2m);
	mpfr_printf("maxABS var3:=+/-%Rf\n",var3m);
	mpfr_printf("maxABS var4:=+/-%Rf\n",var4m);
	mpfr_printf("maxABS var5:=+/-%Rf\n",var5m);
	mpfr_printf("min_exp_Var1:=%d for value:=%Rf ; max_exp_Var1:=%d for value:=%Rf\n",((int)minV1),miV1,((int)maxV1),maV1);
	mpfr_printf("min_exp_Var2:=%d for value:=%Rf ; max_exp_Var2:=%d for value:=%Rf\n",((int)minV2),miV2,((int)maxV2),maV2);
	mpfr_printf("min_exp_Var3:=%d for value:=%Rf ; max_exp_Var3:=%d for value:=%Rf\n",((int)minV3),miV3,((int)maxV3),maV3);
	mpfr_printf("min_exp_Var4:=%d for value:=%Rf ; max_exp_Var4:=%d for value:=%Rf\n",((int)minV4),miV4,((int)maxV4),maV4);
	mpfr_printf("min_exp_Var5:=%d for value:=%Rf ; max_exp_Var5:=%d for value:=%Rf\n",((int)minV5),miV5,((int)maxV5),maV5);
	
	mpfr_printf("max diff 4 var 3 with t on 3 and t on 2 is :=%Rf\n",maxD2);
	mpfr_printf("max diff 4 var 3 with t on 3 and t on 1 is :=%Rf\n",maxD1);
	
	mpfr_printf("Final value with t on 3 bits is:=%Rf\n",acc);
	mpfr_printf("Final value with t on 2 bits is:=%Rf\n",acc2);
	mpfr_printf("Final value with t on 1 bits is:=%Rf\n",acc1);
	mpfr_printf("Final value with t in 0.0111..  for variables 3,4 and t on 4(8) bits for var5 is:=%Rf\n",acc31);
	mpfr_printf("Final value with t on 2(0.25 and 0.75) for variables 3,4 and t on 4(8) bits for var5 is:=%Rf\n",acc4);
	mpfr_printf("Final value with t on 4(8) for variables 3,4 and t on 4(8) bits for var5 is:=%Rf\n",acc21);
	mpfr_printf("Final value with t on 1() for variables 3,4 and t on 4(8) bits for var5 is:=%Rf\n",accn1);
	//mpfr_printf("Final value with t on 1(com) for variables 3,4 and t on 3 bits for var5 is:=%Rf\n",accn2);
	mpfr_printf("Final value with t on 4(1 value) for variables 3,4 and t on 3 bits for var5 is:=%Rf\n",accn3);
	mpfr_printf("Final value with t on 5(0.0111.. ) for variables 3,4 and t on 3 bits for var5 is:=%Rf\n",accn4);
	
	
	mpfr_clears(acc, acc2, acc1, var1, var2, var3, var4, var5, temp1_23, temp2_23, temp3_23, accvar4_1, accvar3_1, accvar5_1, accvar4_2, accvar3_2, accvar5_2,temp1,temp2,temp3,temp4,temp5,temp6,temp7,(mpfr_ptr) 0);
	mpfr_clears(cst1,cst2,cst3,cst4,accvar5,accvar3,accvar4,(mpfr_ptr)0);
	mpfr_clears(var1m,var2m,var3m,var4m,var5m,miV1,maV1,miV2,maV2,miV3,maV3,miV4,maV4,miV5,maV5,(mpfr_ptr)0);
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



mpz_class functionY(int x)
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
	mpfr_clear(yv);
	mpfr_clear(temp);
	mpfr_clear(tauqj);
		mpfr_free_cache();
	return r;
}



mpz_class functionZ(int x)
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
	mpfr_clear(zv);
	mpfr_clear(temp);
	mpfr_clear(tauqj);
		mpfr_free_cache();
	return r;
}


	
mpz_class  functionX(int x)
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

	mpfr_t xv;
	mpfr_init2(xv,1000);
	
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
				
				
		mpfr_mul(xv, xv, temp1, GMP_RNDN);	
		
		mpz_t xvz;
		mpz_init2(xvz,1000);
		mpfr_get_z(xvz,xv,GMP_RNDN);
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
	mpfr_clear(xv);
	mpfr_free_cache();
	mpfr_clear(temp);
	mpfr_clear(tauqj);
	
	mpfr_free_cache();
	return r;
}
