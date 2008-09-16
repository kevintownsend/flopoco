/*
 * A floating-point logarithm for FloPoCo
 * 
 * Author : Florent de Dinechin
 *
 * For a description of the algorithm, see the Arith17 paper. 
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
#include <iostream>
#include <sstream>
#include <fstream>
//#include <mpfrxx.h>
#include <math.h>
#include "math_lib.hpp"
#include "LogArch.hpp"
#include "../LZOC.hpp"

using namespace std;

LogArch::LogArch(int wE, int wF) : 
	wE(wE), wF(wF)
{
	int i;

	gLog = 3; // TODO: almost randomly chosen
	target_prec = wF+((wF+1)>>1) +gLog;

	// First compute the precision of each iteration 

	// Stage 0
	p[0] = 0;
	a[0] = 5; 
	s[0] = wF+2;  
	sfullZ[0] = wF+2;

	p[1] = a[0]-1; 
	sbt[1] = wF+2 ;
	s[1] = wF+2;
	t[1] = 0;
	sfullZ[1] = wF+7;

	// Following stages -- this is true starting from stage 1, although
	// stage 1 needs a specific inverter table
	i=1;
	while(2*p[i] < wF+2){ // ensures 2*p[stages+1]>= wF+2, enough for faithful rounding
		a[i] = 4;
		p[i+1] = p[i] + a[i] - 1;

		// size before truncation
		sbt[i+1] = s[i] +  p[i] + 2;

		if(p[i+1]+sbt[i+1] <= target_prec) 
			{ // No truncation at all
				psize[i] = s[i];
				s[i+1] = sbt[i+1];
				t[i+1] = 0;
			}
		else
			{ // Truncate everybody to targetprec : 
				// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>eps[i]
				// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
				// We need only this many bits of Z[i] to compute it.
				psize[i] = target_prec - 2*p[i];
				if (psize[i]>s[i]) // in the first iterations
				  psize[i]=s[i];
				s[i+1] = target_prec - p[i+1];
				t[i+1] = sbt[i+1] - s[i+1];
			}

 
 
		sfullZ[i+1] =  sfullZ[i] + a[i] + p[i] + 1;
		i++;
	}  


	// Deduce the number of stages
	stages = i-1;

	// MSB of squarer input is p[stages+1];
	// LSB will be target_prec
	// size will be target_prec - p[stages+1]  

	cout<<"LogArch::LogArch: needs "<<stages<<" stages"<<endl;

	// Now allocate the various table objects

	it0 = new FirstInvTable(a[0], a[0]+1);
	lt0 = new FirstLogTable(a[0], target_prec, it0);
	it1 = new SecondInvTable(a[1], p[1]);
	for(i=1; i<=stages; i++) {
		lt[i] = new OtherLogTable(a[i], target_prec - p[i], p[i], i); 
	 }
	int computedG = (int) ceil(log(3*(stages+1))/log(2));

	/* Use FloPoCo operators */
	lzoc = new LZOC(target, wF, intlog2(wF));
	lzoc->set_combinatorial();
} 


LogArch::~LogArch()
{
	delete lzoc;
}

void LogArch::describe(){
	int i;
	cout<<endl;
	for(i=0; i<=stages; i++)
		cout<< i <<" : a=" << a[i] 
				<< "   p="<< p[i]
			//	<< "   sfullZ=" << sfullZ[i]
				<< "   sbt=" << sbt[i]
 	<< "   t=" << t[i]
				<< "   s=" << s[i] 
 	<< "   psize=" << psize[i]
				<< endl;
	cout<<"stages=" <<stages <<"  pfinal="<<p[stages+1]<<"  target_prec="<<  target_prec <<endl;
	cout<<"########### Approx size in LUTs: "<< size_in_LUTs() << endl;
}


// void LogArch::check(){
//   // Check that first stage is OK
//   int tp1 = it0 -> check_accuracy(wF);
//   cout << "p1 computed is "<< tp1 <<endl;
// };


// An auxiliary function because mpz_class doesn't support long double. Strange.
// To remove as soon as it does.
mpz_class mpz_class_(uint64_t x) {
	uint64_t xh,xl;
	mpz_class r;
	xh = x >> 32;
	xl= x - (xh <<32);
	r = (mpz_class((uint32_t)xh) << 32) + mpz_class((uint32_t)xl);
	return r;
}



// 
void showMP(ostream& o, string s, mpz_class x, int p, int size) {
	int j;
	o << s;
	for(j=0; j<p+1; j++) 
		o<<" ";
	printBinNumGMP(o, x, size);
	o<<endl;
}


void LogArch::showFP(ostream& o, double x) {
	uint64_t E, E0,Min;
	mpfr_t xx;
	int i,s; 
	double f, xxx;

	f=frexp(x, &i);   // get mantissa and exponent
	f=f*2; E=i-1;     // frexp returns mantissa between 0.5 and 1, I prefer [1..2]

	//sign
	if (x>=0) s=0; 
	else{
		s=1; 
		f=-f;
	}
	
	E0=(1ull<<(wE-1)) -1; 

	// Throw the input into an mpfr number just right
	mpfr_init2(xx, wF+1);
	mpfr_set_d(xx, f, GMP_RNDN);
	xxx=mpfr_get_d(xx, GMP_RNDN);
	Min=(uint64_t)(xxx*(1ull<<wF)); 

	o << "01"<<s; 
	printBinNum(o, E+E0, wE); 
	printBinNum(o, Min - (1ull<<wF), wF); 
}








// This function is more and more out of sync with the VHDL output.

double LogArch::simulate(double x) 
{
	uint64_t E, E0, Min;
	mpz_class Z[42], A[42], B[42], P[42], fullZ[42], fullInvA[42];
	mpz_class Zsquare, Log1p, Log, Log2, ELog2;
	mpfr_t mantissa_m,x_m,logM_m, LogMCalc_m, LogXCalc_m, LogX_m, two, Log2_m;
	mpz_t zlog2;
	double f, mantissa_d, LogMCalc_d,logM_d, LogX_d;
	int i,j; 
	bool FirstBit;

	mpfr_init2(x_m,wF+1);
	mpfr_set_d(x_m, x, GMP_RNDN);

	mpfr_init2(LogX_m,wF+1);
	mpfr_log(LogX_m, x_m, GMP_RNDN);
	LogX_d= mpfr_get_d(LogX_m,GMP_RNDN);

	cout << "        X: "; showFP(cout,x); cout<<endl;
	cout << "CR result: "; showFP(cout, LogX_d); cout<<endl;

	mpfr_init2(LogMCalc_m,200);

	mpfr_init2(LogXCalc_m,200);

	mpfr_init2(two,2);
	mpfr_set_d(two, 2.0, GMP_RNDN);

	f=frexp(x, &i);   // get mantissa and exponent
	f=f*2; E=i-1;     // frexp returns mantissa between 0.5 and 1, I prefer [1..2]

	E0=(1ull<<(wE-1)) -1; 

	// Throw the input into an mpfr number just right
	mpfr_init2(mantissa_m, wF+1);
	mpfr_set_d(mantissa_m, f, GMP_RNDN);

	mantissa_d=mpfr_get_d(mantissa_m, GMP_RNDN);
	Min=(uint64_t)(mantissa_d*(1ull<<wF)); 

	//************************** First stage***************
	// Take the a[0]+1 first bits of the mantissa, but remove the leading implicit 1

	cout <<endl<<"Input exponent: "<<i-1
			 <<"  exponent field: "; 
	printBinNum(cout, E+E0, wE); cout<<endl ;  
	
	cout <<"Input fraction : "<<f<<"   ";
	printBinNum(cout, Min, wF+1); cout<<endl;

	A[0]= mpz_class_(  (Min>>(wF+1 -(a[0]+1)))  - (1ull<<a[0]));

	fullInvA[0]=it0->function(A[0].get_ui());

	FirstBit= (1==(A[0]>>(a[0]-1)));

	// Here we throw the mantissa into an int of size wF+1
	if(FirstBit) {
		fullZ[0] = mpz_class_(Min);
		mantissa_d=mantissa_d/2;
		E=E+1;
	}
	else{
		cout<<"    ... mantissa smaller than 1.5"<<endl;
		fullZ[0] = mpz_class_(Min <<1);
	}
 
	Z[0] = fullZ[0];
	

	int lzo = 0; 
	mpz_class digit; 
	if(FirstBit)  // Z0  between 0.75 and 1, count the ones
		digit = 1;
	else // Z0 between 1 and 1.5, count the zeroes
		digit=0;
	j=wF;
	while ( digit==((Z[0]>>j)&mpz_class(1)) && (j>=0)) {
			lzo++;
			j--;
	}
	cout<<endl<<"**************LZ="<<lzo<<endl;
	
	//Now we have two paths.

	if (lzo>p[stages+1]) {  // TODO ajouter first bit et test sur E
		cout << "*********** Close path"<<endl;
		showMP(cout, "Z0    ", Z[0], 0, 80);
		Log = Z[0] - (mpz_class(1)<<(s[0]-1));
		showMP(cout, "Log   ", Log, 0, 80);
		Log = Log << (lzo+1+gLog); 
		showMP(cout, "LogPS ", Log, 0, 80);
		Zsquare = (Log>>(p[stages+1])) * (Log>>(p[stages+1]));
		showMP(cout, "Zsq   ", Zsquare, 0, 80);
		Zsquare = Zsquare >> (wF+1+gLog - 2*p[stages+1] );  // truncated mult
		showMP(cout, "Zsq   ", Zsquare, 0, 80);
		Zsquare = Zsquare >> (lzo+2); // Align the LSBs and divide by 2
		showMP(cout, "Zsq   ", Zsquare, 0, 80);
		Log = Log - Zsquare;
		showMP(cout, "Log1p ", Log, 0, 80);

		// convert to double 
		mpfr_set_z(LogMCalc_m, Log.get_mpz_t(), GMP_RNDN);
		mpfr_shift_right(LogMCalc_m, lzo + wF + gLog + 2);
		
	}
	else {
		cout << "*********** Far path"<<endl;
		fullZ[1] = fullZ[0] * fullInvA[0];

		Z[1] = (fullZ[1] - (mpz_class(1)<<(sfullZ[1]-1)));
		// check no overflow TODO if bug
	
		//*************  Last stages : we truncate... 
		for(i=1; i<=stages; i++) {
			
			// The bits of A[i] are the  a[i] first bits of fullZ[i], 
			// smaller than 2^-p[i] if everything went well,
			// after removing the leading implicit 1
			A[i]= Z[i] >> (s[i]-a[i])  ;
			B[i] = Z[i] - (A[i]<<(s[i]-a[i]));
			// B has size (s[i]-a[i]), MSB weights -p[i]-a[i]  
			//                         LSB weights (before truncation) -p[i]-s[i] 

			
			// Full mantissa computation, for debugging
			if (i==1){
				// use the Table object
				fullInvA[i]= (mpz_class(1)<<(p[i]+a[i]+1)) - (it1->function(A[i].get_ui()));
			}
			else{
				if (A[i]==0)
				  fullInvA[i]=mpz_class(1)<<(p[i]+a[i]+1);
				else
				  fullInvA[i]= (mpz_class(1)<<(p[i]+a[i]+1)) - (A[i]<<1) + 1;
			}
			fullZ[i+1] = fullZ[i] * fullInvA[i]; // 
			


		// Actual mantissa computation
			P[i] = Z[i]*A[i];
			// P has size (s[i]+a[i]), MSB weights -2*p[i]     LSB weights -2p[i]-s[i]-a[i]
			
			Z[i+1] = (B[i]<<(p[i]+a[i]+1)) - (P[i]<<1);
			
			// Now add epsilon(1+Z[i])
			if (i==1){
				if(1==A[i]>>(a[i]-1)){ // MSB of A
				Z[i+1] +=  Z[i]<<1 ;
				Z[i+1] += ((mpz_class(1))<<(2*p[i]-a[i]+s[i]+1)) ; // 
				}
				else if (A[i]!=0) {
				  Z[i+1] += Z[i];
				  Z[i+1] += (mpz_class(1))<<(2*p[i]-a[i]+s[i]) ; // 
				}
			} 
			else{ //  if i>=2
				if (A[i]!=0){
				  Z[i+1] += Z[i];
				  Z[i+1] += (mpz_class(1))<<(s[i]+p[i]) ; // 
				}
			}
			

			// Truncation
			Z[i+1] = Z[i+1] >> t[i+1];
		}

		fullSizeSim=true;

		// showMP(cout,"       A0= ",A[0], 0, a[0]); 
		// showMP(cout,"fullInvA0= ", fullInvA[0], 0, it0->wOut ); 
		// showMP(cout,"    fullZ0 : " , fullZ[0], 0, sfullZ[0] ); 
		cout<<"   Z"<<0<<" :"; // un espace de moins pour le 1 implicite
		showMP(cout, "", Z[0], 0, s[0]); 
		if(fullSizeSim) {
			showMP(cout,"  fZ1 :" , fullZ[1], 0, sfullZ[1] );
		} 
		cout<<"   Z"<<1<<" : "; showMP(cout, "", Z[1], p[1], s[1]); 


		for(i=1; i<=stages; i++) {
			
			//       cout << "   **   A"<<i<<"= ";
			//       printBinNumGMP(cout, A[i], a[i]); cout<<endl;
			//       cout << "   **   B"<<i<<":     ";
			//       printBinNumGMP(cout, B[i],  (s[i]-a[i])    ); cout<<endl;
			//       cout<<"  fullInvA"<<i<<" = ";
			//       printBinNumGMP(cout, fullInvA[i],  a[i]+p[i]+2); cout<<endl;
			if(fullSizeSim) {
				cout <<"  fZ"<<i+1<<" : "; 
				if (sfullZ[i+1]>120){
				  printBinNumGMP(cout, fullZ[i+1]>>(sfullZ[i+1]-120) , 120); cout<<"..."<<endl;
				}
				else{
				  printBinNumGMP(cout, fullZ[i+1], sfullZ[i+1]); cout<<endl;
				}
			}
			//     cout <<"        Z"<<i+1<<" : ";
			//     for(j=0; j<p[i+1]+1; j++) cout<<" ";
			//     printBinNum(cout, Z[i+1],  (s[i+1] + t[i+1])    ); cout<< "  (before trunc)"<< endl;
			cout<<"   Z"<<i+1<<" : "; showMP(cout, "", Z[i+1], p[i+1], s[i+1]); 
			      
		}  


		// Now compute the square. Let p=p[stages+1] and s= s[stages+1]
		// Z starts at position p, and is of size s
		// Z^2/2 starts at position 2p+1  and and we need it up to position p+s 
		// Its size is 2s. Therefore we need to shift by (2p+2s+1 -(p+s)) = p+s+1
		// Its size after shifting is thus s-p-1
		
		Zsquare = (Z[stages+1]*Z[stages+1]) >> (s[stages+1]+p[stages+1]+1); // In the VHDL, first truncation of Z[stages+1]
		Log1p = Z[stages+1] - Zsquare;
		cout <<" Z"<<stages+1<<"Sq : ";     showMP(cout, "", Zsquare, 2*p[stages+1]+1, s[stages+1]-p[stages+1]-1); 

		cout <<"LogZ"<<stages+1<<" : ";     showMP(cout, "", Log1p, p[stages+1], s[stages+1]); 
		
		// Reconstruction

	
		// two's complement 
		Log =  lt0->function(A[0].get_ui());

		cout<<"   L"<<0<<" : "; 
		showMP(cout, "", Log , p[0], lt0->wOut); 
//     cout<<"   T0 : ";
//     printBinNumGMP(cout, Log, 90); cout<< endl;

		for(i=1; i<=stages; i++) {
			cout<<"   T"<<i<<" : "; 
			showMP(cout, "",lt[i]->function(A[i].get_ui()) , p[i], lt[i]->wOut); 

//       printBinNumGMP(cout, , 90); cout<< endl;
			Log += lt[i]->function(A[i].get_ui());
		// printBinNumGMP(cout, Log, 90); cout<< endl;
		}

		
		//Log = (Log<<(0*leadingzeroes)) + Log1p;
		Log = Log + Log1p;

		// Conversion to signed (sign extension)
		if(1==Log>>(target_prec -1))
			{
				Log = (mpz_class(1)<<target_prec) - Log;
				Log=-Log;
			} 

		showMP(cout, " LogM : ", Log , 0, target_prec); 

		mpfr_set_z(LogMCalc_m, Log.get_mpz_t(), GMP_RNDN);
		mpfr_shift_right(LogMCalc_m, target_prec );
	}
	
	// At this point, we have the log of the mantissa in LogMCalc_m

	// Add E*Log2
	mpfr_init2(Log2_m, wF+3);
	mpfr_log(Log2_m, two, GMP_RNDN);
	mpfr_shift_left(Log2_m, wF+3);

	mpz_init2(zlog2, wF+3);
	mpfr_get_z(zlog2, Log2_m, GMP_RNDN);
	Log2= mpz_class(zlog2);

	ELog2 = mpz_class((int)E) * Log2; 
	
	showMP(cout, "ELog2 ", ELog2, 0, 100/* wE+wF+3 */); 

	
	// Now conversion


	mpfr_init2(logM_m, wF+1);
	mpfr_log(logM_m, mantissa_m, GMP_RNDN);
	logM_d=mpfr_get_d(logM_m, GMP_RNDN);


	cout<<"precision="<<log(fabs((LogMCalc_d-logM_d)/logM_d)) / log(2)<<endl; 
	return LogMCalc_d;
}

void LogArch::output(std::ostream& o, std::string name, bool pipelined)
{
	int i;
	mpfr_t two;
	mpfr_t log2;
	mpz_t zlog2;

	

	o <<
		"-------------------------------------------------------------------------------\n"
		"-- Left barrel shifter\n"
		"--\n"
		"-- Generics:\n"
		"--   - w : width of the input operand\n"
		"--   - n : number of shifting stages (width of the signal s)\n"
		"--\n"
		"-- Ports:\n"
		"--   - i [in]  : input signal\n"
		"--   - s [in]  : shift by s\n"
		"--   - o [out] :  output\n"
		"--\n"
		"-- Recursive structure with n stages.\n"
		"-------------------------------------------------------------------------------\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"\n"
		"entity lshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"end entity;\n"
		"\n"
		"architecture arch of lshift is\n"
		"  component lshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"  end component;\n"
		"  signal o0 : std_logic_vector(w-1 downto 0);\n"
		"begin\n"
		"\n"
		"  check1: if 2**(n-1)>=w generate\n"
		"    o0 <=      i                      when s(n-1) = '0'\n"
		"          else (w-1 downto 0 => '0');\n"
		"  end generate;\n"
		"  check2: if 2**(n-1)<w generate\n"
		"    o0 <= i   when s(n-1) = '0' else\n"
		"          i(w-2**(n-1)-1 downto 0) & (2**(n-1)-1 downto 0 => '0');\n"
		"  end generate;\n"
		"  \n"
		"  ------------------------------------------------------------- Recursive stage\n"
		"\n"
		"  recursive : if n > 1 generate\n"
		"    shift0 : lshift\n"
		"      generic map ( w => w,\n"
		"                    n => n-1 )\n"
		"      port map (  i => o0,\n"
		"                  s => s(n-2 downto 0),\n"
		"                  o => o   );\n"
		"  end generate;\n"
		"\n"
		"  ----------------------------------------------------------------- Final stage\n"
		"  single : if n = 1 generate\n"
		"    o <= o0;\n"
		"  end generate;\n"
		"\n"
		"end architecture; -------------------------------------------------------------\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"-------------------------------------------------------------------------------\n"
		"-- Right barrel shifter\n"
		"--\n"
		"-- Generics:\n"
		"--   - w : width of the input operand\n"
		"--   - n : number of shifting stages (width of the signal s)\n"
		"--\n"
		"-- Ports:\n"
		"--   - i [in]  : input signal\n"
		"--   - s [in]  : shift by s\n"
		"--   - o [out] :  output\n"
		"--\n"
		"-- Recursive structure with n stages.\n"
		"-------------------------------------------------------------------------------\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"\n"
		"entity rshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"end entity;\n"
		"\n"
		"  architecture arch of rshift is\n"
		"  component rshift is\n"
		"  generic ( w : positive;\n"
		"            n : positive );\n"
		"  port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"         s  : in std_logic_vector(n-1 downto 0);\n"
		"         o  : out std_logic_vector(w-1 downto 0));\n"
		"  end component;\n"
		"  signal o0 : std_logic_vector(w-1 downto 0);\n"
		"begin\n"
		"\n"
		"  check1: if 2**(n-1)>=w generate\n"
		"    o0 <= i   when s(n-1) = '0' else\n"
		"     (w-1 downto 0 => '0');\n"
		"  end generate;\n"
		"  check2: if 2**(n-1)<w generate\n"
		"  o0 <= i   when s(n-1) = '0' else\n"
		"           (w-1 downto w-2**(n-1) => '0')  &  i(w-1 downto 2**(n-1));\n"
		"  end generate;\n"
		"  \n"
		"  ------------------------------------------------------------- Recursive stage\n"
		"\n"
		"  recursive : if n > 1 generate\n"
		"    shift0 : rshift\n"
		"      generic map ( w => w,\n"
		"                    n => n-1 )\n"
		"      port map (  i => o0,\n"
		"                  s => s(n-2 downto 0),\n"
		"                  o => o   );\n"
		"  end generate;\n"
		"\n"
		"  ----------------------------------------------------------------- Final stage\n"
		"  single : if n = 1 generate\n"
		"    o <= o0;\n"
		"  end generate;\n"
		"\n"
		"end architecture; -------------------------------------------------------------\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"-------------------------------------------------------------------------------\n"
		"-- Leading-zero counter and normalization\n"
		"--\n"
		"-- Generics:\n"
		"--   - w : width of the input operand\n"
		"--   - n : number of LZC and shifting stages (width of the signal z)\n"
		"--\n"
		"-- Ports:\n"
		"--   - i [in]  : input signal\n"
		"--   - z [out] : number of leading zeros\n"
		"--   - o [out] : normalized signal\n"
		"--\n"
		"-- Recursive structure with n stages. At most 2^n-1 leading zeros are counted.\n"
		"-------------------------------------------------------------------------------\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"\n"
		"entity lzc_norm is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i : in  std_logic_vector(w-1 downto 0);\n"
		"           z : out std_logic_vector(n-1 downto 0);\n"
		"           o : out std_logic_vector(w-1 downto 0) );\n"
		"end entity;\n"
		"\n"
		"architecture arch of lzc_norm is\n"
		"  component lzc_norm is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i : in  std_logic_vector(w-1 downto 0);\n"
		"           z : out std_logic_vector(n-1 downto 0);\n"
		"           o : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  signal z0 : std_logic;\n"
		"  signal o0 : std_logic_vector(w-1 downto 0);\n"
		"begin ---------------------------------------------- Test 2^(n-1) leading zeros\n"
		"  z0 <= '1' when i(w-1 downto w-2**(n-1)) = (w-1 downto w-2**(n-1) => '0') else\n"
		"        '0';\n"
		"\n"
		"  o0 <= i                                                    when z0 = '0' else\n"
		"        i(w-2**(n-1)-1 downto 0) & (2**(n-1)-1 downto 0 => '0');\n"
		"   z(n-1) <= z0;\n"
		"    ----------------------------------------------------------------- Final stage\n"
		"  single : if n = 1 generate\n"
		"    o <= o0;\n"
		"  end generate;\n"
		"  ------------------------------------------------------------- Recursive stage\n"
		"  recursive : if n > 1 generate\n"
		"    lzc_norm0 : lzc_norm\n"
		"      generic map ( w => w,\n"
		"                    n => n-1 )\n"
		"      port map ( i => o0,\n"
		"                 z => z(n-2 downto 0),\n"
		"                 o => o );\n"
		"  end generate;\n"
		"\n"
		"end architecture; -------------------------------------------------------------\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"\n"
		"library ieee;\n"
		"use ieee.std_logic_1164.all;\n"
		"use ieee.std_logic_arith.all;\n"
		"use ieee.std_logic_unsigned.all;\n"
		"\n"
		"-- Entity fp_log defined here\n"
		" \n";

	// TODO replace with constants
	o << "entity " << name << " is " << endl;
	o << "  generic ( wE : positive := " << wE  << ";" << endl;
	o << "            wF : positive := "  << wF << ");" << endl;
	o << "  port ( X : in  std_logic_vector(2+wE+wF downto 0);" << endl;
	o << "         R : out std_logic_vector(2+wE+wF downto 0)";
	if(pipelined)
		o  << ";" << endl << "         clk : in  std_logic  );" << endl;
	else
		o  << "  );" << endl;
		
	o << "end entity;" << endl;
	o << "architecture arch of "  << name << " is " << endl;

	o <<
		"\n"
		"  component lzoc is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i   : in  std_logic_vector(w-1 downto 0);\n"
		"           ozb : in std_logic;\n"
		"           zo  : out std_logic_vector(n-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  component lshift is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"           s  : in std_logic_vector(n-1 downto 0);\n"
		"           o  : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  component rshift is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i  : in  std_logic_vector(w-1 downto 0);\n"
		"           s  : in std_logic_vector(n-1 downto 0);\n"
		"           o  : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  component lzc_norm is\n"
		"    generic ( w : positive;\n"
		"              n : positive );\n"
		"    port ( i : in  std_logic_vector(w-1 downto 0);\n"
		"           z : out std_logic_vector(n-1 downto 0);\n"
		"           o : out std_logic_vector(w-1 downto 0) );\n"
		"  end component;\n"
		"\n"
		"  -- def of component range_red, and many constants and signals,  here\n";

	o <<   "  component range_red is port ("<<endl;
	o <<   "            Y0 : in  std_logic_vector("<<wF+1<<" downto 0);"<<endl;
	o <<   "            A  : in std_logic_vector("<< a[0]-1 <<" downto 0);"<<endl;
	if(pipelined) 
		o << "          clk  : in std_logic;"<<endl;
	o <<   "            Z  : out std_logic_vector("<< s[stages+1]-1 <<" downto 0);"<<endl;
	o <<   "    almostLog  : out std_logic_vector("<< lt0->wOut-1 <<" downto 0)  );"<<endl;
	o <<   "  end component;"<<endl;
	o <<   "  constant g   : positive := "<<gLog<<";"<<endl;
	o <<   "  constant a0 : positive := "<< a[0] <<";"<<endl;
	o <<   "  constant log2wF : positive := "<< intlog2(wF) <<";"<<endl;
	o <<   "  constant targetprec : positive := "<< target_prec <<";"<<endl;
	o <<   "  constant sfinal : positive := "<< s[stages+1] <<";"<<endl;
	o <<   "  constant pfinal : positive := "<< p[stages+1] <<";"<<endl;
	o <<   "  constant lzc_size : positive := "<< max(intlog2(wF), intlog2(wE+p[stages+1]+1)) << ";" << endl;

	// the maximum shift distance in the "small" path ?
	// o << "  constant shiftvalsize : positive := "<< intlog2(wF + 1 - p[stages+1]) <<";"<<endl;

	// The log2 constant
	mpfr_init2(two, 2);
	mpfr_set_d(two, 2.0, GMP_RNDN);
	mpfr_init2(log2, wF+gLog);
	mpfr_log(log2, two, GMP_RNDN);
	mpfr_shift_left(log2, wF+gLog);
	mpz_init2(zlog2, wF+gLog);
	mpfr_get_z(zlog2, log2, GMP_RNDN);
	o << "  signal log2 : std_logic_vector(wF+g-1 downto 0) := \"";
	printBinPosNumGMP(o, mpz_class(zlog2), wF+gLog);
	o << "\";"<<endl;
	o << "  signal E0offset : std_logic_vector(wE-1 downto 0) := \"";
	printBinPosNumGMP(o, (mpz_class(1)<<(wE-1)) -2 + wE , wE);
	o << "\"; -- E0 + wE "<<endl;
	o << "  signal pfinal_s : std_logic_vector(log2wF -1 downto 0) := \"";
	printBinPosNumGMP(o, mpz_class(p[stages+1]), intlog2(wF));
	o << "\";"<<endl;

	o <<
		"  signal FirstBit : std_logic;\n"
		"  signal Y0 : std_logic_vector(wF+1 downto 0);\n"
		"  signal E  : std_logic_vector(wE-1 downto 0);\n"
		"  signal absE  : std_logic_vector(wE-1 downto 0);\n"
		"  signal absELog2  : std_logic_vector(wF+wE+g-1 downto 0);\n"
		"  signal absELog2_pad, LogF_normal_pad, Log_normal, Log_normal_normd  : std_logic_vector(wE+targetprec-1 downto 0);\n"
		"  signal E_small,ER  : std_logic_vector(wE-1 downto 0);\n"
		"  signal E_normal  : std_logic_vector(lzc_size-1 downto 0);\n"
		"  signal Log_small_normd, Log_g  : std_logic_vector(wF+g-1 downto 0);\n"
		"  signal EFR  : std_logic_vector(wE+wF-1 downto 0);\n"
		"  signal lzo : std_logic_vector(log2wF-1 downto 0);\n"
		"  signal shiftval : std_logic_vector(log2wF downto 0);\n"
		"  signal absZ0, absZ0s : std_logic_vector(wF-pfinal+1 downto 0);\n"
		"  signal Zfinal, Log1p_normal : std_logic_vector(sfinal-1 downto 0);\n"
		"  signal Z2o2_full: std_logic_vector(2*(sfinal-pfinal) -1 downto 0);\n"
		"  signal squarerIn: std_logic_vector(sfinal-pfinal-1 downto 0);\n"
		"  signal Z2o2_small_s, Z2o2: std_logic_vector(sfinal-pfinal downto 0);\n"
		"  signal Log_small, Z_small, Z2o2_small: std_logic_vector(wF+g+1 downto 0);  \n"
		"  signal almostLog, logF_normal : std_logic_vector(targetprec-1 downto 0);\n"
		"  signal E0_sub : std_logic_vector(1 downto 0);\n"
		"  signal sR, small, doRR, ufl, sticky, round: std_logic;\n"
		"begin\n"
		"\n"
		"  FirstBit <=  X(wF-1);\n"
		"  Y0 <=      \"1\"  & X(wF-1 downto 0) & \"0\" when FirstBit = '0'\n"
		"        else \"01\" & X(wF-1 downto 0);\n"
		"\n"
		"  E  <= (X(wE+wF-1 downto wF)) - (\"0\" & (wE-2 downto 1 => '1') & (not FirstBit));\n"
		"\n"
		"  sR <= '0'   when    X(wE+wF-1 downto wF)   =   '0' &(wE-2 downto 0 => '1')  -- binade [1..2)\n"
		"        else not X(wE+wF-1);                -- MSB of exponent\n"
		"\n"
		"  absE <= ((wE-1 downto 0 => '0') - E)   when sR = '1'\n"
		"          else E;\n"
		"\n"
		"  absELog2 <= absE * log2;\n"
		"  \n"
		"  lzoc1 : lzoc\n"
		"    generic map (w => wF,  n => log2wF)\n"
		"    port map (  i => Y0(wF downto 1), ozb => FirstBit,  zo => lzo);\n"
		"\n"
		"  shiftval <= ('0' & lzo) - ('0' & pfinal_s); \n"
		"\n"
		"  doRR <= shiftval(log2wF);             -- sign of the result\n"
		"\n"
		"  small <= '1' when ((E=(wE-1 downto 0 => '0')) and (doRR='0'))\n"
		"          else '0';\n"
		"\n"
		"-- The range reduction instance\n"
		"  rr: range_red\n"
		"     port map ( A => X(wF-1 downto wF-a0), Y0 => Y0,\n";

	if(pipelined)
			o << "                clk=>clk,"<<endl;

	o <<
			"                Z => Zfinal, almostLog => almostLog);\n"
			"  absZ0 <=   Y0(wF-pfinal+1 downto 0)          when (sR='0') else\n"
			"             ((wF-pfinal+1 downto 0 => '0') - Y0(wF-pfinal+1 downto 0));\n"
			"\n"
			"--  absZ0 <=   Y0(wF-pfinal downto 0)   xor (wF-pfinal downto 0 => sR);\n"
			"\n"
			"  lshiftsmall: lshift          \n"
			"    generic map (w => wF-pfinal+2,  n => log2wF)    \n"
			"    port map (  i => absZ0, s => shiftval(log2wF-1 downto 0), o => absZ0s );\n"
			"\n"
			"  -- Z2o2 will be of size sfinal-pfinal, set squarer input size to that\n"
			"  sqintest: if sfinal > wf+2 generate\n"
			"    squarerIn <= Zfinal(sfinal-1 downto pfinal) when doRR='1'\n"
			"                 else (absZ0s &  (sfinal-wF-3 downto 0 => '0'));  \n"
			"  end generate sqintest;\n"
			"  sqintest2: if sfinal <= wf+2 generate\n"
			"    squarerIn <= Zfinal(sfinal-1 downto pfinal) when doRR='1'\n"
			"                 else absZ0s(wF-pfinal+1 downto wf+2-sfinal);  \n"
			"  end generate sqintest2;\n"
			"\n"
			"  -- Z2o2 will be of size sfinal - pfinal -1, set squarer input size to that\n"
			"--  sqintest: if sfinal >= wf+3 generate\n"
			"--    squarerIn <= Zfinal(sfinal-1 downto pfinal+1) when doRR='1'\n"
			"--                 else (absZ0s &  (sfinal-wF-4 downto 0 => '0'));  \n"
			"--  end generate sqintest;\n"
			"--  sqintest2: if sfinal < wf+3 generate\n"
			"--    squarerIn <= Zfinal(sfinal-1 downto pfinal+1) when doRR='1'\n"
			"--                 else absZ0s(wF-pfinal+1 downto wf+3-sfinal);  \n"
			"--  end generate sqintest2;\n"
			" \n"
			"  Z2o2_full <= (squarerIn * squarerIn);\n"
			"  Z2o2 <= Z2o2_full (2*(sfinal-pfinal)-1  downto sfinal-pfinal-1);\n"
			"\n"
			"  Log1p_normal  <=   Zfinal  -  ((sfinal-1 downto sfinal-pfinal-1  => '0') & (Z2o2(sfinal-pfinal downto 2)));\n"
			"\n"
			"  LogF_normal <=   almostLog + ((targetprec-1 downto sfinal => '0') & Log1p_normal);\n"
			"\n"
			"  absELog2_pad <=   absELog2 & (targetprec-wF-g-1 downto 0 => '0');       \n"
			"  LogF_normal_pad <= (wE-1  downto 0 => LogF_normal(targetprec-1))  & LogF_normal;\n"
			"  \n"
			"  Log_normal <=  absELog2_pad  + LogF_normal_pad when sR='0'  \n"
			"                else absELog2_pad - LogF_normal_pad;\n"
			"\n"
			"  lzc_norm_0 : lzc_norm\n"
			"    generic map (w => wE+targetprec, n => lzc_size)\n"
			"    port map (i => Log_normal, z => E_normal, o => Log_normal_normd);\n"
			"\n"
			"\n"
			"  rshiftsmall: rshift\n"
			"    generic map (w => sfinal-pfinal+1,  n => log2wF) \n"
			"    port map (i => Z2o2,\n"
			"              s => shiftval(log2wF-1 downto 0),\n"
			"              o => Z2o2_small_s);\n"
			"\n"
			"  -- send the MSB to position pfinal\n"
			"  Z2o2_small <=  (pfinal-1 downto 0  => '0') & Z2o2_small_s & (wF+g-sfinal downto 0  => '0') ;\n"
			"\n"
			"  -- mantissa will be either Y0-z^2/2  or  -Y0+z^2/2,  depending on sR  \n"
			"\n"
			"  Z_small <= (absZ0s & (pfinal+g-1 downto 0 => '0'));\n"
			"  Log_small  <=       Z_small -  Z2o2_small when (sR='0')\n"
			"                else  Z_small +  Z2o2_small;\n"
			"\n"
			"  -- Possibly subtract 1 or 2 to the exponent, depending on the LZC of Log_small\n"
			"  E0_sub <=      \"11\" when Log_small(wF+g+1) = '1'\n"
			"            else \"10\" when Log_small(wF+g+1 downto wF+g) = \"01\"\n"
			"            else \"01\" ;\n"
			"\n"
			"  E_small <=  \"0\" & (wE-2 downto 2 => '1') & E0_sub\n"
			"               - ((wE-1 downto log2wF => '0') & lzo) ;\n"
			"\n"
			"  Log_small_normd <= Log_small(wF+g+1 downto 2) when Log_small(wF+g+1)='1'\n"
			"             else Log_small(wF+g downto 1)  when Log_small(wF+g)='1'  -- remove the first zero\n"
			"             else Log_small(wF+g-1 downto 0)  ; -- remove two zeroes (extremely rare, 001000000 only)\n"
			"                                               \n"
			"  ER <= E_small when small='1'\n"
			"        else E0offset - ((wE-1 downto lzc_size => '0') & E_normal);\n"
			"  -- works only if wE > lzc_size approx log2wF, OK for usual exp/prec\n"
			"\n"
			"  Log_g  <=  Log_small_normd (wF+g-2 downto 0) & \"0\" when small='1'           -- remove implicit 1\n"
			"        else Log_normal_normd(wE+targetprec-2 downto wE+targetprec-wF-g-1 );  -- remove implicit 1\n"
			"\n"
			"  sticky <= '0' when Log_g(g-2 downto 0) = (g-2 downto 0 => '0') else\n"
			"            '1';\n"
			"  round <= Log_g(g-1) and (Log_g(g) or sticky);\n"
			"\n"
			"  -- use a trick: if round leads to a change of binade, the carry propagation\n"
			"  -- magically updates both mantissa and exponent\n"
			"  EFR <= (ER & Log_g(wF+g-1 downto g)) + ((wE+wF-1 downto 1 => '0') & round); \n"
			"\n"
			"\n"
			"  -- The smallest log will be log(1+2^{-wF}) \\approx 2^{-wF}\n"
			"  -- The smallest representable number is 2^{-2^(wE-1)} \n"
			"  -- Therefore, if \n"
			"--    underflow : if max(wE, log2(wF)+1) > wE generate\n"
			"--      ufl <=      '1' when (eR2(wE0-1) = '1') or (eR = (wE-1 downto 0 => '0'))\n"
			"--             else '0';\n"
			"--    end generate;\n"
			"\n"
			"--    no_underflow : if max(wE, log2(wE+wF)+2) = wE generate\n"
			"      ufl <= '0';\n"
			"--    end generate;\n"
			"\n"
			"  R(wE+wF+2 downto wE+wF) <= \"110\" when ((X(wE+wF+2) and (X(wE+wF+1) or X(wE+wF))) or (X(wE+wF+1) and X(wE+wF))) = '1' else\n"
			"                               \"101\" when X(wE+wF+2 downto wE+wF+1) = \"00\"                                                       else\n"
			"                               \"100\" when X(wE+wF+2 downto wE+wF+1) = \"10\"                                                       else\n"
			"                               \"00\" & sR when (((Log_normal_normd(wE+targetprec-1)='0') and (small='0')) or ( (Log_small_normd (wF+g-1)='0') and (small='1'))) or (ufl = '1') else\n"
			"                               \"01\" & sR;\n"
			"\n"
			"  R(wE+wF-1 downto 0) <=  EFR;\n"
			"\n"
			"end architecture;\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"\n"
			"-------------------------------------------------------------------------------\n"
			"-- Range Reduction box\n"
			"-------------------------------------------------------------------------------\n"
			"library ieee;\n"
			"use ieee.std_logic_1164.all;\n"
			"use ieee.std_logic_arith.all;\n"
			"use ieee.std_logic_unsigned.all;\n"
			"\n"
			"\n";

	//---------------Range reduction entity----------------------------------------------

	o <<   "entity range_red is port ("<<endl;
	o <<   "          Y0 : in  std_logic_vector("<<wF+1<<" downto 0);"<<endl;
	o <<   "          A  : in std_logic_vector("<< a[0]-1 <<" downto 0);"<<endl;
	if(pipelined) 
		o << "          clk  : in std_logic;"<<endl;
	o <<   "          Z  : out std_logic_vector("<< s[stages+1]-1 <<" downto 0);"<<endl;
	o <<   "  almostLog  : out std_logic_vector("<< lt0->wOut-1 <<" downto 0)  );"<<endl;
	o <<   "end entity;"<<endl<<endl;


	o << "architecture arch of  range_red is " <<endl<<endl;
	// All the components for the tables
	ostringstream nameIT0;
	nameIT0 <<"invtable0_" << wE << "_" <<wF;
	it0->outputComponent(o, nameIT0.str());

	ostringstream nameLT0;
	nameLT0 <<"logtable0_" << wE << "_" <<wF;
	lt0->outputComponent(o, nameLT0.str());

	for(i=1; i<=stages; i++) {
		ostringstream name3;
		name3 <<"logtable"<<i<<"_" << wE << "_" <<wF;
		lt[i]->outputComponent(o, name3.str()); 
	}


	for (i=0; i<= stages; i++) {
		o << "   signal       A"<<i<<":  std_logic_vector("<< a[i] - 1  <<" downto 0);"<<endl;
	}

	
	for (i=1; i<= stages; i++)
			o << "   signal       B"<<i<<":  std_logic_vector("<< s[i] - a[i] - 1  <<" downto 0);"<<endl;

	for (i=0; i<= stages+1; i++)
		o << "   signal Z"<<i<<", Z"<<i<<"_d:  std_logic_vector("<< s[i] - 1  <<" downto 0);"<<endl;

	for (i=1; i<= stages; i++) {
		o << "   signal    epsZ"<<i<<":  std_logic_vector("<< s[i]+p[i]+1  <<" downto 0);"<<endl;
	}
	// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>eps[i]
	// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
	// We need only this many bits of Z[i] to compute it.
	for (i=1; i<= stages; i++) {
		o << "   signal      ZM"<<i<<":  std_logic_vector("<< psize[i] - 1  <<" downto 0);"<<endl;
	}

	o << "   signal       P0:  std_logic_vector("<<  it0->wOut + s[0] - 1  <<" downto 0);"<<endl;

	for (i=1; i<= stages; i++) {
		o << "   signal       P"<<i<<":  std_logic_vector("<< psize[i]+a[i] - 1  <<" downto 0);"<<endl;
	}

	o << "   signal       L0:  std_logic_vector("<< lt0->wOut -1 <<" downto 0);"<<endl;

	for (i=1; i<= stages; i++)
			o << "   signal       L"<<i<<":  std_logic_vector("<< lt[i]->wOut -1 <<" downto 0);"<<endl;

	// Note that all the Si have the size of L0: absence of carry out is proven.
	for (i=1; i<= stages+1; i++)
		o << "   signal S"<<i<<", S"<<i<<"_d:  std_logic_vector("<< lt0->wOut-1 <<" downto 0);"<<endl;

	o << "   signal    InvA0:  std_logic_vector("<< a[0]  <<" downto 0);"<<endl;
	//  o << "   " <<endl;



	o << "begin" <<endl;
	//  o << "   A0 <= Fx (wF-1 downto wF-"<<a[0]<<");" <<endl;
	o << "   A0 <= A;"<<endl;
	o << "   it0:"<<nameIT0.str()<<" port map (x=>A0, y=>InvA0);" <<endl; 
	o << "   lt0:"<<nameLT0.str()<<" port map (x=>A0, y=>L0);"<<endl;
#if 0 // two pipeline stages per stage
		if(pipelined) 
			{
				o << "   -- Synchronization barrier 1 " <<endl;
				o << "   process(clk)  begin\n     if clk'event and clk='1' then"<<endl;
				o << "     invA0_d <= invA0;"<<endl;
				o << "     L0_d <= L0;"<<endl;
				o << "     end if;\n   end process;"<<endl<<endl;
			}
		else
			{
				o << "   invA0_d <= invA0;"<<endl;
				o << "   L0_d <= L0;"<<endl;
			}
#endif
		o << "   P0 <= InvA0 * Y0;" <<endl <<endl;
		if(pipelined) 
			{
				o << "   -- Synchronization barrier 1 " <<endl;
				o << "   process(clk)  begin\n     if clk'event and clk='1' then"<<endl;
				o << "     Z1_d <= P0("<< s[1] -1<<" downto 0);"<<endl;
				o << "     S1_d <= L0;"<<endl;
				o << "     end if;\n   end process;"<<endl<<endl;
			}
		else
			{
				o << "   Z1_d <= P0("<< s[1] -1<<" downto 0);"<<endl;
				o << "   S1_d <= L0;"<<endl;
			}


	for (i=1; i<= stages; i++) {

		o <<endl;
			//computation
		o << "   A"<<i<<" <= Z"<<i<<"_d(" << s[i] - 1  <<" downto "<< s[i] - a[i]  << ");"<<endl;
		o << "   B"<<i<<" <= Z"<<i<<"_d(" << s[i] - a[i] - 1  <<" downto 0 );"<<endl;
		o << "   lt"<<i<<":logtable"<<i<<"_"<< wE <<"_"<<wF<<" port map (x=>A"<<i<<", y=>L"<<i<<");"<<endl;
		if(psize[i] == s[i])
			o << "   ZM"<<i<<" <= Z"<<i<< "_d;"<<endl;   
		else
			o << "   ZM"<<i<<" <= Z"<<i<<"_d(" <<s[i]-1   <<" downto "<< s[i]-psize[i]  << ");"<<endl;   
		o << "   P"<<i<<" <= A"<<i<<"*ZM"<<i<<";"<<endl;

#if 0  // two pipeline stages per stage
		if (pipelined) 
			{
				o << "   -- Synchronization barrier "<<i+1 <<endl;
				o << "   process(clk)  begin\n     if clk'event and clk='1' then"<<endl;
				o << "     P"<<i<<"_d <= P" <<i<<";"<<endl;
				o << "     L"<<i<<"_d <= L" <<i<<";"<<endl;
				if(i==1) // delay L0 once more 
				  o << "     S1_d2 <= S1_d;"<<endl;
				else
				  o << "     S"<<i<<"_d2 <=   S"<<i<< "_d;"<<endl;
				o << "     end if;\n   end process;"<<endl<<endl;
			}
		else
			{
				o << "   P"<<i<<"_d <= P" <<i<<";"<<endl;
				o << "   L"<<i<<"_d <= L" <<i<<";"<<endl;
				if(i==1) // delay L0 once more 
				  o << "   S1_d2 <= S1_d;"<<endl;
				else
				  o << "   S"<<i<<"_d2 <=   S"<<i<< "_d;"<<endl<<endl;
			}
#endif

		if(i==1) // special case for the first iteration
			{
				o << "   epsZ"<<i<<" <= ("<<s[i]+p[i]+1<<" downto 0 => '0') "
				     << "     when  A1 = ("<<a[1]-1<<" downto 0 => '0')"<<endl
				     << "       else (\"01\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"_d )"
				     << "  when ((A1("<<a[1]-1<<")='0') and (A1("<<a[1]-2<<" downto 0) /= ("<<a[1]-2<<" downto 0 => '0')))"<<endl
				     << "       else "
				     << "(\"1\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"_d  & \"0\") "
				     << ";"<<endl;
			}
		else 
			{
				o << "   epsZ"<<i<<" <=  ("<< s[i]+p[i]+1<<" downto 0 => '0') "
				     << "     when  A"<<i<<" = ("<<a[i]-1<<" downto 0 => '0')"<<endl
				     << "     else    (\"01\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"_d);"<<endl;
			}

		o << "   Z"<<i+1<<" <=   (\"0\" & B"<<i;
		if (s[i+1] > 1+(s[i]-a[i]))  // need to padd Bi
			o << " & ("<<s[i+1] - 1-(s[i]-a[i]) -1<<" downto 0 => '0') ";    
		o <<")"<<endl 
				 << "         - ( ("<<p[i]-a[i]<<" downto 0 => '0') & P"<<i;
		// either pad, or truncate P
		if(p[i]-a[i]+1  + psize[i]+a[i]  < s[i+1]) // size of leading 0s + size of p 
			 o << " & ("<<s[i+1] - (p[i]-a[i]+1  + psize[i]+a[i]) - 1 <<" downto 0 => '0')";  // Pad
		if(p[i]-a[i]+1  + psize[i]+a[i]  > s[i+1]) 
			//truncate
			o <<"("<< psize[i]+a[i] - 1  <<" downto "<<  p[i]-a[i]+1  + psize[i]+a[i]  - s[i+1] << " )";
		o << "  )"<< endl;

		o << "         + epsZ"<<i << "("<<s[i]+p[i]+1<<" downto "<<s[i]+p[i] +2 - s[i+1]<<")"
				 << ";"<<endl;
			

		o << "   S"<<i+1<<" <=   S"<<i<<"_d + (("<<lt0->wOut-1<<" downto "<<lt[i]->wOut<<" =>'0') & L"<<i<<");"<<endl;

		o << endl;
		if(pipelined) 
			{
				o << "   -- Synchronization barrier "<<i+1 <<endl;
				o << "   process(clk)  begin\n     if clk'event and clk='1' then"<<endl;
				o << "     Z"<<i+1<<"_d <=   Z"<<i+1<< ";"<<endl;
				o << "     S"<<i+1<<"_d <=   S"<<i+1<< ";"<<endl;
				o << "     end if;\n   end process;"<<endl<<endl;
			}
		else
			{
				o << "   Z"<<i+1<<"_d <=   Z"<<i+1<< ";"<<endl;
				o << "   S"<<i+1<<"_d <=   S"<<i+1<< ";"<<endl<<endl;
			}
	}


	o << "   Z <= Z"<<stages+1<<"_d;"<<endl;  
	o << "   almostLog <= S"<<stages+1<<"_d;"<<endl;  

	o << "end architecture;" <<endl<<endl;

	// All the tables 
	ostringstream  name4;
	name4 <<"invtable0_" << wE << "_" <<wF;
	it0->output(o, name4.str());

	ostringstream name5;
	name5 <<"logtable0_" << wE << "_" <<wF;
	lt0->output(o, name5.str());

	for(i=1; i<=stages; i++) {
		ostringstream name6;
		name6<<"logtable"<<i<<"_" << wE << "_" <<wF;
		lt[i]->output(o, name6.str()); 
	}
}

int LogArch::size_in_LUTs() {
	int i,c;
	c=0;

	// first inv table 
	c += it0->size_in_LUTs();
	// first log table 
	c += lt0->size_in_LUTs();
	// first mult + add 
	c += (it0->wOut+1) * s[0];
	// first reconstruction add
	c += wF+gLog; //TODO
	
	// Multiplier for Z^2
	c += wF*wF/4;
	for(i=1; i<=stages; i++) {
		// table
		c += lt[i]->size_in_LUTs(); // plein de colonnes à 0 TODO
		// mult + add
		c += (a[i]+1)*s[i]; // pas tronqué... TODO
		// reconstruction add
		c += wF+gLog; //TODO
	}

	// The E*log2 multiplier
	c += wE*(wF+gLog)/2; // assuming semi-clever constant mult -- can do much better

	// The LZC/barrel shifter for large logs
	c += (wF)*(intlog2(wE)+1);

	// The LO1C
	c += wF;

	// The two lzo barrel shifters 
	c += 2*(wF/2)*(intlog2(wF/2)+1);

	// The two final multiplexers
	c += wE + wF;
	return c;
}
	
