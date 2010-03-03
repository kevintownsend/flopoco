/*
 * Dot Product unit for FloPoCo
  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Author :   Bogdan Pasca

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL License, 2008-2010.
 */

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

#include "DotProduct.hpp"


using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	DotProduct::DotProduct(Target* target, int wE, int wFX, int wFY, int MaxMSBX, int LSBA, int MSBA ):
		Operator(target), wE(wE), wFX(wFX), wFY(wFY), MaxMSBX(MaxMSBX), LSBA(LSBA), MSBA(MSBA)  {
	
		ostringstream name;
		/* Set up the name of the entity */
		name <<"DotProduct_"<<wE<<"_"<<wFX<<"_"<<wFY<<"_"
			  <<(MaxMSBX>=0?"":"M")<<abs(MaxMSBX)<<"_"
			  <<(LSBA>=0?"":"M")<<abs(LSBA)<<"_"
			  <<(MSBA>=0?"":"M")<<abs(MSBA) ;

		setName(name.str());
	
		/* Set up the I/O signals of of the entity */
		addInput("X", 2 + 1 + wE + wFX);
		addInput("Y", 2 + 1 + wE + wFY);
		addOutput("A", MSBA-LSBA+1); //the width of the output represents the accumulator size


		/* This operator is sequential.*/
		setSequential();
 
		sizeAcc_ = MSBA-LSBA+1;
		/* Instantiate one FPMultiplier used to multiply the two inputs fp numbers, X and Y */
		int fpMultiplierResultExponentWidth = wE;
		int fpMultiplierResultFractionWidth = wFX + wFY + 1;
		int fpMultiplierResultNormalization = 0;
		fpMultiplier = new FPMultiplier(target, wE, wFX, wE, wFY, fpMultiplierResultExponentWidth, fpMultiplierResultFractionWidth, fpMultiplierResultNormalization);
		oplist.push_back(fpMultiplier);
  
		/* Instantiate one LongAcc used accumulate the products produced by fpMultiplier */
		longAcc = new LongAcc(target, wE, fpMultiplierResultFractionWidth, MaxMSBX, LSBA, MSBA);
		oplist.push_back(longAcc);
  
		/* Signal declaration */
		addSignal("fpMultiplierResultExponent",wE);
		addSignal("fpMultiplierResultSignificand", fpMultiplierResultFractionWidth + 1);
		addSignal("fpMultiplierResultException", 2);
		addSignal("fpMultiplierResultSign",1);
		addSignal("fpMultiplierResult",2 + 1 + wE + fpMultiplierResultFractionWidth);
   
		/* Set the pipeline depth of this operator */
		setPipelineDepth(fpMultiplier->getPipelineDepth() + longAcc->getPipelineDepth() );
	}

	DotProduct::~DotProduct() {
	}



	void DotProduct::outputVHDL(std::ostream& o, std::string name) {
		licence(o,"Bogdan Pasca (2008)");
		Operator::stdLibs(o);
		outputVHDLEntity(o);
		newArchitecture(o, name);
		fpMultiplier->outputVHDLComponent(o); //output the code of the fpMultiplier component
		longAcc->outputVHDLComponent(o);      //output the code of the longAcc component
		outputVHDLSignalDeclarations(o);
		beginArchitecture(o);

		/* Map the signals on the fp multiplier */
		o << tab << "fp_multiplier: " << fpMultiplier->getName() << endl;
		o << tab << "    port map ( X => X, " << endl;
		o << tab << "               Y => Y, " << endl;
		o << tab << "               ResultExponent => fpMultiplierResultExponent, " << endl;
		o << tab << "               ResultSignificand => fpMultiplierResultSignificand, " << endl;
		o << tab << "               ResultException => fpMultiplierResultException, " << endl;
		o << tab << "               ResultSign => fpMultiplierResultSign, " << endl;
		o << tab << "               clk => clk," << endl;
		o << tab << "               rst => rst " << endl;
		o << tab << "                         );" << endl; 
		o << endl;

		/*rebuild the output under the form of a fp number */
		o << tab << "fpMultiplierResult <= fpMultiplierResultException & "
		  <<"fpMultiplierResultSign & "
		  <<"fpMultiplierResultExponent & "
		  <<"fpMultiplierResultSignificand("<<wFX + wFY+1<<" downto "<< 1 <<");"<<endl;

		/* Map the signals on the long accumulator */
		o << tab << "long_acc: " << longAcc->getName() << endl;
		o << tab << "    port map ( X => fpMultiplierResult, " << endl;
		o << tab << "               A => A, " << endl;
		o << tab << "               clk => clk," << endl;
		o << tab << "               rst => rst " << endl;
		o << tab << "                         );" << endl; 
		o << endl;

		endArchitecture(o);
	} 

	void DotProduct::test_precision(int n) {
		mpfr_t ref_acc, long_acc, fp_acc, x, y, r, d, one, two, msb;
		double sum, error;

		// initialisation
		mpfr_init2(ref_acc, 10000);
		mpfr_init2(long_acc, sizeAcc_+1);
		mpfr_init2(fp_acc, wFX+1);
		mpfr_init2(x, wFX+1);
		mpfr_init2(y, wFY+1);
	
		mpfr_init2(r, 2*(wFX+1) );
		mpfr_init2(d, 10000);
		mpfr_init2(one, 100);
		mpfr_init2(two, 100);
		mpfr_init2(msb, 10);
		mpfr_set_d(one, 1.0, GMP_RNDN);
		mpfr_set_d(two, 2.0, GMP_RNDN);
		mpfr_set_d(msb, (double)(1<<(MSBA+1)), GMP_RNDN); // works for MSBA_<32
		mpfr_mul_2si(msb, one, MSBA+1, GMP_RNDN); // works for MSBA_<32

		cout <<endl;
		mpfr_out_str(stdout, 2, sizeAcc_, msb, GMP_RNDN);
		cout <<endl;

		//cout<<"%-------Acc. of positive numbers--------------- "<<endl;
		mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
		mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
		mpfr_set_d(long_acc, 0.0, GMP_RNDN);
		//put a one in the MSBA_+1 bit will turn all the subsequent additions
		// into fixed-point ones
		mpfr_add(long_acc, long_acc, msb, GMP_RNDN);

		gmp_randstate_t state;
		gmp_randinit_default(state);

		for(int i=0; i<n; i++){
			mpfr_urandomb(d,state); // deprecated function; r is [0,1]
			mpfr_set(x,d,GMP_RNDN);		
			mpfr_urandomb(d,state); // deprecated function; r is [0,1]
			mpfr_set(y,d,GMP_RNDN);		

			mpfr_mul(r,x,y,GMP_RNDN);//    mpfr_add(r, one, r, GMP_RNDN); // r in [1 2[
		
			mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
			mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);		
			mpfr_add(long_acc, long_acc, r, GMP_RNDZ);
			
		}
		cout <<endl;
		mpfr_out_str(stdout, 2, sizeAcc_+1, long_acc, GMP_RNDN);
		cout <<endl;
		mpfr_out_str(stdout, 2, sizeAcc_+1, ref_acc, GMP_RNDN);
		cout <<endl;

		// remove the leading one from long acc

		mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);

		cout <<endl;
		mpfr_out_str(stdout, 2, sizeAcc_+1, long_acc, GMP_RNDN);
		cout <<endl;
		mpfr_out_str(stdout, 2, sizeAcc_+1, ref_acc, GMP_RNDN);
		cout <<endl;

		cout << "*************----***********"<<endl;		
		sum=mpfr_get_d(ref_acc, GMP_RNDN);
		cout  << "% unif[0 1] :sum="<< sum;
		sum=mpfr_get_d(fp_acc, GMP_RNDN);
		cout << "   FPAcc="<< sum;
		sum=mpfr_get_d(long_acc, GMP_RNDN);
		cout << "   LongAcc="<< sum;

		cout <<endl << n << " & ";
		// compute the error for the FP adder
		mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		// cout << " Relative error between fp_acc and ref_acc is "<< error << endl;
		cout << scientific << setprecision(2)  << error << " &1 ";
		// compute the error for the long acc
		mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);

		error=mpfr_get_d(d, GMP_RNDN);
		cout << endl<< " Diff "<< scientific << setprecision(20)<< error << endl;

		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		//  cout << "Relative error between long_acc and ref_acc is "<< error << endl;
		cout << scientific << setprecision(2)  << error << " &2 ";

		/*
		//cout<<"%-------Acc. of positive/negative numbers--------------- "<<endl;

		mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
		mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
		//put a one in the MSBA_+1 bit will turn all the subsequent additions
		// into fixed-point ones
		mpfr_set_d(long_acc, (double)(1<<(MSBA+1)), GMP_RNDN);

		for(int i=0; i<n; i++){
		mpfr_random(r); // deprecated function; r is [0,1]
		mpfr_mul(r, r, two, GMP_RNDN); 
		mpfr_sub(r, r, one, GMP_RNDN); 
		
		mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
		mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);
		mpfr_add(long_acc, long_acc, r, GMP_RNDN);
		}

		// remove the leading one from long acc
		mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);


		// compute the error for the FP adder
		mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		// cout << "Relative error between fp_acc and ref_acc is "<< error << endl;
		cout << scientific << setprecision(2) << error << " & ";

		// compute the error for the long acc
		mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		//  cout << "Relative error between long_acc and ref_acc is "<< error << endl;
		cout << scientific << setprecision(2)  << error << " \\\\ \n     \\hline \n";

		sum=mpfr_get_d(ref_acc, GMP_RNDN);
		cout << "% unif[-1 1] : sum="<< sum;
		sum=mpfr_get_d(fp_acc, GMP_RNDN);
		cout << "   FPAcc="<< sum;
		sum=mpfr_get_d(long_acc, GMP_RNDN);
		cout << "   LongAcc="<< sum;
		cout <<endl;
		*/
	}

}
