/*
 * A long accumulator for FloPoCo
 *
 * Authors : Florent de Dinechin and Bogdan Pasca
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
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "utils.hpp"
#include "Operator.hpp"
#include "LongAcc.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;
	extern int LongAccN;

	LongAcc::LongAcc(Target* target, int wEX, int wFX, int MaxMSBX, int LSBA, int MSBA): 
		Operator(target), 
		wEX_(wEX), wFX_(wFX), MaxMSBX_(MaxMSBX), LSBA_(LSBA), MSBA_(MSBA), AccValue_(0), xOvf(0)
	{
		int i;
		setCopyrightString("Florent de Dinechin, Bogdan Pasca (2008-2009)");		

		//check input constraints, i.e, MaxMSBX <= MSBA, LSBA<MaxMSBx
		if ((MaxMSBX_ > MSBA_)){
			cerr << " LongAcc: Input constraint MaxMSBX <= MSBA not met." << endl;
			exit (EXIT_FAILURE);
		}
		if ((LSBA_ >= MaxMSBX_)){
			cerr << " LongAcc: Input constraint LSBA<MaxMSBx not met:"<<
				" This accumulator would never accumulate a bit." << endl;
			exit (EXIT_FAILURE);
		}

		ostringstream name; 
		name <<"LongAcc_"<<wEX_<<"_"<<wFX_<<"_"
			  <<(MaxMSBX_>=0?"":"M")<<abs(MaxMSBX_)<<"_"
			  <<(LSBA_>=0?"":"M")<<abs(LSBA_)<<"_"
			  <<(MSBA_>=0?"":"M")<<abs(MSBA_) ;
		setName(name.str());

		// This operator is a sequential one
		setSequential();

		// Set up various architectural parameters
		sizeAcc_ = MSBA_-LSBA_+1;

		addFPInput ("X", wEX_,wFX_);
		addInput   ("newDataSet");
		addOutput  ("ready");
		addOutput  ("A", sizeAcc_);  
		addOutput  ("C", sizeAcc_);
		addOutput  ("XOverflow");  
		addOutput  ("XUnderflow");  
		addOutput  ("AccOverflow");  
		
		maxShift_        = MaxMSBX_-LSBA_;              // shift is 0 when the implicit 1 is at LSBA_
		sizeShift_       = intlog2(maxShift_);         // the number of bits needed to control the shifter
		sizeSummand_     = MaxMSBX_-LSBA_+1;         // the size of the summand (the maximum one - when the inplicit 1 is on MaxMSBX)
		sizeShiftedFrac_ = maxShift_ + wFX_+1;   
		E0X_ = (1<<(wEX_-1)) -1;                 // exponent bias
		// Shift is 0 when implicit 1 is on LSBA, that is when EX-bias = LSBA
		// that is, EX-bias-LSBA = 0, EX-(bias + LSBA) = 0
		// We are working in sign-magnitude, so we differentiate 2 cases:
		// 1. when bias+LSBA >=0, that is, sign bit is 0, then shiftValue = Exp - (bias+LSBA) ;  (0)Exp - (0)(bias+LSBA)
		// 2. when bias+LSBA <=0, which means that the sign bit is 1. Because in this case we have a substraction of a negative number,
		// we replace it whit an addition of -(bias+LSBA), that is (0)Exp + (0)(- (bias + LSBA))
		//	
		// If the shift value is negative, then the input is shifted out completely, and a 0 is added to the accumulator 
		// All above computations are performed on wEx+1 bits	
	
		//MaxMSBx is one valid exponent value, that is, 
		//1. after bias is added, value should be >= 0
		//2. after bias is added, representation should still fit on no more than wEX bits
		int biasedMaxMSBX_ = MaxMSBX_ + E0X_;
		if(biasedMaxMSBX_ < 0 || intlog2(biasedMaxMSBX_)>wEX_) {
			cerr<<"ERROR in LongAcc: MaxMSBX_="<<MaxMSBX_<<" is not a valid exponent of X (range "
				 <<(-E0X_)<< " to " << ((1<<wEX_)-1)-E0X_ <<endl;
			exit (EXIT_FAILURE);
		}
		
		int chunkSize_;
		target->suggestSubaddSize(chunkSize_ , sizeAcc_);
		int nbOfChunks = ceil(double(sizeAcc_)/double(chunkSize_));
		int lastChunkSize = ( sizeAcc_ % chunkSize_ == 0 ? chunkSize_  : sizeAcc_ % chunkSize_);

		vhdl << tab << declare("fracX",wFX_+1) << " <=  \"1\" & X" << range(wFX_-1,0) << ";" << endl;
		vhdl << tab << declare("expX" ,wEX_  ) << " <= X" << range(wEX_+wFX_-1,wFX_) << ";" << endl;
		vhdl << tab << declare("signX",1     ) << " <= X" << of(wEX_+wFX_) << ";" << endl;
		vhdl << tab << declare("exnX" ,2     ) << " <= X" << range(wEX_+wFX_+2,wEX_+wFX_+1) << ";" << endl;
	
		vhdl << tab << declare("xOverflowCond",1) << " <= '1' when (( expX > CONV_STD_LOGIC_VECTOR("<<MaxMSBX_ + E0X_<<","<< wEX_<<")) or (exnX >= \"10\")) else '0' ;"<<endl; 
		vhdl << tab << declare("xUnderflowCond",1) << " <= '1' when (expX < CONV_STD_LOGIC_VECTOR("<<LSBA_ + E0X_<<","<<wEX_<<")) else '0' ;" << endl;  
		//determination of the shift value
		int64_t exp_offset = E0X_+LSBA_;
		vhdl << tab << declare("shiftVal",wEX_+1) << " <= (\"0\" & expX) - CONV_STD_LOGIC_VECTOR("<< exp_offset <<","<<  wEX_+1<<");" << endl;

		//input shifter mappings
		shifter_ = new Shifter(target, wFX_+1, maxShift_, Shifter::Left);
		oplist.push_back(shifter_);

		inPortMap   (shifter_, "X", "fracX");
		inPortMapCst(shifter_, "S", use("shiftVal")+range(shifter_->getShiftInWidth() - 1,0));
		outPortMap  (shifter_, "R", "shifted_frac");
		vhdl << instance(shifter_, "LongAccInputShifter");
	
		syncCycleFromSignal("shifted_frac");
		//determine if the input has been shifted out from the accumulator. 
		//in this case the accumulator will be incremented with 0
		vhdl << tab << declare("flushedToZero",1) << " <= '1' when (shiftVal" << of(wEX_)<<"='1' -- negative left shift " << endl;
		vhdl << tab << "                               or exnX=\"00\")" << endl;
		vhdl << tab << "                 else '0';" << endl;
		vhdl << tab << declare("summand", sizeSummand_, true, Signal::registeredWithSyncReset) << "<= " << rangeAssign(sizeSummand_-1,0,"'0'") << " when flushedToZero='1' " << 
			" else shifted_frac" << range(sizeShiftedFrac_-1,wFX_)<<";" << endl;

		vhdl << tab << "-- 2's complement of the summand" << endl;
		// Don't compute 2's complement just yet, just invert the bits and leave the addition of the extra 1 in accumulation, as a carry in bit
		vhdl << tab << declare("summand2c", sizeSummand_, true) << " <= summand when (signX='0' or flushedToZero='1') "<<
			"else not(summand);"<< endl;

		vhdl << tab << "-- extension of the summand to accumulator size" << endl;
		vhdl << tab << declare("ext_summand2c",sizeAcc_,true) << " <= " << (sizeAcc_-1<sizeSummand_?"":rangeAssign(sizeAcc_-1, sizeSummand_, use("signX")+ " and not " + use("flushedToZero")) + " & " ) << "summand2c;" << endl;

		vhdl << tab << "-- accumulation itself" << endl;
		//determine the value of the carry in bit
		vhdl << tab << declare("carryBit_0",1,false) << " <= signX and not flushedToZero;" << endl; 
		
		
		for (i=0; i < nbOfChunks; i++) {
			ostringstream accReg;
			accReg<<"acc_"<<i;

			vhdl << tab << declare(join("acc_",i),(i!=nbOfChunks-1?chunkSize_:lastChunkSize) ,true, Signal::registeredWithSyncReset) << " <= " << 
				join("acc_",i,"_ext")<<range((i!=nbOfChunks-1?chunkSize_-1:lastChunkSize-1),0) << ";" << endl;
			vhdl << tab << declare(join("carryBit_",i+1),1, false, Signal::registeredWithSyncReset) <<"  <= " << join("acc_",i,"_ext")<<of((i!=nbOfChunks-1?chunkSize_:lastChunkSize)) << ";" << endl;
			nextCycle();		
			vhdl << tab << declare(join("acc_",i,"_ext"),(i!=nbOfChunks-1?chunkSize_:lastChunkSize)+1) << " <= ( \"0\" & (" << use(join("acc_",i)) << " and "<<rangeAssign( (i!=nbOfChunks-1?chunkSize_:lastChunkSize)-1,0, "not(newDataSet)") <<")) + " <<
				"( \"0\" & ext_summand2c" << range( (i!=nbOfChunks-1?chunkSize_*(i+1)-1:sizeAcc_-1), chunkSize_*i) << ") + " << 
				use(join("carryBit_",i)) << ";" << endl;
			setCycleFromSignal("carryBit_0");
		}

		nextCycle();
		//compose the acc signal 
		
		vhdl << tab << declare("acc", sizeAcc_) << " <= ";
		for (i=nbOfChunks-1;i>=0;i--) {
			if (i>0)
				vhdl << use(join("acc_",i)) << " & ";
			else
				vhdl << use(join("acc_",i)) << ";" << endl;
		}
	
		vhdl << tab << declare("carry", sizeAcc_) << " <= ";
	
		for (i=nbOfChunks-1;i>=0; i--){
			if (i<nbOfChunks-1){
				vhdl << use( join("carryBit_",i+1) ); 
			}

			if (i<nbOfChunks-1)
				vhdl << " & ";
		
			if (i==nbOfChunks-1)			
				{
					if(nbOfChunks-1==0&&lastChunkSize-1<sizeAcc_)
						vhdl << zg(lastChunkSize,0);
					else			
						vhdl << zg(lastChunkSize-1,0);
				}
			else
				if (i>0)
					vhdl << zg(chunkSize_-1,0);
				else
					vhdl << zg(chunkSize_,0);

			if (i>0)
				vhdl << " & ";

		}
		vhdl << ";" << endl;

		vhdl << tab << declare("xOverflowRegister",1) << " <= "; nextCycle(false); vhdl << "xOverflowRegister or xOverflowCond;"<<endl;
		setCycleFromSignal("acc",false);
		vhdl << tab << declare("xUnderflowRegister",1) << " <= "; nextCycle(false); vhdl << "xUnderflowRegister or xUnderflowCond;"<<endl;
		setCycleFromSignal("acc",false);
		vhdl << tab << declare("accOverflowRegister",1) << " <= "; nextCycle(false); vhdl << "accOverflowRegister or "<<use(join("carryBit_",nbOfChunks)) << ";"<<endl;
		setCycleFromSignal("acc",false);

		vhdl << tab << "A <=   acc;" << endl;
		vhdl << tab << "C <=   carry;" << endl;

		nextCycle();
		//if accumulator overflows this flag will be set to 1 until a maunal reset is performed
		vhdl << tab << "AccOverflow <= accOverflowRegister;"<<endl; 
		//if the input overflows then this flag will be set to 1, and will remain 1 until manual reset is performed
		vhdl << tab << "XOverflow <= xOverflowRegister;"<<endl; 
		vhdl << tab << "XUnderflow <= xUnderflowRegister;"<<endl; 

		setCycleFromSignal("acc");
		vhdl << tab << "ready <= newDataSet;"<<endl;
	}


	LongAcc::~LongAcc() {
	}

	void LongAcc::test_precision(int n) {
		mpfr_t ref_acc, long_acc, fp_acc, r, d, one, two, msb;
		double sum, error;

		// initialisation
		mpfr_init2(ref_acc, 10000);
		mpfr_init2(long_acc, sizeAcc_+1);
		mpfr_init2(fp_acc, wFX_+1);
		mpfr_init2(r, wFX_+1);
		mpfr_init2(d, 100);
		mpfr_init2(one, 100);
		mpfr_init2(two, 100);
		mpfr_init2(msb, 100);
		mpfr_set_d(one, 1.0, GMP_RNDN);
		mpfr_set_d(two, 2.0, GMP_RNDN);
		mpfr_set_d(msb, (double)(1<<(MSBA_+1)), GMP_RNDN); // works for MSBA_<32

		//cout<<"%-------Acc. of positive numbers--------------- "<<endl;
		mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
		mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
		//put a one in the MSBA_+1 bit will turn all the subsequent additions
		// into fixed-point ones
		mpfr_set_d(long_acc, (double)(1<<(MSBA_+1)), GMP_RNDN);


		for(int i=0; i<n; i++){
			//mpfr_random(r); // deprecated function; r is [0,1]
			mpfr_urandomb(r, FloPoCoRandomState::m_state);
			//    mpfr_add(r, one, r, GMP_RNDN); // r in [1 2[
		
			mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
			mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);
			mpfr_add(long_acc, long_acc, r, GMP_RNDN);
			if(mpfr_greaterequal_p(long_acc, msb)) 
				mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);
			
		}

		// remove the leading one from long acc
		if(mpfr_greaterequal_p(long_acc, msb)) 
			mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);
		
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
		cout << scientific << setprecision(2)  << error << " & ";
		// compute the error for the long acc
		mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		//  cout << "Relative error between long_acc and ref_acc is "<< error << endl;
		cout << scientific << setprecision(2)  << error << " & ";


		//cout<<"%-------Acc. of positive/negative numbers--------------- "<<endl;

		mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
		mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
		//put a one in the MSBA_+1 bit will turn all the subsequent additions
		// into fixed-point ones
		mpfr_set_d(long_acc, (double)(1<<(MSBA_+1)), GMP_RNDN);

		for(int i=0; i<n; i++){
			//mpfr_random(r); // deprecated function; r is [0,1]
			mpfr_urandomb(r, FloPoCoRandomState::m_state);
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

	}

	// read the values from a file and accumulate them
	void LongAcc::test_precision2() {

		mpfr_t ref_acc, long_acc, fp_acc, r, d, one, two, msb;
		double dr,  sum, error;

		// initialisation
#define DPFPAcc 1
		mpfr_init2(ref_acc, 10000);
		mpfr_init2(long_acc, sizeAcc_+1);
#if DPFPAcc
		mpfr_init2(fp_acc,52 +1);
#else
		mpfr_init2(fp_acc, wFX_+1);
#endif
		mpfr_init2(r, wFX_+1);
		mpfr_init2(d, 100);
		mpfr_init2(one, 100);
		mpfr_init2(two, 100);
		mpfr_init2(msb, 100);
		mpfr_set_d(one, 1.0, GMP_RNDN);
		mpfr_set_d(two, 2.0, GMP_RNDN);
		mpfr_set_d(msb, (double)(1<<(MSBA_+1)), GMP_RNDN); // works for MSBA_<32

		mpfr_set_d(ref_acc, 0.0, GMP_RNDN);
		mpfr_set_d(fp_acc, 0.0, GMP_RNDN);
		//put a one in the MSBA_+1 bit will turn all the subsequent additions
		// into fixed-point ones
		mpfr_set_d(long_acc, (double)(1<<(MSBA_+1)), GMP_RNDN);


		ifstream myfile ("/home/fdedinec/accbobines.txt");
		int i=0;
		if (myfile.is_open())
			{
				while (! myfile.eof() )
					{
						myfile >> dr;
						mpfr_set_d(r, dr, GMP_RNDN);
						mpfr_add(fp_acc, fp_acc, r, GMP_RNDN);
						mpfr_add(ref_acc, ref_acc, r, GMP_RNDN);
						mpfr_add(long_acc, long_acc, r, GMP_RNDN);
						i++;
						if (i%100000==0) {
							cout<<"i="<<i<<" "<<endl;

							// compute the error for the FP adder
							mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
							mpfr_div(d, d, ref_acc, GMP_RNDN);
							error=mpfr_get_d(d, GMP_RNDN);
							cout << "Relative error between fp_acc and ref_acc is "<< error << endl;
							//cout << scientific << setprecision(2) << error << " & ";

							// remove the leading one from long acc
							mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);
							// compute the error for the long acc
							mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
							mpfr_div(d, d, ref_acc, GMP_RNDN);
							error=mpfr_get_d(d, GMP_RNDN);
							cout << "Relative error between long_acc and ref_acc is "<< error << endl;
							//cout << scientific << setprecision(2)  << error << " \\\\ \n     \\hline \n";

							sum=mpfr_get_d(ref_acc, GMP_RNDN);
							cout << " exact sum="<< sum;
							sum=mpfr_get_d(fp_acc, GMP_RNDN);
							cout << "   FPAcc="<< sum;
							sum=mpfr_get_d(long_acc, GMP_RNDN);
							cout << "   LongAcc="<< sum;
							cout <<endl;
							// add the leading one back
							mpfr_add(long_acc, long_acc, msb, GMP_RNDN);
						}
					}
				myfile.close();
			}
		// remove the leading one from long acc
		mpfr_sub(long_acc, long_acc, msb, GMP_RNDN);


		// compute the error for the FP adder
		mpfr_sub(d, fp_acc, ref_acc, GMP_RNDN);
		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		cout << "Relative error between fp_acc and ref_acc is "<< error << endl;
		//cout << scientific << setprecision(2) << error << " & ";

		// compute the error for the long acc
		mpfr_sub(d, long_acc, ref_acc, GMP_RNDN);
		mpfr_div(d, d, ref_acc, GMP_RNDN);
		error=mpfr_get_d(d, GMP_RNDN);
		cout << "Relative error between long_acc and ref_acc is "<< error << endl;
		//cout << scientific << setprecision(2)  << error << " \\\\ \n     \\hline \n";

		sum=mpfr_get_d(ref_acc, GMP_RNDN);
		cout << " exact sum="<< sum;
		sum=mpfr_get_d(fp_acc, GMP_RNDN);
		cout << "   FPAcc="<< sum;
		sum=mpfr_get_d(long_acc, GMP_RNDN);
		cout << "   LongAcc="<< sum;
		cout <<endl;


		exit(0);

	}




	//FIXME 	: IOs are defined as
   // addFPInput ("X", wEX_,wFX_);
	// addOutput  ("A", sizeAcc_);  
	// addOutput  ("XOverflow");  
	// addOutput  ("XUnderflow");  
	// addOutput  ("AccOverflow");  
	// In the following there is no check on the two last ones

	void LongAcc::fillTestCase(mpz_class a[])
	{
		mpz_class& sX = a[0];
		mpz_class& sA = a[1];
		mpz_class& sXOverflow = a[2];
	
		currentIteration++;

		//the last additionNumberOfChunks_-1 inputs must be 0 in order to propagate the carry bit in the accumulator
		if (currentIteration  > LongAccN - additionNumberOfChunks_ + 1)
			sX = 0;

		//convert this raw number into a FP number accornding to the flopoco FP format
		FPNumber fpX(wEX_, wFX_);
		fpX = sX;
	
		//for now, we do not accept inputs which overflow, negative inputs or inputs having the exception bits INF and NaN
		while ((fpX.getExponentSignalValue()-(intpow2(wEX_-1)-1)>MaxMSBX_)
				 //||(fpX.getSignSignalValue()==1) 
				 || (fpX.getExceptionSignalValue()>1)){
			sX = getLargeRandom(wEX_+wFX_+3);
			fpX = sX;
		}

	
		if ((fpX.getExceptionSignalValue()>1) || (fpX.getExponentSignalValue()-(intpow2(wEX_-1)-1)>MaxMSBX_))
			xOvf = 1;
	
		sXOverflow = xOvf;
	
		if (verbose==1)
			cout<<" i="         << currentIteration
#ifdef  _WIN32	
				 <<",a=" << mpz2string(AccValue_);
#else
		<<",a=" << AccValue_;
#endif

		if (fpX.getExceptionSignalValue()!=0)  {
			if (fpX.getSignSignalValue()==0)
				AccValue_ = AccValue_ + mapFP2Acc(fpX);
			else
				AccValue_ = AccValue_ - mapFP2Acc(fpX);
		}
		if (verbose==1){
			if (fpX.getExceptionSignalValue()!=0)
#ifdef  _WIN32
				cout<< (fpX.getSignSignalValue()==0?" + ":" - ")<< mpz2string(mapFP2Acc(fpX)) << " = "<<mpz2string(AccValue_)<<";"; 
#else 
			cout<< (fpX.getSignSignalValue()==0?" + ":" - ")<< mapFP2Acc(fpX) << " = "<<AccValue_<<";"; 
#endif

#ifdef  _WIN32
			cout<< "Exc="<<mpz2string(fpX.getExceptionSignalValue())
				 <<", S="<<mpz2string(fpX.getSignSignalValue())
				 <<", Exp="<<mpz2string(fpX.getExponentSignalValue()-(intpow2(wEX_-1)-1))
				 <<", Frac="<<mpz2string(fpX.getFractionSignalValue())<<endl;

#else
			cout<< "Exc="<<fpX.getExceptionSignalValue()
				 <<", S="<<fpX.getSignSignalValue()
				 <<", Exp="<<fpX.getExponentSignalValue()-(intpow2(wEX_-1)-1)
				 <<", Frac="<<fpX.getFractionSignalValue()<<endl;
#endif				 
		}
	
		//assign value to sA only at final iteration	
		if (currentIteration==LongAccN)
			sA = sInt2C2(AccValue_, sizeAcc_);
	
 
	}

	mpz_class LongAcc::mapFP2Acc(FPNumber X)
	{
		//get true exponent of X
		mpz_class expX = X.getExponentSignalValue() - ( intpow2(wEX_-1)-1 ); 
	
		int keepBits = -LSBA_ + expX.get_si();
		if (keepBits<0)
			return 0;
		else if (wFX_>keepBits)
			return (X.getFractionSignalValue()/mpz_class(intpow2(wFX_-keepBits)));
		else
			return (X.getFractionSignalValue()*mpz_class(intpow2(keepBits-wFX_)));
	}

	mpz_class LongAcc::sInt2C2(mpz_class X, int width)
	{
		if (X>=0)
			return X;
		else{
			return (intpow2( width)+X);
		}
	}

}
