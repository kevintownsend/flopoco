/*
 * A long accumulator for FloPoCo
 *
 * Author : Florent de Dinechin and Bogdan Pasca
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
#include "utils.hpp"
#include "Operator.hpp"
#include "LongAcc.hpp"

using namespace std;

extern vector<Operator*> oplist;

LongAcc::LongAcc(Target* target, int wEX, int wFX, int MaxMSBX, int LSBA, int MSBA): 
	Operator(target), 
	wEX_(wEX), wFX_(wFX), MaxMSBX_(MaxMSBX), LSBA_(LSBA), MSBA_(MSBA)
{
	int i;

	//check input constraints, i.e, MaxMSBX <= MSBA, LSBA<MaxMSBx
	if ((MaxMSBX_ > MSBA_)||(LSBA_>=MaxMSBX_)){
		cerr << " Input constraints are not met. Please consult FloPoCo help for more details."<<endl;
		exit (EXIT_FAILURE);
	}

	setOperatorName();
	// This operator is a sequential one
	setSequential();

	// Set up various architectural parameters
	sizeAcc_ = MSBA_-LSBA_+1;
		
	maxShift_ = MaxMSBX_-LSBA_; // shift is 0 when the implicit 1 is at LSBA_
	sizeShift_ = intlog2(maxShift_);
	sizeSummand_ = MaxMSBX_-LSBA_+1; 
	sizeShiftedFrac_ = maxShift_ + wFX_+1;
	E0X_ = (1<<(wEX_-1)) -1; // exponent bias
	// The left shift value is MaxMSBX_ - (EX-E0X_) if it is positive
	// If it is negative, it means the summand is shifted out completely: it will be set to 0.
	// So we have to compute (MaxMSBX_+E0X_) - EX
	// The constant should fit on wEX_ bits, because MaxMSBX_ is one possible exponent of X
	// so we have a subtraction of two wEX_-bit numbers.
	int biasedMaxMSBX_ = MaxMSBX_ + E0X_;
	if(biasedMaxMSBX_ < 0 || intlog2(biasedMaxMSBX_)>wEX_) {
		cerr << "ERROR in LongAcc: MaxMSBX_="<<MaxMSBX_<<" is not a valid exponent of X (range " << (-E0X_) << " to " << ((1<<wEX_)-1)-E0X_ <<endl;
		exit (EXIT_FAILURE);
	}

	setOperatorName();	

	// Create an instance of the required input shifter_. 
	shifter_ = new Shifter(target, wFX_+1, maxShift_, Left);
	oplist.push_back(shifter_);

	if(verbose)
		cout << "   shifter_ pipeline depth is "<<shifter_->getPipelineDepth()<<endl;


	addFPInput ("X", wEX_,wFX_);
	addOutput ("A", sizeAcc_);  
	//  add_output ("XOverflow");  
	//  add_output ("XUnderflow");  
	//  add_output ("AccOverflow");  

	// Unregistered signal
	addSignal("exnX", 2);
	addSignal("expX", wEX_);
	addSignal("fracX", wFX_+1);
	addSignal("shifted_frac", sizeShiftedFrac_);
	addSignal("shiftval", wEX_+1); //, "includes a sign bit");
	
	// setup the pipeline 

	if(target->isPipelined()) {
		// TODO here code to pipeline the 2's complement if the target frequency is high, see IntAdder 
		// meanwhile we just do it in 1 level, assuming there is a register at the out of the shifter_.
		c2PipelineDepth_ = 1;
	} 
	else {
		c2PipelineDepth_ = 0;
	}
	
	int suggestedAdditionChunkSize;
	bool status = target->suggestSubaddSize(suggestedAdditionChunkSize ,sizeAcc_);

	if (!status)
		cout << "Warning: the desired frequency is not possible; optimizing to maximum frequency"<<endl;
	
	additionNumberOfChunks_ = (sizeAcc_/suggestedAdditionChunkSize) + (((sizeAcc_%suggestedAdditionChunkSize)==0)?0:1); 
	//rebalance chunks;
	rebalancedAdditionChunkSize_ = (sizeAcc_ / additionNumberOfChunks_)  + (((sizeAcc_%additionNumberOfChunks_)==0)?0:1); 
	rebalancedAdditionLastChunkSize_ = sizeAcc_ - (additionNumberOfChunks_-1)*rebalancedAdditionChunkSize_;

	if (verbose){
		cout << "acumulator size ="<<sizeAcc_<<endl;
		cout << "suggested Addition Chunk Size =  "<<suggestedAdditionChunkSize<<endl; 
		cout << "addition Number Of Chunks =  "<<additionNumberOfChunks_<<endl; 
		cout << "rebalanced Addition Chunk Size =  "<<rebalancedAdditionChunkSize_<<endl; 
		cout << "rebalanced Addition Last Chunk Size =  "<<rebalancedAdditionLastChunkSize_<<endl;
	}
	
	//define the registers which now form the accumulator
	for (i=0;i<additionNumberOfChunks_;i++) {
		ostringstream accReg;
		accReg<<"acc_"<<i;
			if (i==additionNumberOfChunks_-1){
				addRegisteredSignalWithSyncReset(accReg.str(), rebalancedAdditionLastChunkSize_); 
				accReg << "_ext";
				addSignal(accReg.str(), rebalancedAdditionLastChunkSize_ + 1);   
			}
			else{
				addRegisteredSignalWithSyncReset(accReg.str(), rebalancedAdditionChunkSize_);  
				accReg << "_ext";
				addSignal(accReg.str(), rebalancedAdditionChunkSize_ + 1);   
			}
	}
	 
	addRegisteredSignalWithSyncReset("carryIn",1);
	//define the carry propagation registers
	if (additionNumberOfChunks_>1)
		for (i=1;i<additionNumberOfChunks_;i++){
			ostringstream carryBit;
			carryBit<< "carryBit_"<<i;
			addRegisteredSignalWithSyncReset(carryBit.str(), 1);  
		}
			

	// on one side, add delays for the non-complemented signal
	addDelaySignal("summand", sizeSummand_, c2PipelineDepth_); 
	addDelaySignal("flushedToZero", 1, shifter_->getPipelineDepth());
	addDelaySignal("signX", 1, shifter_->getPipelineDepth());

	addSignal("summand2c", sizeSummand_);

	// final pipeline depth is the sum of shifter_ pipeline and 2's
	// complement pipeline, plus 1 for the accumulator itself.
	//TODO when the accumulator is pipelined: replace this 1 with the acc's pipeline depth
	setPipelineDepth(shifter_->getPipelineDepth() + c2PipelineDepth_+  additionNumberOfChunks_ + 1);

	// it is rather stupid to register all the extended bits as they are
	// all equal, but the synthesiser optimises it out
	addRegisteredSignalWithSyncReset("ext_summand2c", sizeAcc_);
	addSignal("acc", sizeAcc_); //,  "includes overflow bit");

	if(verbose)
		cout << tab <<getOperatorName()<< " pipeline depth is " << getPipelineDepth() << " cycles" <<endl;
}

LongAcc::~LongAcc() {
}

void LongAcc::setOperatorName(){
	ostringstream name; 
	name <<"LongAcc_"<<wEX_<<"_"<<wFX_<<"_"
			 <<(MaxMSBX_>=0?"":"M")<<abs(MaxMSBX_)<<"_"
			 <<(LSBA_>=0?"":"M")<<abs(LSBA_)<<"_"
			 <<(MSBA_>=0?"":"M")<<abs(MSBA_) ;
	uniqueName_=name.str();
}


void LongAcc::outputVHDL(ostream& o, string name) {
int i;

	licence(o,"Florent de Dinechin (2007)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	shifter_->outputVHDLComponent(o);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	o << tab << "fracX <= \"1\" & X("<<wFX_-1<<" downto 0);" << endl;
	o << tab << "expX <= X("<<wEX_+wFX_-1<<" downto "<<wFX_<<");" << endl;
	o << tab << "signX <= X("<<wEX_+wFX_<<");" << endl;
	o << tab << "exnX <= X("<<wEX_+wFX_+2<<" downto "<<wEX_+wFX_+1<<");" << endl;
	int64_t exp_offset = E0X_+LSBA_;
	if(exp_offset >=0) {
		o << tab << "shiftval <=  (\"0\"&expX) - \"";  printBinNum(o, exp_offset,  wEX_+1); o << "\";" << endl;
	}
	else {
		o << tab << "shiftval <=  (\"0\"&expX) + \"";  printBinNum(o, -exp_offset,  wEX_+1); o  << "\";" << endl;
	}
	o << endl;
	o << tab << "-- shift of the input into the proper place " << endl;
	o << tab << "input_shifter : " << shifter_->getOperatorName() << endl;
	o << tab << "    port map ( X => fracX, " << endl;
	o << tab << "               S => shiftval("<< shifter_->getShiftAmmount() - 1 <<" downto 0), " << endl;
	o << tab << "               R => shifted_frac";
	if (shifter_->isSequential()) {
		o<<"," << endl;
		o << tab << "             clk => clk," << endl;
		o << tab << "             rst => rst " << endl;
	}
	o << tab << "             );" << endl; 
	o << endl;
 
	o << tab << "flushedToZero <=     '1' when (shiftval("<<wEX_<<")='1' -- negative left shift " << endl;
	o << tab << "                               or exnX=\"00\")" << endl;
	o << tab << "                 else'0';" << endl;
	o << tab << "summand <= ("<<sizeSummand_-1<<" downto 0 => '0')  when "<< getDelaySignalName("flushedToZero", shifter_->getPipelineDepth()) << "='1'  else shifted_frac("<<sizeShiftedFrac_-1<<" downto "<<wFX_<<");" << endl;
	o << endl;

	o << tab << "-- 2's complement of the summand" << endl;
	//Don't compute 2's complement just yet, just invert the bits and leave the addition of the extra 1 in accumulation.
	o << tab << "summand2c <= summand when "<< getDelaySignalName("signX", shifter_->getPipelineDepth()) <<"='0' else not(summand); "<< endl;
	
	o << endl;
	o << tab << "-- extension of the summand to accumulator size" << endl;
	o << tab << "ext_summand2c <= ("<<sizeAcc_-1<<" downto "<<sizeSummand_<<"=> (" 
		<< getDelaySignalName("signX", shifter_->getPipelineDepth()) 
		<< " and not " << getDelaySignalName("flushedToZero", shifter_->getPipelineDepth())
		<<")) & summand2c;" << endl;
	o << endl;
	o << tab << "-- accumulation itself" << endl;
	
	o << tab << "carryIn <= ("<<getDelaySignalName("signX", shifter_->getPipelineDepth())
		<< " and not " << getDelaySignalName("flushedToZero", shifter_->getPipelineDepth() ) << ");"<<endl;  
	
	
	
	for (i=0;i<additionNumberOfChunks_;i++) {
		ostringstream accReg;
		accReg<<"acc_"<<i;
			if ((i==0) & (additionNumberOfChunks_==1) )
				o << tab << "acc_0 <= acc_0_d + ext_summand2c_d + carryIn_d;"<<endl;
			else 
				if ((i==0) & (additionNumberOfChunks_!=1) )
					o << tab << "acc_0 <= acc_0_d + ext_summand2c_d("<<rebalancedAdditionChunkSize_ -1 << " downto  0" <<") + carryIn_d;"<<endl;
				else 
					if (i!=additionNumberOfChunks_-1){
						o<<tab<<"acc_"<<i<<"_ext <= ( \"0\" & acc_"<<i<<"_d ) + ( \"0\" & ext_summand2c_d("<<rebalancedAdditionChunkSize_*(i+1)-1 << " downto "<<rebalancedAdditionChunkSize_*i<<")) + carryBit_"<<i<<"_d;"<<endl;
					
						o<<tab<<"acc_"<<i<<" <= acc_"<<i<<"_ext("<<rebalancedAdditionChunkSize_ -1 <<" downto 0);"<<endl;
						o<<tab<<"carryBit_"<<i<<" <= acc_"<<i<<"_ext("<<rebalancedAdditionChunkSize_<<");"<<endl;
					}  
					else{
						o<<tab<<"acc_"<<i<<"_ext <= ( \"0\" & acc_"<<i<<"_d) + (\"0\" & ext_summand2c_d("<<sizeAcc_-1 << " downto "<<rebalancedAdditionChunkSize_*i<<")) + carryBit_"<<i<<"_d;"<<endl;
						
						o<<tab<<"acc_"<<i<<" <= acc_"<<i<<"_ext("<<rebalancedAdditionLastChunkSize_-1<<" downto 0);"<<endl;
						o<<tab<<"carryBit_"<<i<<" <= acc_"<<i<<"_ext("<<rebalancedAdditionLastChunkSize_<<");"<<endl;
					}
	}
	

//  o << tab << "acc <= ext_summand2c_d   +   acc_d  +  lateAdditionBit;" << endl;
	
		//compose the acc signal 
		
		o << tab << "acc <= ";
		for (i=additionNumberOfChunks_-1;i>=0;i--) {
			
			if (additionNumberOfChunks_==1) 
				o<< "acc_0_d;";
			else
				if (i!=0)
					o<< "acc_"<<i<<"_d & ";
				else
					o<< "acc_"<<i<<"_d;";
		}
	
	o << endl;

	outputVHDLRegisters(o);

	o << tab << "  A <=   acc;" << endl;
	
	//TODO sticky overflow for the accumulator.
	//TODO sticky overflow for the input.
	//TODO sticky underflow for the input.

	endArchitecture(o);
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
		mpfr_random(r); // deprecated function; r is [0,1]
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

}

// read the values from a file and accumulate them
void LongAcc::test_precision2() {

	mpfr_t ref_acc, long_acc, fp_acc, r, d, one, two, msb;
	double dr, dacc, sum, error;

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
