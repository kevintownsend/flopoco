/*
 * A post-normalization unit for the FloPoCo Long Accumulator 
 *
 * Author : Bogdan Pasca
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
#include "LongAcc2FP.hpp"

using namespace std;

extern vector<Operator*> oplist;

LongAcc2FP::LongAcc2FP(Target* target, int MaxMSBX, int LSBA, int MSBA, int wE_out, int wF_out): 
	Operator(target), 
	MaxMSBX(MaxMSBX), LSBA(LSBA), MSBA(MSBA), wE_out(wE_out), wF_out(wF_out)
{
	ostringstream name; 
	int i;
	
	name <<"LongAcc2FP_"
			 <<(MaxMSBX>=0?"":"M")<<abs(MaxMSBX)<<"_"
			 <<(LSBA>=0?"":"M")<<abs(LSBA)<<"_"
			 <<(MSBA>=0?"":"M")<<abs(MSBA)<<"_"
			 <<wE_out<<"_"<<wF_out;
	unique_name=name.str();


	// Set up various architectural parameters
	sizeAcc = MSBA-LSBA+1;	
	
	//instantiate a leading zero/one counter
	wOutLZOC = log2(sizeAcc)+1;
	leadZOCounter = new LZOC(target, sizeAcc, wOutLZOC);
	oplist.push_back(leadZOCounter);

	//instantiate a left shifter
	leftShifter = new Shifter(target,sizeAcc,sizeAcc,Left);
	oplist.push_back(leftShifter);

	// This operator is a sequential one
	set_sequential();
	set_pipeline_depth(leadZOCounter->pipeline_depth() + leftShifter->pipeline_depth());


	//compute the bias value
	expBias = (1<<(wE_out-1)) -1; 

	//inputs and outputs
	add_input  ("A", sizeAcc);
	add_output ("R", 3 + wE_out + wF_out);

	
	add_delay_signal_no_reset("resultSign0",1,leadZOCounter->pipeline_depth() + leftShifter->pipeline_depth());	
	add_delay_signal_bus_no_reset("pipeA",sizeAcc,leadZOCounter->pipeline_depth());
	add_delay_signal_no_reset("expRAdjusted",wE_out,leftShifter->pipeline_depth());
	add_delay_signal_no_reset("excBits",2, leftShifter->pipeline_depth());
	
	
	add_signal("nZO",wOutLZOC);
	add_signal("c2nZO",wOutLZOC+1);
	add_signal("expTest",1);
	add_signal("expAdjValue",wE_out);
	add_signal("excRes",2);
	add_signal("expR",wE_out);
	
	add_signal("resFrac",sizeAcc + sizeAcc);
	
	add_signal("posExpAdj",wOutLZOC+1);
	add_signal("negExpAdj",wOutLZOC+1);
	add_signal("expAdj",wOutLZOC+1);
	
	
	
	add_signal("expAdjustedExt",wOutLZOC+1);
	add_signal("signExpAdjustedExt",1);
	add_signal("modulusExpAdjusted",wOutLZOC);
	add_signal("maxExponent",wOutLZOC);
	add_signal("expOverflow",1);
	add_signal("expUnderflow",1);
}





LongAcc2FP::~LongAcc2FP() {
}


void LongAcc2FP::output_vhdl(ostream& o, string name) {
int i;
	Licence(o,"Bogdan Pasca (2008)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	new_architecture(o,name);
	leadZOCounter->output_vhdl_component(o);
	leftShifter->output_vhdl_component(o);
	output_vhdl_signal_declarations(o);
	begin_architecture(o);
	output_vhdl_registers(o); o<<endl;

	//the sign of the result 
	o<<tab<<"resultSign0 <= A("<<sizeAcc-1<<");"<<endl;
	
	o<<tab<<"pipeA <= A;"<<endl;
	
	//count the number of zeros/ones in order to determine the size of the exponent
	//Leading Zero/One counter 
	o<<tab<< "LZOC_component: " << leadZOCounter->unique_name << endl;
	o<<tab<< "      port map ( I => A, "                      << endl; 
	o<<tab<< "                 OZB => resultSign0, "          << endl; 
	o<<tab<< "                 O => nZO, "                    << endl; 
	o<<tab<< "                 clk => clk, "                  << endl;
	o<<tab<< "                 rst => rst "                   << endl;
	o<<tab<< "               );"                       << endl<< endl;		

	//readjust exponent
	
	//compute the 2's complement value of the number of leading Z/O
	o<<tab<<"c2nZO <= ("<<wOutLZOC<<" downto 0 => '1') - (\"0\" & nZO) + '1';"<<endl;
	
	//determine if the value to be added to the exponent bias is negative or positive.
	o<<tab<<"expTest <= '1' when (CONV_STD_LOGIC_VECTOR("<<MSBA<<","<<wOutLZOC<<")>nZO) else '0';"<<endl;
	
	//compute the exponent operations in both cases => may need pipelining but are generally short additions	
	//o<<tab<<"posExpAdj <= CONV_STD_LOGIC_VECTOR("<<MSBA<<","<<wOutLZOC+1<<")- (\"0\" & nZO); "<<endl;
	o<<tab<<"negExpAdj <= CONV_STD_LOGIC_VECTOR("<<MSBA<<","<<wOutLZOC+1<<") + c2nZO;"<<endl;
	
	//select the actual value which will be added to the expoent bias
	//o<<tab<<"expAdj <= posExpAdj when expTest='1' else"<<endl;
	//o<<tab<<"          negExpAdj;"<<endl;					
	
	o<<tab<<"expAdj <= negExpAdj;"<<endl;
		
	//the exponent bias	
	o<<tab<<"expR <= CONV_STD_LOGIC_VECTOR("<<expBias<<","<<wE_out<<");"<<endl;
	
	
	if (wOutLZOC+1 < wE_out){
		o<<tab<<"expRAdjusted <= expR + (("<<wE_out-1<<" downto "<<wOutLZOC+1<<" => expAdj("<<wOutLZOC<<")) & expAdj);"<<endl;
		o<<tab<<"excBits <=\"01\";"<<endl;
	}
	else
		if (wOutLZOC+1 == wE_out){
			o<<tab<<"expRAdjusted <= expR + expAdj;"<<endl;
			o<<tab<<"excBits <=\"01\";"<<endl;
		}	
		else
		{
			o<<tab<<"expAdjustedExt <= (("<<wOutLZOC<<" downto "<<wE_out<<" =>'0') & expR) + expAdj;"<<endl;
			
			o<<tab<<"signExpAdjustedExt <= expAdjustedExt("<<wOutLZOC<<");"<<endl;
			o<<tab<<"modulusExpAdjusted <= (("<<wOutLZOC-1<<" downto 0 => signExpAdjustedExt) xor expAdjustedExt("<<wOutLZOC-1<<" downto 0)) + signExpAdjustedExt;"<<endl;
				
			o<<tab<<"maxExponent <= CONV_STD_LOGIC_VECTOR("<<(1<<wE_out)-1<<" , "<<wOutLZOC<<");"<<endl;
			
			o<<tab<<"expOverflow <= '1' when (modulusExpAdjusted>maxExponent) and (signExpAdjustedExt='0') else '0';"<<endl;
			o<<tab<<"expUnderflow <= '1' when (modulusExpAdjusted>maxExponent) and (signExpAdjustedExt='1') else '0';"<<endl;
			
			
			o<<tab<<"excBits(1) <= expOverflow and  not(expUnderflow);"<<endl;
			o<<tab<<"excBits(0) <= not(expOverflow) and not(expUnderflow);"<<endl;
					
			o<<tab<<"expRAdjusted<=expAdjustedExt("<<wE_out-1<<" downto 0);"<<endl;
		
		}			
	
		
	// shift 
	o<<tab<< "left_shifter_component: " << leftShifter->unique_name << endl;
	o<<tab<< "      port map ( X => "<<get_delay_signal_name("pipeA",leadZOCounter->pipeline_depth())<<", "<< endl; 
	o<<tab<< "                 S => nZO, " << endl; 
	o<<tab<< "                 R => resFrac, " <<endl; 
	o<<tab<< "                 clk => clk, " << endl;
	o<<tab<< "                 rst => rst " << endl;
	o<<tab<< "               );" << endl<<endl;		
		
	o<<tab<<"excRes <= "<<get_delay_signal_name("excBits", leftShifter->pipeline_depth() ) <<";"<<endl;	
	o<<tab<< "R <= excRes & "
			 <<get_delay_signal_name("resultSign0",leadZOCounter->pipeline_depth())<<" & "
			 <<get_delay_signal_name("expRAdjusted",leftShifter->pipeline_depth())<<"("<<wE_out-1<<" downto 0) & ";

	if (sizeAcc-1-wF_out>=0)		 
			o <<"resFrac("<<sizeAcc-2<<" downto "<<sizeAcc-1-wF_out<<");"<<endl;
	else
			o <<"resFrac("<<sizeAcc-2<<" downto "<<0<<")& CONV_STD_LOGIC_VECTOR(0,"<<-(sizeAcc-1-wF_out)<<");"<<endl;
	
	end_architecture(o);
		
}



