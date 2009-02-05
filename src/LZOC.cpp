/*
 * A leading zero/one counter for FloPoCo
 *
 * Author : Florent de Dinechin, Bogdan Pasca
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
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "LZOC.hpp"

using namespace std;

LZOC::LZOC(Target* target, int wIn, int wOut) :
	Operator(target), wIn_(wIn), wOut_(wOut)  {

	p2wOut_ = 1<<wOut_; // no need for GMP here

	setOperatorName();
	setOperatorType();
	
	//Set up the IO signals
	addInput ("I", wIn_);
	addInput ("OZB");  
	addOutput("O", wOut_);

	//Set up the internal architecture signals
	if (isSequential()){
			addDelaySignal("sozb",1,wOut_-1);
			
			for (int i=wOut_; i>0; i--){
			ostringstream signalName;
			signalName<<"level"<<i;
			if (i!=1)
				addDelaySignal(signalName.str(), (1<<i));
			else
				addSignal(signalName.str(), (1<<i));
			}
	}
	else{
		addSignal("tmpO", wOut_);
	
		for (int i=wOut_; i>0; i--){
			ostringstream signalName;
			signalName<<"level"<<i;
			addSignal(signalName.str(), (1<<i));
		}
	}
	 
	if (isSequential())
	setPipelineDepth(wOut_-1); 

}

LZOC::~LZOC() {}

void LZOC::setOperatorName(){
	ostringstream name; 
	name <<"LZOC_"<<wIn_<<"_"<<wOut_;
	uniqueName_ = name.str();
}

void LZOC::outputVHDL(std::ostream& o, std::string name) {
	licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	outputVHDLSignalDeclarations(o);	
	
	for (int i=wOut_; i>0; i--){
		o << tab << "signal tmpO"<<i<<"  : std_logic_vector("<<wOut_-1<<" downto "<<i-1<<");"<<endl;
		o << tab << "signal tmpO"<<i<<"_d: std_logic_vector("<<wOut_-1<<" downto "<<i-1<<");"<<endl;				
	}
		
	beginArchitecture(o);
	
	if (isSequential()){
		
		outputVHDLRegisters(o); o<<endl;
		
		o << tab << "process(clk, rst)"<<endl;
		o << tab << "  begin"<<endl;
		o << tab << "if clk'event and clk = '1' then"<<endl;
		o << tab << "     if rst = '1' then"<<endl;
		for (int i=wOut_; i>0; i--)
		o << tab << "        tmpO"<<i<<"_d <= ("<<wOut_-1<<" downto "<<i-1<<" => '0');"<<endl;				
		o << tab << "     else"<<endl;
		for (int i=wOut_; i>0; i--)
		o << tab << "        tmpO"<<i<<"_d <= tmpO"<<i<<";"<<endl;				
		o << tab << "     end if;"<<endl;
		o << tab << "  end if;"<<endl;
		o << tab << "end process;"<<endl;
		
		o<<tab<<"sozb <= ozb;"<<endl;
		
		
		// connect first stage to I
		if (p2wOut_==wIn_)
			o << "  level"<<wOut_<<" <=  I ;" << endl;
		else if (p2wOut_>wIn_) // pad input with zeroes/ones function of what we count. If LZC pad with 1, else pad with 
			o << "  level"<<wOut_<<" <=  I & ("<<p2wOut_-wIn_-1<<" downto 0 => not(sozb)) ;" << endl;
		else if (p2wOut_<wIn_)
			o << "  level"<<wOut_<<" <=  I("<<wIn_-1<<" downto "<<wIn_ - p2wOut_ <<") ;" << endl;
		// recursive structure
		for (int i=wOut_; i>1; i--) {
			if (i!=wOut_)
			o << "  tmpO"<<i<<"("<<wOut_-1<<" downto "<<i<<") <= tmpO"<<i+1<<"_d("<<wOut_-1<<" downto "<<i<<");"<<endl;
			
			o << "  tmpO"<<i<<"("<<i-1<<") <= '1' when level"<<i<<"("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<") = ("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<" => "
											<<delaySignal("sozb",wOut_-i)<<") else '0';"<< endl;
			o << "  level"<<i-1<<" <= level"<<i<<"_d("<<(1<<(i-1))-1<<" downto 0) when tmpO"<<i<<"_d("<<i-1<<") = '1'"<< endl
				<< "               else level"<<i<<"_d("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<");" << endl;
		}
		o << "  tmpO"<<1<<"("<<wOut_-1<<" downto "<<1<<") <= tmpO"<<2<<"_d("<<wOut_-1<<" downto "<<1<<");"<<endl;
		o << "  tmpO1(0) <= '1' when level1(1) =  "<<delaySignal("sozb",wOut_-1)<<" else '0';"<< endl;
				
		o << " O <= tmpO1;"<<endl;
	}
	else{
		// connect first stage to I
		cout<<"p2wOut_ is = "<<p2wOut_<<endl;
		if (p2wOut_==wIn_)
			o << "  level"<<wOut_<<" <=  I ;" << endl;
		else if (p2wOut_>wIn_) // pad input with zeroes/ones function of what we count. If LZC pad with 1, else pad with 
			o << "  level"<<wOut_<<" <=  I & ("<<p2wOut_-wIn_-1<<" downto 0 => not(ozb)) ;" << endl;
		else if (p2wOut_<wIn_)
			o << "  level"<<wOut_<<" <=  I("<<wIn_-1<<" downto "<<wIn_ - p2wOut_ <<") ;" << endl;
		// recursive structure
		for (int i=wOut_; i>1; i--) {
			o << "  tmpO("<<i-1<<") <= '1' when level"<<i<<"("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<") = ("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<" => ozb) else '0';"<< endl;
			o << "  level"<<i-1<<" <= level"<<i<<"("<<(1<<(i-1))-1<<" downto 0) when tmpO("<<i-1<<") = '1'"<< endl
				<< "               else level"<<i<<"("<<(1<<i)-1<<" downto "<<(1<<(i-1))<<");" << endl;
		}
		o << "  tmpO(0) <= '1' when level1(1) =  ozb else '0';"<< endl;
		o << " O <= tmpO;"<<endl;
	}
	
	endArchitecture(o);
}

TestIOMap LZOC::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("I"));
	tim.add(*getSignalByName("OZB"));
	tim.add(*getSignalByName("O"));
	return tim;
}

void LZOC::fillTestCase(mpz_class a[])
{
	mpz_class& si   = a[0];
	mpz_class& sozb = a[1];
	mpz_class& so   = a[2];
	
	int j;
	int bit = (sozb == 0) ? 0 : 1;
	for (j = 0; j < wIn_; j++)
	{
		if (mpz_tstbit(si.get_mpz_t(), wIn_ - j - 1) != bit)
			break;
	}
	so = j;
}

