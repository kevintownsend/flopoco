/*
 * The range reduction box for the FP logarithm
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
#include <math.h>
#include "../utils.hpp"
#include "LogRangeRed.hpp"
#include "../FPLog.hpp"

using namespace std;

extern vector<Operator*> oplist;


/** The constructor when LogRangeRed is instantiated from an FPLog
 * (who knows, other situations may arise) 
*/

LogRangeRed :: LogRangeRed(Target* target,
									FPLog *fplog
									): 
	Operator(target), 
	fplog(fplog),
	stages(fplog->stages),
	a(fplog->a), p(fplog->p), psize(fplog->psize), s(fplog->s)

{
	int i;
	ostringstream name; 
	name << fplog->getName() <<"_RangeRed" ;
	setName(name.str());
	setOperatorType();


	// Now allocate the various table objects


	it0 = new FirstInvTable(target, a[0], a[0]+1);
	oplist.push_back(it0);

	lt0 = new FirstLogTable(target, a[0], fplog->target_prec, it0);
	oplist.push_back(lt0);

	// TODO: check this one is useful
	// it1 = new SecondInvTable(target, fplog->a[1], fplog->p[1]);
	// oplist.push_back(it1);

	for(i=1;  i<= stages; i++) {
		lt[i] = new OtherLogTable(target, a[i], fplog->target_prec - p[i], p[i], i); 
		oplist.push_back(lt[i]);
	 }


	// IO declarations

	// Declare input and output signals
	// o <<   "          Y0 : in  std_logic_vector("<<wF+1<<" downto 0);"<<endl;
	// o <<   "          A  : in std_logic_vector("<< a[0]-1 <<" downto 0);"<<endl;
	//	if(isSequential()) 
	// 	o << "          clk  : in std_logic;"<<endl;
	// o <<   "          Z  : out std_logic_vector("<< s[stages+1]-1 <<" downto 0);"<<endl;
	// o <<   "  almostLog  : out std_logic_vector("<< lt0->wOut-1 <<" downto 0)  );"<<endl;


	addInput("Y0", fplog->wF+2);
	addInput("A", a[0]);
	addOutput("Z", s[stages+1]);
	addOutput("almostLog", lt0->wOut);


	// Signal declarations

	for (i=0; i<= stages; i++) {
	// 		o << "   signal       A"<<i<<":  std_logic_vector("<< a[i] - 1  <<" downto 0);"<<endl;
		ostringstream name; 
		name << "A" << i ;
		addSignal(name.str(), a[i]);
	}
	 	
 	for (i=1; i<= stages; i++) {
		// o << "   signal       B"<<i<<":  std_logic_vector("<< s[i] - a[i] - 1  <<" downto 0);"<<endl;
		ostringstream name; 
		name << "B" << i ;
		addSignal(name.str(), s[i] - a[i]);
	}
	
	for (i=0; i<= stages+1; i++) {
		// o << "   signal Z"<<i<<", Z"<<i<<"_d:  std_logic_vector("<< s[i] - 1  <<" downto 0);"<<endl;
		ostringstream name; 
		name << "Z" << i ;
		addSignal(name.str(), s[i]);
	}
	
	for (i=1; i<= stages; i++) {
		// o << "   signal    epsZ"<<i<<":  std_logic_vector("<< s[i]+p[i]+1  <<" downto 0);"<<endl;
		ostringstream name; 
		name << "epsZ" << i ;
		addSignal(name.str(), s[i] + p[i]);
	}

	// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>eps[i]
	// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
	// We need only this many bits of Z[i] to compute it.
	for (i=1; i<= stages; i++) {
	// 		o << "   signal      ZM"<<i<<":  std_logic_vector("<< psize[i] - 1  <<" downto 0);"<<endl;
		ostringstream name; 
		name << "ZM" << i ;
		addSignal(name.str(),  psize[i]);
	}
	
	// 	o << "   signal       P0:  std_logic_vector("<<  it0->wOut + s[0] - 1  <<" downto 0);"<<endl;
	addSignal("P0",  it0->wOut + s[0]);
	for (i=1; i<= stages; i++) {
		// o << "   signal       P"<<i<<":  std_logic_vector("<< psize[i]+a[i] - 1  <<" downto 0);"<<endl;
		ostringstream name; 
		name << "P" << i ;
		addSignal(name.str(),  psize[i] + a[i]);
	}
	
	// 	o << "   signal       L0:  std_logic_vector("<< lt0->wOut -1 <<" downto 0);"<<endl;
	addSignal("L0",  lt0->wOut);
	
	for (i=1; i<= stages; i++) {
	// 			o << "   signal       L"<<i<<":  std_logic_vector("<< lt[i]->wOut -1 <<" downto 0);"<<endl;
		ostringstream name; 
		name << "L" << i ;
		addSignal(name.str(),  lt[i]->wOut);
	}
	
	// Note that all the Si have the size of L0: absence of carry out is proven.
	for (i=1; i<= stages+1; i++) {
		//o << "   signal S"<<i<<", S"<<i<<"_d:  std_logic_vector("<< lt0->wOut-1 <<" downto 0);"<<endl;
		ostringstream name; 
		name << "S" << i ;
		addSignal(name.str(),  lt0->wOut);
	}

	// 	o << "   signal    InvA0:  std_logic_vector("<< a[0]  <<" downto 0);"<<endl;
	addSignal("invA0",  a[0]+1);

} 


LogRangeRed::~LogRangeRed()
{
	delete it0;
	delete lt0;
	//delete it1;
	for(int i=1; i<=stages; i++) 
		delete lt[i];
}





void LogRangeRed::outputVHDL(std::ostream& o, std::string name)
{
	int i;
	licence(o, "F. de Dinechin, C. Klein  (2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	it0->outputVHDLComponent(o);
	lt0->outputVHDLComponent(o);
	for(i=1;  i<= stages; i++) {
		lt[i]->outputVHDLComponent(o);
	 }
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);

	o << tab << "A0 <= A;"<<endl;
	o << tab << "it0: "<< it0->getName() << " port map (x=>A0, y=>InvA0);" <<endl; 
	o << tab << "lt0: "<< lt0->getName() << " port map (x=>A0, y=>L0);"<<endl;
	o << tab << "P0 <= InvA0 * Y0;" <<endl <<endl;
	o << tab << "Z1 <= P0("<< s[1] -1<<" downto 0);"<<endl;
	o << tab << "S1 <= L0;"<<endl;

	for (i=1; i<= stages; i++) {
		o <<endl;
			//computation
		o << tab << "A"<<i<<" <= Z"<<i<<"(" << s[i] - 1  <<" downto "<< s[i] - a[i]  << ");"<<endl;
		o << tab << "B"<<i<<" <= Z"<<i<<"(" << s[i] - a[i] - 1  <<" downto 0 );"<<endl;
		o << tab << "lt"<<i<<": " << lt[i]->getName() << " port map (x=>A"<<i<<", y=>L"<<i<<");"<<endl;
		if(psize[i] == s[i])
			o << tab << "ZM"<<i<<" <= Z"<<i<< ";"<<endl;   
		else
			o << tab << "ZM"<<i<<" <= Z"<<i<<"(" <<s[i]-1   <<" downto "<< s[i]-psize[i]  << ");"<<endl;   
		o << tab << "P"<<i<<" <= A"<<i<<"*ZM"<<i<<";"<<endl;

		if(i==1) { // special case for the first iteration
			o << tab << "epsZ"<<i<<" <= ("<<s[i]+p[i]+1<<" downto 0 => '0') "
			  << tab << "  when  A1 = ("<<a[1]-1<<" downto 0 => '0')"<<endl
			  << tab << "    else (\"01\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<" )"
			  << "  when ((A1("<<a[1]-1<<")='0') and (A1("<<a[1]-2<<" downto 0) /= ("<<a[1]-2<<" downto 0 => '0')))"<<endl
			  << tab << "    else "
			  << "(\"1\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<"  & \"0\") "
			  << ";"<<endl;
		}
		else {
			o << tab << "epsZ"<<i<<" <=  ("<< s[i]+p[i]+1<<" downto 0 => '0') "
			  << tab << "  when  A"<<i<<" = ("<<a[i]-1<<" downto 0 => '0')"<<endl
			  << tab << "  else    (\"01\" & ("<<p[i]-1<<" downto 0 => '0') & Z"<<i<<");"<<endl;
		}

		o << tab << "Z"<<i+1<<" <=   (\"0\" & B"<<i;
		if (s[i+1] > 1+(s[i]-a[i]))  // need to padd Bi
			o << " & ("<<s[i+1] - 1-(s[i]-a[i]) -1<<" downto 0 => '0') ";    
		o <<")"<<endl 
		  << tab << "      - ( ("<<p[i]-a[i]<<" downto 0 => '0') & P"<<i;
		// either pad, or truncate P
		if(p[i]-a[i]+1  + psize[i]+a[i]  < s[i+1]) // size of leading 0s + size of p 
			o << " & ("<<s[i+1] - (p[i]-a[i]+1  + psize[i]+a[i]) - 1 <<" downto 0 => '0')";  // Pad
		if(p[i]-a[i]+1  + psize[i]+a[i]  > s[i+1]) 
			//truncate
			o <<"("<< psize[i]+a[i] - 1  <<" downto "<<  p[i]-a[i]+1  + psize[i]+a[i]  - s[i+1] << " )";
		o << "  )"<< endl;
		
		o << tab << "      + epsZ"<<i << "("<<s[i]+p[i]+1<<" downto "<<s[i]+p[i] +2 - s[i+1]<<")"
		  << ";"<<endl;
		
		
		o << tab << "S"<<i+1<<" <=   S"<<i<<" + (("<<lt0->wOut-1<<" downto "<<lt[i]->wOut<<" =>'0') & L"<<i<<");"<<endl;
	}
	o << "   Z <= Z"<<stages+1<<";"<<endl;  
	o << "   almostLog <= S"<<stages+1<<";"<<endl;  
	endArchitecture(o);

}
