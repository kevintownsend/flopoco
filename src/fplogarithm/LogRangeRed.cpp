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

	setCopyrightString("F. de Dinechin (2008)");

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



	addInput("Y0", fplog->wF+2);
	addInput("A", a[0]);
	addOutput("Z", s[stages+1]);
	addOutput("almostLog", lt0->wOut);


	vhdl << tab << declare("A0", a[0]) << " <= A;"<<endl;

	vhdl << tab << "-- First inv table" << endl;
	inPortMap       (it0, "X", "A0");
	outPortMap      (it0, "Y", "InvA0");
	vhdl << instance(it0, "itO");

	vhdl << tab << "-- First log table" << endl;
	inPortMap       (lt0, "X", "A0");
	outPortMap      (lt0, "Y", "L0");
	vhdl << instance(lt0, "ltO");

	vhdl << tab << declare("P0",  it0->wOut + s[0]) << " <= InvA0 * Y0;" <<endl <<endl;
	vhdl << tab << declare("Z1", s[1]) << " <= P0("<< s[1] -1<<" downto 0);"<<endl;
	vhdl << tab << declare("S1", lt0->wOut) << " <= L0;"<<endl;

	for (i=1; i<= stages; i++) {
		vhdl <<endl;
		//computation
		vhdl << tab << declare(join("A", i), a[i]) <<     " <= " << join("Z",i) << range(s[i]-1,      s[i]-a[i]) << ";"<<endl;
		vhdl << tab << declare(join("B",i), s[i]-a[i]) << " <= " << join("Z",i) << range(s[i]-a[i]-1, 0        ) << ";"<<endl;
		inPortMap       (lt[i], "X", join("A", i));
		outPortMap      (lt[i], "Y", join("L", i));
		vhdl << instance(lt[i], join("lt",i));
		vhdl << tab << declare(join("ZM",i), psize[i]) << " <= ";
		if(psize[i] == s[i])
			vhdl << join("Z",i) << ";"<<endl;   
		else
			vhdl << join("Z",i) << range(s[i]-1, s[i]-psize[i])  << ";" << endl;   
		vhdl << tab << declare(join("P",i),  psize[i] + a[i]) << " <= " << join("A",i) << "*" << join("ZM",i) << ";" << endl;

		vhdl << tab << declare(join("epsZ",i), s[i] + p[i] +2 ) << " <= " ;
		if(i==1) { // special case for the first iteration
			vhdl << 	rangeAssign(s[i]+p[i]+1,  0,  "'0'")  << " when  A1 = " << rangeAssign(a[1]-1, 0,  "'0'") << endl
				  << tab << "    else (\"01\" & "<< rangeAssign(p[i]-1, 0, "'0'") << " & Z"<<i<<" )"
				  << "  when ((A1("<< a[1]-1 << ")='0') and (A1" << range(a[1]-2, 0) << " /= " << rangeAssign(a[1]-2, 0, "'0'") << "))" << endl
			  << tab << "    else " << "(\"1\" & " << rangeAssign(p[i]-1, 0, "'0'") << " & " << join("Z",i) << "  & \"0\") "
			  << ";"<<endl;
		}
		else {
			vhdl << rangeAssign(s[i]+p[i]+1,  0,  "'0'") 
				  << tab << "  when " << join("A",i) << " = " << rangeAssign(a[i]-1,  0,  "'0'") << endl
				  << tab << "  else    (\"01\" & " << rangeAssign(p[i]-1,  0,  "'0'") << " & " << join("Z",i) <<");"<<endl;
		}

		vhdl << tab << declare(join("Z", i+1), s[i+1]) << " <=   (\"0\" & " << join("B",i);
		if (s[i+1] > 1+(s[i]-a[i]))  // need to padd Bi
			vhdl << " & " << rangeAssign(s[i+1] - 1-(s[i]-a[i]) -1,  0 , "'0'");    
		vhdl <<")"<<endl 
			  << tab << "      - ( " << rangeAssign(p[i]-a[i],  0,  "'0'") << " & " << join("P", i);
		// either pad, or truncate P
		if(p[i]-a[i]+1  + psize[i]+a[i]  < s[i+1]) // size of leading 0s + size of p 
			vhdl << " & ("<<s[i+1] - (p[i]-a[i]+1  + psize[i]+a[i]) - 1 <<" downto 0 => '0')";  // Pad
		if(p[i]-a[i]+1  + psize[i]+a[i]  > s[i+1]) 
			//truncate
			vhdl <<"("<< psize[i]+a[i] - 1  <<" downto "<<  p[i]-a[i]+1  + psize[i]+a[i]  - s[i+1] << " )";
		vhdl << "  )"<< endl;
		
		vhdl << tab << "      + " << join("epsZ",i) << range(s[i]+p[i]+1,  s[i]+p[i] +2 - s[i+1]) << ";"<<endl;
		
		vhdl << tab << declare(join("S",i+1),  lt0->wOut) << " <= " 
			  <<  join("S",i)  << " + (" << rangeAssign(lt0->wOut-1, lt[i]->wOut,  "'0'") << " & " << join("L",i) <<");"<<endl;
	}
	vhdl << "   Z <= " << join("Z", stages+1) << ";" << endl;  
	vhdl << "   almostLog <= " << join("S",stages+1) << ";" << endl;  

} 


LogRangeRed::~LogRangeRed()
{
#if 0 // better done on oplist 
	delete it0;
	delete lt0;
	//delete it1;
	for(int i=1; i<=stages; i++) 
		delete lt[i];
#endif
}




