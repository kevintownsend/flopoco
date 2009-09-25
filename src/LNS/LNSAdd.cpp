/*
 * LNS Summation operator : log2(1+2^x)
 *
 * Author : Sylvain Collange
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

#include "LNSAdd.hpp"
#include <sstream>
#include <vector>
#include "../HOTBM.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;	// ???

	LNSAdd::LNSAdd(Target * target, int wE, int wF, int o) :
		Operator(target), wE(wE), wF(wF), order(o)
	{
		t[0] = t[1] = t[2] = 0;

		ostringstream name;
		name<<"LNSAdd_"<< wE <<"_"<< wF; 
		uniqueName_ = name.str();
	
		addInput ("x", wE + wF);
		addOutput("r", wF, 2); // Faithful rounding;
	
		if(wF > 7) {
			addSignal("out_t0", wF - 6);
		
			// Who is going to free this memory??
			// T0 always order 1?...
			t[0] = new HOTBM(target, "log2(1+2^x)", uniqueName_, wF+wE-1, wF-7, o, -(1 << wE), -8, 1 << 7);
			oplist.push_back(t[0]);
		}

		if(wF > 6) {
			addSignal("out_t1", wF - 3);
		
			t[1] = new HOTBM(target, "log2(1+2^x)", uniqueName_, wF+2, wF-4, o, -8, -4, 1 << 4);
			oplist.push_back(t[1]);
		}
		addSignal("out_t2", wF + 1);
	
		t[2] = new HOTBM(target, "log2(1+2^x)", uniqueName_, wF+2, wF, o, -4, 0, 1);
		oplist.push_back(t[2]);

	
	}

	LNSAdd::~LNSAdd()
	{
	}

	void LNSAdd::outputVHDL(std::ostream& o, std::string name)
	{
		licence(o,"Sylvain Collange (2008)");

		//Operator::stdLibs(o);
		o<<"library ieee;\nuse ieee.std_logic_1164.all;"<<endl;
	
		outputVHDLEntity(o);
		newArchitecture(o,name);	
		outputVHDLSignalDeclarations(o);	  

		t[2]->outputVHDLComponent(o);

		if(wF > 6) {
			t[1]->outputVHDLComponent(o);
		}

		if(wF > 7) {
			t[0]->outputVHDLComponent(o);
		}

		beginArchitecture(o);

		if(wF > 7) {
			o <<
				"  inst_t0 : " << t[0]->getName() << endl <<
				"    port map ( x => x(" << (wF+wE-2) << " downto 0),\n"
				"               r => out_t0 );\n";
		}

		if(wF > 6) {
			o <<
				"\n"
				"  inst_t1 : " << t[1]->getName() << endl <<
				"    port map ( x => x(" << (wF+1) << " downto 0),\n"
				"               r => out_t1 );\n";
		}

		o <<
			"\n"
			"  inst_t2 : " << t[2]->getName() << endl <<
			"    port map ( x => x(" << (wF+1) << " downto 0),\n"
			"               r => out_t2 );\n"
			"\n";

	
		o << "  r <= ";
	
		if(wF > 7) {
			o << "(" << (wF-1) << " downto " << wF-6 << " => '0') & out_t0(" << (wF-7) << " downto 0)\n"
			  << "         when x(" << (wF+wE-1) << " downto " << (wF+3) << ") /= (" << (wF+wE-1) << " downto " << (wF+3) << " => '1') else\n       ";
		}


		if(wF > 6) {
			o << "(" << (wF-1) << " downto " << (wF-3) << " => '0') & out_t1(" << (wF-4) << " downto 0)\n"
			  << "         when x(" << (wF+2) << ") /= '1' else\n       ";
		}
		o << "out_t2(" << (wF-1) << " downto 0);\n";
	
		o<< "end architecture;" << endl << endl;
	}


#if 0
	void LNSAdd::fillTestCase(mpz_class a[])
	{
		/* Get inputs / outputs */
		mpz_class &x  = a[0];
		mpz_class &rd = a[1]; // rounded down
		mpz_class &ru = a[2]; // rounded up

		int outSign = 0;

		mpfr_t mpX, mpR;
		mpfr_inits(mpX, mpR, 0);

		/* Convert a random signal to an mpfr_t in [0,1[ */
		mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
		mpfr_div_2si(mpX, mpX, wI, GMP_RNDN);

		/* Compute the function */
		//f.eval(mpR, mpX);

		/* Compute the signal value */
		if (mpfr_signbit(mpR))
			{
				outSign = 1;
				mpfr_abs(mpR, mpR, GMP_RNDN);
			}
		mpfr_mul_2si(mpR, mpR, wO, GMP_RNDN);

		/* NOT A TYPO. HOTBM only guarantees faithful
		 * rounding, so we will round down here,
		 * add both the upper and lower neighbor.
		 */
		mpfr_get_z(rd.get_mpz_t(), mpR, GMP_RNDD);
		ru = rd + 1;

	}

#endif
}
