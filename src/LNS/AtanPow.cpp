/*
 * 4/pi * atan(2^x) function with input domain partitioning
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

#include "AtanPow.hpp"
#include <sstream>
#include <vector>
#include "../HOTBM.hpp"

using namespace std;

namespace flopoco{
	extern vector<Operator*> oplist;	// ???

	AtanPow::AtanPow(Target * target, int wE, int wF, int o, EvaluationMethod method) :
		Operator(target), wE(wE), wF(wF), order(o), method_(method)
	{
		T[0] = T[1] = T[2] = 0;

		ostringstream name;
		name<<"AtanPow_"<< wE <<"_"<< wF << "_" << o; 
		uniqueName_ = name.str();
		setCombinatorial();
		addInput ("x", wE + wF);
		addOutput("r", wF);
	
		if(wF > 7) {
			// Who is going to free this memory??
			T[0] = NewInterpolator(wF+wE-1, wF-7, o, -(1 << wE), -8, 1 << 7);
			oplist.push_back(T[0]);

			vhdl << tab << declare("x0", wF+wE-1) << " <= x(" << (wF+wE-2) << " downto 0);\n";
			inPortMap(T[0], "x", "x0");
			outPortMap(T[0], "r","out_t0");
			vhdl << instance(T[0], "inst_t0");	
		}

		vhdl << tab << declare("x1", wF+2) << " <= x(" << (wF+1) << " downto 0);\n";

		if(wF > 6) {
			T[1] = NewInterpolator(wF+2, wF-4, o, -8, -4, 1 << 4);
			oplist.push_back(T[1]);

			inPortMap(T[1], "x", "x1");
			outPortMap(T[1], "r","out_t1");
			vhdl << instance(T[1], "inst_t1");	
		}

		T[2] = NewInterpolator(wF+2, wF, o, -4, 0, 1);
		oplist.push_back(T[2]);

		inPortMap(T[2], "x", "x1");
		outPortMap(T[2], "r","out_t2");
		vhdl << instance(T[2], "inst_t2");	

		vhdl << tab << "r <= ";
	
		if(wF > 7) {
			vhdl << "(" << (wF-1) << " downto " << wF-6 << " => '0') & out_t0(" << (wF-7) << " downto 0)\n"
			     << "         when x(" << (wF+wE-1) << " downto " << (wF+3) << ") /= (" << (wF+wE-1) << " downto " << (wF+3) << " => '1') else\n       ";
		}


		if(wF > 6) {
			vhdl << "(" << (wF-1) << " downto " << (wF-3) << " => '0') & out_t1(" << (wF-4) << " downto 0)\n"
			     << "         when x(" << (wF+2) << ") /= '1' else\n       ";
		}
		vhdl << "out_t2(" << (wF-1) << " downto 0);\n";

	}
	
	Operator * AtanPow::NewInterpolator(int wI, int wO, int o, double xmin, double xmax, double scale)
	{
		char const * my_func = "4/Pi * atan(2^(x))";
		return NewEvaluator(target_, my_func, uniqueName_, wI, wO, o, xmin, xmax, scale, method_);
	}
	

	AtanPow::~AtanPow()
	{
	}

#if 0
	void AtanPow::outputVHDL(std::ostream& o, std::string name)
	{
#if 0 // commented out for use of deprecated methods
		licence(o,"Sylvain Collange (2008)");

		Operator::stdLibs(o);
		outputVHDLEntity(o);
		newArchitecture(o,name);	
		outputVHDLSignalDeclarations(o);	  

		T[2]->outputVHDLComponent(o);

		if(wF > 6) {
			T[1]->outputVHDLComponent(o);
		}

		if(wF > 7) {
			T[0]->outputVHDLComponent(o);
		}

		beginArchitecture(o);

		if(wF > 7) {
			o <<
			"  inst_t0 : " << T[0]->getName() << endl <<
			"    port map ( x => x(" << (wF+wE-2) << " downto 0),\n"
			"               r => out_t0 );\n";
		}

		if(wF > 6) {
			o <<
			"\n"
			"  inst_t1 : " << T[1]->getName() << endl <<
			"    port map ( x => x(" << (wF+1) << " downto 0),\n"
			"               r => out_t1 );\n";
		}

		o <<
		"\n"
		"  inst_t2 : " << T[2]->getName() << endl <<
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
#endif
	}
#endif
}
