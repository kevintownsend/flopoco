// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
#include "mpfr.h"
#include <gmpxx.h>

// include the header of the Operator
#include "NewCompressorTree.hpp"
#include "PopCount.hpp"

using namespace std;


// personalized parameter
//string NewCompressorTree::operatorInfo = "UserDefinedInfo param0 param1 <options>";


// the length of the result is the same as the one of vops
// the result is truncated on overflow
NewCompressorTree::NewCompressorTree(Target * target, vector<unsigned> vops)
	:Operator(target), w(vops.size()), vert_operands(vops)
{
	{
		// definition of the name of the operator
		ostringstream name;
		name << "NewCompressorTree_" << w;
		vector<unsigned>::const_iterator it = vops.begin();
		for (; it != vops.end(); it++) {
			name << '_' << *it;
		}
		setName(name.str());
	}
	setCopyrightString("Guillaume Sergent 2012");

	for (unsigned i = 0; i < w; i++) {
		addInput (join("X_", i), vops[i]);
	}
	addOutput ("R", w);

	unsigned level = 0, max_height;
	unsigned n = target->lutInputs();
	vector<Operator*> popcounts (n+1, (Operator*) 0);
	for (;;) {
		vector<unsigned>::iterator it;
		for (it = vops.begin(); it != vops.end(); it++) {
			if (max_height < *it)
				max_height = *it;
		}
		bool exit_the_loop = false; // break will only exit the switch
		switch (max_height) {
		case 0:
			if (w)
				vhdl << "R <= \"0\";";
			exit_the_loop = true;
			break;
		case 1:
			// nothing to add
			vhdl << "R <= ";
			for (unsigned i = 0; i < w; i++) {
				if (vops[i]) {
					if (i)
						vhdl << " & ";
					vhdl << "X_" << i << "_level" << level
					     << of(0);
				} else {
					vhdl << "\"0\"";
				}
			}
			vhdl << ";" << endl;
			exit_the_loop = true;
			break;
		case 2:
			// final (binary) addition in the general case
			vhdl << "R <= ";
			for (unsigned i = 0; i < w; i++) {
				if (vops[i]) {
					if (i)
						vhdl << " & ";
					vhdl << "X_" << i << "_level" << level
					     << of(0);
				} else {
					vhdl << "\"0\"";
				}
			}
			vhdl << " + ";
			for (unsigned i = 0; i < w; i++) {
				if (vops[i]) {
					if (i)
						vhdl << " & ";
					vhdl << "X_" << i << "_level" << level
					     << of(1);
				} else {
					vhdl << "\"0\"";
				}
			}
			vhdl << ";" << endl;
			exit_the_loop = true;
			break;
		default:
			vector<unsigned> vops_new (w, 0);
			for (unsigned i = 0; i < w; i++) {
				unsigned inputs = vops[i];
				while (inputs > n) {
					if (!popcounts[n])
						popcounts[n] = new PopCount
							(target, n);
					ostringstream in;
					in << "X_" << i << "_level"
					   << level << "_" << vops[i]-1
					   << "_" << vops[i]-n;
					string out = in.str() + "_popcnt";

					vhdl << declare (in.str(), n)
					     << " <= X_" << i << "_level"
					     << level
					     << range(vops[i]-1,vops[i]-n)
					     << ";" << endl;
					outPortMap (popcounts[n],"R",out);
					inPortMap(popcounts[n],"X",in.str());
					for (int j = 0; j < intlog2(n); j++) {
						if (i+j >= w)
							break;
						ostringstream bit;
						bit << "X_" << i << "_level"
						    << (level+1) << "_bit"
						    << vops_new[i+j];
						vhdl << declare (bit.str())
						     << " <= " << out
						     << of(j) << ";\n";
						vops_new[i+j]++;
					}
					vops[i] -= n;
					inputs = vops[i];
				}
				// if it's compressible
				if (inputs > 2) {
					if (!popcounts[inputs])
						popcounts[inputs] = new
							PopCount (target,
							          inputs);
					ostringstream in;
					in << "X_" << i << "_level"
					   << level << "_" << vops[i]-1 << "_0";
					string out = in.str() + "_popcnt";

					vhdl << declare (in.str(), inputs)
					     << " <= X_" << i << "_level"
					     << level
					     // vops[i] == inputs
					     << range(vops[i]-1, 0)
					     << ";" << endl;
					outPortMap (popcounts[inputs],"R",out);
					inPortMap (popcounts[inputs],"X",
					           in.str());
					for (int j = 0;
					         j < intlog2(inputs);
						 j++) {
						if (i+j >= w)
							break;
						ostringstream bit;
						bit << "X_" << i << "_level"
						    << (level+1) << "_bit"
						    << vops_new[i+j];
						vhdl << declare (bit.str())
						     << " <= " << out
						     << of(j) << ";\n";
						vops_new[i+j]++;
					}
				} else {
					for (unsigned j = 0; j < inputs; j++) {
						ostringstream bit;
						bit << "X_" << i << "_level"
						    << (level+1) << "_bit"
						    << vops_new[i];
						vhdl << declare (bit.str())
						     << " <= X_" << i
						     << "_level" << level
						     << of(i) << ";" << endl;
						vops_new[i]++;
					}
				}
			}
			// and now we constructed every bit of the next level
			// so we are ready to construct the real new inputs
			for (unsigned i = 0; i < w; i++) {
				ostringstream next_lvl;
				next_lvl << "X_" << i << "_level" << (level+1);
				vhdl << declare (next_lvl.str(),vops_new[i])
				     << " <= ";
				// the order of filling of the level vectors
				// doesn't matter: here we fill in reverse
				for (unsigned j = 0; j < vops_new[i]; j++) {
					if (j)
						vhdl << " & ";
					vhdl << "X_" << i << "_level"
					     << (level+1) << "_bit"
					     << j;
				}
				vhdl << ";\n";
			}
			level++;
			vops = vops_new;
			break;
		}
	if (exit_the_loop)
		break;
	}
};

void NewCompressorTree::emulate(TestCase * tc)
{
	std::vector<mpz_class> signals (w, mpz_class(0));
	for (unsigned i = 0; i < w; i++) {
		signals[i] = tc->getInputValue(join("X_",i));
	}
	mpz_class res;
	for (unsigned i = 0; i < w; i++) {
		res += (popcnt (signals[i]) << i);
	}
	tc->addExpectedOutput("R", res);
}


void NewCompressorTree::buildStandardTestCases(TestCaseList * tcl)
{
}

void NewCompressorTree::buildRandomTestCases(TestCaseList * tcl, int n)
{
}

TestCase *NewCompressorTree::buildRandomTestCases(int i)
{
	TestCase *tc = new TestCase(this);
	return tc;
}
