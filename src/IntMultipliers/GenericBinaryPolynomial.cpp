#include <iostream>
#include <sstream>
#include <list>
#include <vector>

#include "gmp.h"
#include "mpfr.h"

#include "GenericBinaryPolynomial.hpp"
#include "FormalBinaryProduct.hpp"
#include "IntCompressorTree.hpp"

//std::string GenericBinaryPolynomial::operatorInfo = "UserDefinedInfo param0 param1 <options>";

using namespace flopoco;

static void print_vhdl_string_of_monomial_option
	(const Option<MonomialOfBits>& o)
{
	if (o.is_empty()) {
		vhdl << "0";
		return;
	}
	MonomialOfBits m = o.get_value();
	bool cont = false;
	size_t i;
	for (i = 0; i < m.data.size(); i++) {
		if (m.data[i]) {
			if (cont)
				vhdl << " and ";
			vhdl << "X" << of (i);
			cont = true;
		}
	}
	if (!cont)
		vhdl << "1";
}

GenericBinaryPolynomial::GenericBinaryPolynomial(Target* target,
                                                 Product p,
						 std::map<std::string,double>
						 	inputDelays)
	:Operator(target,inputDelays), p(p) {

	// don't support zero-sized p
	if (p.data.size())
		throw "Zero-size p not supported";

	ostringstream name;
	name << "GenericBinaryPolynomial_" << p.mon_size << "_" << p.data.size()
	     << "_uid" << Operator::getNewUId();
	setName(name.str());
	setCopyrightString("Guillaume Sergent 2012");

	addInput ("X" , p.mon_size);
	addOutput("R" , p.data.size());

	// here we split (arbitrarily, as of now) p into a sum of vectors
	// that can be easily calculated just by logic ANDs
	// TODO: intelligently sort the list, if necessary
	// (FormalBinaryProduct already implicitly doing sort)
	Product p_copy(p);
	std::list<std::vector<Option<MonomialOfBits> > > p_split;
	for (;;) {
		bool exit_the_loop = true;
		std::vector<Option<MonomialOfBits> > v
			(p.data.size(),Option<MonomialOfBits>());
		std::vector<ProductBit>::iterator it;
		std::vector<Option<MonomialOfBits> >::iterator vit;
		it = p.data.begin(); vit = v.begin();
		// v and p.data have same size
		for (; it != p.data.end(); it++,vit++) {
			if (it->empty()) {
				*vit = Option();
			} else {
				*vit = Option<MonomialOfBits>(it->front());
				it->pop_front();
				exit_the_loop = false;
			}
		}
		// the loop terminates since lists have finite length
		if (exit_the_loop)
			break;
		else
			p_split.push_back (v);
	}
	size_t int_adder_N = 0;
	std::list<std::vector<Option<MonomialOfBits> > >::const_iterator pspit;
	for (pspit = p_split.begin(); pspit != p_split.end(); pspit++) {
		vhdl << tab
		     << declare(join("multiadder_input_",int_adder_N),
				p.data.size())
		     << " <= ";
		bool cont = false;
		std::vector<Option<MonomialOfBits> >::const_iterator vit;
		for (vit = pspit->begin(); vit != pspit->end(); vit++) {
			if (cont)
				vhdl << " & ";
			print_vhdl_string_of_monomial_option (*vit);
		}
		vhdl << ";\n";
		// at the end of the loop int_adder_N will contain
		// the length of the list
		int_adder_N++;
	}
	ima = new IntMultiAdder (target, p.data.size(), int_adder_N,
	                         emptyDelays, false);
	outPortMap (ima, "R", "R");
	for (size_t i = int_adder_N - 1; i >= 0; i++) {
		inPortMap (ima, join("X",i), join("multiadder_input_",i));
	}
	vhdl << instance (ima, "ima");
};

	
void GenericBinaryPolynomial::emulate(TestCase * tc) {
}

void GenericBinaryPolynomial::buildStandardTestCases(TestCaseList * tcl) {
}

void GenericBinaryPolynomial::buildRandomTestCases(TestCaseList *  tcl, int n) {
}

TestCase* GenericBinaryPolynomial::buildRandomTestCases(int i) {
  TestCase* tc = new TestCase(this);
  return tc;
}
