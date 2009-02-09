#include "Operator.hpp"
#include "TestCase.hpp"
#include "FPNumber.hpp"
#include "utils.hpp"





TestCaseList::TestCaseList() { }
TestCaseList::~TestCaseList() { }

void TestCaseList::add(TestCase* tc){
	v.push_back(tc);
}

int TestCaseList::getNumberOfTestCases(){
	return v.size();
}

TestCase* TestCaseList::getTestCase(int i){
	return v[i];
}



/*
A test case is a mapping between I/O signal names and boolean values given as mpz.

The signal names must be those of Operator->iolist_. Whether several
possible output values are possible is stored in the
numberOfPossibleValues_ attribute of the corresponding Signal stored in iolist, and
only there.

The emulate() function of Operator takes a partial test case (mapping
all the inputs) and completes it by mapping the outputs.

*/


TestCase::TestCase(Operator* op) : op_(op){
}

TestCase::~TestCase() {
}


void TestCase::addInput(string s, mpz_class v)
{
	// TODO Check that the signal exists, then check that the value has the right size
	inputs[s] = v;
}


void TestCase::addInput(string name, double x) {
	// get signal size
	Signal* s = op_->getSignalByName(name);
	if(!s->isFP()) {
		throw "TestCase::addInput: Cannot convert a double into non-FP signal";
	} 
	int wE=s->wE();
	int wF=s->wF();
	FPNumber  fpx(wE, wF);
	fpx=x;
	mpz_class mpx = fpx.getSignalValue();

	inputs[name] = mpx;
}


mpz_class TestCase::getInputValue(string s){
	return inputs[s];
}

void TestCase::addExpectedOutput(string s, mpz_class v)
{
	// TODO Check that the signal exists, then check that the value has the right size and that we haven't added too many possible outputs already
	outputs[s].push_back(v);
}

 
string TestCase::getInputVHDL(string prepend)
{
	ostringstream o;

	/* Add comment if there is one */
	if (comment != "")
		o << prepend << "-- " << comment << endl;

	/* Iterate through input signals */
	for (map<string, mpz_class>::iterator it = inputs.begin(); it != inputs.end(); it++)
	{
		string signame = it->first;
		Signal* s = op_->getSignalByName(signame);
		mpz_class v = it->second;
		o << prepend;
		o << signame << " <= " << s->valueToVHDL(v) << "; ";
		o << endl;
	}

	return o.str();
}



string TestCase::getExpectedOutputVHDL(string prepend)
{
	ostringstream o;

	/* Add comment, if there is one */
	if (comment != "")
		o << prepend << "-- " << comment << endl;

	/* Iterate through output signals */
	for (map<string, vector<mpz_class> >::iterator it = outputs.begin(); it != outputs.end(); it++)
	{
		string signame = it->first;
		Signal* s = op_->getSignalByName(signame);
		vector<mpz_class> vs = it->second;
		string expected;

		o << prepend;
		o << "assert false";  // XXX: Too lazy to make an exception for the first value

		/* Iterate through possible output values */
		for (vector<mpz_class>::iterator it = vs.begin(); it != vs.end(); it++)
		{
			mpz_class v = *it;
			o << " or ";
			if (s->isFP())
				o << "fp_equal(" << s->getName() << ",fp" << s->width() << "'("<< s->valueToVHDL(v) << "))";
			else
				o << s->getName() << "=" << s->valueToVHDL(v);
			expected += " " + s->valueToVHDL(v,false);
		}

		o << " report \"Incorrect output value for " << s->getName() << ", expected" << expected << "\" severity ERROR; ";
		o << endl;
	}

	return o.str();
}



void TestCase::addComment(string c) {
	comment = c;
}




