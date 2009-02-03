#include "TestCase.hpp"
#include "utils.hpp"

std::string TestCase::signalValueToVHDL(const Signal& s, mpz_class v, bool quot)
{
	std::string r;

	/* Get base 2 representation */

	r = v.get_str(2);

	/* Some check */
	if (r.size() > s.width())
	{
		std::ostringstream o;
		o << "Error in " <<  __FILE__ << "@" << __LINE__ << ": value (" << r << ") is larger than signal " << s.getName();
		throw o.str();
	}

	/* Do padding */
	while (r.size() < s.width())
		r = "0" + r;

	/* Put apostrophe / quot */
	if (!quot) return r;
	if (s.width() > 1)
		return "\"" + r + "\"";
	else
		return "'" + r + "'";
}

std::string TestCase::signalValueToVHDLHex(const Signal& s, mpz_class v, bool quot)
{
	std::string o;

	/* Get base 16 representation */
	o = v.get_str(16);

	/* Some check */
	/* XXX: Too permissive */
	if (o.size() * 4 > s.width() + 4)
	{
		std::ostringstream o;
		o << "Error in " <<  __FILE__ << "@" << __LINE__ << ": value is larger than signal " << s.getName();
		throw o.str();
	}

	/* Do padding */
	while (o.size() * 4 < s.width())
		o = "0" + o;

	/* Put apostrophe / quot */
	if (!quot) return o;
	if (s.width() > 1)
		return "x\"" + o + "\"";
	else
		return "'" + o + "'";
}

TestCase::TestCase() { }
TestCase::~TestCase() { }

void TestCase::addInput(const Signal& s, mpz_class v)
{
	inputs[s] = v;
}

void TestCase::addExpectedOutput(const Signal& s, mpz_class v)
{
	outputs[s].push_back(v);
}

std::string TestCase::getInputVHDL(std::string prepend)
{
	std::ostringstream o;

	/* Add comment if there is one */
	if (comment != "")
		o << prepend << "-- " << comment << std::endl;

	/* Iterate through input signals */
	for (Inputs::iterator it = inputs.begin(); it != inputs.end(); it++)
	{
		Signal s = it->first;
		mpz_class v = it->second;
		o << prepend;
		o << s.getName() << " <= " << signalValueToVHDL(s, v) << "; ";
		o << std::endl;
	}

	return o.str();
}

std::string TestCase::getExpectedOutputVHDL(std::string prepend)
{
	std::ostringstream o;

	/* Add comment, if there is one */
	if (comment != "")
		o << prepend << "-- " << comment << std::endl;

	/* Iterate through output signals */
	for (Outputs::iterator it = outputs.begin(); it != outputs.end(); it++)
	{
		Signal s = it->first;
		std::vector<mpz_class> vs = it->second;
		std::string expected;

		o << prepend;
		o << "assert false";  // XXX: Too lazy to make an exception for the first value

		/* Iterate through possible output values */
		for (std::vector<mpz_class>::iterator it = vs.begin(); it != vs.end(); it++)
		{
			mpz_class v = *it;
			o << " or ";
			if (s.isFP())
				o << "fp_equal(" << s.getName() << ",fp" << s.width() << "'("<< signalValueToVHDL(s,v) << "))";
			else
				o << s.getName() << "=" << signalValueToVHDL(s,v);
			expected += " " + signalValueToVHDL(s,v,false);
		}

		o << " report \"Incorrect output value for " << s.getName() << ", expected" << expected << "\" severity ERROR; ";
		o << std::endl;
	}

	return o.str();
}

void TestCase::addComment(std::string c)
{
	comment = c;
}

mpz_class TestCase::getInput(const Signal& s)
{
	Inputs::iterator it = inputs.find(s);
	if (it == inputs.end())
		throw std::string("TestCase::getInput: input not found ") + s.getName();
	return it->second;
}

mpz_class TestCase::getOneExpectedOutput(const Signal& s)
{
	std::vector<mpz_class> vs = outputs[s];
	if (vs.size() > 1)
		throw std::string("TestCase::getOneExpectedOutput: multiple expected values for output ") + s.getName();
	if (vs.size() < 1)
		throw 0;	// XXX: Should be a clear exception class.
	return *vs.begin();
}

int TestCase::getExpectedOutputNumber(const Signal& s)
{
	return outputs[s].size();
}

mpz_class TestCase::getExpectedOutput(const Signal& s, int i)
{
	if (i < outputs[s].size())
		return outputs[s][i];
	else
		throw std::string("Expected value index out of range.");
}

TestCaseList::TestCaseList() { }
TestCaseList::~TestCaseList() { }

void TestCaseList::add(TestCase tc)
{
	v.push_back(tc);

	/* Update statistics */
	for (TestCase::Inputs::iterator it = tc.inputs.begin(); it != tc.inputs.end(); it++)
		inputs[it->first] = 0;

	for (TestCase::Outputs::iterator it = tc.outputs.begin(); it != tc.outputs.end(); it++)
	{
		int n = outputs[it->first];
		outputs[it->first] = max(it->second.size(), n);
	}
}

int TestCaseList::getNumberOfTestCases()
{
	return v.size();
}

TestCase TestCaseList::getTestCase(int i)
{
	return v[i];
}

TestCaseList TestCaseList::operator+(TestCaseList second)
{
	TestCaseVector vt;
	TestCaseVector::iterator it;

	TestCaseList tcl;
	tcl.v.reserve(v.size() + second.v.size());
	for (it = v.begin(); it != v.end(); it++)
		tcl.add(*it);
	for (it = second.v.begin(); it != second.v.end(); it++)
		tcl.add(*it);

	return tcl;
}

TestCaseList::Inputs& TestCaseList::getInputMap()
{
	return inputs;
}

TestCaseList::Outputs& TestCaseList::getOutputMap()
{
	return outputs;
}

