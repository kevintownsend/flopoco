/* vim: set tabstop=8 softtabstop=2 shiftwidth=2: */
#include "TestCase.hpp"

// XXX: This might be useful elsewhere. How about moving it to a more general
// place?
/**
 * Converts the value of the signal into a nicely formated VHDL expression,
 * including padding and putting quot or apostrophe.
 * @param s signal (used to determine the width)
 * @param v value
 * @return a VHDL value expression
 */
std::string signalValueToVHDL(Signal s, mpz_class v)
{
  std::string o;

  /* Get base 2 representation */
  o = v.get_str(2);

  /* Some check */
  if (o.size() > s.width())
  {
    std::ostringstream o;
    o << "Error in " <<  __FILE__ << "@" << __LINE__ << ": value is larger than signal " << s.id();
    throw o.str();
  }

  /* Do padding */
  while (o.size() < s.width())
    o = "0" + o;

  /* Put apostrophe / quot */
  if (s.width() > 1)
    return "\"" + o + "\"";
  else
    return "'" + o + "'";
}

TestCase::TestCase() { }
TestCase::~TestCase() { }

void TestCase::addInput(Signal s, mpz_class v)
{
  inputs[s] = v;
}

void TestCase::addExpectedOutput(Signal s, mpz_class v)
{
  outputs[s].insert(v);
}

std::string TestCase::getInputVHDL()
{
  std::ostringstream o;

  for (Inputs::iterator it = inputs.begin(); it != inputs.end(); it++)
  {
    Signal s = it->first;
    mpz_class v = it->second;
    o << s.id() << " <= " << signalValueToVHDL(s, v) << "; ";
  }
  o << "-- " << comment;

  return o.str();
}

std::string TestCase::getExpectedOutputVHDL()
{
  std::ostringstream o;

  /* Iterate through output signals */
  for (Outputs::iterator it = outputs.begin(); it != outputs.end(); it++)
  {
    Signal s = it->first;
    std::set<mpz_class> vs = it->second;

    o << "assert false";  // XXX: Too lazy to make an exception for the first value

    /* Iterate through possible output values */
    for (std::set<mpz_class>::iterator it = vs.begin(); it != vs.end(); it++)
    {
      mpz_class v = *it;
      o << " or " << s.id() << "=" << signalValueToVHDL(s,v);
    }

    o << " report \"Incorrect output value for " << s.id() << "\" severity ERROR; ";
  }
  o << "-- " << comment;

  return o.str();
}

void TestCase::addComment(std::string c)
{
  comment = c;
}

TestCaseList::TestCaseList() { }
TestCaseList::~TestCaseList() { }

void TestCaseList::add(TestCase t)
{
  v.push_back(t);
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

  for (it = v.begin(); it != v.end(); it++)
    vt.push_back(*it);
  for (it = second.v.begin(); it != second.v.end(); it++)
    vt.push_back(*it);

  return TestCaseList(vt);
}

