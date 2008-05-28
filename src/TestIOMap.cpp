#include "TestIOMap.hpp"
#include "Signal.hpp"

TestIOMap::TestIOMap() { }
TestIOMap::~TestIOMap() { }

void TestIOMap::add(Signal s, int maxNumValues)
{
	a.push_back(std::pair<Signal, int>(s, maxNumValues));	
}

const TestIOMap::TestIOMapContainer& TestIOMap::get() const
{
	return a;
}

