#ifndef __TESTCASE_HPP
#define __TESTCASE_HPP

#include <string>
#include <vector>
#include <map>
#include <gmpxx.h>

#include "Signal.hpp"

class TestCase;

/**
 * Represents a list of test cases that an Operator has to pass.
 */
class TestCaseList {
public:
	/**
	 * Creates an empty TestCaseList
	 * @see TestCase
	 */
	TestCaseList();
	
	/** Destructor */
	~TestCaseList();

	/**
	 * Adds a TestCase to this TestCaseList.
	 * @param t TestCase to add
	 */
	void add(TestCase t);

	/* XXX Is this the best way to encapsulate TestCase-es? */
	/**
	 * Get the number of TestCase-es contained in this TestCaseList.
	 * @return The number of TestCase-es
	 */
	int getNumberOfTestCases();

	/* XXX Shouldn't we return by reference to improve speed? */
	/**
	 * Get a specific TestCase.
	 * @param i Which TestCase to return
	 * @return The i-th TestCase.
	 */
	TestCase getTestCase(int i);

	/**
	 * Creates a new TestCaseList, with TestCase-es from both lists.
	 * @param second The second operand of the plus operator.
	 */
	TestCaseList operator+(TestCaseList second);

	/**
	 * Offers an ordering relationship on Signal-s based on their name.
	 * Now we can use Signal-s as keys to map, set etc.
	 */
	struct ltsignal
	{
		bool operator()(const Signal& s1, const Signal& s2) const
		{
			return (s1.getSignalName() < s2.getSignalName());
		}
	};

	/** Map used to store input signals */
	typedef std::map<Signal, int, ltsignal> Inputs;
	
	/** Map used to store output signals */
	typedef std::map<Signal, int, ltsignal> Outputs;

	/**
	 /* Gets a map with the input signals
	 */
	TestCaseList::Inputs& getInputMap();

	/**
	 * Gets a map with the output signals and they maximum number
	 * of expected values.
	 */
	TestCaseList::Outputs& getOutputMap();

private:
	/** Type of vector used to store TestCase-es */
	typedef std::vector<TestCase> TestCaseVector;

	/**
	 * Creates a TestCaseList from a vector of TestCase-es.
	 * Internal constructor used by plus operator.
	 * @param v The vector to use.
	 */
	inline TestCaseList(TestCaseVector v, Inputs inputs, Outputs outputs)
		: v(v), inputs(inputs), outputs(outputs) { }

	/** Stores the TestCase-es */
	TestCaseVector v;

	/** Stores statistics about the TestCase-es */
	Inputs  inputs;  /*< All input signals */
	Outputs outputs; /*< All output signals and their number of expected values */
};

/**
 * Encapsulates data relevant for a test case. The data consists
 * of one value for each input signal, zero, one or more expected
 * values for each output signal. The idea is also to have
 * „good enough” constructors, so that test cases are easily
 * added from all operators.
 * @see TestBench
 * @see Operator
 */
class TestCase {
public:
	/** Creates an empty TestCase */
	TestCase();
	~TestCase();

	/**
	 * Adds an input for this TestCase
	 * @param s The signal which will value the given value
	 * @param v The value which will be assigned to the signal
	 */
	void addInput(const Signal &s, mpz_class v);

	/**
	 * Adds an expected output for this signal
	 * @param s The signal for which to assign an expected output
	 * @param v One possible value which the signal might have
	 */
	void addExpectedOutput(const Signal &s, mpz_class v);

	/**
	 * Adds a comment to the output VHDL. "--" are automatically prepended.
	 * @param c Comment to add.
	 */
	void addComment(std::string c);

	/**
	 * Generates the VHDL code necessary for assigning the input signals.
	 * @param prepend A string to prepend to each line.
	 * @return A multi-line VHDL code.
	 */
	std::string getInputVHDL(std::string prepend = "");

	/**
	 * Generate the VHDL code necessary to assert the expected output
	 * signals.
	 * @param prepend A string to prepend to each line.
	 * @return A single-line VHDL code.
	 */
	std::string getExpectedOutputVHDL(std::string prepend = "");

	/**
	 * Gets the input value stored in the test case for the given
	 * Signal. Throws an exception if no signal is found.
	 * @param s the signal for which the value should be returned.
	 */
	mpz_class getInput(const Signal& s);

	/**
	 * Gets the one single value of for the given output Signal.
	 * Throws an exception if no signal value is found or there
	 * are two expected outputs.
	 */
	mpz_class getOneExpectedOutput(const Signal& s);

	// XXX: Not the best place.
	/**
	 * Converts the value of the signal into a nicely formated VHDL expression,
	 * including padding and putting quot or apostrophe.
	 * @param s signal (used to determine the width)
	 * @param v value
	 * @param quot also put quotes around the value
	 * @return a VHDL value expression
	 */
	static std::string signalValueToVHDL(const Signal& s, mpz_class v, bool quot = true);
	 
	/**
	 * Converts the value of the signal into a nicely formated VHDL expression,
	 * including padding and putting quot or apostrophe. (Hex version.)
	 * @param s signal (used to determine the width)
	 * @param v value
	 * @param quot also put quotes around the value
	 * @return a VHDL value expression in hexa
	 */
	static std::string signalValueToVHDLHex(const Signal& s, mpz_class v, bool quot = true);

	/**
	 * Returns the number of expected output values
	 * @param s signal
	 * @return an integer
	 */
	int getExpectedOutputNumber(const Signal& s);

	/**
	 * Returns one expected output value
	 * @param s signal
	 * @param i zero-indexed expected output to return
	 * @return mpz_class representing a VHDL signal value
	 */
	mpz_class getExpectedOutput(const Signal& s, int i);

private:
	/**
	 * Offers an ordering relationship on Signal-s based on their name.
	 * Now we can use Signal-s as keys to map, set etc.
	 */
	struct ltsignal
	{
		bool operator()(const Signal& s1, const Signal& s2) const
		{
			return (s1.getSignalName() < s2.getSignalName());
		}
	};

	typedef std::map<Signal, mpz_class, ltsignal> Inputs;
	Inputs inputs;

	typedef std::map<Signal, std::vector<mpz_class>, ltsignal> Outputs;
	Outputs outputs;

	std::string comment;

	friend void TestCaseList::add(TestCase);
};

#endif

