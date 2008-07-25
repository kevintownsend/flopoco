#ifndef TESTIOMAP_HPP
#define TESTIOMAP_HPP

#include <vector>

/* Opaque declaration */
class Signal;

/**
 * Stores the inputs and outputs relevant for generating the test case
 * of an operator. For outputs, it also stores the maximum number of values
 * they may store.
 */
class TestIOMap {
public:
	/** Default Constructor */
	TestIOMap();
	
	/** Destructor */
	~TestIOMap();

	/** The data-type which stores the test case input/output map */
	typedef struct std::vector<std::pair<Signal, int> > TestIOMapContainer;

	/**
	 * Adds a new signal to the test input/output map.
	 * @param s the signal to add
	 * @param maxNumValues the maximum number of values that the output signal
	 *                     may have.
	 */
	void add(Signal s, int maxNumValues = 1);

	/**
	 * Returns the already stored input/output map.
	 * @return a TestIOMapContainer reference.
	 */
	const TestIOMapContainer& get() const;

private:
	/** Stores the test input/output map */
	TestIOMapContainer a;
};

#endif

