#ifndef TESTSTATE_HPP
#define TESTSTATE_HPP

#include <ostream>
#include <sstream>
#include <vector>
#include <string>
#include <gmpxx.h>

using namespace std;
namespace flopoco {
	class TestState
	{
	public :
		TestState ( string param );
	
		/**
		 * Test the equality between two TestState
		 **/
		bool equality (  TestState * ts );

		vector < int > vectInt;
		vector < float > vectFloat;
		vector < string > vectString;
		vector < mpz_class > vectMpz;
		vector < bool > vectBool;

		string paramTypes; /** String representing type of parameters of the chosen operator int the right order */
		string toString(); /** get a textual representation of the state, to be passed to the flopoco command line for instance*/
		int counter;   /**< counter to be aware of how many test we have done in this instance*/
	};
}
#endif
