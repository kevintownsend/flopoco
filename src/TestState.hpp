#ifndef TESTSTATE_HPP
#define TESTSTATE_HPP

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
			bool equality ( const TestState * ts ) const;

			vector < int > vectInt;
			vector < float > vectFloat;
			vector < string > vectString;
			vector < mpz_class > vectMpz;
			vector < bool > vectBool;

			string parameterType; /** String representing type of parameters of the chosen operator int the right order */

			int counter;   /**< counter to be aware of how many test we have done in this instance*/
	};
}
#endif
