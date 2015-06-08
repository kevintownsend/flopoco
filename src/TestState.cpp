#include "TestState.hpp"

namespace flopoco {

	TestState::TestState ( string param ) : paramTypes ( param ) {
		counter = 0; /**Initialize the counter of the instance*/

		// integers corresponding to the size of vectors
		int nbInt = 0;
		int nbFloat = 0;
		int nbString = 0;
		int nbMpz = 0;
		int nbBool = 0;

		// parsing the string containing type of parameters for the operator
		// previous integers will then be updated
		unsigned int strSize = 0;
		unsigned int pos = 0;
		while ( strSize <  param.size () ){
			string subParam = param.substr ( pos, 1 );
			pos += 2;
			strSize = pos;

			if ( subParam.compare ( "i" ) == 0 ){
				nbInt += 1;
			}
			else if ( subParam.compare ( "f" ) == 0 ){
				nbFloat += 1;
			}
			else if ( subParam.compare ( "s" ) == 0 ){
				nbString += 1;
			}
			else if ( subParam.compare ( "m" ) == 0 ){
				nbMpz += 1;
			}
			else if ( subParam.compare ( "b" ) == 0 ){
				nbBool += 1;
			}
		}
		// initialize all vectors to the requested size
		vectInt.resize ( nbInt );
		vectFloat.resize ( nbFloat );
		vectString.resize ( nbString );
		vectMpz.resize ( nbMpz );
		vectBool.resize ( nbBool );
	}

	/**
	* Test if the TestState we are working on is equal to the TestState ts
	* return a boolean
	**/
	bool TestState::equality ( TestState * ts ) {
		string s1 =  toString();
		string s2 = ts->toString();
		return (s1.compare(s2) == 0) ;
	}

	string TestState::toString() {
		// list of counter used to know the position inside each vectors
		int counterInt = 0;
		int counterFloat = 0;
		int counterString = 0;
		int counterMpz = 0;
		int counterBool = 0;
		int strSize = paramTypes.size ();
		int currentStrSize = 0;
		int pos=0;
		ostringstream s;
		while ( currentStrSize < strSize ){
			string subParam = paramTypes.substr ( pos, 1 );

			if ( subParam.compare ( "i" ) == 0 ){
				s << " " << vectInt [ counterInt ];
				counterInt++;
			}
			else if ( subParam.compare ( "f" ) == 0 ){
				s << " " << vectFloat [ counterFloat ];
				counterFloat++;
			}
			else if ( subParam.compare ( "s" ) == 0 ){
				s << " " << vectString [ counterString ];
				counterString++;
			}
			else if ( subParam.compare ( "m" ) == 0 ){
				s << " " << vectMpz [ counterMpz ];
				counterMpz++;
			}
			else if ( subParam.compare ( "b" ) == 0 ){
				s << " " << vectBool [ counterBool ];
				counterBool++;
			}

			pos += 2;
			currentStrSize = pos;
		}
		return s.str();
	}


}
