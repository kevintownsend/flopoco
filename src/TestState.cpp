#include "TestState.hpp"

namespace flopoco {
	
	TestState::TestState ( string param ) : parameterType ( param ) {
		counter = 0; /**Initialize the counter of the instance*/

		// integers corresponding to the size of vectors
		int nbInt = 0;
		int nbFloat = 0;
		int nbString = 0;
		int nbMpz = 0;
		int nbBool = 0;

		// parsing the string containing type of parameters for the operator
		// previous integers will then be updated
		int strSize = 0;
		int pos = 0;
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
	bool TestState::equality ( const TestState * ts ) const {
		bool result = false;

		int strSize = this -> parameterType.size ();
		int currentStrSize = 0;
		int pos = 0;
		
		// counters to check the right position in each vectors
		static int counterInt = 0;
		static int counterFloat = 0;
		static int counterString = 0;
		static int counterMpz = 0;
		static int counterBool = 0;
		
		while ( currentStrSize < strSize ){
			string subParam = this -> parameterType.substr ( pos, 1 );

			if ( subParam.compare ( "i" ) == 0 ){
				if ( this -> vectInt [ counterInt ] == ts -> vectInt [ counterInt ] ){
					result = true;
				}
				else{
					return false;
				}
				counterInt++;
			}
			else if ( subParam.compare ( "f" ) == 0 ){
				if ( this -> vectFloat [ counterFloat ] == ts -> vectFloat [ counterFloat ] ){
					result = true;
				}
				else{
					return false;
				}
				counterFloat++;
			}
			else if ( subParam.compare ( "s" ) == 0 ){
				if ( this -> vectString [ counterString ] == ts -> vectString [ counterString ] ){
					result = true;
				}
				else{
					return false;
				}
				counterString++;
			}
			else if ( subParam.compare ( "m" ) == 0 ){
				if ( this -> vectMpz [ counterMpz ] == ts -> vectMpz [ counterMpz ] ){
					result = true;
				}
				else{
					return false;
				}
				counterMpz++;
			}
			else if ( subParam.compare ( "b" ) == 0 ){
				if ( this -> vectBool [ counterBool ] == ts -> vectBool [ counterBool ] ){
					result = true;
				}
				else{
					return false;
				}
				counterBool++;
			}

			pos += 2;
			currentStrSize = pos;
		}
		//reset all counters after the parsing
		counterInt = 0;
		counterFloat = 0;
		counterString = 0;
		counterMpz = 0;
		counterBool = 0;
		return result;
	}
}
