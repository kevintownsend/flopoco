
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <python2.7/Python.h>

#include "FloPoCo.hpp" 

#define srcFileName "main_test"


using namespace std;
using namespace flopoco;



void checkOperators ( string operatorName ){
	// TestState used to get back the TestState modified by nextTest of the selected operator
	TestState * cpyTestState;

	// set a new testState and verify if it has not been already treated for the chosen Operator
	if ( operatorName == "FPConstDiv" ){
		// new called only once when first created
		static TestState * testState = new TestState ( "i i i" );
		FPConstDiv::nextTest ( testState );
		cpyTestState = testState;
	}
	else if ( operatorName == "FPMult" ){
		static TestState * testState = new TestState ( "i i i" );
		FPMult::nextTest ( testState );
		cpyTestState = testState;
	}
	else if ( operatorName == "IntMultiplier"){
		static TestState * testState = new TestState ( "i i i b f b" );
		IntMultiplier::nextTest ( testState );
		cpyTestState = testState;
	}
	else if ( operatorName == "IntConstDiv"){
		static TestState * testState = new TestState ( "i i i" );
		IntConstDiv::nextTest ( testState );
		cpyTestState = testState;
	}
	else if ( operatorName == "FPSquare"){
		static TestState * testState = new TestState ( "i i i" );
		FPSquare::nextTest ( testState );
		cpyTestState = testState;
	}
	else{
		REPORT ( LIST, "Operator " << operatorName << " doesn't exist or hasn't been treated yet." );
	}

	// creation of the python cmd with correct arguments
	string cmdPython = operatorName;
	ostringstream pyOss;
	
	string param = cpyTestState -> parameterType;
	int strSize = param.size ();
	int currentStrSize = 0;
	int pos = 0;

	// list of counter used to know the position inside each vectors
	static int counterInt = 0;
	static int counterFloat = 0;
	static int counterString = 0;
	static int counterMpz = 0;
	static int counterBool = 0;

	while ( currentStrSize < strSize ){
		string subParam = param.substr ( pos, 1 );

		// inserting parameters generated with nextTest () inside the command line for python script
		if ( subParam.compare ( "i" ) == 0 ){
			pyOss << " " << cpyTestState -> vectInt [ counterInt ];
			counterInt++;
		}
		else if ( subParam.compare ( "f" ) == 0 ){
			pyOss << " " << cpyTestState -> vectFloat [ counterFloat ];
			counterFloat++;
		}
		else if ( subParam.compare ( "s" ) == 0 ){
			pyOss << " " << cpyTestState -> vectString [ counterString ];
			counterString++;
		}
		else if ( subParam.compare ( "m" ) == 0 ){
			pyOss << " " << cpyTestState -> vectMpz [ counterMpz ];
			counterMpz++;
		}
		else if ( subParam.compare ( "b" ) == 0 ){
			pyOss << " " << cpyTestState -> vectBool [ counterBool ];
			counterBool++;
		}

		pos += 2;
		currentStrSize = pos;
	}

	// reset of each counter after parsing the string containing types of parameters
	counterInt = 0;
	counterFloat = 0;
	counterString = 0;
	counterMpz = 0;
	counterBool = 0;

	cmdPython += pyOss.str ();
	
	// <<<<<<<<<< 
	// Test of the operator with the Python script 
	// >>>>>>>>>>

	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArg, *pValue;

	Py_Initialize ();

	// Setting which python script we will use
	string pythonPathStr = "./tools/check_operator";
	char* pythonPath = new char [ pythonPathStr.size () ];
	strcpy ( pythonPath, pythonPathStr.c_str () );
	int argcTest = 1;
	PySys_SetArgv ( argcTest, &pythonPath );

	const char* pythonScript = "check_operator";
	pName = PyString_FromString ( pythonScript );
	pModule = PyImport_Import ( pName );
	Py_DECREF ( pName );

	if ( pModule != NULL ){
		// Set which function of the python script we will use
		const char* pythonFunction = "checkOp";
		pFunc = PyObject_GetAttrString ( pModule, pythonFunction );
		if ( pFunc && PyCallable_Check ( pFunc ) ){
			// Set arguments of the function we are using
			pArg = PyTuple_New ( 1 );
			pValue = PyString_FromString ( cmdPython.c_str () );
			if ( !pValue ){
				Py_DECREF ( pArg );
				Py_DECREF ( pModule );
				REPORT ( LIST, "Error when try to get back the pythonCMd" );
			}
			PyTuple_SetItem ( pArg, 0, pValue );
		
			// Call of the function with its arguments	
			pValue = PyObject_CallObject ( pFunc, pArg );
			Py_DECREF ( pArg );
			if ( pValue != NULL ){
				Py_DECREF ( pValue );
			}
			else {
				Py_DECREF ( pFunc );
				Py_DECREF ( pModule );
				REPORT ( LIST, "Python Script Failed" );
			}
		}
		else {
			if ( PyErr_Occurred () ){
				REPORT ( LIST, "Python Error" );
			}
			REPORT ( LIST, "Cannot find Python Function " << pythonFunction );
		}
		Py_XDECREF ( pFunc );
		Py_DECREF ( pModule );
	}
	else {
		REPORT ( LIST, "Python Error with python module" << endl << "Failed to load " << pythonScript );
	}
	Py_Finalize ();
}

void checkList ( string listOperator ){
	// if listOperator is  not determined, set default value
	if ( listOperator == "" ){
		listOperator = "soaktestList.txt";
	}
	//open the file containing the list of Operator to treat
	ifstream file;
	file.open ( listOperator.c_str(), fstream::in );
	int iterator = 0;
	if ( file.is_open() ){
		while ( !file.eof() )
		{
			string operatorChosen;
			//get the operator to treat
			getline ( file, operatorChosen );
			if ( operatorChosen.size() != 0 ){
				//launch the treatment on the selected operator
				checkOperators ( operatorChosen );
				iterator++;
			}
			// infinite loop : when end of file reached -> re-open the file and restart
			if ( file.eof()	&& iterator < 12 ){
				file.close ();
				file.open ( listOperator.c_str (), fstream::in );
			}
		}
		file.close ();
	}
	else{
		REPORT ( LIST, "An error happened at the opening of the file " << listOperator );
	}
}


int main(int argc, char* argv[] )
{
	// default list of operator if not determined
	string listOp = "";
	if ( argc >= 1 ){
		// search for the name of the list of operator if precised
		for (int i = 1; i < argc; i++){		
			string option = argv[i];
			int pos;
			if ( option.find ( "list=" ) != string::npos ){
				pos = option.find ( "list=" );
				listOp = option.substr ( pos + 5, option.size() - 5 );
			}
			else{
				REPORT ( LIST, "Wrong arguments passed, waited : \"list=nameList\"" ) ;
				return 0;
			}
		} 
	}
	checkList ( listOp );

	return 0;
}
