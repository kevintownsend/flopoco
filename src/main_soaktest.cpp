
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <python2.7/Python.h>

#include "FloPoCo.hpp"

// for REPORT and THROWERROR
#define srcFileName "main_test"
#define uniqueName_ "main_test"


using namespace std;
using namespace flopoco;



void checkOperator ( string operatorName ){
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
		REPORT ( LIST, "Operator " << operatorName << " doesn't exist in checkOperator yet." );
	}

	// creation of the python cmd with correct arguments
	string cmdPython = operatorName;
	ostringstream pyOss;

	string params = cpyTestState -> toString();

	cmdPython += params;
	REPORT (LIST, "Testing Operator " << operatorName << " with command: " << cmdPython);

	// <<<<<<<<<<
	// Test of the operator with the Python script
	// >>>>>>>>>>

	PyObject *pName, *pModule, /**pDict,*/ *pFunc;
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





vector<string> readOpList ( string opListFileName ){
	//First read the file
	vector<string> operators;
	string opName;
	ifstream file;
	file.open ( opListFileName.c_str(), fstream::in );
	if ( file.is_open() ){
		while ( !file.eof() )
		{
			//get the operator to treat
			getline ( file, opName );
			if ( opName.size() != 0 ){
				operators.push_back(opName);
				REPORT(LIST, "Going to test " << opName);
			}
		}
		file.close ();

		//launch the treatment on the selected operator
		//		checkOperator ( opName );
	}
	else{
		cerr << "Could not open file " << opListFileName << endl;
	}
	return operators;
}


int main(int argc, char* argv[] )
{
	// default list of operator if not determined
	string listOp = "soaktestList.txt";	//  default value

	if ( argc >= 2 ){
		listOp = argv[1];
	}

	vector<string> opList;
	opList = readOpList ( listOp );

		while(true) {
		for (unsigned int i=0; i<opList.size(); i++)
				checkOperator (opList[i]);
		}
	return 0;
}
