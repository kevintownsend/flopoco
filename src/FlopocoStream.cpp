/**
  FloPoCo Stream for VHDL code including cycle information

  This file is part of the FloPoCo project developed by the Arenaire
  team at Ecole Normale Superieure de Lyon
  
  Authors :   Bogdan Pasca, Nicolas Brunie

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
  All rights reserved. */


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <cstdlib>
#include <gmpxx.h>
#include "utils.hpp"
#include "FlopocoStream.hpp"


using namespace std;

namespace flopoco{

	/** The FlopocoStream class.  */
	FlopocoStream::FlopocoStream(){
		vhdlCode.str("");
		vhdlCodeBuffer.str("");
		currentCycle_ = 0;
	}


	FlopocoStream::~FlopocoStream(){
	}

	template <class paramType> FlopocoStream& FlopocoStream::operator <<(
			FlopocoStream& output, paramType c) {
		output.vhdlCodeBuffer << c;
		return output;
	}
	
	FlopocoStream & FlopocoStream::operator<<(FlopocoStream& output, FlopocoStream fs) {
		output.vhdlCodeBuffer << fs.str();
		return output; 
	}
	
	FlopocoStream& FlopocoStream::operator<<( FlopocoStream& output, UNUSED(ostream& (*f)(ostream& fs)) ){
		output.vhdlCodeBuffer << std::endl;
		return output;
	}
	
	string FlopocoStream::str(){
		flush(currentCycle_);
		return vhdlCode.str();
	}

	string FlopocoStream::str(string UNUSED(s) ){
		vhdlCode.str("");
		vhdlCodeBuffer.str("");
		return "";
	}

	void FlopocoStream::flush(int currentCycle){
		if (! disabledParsing ){
			ostringstream bufferCode;
			if ( vhdlCodeBuffer.str() != string("") ){
				/* do processing if buffer is not empty */

				/* scan buffer sequence and annotate ids */
				bufferCode << annotateIDs( currentCycle );

				/* the newly processed code is appended to the existing one */
				vhdlCode << bufferCode.str();

			}
		}else{
			vhdlCode << vhdlCodeBuffer.str();				
		}
		/* reset buffer */
		vhdlCodeBuffer.str(""); 
	}

	void FlopocoStream::flush(){
		flush ( currentCycle_ );
	}
}
