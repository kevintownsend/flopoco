/*
An integer adder for FloPoCo

It may be pipelined to arbitrary frequency.
Also useful to derive the carry-propagate delays for the subclasses of Target

Authors:  Bogdan Pasca, Florent de Dinechin

This file is part of the FloPoCo project
developed by the Arenaire team at Ecole Normale Superieure de Lyon

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
CeCILL license, 2008-2010.
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntAdderClassical.hpp"

using namespace std;

namespace flopoco {
	extern vector<Operator*> oplist;
	
	IntAdderClassical::IntAdderClassical ( Target* target, int wIn, string name, map<string, double> inputDelays, int optimizeType, bool srl) :
	Operator ( target, inputDelays ), wIn_ ( wIn ) {
		srcFileName="IntAdderClassical";
		setCopyrightString ( "Bogdan Pasca, Florent de Dinechin (2008-2010)" );
		setName( name );
		
		// Set up the IO signals
		addInput  ( "X"  , wIn_, true );
		addInput  ( "Y"  , wIn_, true );
		addInput  ( "Cin", 1 );
		addOutput ( "R"  , wIn_, 1 , true );
		
		objectivePeriod	 = 1 / target->frequency();
		maxInputDelay      = getMaxInputDelays ( inputDelays );
		if ( maxInputDelay > objectivePeriod ) {
			nextCycle();
			maxInputDelay = 0.0;
			maxInputDelay = target->ffDelay() + target->localWireDelay();	
			inputDelays.clear();
		}	
		
		classicalSlackVersion = -1;
		
		switch (optimizeType) {
			case 0:  cost = getLutCostClassical(target,wIn, inputDelays, srl); break;
			case 1:  cost = getRegCostClassical(target,wIn, inputDelays, srl); break;
			case 2:  cost = getSliceCostClassical(target,wIn, inputDelays, srl); break;
			default: cost = getSliceCostClassical(target,wIn, inputDelays, srl); break;
		}
		

		
		vhdl << tab << "--Classical"<<endl;
		if ( isSequential() ) {
			if ( maxInputDelay == (target->localWireDelay() + target->ffDelay()) || classicalSlackVersion==-1 ) {
				/* the non-slack version */
				updateParameters ( target, alpha, beta, k );
				REPORT ( DEBUG, "1) alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			} else {
				if ( classicalSlackVersion == 0 ) {
					/* the slack version that does not buffer the inputs*/
					updateParameters ( target, inputDelays, alpha, beta, gamma, k );
					REPORT ( DEBUG, "alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" k="<<k );
				} else {
					nextCycle(); /* bufferning the inputs */
					REPORT ( DETAILED, "Building architecture for classical version with slack: buffering" );
					updateParameters ( target, alpha, beta, k );
					REPORT ( DEBUG, "alpha="<<alpha<<" beta="<<beta<<" k="<<k );
				}
			}
			
			if ( k>1 ) {
				/* init the array with chunk sizes */
				cSize = new int[k+1];
				if ( (maxInputDelay != 0 ) && (classicalSlackVersion == 0) )
					cSize[0] = gamma;
				else
					cSize[0] = alpha;

				for ( int i=1; i<k-1; i++ )
					cSize[i] = alpha;
				
				cSize[k-1] = beta;
			
				/* the indexes of the chunks */
				cIndex = new int[k];
				cIndex[0] = cSize[0];
				for ( int i=1; i < k; i++ )
					cIndex[i] = cIndex[i-1] + cSize[i];
				
				/* The implementation */
				for ( int i=0; i < k; i++ ) {
					vhdl << tab << declare ( join ( "x",i ), cSize[i],true ) << " <= X" << range ( cIndex[i]-1, ( i>0?cIndex[i-1]:0 ) ) << ";" << endl;
					vhdl << tab << declare ( join ( "y",i ), cSize[i],true ) << " <= Y" << range ( cIndex[i]-1, ( i>0?cIndex[i-1]:0 ) ) << ";" << endl;
				}
				
				for ( int i=0; i < k; i++ ) {
					vhdl << tab << declare ( join ( "sum",i ),cSize[i]+1,true ) << " <= ( \"0\" & "<< join ( "x",i ) << ") + ( \"0\" & "<< join ( "y",i ) << ")  + ";
					if ( i==0 ) vhdl << "Cin";
					else      vhdl << join ( "sum",i-1 ) <<of ( cSize[i-1] );
					vhdl << ";"	<< endl;
					if ( i < k-1 )
						nextCycle(); //////////////////////////////////////////
				}
				
				vhdl << tab << "R <= ";
				for ( int i=k-1; i >= 1; i-- ) {
					vhdl << join ( "sum",i ) <<range ( cSize[i]-1,0 ) << " & ";
				}
				vhdl << "sum0" << range ( cSize[0]-1,0 ) <<";"<<endl;
				/* the output is asociated with the combinatorial delay caused
				by the most-significant bits addition */
				outDelayMap["R"] = target->adderDelay ( cSize[k-1] ) + ( getCurrentCycle() >0?0:getMaxInputDelays ( inputDelays ) );
			} else {
				vhdl << tab << " R <= X + Y + Cin;" << endl;
				if ( (maxInputDelay == 0 ) || (classicalSlackVersion==0) )
					outDelayMap["R"] = target->adderDelay ( wIn ) + getMaxInputDelays ( inputDelays );
				else 
					outDelayMap["R"] = target->adderDelay ( wIn );
			}
	
		} else {
			vhdl << tab << " R <= X + Y + Cin;" << endl;
			outDelayMap["R"] = target->adderDelay ( wIn_ ) + getMaxInputDelays ( inputDelays );
		}
	}
	
	/**************************************************************************/
	IntAdderClassical::~IntAdderClassical() {
	}
	

	/**************************************************************************/
	/*******************             LUTS             *************************/
	/**************************************************************************/
	int IntAdderClassical::getLutCostClassical ( Target* target, int wIn, map<string, double> inputDelays, bool srl ) {
		REPORT ( DEBUG, DEBUG_SEPARATOR );
		if ( getMaxInputDelays ( inputDelays ) == 0 ) {
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters ( target, alpha, beta, k );
			REPORT ( DEBUG, "LUT, Classical, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			if ( srl ) /* Allow hardware shift-registers */
				if ( ( k == 1 ) || ( k == 2 ) )  cost = wIn;
				else                       cost = wIn + ( k-2 ) *alpha; //(k-2)alpha SRLs
					else /* Don't allow hardware shift-registers */
						cost = wIn;
					REPORT ( DETAILED, "Selected: Classical, NO-SLACK with LUT cost " << cost );
				return cost;
		} else {
			int version0, version1;
			int alpha, beta, gamma, k;
			/* Version 0: we try to adapt the architecture to the new slack */
			updateParameters ( target, inputDelays, alpha, beta, gamma, k );
			REPORT ( DEBUG, "LUT, Classical, SLACK, Version 0: alpha="<<alpha<<" beta`="<<beta<<" gamma="<<gamma<<" k="<<k );
			
			if ( k>0 )
				/* an alternative possible splitting */
				if ( srl ) {
					if ( ( k==1 ) || ( k==2 ) )
						version0 = wIn;
					else
						version0 = 4*wIn - 3*alpha - 2*gamma  - beta;
				} else {
					/* NO SRLs */
					version0 = wIn;
				}
				else {
					/* no solution was found, cost is +INF */
					version0 = PINF;
				}
				
				/* Version 1: we buffer the inputs and proceed */
				updateParameters ( target, alpha, beta, k );
				REPORT ( DEBUG, "LUT, Classical, SLACK, Version 1 (buffered inputs): alpha="<<alpha<<" beta="<<beta<<" k="<<k << " p="<<k+1 );
				
				if ( srl ) {
					if ( k==1 ) {
						version1 = wIn;
					} else if ( k == 2 ) {
						version1 = w + 2*beta;
					} else {
						version1 = 4*wIn - 3*alpha - beta;
					}
				} else {
					version1 = wIn;
				}
				
				REPORT ( DEBUG, "LUT, Classical, SLACK, Version 0: " << version0 );
				REPORT ( DEBUG, "LUT, Classical, SLACK, Version 1 (buffered inputs) " << version1 );
				
				if ( version0 <= version1 ) {
					/* for equality version 0 has less pipeline levels */
					classicalSlackVersion = 0;
					REPORT ( DETAILED, "Selected: Classical SLACK version is 0 with LUT cost " << version0 );
					return version0;
				} else {
					classicalSlackVersion = 1;
					REPORT ( DETAILED, "Selected: Classical SLACK version is 1 (buffered inputs) with LUT cost " << version1 );
					return version1;
				}
		}
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit ( -1 );
		return -1;
	}
	

	
	/**************************************************************************/
	/*******************          REGISTERS           *************************/
	/**************************************************************************/
	int IntAdderClassical::getRegCostClassical ( Target* target, int wIn, map<string, double> inputDelays, bool srl ) {
		REPORT ( DEBUG, DEBUG_SEPARATOR );
		if ( getMaxInputDelays ( inputDelays ) == 0 ) {
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters ( target, alpha, beta, k );
			REPORT ( DEBUG, "REG, CLASSICAL, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			if ( k == 1 ) {
				cost = 0;
			} else {
				if ( srl )
					cost = wIn - beta;
				else
					cost = ( ( 3*k*k -7*k+4 ) *alpha/2 + 2* ( k-1 ) *beta + k-1 );
			}
			
			REPORT ( DETAILED, "Selected: Classical, NO-SLACK, with REG cost " << cost );
			return cost;
			
		} else {
			int version0, version1;
			int alpha, beta, gamma, k;
			
			/* Version 0: we try to adapt the architecture to the new slack */
			updateParameters ( target, inputDelays, alpha, beta, gamma, k );
			REPORT ( DEBUG, "REG, Classical, SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" k="<<k );
			
			if ( k>0 ) {
				if ( k==1 )
					version0 = 0;
				
				if ( srl ) {
					if ( k == 2 )
						version0 = wIn + alpha + 1;
					else
						version0 = 3*wIn - 2*gamma  - beta + k - 1;
				} else {
					if ( k == 2 ) {
						version0 = wIn + alpha + 1;
					} else {
						version0 = ( k-1 ) *gamma + 2* ( k-1 ) *beta + 3* ( k*k-3*k+2 ) /2;
					}
				}
				
			} else
				version0 = PINF; //infinity
				
				/* Version 1: we buffer the inputs and proceed */
				updateParameters ( target, alpha, beta, k );
			REPORT ( DEBUG, "REG, Classical, SLACK, Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			
			if ( srl ) {
				if ( k==1 ) {
					version1 = wIn;
				} else if ( k == 2 ) {
					version1 = 3*alpha + 2*beta + 2;
				} else {
					version1 = 3*wIn - beta + k - 1;
				}
			} else {
				if ( k==1 )
					version1 = wIn;
				else
					version1 = wIn + ( ( 3*k*k -7*k+4 ) *alpha/2 + 2* ( k-1 ) *beta + k-1 );
			}
			
			REPORT ( DEBUG, "REG, Classical, SLACK, Version 0: " << version0 );
			REPORT ( DEBUG, "REG, Classical, SALCK, Version 1 (buffered inputs): " << version1 );
			
			if ( version0 <= version1 ) {
				classicalSlackVersion = 0;
				REPORT ( DETAILED, "Selected: Classical SLACK version is 0 with REG cost " << version0 );
				return version0;
			} else {
				classicalSlackVersion = 1;
				REPORT ( DETAILED, "Selected: Classical SLACK version is 1 with REG cost " << version1 );
				return version1;
			}
		}
		
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit ( -1 );
		return -1;
	}
	

	
	/**************************************************************************/
	/*******************           SLICES             *************************/
	/**************************************************************************/
	int IntAdderClassical::getSliceCostClassical ( Target* target, int wIn, map<string, double> inputDelays, bool srl ) {
		REPORT ( DEBUG, DEBUG_SEPARATOR );
		if ( getMaxInputDelays ( inputDelays ) == 0 ) {
			/* no input slack problem */
			int alpha, beta, k, cost;
			updateParameters ( target, alpha, beta, k );
			REPORT ( DEBUG, "SLICE, CLASSICAL, NO-SLACK: alpha="<<alpha<<" beta="<<beta<<" k="<<k );
			
			if ( k == 1 ) {
				cost = int ( ceil ( double ( wIn ) / double ( 2 ) ) );
			} else {
				if ( srl )
					cost = int ( ceil ( double ( wIn + ( k-2 ) *alpha  + k - 1 ) / double ( 2 ) ) );
				else
					cost = int ( ceil ( double ( wIn + ( ( 3* ( k*k-3*k+2 ) ) /2 ) *alpha + 2* ( k-1 ) *beta ) / double ( 2 ) ) );
			}
			REPORT ( DETAILED, "Selected: Classical NO-SLACK, with SLICE cost " << cost );
			return cost;
			
		} else {
			int version0, version1;
			int alpha, beta, gamma, k;
			
			/* Version 0 */
			updateParameters ( target, inputDelays, alpha, beta, gamma, k );
			REPORT ( DEBUG, "SLICE, Classical, SLACK, Version 0: alpha="<<alpha<<" beta="<<beta<<"gamma="<<gamma<< " k="<<k );
			
			if ( k>0 ) {
				if ( k==1 )
					version0= int ( ceil ( double ( wIn ) /double ( 2 ) ) );
				
				if ( srl ) {
					if ( k == 2 ) {
						version0 = int ( ceil ( double ( 3*alpha + gamma + 1 ) /double ( 2 ) ) );
					} else {
						version0 = int ( ceil ( double ( 4*wIn - alpha - 2*gamma -2*beta + k - 1 ) /double ( 2 ) ) );
					}
				} else {
					if ( k == 2 ) {
						version0 = int ( ceil ( double ( gamma + 3*alpha + 1 ) / double ( 2 ) ) );
					} else {
						version0 = int ( ceil ( double ( wIn + ( k-2 ) *gamma + 2* ( k-1 ) *beta + alpha* ( 2*k*k-11*k+10 ) /2 ) /double ( 2 ) ) );
					}
				}
			} else
				version0 = PINF; //infinity
				
				/* Version 1 */
				updateParameters ( target, alpha, beta, k );
			REPORT ( DEBUG, "SLICE, Classical, SLACK Version 1: alpha="<<alpha<<" beta="<<beta<<" k="<<k<<" p="<<k+1 );
			
			if ( k==1 )
				version1= int ( ceil ( double ( wIn ) /double ( 2 ) ) );
			
			if ( srl ) {
				if ( k == 2 ) {
					version1 = int ( ceil ( double ( 3*alpha + 2*beta + 2 ) /double ( 2 ) ) );
				} else {
					version1= int ( ceil ( double ( 4*wIn - alpha - beta ) /double ( 2 ) ) );
				}
			} else {
				if ( k == 2 ) {
					version1 = int ( ceil ( double ( 2*wIn + ( ( 3* ( k*k-3*k+2 ) ) /2 ) *alpha + 2* ( k-1 ) *beta ) / double ( 2 ) ) );
				} else {
					version1 = int ( ceil ( double ( 3*wIn + ( ( 3* ( k*k-3*k+2 ) ) /2 ) *alpha + 2* ( k-1 ) *beta ) / double ( 2 ) ) );
				}
			}
			
			REPORT ( DETAILED, "SLICE, Classical, SLACK, Version 0: " << version0 );
			REPORT ( DETAILED, "SLICE, Classical, SLACK, Version 1: " << version1 );
			
			if ( version0 <= version1 ) {
				classicalSlackVersion = 0;
				REPORT ( DETAILED, "Selected: Classical SLACK version is 0 with SLICE cost " << version0 );
				return version0;
			} else {
				classicalSlackVersion = 1;
				REPORT ( DETAILED, "Selected: Classical SLACK version is 1 with SLICE cost " << version1 );
				return version1;
			}
		}
		
		cerr << "Error in " <<  __FILE__ << "@" << __LINE__;
		exit ( -1 );
		return -1;
	}
	

	
	/**************************************************************************/
	void IntAdderClassical::updateParameters ( Target* target, int &alpha, int &beta, int &k ) {
		
		target->suggestSlackSubaddSize ( alpha , wIn_, target->ffDelay() + target->localWireDelay() ); /* chunk size */
		if ( wIn_ == alpha ) {
			/* addition requires one chunk */
			beta = 0;
			k    = 1;
		} else {
			beta = ( wIn_ % alpha == 0 ? alpha : wIn_ % alpha );
			k    = ( wIn_ % alpha == 0 ? wIn_ / alpha : ceil ( double ( wIn_ ) / double ( alpha ) ) );
		}
	}
	
	/**************************************************************************/
	void IntAdderClassical::updateParameters ( Target* target, map<string, double> inputDelays, int &alpha, int &beta, int &gamma, int &k ) {
		
		int typeOfChunks = 1;
		bool status = target->suggestSlackSubaddSize ( gamma , wIn_, getMaxInputDelays ( inputDelays ) ); // the first chunk size
		REPORT ( DEBUG, "suggestSlackSubaddSize returns gamma="<<gamma<<" with status:"<< ( status?"true":"false" ) );
		
		if ( ! status ) {
			k=-1;
			alpha=0;
			beta=0;
			gamma=0;
		} else
			if ( wIn_ - gamma > 0 ) { //more than 1 chunk
				target->suggestSlackSubaddSize (alpha, wIn_-gamma, target->ffDelay() + target->localWireDelay());
				if ( wIn_-gamma == alpha )
					typeOfChunks++; //only two types of chunks
					else
						typeOfChunks+=2; //three types of chunks
						
						REPORT ( DETAILED, "Types of chunks = " << typeOfChunks );
					
					if ( typeOfChunks==3 )
						beta = ( ( wIn_-gamma ) % alpha == 0 ? alpha : ( wIn_-gamma ) % alpha );
					else
						beta = alpha;
					
					
					if ( typeOfChunks==2 )
						k = 2;
					else
						k = 2 +   int ( ceil ( double ( wIn_ - beta - gamma ) / double ( alpha ) ) );
					
					
			} else {
				alpha = 0;
				beta = 0;
				k=1;
			}
	}
	
	/**************************************************************************/
	void IntAdderClassical::updateParameters ( Target* target, map<string, double> inputDelays, int &alpha, int &beta, int &k ) {
		bool status = target->suggestSlackSubaddSize ( alpha , wIn_,  getMaxInputDelays ( inputDelays ) ); /* chunk size */
		if ( !status ) {
			k=-1;
			alpha=0;
			beta=0;
		} else
			if ( wIn_ == alpha ) {
				/* addition requires one chunk */
				beta = 0;
				k    = 1;
			} else {
				beta = ( wIn_ % alpha == 0 ? alpha : wIn_ % alpha );
				k    = ( wIn_ % alpha == 0 ? wIn_ / alpha : ceil ( double ( wIn_ ) / double ( alpha ) ) );
			}
			
	}
	

	
	
/*	
	void IntAdderClassical::emulate ( TestCase* tc ) {
		mpz_class svX = tc->getInputValue ( "X" );
		mpz_class svY = tc->getInputValue ( "Y" );
		mpz_class svC = tc->getInputValue ( "Cin" );
		
		mpz_class svR = svX + svY + svC;
		// Don't allow overflow
		mpz_clrbit ( svR.get_mpz_t(),wIn_ );
		
		tc->addExpectedOutput ( "R", svR );
	}*/
	
	
}


