/*
 * An integer multiplier for FloPoCo using the Karatsuba method
 *
 * Authors : Bogdan Pasca
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"
#include "Karatsuba.hpp"

using namespace std;
extern vector<Operator*> oplist;

Karatsuba:: Karatsuba(Target* target, int wInX, int wInY) :
	Operator(target), wInX_(wInX), wInY_(wInY), wOut_(wInX_ + wInY_){

	int depth, i, j;
	ostringstream name, nameH, nameL;
	setOperatorType();
	
	/* Set up the IO signals
	  * X and Y have wInX_ and wInY_ bits respectively 
	  * R has wOut_ bits where wOut_ = (wInX_ + wInY_) bits
	  */
	addInput ("X", wInX_);
	addInput ("Y", wInY_);
	addOutput("R", wOut_);

	/* set parameters */

	//for the current version of the code, the variable multiplierYWidth will be discarded, considering that the multipliers are symetrical
	int multiplierYWidth; 
  
  /* - gather in multiplierXWidth_ and multiplierYWidth the maximum multiplier dimmensions so that the required ferquency is reached
	   - the status variable reports if the required frequency can be reached. If the value is false, the multiplierXWidth_ and multiplierYWidth will return the maximum
	     the values for the best reachable frequency */   
	bool status = target->suggestSubmultSize(multiplierXWidth_, multiplierYWidth, wInX_, wInY_); 
	if (!status)
		cout<<"Warning: the requested frequency cannot be reached. "<<endl;

	if (verbose)
	cout<<"The suggested width for the multiplier is: "<<  multiplierXWidth_<<endl;

	if (isSequential()){
		//TODO sequential version
	}
	else{
		//The signal declarations which belong here have been moved to the recursive method which builds the architecture
	}
}

Karatsuba::~Karatsuba() {
}

void Karatsuba::setOperatorName(){
	/* Name Setup procedure
	 *  The name has the format: Karatsuba_wInX__wInY_
	 *  wInX_ = width of the X input
	 *  wInY_ = width of the Y input
	 */  
	ostringstream name;  
	name.str("");;
	name <<"Karatsuba_"<<wInX_<<"_"<<wInY_;
	uniqueName_ = name.str(); 
}

void Karatsuba::outputVHDL(std::ostream& o, std::string name) {
	//temporary stream. When the method BuildCombinationalKaratsuba will return, this stream will contain the contents of the combinational version for the Karatsuba
	// architecture (all lines of code between "begin" and "end architecture"
	ostringstream tempCombinationalKaratsubaStream;
  
	licence(o,"Bogdan Pasca (2008)");
	Operator::stdLibs(o);  
	outputVHDLEntity(o);

 
	if (isSequential()){
		//TODO - the sequential version is due
	}	
	else{
	  // the combinational version for the archtecture code will be outputed in the tempCombinationalKaratsubaStream stream. 
	  // the intermediary signals are also added to the signal lists within this recursive method 
		BuildCombinationalKaratsuba( tempCombinationalKaratsubaStream ,wInX_-1,0,"X","Y","R", 0,"mb");
	}
	
	newArchitecture(o,name);
	//the signal lists now contain all the signals needed by the architecture.
	outputVHDLSignalDeclarations(o);
	
	beginArchitecture(o);
	
	if (isSequential()){
		//TODO - the sequential version is due
	}	
	else{
		//the main stream is feed with the architecture code, wich was previously obtain using the BuildCombinationalKaratsuba method
		o<<tempCombinationalKaratsubaStream.str();
	}
		
	endArchitecture(o);
}

void Karatsuba::BuildCombinationalKaratsuba(std::ostream& o, int L, int R, string lName, string rName, string resultName, int depth, string branch)
{
	int termWidth; // will contain the width of the multiplication terms. This width will be computed as L - R + 1 
	int basisBits; // the width of the determined basis for this level of recursion. The number of bits of the basis are computed as floor(termWidth/2).
	ostringstream t1_LeftTerm, t1_RightTerm, t1, t0, t2, sum_t0_t2, t1_PostSubstractions, sum_shifted_t2_t0;

  // the recursivity ends when the width of the multiplication terms becames <= multiplierXWidth_ 
  if (L-R+1 <= multiplierXWidth_)
    o << resultName << " <= "<< lName <<"("<<L<<" downto "<<R<<") * "<<
                                rName <<"("<<L<<" downto "<<R<<");"<<endl;
  else{ 
		/* parameters setup */
		termWidth = L - R + 1; 
		if (termWidth % 2 == 0)
			basisBits = termWidth / 2;
		else
			basisBits = termWidth / 2 + 1;
    
    /* 
    Let X = x1 * b^k + x0
        Y = y1 * b^k + y0, where x0,x1,y0,y1<b^K
        
        Z = X * Y = (x1 * b^k + x0) * (y1 * b^k + y0) = x1*y1 *b^2k + (x1*y0 + x0*y1)*b^k + x0*y0
        
        Instead of 4 multiplications: x1*y1, x1*y0, x0*y1 and x0*y0, only 3 multiplications are computed.
        Denote by: 
        t0 = x0*y0
        t1 = (x0 + x1)*(y0 + y1)
        t2 = x1*y1
        
        The product Z can be written as:
        Z = t2 * b^2k + (t1 - t0 - t2)*2^k + t0
        
        (x0 + x1)      -> t1_LeftTerm
        (y0 + y1)      -> t1_RightTerm 
        (t1 - t0 - t2) -> t1_PostSubstractions
        t0 + t2        -> sum_t0_t2
        t2 * b^2k + t0 -> sum_shifted_t2_t0
        
    */
    
    
    //Initialize the streams which hold the name of the signals for the current level of recursion    
		t1_LeftTerm             <<"t1_LeftTerm_"          <<depth<<"_"<<branch;
		t1_RightTerm            <<"t1_RightTerm_"         <<depth<<"_"<<branch;
		t1                      <<"t1_"                   <<depth<<"_"<<branch;
		t0                      <<"t0_"                   <<depth<<"_"<<branch;
		t2                      <<"t2_"                   <<depth<<"_"<<branch;
		sum_t0_t2               <<"sum_t0_t2_"            <<depth<<"_"<<branch;
		t1_PostSubstractions    <<"t1_PostSubstractions_" <<depth<<"_"<<branch;
		sum_shifted_t2_t0       <<"sum_shifted_t2_t0_"    <<depth<<"_"<<branch;

		
		//add the signals corresponding to the terms of the t1 multiplication
		addSignal( t1_LeftTerm.str() , basisBits+1); 
		addSignal( t1_RightTerm.str(), basisBits+1);
		
		o<<  t1_LeftTerm.str() <<" <= (\"0"<<  ((termWidth % 2 == 0)?"":"0") <<"\" & "<<lName<<"("<< L            <<" downto "<< L - basisBits +((termWidth % 2 == 0)?1:2)<<")) + "
		                       <<"(\"0\" & "<<lName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;    
		o<< t1_RightTerm.str() <<" <= (\"0"<<  ((termWidth % 2 == 0)?"":"0") <<"\" & "<<rName<<"("<< L            <<" downto "<< L - basisBits +((termWidth % 2 == 0)?1:2)<<")) + "
		                       <<"(\"0\" & "<<rName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;

		addSignal(t1.str(),2*(basisBits+1));
	
		//add the signals for the products t0 and t2
		addSignal(t0.str(),2* basisBits);
		addSignal(t2.str(),2* ((termWidth % 2 == 0)?basisBits:(basisBits-1))  );

		BuildCombinationalKaratsuba(o, basisBits, 0 , t1_LeftTerm.str(), t1_RightTerm.str(), t1.str(), depth+1, branch + "_cb");
		BuildCombinationalKaratsuba(o, L-basisBits+((termWidth % 2 == 0)?1:2)-1, L-2*basisBits+((termWidth % 2 == 0)?1:2), lName , rName, t0.str(), depth+1, branch + "_rb");
		BuildCombinationalKaratsuba(o, L , L-basisBits + ((termWidth % 2 == 0)?1:2) , lName, rName, t2.str(), depth+1, branch + "_lb");


		addSignal(sum_t0_t2.str()           ,2* basisBits+1 );
		addSignal(t1_PostSubstractions.str(),2*(basisBits+1));
		addSignal(sum_shifted_t2_t0.str()   ,2*(L-R+1)      ); 
		  
		o<< sum_t0_t2.str()            <<" <= "<<"( \"0\" & "<<t0.str() << ") + ( \"0\" & "<< t2.str()<<");"<<endl;
		o<< t1_PostSubstractions.str() <<" <= "<< t1.str()              <<  " - "<<"(\"0\" & "<< sum_t0_t2.str()<<")"<<";"<<endl;

		o<< sum_shifted_t2_t0.str()    <<" <= "<< "("<< t2.str()                            <<" & "<< zeroGenerator(2*basisBits,0)<<") + "
		                                       << "("<< zeroGenerator( 2*(L-R+1) - 2*basisBits,0)<<" & "<< t0.str()               <<");"<<endl;
		o<< resultName             <<" <= "<< sum_shifted_t2_t0.str() << " + "
		                                << "("<<zeroGenerator(2*(L-R+1)-(3*basisBits+2),0)<<" & "<< t1_PostSubstractions.str()<<" & "<<zeroGenerator(basisBits,0)<<");"<<endl;   
	}
}                            

TestIOMap Karatsuba::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*getSignalByName("X"));
	tim.add(*getSignalByName("Y"));
	tim.add(*getSignalByName("R"));
	return tim;
}

void Karatsuba::fillTestCase(mpz_class a[])
{
	mpz_class& x = a[0];
	mpz_class& y = a[1];
	mpz_class& r = a[2];
	r = x * y;
}

