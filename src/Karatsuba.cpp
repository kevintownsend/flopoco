/*
 * An integer multiplier for FloPoCo
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

/** 
 * The constructor of the Karatsuba class
 * @param target argument of type Target containing the data for which this operator will be optimized
 * @param wInX integer argument representing the width in bits of the input X 
 * @param wInY integer argument representing the width in bits of the input Y
 **/
Karatsuba:: Karatsuba(Target* target, int wInX, int wInY) :
	Operator(target), wInX(wInX), wInY(wInY), wOut(wInX + wInY){

	int depth, i, j;
	ostringstream name, nameH, nameL;

	//gets information if the multiplier is or not to be pipelined and sets the corresponding operator attribute
	if(target->is_pipelined())
	  set_sequential();
	else
	  set_combinatorial();  

	/* Name Setup procedure
	  *  The name has the format: Karatsuba_wInX_wInY
	  *  wInX = width of the X input
	  *  wInY = width of the Y input
	  */  
	name.str("");;
	name <<"Karatsuba_"<<wInX<<"_"<<wInY;
	unique_name = name.str(); 

	/* Set up the IO signals
	  * X and Y have wInX and wInY bits respectively 
	  * R has wOut bits where wOut = (wInX + WInY) bits
	  */
	add_input ("X", wInX);
	add_input ("Y", wInY);
	add_output("R", wOut);

	/* set parameters */

	//for the current version of the code, the variable multiplierYWidth will be discarded, considering that the multipliers are symetrical
	int multiplierYWidth; 
  
  /* - gather in multiplierXWidth and multiplierYWidth the maximum multiplier dimmensions so that the required ferquency is reached
	   - the status variable reports if the required frequency can be reached. If the value is false, the multiplierXWidth and multiplierYWidth will return the maximum
	     the values for the best reachable frequency */   
	bool status = target->suggest_submult_size(multiplierXWidth, multiplierYWidth, wInX, wInY); 
	if (!status)
		cout<<"Warning: the requested frequency cannot be reached. "<<endl;

	if (verbose)
	cout<<"The suggested width for the multiplier is: "<<  multiplierXWidth<<endl;

	if (is_sequential()){
		//TODO sequential version
	}
	else{
		//The signal declarations which belong here have been moved to the recursive method which builds the architecture
	}
}



/**
 * Karatsuba destructor
 */
Karatsuba::~Karatsuba() {
}

/**
 * Method belonging to the Operator class overloaded by the Karatsuba class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/
void Karatsuba::output_vhdl(std::ostream& o, std::string name) {
	
	//temporary stream. When the method BuildCombinationalKaratsuba will return, this stream will contain the contents of the combinational version for the Karatsuba
	// architecture (all lines of code between "begin" and "end architecture"
	ostringstream tempCombinationalKaratsubaStream;
  
	Licence(o,"Bogdan Pasca (2008)");
	Operator::StdLibs(o);  
	output_vhdl_entity(o);

 
	if (is_sequential()){
		//TODO - the sequential version is due
	}	
	else{
	  // the combinational version for the archtecture code will be outputed in the tempCombinationalKaratsubaStream stream. 
	  // the intermediary signals are also added to the signal lists within this recursive method 
		BuildCombinationalKaratsuba( tempCombinationalKaratsubaStream ,wInX-1,0,"X","Y","R", 0,"mb");
	}
	
	new_architecture(o,name);
	//the signal lists now contain all the signals needed by the architecture.
	output_vhdl_signal_declarations(o);
	
	begin_architecture(o);
	
	if (is_sequential()){
		//TODO - the sequential version is due
	}	
	else{
		//the main stream is feed with the architecture code, wich was previously obtain using the BuildCombinationalKaratsuba method
		o<<tempCombinationalKaratsubaStream.str();
	}
		
	end_architecture(o);
}

/**
 * Method which generates the code for the combinational version of of the Karatsuba architecture
 * @param[in,out] o           the stream where the code generated by the method is returned to
 * @param[in]     L           the left index of the lName and rName arrays which will be considered as the MSB of the inputs
 * @param[in]     R           the right index of the lName and rName arrays which will be considered as the LSB of the inputs. 
 *                            For example, if lName="X" and rName="Y", this method will compute the product X(L downto R)*Y(L downto R)
 * @param[in]     lName       the name of the left term of the product
 * @param[in]     rName       the name of the right term of the product
 * @param[in]     resultName  the name of the product. Adding a signal with this name and the correct size is a responasbility of the calling function.
 * @param[in]     depth       the recursivity depth
 * @param[in]     branch      the branch we are currently on inside the recursivity. 
 *                            The main branch is named mb, the central branch cb, left branch lb and right branch rb
 **/
void Karatsuba::BuildCombinationalKaratsuba(std::ostream& o, int L, int R, string lName, string rName, string resultName, int depth, string branch)
{
// will contain the width of the multiplication terms. This width will be computed as L - R + 1 
int termWidth;
// the width of the determined basis for this level of recursion. 
// The number of bits of the basis are computed as floor(termWidth/2).
int basisBits;


ostringstream t1_LeftTerm, t1_RightTerm, t1, t0, t2, sum_t0_t2, t1_PostSubstractions, sum_shifted_t2_t0;

  // the recursivity ends when the width of the multiplication terms becames <= multiplierXWidth 
  if (L-R+1 <= multiplierXWidth)
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
		add_signal( t1_LeftTerm.str() , basisBits+1); 
		add_signal( t1_RightTerm.str(), basisBits+1);
		
		o<<  t1_LeftTerm.str() <<" <= (\"0"<<  ((termWidth % 2 == 0)?"":"0") <<"\" & "<<lName<<"("<< L            <<" downto "<< L - basisBits +((termWidth % 2 == 0)?1:2)<<")) + "
		                       <<"(\"0\" & "<<lName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;    
		o<< t1_RightTerm.str() <<" <= (\"0"<<  ((termWidth % 2 == 0)?"":"0") <<"\" & "<<rName<<"("<< L            <<" downto "<< L - basisBits +((termWidth % 2 == 0)?1:2)<<")) + "
		                       <<"(\"0\" & "<<rName<<"("<<L - basisBits +((termWidth % 2 == 0)?1:2) -1 <<" downto "<<L -2*basisBits +((termWidth % 2 == 0)?1:2)<<"));"<<endl;

		add_signal(t1.str(),2*(basisBits+1));
	
		//add the signals for the products t0 and t2
		add_signal(t0.str(),2* basisBits);
		add_signal(t2.str(),2* ((termWidth % 2 == 0)?basisBits:(basisBits-1))  );

		BuildCombinationalKaratsuba(o, basisBits, 0 , t1_LeftTerm.str(), t1_RightTerm.str(), t1.str(), depth+1, branch + "_cb");
		BuildCombinationalKaratsuba(o, L-basisBits+((termWidth % 2 == 0)?1:2)-1, L-2*basisBits+((termWidth % 2 == 0)?1:2), lName , rName, t0.str(), depth+1, branch + "_rb");
		BuildCombinationalKaratsuba(o, L , L-basisBits + ((termWidth % 2 == 0)?1:2) , lName, rName, t2.str(), depth+1, branch + "_lb");


		add_signal(sum_t0_t2.str()           ,2* basisBits+1 );
		add_signal(t1_PostSubstractions.str(),2*(basisBits+1));
		add_signal(sum_shifted_t2_t0.str()   ,2*(L-R+1)      ); 
		  
		o<< sum_t0_t2.str()            <<" <= "<<"( \"0\" & "<<t0.str() << ") + ( \"0\" & "<< t2.str()<<");"<<endl;
		o<< t1_PostSubstractions.str() <<" <= "<< t1.str()              <<  " - "<<"(\"0\" & "<< sum_t0_t2.str()<<")"<<";"<<endl;

		o<< sum_shifted_t2_t0.str()    <<" <= "<< "("<< t2.str()                            <<" & "<< zero_generator(2*basisBits,0)<<") + "
		                                       << "("<< zero_generator( 2*(L-R+1) - 2*basisBits,0)<<" & "<< t0.str()               <<");"<<endl;
		o<< resultName             <<" <= "<< sum_shifted_t2_t0.str() << " + "
		                                << "("<<zero_generator(2*(L-R+1)-(3*basisBits+2),0)<<" & "<< t1_PostSubstractions.str()<<" & "<<zero_generator(basisBits,0)<<");"<<endl;   
	}
}                            




/**
 * A zero generator method which takes as input two arguments and returns a string of zeros with quotes as stated by the second argurment
 * @param[in] n        integer argument representing the number of zeros on the output string
 * @param[in] margins  integer argument determining the position of the quotes in the output string. The options are: -2= no quotes; -1=left quote; 0=both quotes 1=right quote
 * @return returns a string of zeros with the corresonding quotes given by margins
 **/
string Karatsuba::zero_generator(int n, int margins)
{
	ostringstream left,full, right, zeros;
	int i;

	//generate the zeros without quotes
	for (i=1; i<=n;i++)
		zeros<<"0";

	//generate the 3 strings with left quote, both and right quote
	left<<"\""<<zeros.str();
	full<<left.str()<<"\"";
	right<<zeros.str()<<"\"";

	//margins determines which of the strings are outputed
	switch(margins){
		case -2: return zeros.str(); break;
		case -1: return left.str();  break;
		case  0: return full.str();  break;
		case  1: return right.str(); break;
		default: return full.str();
	}
}

TestCaseList Karatsuba::generateRandomTestCases(int n)
{
	// TODO
	/* Signals */
	Signal sx = *get_signal_by_name("X");
	Signal sy = *get_signal_by_name("Y");
	Signal sr = *get_signal_by_name("R");

	TestCaseList tcl;	/* XXX: Just like Lyon's Transporation company. :D */
	mpz_class x, y, r;

	for (int i = 0; i < n; i++)	
	{
		x = getLargeRandom(sx.width());
		y = getLargeRandom(sy.width());
		r = x * y;

		TestCase tc;
		tc.addInput(sx, x);
		tc.addInput(sy, y);
		tc.addExpectedOutput(sr, r);
		tc.addComment(x.get_str() + " * " + y.get_str() + " = " + r.get_str());
		tcl.add(tc);
	}
	return tcl;
}

TestCaseList Karatsuba::generateStandardTestCases(int n)
{
	// TODO
	return TestCaseList();
}
