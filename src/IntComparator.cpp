/*
 * An integer comparator unit
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
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntComparator.hpp"

using namespace std;



namespace flopoco{


	IntComparator::IntComparator(Target* target, int wIn, int criteria, bool constant, int constValue, map<string, double> inputDelays) :
		Operator(target, inputDelays), wIn_(wIn), criteria_(criteria) {
	
		// -------- Parameter set up -----------------
		srcFileName = "IntComaprator";
		if (criteria<-2 && criteria>2)
			criteria_=0;
		
		ostringstream name;
		switch(criteria){
			case -2: name << "IntComparator_"<< wIn<<"_"<<"less";   break;
			case -1: name << "IntComparator_"<< wIn<<"_"<<"leq";    break;
			case  0: name << "IntComparator_"<< wIn<<"_"<<"eq";     break;
			case  1: name << "IntComparator_"<< wIn<<"_"<<"geq";    break;
			case  2: name << "IntComparator_"<< wIn<<"_"<<"greater";break;
			default: name << name << "IntComparator_"<< wIn<<"_"<<"eq";
		}
		setName(name.str());
		setCopyrightString("Bogdan Pasca (2010)");

		addInput ("X", wIn_);
		addInput ("Y", wIn_);
		addOutput("R"); 

		switch(criteria_){
			case -2: vhdl << tab << "R <= '1' when (X <  Y) else '0';"<<endl; break;
			case -1: vhdl << tab << "R <= '1' when (X <= Y) else '0';"<<endl; break;
			case  0:{ 
				
				//equality is implemented as tree of luts
				//how many luts we need for the first level?
				int l = target->lutInputs();
				int nLUTsStage = (wIn % l/2 == 0? wIn/ (l/2): wIn/ (l/2) + 1);
				int numberOfLevels = 1;

				REPORT(INFO, "number of luts for first level is " << nLUTsStage);
				
				int sizeFirstStage = (l/2)*nLUTsStage;
				vhdl << tab << declare("eX",sizeFirstStage)<< " <= " << zg(sizeFirstStage - wIn) << " & X;"<<endl;
				vhdl << tab << declare("eY",sizeFirstStage)<< " <= " << zg(sizeFirstStage - wIn) << " & Y;"<<endl;
				
				setCriticalPath(0.0);
				manageCriticalPath(target->localWireDelay() + target->lutDelay());
				//first stage
				for (int i=0; i<nLUTsStage;i++){
					vhdl << tab << declare( join("level0_LUT",i) ) << "<= '1' when eX"<<range((l/2)*(i+1)-1,(l/2)*i)<<"=eY"<<range((l/2)*(i+1)-1,(l/2)*i)<<"else '0';"<<endl;
				}
				
				vhdl << tab << declare("Xlevel0",nLUTsStage) << " <= ";
				for (int i=0; i<nLUTsStage;i++){
					if (i==nLUTsStage-1)
						vhdl << " level0_LUT"<<i<<";"<<endl;
					else
						vhdl << " level0_LUT"<<i<<" & ";
				}
				
				
				while (nLUTsStage != 1){
					int inputBits  = nLUTsStage;
					int outputBits = 0;
					int idx        = 0;
					int eaten;
					
					manageCriticalPath(target->localWireDelay() + target->lutDelay());
					while (inputBits >  0 ) {
						if (inputBits >= l){
							eaten = l;
							vhdl << tab << declare( join("level",numberOfLevels,"_LUT",outputBits) ) 
							<< "<= '1' when Xlevel"<<numberOfLevels-1<<range(idx+l-1,idx)<<" = not("<<zg(l,0)<<") else '0';"<<endl;
						}else{
							eaten = inputBits;
							vhdl << tab << declare( join("level",numberOfLevels,"_LUT",outputBits) )
							<< "<= '1' when Xlevel"<<numberOfLevels-1<<range(idx+eaten-1,idx)<<"=not("<<zg(eaten,0)<< ") else '0';"<<endl;
						}
							
						idx+=eaten;
						outputBits++;
						inputBits-=eaten; 
					}
					
					//create the new X
					vhdl << tab << declare(join("Xlevel",numberOfLevels),outputBits) << " <= ";
					for (int i=0; i<outputBits;i++){
						if (i==outputBits-1)
							vhdl << " level"<<numberOfLevels<<"_LUT"<<i<<";"<<endl;
						else
							vhdl << " level"<<numberOfLevels<<"_LUT"<<i<<" & ";
					}
					
					
					nLUTsStage = outputBits; //reinit for next iteration
					numberOfLevels++;
					
				}
				
				numberOfLevels--;
				REPORT(INFO, "number of levels was " << numberOfLevels);
				
				vhdl << tab << "R <= Xlevel"<<numberOfLevels<<";" <<endl; break;
			}
			case  1: vhdl << tab << "R <= '1' when (X >  Y) else '0';"<<endl; break;
			case  2: vhdl << tab << "R <= '1' when (X = Y) else '0';"<<endl; break;
			default:;
		}
		
		
		}

	IntComparator::~IntComparator() {
	}
	
	void IntComparator::emulate(TestCase* tc){
		mpz_class svX = tc->getInputValue ( "X" );
		mpz_class svY = tc->getInputValue ( "Y" );
		
		mpz_class svR;
		
		switch(criteria_){
			case -2: if ( svX <  svY ) svR = 1; else svR = 0; break;
			case -1: if ( svX <= svY ) svR = 1; else svR = 0; break;
			case  0: if ( svX == svY ) svR = 1; else svR = 0; break;
			case  1: if ( svX >= svY ) svR = 1; else svR = 0; break;
			case  2: if ( svX > svY  ) svR = 1; else svR = 0; break;
			default:; 
		}
		
		tc->addExpectedOutput ( "R", svR );
	}


}
