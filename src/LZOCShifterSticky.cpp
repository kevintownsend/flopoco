/*
 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
 *
 * Authors : Florent de Dinechin, Bogdan Pasca
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

#include "LZOCShifterSticky.hpp"

using namespace std;

//The width of the lzo output will be floor(log2(wIn+1)). It should be accessed as operator->getCountWidth().
// TODO to be optimal for FPAdder, we have to provide a way to disable the sticky computation.

/** 
 * The LZOCShifterSticky constructor
 * @param[in] target the target device for this operator
 * @param[in] wIn the width of the mantissa input
 * @param[in] wOut the width of the mantissa output
 */
LZOCShifterSticky::LZOCShifterSticky(Target* target, int wIn, int wOut, bool compute_sticky, const int countType) :
	Operator(target), wIn(wIn), wOut(wOut), computeSticky(compute_sticky), countType(countType) {

	entityType=(countType<0?generic:specific);
	
	int i, p2i;

	Operator::set_operator_name();
	set_operator_type();


	/* Set up the internal architecture signals */

	/* Terminology: 
	   - level i is the signal on which we will test wether the 2^(i-1) leading bits are zeroes.
	   - i goes from wCount to 0 */
	
	/*  from the end of the pipeline to the beginning, 
	    determine the sizes of all register levels */
	i=0;
	size[0] = wOut; /* size of the result */
	while(size[i]-wOut < wIn){
		i++;
		p2i = 1<<(i-1);
		size[i] = size[i-1] + p2i;
		/* Invariant: size[i] = wOut + 2^i -1 */
	}
	/* the attribute that gives the number of bits of the LZO count */
	wCount = i;
	/* the first stage doesn't need to register zeroes */
	size[wCount] = 1<<wCount;
	
	/* should be identical to : wCount = intlog2(wIn+1); 
	   +1 for the case all zeroes */
	
	/* Set up the IO signals */
		
	si = add_input ("I", wIn);
		
	add_output("Count", wCount);
	add_output("O", wOut);
	/* if we generate a generic LZOC */
	if (entityType==generic) 
		add_input ("OZb"); 
	/* if we require a sticky bit computation */
	if(compute_sticky) 
		add_output("Sticky"); 
	
	
	if(verbose){
		cout <<endl<<"  wCount=" << wCount << "    Level sizes: ";
		for(i=0; i<=wCount; i++) 
			cout << size[i]<<" ";
		cout <<endl;
	}

	double critical_path = 0.0;
	bool registered;

	for (int i=wCount; i>=0; i--){
		ostringstream levelName, leveldName, stickyName;
		levelName << "level"  << i;
		stickyName << "sticky" << i;
		level[i] = levelName.str();
		double stage_delay = 1.2 * target->local_wire_delay() * (1<<i);
		
		if (i==wCount)
		{
			critical_path=stage_delay;
			level_registered[i] = false;
		}	
		else
			if (critical_path + stage_delay > 1/target->frequency()) {
				critical_path=stage_delay; 	// reset critical path
				level_registered[i] = true;
				increment_pipeline_depth();
			}
	}
	
	if (is_sequential())
	{
		for (int i=wCount; i>=0; i--){
			for (int j=0; j<=i; j++)
				if (level_registered[j])
					countDepth[i]++;	
			if (level_registered[i]){
				leveld[i] = level[i] + "_d" ;
				ostringstream levelName;
				levelName << "level"<<i;		
				add_registered_signal(levelName.str(), size[i]);	
			}
			else{
				leveld[i] = level[i];
				ostringstream levelName;
				levelName << "level"<<i;		
				add_signal(levelName.str(), size[i]);	
			}
		}
	
		for (int i=wCount-1; i>=0; i--){
			ostringstream countName;
			countName << "count"<<i;		
			add_delay_signal(countName.str(), 1, countDepth[i]);
		}
			
		if(compute_sticky){
			for (int j=wCount; j>=0; j--){
				ostringstream stickyName;
				ostringstream eqVer;
				ostringstream newSticky;
				eqVer<<"eqVer"<<j;
				newSticky<<"newSticky"<<j;
				stickyName<<"sticky"<<j;
				
				if (level_registered[j])			
					add_registered_signal(stickyName.str());
				else
					add_signal(stickyName.str());
			
				add_signal(eqVer.str());
				add_signal(newSticky.str());	
			}
		}						
			
		if (entityType==generic)
			add_delay_signal_no_reset("sozb",1, pipeline_depth());

	}else
	{
		set_pipeline_depth(0);
		for (int i=wCount; i>=0; i--){
				leveld[i] = level[i];
				ostringstream levelName;
				levelName << "level"<<i;		
				add_signal(levelName.str(), size[i]);
				if(compute_sticky){
					ostringstream stickyName;
					ostringstream eqVer;
					ostringstream newSticky;
					stickyName<<"sticky"<<i;
					eqVer<<"eqVer"<<i;
					newSticky<<"newSticky"<<i;
					add_signal(stickyName.str());
					add_signal(eqVer.str());
					add_signal(newSticky.str());
				}					
		}
		
		for (int i=wCount-1; i>=0; i--){
			ostringstream countName;
			countName << "count"<<i;		
			add_signal(countName.str(), 1);
		}
	
		add_signal("sozb",1);
		
	}
	add_signal("preCount", wCount);	
	
	
		/*	
		if (is_sequential() && i!=wCount) 
			add_registered_signal(level[i], size[i]);
			leveld[i] = level[i] + "_d" ;
			if(compute_sticky)
				add_registered_signal(stickyName.str());
		}
		else {
			add_signal(           level[i], size[i]);
			leveld[i] = level[i];
			if(compute_sticky)
				add_signal(           stickyName.str());
		}
		// The signal that holds the leading zero count
		if(i<wCount) {
			ostringstream name;
			name << "count" << i;
			if (is_sequential()) 
				add_delay_signal_no_reset(name.str(),1, i+1);				
			else 
				add_signal(name.str());
		}
		*/
	
	//}
/*
	if (is_sequential()){
		//set_pipeline_depth(wCount); //was wCount-1
		if (entityType==generic)
		add_delay_signal_no_reset("sozb",1, wCount-1);
	}
*/
}

/** The LZOCShifterSticky destructor
*/
LZOCShifterSticky::~LZOCShifterSticky() {
}


/* Accesor methods */
/* ===========================================================================*/

/** Returns the number of bits of the count
 * @return the number of bits of the count
 */
int LZOCShifterSticky::getCountWidth() const{
	return wCount;
}


/* Overriden methods */
/* ===========================================================================*/

/** Method for setting the operator name
	* @param[in] prefix the prefix that is going to be placed to the default name
	*                   of the operator
	* @param[out] postfix the postfix that is going to be placed to the default
	*                     name of the operator
	*/
void LZOCShifterSticky::set_operator_name(std::string prefix, std::string postfix){
	ostringstream name; 
	ostringstream computationInterface;
	/* The computationInterface refers to the VHDL entity ports. 
	   If computation is generic, then an extra input port is 
	   available on the entity for specifying the count type */ 
	if (countType < 0) 
		computationInterface<<"gen";
	else 
		computationInterface<<"spec";
	
	name<<prefix
	    <<"LZOCShifterSticky_"<<wIn<<"_"<<wOut<<"_"<<computeSticky<<"_" //semantic name
	    <<computationInterface.str()
	    <<postfix;
	
	unique_name=name.str();
}


/**
 * Method belonging to the Operator class overloaded by the LZOCShifterSticky class
 * @param[in,out] o     the stream where the current architecture will be outputed to
 * @param[in]     name  the name of the entity corresponding to the architecture generated in this method
 **/

void LZOCShifterSticky::output_vhdl(std::ostream& o, std::string name) {
	Licence(o,"Florent de Dinechin, Bogdan Pasca (2007)");
	Operator::StdLibs(o);
	output_vhdl_entity(o);
	new_architecture(o,name);
	output_vhdl_signal_declarations(o);	
	begin_architecture(o);
	
	/* connect input to level wCount, possibly with 0 padding */
	o << tab << "level" << wCount<< " <= " ;
	if(wIn == size[wCount]) // no padding needed
		o << "I;" <<endl;
	else 
		o << "I & (" << size[wCount]-wIn-1 << " downto 0 => '0');" << endl ;

	if (entityType==generic)
		o<<tab<<"sozb <= OZb;"<<endl;
	
	if(computeSticky)
		o<<tab<<"sticky" << wCount << " <= '0';"<<endl;

	for  (int i = wCount; i>=1; i--){
		int p2i = 1 << i;
		int p2io2 = 1 << (i-1);

		//=====================================
		o << tab << "count" << i-1 << " <= '1' when " << leveld[i] << "(" << size[i]-1 << " downto " << size[i]-p2io2 << ") = (" << p2io2-1 << " downto 0 => "; 
		if (entityType==generic)
			o<<get_delay_signal_name("sozb",pipeline_depth()-countDepth[i-1]);
		else
			o<<"'"<<countType<<"'";
		o<<" )   else '0';" << endl; 
		//=====================================
		// REM in the following,  size[i]-size[i-1]-1 is almost always equal to p2io2, except in the first stage
		if (i==wCount){
			o << tab << level[i-1] << " <= " ;
			o << " " << "("<<leveld[i] << "(" << size[i]/2 -1 << " downto " << 0 << ") & "<<zero_generator(size[i-1]-size[i]/2, 0)<<" )  when count" << i-1 << "='1'" << endl;
			
			if (size[i-1]-size[i]>0)
				o << tab << tab << tab <<  "else (" << leveld[i] << "(" << size[i]-1 << " downto 0) &  "<<zero_generator(size[i-1]-size[i],0)<<") ;" << endl; 
			else
				o << tab << tab << tab <<  "else " << leveld[i] << "(" << size[i]-1 << " downto "<<size[i]-size[i-1]<<");" << endl; 
		}
		else
		{
			o << tab << level[i-1] << " <= " ;
			o << " " << leveld[i] << "(" << size[i]-1 << " downto " << size[i]-size[i-1] << ")   when count" << i-1 << "='0'" << endl;
			o << tab << tab << tab <<  "else " << leveld[i] << "(" << size[i-1]-1 << " downto 0);" << endl; 
		}
		
		if(computeSticky){
			ostringstream eqVer;
			ostringstream newSticky;
			eqVer << "eqVer"<<i;
			newSticky <<"newSticky"<<i;
			
			o << tab << eqVer.str() <<" <= '1' when " << leveld[i] << "(" << size[i]-size[i-1]-1 << " downto 0) = (" << size[i]-size[i-1]-1 << " downto 0 => '0') else '0';"
			  << endl;
			o << tab << newSticky.str() <<" <= sticky" << i << " or not("<<eqVer.str()<<");"<<endl;  
		
			if (size[i]-size[i-1]-1>=0){
				o << tab << "sticky" << i-1 << " <= sticky"<<i;
				if(is_sequential() && i<wCount && level_registered[i]) 
					o << "_d";
				o <<" when (count" << i-1 << "='1') else "<<newSticky.str();
			}
			else{
				o << tab << "sticky" << i-1 << " <= sticky"<<i<<" when count" << i-1 << "='1'  else sticky" << i;	
				if(is_sequential() && i<wCount && level_registered[i]) 
					o << "_d";
			}
			o << ";" << endl;
		}
			
		o << endl;
	}

	if(computeSticky)
		if (level_registered[0])
			o << tab << "sticky <= sticky0_d;" << endl;
		else
			o << tab << "sticky <= sticky0;" << endl;
			
			
	//o << tab << "O      <= level0_d;" << endl;
	//leveld[i]
	o << tab << "O      <= "<<leveld[0]<<";" << endl;
	o << tab << "preCount  <= ";
	for (int i=wCount-1; i >=0; i--){
		ostringstream name;
		name << "count" << i;
		if (is_sequential())
			o << get_delay_signal_name(name.str(), countDepth[i]);
		else
			o << name.str();
		if (i>0) o << " & ";
	}
	
	o << ";" << endl;

	ostringstream binInputWidth;
	
	printBinNum(binInputWidth, wIn, wCount);	

	o << tab << "Count <= \""<<binInputWidth.str()<<"\" when preCount = ("<<wCount-1<<" downto  0 =>";
//	if (entityType==specific)
		o <<"'1')";
//	else
//		o<<"not("<<get_delay_signal_name("sozb",pipeline_depth())<<"))";
	
	o<<" else preCount;"<<endl;
	 
	if (is_sequential()){		
		output_vhdl_registers(o); o<<endl;
	}

	end_architecture(o);
	
}







TestIOMap LZOCShifterSticky::getTestIOMap()
{
	TestIOMap tim;
	tim.add(*get_signal_by_name("I"));
	//tim.add(*get_signal_by_name("OZb"));
	tim.add(*get_signal_by_name("Count"));
	tim.add(*get_signal_by_name("O"));
	if (computeSticky)
		tim.add(*get_signal_by_name("Sticky"));
	
	return tim;
}





void LZOCShifterSticky::fillTestCase(mpz_class a[])
{

	if (entityType==specific)
	{
		mpz_class& si     = a[0];
		mpz_class& scount = a[1];
		mpz_class& so     = a[2];
		mpz_class& ssticky = a[3];
		
		//if (computeSticky)
		//	ssticky = a[3];
		
		int sticky=0;
	
		/* Count the leading zero/one s */
		int j=(wIn-1);//the index of the MSB of the input
		int bit = (countType == 0) ? 0 : 1; //what are we counting
		for (j = (wIn-1); j >= 0; j--)
				if (mpz_tstbit(si.get_mpz_t(), j) != bit)
				break;
		
		int icount =(wIn-1)-j;
		scount = icount; //the count result
	
		//compute max value on wOut bits
		maxValue=2;
		for (int i=2;i<=wOut;i++)
		maxValue=maxValue*2;
		maxValue--;
	
		mpz_class inputValue;
		inputValue=si;
		sticky=0;

		if (countType==0) /* if we are counting zeros */
			if (inputValue!=0) 
				while (!((inputValue<=maxValue)&&(2*inputValue>maxValue)))
					if (inputValue>maxValue){
						if(mpz_tstbit(inputValue.get_mpz_t(), 0)==1)
							sticky=1;
						inputValue=inputValue/2;
					}
					else
						inputValue=inputValue*2;
			else {}
		else /* if we are counting ones */
		{
			int restOfBits = wIn - icount;
			if (icount>0)
			{
				mpz_class ones=1;
				for (int i=1;i<=icount;i++)
					ones = ones*2;
			
				ones=ones-1;
						
				for (int i=1;i<=restOfBits;i++)
					ones=ones*2;
				inputValue=inputValue-ones; // the input without the leading ones
			}

			if ((wIn<=wOut) || ((wIn>wOut) && (restOfBits<wOut) ))	//shift result in place	
				for (int i=1;i<=(wOut-restOfBits);i++)
					inputValue=inputValue*2;
			else
				for (int i=1;i<=restOfBits-wOut;i++){
					if(mpz_tstbit(inputValue.get_mpz_t(), 0)==1)
							sticky=1;
					inputValue=inputValue/2;
				}
		}
				
		so=inputValue;
		
		
		
		if (computeSticky)
		{
		//cout << "compute sticky !!!";
		ssticky = sticky;
		}	
		
		
		
		if (verbose)
		{
			cout<<"TestCase report:"<<endl;
			cout<<tab<<"Input value:"<<si<<endl;
			cout<<tab<<"Count value:"<<icount<<endl;
			cout << tab <<"Output value:"<<so<<endl;
		}		
			
	}
	
	
}



