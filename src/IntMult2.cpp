#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntMult2.hpp"

using namespace std;
extern vector<Operator*> oplist;

IntMult2:: IntMult2(Target* target, int wInX, int wInY) :
	Operator(target), wInX_(wInX), wInY_(wInY), wOut_(wInX + wInY){
 
	int depth, i, j;

	setOperatorType();
	setOperatorName();
	
	addInput ("X", wInX_);
	addInput ("Y", wInY_);
	addOutput("R", wOut_); /* wOut_ = wInX_ + wInY_ */

	partsX=4;
	partsY=4;
	size=17;

	for (j=1;j<=partsY;j++){
		ostringstream jj; 
		jj<<j;
		for (i=1;i<=partsX;i++){
			ostringstream ii;
			ii<<i;
		
			addDelaySignal("p"+jj.str()+ii.str(),2*size,i-2);
			if (i==1) addDelaySignal("p"+jj.str()+ii.str()+"Low",size,partsX-i);
			else{
				addDelaySignal("s"+jj.str()+ii.str(),2*size+(i==partsX?0:1),0);
				if (i<partsX) addDelaySignal("s"+jj.str()+ii.str()+"Low",size,partsX-i);
			}
		}
		addSignal("sum"+jj.str(),size*partsX);
	}
	
	setPipelineDepth(partsX);

}  

IntMult2::~IntMult2() {
}

void IntMult2::setOperatorName(){	
	/* Name Setup procedure
	 *  The name has the format: IntMult2_wInX__wInY_
	 *  wInX_ = width of the X input
	 *  wInY_ = width of the Y input
	 */  
	ostringstream name;
	name <<"IntMult22_"<<wInX_<<"_"<<wInY_;
	uniqueName_ = name.str(); 
	
}

void IntMult2::outputVHDL(std::ostream& o, std::string name) {
	int i;
	licence(o,"Bogdan Pasca(2008)");
	Operator::stdLibs(o);
	outputVHDLEntity(o);
	newArchitecture(o,name);
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	outputVHDLRegisters(o);
	o<<"-- multiplications and summations"<<endl;
	o<<"process(clk,X,Y)"<<endl;
	o<<"begin"<<endl;
	o<<tab<<"if clk'event and clk='1' then"<<endl;
	for (i=1;i<=partsX;i++){
			o<<tab<<tab<<"p1"<<i<<" <= X("<<i*size-1<<" downto "<<(i-1)*size<<") * Y("<<size-1<<" downto 0);"<<endl;
			if (i>1){
				o<<tab<<tab<<"s1"<<i<<" <= ("<<(i==partsX?"":"\"0\"&")<<"p1"<<i<<delaySignal("",i-2)<<")+"; 
				if (i>2) o<<"s1"<<i-1<<"("<<2*size<<" downto "<<size<<");"<<endl;
				else     o<<"p1"<<i-1<<"("<<2*size-1<<" downto "<<size<<");"<<endl;
			}
	}	
	o<<tab<<"end if;"<<endl;
	o<<"end process;"<<endl;	

	o<<"--assigmnents"<<endl;
	for(i=1;i<=partsX;i++){
		if (i==1) o<<"p1"<<i<<"Low <= p11("<<size-1<<"  downto 0);"<<endl;
		else if (i!=partsX) o<<"s1"<<i<<"Low <= s1"<<i<<"("<<size-1<<" downto 0);"<<endl;
	}
		
	o<<"R<=";
	for (i=partsX;i>=1;i--){
		if (i>1) o<<"s1"<<i<<(i<partsX?"Low":"")<<delaySignal("",partsX-i); else o<<"p1"<<i<<"Low"<<delaySignal("",partsX-i);
		if (i>1) o<<" & "; else o<<";"<<endl;
	}

	endArchitecture(o);
}



void IntMult2::fillTestCase(mpz_class a[])
{
	mpz_class &x = a[0];
	mpz_class &y = a[1];
	mpz_class &r = a[2];

	r = x * y;
}

