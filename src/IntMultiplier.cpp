/*
   An integer multiplier mess for FloPoCo

Authors:  Bogdan Pasca, being cleaned by F de Dinechin, Kinga Illyes and Bogdan Popa

This file is part of the FloPoCo project
developed by the Arenaire team at Ecole Normale Superieure de Lyon

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
2008-2010.
All rights reserved.
*/


/*
   Interface TODO
   support shared bitheap! In this case,
   - do not compress at the end
   - do not add the round bit
   - import the heap LSB
   - export the truncation error
   - ...
   */


/* 
   Who calls whom
   the constructor calls buildLogicOnly or buildTiling
   (maybe these should be unified some day)
   They call buildHeapTiling or buildHeapLogicOnly

*/



/*
   TODO tiling

   - define intermediate data struct (list of multiplier blocks)  
   multiplier block:
   - a bit saying if it should go into a DSP 
   - x and y size
   - x and y position
   - cycle ?
   - pointer to the previous (and the following?) tile if this tile belongs to a supertile

   - a function that explores and builds this structure
   recycle Bogdan's, with optim for large mults
   at least 4 versions: (full, truncated)x(altera, xilinx), please share as much code as possible

   - a function that generates VHDL (adding bits to the bit heap)
   */

/* VHDL variable names:
   X, Y: inputs
   XX,YY: after swap

*/






/* For two's complement arithmetic on n bits, the representable interval is [ -2^{n-1}, 2^{n-1}-1 ]
   so the product lives in the interval [-2^{n-1}*2^{n-1}-1,  2^n]
   The value 2^n can only be obtained as the product of the two minimal negative input values
   (the weird ones, which have no symmetry)
   Example on 3 bits: input interval [-4, 3], output interval [-12, 16] and 16 can only be obtained by -4*-4.
   So the output would be representable on 2n-1 bits in two's complement, if it werent for this weird*weird case.

   So for signed multipliers, we just keep the 2n bits, including one bit used for only for this weird*weird case.
   Current situation is: this must be managed from outside:
   An application that knows that it will not use the full range (e.g. signal processing, poly evaluation) can ignore the MSB bit, 
   but we produce it.



   A big TODO ?

   But for truncated signed multipliers, we should hackingly round down this output to 2^n-1 to avoid carrying around a useless bit.
   This would be a kind of saturated arithmetic I guess.

   Beware, the last bit of Baugh-Wooley tinkering does not need to be added in this case. See the TODO there.

   So the interface must provide a bit that selects this behaviour.

   A possible alternative to Baugh-Wooley that solves it (tried below but it doesn't work, zut alors)
   initially complement (xor) one negative input. This cost nothing, as this xor will be merged in the tables. Big fanout, though.
   then -x=xb+1 so -xy=xb.y+y    
   in case y pos, xy = - ((xb+1)y)  = -(xb.y +y)
   in case x and y neg, xy=(xb+1)(yb+1) = xb.yb +xb +yb +1
   It is enough not to add this lsb 1 to round down the result in this case.
   As this is relevant only to the truncated case, this lsb 1 will indeed be truncated.
   let sx and sy be the signs

   unified equation:

   px = sx xor rx  (on one bit less than x)
   py = sy xor ry

   xy = -1^(sx xor sy)( px.py + px.syb + py.sxb  )   
   (there should be a +sxsy but it is truncated. However, if we add the round bit it will do the same, so the round bit should be sx.sy)
   The final negation is done by complementing again.  


   Note that this only applies to truncated multipliers.

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
#include "IntMultiplier.hpp"
#include "IntAdder.hpp"
using namespace std;

namespace flopoco {

#define vhdl parentOp->vhdl
#define declare parentOp->declare
#define inPortMap parentOp->inPortMap
#define outPortMap parentOp->outPortMap
#define instance parentOp->instance
#define manageCriticalPath parentOp->manageCriticalPath
#define oplist parentOp->getOpListR()

	int IntMultiplier::neededGuardBits(int wX, int wY, int wOut)
	{
		int g;
		if(wX+wY==wOut)
			g=0;
		else {
			unsigned i=0;
			mpz_class ulperror=1;
			while(wX+wY - wOut  > intlog2(ulperror)) {
				// REPORT(DEBUG,"i = "<<i<<"  ulperror "<<ulperror<<"  log:"<< intlog2(ulperror) << "  wOut= "<<wOut<< "  wFull= "<< wX+wY);
				i++;
				ulperror += (i+1)*(mpz_class(1)<<i);
			}
			g=wX+wY-i-wOut;
			// REPORT(DEBUG, "ulp truncation error=" << ulperror << "    g=" << g);
		}
		return g;
	}


	void IntMultiplier::initialize() 
	{
		if(wOut<0 || wXdecl<0 || wYdecl<0) {
			THROWERROR("negative input/output size");
		}

		wFull = wXdecl+wYdecl;

		if(wOut > wFull){
			THROWERROR("wOut=" << wOut << " too large for " << (signedIO?"signed":"unsigned") << " " << wXdecl << "x" << wYdecl <<  " multiplier." );
		}

		if(wOut==0){ 
			wOut=wFull;
		}

		if(wOut<min(wXdecl, wYdecl))
			REPORT(0, "wOut<min(wX, wY): it would probably be more economical to truncate the inputs first, instead of building this multiplier.");

		wTruncated = wFull - wOut;

		g = neededGuardBits(wXdecl, wYdecl, wOut);
		REPORT(DEBUG, "    g=" << g);

		maxWeight = wOut+g;
		weightShift = wFull - maxWeight;  



		// Halve number of cases by making sure wY<=wX:
		// interchange x and y in case wY>wX

		if(wYdecl> wXdecl)
		{
			wX=wYdecl;	 
			wY=wXdecl;	 
			vhdl << tab << declare(addUID("XX"), wX, true) << " <= " << yname << " ;" << endl;	 
			vhdl << tab << declare(addUID("YY"), wY, true) << " <= " << xname << " ;" << endl;	 
		}
		else
		{
			wX=wXdecl;	 
			wY=wYdecl;	 
			vhdl << tab << declare(addUID("XX"), wX, true) << " <= " << xname << ";" << endl;	 
			vhdl << tab << declare(addUID("YY"), wY, true) << " <= " << yname << ";" << endl;	 
		}		

	}



	// The virtual constructor
	IntMultiplier::IntMultiplier (Operator* parentOp_, BitHeap* bitHeap_, Signal* x_, Signal* y_, int wX_, 
			int wY_, int wOut_, int lsbWeight_, bool negate_, bool signedIO_, float ratio_):
		Operator ( parentOp_->getTarget()), 
		wXdecl(wX_), wYdecl(wY_), wOut(wOut_), ratio(ratio_),  maxError(0.0), 
		parentOp(parentOp_), bitHeap(bitHeap_), lsbWeight(lsbWeight_),
		x(x_), y(y_), negate(negate_), signedIO(signedIO_) 
	{

			isOperator=false;

			multiplierUid=parentOp->getNewUId();
			srcFileName="IntMultiplier";
			useDSP = (ratio>=0) && getTarget()->hasHardMultipliers();

			ostringstream name;
			name <<"VirtualIntMultiplier";
			if(useDSP) 
				name << "UsingDSP_";
			else
				name << "LogicOnly_";
			name << wXdecl << "_" << wYdecl <<"_" << wOut << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
			setName ( name.str() );

			xname = x->getName();
			yname = y->getName();

			plotter = new Plotter(bitHeap);

			bitHeap->setPlotter(plotter);
			bitHeap->setSignedIO(signedIO);
			//plotter->setBitHeap(bitHeap);

			initialize();

			fillBitHeap();


			// leave the compression to the parent op
		}







	// The constructor for a stand-alone operator
	IntMultiplier::IntMultiplier (Target* target, int wX_, int wY_, int wOut_, bool signedIO_, float ratio_, map<string, double> inputDelays_):
		Operator ( target, inputDelays_ ), wXdecl(wX_), wYdecl(wY_), wOut(wOut_), ratio(ratio_), maxError(0.0), signedIO(signedIO_) 
	{

			isOperator=true;
			srcFileName="IntMultiplier";
			setCopyrightString ( "Florent de Dinechin, Kinga Illyes, Bogdan Popa, Bogdan Pasca, 2012" );

			// useDSP or not? 
			//useDSP = (ratio>0) && target->hasHardMultipliers();
			useDSP = (ratio>=0)&&target->hasHardMultipliers();


			{
				ostringstream name;
				name <<"IntMultiplier";
				if(useDSP) 
					name << "UsingDSP_";
				else
					name << "LogicOnly_";
				name << wXdecl << "_" << wYdecl <<"_" << wOut << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
				setName ( name.str() );
			}

			parentOp=this;
			multiplierUid=parentOp->getNewUId();
			xname="X";
			yname="Y";

			initialize();
			lsbWeight=g; // g was computed in initialize()

			// Set up the IO signals
			addInput ( xname  , wXdecl, true );
			addInput ( yname  , wYdecl, true );

			// TODO FIXME This 1 should be 2. It breaks TestBench but not TestBenchFile. Fix TestBench first! (check in addExpectedOutput or something)
			addOutput ( "R"  , wOut, 1 , true );

			// Set up the VHDL library style
			if(signedIO)
				useStdLogicSigned();
			else
				useStdLogicUnsigned();

			// The bit heap
			bitHeap = new BitHeap(this, wOut+g);

			plotter = new Plotter(bitHeap);

			bitHeap->setPlotter(plotter);
			bitHeap->setSignedIO(signedIO);




			// initialize the critical path
			setCriticalPath(getMaxInputDelays ( inputDelays_ ));

			fillBitHeap();

			bitHeap -> generateCompressorVHDL();			
			vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << range(wOut+g-1, g) << ";" << endl;
		}





	void  IntMultiplier::fillBitHeap()
	{

		///////////////////////////////////////
		//  architectures for corner cases   //
		///////////////////////////////////////

		// To manage stand-alone cases, we just build a bit-heap of max height one, so the compression will do nothing

		// The really small ones fit in two LUTs and that's as small as it gets  
		if(wX+wY <= target()->lutInputs()+2) {
			vhdl << tab << "-- Ne pouvant me fier à mon raisonnement, j'ai appris par coeur le résultat de toutes les multiplications possibles" << endl;

			SmallMultTable *t = new SmallMultTable( target(), wX, wY, wOut, negate, signedIO, signedIO);
			//useSoftRAM(t);
			//			oplist.push_back(t);
			t->addToGlobalOpList();

			vhdl << tab << declare(addUID("XY"), wX+wY) << " <= "<<addUID("YY")<<" & "<<addUID("XX")<<";"<<endl;

			inPortMap(t, "X", addUID("XY"));
			outPortMap(t, "Y", addUID("RR"));
			vhdl << instance(t, "multTable");

			plotter->addSmallMult(0,0,wX,wY);
			bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);

			for(int w=wOut-1; w>=0; w--)	{ // this is a weight in the multiplier output
				stringstream s;
				s<<addUID("RR")<<of(w);
				bitHeap->addBit(lsbWeight-g + w, s.str()); 
			}		
			return;
		}

		// Multiplication by 1-bit integer is simple
		if ((wY == 1)){
			if (signedIO){
				manageCriticalPath( target()->localWireDelay(wX) + target()->adderDelay(wX+1) );

				vhdl << tab << addUID("R") <<" <= (" << zg(wX+1)  << " - ("<<addUID("XX")<< of(wX-1) 
					<< " & "<<addUID("XX")<<")) when "<<addUID("YY")<<"(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	
			}
			else {
				manageCriticalPath( target()->localWireDelay(wX) + target()->lutDelay() );

				vhdl << tab << addUID("R")<<" <= (\"0\" & "<<addUID("XX")<<") when "<< addUID("YY")<<"(0)='1' else "<<zg(wX+1,0)<<";"<<endl;	

			}
			outDelayMap[addUID("R")] = getCriticalPath();
			return;
		}


		if ((wY == 2)) {
			// Multiplication by 2-bit integer is one addition, which is delegated to BitHeap compression anyway

			vhdl << tab << declare(addUID("R0"),wX+2) << " <= (";
			if (signedIO) 
				vhdl << addUID("XX") << of(wX-1) << " & "<< addUID("XX") << of(wX-1);  
			else  
				vhdl << "\"00\"";
			vhdl <<  " & "<<addUID("XX")<<") when "<<addUID("YY")<<"(0)='1' else "<<zg(wX+2,0)<<";"<<endl;	

			vhdl << tab << declare(addUID("R1i"),wX+2) << " <= ";
			if (signedIO) 
				vhdl << "("<<addUID("XX")<< of(wX-1) << "  &  " <<addUID("XX")<<" & \"0\")";
			else  
				vhdl << "(\"0\" & "<< addUID("XX") <<" & \"0\")";
			vhdl << " when "<<addUID("YY")<<"(1)='1' else " << zg(wX+2,0) << ";"<<endl;	

			vhdl << tab << declare(addUID("R1"),wX+2) << " <= ";
			if (signedIO) 
				vhdl << "not "<<addUID("R1i")<<";" <<endl;
			else  
				vhdl << addUID("R1i")<<";"<<endl;

			for (int w=0; w<wOut+g; w++){
				stringstream s0,s1;
				s0<<addUID("R0")<<of(w+(wX+2-wOut-g));
				bitHeap->addBit(lsbWeight-g + w, s0.str());
				s1<<addUID("R1")<<of(w+(wX+2-wOut-g));
				bitHeap->addBit(lsbWeight-g + w, s1.str());
			}
			// Rounding bit (or carry in bit for signed inputs)
			if(g || signedIO)
				bitHeap->addConstantOneBit(lsbWeight-g);
			// and that's it
			return;

		} 




		// Now getting more and more generic
		if(useDSP) {
			//test if the multiplication fits into one DSP
			//int wxDSP, wyDSP;
			target()->getDSPWidths(wxDSP, wyDSP, signedIO);
			bool testForward, testReverse, testFit;
			testForward     = (wX<=wxDSP)&&(wY<=wyDSP);
			testReverse = (wY<=wxDSP)&&(wX<=wyDSP);
			testFit = testForward || testReverse;


			REPORT(DETAILED,"useDSP");
			if (testFit){
				REPORT(DETAILED,"testfit");

				if( checkThreshold(0,0,wX, wY,wxDSP,wyDSP))
				{	REPORT(INFO,"worthUsingDSP");
					manageCriticalPath(target()->DSPMultiplierDelay());
					/*if (signedIO)
					  {
					  vhdl << tab << declare(addUID("rfull"), wFull+1) << " <= "<<addUID("XX")<<"  *  "
					  << addUID("YY")<<"; -- that's one bit more than needed"<<endl; 
					  vhdl << tab << addUID("R")<<" <= "<< addUID("rfull") <<range(wFull-1, wFull-wOut)<<";"<<endl;	
					  outDelayMap[addUID("R")] = getCriticalPath();
					  }
					  else //sign extension is necessary for using use ieee.std_logic_signed.all; 
					  {	// for correct inference of Xilinx DSP functions

					//vhdl << tab << declare(addUID("rfull"), wX + wY + 2) << " <= (\"0\" & "<<addUID("XX")<<
					//") * (\"0\" &"<<addUID("YY")<<");"<<endl;
					*/

					int topx=wX-wxDSP;
					int topy=wY-wyDSP;

					REPORT(DETAILED,"wxDSSSSPPP=="<<wxDSP);

					stringstream inx,iny;
					inx<<addUID("XX");
					iny<<addUID("YY");

					MultiplierBlock* m = new MultiplierBlock(wxDSP,wyDSP,topx, topy,
							inx.str(),iny.str(),weightShift + lsbWeight-g);
					m->setNext(NULL);		
					m->setPrevious(NULL);			
					localSplitVector.push_back(m);
					bitHeap->addMultiplierBlock(m);
					bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);
					//	}	

				}


				else {
					// For this target and this size, better do a logic-only implementation
					REPORT(DETAILED,"before buildlogic only");
					buildLogicOnly();
				}
			}
			else {

				buildTiling();
				//buildLogicOnly();
				//buildHeapLogicOnly(0,0,50,25,2);

			}
		} 



		else {// This target has no DSP, going for a logic-only implementation

			buildLogicOnly();
		}

		bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);
	}






	/**************************************************************************/
	void IntMultiplier::buildLogicOnly() 
	{
		buildHeapLogicOnly(0,0,wX,wY);
		//adding the round bit
		if(g>0) {
			int weight=wFull-wOut-1- weightShift;
			if(negate)
				bitHeap->subConstantOneBit(lsbWeight-g + weight);
			else
				bitHeap->addConstantOneBit(lsbWeight-g + weight);
		}
	}



	/**************************************************************************/
	void IntMultiplier::buildTiling() 
	{


		if((!signedIO)&&((wX==41)&&(wY==41)&&(wFull-wOut-g==0)))
			buildFancyTiling();
		else
			buildHeapTiling();
		//adding the round bit
		//bitHeap->getPlotter()->plotMultiplierConfiguration(multiplierUid, localSplitVector, wX, wY, wOut, g);
		if(g>0) {
			int weight=wFull-wOut-1- weightShift;
			if(negate)
				bitHeap->subConstantOneBit(lsbWeight-g + weight);
			else
				bitHeap->addConstantOneBit(lsbWeight-g + weight);
		}


	}


	//the fancy tiling is used only for a hardwired case 41 41 82 
	/***********************************************************************/
	void IntMultiplier::buildFancyTiling()
	{
		//THROWERROR("fancy tiling not implemented yet");

		stringstream inx,iny;
		inx<<addUID("XX");
		iny<<addUID("YY");

		//REPORT(INFO,"chr DSP at "<<topx<<" "<<topy<<" width= "<<widthX<<" height= "<<widthY);		
		int widthX, widthY,topx,topy;

		//topright dsp;
		widthX=wxDSP;
		widthY=wyDSP;
		topx=0;
		topy=0;
		MultiplierBlock* m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(),weightShift + lsbWeight-g);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);


		//topleft dsp
		widthX=wyDSP;
		widthY=wxDSP;
		topx=wxDSP;
		topy=0;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(),weightShift + lsbWeight-g);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//bottomleft dsp
		widthX=wxDSP;
		widthY=wyDSP;
		topx=wyDSP;
		topy=wxDSP;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(),weightShift + lsbWeight-g);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//bottomright dsp
		widthX=wyDSP;
		widthY=wxDSP;
		topx=0;
		topy=wyDSP;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(),weightShift + lsbWeight-g);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//logic

		buildHeapLogicOnly(wyDSP,wyDSP,wxDSP,wxDSP,parentOp->getNewUId());

	}


	/**************************************************************************/
	void IntMultiplier::buildHeapLogicOnly(int topX, int topY, int botX, int botY,int blockUid)
	{

		REPORT(INFO,"topX"<<topX<<" topY"<<topY<<" botX"<<botX<<" botY"<<botY);
		Target *target=getTarget();
		if(blockUid==-1)
			blockUid++;    /// ???????????????

		vhdl << tab << "-- code generated by IntMultiplier::buildHeapLogicOnly()"<< endl;


		int dx, dy;				// Here we need to split in small sub-multiplications
		int li=target->lutInputs();

		dx = li>>1;
		dy = li - dx; 
		REPORT(DEBUG, "dx="<< dx << "  dy=" << dy );

		int wXX=wX;
		int wYY=wY;

		int wX=botX-topX;
		int wY=botY-topY;
		int chunksX =  int(ceil( ( double(wX) / (double) dx) ));
		int chunksY =  int(ceil( ( 1+ double(wY-dy) / (double) dy) ));
		int sizeXPadded=dx*chunksX; 
		int sizeYPadded=dy*chunksY;
		int padX=sizeXPadded-wX;
		int padY=sizeYPadded-wY;

		REPORT(DEBUG, "X split in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
		REPORT(DEBUG, " sizeXPadded="<< sizeXPadded << "  sizeYPadded="<< sizeYPadded);
		if (chunksX + chunksY > 2) { //we do more than 1 subproduct // FIXME where is the else?

			// Padding X to the right
			vhdl << tab << declare(addUID("Xp", blockUid), sizeXPadded) << " <= ";
			vhdl << addUID("XX") << range(botX-1,topX) << " & "<<zg(padX)<<";"<<endl;
			REPORT(DETAILED,addUID("XX") << range(botX-1,topX) << " & "<<zg(padX)<<";");

			// Padding Y to the right
			vhdl << tab << declare(addUID("Yp",blockUid), sizeYPadded)<<" <= ";
			vhdl << addUID("YY") << range(botY-1,topY) << " & "<<zg(padY)<<";"<<endl;

			REPORT(DETAILED,addUID("YY") << range(botY-1,topY) << " & "<<zg(padY)<<";");
			//SPLITTINGS
			for (int k=0; k<chunksX ; k++)
			{
				vhdl << tab << declare(join(addUID("x",blockUid),"_",k),dx) << " <= "<< addUID("Xp",blockUid) << range((k+1)*dx-1,k*dx)<<";"<<endl;
				REPORT(DETAILED,join(addUID("x",blockUid),"_",k)<<" <= "<< addUID("Xp",blockUid) << range((k+1)*dx-1,k*dx)<<";");
			}	
			for (int k=0; k<chunksY ; k++)
			{
				vhdl << tab << declare(join(addUID("y",blockUid),"_",k),dy) << " <= " << addUID("Yp",blockUid) << range((k+1)*dy-1, k*dy)<<";"<<endl;
				REPORT(DETAILED,join(addUID("y",blockUid),"_",k)<<" <= "<< addUID("Yp",blockUid) << range((k+1)*dy-1,k*dy)<<";");
			}	


			REPORT(DEBUG, "maxWeight=" << maxWeight <<  "    weightShift=" << weightShift);
			SmallMultTable *tUU, *tSU, *tUS, *tSS;


			tUU = new SmallMultTable( target, dx, dy, dx+dy, negate, false, false);
			//useSoftRAM(tUU);
			//oplist.push_back(tUU);
			tUU->addToGlobalOpList();

			if(signedIO) { // need for 4 different tables

				tSU = new SmallMultTable( target, dx, dy, dx+dy, negate, true, false );
				//useSoftRAM(tSU);
				// oplist.push_back(tSU);
				tSU->addToGlobalOpList();

				tUS = new SmallMultTable( target, dx, dy, dx+dy, negate, false, true );
				//useSoftRAM(tUS);
				// oplist.push_back(tUS);
				tUS->addToGlobalOpList();


				tSS = new SmallMultTable( target, dx, dy, dx+dy, negate, true, true );
				//useSoftRAM(tSS);
				//oplist.push_back(tSS);
				tSS->addToGlobalOpList();

			}

			setCycle(0); // TODO FIXME for the virtual multiplier case where inputs can arrive later
			setCriticalPath(initialCP);
			// SmallMultTable is built to cost this:
			manageCriticalPath( getTarget()->localWireDelay(chunksX) + getTarget()->lutDelay() ) ;  
			for (int iy=0; iy<chunksY; iy++){

				vhdl << tab << "-- Partial product row number " << iy << endl;

				for (int ix=0; ix<chunksX; ix++) { 


					SmallMultTable *t;
					if (!signedIO) 
					{
						t=tUU;
					}

					else 
					{ // 4  cases 
						if( ((ix==chunksX-1)&&(botX==wXX)) && ((iy==chunksY-1)&&(botY==wYY) ))
							t=tSS;
						else if ((ix==chunksX-1)&&(botX==wXX)) 
							t=tSU;
						else if ((iy==chunksY-1)&&(botY==wYY))
							t=tUS;
						else
							t=tUU; 
					}



					if(dx*(ix+1)+dy*(iy+1)+topX+topY-padX-padY>wFull-wOut-g)
					{



						plotter->addSmallMult(dx*(ix)+topX-padX, dy*(iy)+topY-padY,dx,dy);

						REPORT(DETAILED,XY(ix,iy,blockUid)<<" <= " << addUID("y",blockUid) <<"_"<< iy << " & " << addUID("x",blockUid) <<"_"<< ix << ";");
						vhdl << tab << declare(XY(ix,iy,blockUid), dx+dy) 
							<< " <= " << addUID("y",blockUid) <<"_"<< iy << " & " << addUID("x",blockUid) <<"_"<< ix << ";"<<endl;

						inPortMap(t, "X", XY(ix,iy,blockUid));
						outPortMap(t, "Y", PP(ix,iy,blockUid));
						vhdl << instance(t, PPTbl(ix,iy,blockUid));

						vhdl << tab << "-- Adding the relevant bits to the heap of bits" << endl;

						// Two's complement management
						// There are really 2 cases:
						// If the result is known positive, ie if tUU and !negate, nothing to do
						// If the result is in two's complement  
						//    sign extend by adding ones on weights >= the MSB of the table, so its sign is propagated.
						//    Also need to complement the sign bit
						// Note that even when negate and tUU, the result is known negative, but it may be zero, so its actual sign is not known statically


						bool resultSigned =false;  
						if(negate ||(t==tSS) || (t==tUS) || (t==tSU)) 
							resultSigned = true ;

						int maxK=dx+dy; 
						int minK=0;
						//if(ix == chunksX-1)
						if ((ix == 0))
							minK+=padX;
						if((iy == 0))
							//maxK-=padY;
							minK+=padY;
						REPORT(INFO,"The bits will be added from: mink"<<minK<<"	 maxk"<<maxK);
						REPORT(DEBUG,  "ix=" << ix << " iy=" << iy << "  maxK=" << maxK  << "  negate=" << negate <<  "  resultSigned="  << resultSigned );
						for (int k=minK; k<maxK; k++) {
							ostringstream s;
							s << PP(ix,iy,blockUid) << of(k); // right hand side

							int weight = ix*dx+iy*dy+k-weightShift+topX+topY-padX-padY;

							if(weight>=0) {// otherwise these bits deserve to be truncated
								if(resultSigned && (k==maxK-1)) { 
									ostringstream nots;
									nots << "not " << s.str(); 
									bitHeap->addBit(lsbWeight-g + weight, nots.str());
									REPORT(INFO, "adding neg bit " << nots.str() << " at weight " << weight); 
									REPORT(INFO,  "  adding constant ones from weight "<< weight << " to "<< bitHeap->getMaxWeight());
									for (unsigned w=weight; w<bitHeap->getMaxWeight(); w++)
										bitHeap->addConstantOneBit(lsbWeight-g + w);
								}
								else { // just add the bit
									REPORT(INFO, "adding bit " << s.str() << " at weight " << weight); 
									bitHeap->addBit(lsbWeight-g + weight, s.str());
								}
							}
						}

						vhdl << endl;

					}
				}
			}				



		}

	}



			void IntMultiplier::addExtraDSPs(int topX, int topY, int botX, int botY,int wxDSP,int wyDSP)
		{
		
				
				int widthX=wxDSP;
				int widthY=wyDSP;
				int botx=botX;
			
				if(signedIO)
				{
					if(botx!=wX)
						widthX--;
					if(botY!=wY)	
						widthY--;
				}
			
		
				int height=botY-topY;
				int width=botX-topX;
				int dspArea=widthX*widthY;
				bool was=false; // tells if the while loop was executed or not
				int topx=topX;
				int topy=topY;
				int dsp=0;//number of used DSPs
				int usedWidth=0;//the length of the used width on x
			
			
			
			
				//if the width is larger then a dsp width, than we have to checkTreshHold the good coordinates for the dsp
				if (width>widthX)
					topx=botx-widthX;
		
				while (width>widthX)
				{	
					REPORT(DETAILED,"width greater than DSPwidth");
					//we need to split the block
					was=true;
					float blockArea=widthX*height; //the area of the block that will be analyzed
				
					if(checkThreshold(topx,topy,botx,botY,widthX,widthY))
					{  
				
						if(height<widthY)
							topy=topY-(widthY-height);
						stringstream inx,iny;
						inx<<addUID("XX");
						iny<<addUID("YY");
								
						//REPORT(INFO,"chr DSP at "<<topx<<" "<<topy<<" width= "<<widthX<<" height= "<<widthY);							
						MultiplierBlock* m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(),weightShift + lsbWeight-g);
						m->setNext(NULL);		
						m->setPrevious(NULL);			
						localSplitVector.push_back(m);
						bitHeap->addMultiplierBlock(m);
					}
					else
					{
					
						if((topx<botX-dsp*widthX-1))
							buildHeapLogicOnly(topx, topY,(botX-usedWidth),botY,parentOp->getNewUId());	
						//REPORT(INFO,"Logic computed xt="<<topx<<" yt="<<topY<<"  xb="<<(botX-dsp*widthX)<<" yb="<<botY); 	
					}
				
					//updating the coordinates for the next block
					dsp++;
					botx=topx-1;
					topx=topx-widthX;
					width=width-widthX;
					usedWidth+=widthX;
					widthX=wxDSP;
					widthY=wyDSP;
						
					if(signedIO)
					{
						if(botx!=wX)
							widthX--;
						if(botY!=wY)	
							widthY--;
					}	
					
					if(width<=widthX)
						topx=topX;
					
				}
				//now the remaining block is smaller (for sure!) than the DSP width
				//we do the same area / ratio checking again	
				if(checkThreshold(topx,topy,botx,botY,widthX,widthY))
				{
			
					//computing the coordinates
					if(was)
						topx=topx-(widthX-(botx-topx+1));
					else 
						topx=topX-(widthX-(botX-topX+1))-1; 
							
					if(height<widthY)
						topy=topY-(widthY-height);
							
					stringstream inx,iny;
								
					MultiplierBlock* m = new MultiplierBlock(widthX,widthY,topx,topy,addUID("XX"),addUID("YY"),weightShift + lsbWeight-g);
					m->setNext(NULL);		
					m->setPrevious(NULL);			
					localSplitVector.push_back(m);
					bitHeap->addMultiplierBlock(m);
				
				}
				else
				{
					//build logic only is at least 1 bit to compute
					if((topx<botX-dsp*widthX))
						buildHeapLogicOnly(topx,topY,botX-usedWidth,botY,parentOp->getNewUId());
							
				}
				
		
		}	
	
	



		bool IntMultiplier::checkThreshold(int topX, int topY, int botX, int botY,int wxDSP,int wyDSP)
		{
			int widthX=botX-topX;
			int widthY=botY-topY;
			double blockArea;
			double triangle=0.0;
			double dspArea=wxDSP*wyDSP;
		
			//drop the unused area
			if(botX+topY<wFull-wOut-g)
			{
				int x=botX;
				int y=topY;
				while(x+y<wFull-wOut-g)
					y++;
			
				topY=y;
			}
		
			//***the triangle is the area which will be lost from the block area***//
			//**computing the triangle's edge (45degree => area=(edge^2)/2) ***//		
			int x=topX;
			int y=topY;
			while(x+y<wFull-wOut-g)
			x++;
			//**************************************************//
			
			
			//computing the triangle's area
			if(topX+topY<=wFull-wOut-g)
				triangle=((x-topX)*(x-topX))/2;
			else
				triangle=0.0;
			
			
			//the final area which is used
			blockArea=widthX*widthY-triangle;		
		
			//checking according to ratio/area
			if(blockArea>=(1.0-ratio)*dspArea)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		
					
		void IntMultiplier::buildHeapTiling()
		{
			
			int botx=wX;
			int boty=wY;
			int topx=botx-wxDSP;
			int topy=boty-wyDSP;
			int widthX=wxDSP;
			int widthY=wyDSP;
			int ok;
			while(boty>0)
			{	ok=0;
				botx=wX;
				while((botx>0)&&(ok==0))
				{	
					widthX=wxDSP;
					widthY=wyDSP;
					
					if(signedIO)
					{
						if(boty!=wY)
							widthY--;
						if(botx!=wX)
							widthX--;	
					}
				
					
				
					topx=botx-widthX;
					topy=boty-widthY;
				
					if(topx+topy>=wFull-wOut-g)
					{
						MultiplierBlock* m = new MultiplierBlock(widthX,widthY,topx,topy,addUID("XX"),addUID("YY"),weightShift + lsbWeight-g);
						m->setNext(NULL);		
						m->setPrevious(NULL);			
						localSplitVector.push_back(m);
						bitHeap->addMultiplierBlock(m);
						botx=botx-widthX;
						ok=0;
						//REPORT(INFO,"DSP at "<<topx<<" "<<topy<<" width= "<<widthX<<" height= "<<widthY);
					}
					
					else
					{
						ok=1;
						
					}
					
					
				
				}
				
				//updating the dsp's width and height
				widthX=wxDSP;
				widthY=wyDSP;
					
				if(signedIO)
				{
					if(boty!=wY)
						widthY--;
					if(botx!=wX)
						widthX--;	
				}
				
				
				//determination of the x coordinate
				if(topy<0)
					topy=0;
				int y=boty;
				int x=wX;
				while((x+y>wFull-wOut-g) && (x>0))
					x--;
				
				
				//call the function only if at least 1 bit remaining
				if((botx>0))
					addExtraDSPs(x,topy, botx, boty, wxDSP, wyDSP); 
				
				boty=boty-widthY;
			}
		}



	IntMultiplier::~IntMultiplier() {
	}


	//signal name construction

	string IntMultiplier::addUID(string name, int blockUID)
	{
		ostringstream s;
		s << name << "_m" << multiplierUid;
		if (blockUID!=-1) 
			s << "b"<< blockUID;
		return s.str() ;
	};

	string IntMultiplier::PP(int i, int j, int uid ) {
		std::ostringstream p;		
		if(uid==-1) 
			p << "PP" <<  "_X" << i << "Y" << j;
		else
			p << "PP" <<uid<<"X" << i << "Y" << j;
		return  addUID(p.str());
	};

	string IntMultiplier::PPTbl(int i, int j, int uid ) {
		std::ostringstream p;		
		if(uid==-1) 
			p << addUID("PP") <<  "_X" << i << "Y" << j << "_Tbl";
		else
			p << addUID("PP") <<"_"<<uid<<"X" << i << "Y" << j << "_Tbl";
		return p.str();
	};


	string IntMultiplier::XY( int i, int j,int uid) {
		std::ostringstream p;		
		if(uid==-1) 
			p  << "Y" << j<< "X" << i;
		else
			p  << "Y" << j<< "X" << i<<"_"<<uid;
		return  addUID(p.str());	
	};






	IntMultiplier::SmallMultTable::SmallMultTable(Target* target, int dx, int dy, int wO, bool negate, bool  signedX, bool  signedY ) : 
		Table(target, dx+dy, wO, 0, -1, true), // logic table
		dx(dx), dy(dy), negate(negate), signedX(signedX), signedY(signedY) {
			ostringstream name; 
			srcFileName="LogicIntMultiplier::SmallMultTable";
			// No getUid() in the name: this kind of table should be added to the globalOpList 
			name <<"SmallMultTable"<< (negate?"M":"P") << dy << "x" << dx << "r" << wO << (signedX?"Xs":"Xu") << (signedY?"Ys":"Yu");
			setName(name.str());				
		};


	mpz_class IntMultiplier::SmallMultTable::function(int yx){
		mpz_class r;
		int y = yx>>dx;
		int x = yx -(y<<dx);
		int wF=dx+dy;

		if(signedX){
			if ( x >= (1 << (dx-1)))
				x -= (1 << dx);
		}
		if(signedY){
			if ( y >= (1 << (dy-1)))
				y -= (1 << dy);
		}
		// if(negate && signedX && signedY) cerr << "  y=" << y << "  x=" << x;
		r = x * y;
		// if(negate && signedX && signedY) cerr << "  r=" << r;
		if(negate)
			r=-r;
		// if(negate && signedX && signedY) cerr << "  -r=" << r;
		if ( r < 0)
			r += mpz_class(1) << wF; 
		// if(negate && signedX && signedY) cerr << "  r2C=" << r;

		if(wOut<wF){ // wOut is that of Table
			// round to nearest, but not to nearest even
			int tr=wF-wOut; // number of truncated bits 
			// adding the round bit at half-ulp position
			r += (mpz_class(1) << (tr-1));
			r = r >> tr;
		}

		// if(negate && signedX && signedY) cerr << "  rfinal=" << r << endl;

		return r;

	};



	void IntMultiplier::emulate ( TestCase* tc ) {
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svR;

		if (! signedIO){
			svR = svX * svY;
		}

		else{ // Manage signed digits
			mpz_class big1 = (mpz_class(1) << (wXdecl));
			mpz_class big1P = (mpz_class(1) << (wXdecl-1));
			mpz_class big2 = (mpz_class(1) << (wYdecl));
			mpz_class big2P = (mpz_class(1) << (wYdecl-1));

			if ( svX >= big1P)
				svX -= big1;

			if ( svY >= big2P)
				svY -= big2;

			svR = svX * svY;
		}
		if(negate)
			svR = -svR;

		// manage two's complement at output
		if ( svR < 0){
			svR += (mpz_class(1) << wFull); 
		}
		if(wTruncated==0) 
			tc->addExpectedOutput("R", svR);
		else {
			// there is truncation, so this mult should be faithful
			svR = svR >> wTruncated;
			tc->addExpectedOutput("R", svR);
			svR++;
			svR &= (mpz_class(1) << (wOut)) -1;
			tc->addExpectedOutput("R", svR);
		}
	}



	void IntMultiplier::buildStandardTestCases(TestCaseList* tcl)
	{
		TestCase *tc;

		mpz_class x, y;

		// 1*1
		x = mpz_class(1); 
		y = mpz_class(1); 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1 * -1
		x = (mpz_class(1) << wXdecl) -1; 
		y = (mpz_class(1) << wYdecl) -1; 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// The product of the two max negative values overflows the signed multiplier
		x = mpz_class(1) << (wXdecl -1); 
		y = mpz_class(1) << (wYdecl -1); 
		tc = new TestCase(this); 
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}




}
