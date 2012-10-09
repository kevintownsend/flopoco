
/*
  An integer multiplier mess for FloPoCo

  TODO in the virtual multiplier case, manage properly the case when the initial instant is not  0

  Authors:  Bogdan Pasca, but we (F de Dinechin, Kinga Illyes and Bogdan Popa) spent two months getting rid of the last bits of his code

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2012.
  All rights reserved.
*/


/* 
   Who calls whom
   the constructor calls buildLogicOnly or buildTiling
   (maybe these should be unified some day)
   They call buildXilinxTiling or buildAlteraTiling or buildHeapLogicOnly

*/



/* VHDL variable names:
   X, Y: inputs
   XX,YY: after swap

*/






/* For two's complement arithmetic on n bits, the representable interval is [ -2^{n-1}, 2^{n-1}-1 ]
   so the product lives in the interval [-2^{n-1}*2^{n-1}-1,  2^n]
   The value 2^n can only be obtained as the product of the two minimal negative input values
   (the weird ones, which have no opposite)
   Example on 3 bits: input interval [-4, 3], output interval [-12, 16] and 16 can only be obtained by -4*-4.
   So the output would be representable on 2n-1 bits in two's complement, if it werent for this weird*weird case.

   So even for signed multipliers, we just keep the 2n bits, including one bit used for only for this weird*weird case.
   Current situation is: if you don't like this you manage it from outside:
   An application that knows that it will not use the full range (e.g. signal processing, poly evaluation) can ignore the MSB bit, 
   but we produce it.



   A big TODO ?

   But for truncated signed multipliers, maybe we could hackingly round down this output to 2^n-1 to avoid carrying around a useless bit.
   This would be a kind of saturated arithmetic I guess.


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
#include "Targets/StratixII.hpp"
#include "Targets/StratixIII.hpp"
#include "Targets/StratixIV.hpp"
#include "Plotter.hpp"

using namespace std;

namespace flopoco {

#define vhdl parentOp->vhdl
#define declare parentOp->declare
#define inPortMap parentOp->inPortMap
#define outPortMap parentOp->outPortMap
#define instance parentOp->instance
#define manageCriticalPath parentOp->manageCriticalPath
#define getCriticalPath parentOp->getCriticalPath
#define setCycle parentOp->setCycle
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

		weightShift = wFull - (wOut+g); 
		REPORT(DEBUG,   "weightShift=" << weightShift);


		// Halve number of cases by making sure wY<=wX:
		// interchange x and y in case wY>wX
		// After which we negate y (the smaller) by 1/ complementing it and 2/  adding it back to the bit heap

		string newxname, newyname;
		if(wYdecl> wXdecl) {
			wX=wYdecl;	 
			wY=wXdecl;	 
			newxname=yname;
			newyname=xname;
		}
		else {
			wX=wXdecl;	 
			wY=wYdecl;
			newxname=xname;
			newyname=yname;
		}

		// The larger of the two 
		vhdl << tab << declare(addUID("XX"), wX, true) << " <= " << newxname << " ;" << endl;	 

		// possibly negate the smaller
		if(!negate)
			vhdl << tab << declare(addUID("YY"), wY, true) << " <= " << newyname << " ;" << endl;	 
		else {	
			vhdl << tab << "-- we compute -xy as x(not(y)+1)" << endl;
			vhdl << tab << declare(addUID("YY"), wY, true) << " <= not " << newyname << " ;" << endl;	 
		}
	}





	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// The virtual constructor 
	IntMultiplier::IntMultiplier (Operator* parentOp_, BitHeap* bitHeap_, Signal* x_, Signal* y_, int wX_, 
	                              int wY_, int wOut_, int lsbWeight_, bool negate_, bool signedIO_, float ratio_):
		Operator (  parentOp_->getTarget()), 
		wxDSP(0), wyDSP(0), wXdecl(wX_), wYdecl(wY_), wX(0), wY(0), wOut(wOut_), ratio(ratio_),  maxError(0.0), 
		parentOp(parentOp_), bitHeap(bitHeap_), lsbWeight(lsbWeight_),
		x(x_), y(y_), negate(negate_), signedIO(signedIO_) 
	{

		multiplierUid=parentOp->getNewUId();
		srcFileName="IntMultiplier";
		useDSP = (ratio>=0) &&  parentOp->getTarget()->hasHardMultipliers();

		ostringstream name;
		name <<"VirtualIntMultiplier";
		if(useDSP) 
			name << "UsingDSP_";
		else
			name << "LogicOnly_";
		name << wXdecl << "_" << wYdecl <<"_" << wOut << "_" << (signedIO?"signed":"unsigned") << "_uid"<<Operator::getNewUId();
		setName ( name.str() );

		REPORT(DEBUG, "Building " << name.str() )

			xname = x->getName();
		yname = y->getName();

		// // we create a separate Plotter for this mult. Is it useful?
		// plotter = new Plotter(bitHeap);
		//plotter->setBitHeap(bitHeap);

		bitHeap->setSignedIO(signedIO);
		

		initialize();

		fillBitHeap();


		// leave the compression to the parent op
	}







	// The constructor for a stand-alone operator
	IntMultiplier::IntMultiplier (Target* target, int wX_, int wY_, int wOut_, bool signedIO_, float ratio_, map<string, double> inputDelays_, bool enableSuperTiles_):
		Operator ( target, inputDelays_ ), wxDSP(0), wyDSP(0), wXdecl(wX_), wYdecl(wY_), wX(0), wY(0), wOut(wOut_), wFull(0), ratio(ratio_),  maxError(0.0), negate(false), signedIO(signedIO_),enableSuperTiles(enableSuperTiles_) 
	{
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
			REPORT(DEBUG, "Building " << name.str() )
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
		bitHeap = new BitHeap(this, wOut+g, enableSuperTiles);

		// TODO CHECK ??? A bit heap is sign-agnostic.
		bitHeap->setSignedIO(signedIO);



		// initialize the critical path
		setCriticalPath(getMaxInputDelays ( inputDelays_ ));

		fillBitHeap();

		bitHeap -> generateCompressorVHDL();			
		vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << range(wOut+g-1, g) << ";" << endl;
	}






	void  IntMultiplier::fillBitHeap()
	{

		Plotter* plotter= bitHeap->getPlotter();
		///////////////////////////////////////
		//  architectures for corner cases   //
		///////////////////////////////////////

		// TODO Support negate in all the corner cases
		// To manage stand-alone cases, we just build a bit-heap of max height one, so the compression will do nothing
		// The really small ones fit in two LUTs and that's as small as it gets  
		if(wX+wY <=  parentOp->getTarget()->lutInputs()+2) {
			vhdl << tab << "-- Ne pouvant me fier à mon raisonnement, j'ai appris par coeur le résultat de toutes les multiplications possibles" << endl;

			SmallMultTable *t = new SmallMultTable(  parentOp->getTarget(), wX, wY, wOut, negate, signedIO, signedIO);
			//useSoftRAM(t);
			//			oplist.push_back(t);
			t->addToGlobalOpList();

			vhdl << tab << declare(addUID("XY"), wX+wY) << " <= "<<addUID("YY")<<" & "<<addUID("XX")<<";"<<endl;

			inPortMap(t, "X", addUID("XY"));
			outPortMap(t, "Y", addUID("RR"));
			vhdl << instance(t, "multTable");

			plotter->addSmallMult(0,0,wX,wY);
			plotter->plotMultiplierConfiguration(getName(), localSplitVector, wX, wY, wOut, g);

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
				manageCriticalPath(  parentOp->getTarget()->localWireDelay(wX) +  parentOp->getTarget()->adderDelay(wX+1) );

				vhdl << tab << addUID("R") <<" <= (" << zg(wX+1)  << " - ("<<addUID("XX")<< of(wX-1) 
				     << " & "<<addUID("XX")<<")) when "<<addUID("YY")<<"(0)='1' else "<< zg(wX+1,0)<<";"<<endl;	
			}
			else {
				manageCriticalPath(  parentOp->getTarget()->localWireDelay(wX) +  parentOp->getTarget()->lutDelay() );

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

		// TODO HERE

#if 0
		// Finish the negation of the smaller input by adding X (since -yx=not(y)x +x)
		setCycle(0); // TODO FIXME for the virtual multiplier case where inputs can arrive later
		setCriticalPath(initialCP);

		// TODO FIXME need to sign extend
		if(negate) {
			for(int i=0; i<wX; i++) {
				int w = lsbWeight + i-weightShift;
				if(w>=0) {
					ostringstream rhs;
					rhs << addUID("XX") << of(i);
					bitHeap->addBit(w, rhs.str());
				}
			}
		}
#endif

		if(useDSP) 
			{
				REPORT(DETAILED,"useDSP");
				 parentOp->getTarget()->getDSPWidths(wxDSP, wyDSP, signedIO);
				buildTiling();
			}

		else {// This target has no DSP, going for a logic-only implementation	
			buildLogicOnly();
		}


		// TODO weight=-1 here
		//add the round bit
		if(g>0) {
			int weight = lsbWeight-1;
			if(negate)
				bitHeap->subConstantOneBit(weight);
			else
				bitHeap->addConstantOneBit(weight);
		}

		
		plotter->plotMultiplierConfiguration(getName(), localSplitVector, wX, wY, wOut, g);
	}
	



	/**************************************************************************/
	void IntMultiplier::buildLogicOnly() 
	{
		buildHeapLogicOnly(0,0,wX,wY);
	}



	/**************************************************************************/
	void IntMultiplier::buildTiling() 
	{
		int* multiplierWidth;
		int size;	
		if( parentOp->getTarget()->getVendor()=="Altera")
			{	
				if ( parentOp->getTarget()->getID()=="StratixII")
					{
						StratixII* t = (StratixII*) parentOp->getTarget();
						multiplierWidth = t->getDSPMultiplierWidths();
						size = t->getNrDSPMultiplier();	

					}
				else
					if( parentOp->getTarget()->getID()=="StratixIII")
						{
							StratixIII* t = (StratixIII*) parentOp->getTarget();
							multiplierWidth = t->getDSPMultiplierWidths();
							size = t->getNrDSPMultiplier();	

						}
					else
						if( parentOp->getTarget()->getID()=="StratixIV")
							{
								StratixIV* t = (StratixIV*) parentOp->getTarget();
								multiplierWidth = t->getDSPMultiplierWidths();
								size = t->getNrDSPMultiplier();	

							}
				//add Altera boards here
						else
							{
								StratixII* t = (StratixII*) parentOp->getTarget();
								multiplierWidth = t->getDSPMultiplierWidths();
								size = t->getNrDSPMultiplier();	
							}


				for(int i=0; i<size; i++)
					multWidths.push_back(multiplierWidth[i]);

				buildAlteraTiling(0, 0 ,wX ,wY ,0 );




			}

		else  // Xilinx here
			{
				if((!signedIO)&&((wX==41)&&(wY==41)&&(wFull-wOut-g==0)))
					buildFancy41x41Tiling();
				else
					buildXilinxTiling();
			}
	}


	//the fancy tiling is used only for a hardwired case 41 41 82 
	/***********************************************************************/
	void IntMultiplier::buildFancy41x41Tiling()
	{
		//THROWERROR("fancy tiling not implemented yet");

		stringstream inx,iny;
		inx<<addUID("XX");
		iny<<addUID("YY");


		int widthX, widthY,topx,topy;

		//topright dsp;
		widthX=wxDSP;
		widthY=wyDSP;
		topx=0;
		topy=0;
		MultiplierBlock* m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeight);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);


		//topleft dsp
		widthX=wyDSP;
		widthY=wxDSP;
		topx=wxDSP;
		topy=0;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeight);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//bottomleft dsp
		widthX=wxDSP;
		widthY=wyDSP;
		topx=wyDSP;
		topy=wxDSP;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeight);
		m->setNext(NULL);		
		m->setPrevious(NULL);			
		localSplitVector.push_back(m);
		bitHeap->addMultiplierBlock(m);

		//bottomright dsp
		widthX=wyDSP;
		widthY=wxDSP;
		topx=0;
		topy=wyDSP;
		m = new MultiplierBlock(widthX,widthY,topx,topy,inx.str(),iny.str(), lsbWeight);
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

		REPORT(DETAILED,"buildheaplogiconly called for "<<topX<<" "<<topY<<" "<<botX<<" "<<botY);
		Target *target= parentOp->getTarget();
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
		REPORT(DEBUG, "sizeXpadded"<<sizeXPadded);	

		int padX=sizeXPadded-wX;
		int padY=sizeYPadded-wY;

		REPORT(DEBUG, "X split in "<< chunksX << " chunks and Y in " << chunksY << " chunks; ");
		REPORT(DEBUG, " sizeXPadded="<< sizeXPadded << "  sizeYPadded="<< sizeYPadded);
		if (chunksX + chunksY >= 2) { //we do more than 1 subproduct // FIXME where is the else? 

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

			SmallMultTable *tUU, *tSU, *tUS, *tSS;


			// The output size of tUU needs to be one bit larger in case of negate:
			// example for dx=3, dy=3, -0*0 = 0 (positive), -7*7=-49 (negative but sign bit at weight 6)
			tUU = new SmallMultTable( target, dx, dy, dx+dy, false, false, false);

			//useSoftRAM(tUU);
			//oplist.push_back(tUU);
			tUU->addToGlobalOpList();

			if(signedIO) { // need for 4 different tables

				tSU = new SmallMultTable( target, dx, dy, dx+dy, false, true, false );
				//useSoftRAM(tSU);
				// oplist.push_back(tSU);
				tSU->addToGlobalOpList();

				tUS = new SmallMultTable( target, dx, dy, dx+dy, false, false, true );
				//useSoftRAM(tUS);
				// oplist.push_back(tUS);
				tUS->addToGlobalOpList();


				tSS = new SmallMultTable( target, dx, dy, dx+dy, false, true, true );
				//useSoftRAM(tSS);
				//oplist.push_back(tSS);
				tSS->addToGlobalOpList();

			}

			setCycle(0); // TODO FIXME for the virtual multiplier case where inputs can arrive later
			setCriticalPath(initialCP);
			// SmallMultTable is built to cost this:
			manageCriticalPath(  parentOp->getTarget()->localWireDelay(chunksX) +  parentOp->getTarget()->lutDelay() ) ;  
			for (int iy=0; iy<chunksY; iy++){

				vhdl << tab << "-- Partial product row number " << iy << endl;

				for (int ix=0; ix<chunksX; ix++) {

					SmallMultTable *t;
					if (!signedIO) 
						t=tUU;

					else 		{ // 4  cases 
						if( ((ix==chunksX-1)&&(botX==wXX)) && ((iy==chunksY-1)&&(botY==wYY) ))
							t=tSS;
						else if ((ix==chunksX-1)&&(botX==wXX)) 
							t=tSU;
						else if ((iy==chunksY-1)&&(botY==wYY))
							t=tUS;
						else
							t=tUU; 
					}




					//smallMultTable needed only if is on the left size of the truncation 	
					if(dx*(ix+1)+dy*(iy+1)+topX+topY-padX-padY>wFull-wOut-g)
						{
							bitHeap->getPlotter()->addSmallMult(dx*(ix)+topX-padX, dy*(iy)+topY-padY,dx,dy);

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


							// The following comments are obsolete since we negate X at the beginning of the operator:
							// Note that even when negate and tUU, the result is known negative, but it may be zero, so its actual sign is not known statically
							// Note also that in this case (negate and tUU), the result overflows the dx+dy two's complement format.
							// This is why tUU is never negated above, and we negate it below. A bit messy, but probably the most resource-efficient


							bool resultSigned =false;  
							if((t==tSS) || (t==tUS) || (t==tSU)) 
								resultSigned = true ;

							int maxK=t->wOut; // or, dx+dy + (t==tUU && negate?1:0); 
							int minK=0;
							//if(ix == chunksX-1)
							if ((ix == 0))
								minK+=padX;
							if((iy == 0))
								//maxK-=padY;
								minK+=padY;
							REPORT(DEBUG,"The bits will be added from mink="<<minK<<"	to maxk="<<maxK);
							REPORT(DEBUG,  "ix=" << ix << " iy=" << iy << "  maxK=" << maxK  << "  negate=" << negate <<  "  resultSigned="  << resultSigned );
							for (int k=minK; k<maxK; k++) {
								ostringstream s, nots;
								s << PP(ix,iy,blockUid) << of(k); // right hand side
								nots << "not " << s.str(); 

								int weight =  ix*dx+iy*dy+k  +topX+topY-padX-padY  - (wFull-wOut) + lsbWeight;

								if(weight>=0) {// otherwise these bits are just truncated
									if(resultSigned && (k==maxK-1)) { 
										// This is a sign bit: sign-extend it by 1/ complementing and 2/ adding constant 1s all the way up to maxweight
										REPORT(DEBUG, "adding neg bit " << nots.str() << " at weight " << weight); 
										bitHeap->addBit(weight, nots.str());
										REPORT(DEBUG,  "  adding constant ones from weight "<< weight << " to "<< bitHeap->getMaxWeight());
										for (unsigned w=weight; w<bitHeap->getMaxWeight(); w++)
											bitHeap->addConstantOneBit(w);
									}
									else { // just add the bit
										bitHeap->addBit(weight, s.str());
									}
								}
							}

							vhdl << endl;
						}
				}
			}				



		}

	}

	/** checks how many DSPs will be used in case of a tiling **/
	int IntMultiplier::checkTiling(int wxDSP, int wyDSP, int& horDSP, int& verDSP)
	{
		int widthOnX=wX;
		int widthOnY=wY;
		int horDS=0;
		int verDS=0;

		//**** how many dsps will be vertical*******************************/
		int hor=0;

		//if the multiplication is signed, the first DSP will have different size, will be bigger
		if( widthOnX>=wxDSP)
			{
				hor++;
				widthOnX-=wxDSP;
			}

		if(signedIO)	
			wxDSP--;

		//how many DSPs fits on the remaining part, without the first one
		horDS=int(ceil ( (double(widthOnX) / (double) wxDSP)))+hor;
		/***********************************************************************/


		//*** how many dsps will be horizontal**********************************/
		int ver=0;

		if( widthOnY>=wyDSP)
			{
				ver++;
				widthOnY-=wyDSP;
			}


		if(signedIO)	
			wyDSP--;

		verDS=int(ceil ( (double(widthOnY) / (double) wyDSP)))+ver;
		//***********************************************************************/

		horDSP=horDS;
		verDSP=verDS;

		return verDS*horDS;
	}


	//** checks against the ratio the given block and adds a DSP or logic**/
	void IntMultiplier::addExtraDSPs(int topX, int topY, int botx, int boty,int wxDSP,int wyDSP)
	{

		REPORT(DEBUG, "in addExtraDSPs: topX=" << topX << " topY=" << topY << " botX=" << botx << " botY=" << boty);
		int topx=topX,topy=topY;
		//if the block is on the margins of the multipliers, then the coordinates have to be reduced.
		if(topX<0)
			topx=0;
		else
			topx=topX;	

		if(topY<0)
			topy=0;	
		else
			topy=topY;

		//if the truncation line splits the block, then the used block is smaller, the coordinates need to be updated
		if((botx+boty>wFull-wOut-g)&&(topx+topy<wFull-wOut-g))
			{
				int x=topx;
				int y=boty;
				while((x+y<wFull-wOut-g)&&(x<botx))
					{
						x++;
						topx=x;
					}

				x=botx;
				y=topy;
				while((x+y<wFull-wOut-g)&&(y<boty))
					{
						y++;
						topy=y;	
					}	
			}	

		//REPORT(INFO,"in addextradsps after truncation line sets "<< topx <<" "<<topy<<" "<<botx<<" "<<boty);

		//now is the checking against the ratio
		if(checkThreshold(topx,topy,botx,boty,wxDSP,wyDSP))
			{  
				//worth using DSP
				stringstream inx,iny;
				inx<<addUID("XX");
				iny<<addUID("YY");

				topx=botx-wxDSP;
				topy=boty-wyDSP;

				MultiplierBlock* m;
				m = new MultiplierBlock(wxDSP,wyDSP,topx,topy,inx.str(),iny.str(),weightShift);
				m->setNext(NULL);		
				m->setPrevious(NULL);			
				localSplitVector.push_back(m);
				bitHeap->addMultiplierBlock(m);
			}
		else
			{
				//build logic	
				if((topx<botx)&&(topy<boty))
					buildHeapLogicOnly(topx,topy,botx,boty,parentOp->getNewUId());	
			}
	}	



	/** checks the area usage of 1 dsp according to a given block and ratio(threshold)**/
	/** ratio(threshold) = says what percentage of 1 DSP area is allowed to be lost**/
	bool IntMultiplier::checkThreshold(int topX, int topY, int botX, int botY,int wxDSP,int wyDSP)
	{
	
		REPORT(DEBUG, "in checktreshhold "<<topX<<" "<<topY<<" "<<botX<<" "<<botY);
		int widthX=(botX-topX);
		int widthY=(botY-topY);
		double blockArea;
		double triangle=0.0;
		double dspArea=wxDSP*wyDSP;

		//** if the truncation line splits the block, we need to subtract the area of the lost corner**/		
		//***the triangle is the area which will be lost from the block area***//
		//**computing the triangle's edge (45degree => area=(edge^2)/2) ***//		
		int x=topX;
		int y=topY;
		
		while(x+y<wFull-wOut-g)
			x++;
		REPORT(DEBUG, "in checkThreshold: blocArea full "<<widthX*widthY << " x= "<<x<<" topX= "<<topX);
		//computing the triangle's area
		if(topX+topY<=wFull-wOut-g)
			triangle=((x-topX)*(x-topX))/2;
		else
			triangle=0.0;
		//REPORT(INFO, "triangle=" << triangle);
		//the final area which is used
		blockArea=widthX*widthY-triangle;
		REPORT(DEBUG, "* blockArea=" << blockArea);
		//if(( parentOp->getTarget()->getVendor()=="Altera") && 
		//	((wxDSP >=widthX) || (wyDSP >= widthY)))
		//	blockArea -= (wxDSP*wyDSP)  (widthX*widthY);		
		//REPORT(INFO, "**blockArea=" << blockArea);
		//checking according to ratio/area

		if(blockArea>=(1.0-ratio)*dspArea)
				return true;
		else
				return false;
	}



	void IntMultiplier::buildAlteraTiling(int blockTopX, int blockTopY, int blockBottomX, 
	                                      int blockBottomY, int multIndex)
	{
		REPORT(DEBUG, "in Altera tiling");
		int dspSizeX,dspSizeY;

		int size = multWidths.size();

		dspSizeX = multWidths[size-multIndex-1];

		//FIXME
		if(dspSizeX==12)
			dspSizeX=9;
		int sizeLimit;
		if(signedIO)
			{
				sizeLimit=18;
				dspSizeX=18;
			}
		else
			sizeLimit=9;


		dspSizeY=dspSizeX;
		int dsX=dspSizeX;
		int dsY=dspSizeY;

		int width=blockBottomX-blockTopX+1;
		int height=blockBottomY-blockTopY+1;

		int verticalDSP=height/dspSizeY + 1;
		int horizontalDSP=width/dspSizeX + 1;
		//int restY=height-verticalDSP*dspSizeY;
		int botX=blockBottomX;
		int botY=blockBottomY;
		int topX=botX-dspSizeX;
		int topY=botY-dspSizeY;

		REPORT(DEBUG,"verticaldsps= "<<verticalDSP<<" hordsps= "<<horizontalDSP);

		for(int i=0;i<verticalDSP;i++)
			{

				for(int j=0;j<horizontalDSP;j++)
					{
						if (topX<0)
							topX=0;

						if (topY<0)
							topY=0;

						if((signedIO)&&(botX!=wX))
							dspSizeX=dsX-1;
						else
							dspSizeX=dsX;	

						if((signedIO)&&(botY!=wY))
							dspSizeY=dsY-1;
						else
							dspSizeY=dsY;		

						//REPORT(INFO, "in for    dspSizeX=" << dspSizeX);

						if(dspSizeX > sizeLimit)
							{
								//REPORT(INFO, "MAXIMUM CHECKS " << dspSizeX <<" "<< dspSizeY << " " << botX-botY<<" "<<topX-topY);		

								if(checkThreshold(topX, topY, botX, botY, dspSizeX, dspSizeY)
								   && (min(botY-topY+1,botX-topX+1)>=max(dspSizeX,dspSizeY)))
									{	//REPORT(INFO,"22222222222222222222222222222222222222");
										if((topX<botX)&&(topX<botX))		
											addExtraDSPs(topX, topY, botX, botY, dspSizeX, dspSizeY);
									}
								else
									{
										buildAlteraTiling(topX, topY, botX, botY, multIndex+1);
									}
							}
						else
							{
					
								if((topX<botX)&&(topY<botY))
									{
										REPORT(DEBUG, " in BuildAlteraTiling, at " << dspSizeX << " =>  " << topX << " " << topY  << " " << botX << " " << botY << " "  << botX-topX << " " << botY - topY );						
										addExtraDSPs(topX, topY, botX, botY, dspSizeX, dspSizeY);
									}
							}

						botX = topX;

						if((signedIO)&&(botX!=wX))
							dspSizeX=dsX-1;
						else
							dspSizeX=dsX;	

						if((signedIO)&&(botY!=wY))
							dspSizeY=dsY-1;
						else
							dspSizeY=dsY;

						topX = topX-dspSizeX;

					}

				botY = topY;

				if((signedIO)&&(botX!=wX))
					dspSizeX=dsX-1;
				else
					dspSizeX=dsX;	

				if((signedIO)&&(botY!=wY))
					dspSizeY=dsY-1;
				else
					dspSizeY=dsY;

				topY = topY-dspSizeY;
				botX = blockBottomX;

				if((signedIO)&&(botX!=wX))
					dspSizeX=dsX-1;
				else
					dspSizeX=dsX;	

				if((signedIO)&&(botY!=wY))
					dspSizeY=dsY-1;
				else
					dspSizeY=dsY;

				topX = botX-dspSizeX;


			}

	}



	void IntMultiplier::buildXilinxTiling()
	{

		int widthXX,widthX;//local wxDSP
		int widthYY,widthY;//local wyDSP
		int hor1,hor2,ver1,ver2;	
		int horizontalDSP,verticalDSP;
		int nrDSPvertical=checkTiling(wyDSP,wxDSP,hor1,ver1); //number of DSPs used in the case of vertical tiling
		int nrDSPhorizontal=checkTiling(wxDSP,wyDSP,hor2,ver2);//number of DSPs used in case of horizontal tiling
		int botx=wX;
		int boty=wY;
		int topx,topy;

		//decides if a horizontal tiling will be used or a vertical one
		if(nrDSPvertical<nrDSPhorizontal)
			{
				widthXX=wyDSP;
				widthYY=wxDSP;
				horizontalDSP=hor1;
				verticalDSP=ver1;


			}
		else
			{
				widthXX=wxDSP;
				widthYY=wyDSP;
				horizontalDSP=hor2;
				verticalDSP=ver2;
			}


		//applying the tiles
		for(int i=0;i<verticalDSP;i++)
			{	
				//managing the size of a DSP according to its position if signed
				if((signedIO)&&(i!=0))
					widthY=widthYY-1;
				else
					widthY=widthYY;	

				topy=boty-widthY;
				botx=wX;

				for(int j=0;j<horizontalDSP;j++)
					{
						//managing the size of a DSP according to its position if signed
						if((signedIO)&&(j!=0))
							widthX=widthXX-1;
						else
							widthX=widthXX;

						topx=botx-widthX;

						if(botx+boty>wFull-wOut-g)
							addExtraDSPs(topx,topy,botx,boty,widthX,widthY);
						botx=botx-widthX;			
					}


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
		dx(dx), dy(dy), wO(wO), negate(negate), signedX(signedX), signedY(signedY) {
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
		// if(negate && !signedX && !signedY) cerr << "  y=" << y << "  x=" << x;
		r = x * y;
		//  if(negate && !signedX && !signedY) cerr << "  r=" << r;
		if(negate)
			r=-r;
		// if(negate && signedX && signedY) cerr << "  -r=" << r;
		if ( r < 0)
			r += mpz_class(1) << wOut; 
		// if(negate && signedX && signedY) cerr << "  r2C=" << r;

		if(wOut<wF){ // wOut is that of Table
			// round to nearest, but not to nearest even
			int tr=wF-wOut; // number of truncated bits 
			// adding the round bit at half-ulp position
			r += (mpz_class(1) << (tr-1));
			r = r >> tr;
		}

		// if(negate && !signedX && !signedY) cerr << "  rfinal=" << r << endl;

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
