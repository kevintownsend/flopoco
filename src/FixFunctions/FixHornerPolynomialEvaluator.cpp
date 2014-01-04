/*
  
  This file is part of the FloPoCo project
  initiated by the Aric team at Ecole Normale Superieure de Lyon
  and developed by the Socrate team at Institut National des Sciences Appliquées de Lyon
  
  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,  
  2008-2014.
  All rights reserved.

*/
#include <iostream>

#include "FixHornerPolynomialEvaluator.hpp"
#include "FixMultAdd.hpp"

using namespace std;



#define LARGE_PREC 1000 // 1000 bits should be enough for everybody

namespace flopoco{

  FixHornerPolynomialEvaluator::FixHornerPolynomialEvaluator(Target* target, int wX_, int degree_, int LSB_, vector<int> MSB_, bool signedCoeffs_, map<string, double> inputDelays)
    : Operator(target), wX(wX_), degree(degree_), LSB(LSB_), MSB(MSB_), signedCoeffs(signedCoeffs_)
  { 
    
    /* Generate unique name */
    {
      std::ostringstream o;
      o << "FixHornerPolynomialEvaluator_" << wX << "_" << LSB;
      for(int i=0; i<=degree; i++)
	o << "_" << MSB[i];
      o << "_" << (signedCoeffs?1:0);
      if(target->isPipelined()) 
	o << target->frequencyMHz() ;
      else
	o << "comb";
      uniqueName_ = o.str();
    }
    
    setCopyrightString("F. de Dinechin (2014)");
    srcFileName="FixHornerPolynomialEvaluator";


    // declaring inputs
    addInput("X", wX);

    // declaring outputs
    addOutput("R", MSB[0]-LSB+1);


    // First analyse the polynomial coefficients to determine the intermediate sizes
    // As this class doesn't know the coefficients themselves, this is a worst-case analysis
    // This is a simplified version of the computation in the ASAP 2010 paper, because x is in [0,1)
    /* 
       
     */
    for (int i=0; i<degree; i++) {
    }

    // Now assemble faithful FixMultAdd operators

    for (int i=0; i<degree; i++) {
#if 0
      int wInX=wX + (signedCoeffs? 1:0);
      int wInPii=;
	FixMultAdd * ma = new  FixMultAdd(target, 
					  wInX, /* +1 if X is going to be considered signed */
					  ,
					  wA,
					  wR , 
					  msbP, /* msbP */
					  lsbA,
					  true /*signedIO*/, 
					  thresholdForDSP);
	//cout << " Fin du constr" << endl;
	oplist.push_back(ma);
	
	inPortMap ( ma, "X", join("yT",i));
	inPortMap ( ma, "Y", join("sigmaP",i-1));
	inPortMap ( ma, "A", join("a",degree_-i));
	outPortMap( ma, "R", join("sigmaP",i));
	vhdl << instance (ma, join("MultAdd_",i) );
	syncCycleFromSignal( join("sigmaP",i) );
	setCriticalPath( ma->getOutputDelay("R") );
#endif
      }

  }


}
