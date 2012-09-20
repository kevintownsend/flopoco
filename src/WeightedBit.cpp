/*
  A class to manage heaps of weighted bits in FloPoCo
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2012.
  All rights reserved.

*/

#include "WeightedBit.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>	

#include "utils.hpp"
#include <vector>
#include <list>

using namespace std;


namespace flopoco
{

	WeightedBit::WeightedBit(int guid, int uid, int weight_, int type_, int cycle_, double criticalPath_) : 
		cycle(cycle_),  criticalPath(criticalPath_), weight(weight_), type(type_) 
	{
		std::ostringstream p;

		p  << "heap_bh" << guid << "_w" << weight << "_" << uid;
		name=p.str();

		
	}


	WeightedBit::WeightedBit(string name_, int uid_, int weight_, int type_, 
			int cycle_, double criticalPath_) :
	cycle(cycle_), criticalPath(criticalPath_),	weight(weight_), type(type_) 
	{
		name=name_;
		uid=uid_;
		
	}


	
	/* which bit was defined earlier */
	bool WeightedBit::operator< (const WeightedBit& b){
		if ((cycle<b.cycle) || (cycle==b.cycle && criticalPath<b.criticalPath)) 
			return true;
		else
			return false;
	} 

	bool WeightedBit::operator<= (const WeightedBit& b){
		if ((cycle<b.cycle) || (cycle==b.cycle && criticalPath<=b.criticalPath)) 
			return true;
		else
			return false;
	} 
	
	bool WeightedBit::operator== (const WeightedBit& b){
		if ((cycle==b.cycle) && (criticalPath==b.criticalPath)) 
			return true;
		else
			return false;
	} 
	
	double WeightedBit::getCriticalPath(int atCycle)
	{
#if 0
		if (cycle>atCycle){
			THROWERROR("For bit " << name << ", getCriticalPath() called for cycle "<<atCycle<< " but this bit is defined only at cycle "<< cycle);
		}
#endif
		if (cycle==atCycle)
			return criticalPath;
		if (cycle<atCycle)
			return 0.0;
		
		return 0.0;  //because it returned no value on this branch
	}

	int WeightedBit::computeStage(int stagesPerCycle, double elementaryTime)
	{
		return (getCycle()*stagesPerCycle + getCriticalPath(getCycle())/elementaryTime);        
	}

	
}
