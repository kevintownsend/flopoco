/*
Tools for performing floorplaning

Author : Matei Istoan, Florent de Dinechin

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
2012.
  All rights reserved.

*/


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "Operator.hpp"
#include "utils.hpp"
#include "FloorplanningHelper.hpp"


namespace flopoco{
	

	FloorplanningHelper::FloorplanningHelper(Target* target_, Operator* op_){
		
		target = target_;
		parentOp = op_;
	}
	
	void FloorplanningHelper::addToFlpComponentList(std::string name){
		
		flComponentList.push_back(name);
	}
	
	void FloorplanningHelper::addToInstanceNames(std::string componentName, std::string instanceName){
		
		flInstanceNames[componentName] = instanceName;
	}
	
	void FloorplanningHelper::initFloorplanning(double ratio){
		
		floorplanningRatio = ratio;
		
		virtualModuleId = 0;
		
		prevEstimatedCountFF 		 = 0;
		prevEstimatedCountLUT 		 = 0;
		prevEstimatedCountMemory 	 = 0;
		prevEstimatedCountMultiplier = 0;
	}
	
	std::string FloorplanningHelper::manageFloorplan(){
		//create a new virtual module for the resources that are between
		//	two modules, and add it to the corresponding lists
		std::ostringstream result;
		Operator* virtualComponent = new Operator(target); 
		std::string moduleName = join("virtual_module_", virtualModuleId);
		
		result << "";
		
		virtualModuleId++;
		
		virtualComponent->setName(moduleName);
		virtualComponent->reHelper->estimatedCountLUT 			= parentOp->reHelper->estimatedCountLUT 	   - prevEstimatedCountLUT;
		virtualComponent->reHelper->estimatedCountFF  			= parentOp->reHelper->estimatedCountFF 		   - prevEstimatedCountFF;
		virtualComponent->reHelper->estimatedCountMultiplier 	= parentOp->reHelper->estimatedCountMultiplier - prevEstimatedCountMultiplier;
		virtualComponent->reHelper->estimatedCountMemory 		= parentOp->reHelper->estimatedCountMemory 	   - prevEstimatedCountMemory;
		
		if((virtualComponent->reHelper->estimatedCountLUT!=0) || (virtualComponent->reHelper->estimatedCountFF!=0) 
				|| (virtualComponent->reHelper->estimatedCountMultiplier!=0) || (virtualComponent->reHelper->estimatedCountMemory!=0)){
			flComponentList.push_back(moduleName);
			flVirtualComponentList[moduleName] = virtualComponent;
			
			result << "Created virtual module - " << moduleName << endl;
			result << tab << "Added " << parentOp->reHelper->estimatedCountLUT		  - prevEstimatedCountLUT << " function generators" << endl;
			result << tab << "Added " << parentOp->reHelper->estimatedCountFF		  - prevEstimatedCountFF << " registers" << endl;
			result << tab << "Added " << parentOp->reHelper->estimatedCountMultiplier - prevEstimatedCountMultiplier << " multipliers" << endl;
			result << tab << "Added " << parentOp->reHelper->estimatedCountMemory	  - prevEstimatedCountMemory << " memories" << endl;
			
			prevEstimatedCountLUT 			= parentOp->reHelper->estimatedCountLUT;
			prevEstimatedCountFF 			= parentOp->reHelper->estimatedCountFF;
			prevEstimatedCountMultiplier 	= parentOp->reHelper->estimatedCountMultiplier;
			prevEstimatedCountMemory 		= parentOp->reHelper->estimatedCountMemory;
		}
		
		return result.str();
	}
	
	std::string FloorplanningHelper::addPlacementConstraint(std::string source, std::string sink, int type){
		std::ostringstream result;
		constraintType newConstraint;
		map<string, Operator*> subComponents = parentOp->subComponents_;
		
		//check to see if the type of the constraint is valid
		if(!((type==TO_LEFT_OF) || (type==TO_RIGHT_OF) || (type==ABOVE) || (type==UNDER) || 
				(type==TO_LEFT_OF_WITH_EXTRA) || (type==TO_RIGHT_OF_WITH_EXTRA) || (type==ABOVE_WITH_EXTRA) || (type==UNDER_WITH_EXTRA))){
			cerr << "Error: Trying to add a placement constraint of undefined type." << endl;
			exit(1);
		}
		
		//check if the source is a valid sub-component name
		if(subComponents.find(source)==subComponents.end()){
			cerr << "Error: source sub-component " << source << " was not found" << endl;
			exit(1);
		}
		
		//check if the sink is a valid sub-component name
		if(subComponents.find(sink)==subComponents.end()){
			cerr << "Error: sink sub-component " << sink << " was not found" << endl;
			exit(1);
		}
		
		newConstraint.type 		= PLACEMENT;
		newConstraint.source 	= source;
		newConstraint.sink 		= sink;
		newConstraint.value 	= type;
		
		flConstraintList.push_back(newConstraint);
		
		result << "Created new placement constraint" << endl;
		result << tab << "Sink module " << sink
				<< " is " << ((type==TO_LEFT_OF) ? "to the left of" 							: 
							(type==TO_RIGHT_OF) ? "to the right of" 							: 
							(type==ABOVE) ? "above" 											: 
							(type==UNDER) ? "under" 											: 
							(type==TO_LEFT_OF_WITH_EXTRA) ? "to the left of (with glue logic)"  : 
							(type==TO_RIGHT_OF_WITH_EXTRA) ? "to the right of (with glue logic)": 
							(type==ABOVE_WITH_EXTRA) ? "above (with glue logic)" 				: "under (with glue logic)")
				<< " source module " << source << endl;
		
		return result.str();
	}
	
	std::string FloorplanningHelper::addConnectivityConstraint(std::string source, std::string sink, int nrWires){
		std::ostringstream result;
		constraintType newConstraint;
		map<string, Operator*> subComponents = parentOp->subComponents_;
		
		//no non-positive values allowed for the number of wires
		if(nrWires<1){
			cerr << "Error: trying to add an invalid number of wires:" << nrWires << endl;
			exit(1);
		}
		
		//check if the source is a valid sub-component name
		if(subComponents.find(source)==subComponents.end()){
			cerr << "Error: source sub-component " << source << " was not found" << endl;
			exit(1);
		}
		
		//check if the sink is a valid sub-component name
		if(subComponents.find(sink)==subComponents.end()){
			cerr << "Error: sink sub-component " << sink << " was not found" << endl;
			exit(1);
		}
		
		newConstraint.type 		= CONNECTIVITY;
		newConstraint.source 	= source;
		newConstraint.sink 		= sink;
		newConstraint.value 	= nrWires;
		
		flConstraintList.push_back(newConstraint);
		
		result << "Created new connectivity constraint" << endl;
		result << tab << "Source module " << source 
				<< "is connected to sink module " << sink << " by " << nrWires << " wires" << endl;
		
		return result.str();
	}
	
	std::string FloorplanningHelper::processPlacementConstraints(){
		std::ostringstream result;
		
		//process each constraint, one by one, and update the virtual grid
		for(unsigned int i=0; i<flConstraintList.size(); i++){
			constraintType newConstraint = flConstraintList[i];
			coordinateType sourceCoordinate, sinkCoordinate, newCoordinate;
			bool positionOccupied = false;
			
			result << "Processed placement constraint" << endl;
			
			//get coordinates of source and sink modules
			sourceCoordinate = flComponentCordVirtual[newConstraint.source];
			sinkCoordinate = flComponentCordVirtual[newConstraint.sink];
			if(newConstraint.type==PLACEMENT){
				if(newConstraint.value==TO_LEFT_OF || newConstraint.value==TO_LEFT_OF_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x - 1;
					newCoordinate.y = sourceCoordinate.y;
					
					if(sinkCoordinate.x<sourceCoordinate.x && sinkCoordinate.y==sourceCoordinate.y)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						unsigned int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.x -= 1;
						}
					}while(positionOccupied == true);
					
					if(newConstraint.value==TO_LEFT_OF_WITH_EXTRA){
						//look for the virtual module just before the sink module, if it exists
						// if it does, then place it between the source and the sink
						std::string prevModuleName = "";
						
						for(unsigned int i=0; i<flComponentList.size(); i++){
							if(newConstraint.sink == flComponentList[i])
								break;
							else
								prevModuleName = flComponentList[i];
						}
						if(prevModuleName.find("virtual_module_") != string::npos){
							flComponentCordVirtual[prevModuleName] = newCoordinate;
							newCoordinate.x -= 1;
						}
					}
				}else if(newConstraint.value==TO_RIGHT_OF || newConstraint.value==TO_RIGHT_OF_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x + 1;
					newCoordinate.y = sourceCoordinate.y;
					
					if(sinkCoordinate.x>sourceCoordinate.x && sinkCoordinate.y==sourceCoordinate.y)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						unsigned int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.x += 1;
						}
					}while(positionOccupied == true);
					
					if(newConstraint.value==TO_RIGHT_OF_WITH_EXTRA){
						//look for the virtual module just before the sink module, if it exists
						// if it does, then place it between the source and the sink
						std::string prevModuleName = "";
						
						for(unsigned int i=0; i<flComponentList.size(); i++){
							if(newConstraint.sink == flComponentList[i])
								break;
							else
								prevModuleName = flComponentList[i];
						}
						if(prevModuleName.find("virtual_module_") != string::npos){
							flComponentCordVirtual[prevModuleName] = newCoordinate;
							newCoordinate.x += 1;
						}
					}
				}else if(newConstraint.value==ABOVE || newConstraint.value==ABOVE_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x;
					newCoordinate.y = sourceCoordinate.y - 1;
					
					if(sinkCoordinate.y<sourceCoordinate.y && sinkCoordinate.x==sourceCoordinate.x)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						unsigned int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.y -= 1;
						}
					}while(positionOccupied == true);
					
					if(newConstraint.value==ABOVE_WITH_EXTRA){
						//look for the virtual module just before the sink module, if it exists
						// if it does, then place it between the source and the sink
						std::string prevModuleName = "";
						
						for(unsigned int i=0; i<flComponentList.size(); i++){
							if(newConstraint.sink == flComponentList[i])
								break;
							else
								prevModuleName = flComponentList[i];
						}
						if(prevModuleName.find("virtual_module_") != string::npos){
							flComponentCordVirtual[prevModuleName] = newCoordinate;
							newCoordinate.y -= 1;
						}
					}
				}else if(newConstraint.value==UNDER || newConstraint.value==UNDER_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x;
					newCoordinate.y = sourceCoordinate.y + 1;
					
					if(sinkCoordinate.y>sourceCoordinate.y && sinkCoordinate.x==sourceCoordinate.x)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						unsigned int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.y += 1;
						}
					}while(positionOccupied == true);
					
					if(newConstraint.value==UNDER_WITH_EXTRA){
						//look for the virtual module just before the sink module, if it exists
						// if it does, then place it between the source and the sink
						std::string prevModuleName = "";
						
						for(unsigned int i=0; i<flComponentList.size(); i++){
							if(newConstraint.sink == flComponentList[i])
								break;
							else
								prevModuleName = flComponentList[i];
						}
						if(prevModuleName.find("virtual_module_") != string::npos){
							flComponentCordVirtual[prevModuleName] = newCoordinate;
							newCoordinate.y += 1;
						}
					}
				}
				
				//change the coordinates of the sink module
				flComponentCordVirtual[newConstraint.sink] = newCoordinate;
				
				result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
				result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
				result << tab << tab 
						<< ((newConstraint.value==TO_LEFT_OF  || newConstraint.value==TO_LEFT_OF_WITH_EXTRA)  ? " moving under module " :
						   (newConstraint.value==TO_RIGHT_OF || newConstraint.value==TO_RIGHT_OF_WITH_EXTRA) ? " moving under module " : 
						   (newConstraint.value==ABOVE 		 || newConstraint.value==ABOVE_WITH_EXTRA) 		 ? " moving under module " : 
						   (newConstraint.value==UNDER 		 || newConstraint.value==UNDER_WITH_EXTRA) 		 ? " moving under module " : "")
						<< newConstraint.source << endl;
				result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
			}
		}
		
		//re-normalize the virtual grid
		//the constraints might have created wholes in the grid, 
		//	leading to inefficiency and unused space
		
		//scan the components' y coordinates; if there are negative 
		//	coordinates, shift all coordinates downwards, so as to
		//	have all coordinates non-negative
		int minLine = 0;
		
		for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
			if((it->second).y < minLine)
				minLine = (it->second).y;
		}
		if(minLine < 0){
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				(it->second).y += abs(minLine);
			}
		}
		
		//scan the components' y coordinates; if there are empty lines,
		//	move the rest of the lines up
		bool hasEmptyLines = false;
		int maxLine = 0;
		
		//compute the lowest line
		for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
			coordinateType tempCoordinate = it->second;
			if(tempCoordinate.y>maxLine)
				maxLine = tempCoordinate.y;
		}
		
		//as long as the closest line to the current is more than a line 
		//	away, move the line upwards so as to fill the gaps
		do{
			int minGapSize = flComponentCordVirtual.size(), minGapLine = 0;
			
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				coordinateType tempCoordinate = it->second;
				coordinateType closestCoordinate = it->second;
				
				if(tempCoordinate.y == maxLine)
					continue;
				
				closestCoordinate.y +=flComponentCordVirtual.size();
				
				for(map<string, coordinateType>::iterator it2 = flComponentCordVirtual.begin(); it2 != flComponentCordVirtual.end(); it2++){
					coordinateType tempCoordinate2 = it2->second;
					
					if((tempCoordinate2.y > tempCoordinate.y) && (tempCoordinate2.y < closestCoordinate.y) && (it->first != it2->first))
						closestCoordinate = tempCoordinate2;
				}
				
				if((closestCoordinate.y-tempCoordinate.y < minGapSize) && (closestCoordinate.y-tempCoordinate.y > 1)){
					minGapSize = closestCoordinate.y-tempCoordinate.y;
					minGapLine = tempCoordinate.y;
				}
			}
			
			if(minGapSize>1 && minGapSize!=flComponentCordVirtual.size()){
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).y >= minGapLine)
						flComponentCordVirtual[(it->first)].y -= (minGapSize-1);
				}
			}else{
				hasEmptyLines = false;
			}
		}while(hasEmptyLines == true);
		
		//scan the components' x coordinates; if there are negative 
		//	coordinates, shift all coordinates to the right, so as to
		//	have all coordinates non-negative
		int minColumn = 0;
		
		for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
			if((it->second).x < minColumn)
				minColumn = (it->second).x;
		}
		if(minColumn < 0){
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				(it->second).x += abs(minColumn);
			}
		}
		
		return result.str();
	}
	
	std::string FloorplanningHelper::processConnectivityConstraints(){
		//
		 // currently, as the user decides where each module goes, without
		 // the automation of the process, the relevance of the number of
		 // links between two modules is questionable.
		 // futrher modifications to follow.
		 //
		 return "";
	}
	
	std::string FloorplanningHelper::createVirtualGrid(){
		std::ostringstream result;
		
		result << "Created virtual grid of sub-components" << endl;
		
		//create the virtual component grid, and place the modules one
		//	under the other, in the order that they are instantiated
		//virtual modules, for the glue logic, are also considered as
		//	modules
		for(unsigned int i=0; i<flComponentList.size(); i++){
			std::string componentName = flComponentList[i];
			coordinateType newCoordinate;
			
			newCoordinate.x = 0;
			newCoordinate.y = i;
			
			flComponentCordVirtual[componentName] = newCoordinate;
			
			result << tab << "Sub-component " << componentName << " placed by default at x=0 and y=" << i << endl;
		}
		
		return result.str();
	}
	
	/*
	bool FloorplanningHelper::createHorrizontalPlacement(){
		
	}
	
	bool FloorplanningHelper::createVerticalPlacement(){
		
	}
	*/
	
	std::string FloorplanningHelper::createPlacementGrid(){
		ostringstream result;
		map<string, Operator*> subComponents = parentOp->subComponents_;
		
		result << "Created real grid of components" << endl;
		
		//different processing techniques and efforst for modules that
		//	contain and that don't contain DSP and RAM blocks
		if(parentOp->reHelper->estimatedCountMultiplier!=0 || parentOp->reHelper->estimatedCountMultiplier!=0 
				|| parentOp->reHelper->estimatedCountMemory!=0){
			int nrLines = 0, nrColumns = 0;
			int maxDSPperLine = 0, maxDSPperColumn = 0, maxRAMperLine = 0, maxRAMperColumn = 0;
			vector< vector<string> > subcomponentMatrixLines, subcomponentMatrixColumns;
			//FEATURE NOT YET WORKING
			//bool invertAxes = false, breakLines = false;
			//int virtualDSPperColumn, virtualDSPperRow, virtualRAMperColumn, virtualRAMperRow;
			int maxComponentHeight = 0;
			
			result << tab << "Creating the grid of real coordintes with DSP or/and RAM requirements" << endl;
			
			//determine the maximum number of DSPs in one line of the floorplan
			//	if the nb. is larger the number of lines in the chip, throw warning to the user,
			//	then check to see the maximum number of DSPs in one column
			//-> if viable, invert axes and continue with the floorplan
			//		if not, break the lines longer than the capacity
			
			result << tab << "Created a (almost) 2D version of the virtual component grid, organized by lines" << endl;
			
			//construct a (almost) 2D version of the virtual grid, organized by lines
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				if((it->second).y > nrLines)
					nrLines = (it->second).y;
			}
			for(int i=0; i<=nrLines; i++){
				vector<string> tempLevelList;
				int nrDSP = 0, nrRAM = 0;
				
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).y == i){
						tempLevelList.push_back(it->first);
						
						Operator* tempOperator = subComponents[it->first];
						nrDSP += (tempOperator->reHelper->estimatedCountMultiplier > 0) 	? 1 : 0;
						nrRAM += (tempOperator->reHelper->estimatedCountMemory > 0) 		? 1 : 0;
					}
				}
				subcomponentMatrixLines.push_back(tempLevelList);
				
				if(nrDSP > maxDSPperLine)
					maxDSPperLine = nrDSP;
				if(nrRAM > maxRAMperLine)
					maxRAMperLine = nrRAM;
			}
			
			result << tab << "Created a (almost) 2D version of the virtual component grid, organized by columns" << endl;
			
			//construct a (almost) 2D version of the virtual grid, organized by columns
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				if((it->second).y > nrLines)
					nrColumns = (it->second).x;
			}
			for(int i=0; i<=nrColumns; i++){
				vector<string> tempLevelList;
				int nrDSP = 0, nrRAM = 0;
				
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).x == i){
						tempLevelList.push_back(it->first);
						
						Operator* tempOperator = subComponents[it->first];
						nrDSP += tempOperator->reHelper->estimatedCountMultiplier;
						nrRAM += tempOperator->reHelper->estimatedCountMemory;
					}
				}
				subcomponentMatrixColumns.push_back(tempLevelList);
				
				if(nrDSP > maxDSPperColumn)
					maxDSPperColumn = nrDSP;
				if(nrRAM > maxRAMperColumn)
					maxRAMperColumn = nrRAM;
			}
			
			result << tab << "Testing if columns and lines should/shouldn't be inverted" << endl;
			
			//FEATURE UNDER TEST - NOT YET USED
			//determine if the axes should be inverted and if the lines
			//	should be broken in two
			//if(maxDSPperLine > (target->multiplierPosition).size()){
				//if(maxDSPperColumn < (target->multiplierPosition).size())
					//invertAxes = true;
				//else
					//breakLines = true;
			//}
			//if(maxRAMperLine > (target->memoryPosition).size()){
				//if(maxRAMperColumn < (target->memoryPosition).size())
					//invertAxes = true;
				//else
					//breakLines = true;
			//}
			
			//if(invertAxes == true){
				//for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					//int temp = (it->second).x;
					//(it->second).x = (it->second).y;
					//(it->second).y = temp;
				//}
				
				//vector< vector<string> > tempMatrix;
				//tempMatrix = subcomponentMatrixColumns;
				//subcomponentMatrixColumns = subcomponentMatrixLines;
				//subcomponentMatrixLines = tempMatrix;
				
				//virtualDSPperColumn = (target_->multiplierPosition).size();
				//virtualDSPperRow = target_->dspPerColumn;
				
				//virtualRAMperColumn = (target_->memoryPosition).size();
				//virtualRAMperRow = target_->ramPerColumn;
				
				//int temp = maxDSPperColumn;
				//maxDSPperColumn = maxDSPperLine;
				//maxDSPperLine = temp;
				//temp = maxRAMperColumn;
				//maxRAMperColumn = maxRAMperLine;
				//maxRAMperLine = temp;
			//}else{
				//virtualDSPperColumn = target_->dspPerColumn;
				//virtualDSPperRow = (target_->multiplierPosition).size();
				
				//virtualRAMperColumn = target_->ramPerColumn;
				//virtualRAMperRow = (target_->memoryPosition).size();
			//}
			//if(breakLines == true){
				////should break the lines, but this may lead to considerable
				////	differences from the original floorplan
				////  so an Error is signaled to the user to reconsider the strategy
				//cerr << "Warning: the current floorplanning strategy, on the target FPGA, cannot be performed." 
						//<< " Please reconsider. No floorplan has been generated." << endl;
			//}
			//
			
			//determine the maximum height of a bounding box for a 
			//	sub-component based on the largets chain of DSPs or RAMs 
			//	in any of the sub-components
			for(map<string, Operator*>::iterator it = subComponents.begin(); it != subComponents.end(); it++){
				Operator* tempOperator = it->second;
				int multiplierHeightRatio, ramHeightRatio;
				
				multiplierHeightRatio = tempOperator->reHelper->estimatedCountMultiplier/target->dspHeightInLUT;
				ramHeightRatio 		  = tempOperator->reHelper->estimatedCountMemory/target->ramHeightInLUT;
				if(multiplierHeightRatio > maxComponentHeight)
					maxComponentHeight = multiplierHeightRatio;
				if(ramHeightRatio > maxComponentHeight)
					maxComponentHeight = ramHeightRatio;
			}
			
			result << tab << "Ordering the elements of the lines in the order that they appear in the virtual grid" << endl;
			
			//TODO: version not working properly; take ordering from logic-only version
			
			//order the lines of the 2D virtual placement matrix, in
			//	increasing order of their x coordinate
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> tempList = subcomponentMatrixLines[i];
				
				for(unsigned int j=0; j<tempList.size(); j++)
					for(unsigned int k=0; k<tempList.size(); k++){
						if((flComponentCordVirtual[tempList[j]]).x > (flComponentCordVirtual[tempList[k]]).x){
							string tempString = tempList[j];
							tempList[j] = tempList[k];
							tempList[k] = tempString;
						}
					}
					
				subcomponentMatrixLines[i] = tempList;
			}
			
			result << tab << "Creating the grid of real coordintes" << endl;
			
			//go through the 2D virtual placement list and generate the
			//	real coordinate list
			int prevCoordX, prevCoordY;
			
			result << tab << "Creating the initial placement in the real coordinate grid" << endl;
			
			//create an initial placement, then move the components around
			//	so as to satisfy their resource requirements
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> matrixLine = subcomponentMatrixLines[i];
				
				prevCoordX = 0;
				prevCoordY = i * maxComponentHeight * (1.0/floorplanningRatio);
				for(unsigned int j=0; j<matrixLine.size(); j++){
					int lutWidth, ffWidth, componentWidth;
					coordinateType componentLocation, componentDimension;
					Operator* currentComponent = subComponents[matrixLine[j]];
					
					lutWidth = (currentComponent->reHelper->estimatedCountLUT/target->lutPerSlice)/maxComponentHeight;		//width in slices
					ffWidth = (currentComponent->reHelper->estimatedCountFF/target->ffPerSlice)/maxComponentHeight;		//width in slices
					componentWidth = (lutWidth>ffWidth) ? lutWidth : ffWidth;
					
					componentLocation.x = prevCoordX + 1;
					componentLocation.y = prevCoordY;
					componentDimension.x = componentWidth * (1.0/floorplanningRatio);
					componentDimension.y = maxComponentHeight * (1.0/floorplanningRatio);
					
					flComponentCordReal[matrixLine[j]] = componentLocation;
					flComponentDimension[matrixLine[j]] = componentDimension;
					
					prevCoordX = prevCoordX + 1 + componentDimension.x;
				}
			}
			
			result << tab << "Adjusting the component grid so that components with RAM and DSP requirement are properly placed" << endl;
			
			//shift the components to the left so that they meet the
			//	DSP and RAM requirements
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> matrixLine = subcomponentMatrixLines[i];
				
				for(unsigned int j=0; j<matrixLine.size(); j++){
					bool dspSatisfied = false, ramSatisfied = false;
						
					for(unsigned int k=0; k<target->multiplierPosition.size(); k++){
						if(target->multiplierPosition[k] == flComponentCordReal[matrixLine[j]].x+1){
							dspSatisfied = true;
							break;
						}
					}
					for(unsigned int k=0; k<target->memoryPosition.size(); k++){
						if(target->memoryPosition[k] == flComponentCordReal[matrixLine[j]].x+1){
							dspSatisfied = true;
							break;
						}
					}
					
					if(dspSatisfied && ramSatisfied)
						continue;
					else if(!dspSatisfied && !ramSatisfied){
						int closestDSPColumn = flComponentCordReal[matrixLine[j]].x, closestRAMColumn = flComponentCordReal[matrixLine[j]].x;
						int closestMove, dimensionIncrease;
						
						//decide which is closest, DSP or RAM, then shift to that column, 
						//	and then extend the width to accomodate the other resource
						for(unsigned int k=0; k<target->multiplierPosition.size(); k++){
							if(target->multiplierPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestDSPColumn = target->multiplierPosition[k]+1;
								break;
							}
						}
						for(unsigned int k=0; k<target->memoryPosition.size(); k++){
							if(target->memoryPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestRAMColumn = target->memoryPosition[k]+1;
								break;
							}
						}
						
						closestMove = (closestDSPColumn>closestRAMColumn) ? closestDSPColumn : closestRAMColumn;
						dimensionIncrease += abs(closestDSPColumn - closestRAMColumn);
						//if needed, shift all the components to the left
						if((closestMove - flComponentCordReal[matrixLine[j]].x) > 0){
							for(unsigned int k=j; k<matrixLine.size(); k++){
								flComponentCordReal[matrixLine[k]].x += dimensionIncrease + 
																			(closestMove - flComponentCordReal[matrixLine[j]].x);
							}
						}
					}else if(!dspSatisfied){
						int closestDSPColumn = flComponentCordReal[matrixLine[j]].x;
						
						for(unsigned int k=0; k<target->multiplierPosition.size(); k++){
							if(target->multiplierPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestDSPColumn = target->multiplierPosition[k]+1;
								break;
							}
						}
						
						//if needed, shift all the components to the left
						int shiftSize = closestDSPColumn - flComponentCordReal[matrixLine[j]].x;
						if(shiftSize > 0){
							for(unsigned int k=j; k<matrixLine.size(); k++){
								flComponentCordReal[matrixLine[k]].x += shiftSize;
							}
						}
					}else if(!ramSatisfied){
						int closestRAMColumn = flComponentCordReal[matrixLine[j]].x;
						
						for(unsigned int k=0; k<target->memoryPosition.size(); k++){
							if(target->memoryPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestRAMColumn = target->memoryPosition[k]+1;
								break;
							}
						}
						
						//if needed, shift all the components to the left
						int shiftSize = closestRAMColumn - flComponentCordReal[matrixLine[j]].x;
						if(shiftSize > 0){
							for(unsigned int k=j; k<matrixLine.size(); k++){
								flComponentCordReal[matrixLine[k]].x += shiftSize;
							}
						}
					}
				}
			}
			
		}else{
			//creating the placement for operators that don not have DSPs and RAMs
			result << tab << "Creating the grid of real coordintes" << endl;
			
			int nrLines = 0;
			vector< vector<string> > subcomponentMatrixLines;
			int maxComponentHeight = 0;
			
			result << tab << "Construct a (almost) 2D version of the virtual grid, organized by lines" << endl;
			
			//construct a (almost) 2D version of the virtual grid, organized by lines
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				if((it->second).y > nrLines)
					nrLines = (it->second).y;
			}
			for(int i=0; i<=nrLines; i++){
				vector<string> tempLevelList;
				
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).y == i){
						tempLevelList.push_back(it->first);
					}
				}
				subcomponentMatrixLines.push_back(tempLevelList);
			}
			
			//determine the maximum height of a bounding box of any of 
			//	the sub-components; the goal is to have balanced dimensions 
			//	for the boxes, so make the largest box as close to being 
			//	square -> hence the square root
			for(unsigned int i=0; i<flComponentList.size(); i++){
				Operator* tempOperator;
				int operatorWidth; 
				
				if(subComponents.find(flComponentList[i]) != subComponents.end()){
					tempOperator = subComponents[flComponentList[i]];
				}else{
					tempOperator = flVirtualComponentList[flComponentList[i]];
				}
				
				operatorWidth = ceil(sqrt((tempOperator->reHelper->estimatedCountLUT)/(double)target->lutPerSlice));
				if(operatorWidth>maxComponentHeight)
					maxComponentHeight = operatorWidth;
				operatorWidth = ceil(sqrt((tempOperator->reHelper->estimatedCountFF)/(double)target->ffPerSlice));
				if(operatorWidth>maxComponentHeight)
					maxComponentHeight = operatorWidth;
			}
			
			//if the modules chain vertically and they run out of space, expand design horizontally
			while(ceil(nrLines*maxComponentHeight*sqrt(1.0/floorplanningRatio)) > target->topSliceY){
				maxComponentHeight--;
				
				//design cannot be floorplanned with the current constraints
				if(maxComponentHeight == 0){
					cerr << "Error: the design cannot be floorplanned with the current constraints. Please reconsider and re-run." << endl;
					cerr << tab << "number of lines= " << nrLines << " maximum component height=" << maxComponentHeight << " ratio=" << 1.0/floorplanningRatio << endl;
					cerr << tab << "maximum allowed height=" << target->topSliceY << endl;
					exit(1);
				}
			}
			
			result << tab << "Order the lines of the 2D virtual placement matrix" << endl;
			
			//order the lines of the 2D virtual placement matrix, in the
			//	increasing order of their x coordinates
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> tempList = subcomponentMatrixLines[i];
				
				for(unsigned int j=0; j<tempList.size()-1; j++)
					for(unsigned int k=j+1; k<tempList.size(); k++){
						if((flComponentCordVirtual[tempList[j]]).x > (flComponentCordVirtual[tempList[k]]).x){
							coordinateType tempCoord = flComponentCordVirtual[tempList[j]];
							flComponentCordVirtual[tempList[j]] = flComponentCordVirtual[tempList[k]];
							flComponentCordVirtual[tempList[k]] = tempCoord;
							
							string tempString = tempList[j];
							tempList[j] = tempList[k];
							tempList[k] = tempString;
						}
					}
					
				subcomponentMatrixLines[i] = tempList;
			}
			
			result << tab << "Generate the real coordinate grid" << endl;
			
			//go through the 2D virtual placement list and generate the
			//	real coordinate list
			int prevCoordX, prevCoordY;
			
			//place each line at a time; lines with lower y coordinates 
			//	are processed first, and inside the lines the elements 
			// are processed in the increasing order of their x coordinates
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> matrixLine = subcomponentMatrixLines[i];
				
				prevCoordX = 0;
				prevCoordY = ceil(i * maxComponentHeight * sqrt(1.0/floorplanningRatio));
				for(unsigned int j=0; j<matrixLine.size(); j++){
					int lutWidth, ffWidth, componentWidth;
					coordinateType componentLocation, componentDimension;
					Operator* currentComponent;
					
					if(subComponents.find(matrixLine[j]) != subComponents.end()){
						currentComponent = subComponents[matrixLine[j]];
						
						lutWidth = ceil(((double)(currentComponent->reHelper->estimatedCountLUT)/(target->lutPerSlice))/maxComponentHeight);		//width in slices
						ffWidth = ceil(((double)(currentComponent->reHelper->estimatedCountFF)/(target->ffPerSlice))/maxComponentHeight);		//width in slices
						componentWidth = (lutWidth>ffWidth) ? lutWidth : ffWidth;
						
						componentLocation.x = prevCoordX + 1;
						componentLocation.y = prevCoordY;
						componentDimension.x = ceil(componentWidth * sqrt(1.0/(double)floorplanningRatio));
						componentDimension.y = ceil(maxComponentHeight * sqrt(1.0/(double)floorplanningRatio));
						
						flComponentCordReal[matrixLine[j]] = componentLocation;
						flComponentDimension[matrixLine[j]] = componentDimension;
						
						prevCoordX = prevCoordX + 1 + componentDimension.x;
					}
				}
			}
		}
		
		return result.str();
	}
	
	std::string FloorplanningHelper::createConstraintsFile(){
		ofstream file;
		ostringstream result;
		
		result << "Created output constraints file" << endl;
		
		//create the physical file
		if(target->getVendor() == "Xilinx")
			file.open("flopoco.ucf");
		else if(target->getVendor() == "Altera"){
		
		}
		
		result << tab << "Adding constraints" << endl;
		result << tab << "Added constraints to contain the entire operator" << endl;
		file << createPlacementForComponent("root");
		
		//create the placement constraints for each sub-component at a time
		//this time, the boxes for the glue logic aren't gien any bounds, 
		//	but they will have the empty space between the modules, which 
		//	should fit them
		for(unsigned int i=0; i<flComponentList.size(); i++){
			result << tab << "Added constraints for sub-component " << flComponentList[i] << endl;
			file << createPlacementForComponent(flComponentList[i]);
		}
		
		file.close();
		
		if(target->getVendor() == "Xilinx"){
			cerr << "***Floorplan written to \'flopoco.ucf\' constraints file" << endl;
		}else if(target->getVendor() == "Altera"){
			cerr << "ERROR: Floorplanning feature not yet implemented for Altera target" << endl;
		}else{
			cerr << "ERROR: Floorplanning not yet implemented for this target" << endl;
		}
		
		return result.str();
	}
	
	std::string FloorplanningHelper::createPlacementForComponent(std::string moduleName){
		ostringstream constraintString;
		map<string, Operator*> subComponents = parentOp->subComponents_;
		
		constraintString << "";
		
		//create the constraint for the whole operator
		if(moduleName == "root"){
			
			int maxX, maxY, minX, minY;
			map<string, coordinateType>::iterator it = flComponentCordReal.begin();
			
			minX = (it->second).x;
			minY = (it->second).y;
			
			maxX = minX + (flComponentDimension[it->first]).x;
			maxY = minY + (flComponentDimension[it->first]).y;
			
			it++;
			while(it != flComponentCordReal.end()){
				if((it->second).x < minX){
					minX = (it->second).x;
				}
				if((it->second).y < minY){
					minY = (it->second).y;
				}
				if((it->second).x+(flComponentDimension[it->first]).x > maxX){
					maxX = (it->second).x+(flComponentDimension[it->first]).x;
				}
				if((it->second).y+(flComponentDimension[it->first]).y > maxY){
					maxY = (it->second).y+(flComponentDimension[it->first]).y;
				}
				it++;
			}
			
			constraintString << "INST \"*\" AREA_GROUP=\"pblock_root\";" << endl;
			constraintString << "AREA_GROUP \"pblock_root\" RANGE=SLICE_X" 
				<< minX << "Y" << minY << ":SLICE_X" << maxX << "Y" << maxY << ";" << endl;
			constraintString << endl;
				
			return constraintString.str();
		}
		
		if(flVirtualComponentList.find(moduleName) != flVirtualComponentList.end())
			return constraintString.str();
		
		string instanceName = flInstanceNames[moduleName];
		
		if(target->getVendor() == "Xilinx"){
			//create the constraint
			constraintString << "INST \"" << instanceName << "\" AREA_GROUP=\"pblock_" << instanceName << "\";" << endl;
			//add constraints for function generators and registers
			constraintString << "AREA_GROUP \"pblock_" << instanceName 
				<< "\" RANGE=SLICE_X" << (flComponentCordReal[moduleName]).x << "Y" << (flComponentCordReal[moduleName]).y
				<< ":SLICE_X" << (flComponentCordReal[moduleName]).x + (flComponentDimension[moduleName]).x
				<< "Y" << (flComponentCordReal[moduleName]).y + (flComponentDimension[moduleName]).y << ";" << endl;
			constraintString << "AREA_GROUP \"pblock_" << instanceName << "\" GROUP=OPEN;" << endl;
			constraintString << "AREA_GROUP \"pblock_" << instanceName << "\" PLACE=OPEN;" << endl;
			//add constraints for DSPs
			if((subComponents[moduleName])->reHelper->estimatedCountMultiplier != 0){
				vector<int> dspPositions;
				int dspInColumn = (flComponentDimension[moduleName]).y/(target->dspHeightInLUT);
				
				for(unsigned int i=0; i<target->multiplierPosition.size(); i++){
					int currentDSPColumn = target->multiplierPosition[i];
					
					if(currentDSPColumn == (flComponentCordReal[moduleName]).x-1)
						dspPositions.push_back(currentDSPColumn);
					else if((currentDSPColumn>=(flComponentCordReal[moduleName]).x) && (currentDSPColumn<=(flComponentCordReal[moduleName]).x+(flComponentDimension[moduleName]).x))
						dspPositions.push_back(currentDSPColumn);
				}
				
				for(unsigned int i=0; i<dspPositions.size(); i++){
					int currentDSPColumn = dspPositions[i];
					
					constraintString << "AREA_GROUP \"pblock_" << flInstanceNames[moduleName] 
						<< "\" RANGE=DSP48_X" << currentDSPColumn << "Y" << (flComponentCordReal[moduleName]).y/target->dspHeightInLUT
						<< ":DSP48_X" << currentDSPColumn
						<< "Y" << (flComponentCordReal[moduleName]).y/target->dspHeightInLUT + dspInColumn << ";" << endl;
				}
			}
			//add constraints for RAMs
			if((subComponents[moduleName])->reHelper->estimatedCountMemory != 0){
				vector<int> ramPositions;
				int ramInColumn = (flComponentDimension[moduleName]).y/target->dspHeightInLUT;
				
				for(unsigned int i=0; i<target->multiplierPosition.size(); i++){
					int currentRAMColumn = target->multiplierPosition[i];
					
					if(currentRAMColumn == (flComponentCordReal[moduleName]).x-1)
						ramPositions.push_back(currentRAMColumn);
					else if((currentRAMColumn>=(flComponentCordReal[moduleName]).x) && (currentRAMColumn<=(flComponentCordReal[moduleName]).x+(flComponentDimension[moduleName]).x))
						ramPositions.push_back(currentRAMColumn);
				}
				
				for(unsigned int i=0; i<ramPositions.size(); i++){
					int currentRAMColumn = ramPositions[i];
					
					if(target->getID() == "Virtex4" || target->getID() == "Spartan3"){
					constraintString << "AREA_GROUP \"pblock_" << flInstanceNames[moduleName] 
						<< "\" RANGE=RAMB16_X" << currentRAMColumn << "Y" << (flComponentCordReal[moduleName]).y/target->ramHeightInLUT
						<< ":RAMB16_X" << currentRAMColumn
						<< "Y" << (flComponentCordReal[moduleName]).y/target->ramHeightInLUT + ramInColumn << ";" << endl;
					}else{
					constraintString << "AREA_GROUP \"pblock_" << flInstanceNames[moduleName] 
						<< "\" RANGE=RAMB36_X" << currentRAMColumn << "Y" << (flComponentCordReal[moduleName]).y/target->ramHeightInLUT
						<< ":RAMB36_X" << currentRAMColumn
						<< "Y" << (flComponentCordReal[moduleName]).y/target->ramHeightInLUT + ramInColumn << ";" << endl;	
					}
				}
			}
			//end of constraint
			constraintString << endl;
		}else if(target->getVendor() == "Altera"){
			//add constraints for function generators
			
			//add constraints for registers
			
			//add constraints for DSPs
			
			//add constraints for RAMs
			
		}
		
		return constraintString.str();
	}
	
	std::string FloorplanningHelper::createFloorplan(){
		ostringstream result;
		
		cerr << "=========================================================================" << endl;
		cerr << "*                          Floorplanning                                *" << endl;
		cerr << "=========================================================================" << endl;
		cerr << "Starting the creation of the floorplan for operator " << parentOp->getName() << endl;
		cerr << "***Triggered creation of the virtual arrangement of the sub-components" << endl;
		result << createVirtualGrid();
		cerr << "***Triggered processing of placement constraints" << endl;
		result << processPlacementConstraints();
		cerr << "***Triggered processing of connectivity constraints" << endl;
		result << processConnectivityConstraints();
		cerr << "***Triggered creation of the actual arrangement of the sub-components" << endl;
		result << createPlacementGrid();
		cerr << "***Triggered creation of the constraints file" << endl;
		result << createConstraintsFile();
		cerr << "Finished creating the floorplan for operator " << parentOp->getName() << endl;
		cerr << "=========================================================================" << endl;
		
		return result.str();
	}
}
