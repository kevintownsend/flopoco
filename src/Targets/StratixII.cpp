/*
 * A model of Stratix II FPGA 
 *
 * Author : Florent de Dinechin
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
#include "StratixII.hpp"
#include <iostream>
#include <sstream>
#include "../utils.hpp"

double StratixII::adder_delay(int size) {
  return _lut_delay  +  size * _fastcarry_delay; 
};

double StratixII::carry_propagate_delay() {
  return  _fastcarry_delay; 
};

double StratixII::local_wire_delay(){
  return _lut2lut_delay;
};

double StratixII::distant_wire_delay(int n){
  return n*_elem_wire_delay;
};

double StratixII::lut_delay(){
  return _lut_delay;
};


bool StratixII::suggest_submult_size(int &x, int &y, int wInX, int wInY){
int i;

	if (get_use_hard_multipliers()){
		if ((wInX<=17) && (wInY<=17))	{
			x = max(wInX, wInY);
			y = x;
			 if (frequency()>600000000)
			 	return false;
			 else
			 	return true;
		}else{
			int f1=(wInX % 17 ==0)?0:1;
			int f2=(wInY % 17 ==0)?0:1;
			int k=wInX/17+wInY/17 + f1+ f2;
			x = 17;	y = 17;
					
			if (k<=4)
				if (frequency()<=400000000)
					return true;
				else 
					return false; 
			else{
				double freq;
				freq = 11.2 + 1560/k;
				if (frequency()<=freq*1000000)			
					return true;
				else 
					return false;
			} 
			
			
		}
	}else{
		int f1=(wInX % 17 ==0)?0:1;
		int f2=(wInY % 17 ==0)?0:1;
		int k=wInX/17+wInY/17 + f1+ f2;
		double freq;
		
		if ((max(wInX,wInY)<=4)&&(max(wInX,wInY)>=2))
		{
			freq = 669-13* (max(wInX,wInY));
			
			x=wInX;
			y=wInY;
			if (frequency()<=freq*1000000)			
				return true;
			else 
				return false;
		} else if ((max(wInX,wInY)<=15)&&(max(wInX,wInY)>=5)){
			freq = 411-9*(max(wInX,wInY));
			if (frequency()<=freq*1000000){			
				x=wInX;
				y=wInY;			
				return true;
			}
			else{
				freq = 121.3+549/2+988/4;
				if (frequency()>freq*1000000){
					x=2;
					y=2;
					return false;
				}else{
					int i=2;
					while ( ((121.3+549/i+988/(i*i))*1000000) > frequency())
						i++;
					
					x=i-1;
					y=i-1;
					return true;
				}
			}
		} else if (max(wInX,wInY)>15){
			freq = 80.2+34037/(max(wInX,wInY)*max(wInX,wInY));
			if (frequency()<=freq*1000000){
				x=wInX;
				y=wInY;			
				return true;
			}
			else{
				freq = 121.3+549/2+988/4;
				if (frequency()>freq*1000000){
					x=2;
					y=2;
					return false;
				}else{
					int i=2;
					while ( ((121.3+549/i+988/(i*i))*1000000) > frequency())
						i++;
					
					x=i-1;
					y=i-1;
					return true;
				}
			}	
		}
	
		
	}
};	 
	 
	 
bool StratixII::suggest_subadd_size(int &x, int wIn){

	int chunk_size = (int)floor( (1./frequency() - lut_delay()) / carry_propagate_delay()); // 1 if no need for pipeline
	
	x = chunk_size;		
	
	if (x>0) return true;
	else {
		x=1;		
		return false;
	} 
};
  
