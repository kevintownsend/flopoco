/*
 * An integer multiplier for FloPoCo
 *
 * Authors : Bogdan Pasca and Florent de Dinechin
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
#include "IntMultiplier.hpp"

using namespace std;


IntMultiplier:: IntMultiplier(Target* target, int wInX, int wInY) :
   Operator(target), wInX(wInX), wInY(wInY), wOut(wInX + wInY), reverse(false){
 
   int depth, i, j;
   ostringstream name, nameH, nameL;

   //gets information if the multiplier is or not to be pipelined
   if(target->is_pipelined())
      set_sequential();
   else
      set_combinatorial();  

   //============================ Name Setup procedure==========================
   //  The name has the format: IntMultiplier_wInX_wInY
   //  wInX = width of the X input
   //  wInY = width of the Y input
   name.str("");;
   name <<"IntMultiplier_"<<wInX<<"_"<<wInY;
   unique_name = name.str(); //attribute assignment
   //===========================================================================   

   //============================ Set up the IO signals=========================
   //  X and Y have wInX and wInY respectively 
   //  R has wOut bits where wOut = (wInX + WInY) bits
   add_input ("X", wInX);
   add_input ("Y", wInY);
   add_output("R", wOut);
   //===========================================================================
  
   if(is_sequential()) 
   {
      multiplier_width = target->mult_x_inputs()-1;

      // decide in how many parts should the numbers be split depending on the
      // width of the target multipliers
      partsX = int(   ceil( double(wInX) / double(multiplier_width)));
      partsY = int(   ceil( double(wInY) / double(multiplier_width)));

      //if X's width is smaller then Y, make Y=X and X=Y    
      if (partsX < partsY){
         i=partsX;  
         partsX=partsY; 
         partsY=i;
         reverse=true;
         // number of zeros that the numbers must be padded with, so that their 
         // length becomes a multiple of the target multiplier width
         number_of_zerosX = partsY * multiplier_width - wInX;
         number_of_zerosY = partsX * multiplier_width - wInY;
      }
      else{
         number_of_zerosX = partsX * multiplier_width - wInX;
         number_of_zerosY = partsY * multiplier_width - wInY;
      }    

      // signals that will contain the inputs X and Y but padded with 0 so that 
      // their width becomes a multiple of the target multiplier width
      add_signal("i_X", partsX * multiplier_width);	
      add_signal("i_Y", partsY * multiplier_width);

      //declare the signals needed to split X 
      for (i=1;i<=partsX;i++){
      	name.str("");
      	name <<"X_"<<i;
      	add_signal(name.str(), multiplier_width); 
      } 	

      //declare the signals needed to split Y 
      for (i=1;i<=partsY;i++){
      	name.str("");
      	name <<"Y_"<<i;
      	add_signal(name.str(), multiplier_width); 
      }	

      //declare the registers needed to store the partial products
      for (i=1;i<=partsY;i++)
         for (j=1;j<=partsX;j++){
      		name.str("");;
      		name <<"Y_"<<i<<"_X_"<<j;
      		add_registered_signal_with_synch_reset(name.str(), 2 * multiplier_width);	
      	}	 

      //declare the high/low parts decomposition signals
      for (i=1;i<=partsY;i++)
      	for (j=1;j<=partsX;j++)	{
            nameH.str("");; nameL.str("");;
            nameH <<"Y_"<<i<<"_X_"<<j<<"_H";
            nameL <<"Y_"<<i<<"_X_"<<j<<"_L";    	
            add_signal(nameH.str(), multiplier_width);	
            add_signal(nameL.str(), multiplier_width);	
      	}	 

      //signals to regroup all the high parts and all the low parts
      for (i=1;i<=partsY;i++)  {	
      	nameH.str("");; nameL.str("");;
      	nameH <<"High_"<<i;
         nameL <<"Low_"<<i;    	
         add_signal(nameH.str(), partsX * multiplier_width);	
      	add_signal(nameL.str(), partsX * multiplier_width);	
      }	  
      
      if (partsY>1) {	
      //declare the signals which form the two trees (high and low)
         for (i=1;i<=partsY;i++)	{
             name.str("");;
             name<<"Buffer_H_"<<i;
             add_registered_signal_with_synch_reset(name.str(), partsX * multiplier_width);	
         }
         
         //for the low and high tree
         for (j=0; j<partsY;j++)	
            for (i=1;i<=partsY-j;i++)	{	
               nameH.str("");; nameL.str("");;
               nameL <<"L_Level_"<<j<<"_Reg_"<<i;	
               nameH <<"H_Level_"<<j<<"_Reg_"<<i;	
               add_registered_signal_with_synch_reset(nameL.str(), partsX * multiplier_width);	
               add_registered_signal_with_synch_reset(nameH.str(), partsX * multiplier_width);	
            }	

         //one more registers for the Low part for synchronization 
         add_registered_signal_with_synch_reset("Low1", partsX * multiplier_width);	
      		
         //one more signal for the high part
         add_signal("High1", partsX * multiplier_width);	

         //partsY! registers to hold and compute the low part of the result in the 
         //pipeline levels
      	for (j=1; j<=partsY;j++)	
      		for (i=1;i<=j;i++){	
      	      name.str("");;
      	      name<<"PartialBits_Level_"<<j<<"_Reg_"<<i;
      	      add_registered_signal_with_synch_reset(name.str(), multiplier_width + 1 );	
          }

         for (j=1; j<=partsX;j++)	
            for (i=1;i<=partsX;i++){	
               name.str("");;
               name<<"Last_Addition_Level_"<<j<<"_Reg_"<<i;
               add_registered_signal_with_synch_reset(name.str(), multiplier_width + 1 );	
         }       

         for (j=1; j<=partsX;j++){	    
            name.str("");;
            name<<"PartialBits_Reg_"<<j;
            add_registered_signal_with_synch_reset(name.str(), partsY*multiplier_width);	
         }       
      	
         add_signal("temp_result",  partsX * multiplier_width );
         add_signal("partial_bits", partsY * multiplier_width );
         add_signal("full_result",  partsX * multiplier_width + partsY * multiplier_width);

      }//end if partsY>1
      else
    	   add_signal("temp_result", (partsX + 1) * multiplier_width );
     
      
      // Set up the pipeline
   	//========================================================================
   	if (partsY==1)
	      depth=1;
	   else
	      depth= 2 + partsY + partsX;
    
      set_pipeline_depth(depth);
  
   }//end if sequential
}


IntMultiplier::~IntMultiplier() {
}


void IntMultiplier::output_vhdl(std::ostream& o, std::string name) {
   Licence(o,"Bogdan Pasca and Florent de Dinechin (2008)");
   Operator::StdLibs(o);
   output_vhdl_entity(o);

   int i, j, k;
   ostringstream zeros, zerosX, zerosY, names, nameH, nameL, concatH, concatL,
                 first_summand, second_summand, the_bits;

   o << "architecture arch of " << name  << " is" << endl;
   output_vhdl_signal_declarations(o);
   o << "begin" << endl;
 
   if (is_sequential()) { 
	   output_vhdl_registers(o); 
        
   	// iX and iY have the length multiple of target->mult_x_inputs
      // generate the strings of zeros needed to pad X and Y
  	   zerosX.str(""); zerosY.str("");
   	for (i=0;i<number_of_zerosX;i++)
   		zerosX<<"0";  
	   for (i=0;i<number_of_zerosY;i++)
  		   zerosY<<"0";  
  	
      //pad X and Y
      o <<endl;
      if (reverse==false){
         o<<tab<< "i_X <= \""<< zerosX.str()<<"\" & X;"<<endl;
  	      o<<tab<< "i_Y <= \""<< zerosY.str()<<"\" & Y;"<<endl<<endl<<endl;
  	   }
  	   else{//reverse the inputs
         o<<tab<< "i_Y <= \""<< zerosX.str()<<"\" & X;"<<endl;
  	      o<<tab<< "i_X <= \""<< zerosY.str()<<"\" & Y;"<<endl<<endl<<endl;
      }
     
      //split X 
      for (i=1;i<=partsX;i++)	{
         names.str("");;
    	   names <<"X_"<<i;
 		   o<<tab<< names.str() << " <= i_X("<< i * multiplier_width -1 <<" downto "<< (i-1) * multiplier_width<<");"<<endl;		
	   }	
	  
      o<<endl;
      //split Y 
      for (i=1;i<=partsY;i++) {
         names.str("");;
         names <<"Y_"<<i;
 		   o<<tab<< names.str() << " <= i_Y("<< i * multiplier_width -1 <<" downto "<< (i-1) * multiplier_width<<");"<<endl;			
      }	

	   o<<endl;
      //make all multiplications concurrently
     	for (i=1;i<=partsY;i++)
     		for (j=1;j<=partsX;j++){  
     			names.str("");;
     			names <<"Y_"<<i<<"_X_"<<j;
     			o<<tab<< names.str() << " <= " << "Y_"<< i <<" * " << "X_"<<j << ";"<<endl;	
     		}	 

      o<<endl;
      //decompose the results into high and low part
      for (i=1;i<=partsY;i++)
      	for (j=1;j<=partsX;j++){
      		names.str("");; nameH.str("");; nameL.str("");;
      		names <<"Y_"<<i<<"_X_"<<j<<"_d";
   	      nameH <<"Y_"<<i<<"_X_"<<j<<"_H";
   	      nameL <<"Y_"<<i<<"_X_"<<j<<"_L";
      		o<<tab<< nameH.str() << " <= " << names.str() <<"("<< 2 * multiplier_width - 1 <<" downto " << multiplier_width <<");"<<endl;
   	      o<<tab<< nameL.str() << " <= " << names.str() <<"("<< multiplier_width - 1 <<" downto " << 0 <<");"<<endl;			
      	}	 

	   o<<endl<<endl;
	   //Regroup all the high parts together and all the low parts together
	   for (i=1;i<=partsY;i++)	{	
    		nameH.str("");; nameL.str("");; concatH.str("");; concatL.str("");;
    		nameH <<"High_"<<i;
		   nameL <<"Low_"<<i;    	
			//create the concatenation strings for low and high parts
		   for (j=partsX;j>=1;j--){
			   if (j>1){   
				   concatL<<"Y_"<<i<<"_X_"<<j<<"_L & ";
			      concatH<<"Y_"<<i<<"_X_"<<j<<"_H & ";	
			   }
			   else{	
				  concatL<<"Y_"<<i<<"_X_"<<j<<"_L";
				  concatH<<"Y_"<<i<<"_X_"<<j<<"_H";
			   }			
		   }	
		   
		   o<<tab<< nameL.str() << "  <= " << concatL.str() << ";"<<endl;	
		   o<<tab<< nameH.str() << " <= " << concatH.str() << ";"<<endl;
	   }

	   o<<endl<<endl;
      //pass to the adder tree
      if (partsY==1) {
		   //add zeros to the left of Low_  & add zeros right of the High_
	      zeros.str("");;
		   for (i=0;i<multiplier_width;i++)
		  	   zeros<<"0";
      
         o<<tab<<"temp_result <= (\""<< zeros.str()<< "\" & Low_1) + (High_1 & \"" << zeros.str()<<"\");"<<endl;
		   //output the result
 		   o<<tab<<"R <= temp_result("<<(partsX+ partsY)*multiplier_width-1 -number_of_zerosX - number_of_zerosY<<" downto 0);"<<endl;
      }
      else
      {
  	      //Connect the low part tohether 
         for (i=0;i<partsY-1;i++) {
      	   for (j=1;j<=partsY-i;j++) {
      	   	if (j==1) {
      	   		//concatenate zeros in front of the first register, add it with the
      	   		//second register and place the result in the next level's first register
      	   		zeros.str("");;
      	   		for (k=0;k<multiplier_width;k++)
      	   			zeros<<"0";
      	
      	   		first_summand.str("");; second_summand.str("");;
      	   		first_summand<<"(\""<< zeros.str()<<"\" & L_Level_"<<i<<"_Reg_"<<j<<"_d("<<partsX*multiplier_width-1<<" downto "<<multiplier_width<<"))"; 
                  second_summand<<"L_Level_"<<i<<"_Reg_"<<j+1<<"_d";				
                  o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j<<"<="<< first_summand.str() << " + " << second_summand.str()<<";"<< endl;
  		         }
  		         else
  			         if (j==2) {}	 //do nothing (the second register on the level is always processed with the first
  			            else //for the rest of the registers just propagate their contents to the next level
  			               o<<tab<<"L_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<<"L_Level_"<<i<<"_Reg_"<<j<<"_d"<<";"<<endl;
  	         }
  	      }

         //Connect the HIGH part tohether 
         for (i=0;i<partsY-1;i++){
      	   for (j=1;j<=partsY-i;j++){
      	   	if (j==1){
      	   		//concatenate zeros in front of the first register, add it with the
      	   		//second register and place the result in the next level's first register
      	   		zeros.str("");;
      	   		for (k=0;k<multiplier_width;k++)
      	   			zeros<<"0";
      	     
      	     		first_summand.str("");; second_summand.str("");;
                  first_summand<<"(\""<< zeros.str()<<"\" & H_Level_"<<i<<"_Reg_"<<j<<"_d("<<partsX*multiplier_width-1<<" downto "<<multiplier_width<<"))"; 
                  second_summand<<"H_Level_"<<i<<"_Reg_"<<j+1<<"_d";				
                  o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j<<"<="<< first_summand.str() << " + " << second_summand.str()<<";"<< endl;
      		   }
      		   else
      			   if (j==2) {} //do nothing (the second register on the level is always processed with the first
      			      else //for the rest of the registers just propagate their contents to the next level
      			         o<<tab<<"H_Level_"<<i+1<<"_Reg_"<<j-1<<"<="<<"H_Level_"<<i<<"_Reg_"<<j<<"_d;"<<endl;
            }
         }

         o<<endl<<endl;
         //connect the low part of the multiplication results to the low tree
         for(i=1;i<=partsY;i++)
  	         o<<tab<<"L_Level_0_Reg_"<<i<<" <= Low_"<<i<<";"<<endl;

         o<<endl<<endl;
         //the high part of the multiplication needs to be delayed one clock cycle before
         //it will be inserted into the adder tree
         for(i=1;i<=partsY;i++)	
  	         o<<tab<<"Buffer_H_"<< i << " <= High_" <<i<<";"<<endl;

         o<<endl<<endl;
         //the buffers need to be connected to the high tree
         for(i=1;i<=partsY;i++)
  	         o<<tab<<"H_Level_0_Reg_"<<i<<" <= Buffer_H_"<<i<<"_d;"<<endl;

         o<<endl<<endl;
         o<<tab<<"High1 <= H_Level_"<<partsY-1<<"_Reg_1_d"<<";"<<endl;

         zeros.str("");
         for (i=0;i<multiplier_width;i++)
            zeros<<"0";

         o<<endl<<endl;
         //the zeros + [the high bits of the last register of the low part] 
         o<<tab<<"Low1 <= \""<<zeros.str()<<"\" & L_Level_"<<partsY-1<<"_Reg_1_d("<<partsX * multiplier_width-1<<" downto "<<multiplier_width<<")"<< ";"<<endl;

         o<<endl<<endl;
         //connect the partsY! registers which compute the low part of the result 
    	   for (j=1; j<=partsY;j++)	
    		   for (i=1;i<=j;i++){	
  		         if ((j==1) and (i==1) )
  			         o<<tab<<"PartialBits_Level_1_Reg_1 <= \"0\" & L_Level_0_Reg_1_d("<< multiplier_width -1 <<" downto 0);"<< endl;
  		         else{
  			         //just propagate the registers from 1 to j-1
  			         if (i<j)									
  			            o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i <<"<= "<<"PartialBits_Level_"<<j-1<<"_Reg_"<< i <<"_d;"<<endl;
  			         else //if i=j compute the jth register
  			            o<<tab<<"PartialBits_Level_"<<j<<"_Reg_"<<i <<"<= "<<"(\"0\" & L_Level_"<<j-1<<"_Reg_1_d("<<multiplier_width-1<<" downto 0)) + "<<"(\"0\" & H_Level_"<<j-2<<"_Reg_1_d("<<multiplier_width-1<<" downto 0)) + "<< "PartialBits_Level_"<<j-1<<"_Reg_"<<i-1<<"_d("<< multiplier_width  <<") ;"<<endl; 
  		         }	
  	         }

         o<<endl<<endl;
     
         //connect first parts of these registers
         for (j=1; j<=partsX;j++){
            if (j==1)    
               o<<tab<<"Last_Addition_Level_1_Reg_"<<j<<" <= (\"0\" & Low1_d("<<j*multiplier_width -1<<" downto "<<(j-1)*multiplier_width<<")) + ( \"0\" & High1("<<j*multiplier_width -1<<" downto "<<(j-1)*multiplier_width<<")) + PartialBits_Level_"<< partsY <<"_Reg_"<< partsY<<"_d("<< multiplier_width <<");"<<endl;        
            else
               o<<tab<<"Last_Addition_Level_1_Reg_"<<j<<" <= (\"0\" & Low1_d("<<j*multiplier_width -1<<" downto "<<(j-1)*multiplier_width<<")) + (\"0\" & High1("<<j*multiplier_width -1<<" downto "<<(j-1)*multiplier_width<<"));"<<endl  ;     
         }	
      
         //put tohether the pipeline
         for (j=1; j<=partsX-1;j++)	
    		   for (i=1;i<=partsX;i++)
			      if ((i<=j)||(i>j+1))
				      o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d;"<<endl;
				   else
			         o<<tab<<"Last_Addition_Level_"<<j+1<<"_Reg_"<<i<<" <= Last_Addition_Level_"<<j<<"_Reg_"<<i<<"_d + Last_Addition_Level_"<<j<<"_Reg_"<<i-1<<"_d("<<multiplier_width<<") ;"<<endl; 
			    
         //put the result back together
         zeros.str("");
         for (i=partsX;i>0;i--) 
            if (i!=1)
               zeros<<"Last_Addition_Level_"<<partsX<<"_Reg_"<<i<<"_d("<<multiplier_width -1<<" downto 0) &" ;
            else
               zeros<<"Last_Addition_Level_"<<partsX<<"_Reg_"<<i<<"_d("<<multiplier_width -1<<" downto 0);" ;
               
         o<<tab<<"temp_result <= "<<zeros.str()<<endl;

         //make the string that will gather all partial bits from the last level
         the_bits.str("");;	
         for (i=partsY;i>0;i--)
  	         if (i!=1) 	
  	            the_bits << "PartialBits_Level_"<< partsY <<"_Reg_"<<i<<"_d("<< multiplier_width-1 << " downto 0" <<")"<<" & ";	
  	         else
  	            the_bits << "PartialBits_Level_"<< partsY <<"_Reg_"<<i<<"_d("<< multiplier_width-1 << " downto 0" <<")";		

   
         //setup the pipeline part that will carry the low part of the final result
         o<<tab<<"PartialBits_Reg_1 <= "<<the_bits.str()<<";"<<endl;
         for (i=2;i<=partsX;i++)
            o<<tab<<"PartialBits_Reg_"<<i<<" <= PartialBits_Reg_"<<i-1<<"_d;"<<endl;
      
      
         o<<tab<<"partial_bits <= PartialBits_Reg_"<<partsX<<"_d;"<<endl;
         o<<tab<<"full_result <= temp_result & partial_bits;"<<endl; 
         o<<tab<<"R <= full_result("<< (partsX+ partsY)*multiplier_width-1 -number_of_zerosX - number_of_zerosY<<" downto 0);"<<endl;
        
         o<<endl<<endl;
      }	
   }//endl if sequential
   else //the combinational version
	  o << tab<<"R <= X * Y ;" << endl;
  
   
   
   o << "end architecture;" << endl << endl;
}




