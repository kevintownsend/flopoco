/*
  TaMaDiCore evaluation of worst cases for rounding functions without 
  overheating the planet

  This file is part of the FloPoCo project 
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author :   Bogdan Pasca

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  CeCILL license, 2008-2010.

 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <cstdlib>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"
#include "../IntAddition/IntAdderSpecific.hpp"

#include "TaMaDiCore.hpp"
//#define IAS

using namespace std;

namespace flopoco{

	extern vector<Operator*> oplist;

	TaMaDiCore::TaMaDiCore(Target* target, int wp, int d, int iterations, int wIntervalID, int version, int pathWidth):
	Operator(target), wp(wp),d(d), iterations(iterations), wIntervalID(wIntervalID) 
	{
		srcFileName="TaMaDiCore";
		ostringstream name;

		name <<"TaMaDiCore_wp"<<wp<<"_interations"<<iterations<<"_degree"<<d<<"_uid"<<getNewUId();
		setName(name.str()); 

		setCopyrightString("Bogdan Pasca (2011)");		

		/* Set up the I/O signals of of the entity */
		addInput   ("SerialPin", wp);
		addInput   ("Initialize");
		addInput   ("IntervalID",wIntervalID);
		addInput   ("CE");

		int counterWidth = intlog2(iterations);
		addOutput ("potentialUlpNumber", counterWidth);
		addOutput ("potentialInterval", wIntervalID);
		addOutput ("potentialOutput");

		addOutput ("finished"); 
		 
		/* computing core description */

		/* describe the counter */
		int countSize = intlog2(d+1);
		
		vhdl << tab << declare( "counter_ce") << " <= (Initialize and counter_is_zero) or ((not counter_is_zero) and (not finished_signal));"<<endl;
		vhdl << tab << declare( "ulp_counter_reset") << "<= Initialize;" <<endl;

		declare("c",wp);
		vhdl << tab << " process(clk, rst, c_ce, CE) " << endl;
		vhdl << tab << "      begin" << endl;
		vhdl << tab << "         if rst = '1' then" << endl;
		vhdl << tab << "               c <=  (others => '0');" << endl;
		vhdl << tab << "         elsif clk'event and clk = '1' then" << endl;
		vhdl << tab << "            if c_ce='1' and CE='1' then " << endl;
		vhdl << tab << "               c <=  SerialPin;" << endl;
		vhdl << tab << "         	end if;" << endl;
		vhdl << tab << "         end if;" << endl;
		vhdl << tab << "      end process;" << endl;

		vhdl << tab << " process(clk, rst, counter_ce, CE) " << endl;
		vhdl << tab << "      begin" << endl;
		vhdl << tab << "         if rst = '1' then" << endl;
		vhdl << tab << "               counter_reg_d1 <=  (others => '0');" << endl;
		vhdl << tab << "         elsif clk'event and clk = '1' then" << endl;
		vhdl << tab << "            if counter_ce='1' and CE='1' then " << endl;
		vhdl << tab << "               counter_reg_d1 <=  counter_reg;" << endl;
		vhdl << tab << "               finished_signal_d1 <=  finished_signal;" << endl;
		vhdl << tab << "         	end if;" << endl;
		vhdl << tab << "         end if;" << endl;
		vhdl << tab << "      end process;" << endl;

		vhdl << tab << " process(clk, rst, ulp_counter_reset,CE) " << endl;
		vhdl << tab << "      begin" << endl;
		vhdl << tab << "         if rst = '1' or ulp_counter_reset='1' then" << endl;
		vhdl << tab << "               ulp_counter_d1 <=  (others => '0');" << endl;
		vhdl << tab << "         elsif clk'event and clk = '1' then" << endl;
		vhdl << tab << "            if CE='1' then"<<endl;
		vhdl << tab << "               ulp_counter_d1 <=  ulp_counter;" << endl;
		for (int i=1;i<=d;i++)
		vhdl << tab << "               ulp_counter_d"<<i+1<<" <=  ulp_counter_d"<<i<<";" << endl;
		vhdl << tab << "            end if; "<<endl;
		vhdl << tab << "         end if;" << endl;
		vhdl << tab << "      end process;" << endl;

		vhdl << tab << " process(clk, rst, Initialize) " << endl;
		vhdl << tab << "      begin" << endl;
		vhdl << tab << "         if rst = '1' then" << endl;
		vhdl << tab << "               iid <=  (others => '0');" << endl;
		vhdl << tab << "         elsif clk'event and clk = '1' then" << endl;
		vhdl << tab << "            if Initialize='1' then"<<endl;
		vhdl << tab << "               iid <=  IntervalID;" << endl;
		vhdl << tab << "            end if; "<<endl;
		vhdl << tab << "         end if;" << endl;
		vhdl << tab << "      end process;" << endl;
		
		//TODO Add signal declarations with reset and clock enable options
		setCycle(1,false);
		for (int i=1;i<=d;i++)
			declare(join("ulp_counter_d",i+1),counterWidth);
		setCycle(0,false);

		setCycle(1,false);
		declare( "counter_reg_d1", countSize, true);
		declare( "ulp_counter_d1", counterWidth, true);
		declare("iid",wIntervalID);
		setCycle(0,false);
		
		vhdl << tab << declare( "counter_reg", countSize, true /*is bus*/) << " <= 	counter_reg_d1 + '1'; "<< endl;
		vhdl << tab << declare( "ulp_counter", counterWidth) << " <= ulp_counter_d1 + '1';" << endl;
		setCycle(1,false);
		vhdl << tab << declare("counter_is_zero") << "<= '1' when counter_reg_d1 = "<<zg(countSize,0) << " else '0';" << endl;
		setCycle(0,false);

		/*signal Decoder */
		vhdl << tab << declare("c_ce") << "<= '1' when (counter_reg_d1=CONV_STD_LOGIC_VECTOR("<<1<<","<<countSize<<")) else '0';"<<endl;
		for (int i=0; i<d; i++)
			vhdl << tab << declare( join("select",i) ) << "<= '1' when (counter_reg_d1=CONV_STD_LOGIC_VECTOR("<<i+2<<","<<countSize<<")) else '0';"<<endl;

		/* register interval id */
//		vhdl << tab << declare("c",wp) << " <= (others =>'1');"<<endl;

#ifdef IAS
		IntAdderSpecific *ia = new IntAdderSpecific(target, wp);
		oplist.push_back(ia);
#endif
		setCycle(1,false);
		for(int i=0; i<d; i++)
			vhdl << tab << declare(join("regSelect",i)) << " <= " << join("select",i) << ";" << endl;
		setCycle(0,false);

		for(int i=0; i<d; i++){
			setCycle(0);
			vhdl << tab << declare( join("mux",i), wp) << " <= SerialPin when " << join("regSelect",i) << " = '1' else ";
			setCycle(1,false);
			vhdl << join("adder",i) << ";" << endl;
			setCycle(0,false);


#ifndef IAS			
			if (i==0){
				vhdl << tab << declare( join("adder",i), wp ) << " <= " << join("mux",i) << " + c;" << endl;  	
			}else{
				vhdl << tab << declare( join("adder",i), wp ) << " <= " << join("mux",i) << " + "; 
				setCycle(1,false);
				vhdl << join("adder",i-1) <<";"<< endl;  	
				setCycle(0,false);
			}							
#else
		setCycle(1,false);
		inPortMapCst(ia, "X", join("mux",i));
		setCycle(0,false);

		if (i==0)
			inPortMap(ia, "Y", "c");
		else{
			setCycle(1,false);
			inPortMap(ia, "Y", join("adder",i-1));
			setCycle(0,false);
		}
		inPortMapCst(ia, "Cin", "'0'");
		outPortMap (ia, "R", join("adder",i));
		
		vhdl << tab << instance(ia, join("ia",i))<<endl;
#endif
		}
		
		setCycle(1);
		/* one comparator */
		vhdl << tab << "potentialUlpNumber <= ulp_counter_d1;" << endl; /*permanent*/
		setCycle(0);
		vhdl << tab << "potentialInterval <= iid;"<<endl;
		setCycle(1);
		declare("finished_signal_d1");

if (target->getVendor()=="Altera"){
		declare("potentialOutput_1");	
		vhdl << tab << "LPM_COMPARE_component"<<getNewUId()<<" : LPM_COMPARE"<<endl;
		vhdl << tab << "GENERIC MAP ("<<endl;
		vhdl << tab << "lpm_hint => \"ONE_INPUT_IS_CONSTANT=YES\","<<endl;
		vhdl << tab << "lpm_pipeline => 2,"<<endl;
		vhdl << tab << "lpm_representation => \"UNSIGNED\","<<endl;
		vhdl << tab << "lpm_type => \"LPM_COMPARE\","<<endl;
		vhdl << tab << "lpm_width => "<<wp<<endl;;
		vhdl << tab << "	)"<<endl;
		vhdl << tab << "PORT MAP ("<<endl;
		vhdl << tab << "clock => clk,"<<endl;
		vhdl << tab << "dataa => "<<join("adder",d-1)<<","<<endl;
		vhdl << tab << "datab => "<<"\"1"<<zg(wp-1,1) <<","<<endl;
		vhdl << tab << "AeB => potentialOutput_1"<<endl;
		vhdl << tab << "	);"<<endl;

		declare("potentialOutput_2");	
		vhdl << tab << "LPM_COMPARE_component"<<getNewUId()<<" : LPM_COMPARE"<<endl;
		vhdl << tab << "GENERIC MAP ("<<endl;
		vhdl << tab << "lpm_hint => \"ONE_INPUT_IS_CONSTANT=YES\","<<endl;
		vhdl << tab << "lpm_pipeline => 2,"<<endl;
		vhdl << tab << "lpm_representation => \"UNSIGNED\","<<endl;
		vhdl << tab << "lpm_type => \"LPM_COMPARE\","<<endl;
		vhdl << tab << "lpm_width => "<<wp<<endl;;
		vhdl << tab << "	)"<<endl;
		vhdl << tab << "PORT MAP ("<<endl;
		vhdl << tab << "clock => clk,"<<endl;
		vhdl << tab << "dataa => "<<join("adder",d-1)<<","<<endl;
		vhdl << tab << "datab => "<<"\"0"<<og(wp-1,1) <<","<<endl;
		vhdl << tab << "AeB => potentialOutput_2"<<endl;
		vhdl << tab << "	);"<<endl;

		vhdl << tab << "potentialOutput <= potentialOutput_1 or potentialOutput_2;"<<endl;
}else{
		vhdl << tab << "potentialOutput <= '1' when ("<<join("adder",d-1)<<"="<<"\"1"<<zg(wp-1,1)<<" or "<<join("adder",d-1)<<"="<<"\"0"<<og(wp-1,1) <<") else '0';"<<endl;
}

		vhdl << tab << declare("finished_signal") << " <= '1' when ulp_counter_d"<<d<<"=CONV_STD_LOGIC_VECTOR("<<iterations<<","<<counterWidth<<") else '0';"<<endl;
		vhdl << tab << "finished <= finished_signal_d1;"<<endl;
	}

	TaMaDiCore::~TaMaDiCore() {
	}
	
	void TaMaDiCore::outputVHDL(std::ostream& o, std::string name) {
		ostringstream signame;
		licence(o);
		pipelineInfo(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.std_logic_arith.all;" << endl;
		o << "use ieee.std_logic_unsigned.all;" << endl;

		o << "library work;" << endl;
		if (target_->getVendor() == "Xilinx"){
			o << "library UNISIM;"<<endl;
			o << "use UNISIM.VComponents.all;"<<endl;
		}else if(target_->getVendor() == "Altera"){
			o << "LIBRARY lpm;"<<endl;
			o << "USE lpm.all;"<<endl;
		}
		outputVHDLEntity(o);
		newArchitecture(o,name);
		if (target_->getVendor() == "Altera"){
			o << "	COMPONENT lpm_add_sub "<<endl;
			o << "	GENERIC ("<<endl;
			o << "		lpm_direction		: STRING;"<<endl;
			o << "		lpm_hint		: STRING;"<<endl;
			o << "		lpm_representation		: STRING;"<<endl;
			o << "		lpm_type		: STRING;"<<endl;
			o << "		lpm_width		: NATURAL"<<endl;
			o << "	);"<<endl;
			o << "	PORT ("<<endl;
			o << "			cin	: IN STD_LOGIC ;"<<endl;
			o << "			datab	: IN STD_LOGIC_VECTOR (lpm_width-1 DOWNTO 0);"<<endl;
			o << "			cout	: OUT STD_LOGIC ;"<<endl;
			o << "			dataa	: IN STD_LOGIC_VECTOR (lpm_width-1 DOWNTO 0);"<<endl;
			o << "			result	: OUT STD_LOGIC_VECTOR (lpm_width-1 DOWNTO 0)"<<endl;
			o << "	);"<<endl;
			o << "	END COMPONENT;"<<endl;			

			o << "COMPONENT lpm_compare"<<endl;
			o << "GENERIC ("<<endl;
			o << "	lpm_hint		: STRING;"<<endl;
			o << "	lpm_pipeline		: NATURAL;"<<endl;
			o << "	lpm_representation		: STRING;"<<endl;
			o << "	lpm_type		: STRING;"<<endl;
			o << "	lpm_width		: NATURAL"<<endl;
			o << ");"<<endl;
			o << "PORT ("<<endl;
			o << "		AeB	: OUT STD_LOGIC ;"<<endl;
			o << "		clock	: IN STD_LOGIC ;"<<endl;
			o << "		dataa	: IN STD_LOGIC_VECTOR (lpm_width-1 DOWNTO 0);"<<endl;
			o << "		datab	: IN STD_LOGIC_VECTOR (lpm_width-1 DOWNTO 0)"<<endl;
			o << ");"<<endl;
			o << "END COMPONENT;"<<endl;
		
		}

		o << buildVHDLComponentDeclarations();	
		o << buildVHDLSignalDeclarations();
		beginArchitecture(o);		
		o<<buildVHDLRegisters();
		o << vhdl.str();
		endArchitecture(o);
	}
	
	

}



