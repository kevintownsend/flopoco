/*
 * Floating Point Adder for FloPoCo
 *
 * Author :  Radu Tudoran, Bogdan Pasca
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
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>


#include <gmpxx.h>
#include "../../utils.hpp"
#include "../../Operator.hpp"

#include "CoilInductance.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <math.h>
#include <locale>

#include <stdio.h>
#include <mpfr.h>

using namespace std;
using std::ifstream;
extern vector<Operator*> oplist;

#define DEBUGVHDL 0


CoilInductance::CoilInductance(Target* target, int LSBI, int MSBI, int MaxMSBO,int LSBO, int MSBO, char *filepath) :
	Operator(target), MSBI(MSBI), LSBI(LSBI), MaxMSBO(MaxMSBO),LSBO(LSBO), MSBO(MSBO) ,filepath(filepath){
	
	if ((MSBI < LSBI)){
		cerr << 
			" CoilInductance: Input constraint LSBI <= MSBI not met."<<endl;
		exit (EXIT_FAILURE);
	}
	
	if ((MSBO < LSBO)){
		cerr << 
			" CoilInductance: Input constraint LSBO <= MSBO not met."<<endl;
		exit (EXIT_FAILURE);
	}
	
	wE=8;
	wF=23;
	
	ostringstream name; 
	name <<"CoilInductance_"<<abs(LSBI)<<"_"<<abs(MSBI)<<"_"<<abs(LSBO)<<"_"<<abs(MSBO);
	setName(name.str());
	
	setCopyrightString("Bogdan Pasca, Radu Tudoran (2009)");
	
	inputWidth= MSBI-LSBI;
	outputWidth= MSBO-LSBO;
	integratorWidth=3;
	inputWidthSegments=2*(inputWidth-1)+1;
	
	addrWidth =addressLength();
	
	double tempFreq=target->frequency();
	long internFreq = 200000000;
	
	internFreq = tempFreq;
	
	addOutput("O",outputWidth);
	//addFPOutput("O",8,23);
	
	//Counters for addressing the memories and for frequency division
	
	//~ vhdl<<endl;
	//~ vhdl<<tab<<"process (clk)"<<endl;
	//~ vhdl<<tab<<"variable count:std_logic_vector"<<range(integratorWidth,0)<<":=(others=>'0');"<<endl;
	//~ vhdl<<tab<<"begin"<<endl<<tab<<tab<<"if clk'event and clk = '1' then"<<endl;
	//~ vhdl<<tab<<tab<<tab<<"if count < "<<intpow2(integratorWidth)<<" then "<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<"count:=count+'1';"<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<declare("out_clk11",1)<<"<= '0';"<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<declare("out_rst",1)<<"<= '0'; "<<endl;
	//~ vhdl<<tab<<tab<<tab<<"elsif count = "<<intpow2(integratorWidth)<<" then "<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<"count:=count+'1';"<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<use("out_clk11")<<"<= '1'; "<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<use("out_rst")<<"<= '0'; "<<endl;
	//~ vhdl<<tab<<tab<<tab<<"else "<<endl<<tab<<tab<<tab<<tab<<"count:= CONV_STD_LOGIC_VECTOR(0,"<<integratorWidth<<");"<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<use("out_clk11")<<"<= '0'; "<<endl;
	//~ vhdl<<tab<<tab<<tab<<tab<<use("out_rst")<<"<= '1'; "<<endl;
	//~ vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	//~ vhdl<<tab<<tab<<"end if;"<<endl;
	//~ vhdl<<tab<<declare("out_clk2",1)<<"<="<<"clk;"<<endl;
	//~ vhdl<<tab<<declare("signal_tp",integratorWidth+1)<<" <= '0' & "<<" count"<<range(integratorWidth-1,0)<<";"<<endl;
	//~ vhdl<<tab<<"end process;"<<endl;
	//~ vhdl<<endl;
	
	vhdl<<tab<<declare("out_clk1",1)<<" <= clk;"<<endl;
	vhdl<<tab<<declare("out_rst",1)<<" <= '0';"<<endl;
	//vhdl<<tab<<declare("out_clk1",1)<<" <= "<<use("out_clk11")<<" and "<<" clk ;"<<endl;
	
	//cout<<"adresa:="<<addrWidth;
	
	vhdl<<tab<<"process(out_clk1)"<<endl;
	vhdl<<tab<<"variable c1:std_logic_vector(5 downto 0):=(others=> '0');"<<endl;
	vhdl<<tab<<"variable temp:std_logic_vector(1 downto 0):= \"01\";"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk1'event and out_clk1 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if c1 = 0 and temp = \"00\" then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" c1:=(others=>'0');"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" temp:=\"01\";"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<declare("selectionVal1",1)<<" <= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<tab<<declare("out_clk3",1)<<"<= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<"elsif temp= \"01\" then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("out_clk3")<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("selectionVal1")<<" <= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" c1:=c1 + '1';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" temp:=\"10\";"<<endl;
	
	vhdl<<tab<<tab<<tab<<"elsif temp=\"10\" then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("out_clk3")<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("selectionVal1")<<" <= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" c1:=c1 + '1';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" temp:=\"11\";"<<endl;
	
	
	vhdl<<tab<<tab<<tab<<"else"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("out_clk3")<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"temp:=\"00\";"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"c1:=c1+'1';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("selectionVal1")<<" <= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<declare("lowAddress1",6)<<" <= c1;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<"process(out_clk3)"<<endl;
	vhdl<<tab<<"variable c1:std_logic_vector("<<addrWidth - 6 -1<<" downto 0):=(others=> '0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk3'event and out_clk3 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"c1:=c1+'1';"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<declare("highAddress1",addrWidth - 6)<<" <= c1;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	vhdl<<tab<<declare("Address1",addrWidth)<<" <= "<<use("highAddress1")<<" & "<<use("lowAddress1")<<";"<<endl<<endl;
	
	vhdl<<tab<<"process(out_clk3)"<<endl;
	vhdl<<tab<<"variable c2:std_logic_vector(5 downto 0):=(others=> '0');"<<endl;
	vhdl<<tab<<"variable temp:std_logic_vector(1 downto 0):= \"01\";"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk3'event and out_clk3 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"if c2 = 0 and temp = \"00\" then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"c2:=(others=>'0');"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" temp:=\"01\";"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<declare("selectionVal2",1)<<" <= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<tab<<declare("out_clk4",1)<<"<= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<"elsif temp= \"01\" then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("out_clk4")<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("selectionVal2")<<" <= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" c2:=c2 + '1';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" temp:=\"10\";"<<endl;
	
	
	vhdl<<tab<<tab<<tab<<"elsif temp= \"10\" then"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("out_clk4")<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("selectionVal2")<<" <= '1'; "<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" c2:=c2 + '1';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<" temp:=\"11\";"<<endl;
	
	
	
	vhdl<<tab<<tab<<tab<<"else"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("out_clk4")<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"temp:=\"00\";"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<"c2:=c2+'1';"<<endl;
	vhdl<<tab<<tab<<tab<<tab<<use("selectionVal2")<<" <= '0'; "<<endl;
	vhdl<<tab<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<declare("lowAddress2",6)<<" <= c2;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<"process(out_clk4)"<<endl;
	vhdl<<tab<<"variable c2:std_logic_vector("<<addrWidth - 6 -1<<" downto 0):=(others=> '0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if out_clk4'event and out_clk4 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"c2:=c2+'1';"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<declare("highAddress2",addrWidth - 6)<<" <= c2;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	vhdl<<tab<<declare("Address2",addrWidth)<<" <= "<<use("highAddress2")<<" & "<<use("lowAddress2")<<";"<<endl<<endl;
	
	
	
	vhdl<<tab<<"process(Address1,Address2)"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<" if Address1 - Address2 > 1 or Address1 - Address2 < -1 then"<<endl;
	vhdl<<tab<<tab<<tab<<declare("closeAddr",1)<<" <= '0';"<<endl;
	vhdl<<tab<<tab<<"else"<<endl;
	vhdl<<tab<<tab<<tab<<use("closeAddr")<<" <= '1';"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<declare("selectionVal",1)<<" <= "<<use("selectionVal2")<<" or "<<use("selectionVal1")<<" or "<<use("closeAddr")<<";"<<endl;
	
	vhdl<<tab<<"process(out_clk1)"<<endl;
	vhdl<<tab<<"variable temp:std_logic:='0';"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<"if out_clk1'event and out_clk1 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"temp:="<<use("selectionVal")<<";"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<declare("selectionValD",1)<<"<= temp;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	vhdl<<endl;
	
	vhdl<<tab<<declare("selection4Pipeline",1)<<" <= "<<use("selectionVal")<<" or "<<use("selectionValD")<<";"<<endl;
	
	
	//Memories instantiation
	
	 
	
	//cout<<"LSBI:++++ "<<LSBI<<endl;
	memCoordX = new CoordinatesTableX(target,13,LSBI, MSBI,filepath); // perhaps adding a new parameter for number of addresses -> suggestion NO
	memCoordX->changeName(getName()+"memCoordX");	
	oplist.push_back(memCoordX);
	inPortMap  (memCoordX, "X1", "Address1");
	inPortMap  (memCoordX, "X2", "Address2");
	outPortMap (memCoordX, "Y1","dataX1");
	outPortMap (memCoordX, "Y2","dataX2");
	vhdl << instance(memCoordX, "memCoordX");
			
	
	memCoordY = new CoordinatesTableY(target,13,LSBI, MSBI,filepath) ;// perhaps adding a new parameter for number of addresses -> suggestion NO
	memCoordY->changeName(getName()+"memCoordY");	
	oplist.push_back(memCoordY);
	inPortMap  (memCoordY, "X1", "Address1");
	inPortMap  (memCoordY, "X2", "Address2");
	outPortMap (memCoordY, "Y1","dataY1");
	outPortMap (memCoordY, "Y2","dataY2");
	vhdl << instance(memCoordY, "memCoordY");
	
	memCoordZ = new CoordinatesTableZ(target,13,LSBI, MSBI,filepath) ;// perhaps adding a new parameter for number of addresses -> suggestion NO
	memCoordZ->changeName(getName()+"memCoordZ");	
	oplist.push_back(memCoordZ);
	inPortMap  (memCoordZ, "X1", "Address1");
	inPortMap  (memCoordZ, "X2", "Address2");
	outPortMap (memCoordZ, "Y1","dataZ1");
	outPortMap (memCoordZ, "Y2","dataZ2");
	vhdl << instance(memCoordZ, "memCoordZ");
	
	vhdl<<endl;
	vhdl<<tab<<"process(out_clk1)"<<endl;
	vhdl<<tab<<"variable tempX: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempY: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempZ: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<"if  out_clk1'event and out_clk1 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"tempX:="<<use("dataX1")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempY:="<<use("dataY1")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempZ:="<<use("dataZ1")<<";"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("dataD1X1",inputWidth)<<" <= tempX;"<<endl;
	vhdl<<tab<<tab<<declare("dataD1Y1",inputWidth)<<" <= tempY;"<<endl;
	vhdl<<tab<<tab<<declare("dataD1Z1",inputWidth)<<" <= tempZ;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
	
	vhdl<<endl;
	vhdl<<tab<<"process(out_clk1)"<<endl;
	vhdl<<tab<<"variable tempX: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempY: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempZ: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<"if  out_clk1'event and out_clk1 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"tempX:="<<use("dataD1X1")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempY:="<<use("dataD1Y1")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempZ:="<<use("dataD1Z1")<<";"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("dataD2X1",inputWidth)<<" <= tempX;"<<endl;
	vhdl<<tab<<tab<<declare("dataD2Y1",inputWidth)<<" <= tempY;"<<endl;
	vhdl<<tab<<tab<<declare("dataD2Z1",inputWidth)<<" <= tempZ;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
	
	vhdl<<endl;
	vhdl<<tab<<"process(out_clk3)"<<endl;
	vhdl<<tab<<"variable tempX: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempY: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempZ: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<"if  out_clk3'event and out_clk3 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"tempX:="<<use("dataX2")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempY:="<<use("dataY2")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempZ:="<<use("dataZ2")<<";"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("dataD1X2",inputWidth)<<" <= tempX;"<<endl;
	vhdl<<tab<<tab<<declare("dataD1Y2",inputWidth)<<" <= tempY;"<<endl;
	vhdl<<tab<<tab<<declare("dataD1Z2",inputWidth)<<" <= tempZ;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
	
	vhdl<<endl;
	vhdl<<tab<<"process(out_clk3)"<<endl;
	vhdl<<tab<<"variable tempX: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempY: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"variable tempZ: std_logic_vector("<<inputWidth-1<<" downto 0):=(others=>'0');"<<endl;
	vhdl<<tab<<"begin"<<endl;
	vhdl<<tab<<tab<<"if  out_clk3'event and out_clk3 = '1' then"<<endl;
	vhdl<<tab<<tab<<tab<<"tempX:="<<use("dataD1X2")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempY:="<<use("dataD1Y2")<<";"<<endl;
	vhdl<<tab<<tab<<tab<<"tempZ:="<<use("dataD1Z2")<<";"<<endl;
	vhdl<<tab<<tab<<"end if;"<<endl;
	vhdl<<tab<<tab<<declare("dataD2X2",inputWidth)<<" <= tempX;"<<endl;
	vhdl<<tab<<tab<<declare("dataD2Y2",inputWidth)<<" <= tempY;"<<endl;
	vhdl<<tab<<tab<<declare("dataD2Z2",inputWidth)<<" <= tempZ;"<<endl;
	vhdl<<tab<<"end process;"<<endl;
	
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD2X1s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD2X1")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD2Y1s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD2Y1")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD2Z1s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD2Z1")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD1X1s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD1X1")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD1Y1s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD1Y1")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD1Z1s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD1Z1")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD2X2s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD2X2")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD2Y2s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD2Y2")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD2Z2s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD2Z2")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD1X2s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD1X2")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD1Y2s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD1Y2")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	vhdl<<tab<<"with "<<use("selection4Pipeline")<<" select "<<declare("dataD1Z2s",inputWidth)<<" <= "<<endl;
	vhdl<<tab<<tab<<use("dataD1Z2")<<" when '0', "<<endl;
	vhdl<<tab<<tab<<"(others=> '0' ) when others;"<<endl;
	vhdl<<endl;
	
	
	vhdl<<tab<<declare("signal_x0",inputWidth)<<"<= "<<use("dataD2X1s")<<";"<<endl;
	vhdl<<tab<<declare("signal_x1",inputWidth)<<"<= "<<use("dataD1X1s")<<";"<<endl;
	vhdl<<tab<<declare("signal_x2",inputWidth)<<"<= "<<use("dataD2X2s")<<";"<<endl;
	vhdl<<tab<<declare("signal_x3",inputWidth)<<"<= "<<use("dataD1X2s")<<";"<<endl;
	
	vhdl<<tab<<declare("signal_y0",inputWidth)<<"<= "<<use("dataD2Y1s")<<";"<<endl;
	vhdl<<tab<<declare("signal_y1",inputWidth)<<"<= "<<use("dataD1Y1s")<<";"<<endl;
	vhdl<<tab<<declare("signal_y2",inputWidth)<<"<= "<<use("dataD2Y2s")<<";"<<endl;
	vhdl<<tab<<declare("signal_y3",inputWidth)<<"<= "<<use("dataD1Y2s")<<";"<<endl;
	
	vhdl<<tab<<declare("signal_z0",inputWidth)<<"<= "<<use("dataD2Z1s")<<";"<<endl;
	vhdl<<tab<<declare("signal_z1",inputWidth)<<"<= "<<use("dataD1Z1s")<<";"<<endl;
	vhdl<<tab<<declare("signal_z2",inputWidth)<<"<= "<<use("dataD2Z2s")<<";"<<endl;
	vhdl<<tab<<declare("signal_z3",inputWidth)<<"<= "<<use("dataD1Z2s")<<";"<<endl;
	
	
	
	
	//Computing the segments x1-x0 x3-x2 y1-y0 y3-y2 z1-z0 z3-z2
	
	//syncCycleFromSignal("????"); sincronization with memories
	
	//~ addInput("signal_x0",inputWidth);
	//~ addInput("signal_x1",inputWidth);
	//~ addInput("signal_x2",inputWidth);
	//~ addInput("signal_x3",inputWidth);
	
	//~ addInput("signal_y0",inputWidth);
	//~ addInput("signal_y1",inputWidth);
	//~ addInput("signal_y2",inputWidth);
	//~ addInput("signal_y3",inputWidth);
	
	//~ addInput("signal_z0",inputWidth);
	//~ addInput("signal_z1",inputWidth);
	//~ addInput("signal_z2",inputWidth);
	//~ addInput("signal_z3",inputWidth);
	
	/*
	vhdl<<tab<<declare("signal_x0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_x1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_x2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_x3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	
	vhdl<<tab<<declare("signal_y0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_y1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_y2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_y3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	
	vhdl<<tab<<declare("signal_z0",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;  //possible need to add 2 bits for the special bits ; modified to take the appropriate value read from memory
	vhdl<<tab<<declare("signal_z1",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_z2",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	vhdl<<tab<<declare("signal_z3",inputWidth)<<"<= "<<"(others=>'0')"<<";"<<endl;
	*/
	
		
	
	//performing x1-x0
	
	setCycleFromSignal("signal_x0");
	
	vhdl<<tab<<declare("segmentX1mX0",inputWidth)<<" <= "<<use("signal_x1")<<" - "<<use("signal_x0")<<";"<<endl;
	
	//~ //performing x0-x1
	
	//~ setCycleFromSignal("signal_x1");
	
	//~ vhdl<<tab<<declare("segmentX0mX1",inputWidth)<<" <= "<<use("signal_x0")<<" - "<<use("signal_x1")<<";"<<endl;
	
	//performing x0-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("segmentX0mX2",inputWidth)<<" <= "<<use("signal_x0")<<" - "<<use("signal_x2")<<";"<<endl;
	
	
	//performing x3-x2
	
	setCycleFromSignal("signal_x2");
	
	vhdl<<tab<<declare("segmentX3mX2",inputWidth)<<" <= "<<use("signal_x3")<<" - "<<use("signal_x2")<<";"<<endl;
	
	
	//performing x3-x0
	
	setCycleFromSignal("signal_x0");
	
	vhdl<<tab<<declare("segmentX3mX0",inputWidth)<<" <= "<<use("signal_x3")<<" - "<<use("signal_x0")<<";"<<endl;
	
	
	//performing y1-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("segmentY1mY0",inputWidth)<<" <= "<<use("signal_y1")<<" - "<<use("signal_y0")<<";"<<endl;
	
	//~ //performing y0-y1
	
	//~ setCycleFromSignal("signal_y1");
	
	//~ vhdl<<tab<<declare("segmentY0mY1",inputWidth)<<" <= "<<use("signal_y0")<<" - "<<use("signal_y1")<<";"<<endl;
	
	//performing y0-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("segmentY0mY2",inputWidth)<<" <= "<<use("signal_y0")<<" - "<<use("signal_y2")<<";"<<endl;
	
	
	
	//performing y3-y2
	
	setCycleFromSignal("signal_y2");
	
	vhdl<<tab<<declare("segmentY3mY2",inputWidth)<<" <= "<<use("signal_y3")<<" - "<<use("signal_y2")<<";"<<endl;
	
	
	//performing y3-y0
	
	setCycleFromSignal("signal_y0");
	
	vhdl<<tab<<declare("segmentY3mY0",inputWidth)<<" <= "<<use("signal_y3")<<" - "<<use("signal_y0")<<";"<<endl;
	
	
	//performing z1-z0
	
	setCycleFromSignal("signal_z0");

	vhdl<<tab<<declare("segmentZ1mZ0",inputWidth)<<" <= "<<use("signal_z1")<<" - "<<use("signal_z0")<<";"<<endl;

		
	//~ //performing z0-z1
	
	//~ setCycleFromSignal("signal_z1");
	
	//~ vhdl<<tab<<declare("segmentZ0mZ1",inputWidth)<<" <= "<<use("signal_z0")<<" - "<<use("signal_z1")<<";"<<endl;
		
	
	//performing z0-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("segmentZ0mZ2",inputWidth)<<" <= "<<use("signal_z0")<<" - "<<use("signal_z2")<<";"<<endl;
	
	
	//performing z3-z2
	
	setCycleFromSignal("signal_z2");
	
	vhdl<<tab<<declare("segmentZ3mZ2",inputWidth)<<" <= "<<use("signal_z3")<<" - "<<use("signal_z2")<<";"<<endl;
	
	
	//performing z3-z0
	
	setCycleFromSignal("signal_z0");
	
	vhdl<<tab<<declare("segmentZ3mZ0",inputWidth)<<" <= "<<use("signal_z3")<<" - "<<use("signal_z0")<<";"<<endl;
	
	nextCycle();
	
	//referenceCycle1 = getCurrentCycle();
	
	//computing the var1 (x1-x0)*(x3-x2)+(y1-y0)*(y3-y2)+(z1-z0)*(z3-z2)
	
	
		//(x1-x0)*(x3-x2)
	//setCycleFromSignal("signX1mX0v1");
	//setCycle(referenceCycle1);
	
	
	vhdl<<tab<<declare("convertedSegmentXv1",inputWidthSegments)<<" <= "<<use("segmentX1mX0")<<" * "<<use("segmentX3mX2")<<";"<<endl;
	
			
	//(y1-y0)*(y3-y2)
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle1);
	
	
	vhdl<<tab<<declare("convertedSegmentYv1",inputWidthSegments)<<" <= "<<use("segmentY1mY0")<<" * "<<use("segmentY3mY2")<<";"<<endl;
	
	
			
		//(z1-z0)*(z3-z2)
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle1);
	
	
	vhdl<<tab<<declare("convertedSegmentZv1",inputWidthSegments)<<" <= "<<use("segmentZ1mZ0")<<" * "<<use("segmentZ3mZ2")<<";"<<endl;
	
	
		
	nextCycle();
	
	
		//ading the 3 results for segments
	
	vhdl<<tab<<declare("convertedSegmentYv1temp",inputWidthSegments)<<" <= "<<use("convertedSegmentYv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentZv1temp",inputWidthSegments)<<" <= "<<use("convertedSegmentZv1")<<";"<<endl;
	vhdl<<tab<<declare("convertedSegmentXv1temp",inputWidthSegments)<<" <= "<<use("convertedSegmentXv1")<<";"<<endl;
	
	target->setFrequency(((double) internFreq) );
	
	adder4var = new IntNAdder(target,inputWidthSegments,3);
	target->setFrequency(tempFreq);
	
	adder4var->changeName(getName()+"adder4var");	//aici
	oplist.push_back(adder4var);
	inPortMap  (adder4var, "X0", use("convertedSegmentYv1temp"));
	inPortMap  (adder4var, "X1", use("convertedSegmentZv1temp"));
	inPortMap  (adder4var, "X2", use("convertedSegmentXv1temp") );
	inPortMapCst(adder4var,"Cin","'0'");
	outPortMap (adder4var, "R","result4Var1");
	vhdl << instance(adder4var, "adder4var1");
	
	syncCycleFromSignal("result4Var1");
	
	//int signofLSBI=LSBI>=0?(+1):(-1);
	int signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2<<endl;
	
	if(tempFreq>200000000)
		nextCycle();
	
	
	vhdl<<tab<<declare("result4Var1temp",inputWidthSegments)<<" <= "<<use("result4Var1")<<";"<<endl;
	
	target->setFrequency(((double) internFreq) );
	
	convert2FP = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2,1,wE,wF);
	target->setFrequency(tempFreq);
	
	convert2FP->changeName(getName()+"convert2FPv");	//aici
	oplist.push_back(convert2FP);
	inPortMap  (convert2FP, "I", use("result4Var1temp"));
	outPortMap (convert2FP, "O","Var1");
	vhdl << instance(convert2FP, "convert2FPv1");
	
	syncCycleFromSignal("Var1");
	
		
	//computing the var2 sqrt( (x3-x2)^2 +  (y3-y2)^2 +(z3-z2)^2 )
	
		//computing (x3-x2)^2
	
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("sqrX3mX2s",inputWidthSegments)<<" <= "<<use("segmentX3mX2")<<" * "<<use("segmentX3mX2")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrX3mX2ns",(inputWidth-1)*2)<<" <= "<<use("sqrX3mX2s")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
		
	
		//computing (y3-y2)^2
	
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("sqrY3mY2s",inputWidthSegments)<<" <= "<<use("segmentY3mY2")<<" * "<<use("segmentY3mY2")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrY3mY2ns",(inputWidth-1)*2)<<" <= "<<use("sqrY3mY2s")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	
		//computing (z3-z2)^2
	
	setCycleFromSignal("convertedSegmentXv1");
	//setCycle(referenceCycle2);
	
	vhdl<<tab<<declare("sqrZ3mZ2s",inputWidthSegments)<<" <= "<<use("segmentZ3mZ2")<<" * "<<use("segmentZ3mZ2")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrZ3mZ2ns",(inputWidth-1)*2)<<" <= "<<use("sqrZ3mZ2s")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrX3mX2nstemp",(inputWidth-1)*2)<<" <= "<<use("sqrX3mX2ns")<<";"<<endl;
	vhdl<<tab<<declare("sqrY3mY2nstemp",(inputWidth-1)*2)<<" <= "<<use("sqrY3mY2ns")<<";"<<endl;
	vhdl<<tab<<declare("sqrZ3mZ2nstemp",(inputWidth-1)*2)<<" <= "<<use("sqrZ3mZ2ns")<<";"<<endl;
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	target->setFrequency(((double) internFreq) );
	adder4SQRTv = new IntNAdder(target,(inputWidth-1)*2,3);
	target->setFrequency(tempFreq);
	
	adder4SQRTv->changeName(getName()+"adder4SQRTv");	//aici
	oplist.push_back(adder4SQRTv);
	inPortMap  (adder4SQRTv, "X0", use("sqrY3mY2nstemp"));
	inPortMap  (adder4SQRTv, "X1", use("sqrZ3mZ2nstemp"));
	inPortMap  (adder4SQRTv, "X2", use("sqrX3mX2nstemp") );
	inPortMapCst(adder4SQRTv,"Cin","'0'");
	outPortMap (adder4SQRTv, "R","result4SQRTv2");
	vhdl << instance(adder4SQRTv, "adder4SQRTv2");
	
	syncCycleFromSignal("result4SQRTv2");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	target->setFrequency(((double) internFreq) );
	
	convert2FP4sqrtv = new Fix2FP(target,(LSBI)*2,signofMSBI*(abs(MSBI)-1)*2-1,0,wE,wF);
	target->setFrequency(tempFreq);
	
	convert2FP4sqrtv->changeName(getName()+"convert2FP4sqrtv");	//aici
	oplist.push_back(convert2FP4sqrtv);
	inPortMap  (convert2FP4sqrtv, "I", use("result4SQRTv2"));
	outPortMap (convert2FP4sqrtv, "O","fpSQRVar2");
	vhdl << instance(convert2FP4sqrtv, "convert2FP4sqrtv2");
	
	syncCycleFromSignal("fpSQRVar2");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar2temp",wE+wF+3)<<" <= "<<use("fpSQRVar2")<<";"<<endl;
	
	target->setFrequency(((double) internFreq) );
	sqrt4var = new  FPSqrt(target, wE, wF, 1, 0);
	target->setFrequency(tempFreq);
	sqrt4var->changeName(getName()+"sqrt4var");	//aici
	oplist.push_back(sqrt4var);
	inPortMap  (sqrt4var, "X", use("fpSQRVar2temp"));
	outPortMap (sqrt4var, "R","Var2");
	vhdl << instance(sqrt4var, "sqrt4var2");
	
	syncCycleFromSignal("Var2");
	
	
	//computing the var3 sqrt(((x3-x0)+(x0-x1)*t)^2+((y3-y0)+(y0-y1)*t)^2+((z3-z0)+(z0-z1)*t)^2)
	
	
	
	

	//(x3-x0)+(x0-x1)*t
	setCycleFromSignal("convertedSegmentXv1");
	
	target->setFrequency(((double) internFreq) );
	cstMult = new  IntConstMult(target, inputWidth, 14349);
	target->setFrequency(tempFreq);
	cstMult->changeName(getName()+"cstMult");	
	oplist.push_back(cstMult);
	inPortMap  (cstMult, "inX", "segmentX1mX0");
	outPortMap (cstMult, "R","multTXvar3");
	vhdl << instance(cstMult, "cstMult3X");
	 
	syncCycleFromSignal("multTXvar3");
	
	
	vhdl<<tab<<declare("partialConvertedProductSXv3",inputWidth)<<" <= "<<use("multTXvar3")<<of(inputWidth-1 + 14)<<" & "<<use("multTXvar3")<<range(inputWidth-1 + 14,15)<<";"<<endl;
	
	vhdl<<tab<<declare("sumXPartv3",inputWidth)<<" <= "<<use("segmentX3mX0")<<" - "<<use("partialConvertedProductSXv3")<<";"<<endl;

	//nextCycle();
	
	vhdl<<tab<<declare("sqrXsv3",inputWidthSegments)<<" <= "<<use("sumXPartv3")<<" * "<<use("sumXPartv3")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrXnsv3",(inputWidth-1)*2)<<" <= "<<use("sqrXsv3")<<range(inputWidthSegments-2,0)<<";"<<endl;

	

		//(y3-y0)+(y0-y1)*t
	setCycleFromSignal("convertedSegmentXv1");
	
	
	inPortMap  (cstMult, "inX", "segmentY1mY0");
	outPortMap (cstMult, "R","multTYvar3");
	vhdl << instance(cstMult, "cstMult3Y");
	 
	syncCycleFromSignal("multTYvar3");
	
	
	vhdl<<tab<<declare("partialConvertedProductSYv3",inputWidth)<<" <= "<<use("multTYvar3")<<of(inputWidth-1 + 14)<<" & "<<use("multTYvar3")<<range(inputWidth-1 + 14,15)<<";"<<endl;
	
	vhdl<<tab<<declare("sumYPartv3",inputWidth)<<" <= "<<use("segmentY3mY0")<<" - "<<use("partialConvertedProductSYv3")<<";"<<endl;
	
	//nextCycle();
	
	vhdl<<tab<<declare("sqrYsv3",inputWidthSegments)<<" <= "<<use("sumYPartv3")<<" * "<<use("sumYPartv3")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrYnsv3",(inputWidth-1)*2)<<" <= "<<use("sqrYsv3")<<range(inputWidthSegments-2,0)<<";"<<endl;

	

	//(z3-z0)+(z0-z1)*t
	setCycleFromSignal("convertedSegmentXv1");
	
	
	inPortMap  (cstMult, "inX", "segmentZ1mZ0");
	outPortMap (cstMult, "R","multTZvar3");
	vhdl << instance(cstMult, "cstMult3Z");
	 
	syncCycleFromSignal("multTZvar3");
	
	
	vhdl<<tab<<declare("partialConvertedProductSZv3",inputWidth)<<" <= "<<use("multTZvar3")<<of(inputWidth-1 + 14)<<" & "<<use("multTZvar3")<<range(inputWidth-1 + 14,15)<<";"<<endl;
	
	vhdl<<tab<<declare("sumZPartv3",inputWidth)<<" <= "<<use("segmentZ3mZ0")<<" - "<<use("partialConvertedProductSZv3")<<";"<<endl;
	
	//nextCycle();
	
	vhdl<<tab<<declare("sqrZsv3",inputWidthSegments)<<" <= "<<use("sumZPartv3")<<" * "<<use("sumZPartv3")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrZnsv3",(inputWidth-1)*2)<<" <= "<<use("sqrZsv3")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
			
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrXnsv3temp",(inputWidth-1)*2)<<" <= "<<use("sqrXnsv3")<<";"<<endl;
	vhdl<<tab<<declare("sqrYnsv3temp",(inputWidth-1)*2)<<" <= "<<use("sqrYnsv3")<<";"<<endl;
	vhdl<<tab<<declare("sqrZnsv3temp",(inputWidth-1)*2)<<" <= "<<use("sqrZnsv3")<<";"<<endl;
	
	
	
	inPortMap  (adder4SQRTv, "X0", use("sqrYnsv3temp"));
	inPortMap  (adder4SQRTv, "X1", use("sqrZnsv3temp"));
	inPortMap  (adder4SQRTv, "X2", use("sqrXnsv3temp") );
	inPortMapCst(adder4SQRTv,"Cin","'0'");
	outPortMap (adder4SQRTv, "R","result4SQRTv3");
	vhdl << instance(adder4SQRTv, "adder4SQRTv3");
	
	syncCycleFromSignal("result4SQRTv3");
	
	
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	
	//cout<<"new Lsbi:= "<<(LSBI)*2<<"new Msbi:= "<<signofMSBI*(abs(MSBI)-1)*2-1<<endl;
	
	inPortMap  (convert2FP4sqrtv, "I", use("result4SQRTv3"));
	outPortMap (convert2FP4sqrtv, "O","fpSQRVar3");
	vhdl << instance(convert2FP4sqrtv, "convert2FP4sqrtv3");
	
	syncCycleFromSignal("fpSQRVar3");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar3temp",wE+wF+3)<<" <= "<<use("fpSQRVar3")<<";"<<endl;
	
	inPortMap  (sqrt4var, "X", use("fpSQRVar3temp"));
	outPortMap (sqrt4var, "R","Var3");
	vhdl << instance(sqrt4var, "sqrt4var3");
	
	syncCycleFromSignal("Var3");
	

	
	vhdl<<tab<<declare("expVar3shift",wE)<<" <= "<<use("Var3")<<" - "<<" CONV_STD_LOGIC_VECTOR(3,"<<wE<<");"<<endl;
	
	vhdl<<tab<<declare("accVar3",wF+wE+3)<<" <= "<<use("Var3")<<range(wE+wF+3-1,wE+wF)<<" & "<<use("expVar3shift")<<" & "<<use("Var3")<<range(wF-1,0)<<";"<<endl;
	
	
	
		
	
	
	
	//computing the var4 sqrt(((x0-x2)+(x1-x0)*t)^2+((y0-y2)+(y1-y0)*t)^2+((z0-z2)+(z1-z0)*t)^2)
	
	
	
	
	
	
	//((x0-x2)+(x1-x0)*t)^2
	setCycleFromSignal("multTXvar3");
	
	
	vhdl<<tab<<declare("partialConvertedProductSXv4",inputWidth)<<" <= "<<use("multTXvar3")<<of(inputWidth-1 + 14)<<" & "<<use("multTXvar3")<<range(inputWidth-1 + 14,15)<<";"<<endl;
	vhdl<<tab<<declare("sumXv4",inputWidth)<<" <= "<<use("partialConvertedProductSXv4")<<" + "<<use("segmentX0mX2")<<";"<<endl;

	//nextCycle();
	
	vhdl<<tab<<declare("sqrXsv4",inputWidthSegments)<<" <= "<<use("sumXv4")<<" * "<<use("sumXv4")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrXnsv4",(inputWidth-1)*2)<<" <= "<<use("sqrXsv4")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	

	//((y0-y2)+(y1-y0)*t)^2
	setCycleFromSignal("multTYvar3");
	
	vhdl<<tab<<declare("partialConvertedProductSYv4",inputWidth)<<" <= "<<use("multTYvar3")<<of(inputWidth-1 + 14)<<" & "<<use("multTYvar3")<<range(inputWidth-1 + 14,15)<<";"<<endl;
	vhdl<<tab<<declare("sumYv4",inputWidth)<<" <= "<<use("partialConvertedProductSYv4")<<" + "<<use("segmentY0mY2")<<";"<<endl;
	
	//nextCycle();
	
	vhdl<<tab<<declare("sqrYsv4",inputWidthSegments)<<" <= "<<use("sumYv4")<<" * "<<use("sumYv4")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrYnsv4",(inputWidth-1)*2)<<" <= "<<use("sqrYsv4")<<range(inputWidthSegments-2,0)<<";"<<endl;	
	
		
		
	//((z0-z2)+(z1-z0)*t)^2
	setCycleFromSignal("multTZvar3");
	
	vhdl<<tab<<declare("partialConvertedProductSZv4",inputWidth)<<" <= "<<use("multTZvar3")<<of(inputWidth-1 + 14)<<" & "<<use("multTZvar3")<<range(inputWidth-1 + 14,15)<<";"<<endl;
	vhdl<<tab<<declare("sumZv4",inputWidth)<<" <= "<<use("partialConvertedProductSZv4")<<" + "<<use("segmentZ0mZ2")<<";"<<endl;
	
	//nextCycle();
	
	vhdl<<tab<<declare("sqrZsv4",inputWidthSegments)<<" <= "<<use("sumZv4")<<" * "<<use("sumZv4")<<";"<<endl;
	
	vhdl<<tab<<declare("sqrZnsv4",(inputWidth-1)*2)<<" <= "<<use("sqrZsv4")<<range(inputWidthSegments-2,0)<<";"<<endl;
	
	
	
	
	//ading the 3 results for segments to be feed to sqrt
	
	nextCycle();
	
	vhdl<<tab<<declare("sqrXnsv4temp",(inputWidth-1)*2)<<" <= "<<use("sqrXnsv4")<<";"<<endl;
	vhdl<<tab<<declare("sqrYnsv4temp",(inputWidth-1)*2)<<" <= "<<use("sqrYnsv4")<<";"<<endl;
	vhdl<<tab<<declare("sqrZnsv4temp",(inputWidth-1)*2)<<" <= "<<use("sqrZnsv4")<<";"<<endl;
	
	
	inPortMap  (adder4SQRTv, "X0", use("sqrYnsv4temp"));
	inPortMap  (adder4SQRTv, "X1", use("sqrZnsv4temp"));
	inPortMap  (adder4SQRTv, "X2", use("sqrXnsv4temp") );
	inPortMapCst(adder4SQRTv,"Cin","'0'");
	outPortMap (adder4SQRTv, "R","result4SQRTv4");
	vhdl << instance(adder4SQRTv, "adder4SQRTv4");
	
	syncCycleFromSignal("result4SQRTv4");
	
	
	signofMSBI=MSBI>=0?(+1):(-1);
	


	inPortMap  (convert2FP4sqrtv, "I", use("result4SQRTv4"));
	outPortMap (convert2FP4sqrtv, "O","fpSQRVar4");
	vhdl << instance(convert2FP4sqrtv, "convert2FP4sqrtv4");
	
	syncCycleFromSignal("fpSQRVar4");
	
	if(target->frequencyMHz()>=250)
		nextCycle();
	
	vhdl<<tab<<declare("fpSQRVar4temp",wE+wF+3)<<" <= "<<use("fpSQRVar4")<<";"<<endl;
	
	inPortMap  (sqrt4var, "X", use("fpSQRVar4temp"));
	outPortMap (sqrt4var, "R","Var4");
	vhdl << instance(sqrt4var, "sqrt4var4");
	
	syncCycleFromSignal("Var4");



	vhdl<<tab<<declare("expVar4shift",wE)<<" <= "<<use("Var4")<<" - "<<" CONV_STD_LOGIC_VECTOR(3,"<<wE<<");"<<endl;
	
	vhdl<<tab<<declare("accVar4",wF+wE+3)<<" <= "<<use("Var4")<<range(wE+wF+3-1,wE+wF)<<" & "<<use("expVar4shift")<<" & "<<use("Var4")<<range(wF-1,0)<<";"<<endl;
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//computing the var5 ( ((x0-x2)+(x1-x0)*t)*(x3-x2) +((y0-y2)+(y1-y0)*t)*(y3-y2)+((z0-z2)+(z1-z0)*t)*(z3-z2))
	
	
		//((x0-x2)+(x1-x0)*t)*(x3-x2)
	setCycleFromSignal("convertedSegmentXv1");
	
	vhdl<<tab<<declare("shiftX1var5",inputWidth)<<" <= "<<use("segmentX1mX0")<<of(inputWidth-1)<<" & "<<use("segmentX1mX0")<<range(inputWidth-1,1)<<";"<<endl;
	vhdl<<tab<<declare("shiftX4var5",inputWidth)<<" <= "<<use("segmentX1mX0")<<of(inputWidth-1)<<" & "<<use("segmentX1mX0")<<of(inputWidth-1)<<" & "<<use("segmentX1mX0")<<of(inputWidth-1)<<" & "<<use("segmentX1mX0")<<of(inputWidth-1)<<" & "<<use("segmentX1mX0")<<range(inputWidth-1,4)<<";"<<endl;
	vhdl<<tab<<declare("shiftX4var5neg",inputWidth)<<" <= not("<<use("shiftX4var5")<<");"<<endl;
	
	
	target->setFrequency(((double) internFreq) );
	
	adder4coordvar5 = new IntNAdder(target,inputWidth,3);
	target->setFrequency(tempFreq);
	
	adder4coordvar5->changeName(getName()+"adder4coordvar5");	
	oplist.push_back(adder4coordvar5);
	inPortMap  (adder4coordvar5, "X0", "shiftX1var5");
	inPortMap  (adder4coordvar5, "X1", "shiftX4var5neg");
	inPortMap  (adder4coordvar5, "X2", "segmentX0mX2" );
	inPortMapCst(adder4coordvar5,"Cin","'1'");
	outPortMap (adder4coordvar5, "R","sumX4var5");
	vhdl << instance(adder4coordvar5, "adder4coordvar5X");
	
	syncCycleFromSignal("sumX4var5");
	
	
	
	//nextCycle();
	
	vhdl<<tab<<declare("convertedProductSXv5",inputWidthSegments )<<" <= "<<use("segmentX3mX2")<<" * "<<use("sumX4var5")<<";"<<endl;
	
	
	
		//((y0-y2)+(y1-y0)*t)*(y3-y2)
	setCycleFromSignal("convertedSegmentXv1");
	
	vhdl<<tab<<declare("shiftY1var5",inputWidth)<<" <= "<<use("segmentY1mY0")<<of(inputWidth-1)<<" & "<<use("segmentY1mY0")<<range(inputWidth-1,1)<<";"<<endl;
	vhdl<<tab<<declare("shiftY4var5",inputWidth)<<" <= "<<use("segmentY1mY0")<<of(inputWidth-1)<<" & "<<use("segmentY1mY0")<<of(inputWidth-1)<<" & "<<use("segmentY1mY0")<<of(inputWidth-1)<<" & "<<use("segmentY1mY0")<<of(inputWidth-1)<<" & "<<use("segmentY1mY0")<<range(inputWidth-1,4)<<";"<<endl;
	vhdl<<tab<<declare("shiftY4var5neg",inputWidth)<<" <= not("<<use("shiftY4var5")<<");"<<endl;
	
	
	inPortMap  (adder4coordvar5, "X0", "shiftY1var5");
	inPortMap  (adder4coordvar5, "X1", "shiftY4var5neg");
	inPortMap  (adder4coordvar5, "X2", "segmentY0mY2" );
	inPortMapCst(adder4coordvar5,"Cin","'1'");
	outPortMap (adder4coordvar5, "R","sumY4var5");
	vhdl << instance(adder4coordvar5, "adder4coordvar5Y");
	
	syncCycleFromSignal("sumY4var5");
	
	
	
	//nextCycle();
	
	vhdl<<tab<<declare("convertedProductSYv5",inputWidthSegments)<<" <= "<<use("segmentY3mY2")<<" * "<<use("sumY4var5")<<";"<<endl;
	
	
		//((z0-z2)+(z1-z0)*t)*(z3-z2)
	setCycleFromSignal("convertedSegmentXv1");
	
	
	vhdl<<tab<<declare("shiftZ1var5",inputWidth)<<" <= "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<" & "<<use("segmentZ1mZ0")<<range(inputWidth-1,1)<<";"<<endl;
	vhdl<<tab<<declare("shiftZ4var5",inputWidth)<<" <= "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<" & "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<" & "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<" & "<<use("segmentZ1mZ0")<<of(inputWidth-1)<<" & "<<use("segmentZ1mZ0")<<range(inputWidth-1,4)<<";"<<endl;
	vhdl<<tab<<declare("shiftZ4var5neg",inputWidth)<<" <= not("<<use("shiftZ4var5")<<");"<<endl;
	
	
	inPortMap  (adder4coordvar5, "X0", "shiftZ1var5");
	inPortMap  (adder4coordvar5, "X1", "shiftZ4var5neg");
	inPortMap  (adder4coordvar5, "X2", "segmentZ0mZ2" );
	inPortMapCst(adder4coordvar5,"Cin","'1'");
	outPortMap (adder4coordvar5, "R","sumZ4var5");
	vhdl << instance(adder4coordvar5, "adder4coordvar5Z");
	
	syncCycleFromSignal("sumZ4var5");
	
	
	
	//nextCycle();
	
	vhdl<<tab<<declare("convertedProductSZv5",inputWidthSegments )<<" <= "<<use("segmentZ3mZ2")<<" * "<<use("sumZ4var5")<<";"<<endl;
	
	
	
	
	
	
	//Adding the 3 results and the caries that are neded to finalize the transformation of the numbers in 2's complement
	
	nextCycle();
	
	vhdl<<tab<<declare("partialConvertedProductSYv5temp",inputWidthSegments)<<" <= "<<use("convertedProductSYv5")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSZv5temp",inputWidthSegments)<<" <= "<<use("convertedProductSZv5")<<";"<<endl;
	vhdl<<tab<<declare("partialConvertedProductSXv5temp",inputWidthSegments)<<" <= "<<use("convertedProductSXv5")<<";"<<endl;
	
	
	
	
	inPortMap  (adder4var, "X0", use("partialConvertedProductSYv5temp"));
	inPortMap  (adder4var, "X1", use("partialConvertedProductSZv5temp"));
	inPortMap  (adder4var, "X2", use("partialConvertedProductSXv5temp") );
	inPortMapCst(adder4var,"Cin","'0'");
	outPortMap (adder4var, "R","result4Var5");
	vhdl << instance(adder4var, "adder4var5");
	
	syncCycleFromSignal("result4Var5");
	
	//vhdl<<tab<<declare("result4Var5",inputWidthSegments)<<" <= "<<use("result4Var5noshift")<<of(inputWidthSegments-1)<<" & "<<use("result4Var5noshift")<<of(inputWidthSegments-1)<<" & "<<use("result4Var5noshift")<<of(inputWidthSegments-1)<<" & "<<use("result4Var5noshift")<<range(inputWidthSegments-1,3)<<";"<<endl;
	
	
	inPortMap  (convert2FP, "I", use("result4Var5"));
	outPortMap (convert2FP, "O","accVar5");
	vhdl << instance(convert2FP, "convert2FPv5");
	
	syncCycleFromSignal("accVar5");
	
	
	
	
	
	
	
	
	
	
	//de aici in jos toate componentele ar trebui sa functioneze pe clockul out_clk1 (este o diviziune cu 9+1 a clk)
	
	setCycleFromSignal("accVar4");
	
	
	vhdl<<tab<<declare("Var1temp",wE+wF+3)<<" <= "<<use("Var1")<<";"<<endl;
	vhdl<<tab<<declare("Var2temp1",wE+wF+3)<<" <= "<<use("Var2")<<";"<<endl;
	
	target->setFrequency(((double) internFreq) );
	div4Log =new FPDiv(target, wE, wF);
	target->setFrequency(tempFreq);
	div4Log->changeName(getName()+"div4Acc");
	oplist.push_back(div4Log);
	inPortMap  (div4Log, "X", use("Var1temp"));
	inPortMap (div4Log, "Y",use("Var2temp1"));
	outPortMap (div4Log, "R","var1divvar2");
	vhdl << instance(div4Log, "var1divvar24acc","out_clk1","rst");
	//vhdl << instance(div4Log, "var1divvar24acc");
	
	syncCycleFromSignal("var1divvar2");
	
	
	setCycleFromSignal("accVar4");
		
	
	vhdl<<tab<<declare("Var2temp",wE+wF+3)<<" <= "<<use("Var2")<<";"<<endl;
	vhdl<<tab<<declare("accVar3temp",wE+wF+3)<<" <= "<<use("accVar3")<<";"<<endl;
	
	
	target->setFrequency(((double) internFreq) );
	acc4var = new FPAdder(target, wE, wF, wE, wF, wE, wF);	
	target->setFrequency(tempFreq);
	
	acc4var->changeName(getName()+"accumulator4var");	
	oplist.push_back(acc4var);
	
	inPortMap  (acc4var, "X", use("Var2temp"));
	inPortMap (acc4var, "Y",use("accVar3temp"));
	outPortMap (acc4var, "R","var3pvar2");
	vhdl << instance(acc4var, "var3plusvar2");
	
	syncCycleFromSignal("var3pvar2");
	
	nextCycle();
	
	setCycleFromSignal("accVar4");
	
	vhdl<<tab<<declare("accVar5temp",wE+wF+3)<<" <= "<<use("accVar5")<<";"<<endl;
	vhdl<<tab<<declare("Var2temp2",wE+wF+3)<<" <= "<<use("Var2")<<";"<<endl;
	
	inPortMap  (div4Log, "X", use("accVar5temp"));
	inPortMap (div4Log, "Y",use("Var2temp2"));
	outPortMap (div4Log, "R","var5divvar2");
	vhdl << instance(div4Log, "var5divvar24log","out_clk1","rst");
	//vhdl << instance(div4Log, "var5divvar24log");
	
	syncCycleFromSignal("var5divvar2");
	
	vhdl<<tab<<declare("var3pvar2temp",wE+wF+3)<<" <= "<<use("var3pvar2")<<";"<<endl;
	vhdl<<tab<<declare("minusvar5divvar2",wE+wF+3)<<" <= "<<use("var5divvar2")<<range(wE+wF+2,wE+wF+1)<<" & "<<"( not ("<<use("var5divvar2")<<of(wE+wF)<<")) & "<<use("var5divvar2")<<range(wE+wF-1,0)<<";"<<endl;	
	
	inPortMap  (acc4var, "X", use("var3pvar2temp"));
	inPortMap (acc4var, "Y",use("minusvar5divvar2"));
	outPortMap (acc4var, "R","numerator4Log");
	vhdl << instance(acc4var, "numerator4LogAdder");
	
	syncCycleFromSignal("numerator4Log");
	
	
	setCycleFromSignal("var5divvar2");
	
	vhdl<<tab<<declare("accVar4temp",wE+wF+3)<<" <= "<<use("accVar4")<<";"<<endl;
	
	
	inPortMap  (acc4var, "X", use("accVar4temp"));
	inPortMap (acc4var, "Y",use("minusvar5divvar2"));
	outPortMap (acc4var, "R","denominator4Log");
	vhdl << instance(acc4var, "denominator4LogAdder");
	
	syncCycleFromSignal("denominator4Log");
	
	
	inPortMap  (div4Log, "X", use("numerator4Log"));
	inPortMap (div4Log, "Y",use("denominator4Log"));
	outPortMap (div4Log, "R","result4Log");
	vhdl << instance(div4Log, "div4log","out_clk1","rst");
	//vhdl << instance(div4Log, "div4log");
	
	syncCycleFromSignal("result4Log");
	
	

	
	target->setFrequency(((double) internFreq) );
	log4Acc = new FPLog(target, wE, wF);
	target->setFrequency(tempFreq);
	log4Acc->changeName(getName()+"log4Acc");
	oplist.push_back(log4Acc);
	inPortMap  (log4Acc, "X", use("result4Log"));
	outPortMap (log4Acc, "R","resultLog");
	vhdl << instance(log4Acc, "log4acc","out_clk1","rst");
	//vhdl << instance(log4Acc, "log4acc");
	
	
	syncCycleFromSignal("resultLog");
	
	vhdl<<tab<<declare("var1divvar2temp",wE+wF+3)<<" <= "<<use("var1divvar2")<<";"<<endl;
	
	//target->setNotPipelined();
	target->setFrequency(((double) internFreq) );
	mult4Acc = new FPMultiplier(target, wE, wF, wE, wF, wE, wF, 1);
	target->setFrequency(tempFreq);
	//target->setPipelined();
	mult4Acc->changeName(getName()+"mult4Acc");
	oplist.push_back(mult4Acc);
	inPortMap  (mult4Acc, "X", use("var1divvar2temp"));
	inPortMap (mult4Acc, "Y",use("resultLog"));
	outPortMap (mult4Acc, "R","value4LongAcc");
	vhdl << instance(mult4Acc, "mult4acc");
	
	
	syncCycleFromSignal("value4LongAcc");
	
	nextCycle();
	
	vhdl<<tab<<declare("value4LongAcctemp",wE+wF+3)<<" <= "<<use("value4LongAcc")<<";"<<endl;
	
	target->setFrequency(((double) internFreq) );
	finalAcc = new LongAcc(target, wE, wF, MaxMSBO, LSBO, MSBO);
	target->setFrequency(tempFreq);
	finalAcc->changeName(getName()+"finalAcc");
	oplist.push_back(finalAcc);
	inPortMap  (finalAcc, "X", use("value4LongAcctemp"));
	outPortMap (finalAcc, "A","finalResult");
	outPortMap (finalAcc, "XOverflow","XOverflow");
	outPortMap (finalAcc, "XUnderflow","XUnderflow");
	outPortMap (finalAcc, "AccOverflow","AccOverflow");
	vhdl << instance(finalAcc, "finalAcc","out_clk1","rst");
	//vhdl << instance(finalAcc, "finalAcc");
		
	syncCycleFromSignal("finalResult");
		
	
	vhdl<<tab<<"O <= "<<use("finalResult")<<range(outputWidth-1,0)<<";"<<endl;
	
	//vhdl<<tab<<"O<="<<use("Var1")<<range(outputWidth-1,0)<<" or "<<use("accVar4")<<range(outputWidth-1,0)<<" or "<<use("numerator4Log")<<range(outputWidth-1,0)<<";"<<endl;

	}

void CoilInductance::outputVHDL(std::ostream& o, std::string name) {
  
	licence(o);
	o << "library ieee; " << endl;
	o << "use ieee.std_logic_1164.all;" << endl;
	o << "use ieee.std_logic_arith.all;" << endl;
	o << "use ieee.std_logic_signed.all;" << endl;
	o << "library work;" << endl;
	outputVHDLEntity(o);
	newArchitecture(o,name);
	o << buildVHDLComponentDeclarations();	
	o << buildVHDLSignalDeclarations();
	beginArchitecture(o);		
	o<<buildVHDLRegisters();
	o << vhdl.str();
	endArchitecture(o);

}
	
	
CoilInductance::~CoilInductance() {
}

void CoilInductance::emulate(TestCase * tc)
{
}

void CoilInductance::buildStandardTestCases(TestCaseList* tcl){
	
}

int CoilInductance::addressLength()
{
	
	int L;
	int NrSVert;
	int **config;
	int nrTurns;
	float out_radius;
	float turn_radius;
	float insulation;
	int N;
	
ifstream indata; // indata is like cin
	indata.open(filepath);
	
	
	indata>>L>>NrSVert;
	   
	 config = new int*[L];
	   for(int j=0;j<L;j++)
		config[j]= new int[NrSVert];
	 

	   nrTurns=0;
	   for(int j=0;j<L;j++)
		for(int i=0;i<NrSVert;i++)
			{indata>>config[j][i];
			nrTurns+=config[j][i];
			}
	
	indata>>out_radius>>turn_radius>>insulation>>N;

	indata.close();
	
	return ((int) intlog2(N*nrTurns));
}

