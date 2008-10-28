// Copyright blabla

#include "Cotran.hpp"
#include "../utils.hpp"
#include <cmath>
#include "LNSAdd.hpp"

using namespace std;
extern vector<Operator*> oplist;	// ???

Cotran::Cotran(Target * target, int wE, int wF, int j, int wECotran, int o) :
	Operator(target), wE(wE), wF(wF), j(j), wECotran(wECotran)
{
	if(j < 0)
		select_j();
	if(wECotran < 0)
		wECotran = wE;
	
	ostringstream name;
	name<<"Cotran_"<< wE <<"_"<< wF << "_" << j << "_" << wECotran; 
	uniqueName_ = name.str();
	
	name.str("");
	name<<"Cotran_F1_"<< wF <<"_"<< j;
	f1_name  = name.str();
	
	name.str("");
	name<<"Cotran_F2_"<< wF <<"_"<< j;
	f2_name  = name.str();

	name.str("");
	name<<"Cotran_F3_"<< wF <<"_"<< j;
	f3_name  = name.str();
	
//	wEssZero = min(wE, log2(wF+1)+1)+wF;
	wEssZero = min(wE, intlog2(wF)) + wF;
	
	// Translated back from VHDL...
	DBMaxInput = (int)rint(double(1 << wF) * log((1.0 - (pow(2.0, (-pow(2.0, -wF)))))) / log(2.));;
	
	addInput("Z", wE + wF + 1);
	addInput("IsSub");

	addOutput("SBDB", wE + wF + 1);
	
	// Output signals manually...

	// Delayed construction, after computation of j...
	f1 = auto_ptr<CotranF1Table>(new CotranF1Table(wF, j, wECotran));
	f2 = auto_ptr<CotranF2Table>(new CotranF2Table(wF, j));
	f3 = auto_ptr<CotranF3Table>(new CotranF3Table(wF, j));
	
	sb = new LNSAdd(target, wE, wF, o);	// TODO: + guard bits!
	
	oplist.push_back(sb);
}

Cotran::~Cotran()
{
}

void Cotran::select_j()
{
	j = wF / 2 + 2;
}

void Cotran::outputVHDL(std::ostream& o, std::string name)
{
	licence(o,"Sylvain Collange, Jesus Garcia (2005-2008)");
	
	// Output tables
	f1->output(o, f1_name);
	f2->output(o, f2_name);
	f3->output(o, f3_name);

	//Operator::StdLibs(o);
	o
	<< "library ieee;\n"
	<< "use ieee.std_logic_1164.all;\n"
	<< "use ieee.numeric_std.all;\n"
	<< "\n";
	
	outputVHDLEntity(o);
	newArchitecture(o,name);	

	//output_vhdl_signal_declarations(o);
	o
	<< tab << "constant F : positive := " << wF << ";\n"
	<< tab << "constant K : positive := " << wE << ";\n"
	<< tab << "constant wEssZero : positive := " << wEssZero << ";\n"
//	<< tab << "constant wCotran : positive := " << min(wEssZero, wF + wECotran) << ";\n"
	<< tab << "constant wBreak : positive := " << j << ";\n"
	<< tab << "constant DB_Max_Input : std_logic_vector := std_logic_vector(to_signed(" << DBMaxInput << ", K+F+1));\n";

	
	o
	<< tab << "signal SelMuxC, SelMuxB : std_logic_vector (1 downto 0);\n"
	<< tab << "signal Special, IsEssZero : std_logic;\n"
	<< tab << "signal Zh     : std_logic_vector (" << f1->wIn + j - 1 << " downto wBreak);\n"
	<< tab << "signal Zl     : std_logic_vector (wBreak-1 downto 0);\n"
	<< tab << "signal F1_v   : std_logic_vector (K+F downto 0);\n"
	<< tab << "signal F2_v   : std_logic_vector (K+F downto 0);\n"
	<< tab << "signal Zsum   : std_logic_vector (K+F downto 0);\n"
	<< tab << "signal Zdif   : std_logic_vector (K+F downto 0);\n"
	<< tab << "signal ZEssZero : std_logic_vector (K+F downto 0);\n"
	<< tab << "signal Zfinal : std_logic_vector (K+F downto 0);\n"
	<< tab << "signal DB0    : std_logic_vector (K+F downto 0);\n"
//	<< tab << "signal R0     : std_logic_vector (K+F downto 0);\n";
	<< tab << "signal SB0 : std_logic_vector(F-1 downto 0);\n"
	<< tab << "signal SBArg : std_logic_vector(K+F downto 0);\n"
	<< tab << "signal SBPos : std_logic_vector(F downto 0);\n"
	<< tab << "signal SB1 : std_logic_vector(F downto 0);\n"
	<< tab << "signal SBArgSign : std_logic;\n"
	<< tab << "signal SBEssZero : std_logic;\n";


	f1->outputComponent(o, f1_name);
	f2->outputComponent(o, f2_name);
	f3->outputComponent(o, f3_name);
	
	sb->outputVHDLComponent(o);

	beginArchitecture(o);

	o
	<< tab << "Zh <= Z(" << f1->wIn + j - 1 << " downto wBreak);\n"
	<< tab << "Zl <= Z(wBreak-1 downto 0);\n"
	<< tab << "\n"
	<< tab << "IsEssZero <= '1' when Z(K+F downto wBreak) < DB_Max_Input(K+F downto wBreak) else '0';\n"
	<< tab << "\n"
	<< tab << "Special <=\n"
	<< tab << "  '1' when Zh = (" << f1->wIn + j - 1 << " downto wBreak => '1') else\n"
	<< tab << "  '0';\n"
	<< tab << "\n"
//	<< tab << "with IsEssZero or Special select\n"
//	<< tab << "	R0 <=	RMin when '0',\n"
//	<< tab << "			RMax when others;\n"
//	<< tab << "\n"

	<< tab << "f1 : " << f1_name << "\n"
	<< tab << "  port map (\n"
	<< tab << "    x => Zh,\n"
	<< tab << "    y => F1_v(" << f1->wOut - 1 << " downto 0));\n"
	<< tab << "\n";
	
	if(wE + wF + 1 > f1->wOut) {
		o << tab << "F1_v(K+F downto " << f1->wOut << " ) <= (K+F downto " << f1->wOut << " => F1_v(" << f1->wOut - 1 << "));\n";
	}
	
	o
	<< tab << "f2 : " << f2_name << "\n"
	<< tab << "  port map (\n"
	<< tab << "    x => Zl,\n"
	<< tab << "    y => F2_v(" << f2->wOut - 1 << " downto 0));\n"
	<< tab << "\n";
	if(wE + wF + 1 > f2->wOut) {
		o << tab << "F2_v(K+F downto " << f2->wOut << " ) <= (K+F downto " << f2->wOut << " => F2_v(" << f2->wOut - 1 << "));\n";
	}

	o
	<< tab << "f3 : " << f3_name << "\n"
	<< tab << "  port map (\n"
	<< tab << "    x => Z(" << j+1 << " downto 0),\n"
	<< tab << "    y => SBPos);\n";

	o
	<< tab << "\n"
	<< tab << "SelMuxB <= IsSub & (IsEssZero or Special);\n"
	<< tab << "\n"
	<< tab << "SelMuxC <= (IsEssZero or (not IsSub)) & Special;\n";
	

	o	
	<< tab << "Zdif <= std_logic_vector(signed(F2_v) - signed(F1_v) - signed(Z));\n"
	<< tab << "Zsum  <= Z;\n"
	<< tab << "ZEssZero <= '1' & (K+F-1 downto 0 => '0');\n"
	<< tab << "with SelMuxB select                        -- MUX for the address to sb()\n"
	<< tab << "  Zfinal <=\n"
	<< tab << "  Zdif  when \"10\",\n"
	<< tab << "  ZEssZero when \"11\",\n"
	<< tab << "  Zsum  when others;\n"
	<< tab << "\n"
	<< tab << "SBArg <= Zfinal;\n"
	<< tab << "\n"
	<< tab << "with SelMuxC select\n"
	<< tab << "  DB0 <=\n"
	<< tab << "  F1_v                   when \"00\",   -- subtraction, not special case\n"
	<< tab << "  F2_v                   when \"01\",   -- subtraction, special case (not EZ)\n"
	<< tab << "  (F+K downto 0 => '0') when others;  -- addition, or subtraction and EZ\n";
//	<< tab << "DBR <= std_logic_vector(signed(DB0(K+F) & DB0) + signed(R0(K+F) & R0));\n";

	o
	<< tab << "sb : " << sb->getOperatorName() << "\n"
	<< tab << "  port map (\n"
	<< tab << "    x => SBArg(" << sb->wE + sb->wF - 1 << " downto 0),\n"
	<< tab << "    r => SB0);\n";

	o
	<< tab << "SBEssZero <= 	'1' when SBArg(K+F downto wEssZero) /= (K+F downto wEssZero => '1') else\n"
	<< tab << "             	'0';\n"
	<< tab << "\n"
	<< tab << "SBArgSign <= SBArg(K+F);\n"
	<< tab << "\n"
	<< tab << "-- special cases : essential zero, zero, positive\n"
	<< tab << "SB1 <=	SBPos						when SBArgSign = '0' else\n"
	<< tab << "		(F downto 0 => '0')		when SBEssZero = '1' else\n"
	<< tab << "		'0' & SB0;\n";

	o
	<< tab << "SBDB <= std_logic_vector(signed(DB0) + signed('0' & SB1));\n";
	
	o<< "end architecture;" << endl << endl;
}

