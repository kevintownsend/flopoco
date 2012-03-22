
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>


#include "../FloPoCo.hpp" 
 
#include "CordicSinCos.hpp"


using namespace std;
using namespace flopoco;


void addOperator(vector<Operator*> &oplist, Operator *op) {
	
	oplist.push_back(op);
}

Operator* createTestBench(vector<Operator*> &oplist, Operator* op, Target* target, int n, bool tofile){
	Operator* tb = new TestBench(target, op, n, tofile);
	
	cerr << "> TestBench for " << op->getName()<<endl;
	addOperator(oplist, tb);
	cerr << "To run the simulation using ModelSim, type the following in 'vsim -c':" <<endl;
	cerr << tab << "vdel -all -lib work" <<endl;
	cerr << tab << "vlib work" <<endl;
	cerr << tab << "vcom " << "CordicSinCos.vhdl" <<endl;
	cerr << tab << "vsim " << tb->getName() <<endl;
	cerr << tab << "add wave -r *" <<endl;
	cerr << tab << "run " << ((TestBench*)tb)->getSimulationTime() <<"ns" << endl;
	cerr << "To run the simulation using gHDL, type the following in a shell prompt:" <<endl;
	cerr << tab << "ghdl -a --ieee=syntbsys -fexplicit "<< "FPAddSub.vhdl" <<endl;
	cerr << tab << "ghdl -e --ieee=syntbsys -fexplicit " << tb->getName() <<endl;
	cerr << tab << "ghdl -r --ieee=syntbsys " << tb->getName() << " --vcd=" << tb->getName() << ".vcd" <<endl;
	cerr << tab << "gtkwave " << tb->getName() << ".vcd" << endl;	
	
	return tb;
}




int main(int argc, char* argv[] )
{
	Target* target;
	std::string filename;
	int frequency, tests;
	
	vector<Operator*> oplist;

	int w = 32;
	
	CordicSinCos* op;
	
	try{
		cout << "Processing arguments...\n";
		if(argc != 4){
			cout << "Please specify (1)the name of the output file, (2)the target frequency(0 means no pipeline) (3)size of testbench file." << endl;
			exit(1);
		}else{
			filename = argv[1];
			filename.append(".vhdl");
			frequency = atoi(argv[2]);
			tests = atoi(argv[3]);
		}
		
		cout << "Creating target...\n";
		target = new Virtex4();
		if(frequency == 0)
			target->setNotPipelined();
		if (frequency>1 && frequency<10000) {
			target->setFrequency(1e6*(double)frequency);
		}
		cout << "Creating unit under test...\n";
		op = new CordicSinCos(target, w);
		addOperator(oplist, op);
		if(tests != 0){
			cout << "Creating Testbench...\n";
			createTestBench(oplist, op, target, tests, true);
		}
	}catch(std::string str){
		cout << str << endl;
		exit(1);
	}
	
	cout << "Writing to file...\n";
	ofstream file;
	file.open(filename.c_str(), ios::out);
	
	Operator::outputVHDLToFile(oplist, file);

	file.close();
	
	cout << "Finished.\n";
	cerr << endl<<"Final report:"<<endl;
	op->outputFinalReport(0);
	
	return 0;
}



