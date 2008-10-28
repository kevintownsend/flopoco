#ifndef COTRAN_TABLES_HPP
#define COTRAN_TABLES_HPP

#include "../Table.hpp"
#include "../Operator.hpp"
#include <cmath>
#include <gmpxx.h>
#include <memory>

using namespace std;

struct TableOp : Operator
{
	// Acquires ownership of table
	TableOp(Target * target, Table /*const*/ * table) :
		Operator(target),
		table(table)
	{
		// TODO: use better name
		ostringstream name;
		name<<"Table_"<< table->wIn <<"_"<< table->wOut; 
		uniqueName_ = name.str();

	}

	virtual void outputVHDL(std::ostream& o, std::string name) {
		table->output(o, name);
	}

	virtual void outputVHDLComponent(std::ostream& o, std::string name) {
		table->outputComponent(o, name);
	}
		
private:
	auto_ptr<Table> table;
};

struct FXTable : Table
{
	FXTable(int wIn, int wOut, int minIn=0, int maxIn=-1) :
		Table(wIn, wOut, minIn, maxIn) {}
		
	virtual  mpz_class double2output(double x);
	virtual double input2double(int x);

};

struct CotranF1Table : Table
{
	CotranF1Table(int wF, int j, int wE);
		
	virtual ~CotranF1Table() {}

	virtual mpz_class function(int x);

private:
	static int esszero(int wF, int j, int wE);
	static int addrLen(int wF, int j, int wE);
	static int dataLen(int wF, int j, int wE);
	
	int wF;
	int j;
	int k;
	double dh;
	int wE;

};

struct CotranF2Table : Table
{
	CotranF2Table(int wF, int j);
		
	virtual ~CotranF2Table() {}

	virtual mpz_class function(int x);

private:
	static int addrLen(int wF, int j);
	static int dataLen(int wF, int j);
	
	int wF;
	int j;
	int k;
	double dh;

};

struct CotranF3Table : Table
{
	CotranF3Table(int wF, int j);
		
	virtual ~CotranF3Table() {}

	virtual mpz_class function(int x);

private:
	static int addrLen(int wF, int j);
	static int dataLen(int wF, int j);
	static int begin(int wF, int j);
	static int end(int wF, int j);
	double SbArg(int z);
	
	int wF;
	int j;
	int k;
	double dh;

};



#endif
