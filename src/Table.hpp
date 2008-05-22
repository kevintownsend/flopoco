#ifndef TABLE_HPP
#define TABLE_HPP
#include <gmpxx.h>

// A basic hardware table, with its VHDL output. 

// If the input to your table are negative, etc, or if you want to
// define errors, or... then derive a class from this one.


class Table
{
 public:

	int wIn;

	int wOut;

	int minIn; // usually 0

	int maxIn; // usually 2^wIn-1
	
	Table(int _wIn, int _wOut, int _minIn=0, int _maxIn=-1) : 
		wIn(_wIn), wOut(_wOut), minIn(_minIn), maxIn(_maxIn)
	{
		if(maxIn==-1) maxIn=(1<<wIn)-1;
		if(minIn<0) {
			cerr<<"ERROR in Table::Table, minIn<0\n";
			exit(EXIT_FAILURE);
		}
		if(maxIn>=(1<<wIn)) {
			cerr<<"ERROR in Table::Table, maxIn too large\n";
			exit(EXIT_FAILURE);
		}
		if((minIn==0) && (maxIn==(1<<wIn)-1)) 
			full=true;
		else
			full=false;
	}

	virtual ~Table() {};

	// An integer function of [0, 2^wIn-1]  in [0, 2^wOut-1] 
	// It may be defined only on [minIn, int maxIn]

	virtual mpz_class function(int x) =0;


	void outputComponent(ostream& o, string name, int component_id);
	void outputComponent(ostream& o, string name);


	void output(ostream& o, string name, int component_id);
	void output(ostream& o, string name);

	// This function should be overridden by an implementation of Table
	virtual int    double2input(double x);

	// This function should be overridden by an implementation of Table
	virtual double input2double(int x);

	// This function should be overridden by an implementation of Table
	virtual  mpz_class double2output(double x);

	// This function should be overridden by an implementation of Table
	virtual double output2double(mpz_class x);
	
	int size_in_LUTs();
 private:
	bool full; // true if there is no "don't care" inputs, i.e. minIn=0 and maxIn=2^wIn-1
	void output_lines(ostream& o);
};


#endif
