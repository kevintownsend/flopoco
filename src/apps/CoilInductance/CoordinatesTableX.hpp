#include "../../Table.hpp"


class CoordinatesTableX : public Table
{
	public:
	
	CoordinatesTableX(Target* target, int wIn, int LSBI,int MSBI);
	
	~CoordinatesTableX();
  	
	mpz_class function(int x);
	
	
	int    double2input(double x);
	
	double input2double(int x);
	
	mpz_class double2output(double x);
	
	double output2double(mpz_class x);
	
	double maxMulOut;
	double minMulOut;
	
private:

	/** The MSB for the input */
	int MSBI; 
	/** The LSB for the input */
	int LSBI; 
	/** The width of the exponent for the output R */
	int adrWidth; 
	/**		*/
	int wOutm;

};