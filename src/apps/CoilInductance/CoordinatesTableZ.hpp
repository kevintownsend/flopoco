#include "../../DualTable.hpp"


class CoordinatesTableZ : public DualTable
{
	public:
	
	CoordinatesTableZ(Target* target, int wIn, int LSBI,int MSBI,char *filename);
	
	~CoordinatesTableZ();
  	
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
	/** configuration file path from which to initialize the parameters of the coil*/
	char *filepath;
	/** number of coils	*/
	int L;
	/** number of vertical stages	*/
	int NrSVert;
	/** configuration */
	//int  config [L][NrSVert];
	int **config;
	/** number of turns */
	int nrTurns;
	/** outer radius */
	float out_radius;
	/* radius of a turn */
	float turn_radius;
	/**	insulation	*/
	float insulation;
	/** number of points to divide a turn */
	int N;
	

	void initParameters();
	void readParams();
	
	mpfr_t out_rad_m,turn_rad_m,insu_m,nr_m;
	mpfr_t fi;
	mpfr_t h;
	mpfr_t displacementZ[5],alfa[5];//,displacementZ[5];
};

