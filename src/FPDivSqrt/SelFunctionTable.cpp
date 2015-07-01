#include <gmpxx.h>
#include "SelFunctionTable.hpp"
#include <math.h>

using namespace std;
namespace flopoco
{
	SelFunctionTable::SelFunctionTable(Target* target, float dmin, float dmax, int nbd, int nbw, int digit, int base, int wIn, int wOut) : Table(target, wIn, wOut),
																														Dmin(dmin),
																														Dmax(dmax),
																														nbBitD(nbd),
																														nbBitW(nbw),
																														digitSet(digit),
																														radix(base)
	{
		ro = ((float)digitSet)/(radix-1);
		ostringstream name;
		srcFileName="SelFunctionTable";
		name << "SelFunctionTable";
		setName(name.str());
	}

	SelFunctionTable::~SelFunctionTable() {}

	mpz_class SelFunctionTable::function(int x)
	{
		mpz_class result;

		//TODO : Parameterized the table

        int w = (x>>nbBitD);	 //separating w and d (together in x)
        int d = x - (w<<nbBitD);
        int wneg = (w >> (nbBitW-1)) << (nbBitW-1); //getting the sign bit
        w -= 2*wneg; //switching from two's complement to decimal representation, more convenient for upcoming computation

//		if((w>=16 && w<=19) || (w==20 && d<=1) || (w==21 && d==0))
//		{
//			result = mpz_class("01001", 2);
//		}
//		else if ((w==20 && d>=2) || (w==21 && d>=1) || (w==22 && d<=1))
//		{
//			result = mpz_class("11010", 2);
//		}
//		else if ((w==22 && d>=2) || (w==23) || (w==24 && d<=1))
//		{
//			result = mpz_class("11011", 2);
//		}
//		else if ((w==24 && d>=2) || (w==25) || (w==26 && d==0))
//		{
//			result = mpz_class("11100", 2);
//		}
//		else if ((w==26 && d>=1) || (w==27))
//		{
//			result = mpz_class("11101", 2);
//		}
//		else if (w==28)
//		{
//			result = mpz_class("11110", 2);
//		}
//		else if (w==29 || w==30)
//		{
//			result = mpz_class("11111", 2);
//		}
//		else if (w==0 || w==31)
//		{
//			result = mpz_class("00000", 2);
//		}
//		else if (w==1 || w==2)
//		{
//			result = mpz_class("00001", 2);
//		}
//		else if (w==3)
//		{
//			result = mpz_class("00010", 2);
//		}
//		else if ((w==5 && d>=1) || (w==4))
//		{
//			result = mpz_class("00011", 2);
//		}
//		else if ((w==5 && d==0) || (w==6) || (w==7 && d>=2))
//		{
//			result = mpz_class("00100", 2);
//		}
//		else if ((w==7 && d<=1) || (w==8) || (w==9 && d>=2))
//		{
//			result = mpz_class("00101", 2);
//		}
//		else if ((w==9 && d<=1) || (w==10 && d>=1) || (w==11 && d>=2))
//		{
//			result = mpz_class("00110", 2);
//		}
//		else if ((w==10 && d==0) || (w==11 && d<=1) || (w>=12 && w<=15))
//		{
//			result = mpz_class("10111", 2);
//		}

		int decimalResult;
		int nbBitK = ceil(log2(digitSet)+1); //Nb of bit for the entire part of w
		float realw = w/pow(2, nbBitW-nbBitK);

		float hPitch = pow(2, -nbBitD)*(Dmax-Dmin);
		float vPitch = pow(2, -nbBitW+nbBitK);

		float uSlope, lSlope;
		int uCorrecter, lCorrecter;


		for(int k = -digitSet; k <= digitSet; k++)
		{
			uSlope = k+ro;
			lSlope = k-ro;

			uCorrecter = (uSlope<0 ? 1 : 0);
			lCorrecter = (lSlope<0 ? 0 : 1);

			float wMax = ((d+uCorrecter)*hPitch + Dmin) * uSlope;
			float wMin = ((d+lCorrecter)*hPitch + Dmin) * lSlope;

			if((realw+vPitch <= wMax && realw >= wMin) || (k == digitSet && realw >= wMin) || (k == -digitSet && realw+vPitch <= wMax))
			{
				//cout << realw << " " << d << " qui donne qi = "  << k << endl;
				//cout << " parce que evidemment on est entre " << wMin << " et " << wMax << endl;
				decimalResult = k;
				break;
			}
		}


		if(decimalResult < 0)
		{
			decimalResult+=(pow(2, nbBitK)); //switch to two's complement
			if(decimalResult != 9 && radix > 4)
				decimalResult+=pow(2, nbBitK); //qa is negative
		}
		else if(decimalResult == 7)
		{
			decimalResult+=pow(2, nbBitK); //qa is negative
		}

		result = mpz_class(decimalResult);
		return result;
	}
}
