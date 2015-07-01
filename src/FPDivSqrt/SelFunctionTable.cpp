#include <gmpxx.h>
#include "SelFunctionTable.hpp"

namespace flopoco
{
	SelFunctionTable::SelFunctionTable(Target* target) : Table(target, 7, 5)
	{
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

        int w = (x>>2);
        int d = x - (w<<2);

		if((w>=16 && w<=19) || (w==20 && d<=1) || (w==21 && d==0))
		{
			result = mpz_class("01001", 2);
		}
		else if ((w==20 && d>=2) || (w==21 && d>=1) || (w==22 && d<=1))
		{
			result = mpz_class("11010", 2);
		}
		else if ((w==22 && d>=2) || (w==23) || (w==24 && d<=1))
		{
			result = mpz_class("11011", 2);
		}
		else if ((w==24 && d>=2) || (w==25) || (w==26 && d==0))
		{
			result = mpz_class("11100", 2);
		}
		else if ((w==26 && d>=1) || (w==27))
		{
			result = mpz_class("11101", 2);
		}
		else if (w==28)
		{
			result = mpz_class("11110", 2);
		}
		else if (w==29 || w==30)
		{
			result = mpz_class("11111", 2);
		}
		else if (w==0 || w==31)
		{
			result = mpz_class("00000", 2);
		}
		else if (w==1 || w==2)
		{
			result = mpz_class("00001", 2);
		}
		else if (w==3)
		{
			result = mpz_class("00010", 2);
		}
		else if ((w==5 && d>=1) || (w==4))
		{
			result = mpz_class("00011", 2);
		}
		else if ((w==5 && d==0) || (w==6) || (w==7 && d>=2))
		{
			result = mpz_class("00100", 2);
		}
		else if ((w==7 && d<=1) || (w==8) || (w==9 && d>=2))
		{
			result = mpz_class("00101", 2);
		}
		else if ((w==9 && d<=1) || (w==10 && d>=1) || (w==11 && d>=2))
		{
			result = mpz_class("00110", 2);
		}
		else if ((w==10 && d==0) || (w==11 && d<=1) || (w>=12 && w<=15))
		{
			result = mpz_class("10111", 2);
		}

		return result;
	}
}
