/*
 * Helper table for LNS Cotransformation
 *
 * Author : Sylvain Collange
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

#include "CotranTables.hpp"
#include <cmath>
#include <climits>

using namespace std;

mpz_class FXTable::double2output(double x)
{
	// Result in [0, 1[
	// Map to [0, 2^wOut[
	if(x < 0 || x >= 1) {
		throw "FXTable::double2output : argument out of range";
	}
	
	x = ldexp(x, wOut);
	// which rounding??
	return mpz_class(x);
}

double FXTable::input2double(int x)
{
	// x in [0, 2^wIn[
	// Map to [0, 1[
	return ldexp((double)x, -wIn);
}


////////////////////////////////////////////////

double db2(double z)
{
	// Assumes target precision is small enough for double precision to provide correct rounding...
	return log(abs(1. - pow(2, z)))/log(2.);
}

int CotranF1Table::addrLen(int wF, int j, int wE)
{
	int k = min(int(log(wF+1)/log(2)+1), wE);
	return wF + k - j;
}

int CotranF1Table::dataLen(int wF, int j, int wE)
{
	int k2 = int(log(wF+1)/log(2)+1);
	return wF + k2 + 1;
}

int CotranF1Table::esszero(int wF, int j, int wE)
{
	double dh = 1. / (1 << (wF - j));
	int k = min(int(log(wF+1)/log(2)+1), wE);
	// negative number
	double esszero1 = log(1. - pow(2., -(pow(2., -wF)))) / log(2.0) - dh;
//	cout << "esszero = " << esszero1 << endl;
	
	// 2's-complement
	return max(0, int(floor((1.+esszero1 / (1 << k)) * (1 << addrLen(wF, j, wE)))));
}

CotranF1Table::CotranF1Table(int wF, int j, int wE) :
	Table(addrLen(wF, j, wE),
	        dataLen(wF, j, wE),
	        esszero(wF, j, wE),
	        (1 << addrLen(wF, j, wE)) - 2),
	wF(wF), j(j), k(int(log(wF+1)/log(2)+1)), wE(wE)
{
	dh = 1. / (1 << (wF - j));
}


mpz_class CotranF1Table::function(int x)
{
	double y = ((double(x) / (1 << wIn)) - 1) * (1 << k);
	
		
	int decval = (int)rint(db2(-y-dh) * (1 << wF));

	// keep only low-order bits (2's-complement with wOut bits)
	decval &= ((1 << wOut) - 1);
//	cout << "x=" << x << ", y=" << y << ", val=" << decval << endl;
	return mpz_class(decval);
}


/////////////////////////////////////////////////////

int CotranF2Table::addrLen(int wF, int j)
{
	return j;
}

int CotranF2Table::dataLen(int wF, int j)
{
	return wF + int(log(wF+1)/log(2)+1) + 1;
}


CotranF2Table::CotranF2Table(int wF, int j) :
	Table(addrLen(wF, j),
	        dataLen(wF, j)),
	wF(wF), j(j), k(int(log(wF+1)/log(2)+1))
{
	dh = 1. / (1 << (wF - j));
}


mpz_class CotranF2Table::function(int x)
{
	double y = double(x) / (1 << wF);
	
	int decval = (int)rint(db2(y-dh) * (1 << wF));

	// keep only low-order bits (2's-complement with wOut bits)
	decval &= ((1 << wOut) - 1);
//	cout << "x=" << x << ", y=" << y << ", val=" << decval << endl;
	return mpz_class(decval);
}

///////////////////////////////////////////////////


int CotranF3Table::addrLen(int wF, int j)
{
	return j + 2;
}

int CotranF3Table::dataLen(int wF, int j)
{
	return wF + 1;
}

int two_compl(int i, int w)
{
	return i & ((1 << w) - 1);
}

int sign_ext(int i, int w)
{
	int shift_amount = (sizeof(int) * CHAR_BIT - w);
	return (i << shift_amount) >> shift_amount; 
}

int CotranF3Table::begin(int wF, int j)
{
	double dh = 1. / (1 << (wF - j));
	return two_compl(floor(-2 * (dh * (1 << wF))), j + 2);
}

int CotranF3Table::end(int wF, int j)
{
	double dh = 1. / (1 << (wF - j));
	return two_compl(ceil(-log(2 * pow(2., dh) - 1)/log(2) * (1 << wF)), j + 2);
}

CotranF3Table::CotranF3Table(int wF, int j) :
	Table(addrLen(wF, j),
	      dataLen(wF, j),
	      begin(wF, j),
	      end(wF, j)),
	wF(wF), j(j), k(int(log(wF+1)/log(2)+1))
{
	dh = 1. / (1 << (wF - j));
}

double sb2(double z)
{
	return log(1 + pow(2., z))/log(2);
}

double CotranF3Table::SbArg(int z)
{
	// Direct translation from a Perl script from a few years ago.
	// No more idea how it works, but supposed to work.
	// TODO: make sure it actually works
	double zh = (floor((double)z / (1 << j))) * (1 << j) / (1 << wF);
	double zl = (double)(z % (1 << j)) / (1 << wF);
	double f1 = db2(-zh - dh);
	double f2 = db2(zl - dh);
	return f2 - ((double)z / (1 << wF) + f1);
}

mpz_class CotranF3Table::function(int x)
{
	int y = sign_ext(x, wIn);
	double sbarg = SbArg(y);
	int decval = (int)rint(sb2(sbarg) * (1 << wF));

	// keep only low-order bits (2's-complement with wOut bits)
//	decval &= ((1 << wOut) - 1);
//	cout << "x=" << x << ", val=" << decval << endl;
	return mpz_class(decval);
}

