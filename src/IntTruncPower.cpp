#include "IntTruncPower.hpp"
#include <iostream>
#include <map>

using std::cout;
using std::endl;

MonomialOfBits MonomialOfBits::operator* (const MonomialOfBits& rhs) const
{
	if (this->n != rhs.n)
		throw "sizes don't match (MonomialOfBits::operator+)";
	MonomialOfBits res (*this);
	std::vector<bool>::iterator i = res.data.begin();
	std::vector<bool>::const_iterator j = rhs.data.begin();
	for (;j != rhs.data.end (); i++,j++) {
		*i = *i || *j;
	}
	return res;
}
bool MonomialOfBits::operator< (const MonomialOfBits& rhs) const
{
	if (this->n != rhs.n)
		throw "sizes don't match (MonomialOfBits::operator<)";
	std::vector<bool>::const_iterator
		i = this->data.begin(), j = rhs.data.begin();
	for (; i != this->data.end(); i++, j++) {
		if (*i == *j) continue;
		return (*i < *j);
	}
	// these are equal
	return false;
}
std::ostream& operator<< (std::ostream& o, const MonomialOfBits& m)
{
	bool cont = false;
	size_t i;
	for (i = 0; i < m.data.size(); i++) {
		if (m.data[i]) {
			if (cont)
				o << '*';
			o << "x_" << i;
			cont = true;
		}
	}
	// and if there is no x_i at all...
	if (!cont) {
		o << '1';
	}
	return o;
}

// the ProductBitIR must have a 0 or 1 coeff (see also the throw)
ProductBit::ProductBit (const ProductBitIR& rhs)
	:data (std::list<MonomialOfBits>())
{
	std::map<MonomialOfBits, int>::const_iterator it = rhs.data.begin();
	for (; it != rhs.data.end(); it++) {
		switch (it->second) {
			case 0:
				break;
			case 1:
				data.push_back (it->first);
				break;
			default:
				throw "The data must be from a ProductIR"
				      "on which simplify() has been called\n";
		}
	}
}

ProductBitIR::ProductBitIR (const ProductBit& rhs)
	:data (std::map<MonomialOfBits,int>())
{
	std::list<MonomialOfBits>::const_iterator it
		= rhs.data.begin();
	for (; it != rhs.data.end(); it++) {
		addToCoeff (*it, 1);
	}
}
int ProductBitIR::getTimes (const MonomialOfBits& e) const
{
	std::map<MonomialOfBits, int>::const_iterator it
		= data.find (e);
	if (it == data.end()) {
		return 0;
	}
	return it->second;
}
void ProductBitIR::addToCoeff (const MonomialOfBits& e, int coeff)
{
	std::map<MonomialOfBits, int>::iterator it = data.find (e);
	if (it == data.end()) {
		data[e] = coeff;
	} else {
		it->second += coeff;
	}
}
ProductBitIR& ProductBitIR::operator+= (const ProductBitIR& rhs)
{
	std::map<MonomialOfBits, int>::const_iterator it
		= rhs.data.begin();
	for (; it != rhs.data.end(); it++) {
		this->addToCoeff (it->first, it->second);
	}
	return *this;
}
ProductBitIR ProductBitIR::operator* (const ProductBitIR& rhs) const
{
	ProductBitIR res;
	std::map<MonomialOfBits, int>::const_iterator i
		= this->data.begin();
	std::map<MonomialOfBits, int>::const_iterator j;
	for (; i != this->data.end(); i++) {
		j = rhs.data.begin();
		for (; j != rhs.data.end(); j++) {
			MonomialOfBits prod = i->first * j->first;
			res.addToCoeff (prod, i->second * j->second);
		}
	}
	return res;
}
std::ostream& operator<< (std::ostream& o, const ProductBitIR& pbi)
{
	bool cont = false;
	std::map<MonomialOfBits, int>::const_iterator it =
		pbi.data.begin();
	for (; it != pbi.data.end(); it++) {
		if (it->second) {
			if (cont)
				o << " + ";
			// can print "+ -" in some cases
			if (it->second != 1) {
				o << it->second;
				o << ' ';
			}
			o << it->first;
			cont = true;
		}
	}
	// and if there's no term...
	if (!cont)
		o << '0';
	return o;
}
ProductIR::ProductIR (const Product& rhs)
	:msb(rhs.msb),
	data(std::vector<ProductBitIR> (rhs.data.size(),ProductBitIR()))
{
	std::vector<ProductBit>::const_iterator i = rhs.data.begin();
	std::vector<ProductBitIR>::iterator j = data.begin();
	for (; i != rhs.data.end(); i++, j++) {
		*j = ProductBitIR (*i);
	}
}
//src range must include dst range
ProductIR& ProductIR::operator+= (const ProductIR& rhs)
{
	if (this->msb < rhs.msb)
		throw "msb mismatch (ProductIR::operator+=)";
	if (this->msb - this->data.size() > rhs.msb - rhs.data.size())
		throw "lsb mismatch (ProductIR::operator+=)";
	// we use reverse iterators because we know MSBs more readily
	// than LSBs
	std::vector<ProductBitIR>::const_reverse_iterator j
		= rhs.data.rbegin();
	std::vector<ProductBitIR>::reverse_iterator i = this->data.rbegin();
	//align the 2 iterators
	i += (this->msb - rhs.msb);
	for (; j != rhs.data.rend(); i++,j++) {
		*i += *j;
	}
	return *this;
}
ProductIR& ProductIR::operator>>= (int n)
{
	if (n < 0)
		throw "invalid argument (ProductIR::operator>>=)";
	if (n == 0)
		return *this;
	// arithmetically it is a >>=, so since we use an el vector
	// we have to do a <<= in memory
	std::vector<ProductBitIR>::iterator i = this->data.begin();
	std::vector<ProductBitIR>::iterator j = i + n;
	for (; j != this->data.end(); i++,j++) {
		*i = *j;
	}
	for (; i != this->data.end(); i++) {
		*i = ProductBitIR();
	}
	return *this;
}
ProductIR& ProductIR::operator*= (const ProductBitIR& rhs)
{
	std::vector<ProductBitIR>::iterator i = this->data.begin();
	for (; i != this->data.end(); i++) {
		//no ProductBitIR::operator*=
		*i = *i * rhs;
	}
	return *this;
}
std::ostream& operator<< (std::ostream& o, const ProductIR& pi)
{
	int i = pi.msb;
	// we begin by msb so we have to use a reverse iterator
	std::vector<ProductBitIR>::const_reverse_iterator it
		= pi.data.rbegin();
	bool cont = false;
	for (; it != pi.data.rend(); it++, i--) {
		if (cont)
			o << " + ";
		o << '(' << (*it) << ") * 1b" << i;
		cont = true;
	}
	return o;
}
ProductIR ProductIR::operator* (const ProductIR& rhs) const
{
	size_t n = this->data.size(), m = rhs.data.size();
	ProductIR res (n+m-1, this->msb + rhs.msb);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			// since product has exactly the right size
			// it doesn't matter if we make the product
			// in counter-endian mode
			res.data[i+j] += (this->data[i] * rhs.data[j]);
		}
	}
	return res;
}
void ProductIR::simplify (void)
{
	// we begin by lsb: non-rev iterator
	std::vector<ProductBitIR>::iterator i = data.begin(), i_tmp;
	int coeff; // to avoid messing beetween goto and declarations
	for (; i != data.end(); i++) {
		std::map<MonomialOfBits, int>::iterator j = i->data.begin();
		for (; j != i->data.end(); j++) {
			loop_begin:
			coeff = j->second;
			if (coeff < 0) {
				throw "negative coefficient, shouldn't"
				      "happen (ProductIR::simplify)";
			}
			if (coeff > 1) {
				i_tmp = i+1;
				if (i_tmp == data.end()) {
					msb++;
					data.push_back (ProductBitIR());
					// now i and i_tmp are both invalid
					// so regenerate them
					i_tmp = data.end() - 1;
					i = i_tmp - 1;
					// apparently j is also destroyed
					j = i->data.begin();
					// not sure continue; does what I want
					if (j == i->data.end())
						break;
					else
						goto loop_begin;
				}
				i_tmp->addToCoeff (j->first, coeff >> 1);
				j->second = coeff & 0x1;
			}
		}
	}
}
// must call ProductIR::simplify() on rhs before
Product::Product (const ProductIR& rhs)
	:data(std::vector<ProductBit>(rhs.data.size(), ProductBit()))
	 ,msb (rhs.msb)
{
	std::vector<ProductBitIR>::const_iterator i = rhs.data.begin();
	std::vector<ProductBit>::iterator j = data.begin();
	for (; i != rhs.data.end(); i++, j++) {
		*j = ProductBit (*i);
	}
}

#include <cstdlib>
#include <iostream>
int main ()
{
	MonomialOfBits m (3);
	ProductIR a(5, 0);
	for (int i0 = 0; i0 < 5; i0++) {
		for (int i = 0; i < 3; i++)
			m.data[i] = ((std::rand() & 0x1) == 0);
		std::map<MonomialOfBits, int> map;
		map[m] = 1;
		ProductBitIR pbi;
		pbi.data = map;
		a.data[i0] = pbi;
	}
	//std::cout << m << std::endl << m*m << std::endl;
	//std::cout << pbi << std::endl << pbi*pbi << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "a*a = " << a*a << std::endl;
	ProductIR b = (a*a)*a;
	//b.simplify();
	std::cout << "(a*a)*a = " << b << std::endl;
	b = a*(a*a);
	//b.simplify();
	std::cout << "a*(a*a) = " << b << std::endl;
	return 0;
}

