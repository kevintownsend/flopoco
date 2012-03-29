#include "IntTruncPower.hpp"
#include <iostream>

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

// too long a typename
using std::tr1::unordered_map;

// first we need to supply a hash function
template<class MonomialOfBits>
size_t std::tr1::hash<MonomialOfBits>::operator() (MonomialOfBits m) const
{
	return std::tr1::hash<std::vector<bool> >::hash (m.data);
}
ProductBitIR::ProductBitIR (const ProductBit& rhs)
	:data (unordered_map<MonomialOfBits,int>())
{
	std::list<MonomialOfBits>::const_iterator it
		= rhs.data.begin();
	for (; it != rhs.data.end(); it++) {
		data[*it] = getTimes (*it) + 1;
	}
}
int ProductBitIR::getTimes (const MonomialOfBits& e) const
{
	unordered_map<MonomialOfBits, int>::const_iterator it
		= data.find (e);
	if (it == data.end()) {
		return 0;
	}
	return it->second;
}
ProductBitIR& ProductBitIR::operator+= (const ProductBitIR& rhs)
{
	unordered_map<MonomialOfBits, int>::const_iterator it
		= rhs.data.begin();
	for (; it != rhs.data.end(); it++) {
		this->data[it->first] =
			this->getTimes(it->first) + it->second;
	}
	return *this;
}
ProductBitIR ProductBitIR::operator* (const ProductBitIR& rhs)
{
	ProductBitIR res;
	unordered_map<MonomialOfBits, int>::const_iterator i
		= this->data.begin();
	unordered_map<MonomialOfBits, int>::const_iterator j
		= res.data.begin();
	for (; i != this->data.end(); i++) {
		for (; j != res.data.end(); j++) {
			MonomialOfBits prod = i->first * j->first;
			res.data[prod] =
				res.getTimes(prod) + (i->second * j->second);
		}
	}
	return res;
}
std::ostream& operator<< (std::ostream& o, const ProductBitIR& pbi)
{
	bool cont = false;
	std::tr1::unordered_map<MonomialOfBits, int>::const_iterator it =
		pbi.data.begin();
	for (; it != pbi.data.end(); it++) {
		if (cont)
			o << " + ";
		if (it->second) {
			// can print "+ -" in some cases
			o << it->second;
			o << ' ';
			o << it->first;
			cont = true;
		}
	}
	// and if there's no term...
	if (!cont)
		o << '0';
	return o;
}
//src range must include dst range
ProductIR& ProductIR::operator+= (const ProductIR& rhs)
{
	if (this->msb < rhs.msb)
		throw "msb mismatch (ProductIR::operator+=)";
	if (this->msb - this->data.size() > rhs.msb - rhs.data.size())
		throw "lsb mismatch (ProductIR::operator+=)";
	std::vector<ProductBitIR>::const_iterator j = rhs.data.begin();
	std::vector<ProductBitIR>::iterator i = this->data.begin();
	//align the 2 iterators
	i += (this->msb - rhs.msb);
	for (; j != rhs.data.end(); i++,j++) {
		*i += *j;
	}
	return *this;
}
ProductIR& ProductIR::operator>>= (int n)
{
	if (n < 0)
		throw "invalid argument (ProductIR::operator>>=)";
	std::vector<ProductBitIR>::iterator i = this->data.begin();
	std::vector<ProductBitIR>::const_iterator j =
		this->data.begin() + n;
	for (; j != this->data.end(); i++,j++) {
		*i = *j;
	}
	i = this->data.begin();
	j = this->data.begin() + n;
	for (; i != j; i++) {
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
ProductIR ProductIR::operator* (const ProductIR& rhs)
{
	size_t n = this->data.size(), m = rhs.data.size();
	ProductIR res (n+m, this->msb + rhs.msb);
	for (int i = 0; i < n; i++) {
		ProductIR tmp (n+m, this->msb + rhs.msb);
		tmp += rhs;
		tmp >>= i;
		tmp *= this->data[i];
		res += tmp;
	}
	return res;
}
std::ostream& operator<< (std::ostream& o, const ProductIR& pi)
{
	int i = pi.msb;
	std::vector<ProductBitIR>::const_iterator it = pi.data.begin();
	bool cont = false;
	for (; it != pi.data.end(); it++, i--) {
		if (cont)
			o << " + ";
		o << '(' << (*it) << ") * 1b" << i;
		cont = true;
	}
	return o;
}

#include <cstdlib>
#include <iostream>
int main ()
{
	MonomialOfBits m (10);
	for (int i = 0; i < 10; i++)
		m.data[i] = ((std::rand() & 0x1) == 0);
	std::tr1::unordered_map<MonomialOfBits, int> map;
	map[m] = 1;
	ProductIR a(5, 0);
	ProductBitIR pbi;
	pbi.data = map;
	a.data[0] = pbi;
	ProductIR b(a*a);
	std::cout << pbi << std::endl << pbi*pbi << std::endl;
	std::cout << a << std::endl;
	std::cout << b << std::endl;
	return 0;
}

