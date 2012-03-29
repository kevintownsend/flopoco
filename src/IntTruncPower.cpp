#include "IntTruncPower.hpp"
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
// too long a typename
using std::tr1::unordered_map;
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
	std::vector<ProductBitIR>::iterator j = i + n;
	for (; j != this->data.end(); i++,j++) {
		*i = *j;
	}
	i = this->data.begin();
	j = i+n;
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
//incomplete
ProductIR ProductIR::operator* (const ProductIR& rhs)
{
	size_t n = this->data.size(), m = rhs.data.size();
	ProductIR res (n+m, this->msb + rhs.msb);
	ProductIR tmp (n+m, this->msb + rhs.msb);
	for (int i = 0; i < n; i++) {
	}
	return res;
}

